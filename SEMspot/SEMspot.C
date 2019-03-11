/*---------------------------------------------------------------------------*\
License
    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Author: Alex Skillen.       alex.skillen@stfc.ac.uk

Reference. Accuracy and Efficiency Improvements in Synthetic Eddy Methods. A. Skillen A. Revell T. Craft IJHFF (2016) DOI: 10.1016/j.ijheatfluidflow.2016.09.008
\*---------------------------------------------------------------------------*/

#include "SEMspot.H"
#include "turbulentInletSEMFvPatchField.H" 

namespace Foam
{



SEMspot::SEMspot
(
    SEMBase* pf 
)
:
    pf_(pf)
{
}


SEMspot::SEMspot
(
    SEMspot* spt,
    SEMBase* pf 
)
:
    pf_(pf),
    origin_(spt->origin_),
    sigma_(spt->sigma_),
    epsilon_(spt->epsilon_),
    u_(spt->u_),
    nearest_(spt->nearest_),
    donorProcN_(spt->donorProcN_),
    residenceTime_(spt->residenceTime_)
{
}   
//--

void SEMspot::initialise(const bool setToFace )
{
    if( Pstream::master() )
    { 
        pf_->ranGen().randomise01( origin_ );
        origin_ = cmptMultiply( origin_, pf_->semBox().span() );
        origin_ += pf_->semBox().min();
    }
    Pstream::scatter( origin_ );

    interpolatePatchToSpot(); 

    vector nn( -u_ );
    nn /= mag(nn);
    u_=mag(u_)*-nn;

    this->projectBack( nn, setToFace );

    scalar projectedDist = 
    ( donorProcN_ == Pstream::myProcNo() ) ?
    mag( (pf_->patch().Cf()[nearest_] - origin_)&nn ):
    GREAT;

    reduce( projectedDist, minOp<scalar>() );

    if( setToFace )
    {      
        residenceTime_ = pf_->db().time().value() + 2.0*projectedDist/mag(u_);
    }
    else
    {
        residenceTime_ = pf_->db().time().value() + (projectedDist + mag((sigma_&nn)) )/mag(u_);
    }

    if( Pstream::master() )
    {
        //also randomise the sign of the spot (i.e. epsilon_ = 1 or -1)
        pf_->ranGen().randomise01<vector>( epsilon_ );
        for( int i=0; i<3; i++ ) 
        {
            epsilon_[i] = round(epsilon_[i]);
        }

        epsilon_ *= 2;
        epsilon_ -= pTraits<vector>::one;
    }

    Pstream::scatter( epsilon_ );
    Pstream::scatter( origin_ );
    Pstream::scatter( residenceTime_ );
}

//--



void SEMspot::interpolatePatchToSpot()
{
    //get face centres
    const vectorField& fc = pf_->patch().Cf();


    //find nearest face centre to 'this'
    List<scalar> distCur( Pstream::nProcs(), 0.0 );
    List<label>  localNearest( Pstream::nProcs(), 0 );

    distCur[Pstream::myProcNo()] = GREAT;

    forAll( *pf_, facei )
    {
        scalar distance = mag( fc[facei] - origin_ );
        if( distance < distCur[Pstream::myProcNo()] )
        {
            distCur[Pstream::myProcNo()] = distance;
            localNearest[Pstream::myProcNo()] = facei;
        }
    } 

    reduce( localNearest, sumOp<List<label> >() );
    reduce( distCur, sumOp<List<scalar> >() );

    donorProcN_ = findMin( distCur );

    nearest_ = localNearest[Pstream::myProcNo()]; 

    //get bulk velocity of patch as a vector
    scalar patchTotArea = 0.0;
    vector UB = vector(0.0,0.0,0.0); 

    forAll(*pf_, facei)
    {
        UB += pf_->meanField()[facei] * pf_->patch().magSf()[facei];
        patchTotArea += pf_->patch().magSf()[facei];
    }

    reduce( UB, sumOp<vector>() );
    reduce( patchTotArea, sumOp<scalar>() );

    UB /= max( patchTotArea, SMALL);
   
    u_ = UB;

    
    if( donorProcN_ != Pstream::myProcNo() )
    {
         sigma_ = pTraits<vector>::zero;
    }
    else
    {
        sigma_ = pf_->sigma()[nearest_];
    }
 
    reduce( sigma_, maxOp<vector>() );

}

void SEMspot::projectBack( vector nn, bool setToFace )
{
    //project Back from patch
    if( donorProcN_ == Pstream::myProcNo() )
    {
        //set onto patch
        scalar distance=mag( pf_->patch().Cf()[nearest_] - origin_ );
        while(1)
        {
            scalar dist2 = mag( pf_->patch().Cf()[nearest_] - (origin_-0.1*mag(sigma_&nn)*nn) );
            if( dist2 < distance )
            {
                distance = dist2;
                origin_ -= 0.1*mag(sigma_&nn)*nn;
            }
            else
            {
                break;
            }
        }

        if( setToFace )
        {
            origin_ += mag(sigma_&nn)*nn;
        }
        else
        {
            scalar randNum;
            pf_->ranGen().randomise01( randNum );
            randNum*=2.0;
            randNum-=1.0;

            origin_ += randNum*(sigma_&nn)*nn;
        }
    }
    else
    {
        origin_ = pTraits<vector>::one * GREAT;
    }

    reduce( origin_, minOp<vector>() );

}


}
