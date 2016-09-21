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

#include "turbulentInletSEMFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "SEMspot.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


turbulentInletSEMFvPatchField::turbulentInletSEMFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    SEMBase(p, iF),
    a_(p.size()),
    fMean2_(p.size())
{
}



turbulentInletSEMFvPatchField::turbulentInletSEMFvPatchField
(
    const turbulentInletSEMFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    SEMBase(ptf, p, iF, mapper),
    a_(ptf.a_, mapper),
    fMean2_(ptf.fMean2_,mapper)
{
}



turbulentInletSEMFvPatchField::turbulentInletSEMFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    SEMBase(p, iF, dict),
    a_(p.size()),
    fMean2_(p.size())
{
    if (dict.found("fMean2"))
    {
        fMean2_ = Field<scalar>("fMean2", dict, p.size());
    }
    else
    {
        fMean2_ = Field<scalar>(p.size(), 1.0);
    }

    this->initilise();
}



turbulentInletSEMFvPatchField::turbulentInletSEMFvPatchField
(
    const turbulentInletSEMFvPatchField& ptf
)
:
    SEMBase(ptf),
    a_(ptf.a_),
    fMean2_(ptf.fMean2_)
{
    for( int i=0; i<ptf.spot_.size(); i++ )
    {
        spot_.append( new SEMspot( ptf.spot_[i], this ) );
    }
}



turbulentInletSEMFvPatchField::turbulentInletSEMFvPatchField
(
    const turbulentInletSEMFvPatchField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    SEMBase(ptf, iF),
    a_(ptf.a_),
    fMean2_(ptf.fMean2_)
{
    for( int i=0; i<ptf.spot_.size(); i++ )
    {
        spot_.append( new SEMspot( ptf.spot_[i], this ) );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void turbulentInletSEMFvPatchField::allocateSpots()
{
    int n = numEddies();

    //generate spots
    for( int i=0; i<n; i++ ) 
    {
        spot_.append( new SEMspot(this) );
        spot_.last()->initialise(false);
    }
}


void turbulentInletSEMFvPatchField::initilise()
{
    SEMBase::initilise();

    //allocate spots
    this->allocateSpots();

    //set a_
    a_ = pTraits<tensor>::zero;

    forAll(*this, facei)
    {
        a_[facei].xx() = sqrt( max( RIn_[facei].xx(), SMALL ) );
        a_[facei].yx() = RIn_[facei].xy()/max( a_[facei].xx(), SMALL );
        a_[facei].yy() = sqrt( max( RIn_[facei].yy() - pow( a_[facei].yx(), 2), SMALL ));
        a_[facei].zx() = RIn_[facei].xz()/max( a_[facei].xx(), SMALL );
        a_[facei].zy() = (RIn_[facei].yz() - a_[facei].yx()*a_[facei].zx())/max( a_[facei].yy(), SMALL );
        a_[facei].zz() = 
            sqrt( max( RIn_[facei].zz() - pow( a_[facei].zx(), 2) - pow( a_[facei].zy(), 2), SMALL));
    }
}

//---
scalar turbulentInletSEMFvPatchField::f( vector& distance, vector& sigma )
{
    return 
    ( 
        ( 1.0 - mag( distance.x() ) / sigma.x() ) *
        ( 1.0 - mag( distance.y() ) / sigma.y() ) *
        ( 1.0 - mag( distance.z() ) / sigma.z() )
    ) * pow(1.5,1.5);
}





//---
void turbulentInletSEMFvPatchField::updateU()
{
    Field<vector>& patchField = *this;
    const vectorField& fc = this->patch().Cf();

    Field<vector> fluctuatingField( this->size() );
    
    patchField = vector( 0.0, 0.0, 0.0 );
    fluctuatingField = vector( 0.0, 0.0, 0.0 );


    forAll(fluctuatingField, facei)
    {
        scalar sum_f2=0.0;

        for( int i=0; i<spot_.size(); i++ )
        {

            vector distance = fc[facei] - spot_[i]->origin();

            if
            ( 
                mag( distance.x() ) / spot_[i]->sigma().x() < 1.0
             && mag( distance.y() ) / spot_[i]->sigma().y() < 1.0
             && mag( distance.z() ) / spot_[i]->sigma().z() < 1.0
            )
            {

                scalar shapeFun = f( distance, spot_[i]->sigma() );

                sum_f2 += pow( shapeFun, 2 );

                fluctuatingField[facei] += spot_[i]->epsilon()*shapeFun;
            }
        }
       
        scalar dt = this->db().time().deltaTValue();

        scalar alpha = min( max( dt/avgWindow_, 0.0 ), 1.0 );

        fMean2_[facei] = alpha * sum_f2 + (1.0-alpha)*fMean2_[facei]; 

        fluctuatingField[facei] /= pow( max(fMean2_[facei], SMALL), 0.5 ) ;
        

        patchField[facei] += a_[facei] & fluctuatingField[facei];
    }

    patchField += meanField_;

}


void turbulentInletSEMFvPatchField::write(Ostream& os) const
{
    SEMBase::write(os);
    fMean2_.writeEntry("fMean2", os);
}
    


turbulentInletSEMFvPatchField::~turbulentInletSEMFvPatchField()
{
    for( int i=0; i<spot_.size(); i++ )
    {
        delete spot_[i];
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{   
    makePatchTypeField
    (   
        fvPatchVectorField,
        turbulentInletSEMFvPatchField
    );
}

// ************************************************************************* //
