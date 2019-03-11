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

#include "SEMBase.H"
#include "addToRunTimeSelectionTable.H"
#include "polyBoundaryMeshEntries.H"
#include "SEMspot.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTemplateTypeNameAndDebugWithName(IOList<symmTensor>, "symmTensorList", 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


SEMBase::SEMBase
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    ranGen_(label(0)),
    meanField_(p.size()),
    curTimeIndex_(-1),
    RIn_(p.size()),
    sigma_(p.size())
{
}



SEMBase::SEMBase
(
    const SEMBase& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    ranGen_(label(0)),
    meanField_(ptf.meanField_, mapper),
    curTimeIndex_(-1),
    semBox_(ptf.semBox_),
    RIn_(ptf.RIn_, mapper),
    sigma_(ptf.sigma_, mapper),
    avgWindow_(ptf.avgWindow_)
{
}



SEMBase::SEMBase
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF),
    ranGen_(label(0)),
    meanField_(p.size()),
    curTimeIndex_(-1),
    semBox_(p.Cf()), 
    RIn_(p.size())
{
    fixedValueFvPatchField<vector>::operator==(vector(0.0, 0.0, 0.0));

    meanField_ = Field<vector>("UIn", dict, p.size());
    RIn_ = Field<symmTensor>("RIn", dict, p.size()); 
    sigma_ = Field<vector>("sigma", dict, p.size());
}



SEMBase::SEMBase
(
    const SEMBase& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    ranGen_(ptf.ranGen_),
    meanField_(ptf.meanField_),
    curTimeIndex_(-1),
    semBox_(ptf.semBox_), 
    RIn_(ptf.RIn_),
    sigma_(ptf.sigma_), 
    spot_(ptf.spot_), 
    avgWindow_(ptf.avgWindow_)
{
}



SEMBase::SEMBase
(
    const SEMBase& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    ranGen_(ptf.ranGen_),
    meanField_(ptf.meanField_),
    curTimeIndex_(ptf.curTimeIndex_),
    semBox_(ptf.semBox_),
    UBulk_(ptf.UBulk_),
    RIn_(ptf.RIn_),
    sigma_(ptf.sigma_), 
    avgWindow_(ptf.avgWindow_)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void SEMBase::initilise()
{
    //compute bulk velocity of patch
    scalar patchTotArea( 0.0 );
    UBulk_ = 0.0;

    forAll(*this, facei)
    {
        UBulk_ += ( meanField_[facei] & this->patch().nf()()[facei] ) * this->patch().magSf()[facei];
        patchTotArea += this->patch().magSf()[facei];
    }

    reduce( UBulk_, sumOp<scalar>() );
    reduce( patchTotArea, sumOp<scalar>() );

    UBulk_ /= max( patchTotArea, SMALL);


    forAll(*this, facei)
    {
        vector inflate( sigma_[facei] - ( this->patch().nf()()[facei] & sigma_[facei] ) * this->patch().nf()()[facei] );

        for( int k=0; k<3; k++ )
        {
            semBox_.min()[k] = min( semBox_.min()[k], this->patch().Cf()[facei][k] - inflate[k] );
            semBox_.max()[k] = max( semBox_.max()[k], this->patch().Cf()[facei][k] + inflate[k] );
        }
    }
    
    reduce( semBox_.min(), minOp<vector>() );
    reduce( semBox_.max(), maxOp<vector>() );


    //set averaging window size
    avgWindow_ = cmptMax( max(sigma_) )/mag(UBulk_) * 5.0;
    reduce( avgWindow_, maxOp<scalar>() );
}

int SEMBase::numEddies()
{
    int numSpots = 0;
    const int maxSpots=10000;

    scalar minLength = GREAT;
    
    forAll(*this, facei)
    {
        minLength = min( minLength, cmptMin(sigma_[facei]));
    }

    reduce( minLength, minOp<scalar>() );

    scalar patchTotArea( gSum(this->patch().magSf()) );
    numSpots = 5.0 * patchTotArea / (2.0/3.0*3.14 * pow(minLength,2) + SMALL);

    reduce( numSpots, maxOp<scalar>() );
    label requested = numSpots;
    numSpots = min( numSpots, maxSpots );
    Info<<patchTotArea<<" "<<minLength<<endl;
    if( numSpots == maxSpots )
    {
        WarningIn("SEMBase::numEddies()")
            << "Maximum number of eddies reached\n"
            << "\tRequested: "
            << requested
            << endl
            << "\tLimited to: "
            << maxSpots 
            << endl;
    }
    else
    {
        Info<< "Number of eddies used: " 
            << numSpots
            << endl;
    }

    return numSpots; 
}



//---
void SEMBase::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<vector>::autoMap(m);
    meanField_.autoMap(m);
    RIn_.autoMap(m);
    sigma_.autoMap(m);
}


//---
void SEMBase::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<vector>::rmap(ptf, addr);

    const SEMBase& tiptf =
        refCast<const SEMBase >(ptf);

    meanField_.rmap(tiptf.meanField_, addr);
    RIn_.rmap(tiptf.RIn_, addr);
    sigma_.rmap(tiptf.sigma_, addr);
}



//---
void SEMBase::advectPoints()
{
    //advect random points
    scalar dt = this->db().time().deltaTValue();

    for( int i=0; i<spot_.size(); i++ ) 
    {
        spot_[i]->origin() += spot_[i]->u()*dt;

        bool regen =  this->db().time().value() > spot_[i]->residenceTime();
        reduce( regen, orOp<bool>() );

        if(regen)
        {
            scalar origSize = mag( spot_[i]->sigma() );
            do
            {
                spot_[i]->initialise(true);
            } while( mag( spot_[i]->sigma() ) > 1.1*origSize || mag( spot_[i]->sigma() ) < 0.9*origSize );
           
        }
    }
}





void SEMBase::correctMass() 
{
    //correct bulk velocity of patch
    scalar patchTotArea = 0.0;
    scalar Uc = 0.0;

    forAll(*this, facei)
    {
        Uc += ( (*this)[facei] & this->patch().nf()()[facei] ) * this->patch().magSf()[facei];
        patchTotArea += this->patch().magSf()[facei];
    }

    reduce( Uc, sumOp<scalar>() );
    reduce( patchTotArea, sumOp<scalar>() );


    Uc /= max( patchTotArea, SMALL);

    forAll(*this, facei)
    { 
        (*this)[facei] *= UBulk_/Uc;
    } 
}



//---
void SEMBase::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    
    if (curTimeIndex_ != this->db().time().timeIndex())
    {

        this->advectPoints();            
   
        this->updateU(); 
        
        this->correctMass(); 
        
        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchField<vector>::updateCoeffs();

}






void SEMBase::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    this->writeEntry("value", os);
    meanField_.writeEntry("UIn", os);
    RIn_.writeEntry("RIn", os);
    sigma_.writeEntry("sigma", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
