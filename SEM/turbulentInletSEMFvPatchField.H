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
#ifndef turbulentInletSEMFvPatchField_H
#define turbulentInletSEMFvPatchField_H

#include "Random.H"
#include "SEMBase.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class turbulentInletSEMFvPatch Declaration
\*---------------------------------------------------------------------------*/

class turbulentInletSEMFvPatchField
:
    public SEMBase
{
    // Protected data

protected:
        Field<tensor> a_;
        Field<scalar> fMean2_;



public:

    //- Runtime type information
    TypeName("turbulentInletSEM");

    // Constructors

        //- Construct from patch and internal field
        turbulentInletSEMFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentInletSEMFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given turbulentInletSEMFvPatchField
        //  onto a new patch
        turbulentInletSEMFvPatchField
        (
            const turbulentInletSEMFvPatchField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentInletSEMFvPatchField
        (
            const turbulentInletSEMFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new turbulentInletSEMFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        turbulentInletSEMFvPatchField
        (
            const turbulentInletSEMFvPatchField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new turbulentInletSEMFvPatchField(*this, iF)
            );
        }

        ~turbulentInletSEMFvPatchField();
    // Member functions
protected:
        virtual void updateU();
        virtual void initilise() ;
        virtual void allocateSpots();  
        virtual void write(Ostream&) const;

private:
        scalar f( vector&, vector& );
       
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
