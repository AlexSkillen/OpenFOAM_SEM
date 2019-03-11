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

#ifndef SEMBase_H
#define SEMBase_H

#include "Random.H"
#include "fixedValueFvPatchFields.H"
#include "primitivePatchInterpolation.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class SEMspot;

/*---------------------------------------------------------------------------*\
                     Class turbulentInletSEMFvPatch Declaration
\*---------------------------------------------------------------------------*/

class SEMBase
:
    public fixedValueFvPatchField<vector>
{
    // Protected data

protected:
        Random ranGen_;
        Field<vector> meanField_;
        label curTimeIndex_; 
        boundBox semBox_;
        scalar UBulk_;
        Field<symmTensor> RIn_;        
        Field<vector> sigma_; 
        DynamicList<SEMspot*> spot_;
        scalar avgWindow_;

public:
    // Constructors

        //- Construct from patch and internal field
        SEMBase
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        SEMBase
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given SEMBase
        //  onto a new patch
        SEMBase
        (
            const SEMBase&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        SEMBase
        (
            const SEMBase&
        );

        //- Construct as copy setting internal field reference
        SEMBase
        (
            const SEMBase&,
            const DimensionedField<vector, volMesh>&
        );


    // Member functions

        // Access


            const Field<vector>& meanField() const
            {
                return meanField_;
            }

            Field<vector>& meanField()
            {
                return meanField_;
            }

            Field<vector>& sigma()
            {
                return sigma_;
            }

            Field<symmTensor>& RIn()
            {
                return RIn_;
            }

            Random& ranGen()
            {
                return ranGen_;
            }

            boundBox semBox()
            {
                return semBox_;
            }

            List<SEMspot*>& spot()
            {
                return spot_;
            }

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<vector>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();
         
            //- Write
            virtual void write(Ostream&) const;




protected:
            virtual void updateU()=0;
            virtual void advectPoints(); 
            virtual void correctMass();
            virtual void initilise();
            virtual int numEddies();
            virtual void allocateSpots()=0;

      
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
