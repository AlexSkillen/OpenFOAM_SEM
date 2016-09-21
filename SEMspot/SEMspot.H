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

#ifndef SEMspot_H
#define SEMspot_H

#include "fixedValueFvPatchFields.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class SEMBase; 



class SEMspot
{

    SEMBase* pf_;

protected:
    vector origin_;
    vector sigma_;
    vector epsilon_;
    vector u_;
    label nearest_;
    label donorProcN_;
    scalar residenceTime_;

    virtual void interpolatePatchToSpot();
    virtual void projectBack( vector nn, bool );


public:
    
    //constructor
    SEMspot(SEMBase*);  
    SEMspot(SEMspot*, SEMBase*);
    SEMspot() {} 

 
    //- Return the spot origin
    const vector& origin() const
    {
        return origin_;
    }

    //- Return reference to the spot origin to allow adjustment
    vector& origin()
    {
        return origin_;
    }

    //- Return the spot dimensions
    const vector& sigma() const
    {
        return sigma_;
    }

    //- Return reference to the spot dimensions to allow adjustment
    vector& sigma()
    {
        return sigma_;
    }
    
    
    vector& epsilon()
    {
        return epsilon_;
    }

    vector& u()
    {
        return u_;
    }

    scalar residenceTime()
    {
        return residenceTime_;
    }

    label& nearest()
    {  
        return nearest_;
    }

    label& donorProcN()
    {
        return donorProcN_;
    }
 
    void initialise( const bool setToFace=true);

    virtual ~SEMspot() {}
};

}

#endif
