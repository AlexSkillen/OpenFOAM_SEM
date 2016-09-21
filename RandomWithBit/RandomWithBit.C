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

#include "RandomWithBit.H"


namespace Foam
{


RandomWithBit::RandomWithBit(label seed)
:
    Random(seed)
{
}




void RandomWithBit::bitRandomise( scalar& s )
{
    s=bit();
}

void RandomWithBit::bitRandomise( vector& v )
{
    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        v.component(cmpt) = bit();
    }
}

void RandomWithBit::bitRandomise( sphericalTensor& st )
{
    st.ii() = bit();
}


void RandomWithBit::bitRandomise( tensor& t )
{
    for (direction cmpt=0; cmpt<tensor::nComponents; cmpt++)
    {
        t.component(cmpt) = bit();
    }
}

void RandomWithBit::bitRandomise( symmTensor& st )
{
    for (direction cmpt=0; cmpt<symmTensor::nComponents; cmpt++)
    {
        st.component(cmpt) = bit();
    }
}


}
