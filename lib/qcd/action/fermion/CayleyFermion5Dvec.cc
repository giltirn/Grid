    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/CayleyFermion5D.cc

    Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
Author: paboyle <paboyle@ph.ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */


#include <Grid/qcd/action/fermion/CayleyFermion5DvecMethodImpl.h>


namespace Grid {
namespace QCD {  

INSTANTIATE_DPERP(DomainWallVec5dImplD);
INSTANTIATE_DPERP(DomainWallVec5dImplF);
INSTANTIATE_DPERP(ZDomainWallVec5dImplD);
INSTANTIATE_DPERP(ZDomainWallVec5dImplF);

INSTANTIATE_DPERP(DomainWallVec5dImplDF);
INSTANTIATE_DPERP(DomainWallVec5dImplFH);
INSTANTIATE_DPERP(ZDomainWallVec5dImplDF);
INSTANTIATE_DPERP(ZDomainWallVec5dImplFH);

template void CayleyFermion5D<DomainWallVec5dImplF>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);
template void CayleyFermion5D<DomainWallVec5dImplD>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);
template void CayleyFermion5D<ZDomainWallVec5dImplF>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);
template void CayleyFermion5D<ZDomainWallVec5dImplD>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);

template void CayleyFermion5D<DomainWallVec5dImplFH>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);
template void CayleyFermion5D<DomainWallVec5dImplDF>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);
template void CayleyFermion5D<ZDomainWallVec5dImplFH>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);
template void CayleyFermion5D<ZDomainWallVec5dImplDF>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);

}}
