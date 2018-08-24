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

#include <Grid/qcd/action/fermion/CayleyFermion5DdenseMethodImpl.h>


namespace Grid {
namespace QCD {

#ifdef CAYLEY_DPERP_DENSE
INSTANTIATE_DPERP(GparityWilsonImplF);
INSTANTIATE_DPERP(GparityWilsonImplD);
INSTANTIATE_DPERP(WilsonImplF);
INSTANTIATE_DPERP(WilsonImplD);
INSTANTIATE_DPERP(ZWilsonImplF);
INSTANTIATE_DPERP(ZWilsonImplD);

template void CayleyFermion5D<GparityWilsonImplF>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);
template void CayleyFermion5D<GparityWilsonImplD>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);
template void CayleyFermion5D<WilsonImplF>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);
template void CayleyFermion5D<WilsonImplD>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);
template void CayleyFermion5D<ZWilsonImplF>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);
template void CayleyFermion5D<ZWilsonImplD>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);

INSTANTIATE_DPERP(GparityWilsonImplFH);
INSTANTIATE_DPERP(GparityWilsonImplDF);
INSTANTIATE_DPERP(WilsonImplFH);
INSTANTIATE_DPERP(WilsonImplDF);
INSTANTIATE_DPERP(ZWilsonImplFH);
INSTANTIATE_DPERP(ZWilsonImplDF);

template void CayleyFermion5D<GparityWilsonImplFH>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);
template void CayleyFermion5D<GparityWilsonImplDF>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);
template void CayleyFermion5D<WilsonImplFH>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);
template void CayleyFermion5D<WilsonImplDF>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);
template void CayleyFermion5D<ZWilsonImplFH>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);
template void CayleyFermion5D<ZWilsonImplDF>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);
#endif

}}
