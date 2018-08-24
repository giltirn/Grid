/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/DomainWallEOFAFermionvec.cc

Copyright (C) 2017

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: David Murphy <dmurphy@phys.columbia.edu>

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

#include <Grid/qcd/action/fermion/DomainWallEOFAFermionvecMethodImpl.h>

namespace Grid {
namespace QCD {

    #ifdef DOMAIN_WALL_EOFA_DPERP_VEC

        INSTANTIATE_DPERP_DWF_EOFA(DomainWallVec5dImplD);
        INSTANTIATE_DPERP_DWF_EOFA(DomainWallVec5dImplF);
        INSTANTIATE_DPERP_DWF_EOFA(ZDomainWallVec5dImplD);
        INSTANTIATE_DPERP_DWF_EOFA(ZDomainWallVec5dImplF);

        INSTANTIATE_DPERP_DWF_EOFA(DomainWallVec5dImplDF);
        INSTANTIATE_DPERP_DWF_EOFA(DomainWallVec5dImplFH);
        INSTANTIATE_DPERP_DWF_EOFA(ZDomainWallVec5dImplDF);
        INSTANTIATE_DPERP_DWF_EOFA(ZDomainWallVec5dImplFH);

        template void DomainWallEOFAFermion<DomainWallVec5dImplF>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);
        template void DomainWallEOFAFermion<DomainWallVec5dImplD>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);
        template void DomainWallEOFAFermion<ZDomainWallVec5dImplF>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);
        template void DomainWallEOFAFermion<ZDomainWallVec5dImplD>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);

        template void DomainWallEOFAFermion<DomainWallVec5dImplFH>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);
        template void DomainWallEOFAFermion<DomainWallVec5dImplDF>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);
        template void DomainWallEOFAFermion<ZDomainWallVec5dImplFH>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);
        template void DomainWallEOFAFermion<ZDomainWallVec5dImplDF>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);

    #endif

}}
