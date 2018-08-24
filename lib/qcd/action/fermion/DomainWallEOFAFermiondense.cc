/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/DomainWallEOFAFermiondense.cc

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

#include <Grid/qcd/action/fermion/DomainWallEOFAFermiondenseMethodImpl.h>

namespace Grid {
namespace QCD {

    #ifdef DOMAIN_WALL_EOFA_DPERP_DENSE

        INSTANTIATE_DPERP_DWF_EOFA(GparityWilsonImplF);
        INSTANTIATE_DPERP_DWF_EOFA(GparityWilsonImplD);
        INSTANTIATE_DPERP_DWF_EOFA(WilsonImplF);
        INSTANTIATE_DPERP_DWF_EOFA(WilsonImplD);
        INSTANTIATE_DPERP_DWF_EOFA(ZWilsonImplF);
        INSTANTIATE_DPERP_DWF_EOFA(ZWilsonImplD);

        template void DomainWallEOFAFermion<GparityWilsonImplF>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);
        template void DomainWallEOFAFermion<GparityWilsonImplD>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);
        template void DomainWallEOFAFermion<WilsonImplF>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);
        template void DomainWallEOFAFermion<WilsonImplD>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);
        template void DomainWallEOFAFermion<ZWilsonImplF>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);
        template void DomainWallEOFAFermion<ZWilsonImplD>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);

        INSTANTIATE_DPERP_DWF_EOFA(GparityWilsonImplFH);
        INSTANTIATE_DPERP_DWF_EOFA(GparityWilsonImplDF);
        INSTANTIATE_DPERP_DWF_EOFA(WilsonImplFH);
        INSTANTIATE_DPERP_DWF_EOFA(WilsonImplDF);
        INSTANTIATE_DPERP_DWF_EOFA(ZWilsonImplFH);
        INSTANTIATE_DPERP_DWF_EOFA(ZWilsonImplDF);

        template void DomainWallEOFAFermion<GparityWilsonImplFH>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);
        template void DomainWallEOFAFermion<GparityWilsonImplDF>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);
        template void DomainWallEOFAFermion<WilsonImplFH>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);
        template void DomainWallEOFAFermion<WilsonImplDF>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);
        template void DomainWallEOFAFermion<ZWilsonImplFH>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);
        template void DomainWallEOFAFermion<ZWilsonImplDF>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);

    #endif

}}
