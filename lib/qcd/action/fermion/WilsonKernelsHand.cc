    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/WilsonKernelsHand.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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
#include <Grid/qcd/action/fermion/WilsonKernelsHandMethodImpl.h>

namespace Grid {
namespace QCD {

INSTANTIATE_THEM(WilsonImplF);
INSTANTIATE_THEM(WilsonImplD);
INSTANTIATE_THEM(ZWilsonImplF);
INSTANTIATE_THEM(ZWilsonImplD);
INSTANTIATE_THEM(DomainWallVec5dImplF);
INSTANTIATE_THEM(DomainWallVec5dImplD);
INSTANTIATE_THEM(ZDomainWallVec5dImplF);
INSTANTIATE_THEM(ZDomainWallVec5dImplD);
INSTANTIATE_THEM(WilsonImplFH);
INSTANTIATE_THEM(WilsonImplDF);
INSTANTIATE_THEM(ZWilsonImplFH);
INSTANTIATE_THEM(ZWilsonImplDF);
INSTANTIATE_THEM(DomainWallVec5dImplFH);
INSTANTIATE_THEM(DomainWallVec5dImplDF);
INSTANTIATE_THEM(ZDomainWallVec5dImplFH);
INSTANTIATE_THEM(ZDomainWallVec5dImplDF);
INSTANTIATE_THEM(WilsonTwoIndexAntiSymmetricImplF);
INSTANTIATE_THEM(WilsonTwoIndexAntiSymmetricImplD);

}}
