// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================


/**
 * @file HGTHDMIIMSSMBC_particle_aliases.hpp
 *
 * This file was generated with FlexibleSUSY @FlexibleSUSYVersion@ and SARAH @SARAHVersion@ .
 */

#include "cxx_qft/HGTHDMIIMSSMBC_fields.hpp"

namespace flexiblesusy {

using Higgs = HGTHDMIIMSSMBC_cxx_diagrams::fields::hh;
using PseudoscalarHiggs = HGTHDMIIMSSMBC_cxx_diagrams::fields::Ah;
using Hm = HGTHDMIIMSSMBC_cxx_diagrams::fields::Hm;
using Hp = typename HGTHDMIIMSSMBC_cxx_diagrams::fields::conj<HGTHDMIIMSSMBC_cxx_diagrams::fields::Hm>::type;
using WmBoson = HGTHDMIIMSSMBC_cxx_diagrams::fields::VWm;
using WpBoson = typename HGTHDMIIMSSMBC_cxx_diagrams::fields::conj<HGTHDMIIMSSMBC_cxx_diagrams::fields::VWm>::type;
using Photon = HGTHDMIIMSSMBC_cxx_diagrams::fields::VP;
using ZBoson = HGTHDMIIMSSMBC_cxx_diagrams::fields::VZ;
using Gluon = HGTHDMIIMSSMBC_cxx_diagrams::fields::VG;
using ChargedLepton = HGTHDMIIMSSMBC_cxx_diagrams::fields::Fe;
using Neutrino = HGTHDMIIMSSMBC_cxx_diagrams::fields::Fv;
using UpTypeQuark = HGTHDMIIMSSMBC_cxx_diagrams::fields::Fu;
using DownTypeQuark = HGTHDMIIMSSMBC_cxx_diagrams::fields::Fd;

}
