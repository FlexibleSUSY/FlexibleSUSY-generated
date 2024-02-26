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
 * @file THDMII_particle_aliases.hpp
 *
 * This file was generated with FlexibleSUSY @FlexibleSUSYVersion@ and SARAH @SARAHVersion@ .
 */

#include "cxx_qft/THDMII_fields.hpp"

namespace flexiblesusy {

using Higgs = THDMII_cxx_diagrams::fields::hh;
using PseudoscalarHiggs = THDMII_cxx_diagrams::fields::Ah;
using Hm = THDMII_cxx_diagrams::fields::Hm;
using Hp = typename THDMII_cxx_diagrams::fields::conj<THDMII_cxx_diagrams::fields::Hm>::type;
using WmBoson = THDMII_cxx_diagrams::fields::VWm;
using WpBoson = typename THDMII_cxx_diagrams::fields::conj<THDMII_cxx_diagrams::fields::VWm>::type;
using Photon = THDMII_cxx_diagrams::fields::VP;
using ZBoson = THDMII_cxx_diagrams::fields::VZ;
using Gluon = THDMII_cxx_diagrams::fields::VG;
using ChargedLepton = THDMII_cxx_diagrams::fields::Fe;
using Neutrino = THDMII_cxx_diagrams::fields::Fv;
using UpTypeQuark = THDMII_cxx_diagrams::fields::Fu;
using DownTypeQuark = THDMII_cxx_diagrams::fields::Fd;

}
