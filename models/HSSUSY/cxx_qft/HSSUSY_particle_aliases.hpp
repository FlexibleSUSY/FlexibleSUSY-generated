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
 * @file HSSUSY_particle_aliases.hpp
 *
 * This file was generated with FlexibleSUSY @FlexibleSUSYVersion@ and SARAH @SARAHVersion@ .
 */

#include "cxx_qft/HSSUSY_fields.hpp"

namespace flexiblesusy {

using Higgs = HSSUSY_cxx_diagrams::fields::hh;
using WpBoson = HSSUSY_cxx_diagrams::fields::VWp;
using WmBoson = typename HSSUSY_cxx_diagrams::fields::conj<HSSUSY_cxx_diagrams::fields::VWp>::type;
using Photon = HSSUSY_cxx_diagrams::fields::VP;
using ZBoson = HSSUSY_cxx_diagrams::fields::VZ;
using Gluon = HSSUSY_cxx_diagrams::fields::VG;
using ChargedLepton = HSSUSY_cxx_diagrams::fields::Fe;
using Neutrino = HSSUSY_cxx_diagrams::fields::Fv;
using UpTypeQuark = HSSUSY_cxx_diagrams::fields::Fu;
using DownTypeQuark = HSSUSY_cxx_diagrams::fields::Fd;

}
