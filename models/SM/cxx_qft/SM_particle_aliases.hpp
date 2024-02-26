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
 * @file SM_particle_aliases.hpp
 *
 * This file was generated with FlexibleSUSY @FlexibleSUSYVersion@ and SARAH @SARAHVersion@ .
 */

#include "cxx_qft/SM_fields.hpp"

namespace flexiblesusy {

using Higgs = SM_cxx_diagrams::fields::hh;
using WpBoson = SM_cxx_diagrams::fields::VWp;
using WmBoson = typename SM_cxx_diagrams::fields::conj<SM_cxx_diagrams::fields::VWp>::type;
using Photon = SM_cxx_diagrams::fields::VP;
using ZBoson = SM_cxx_diagrams::fields::VZ;
using Gluon = SM_cxx_diagrams::fields::VG;
using ChargedLepton = SM_cxx_diagrams::fields::Fe;
using Neutrino = SM_cxx_diagrams::fields::Fv;
using UpTypeQuark = SM_cxx_diagrams::fields::Fu;
using DownTypeQuark = SM_cxx_diagrams::fields::Fd;

}
