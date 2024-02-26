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
 * @file UMSSM_particle_aliases.hpp
 *
 * This file was generated with FlexibleSUSY @FlexibleSUSYVersion@ and SARAH @SARAHVersion@ .
 */

#include "cxx_qft/UMSSM_fields.hpp"

namespace flexiblesusy {

using Higgs = UMSSM_cxx_diagrams::fields::hh;
using PseudoscalarHiggs = UMSSM_cxx_diagrams::fields::Ah;
using Hm = UMSSM_cxx_diagrams::fields::Hpm;
using Hp = typename UMSSM_cxx_diagrams::fields::conj<UMSSM_cxx_diagrams::fields::Hpm>::type;
using WmBoson = UMSSM_cxx_diagrams::fields::VWm;
using WpBoson = typename UMSSM_cxx_diagrams::fields::conj<UMSSM_cxx_diagrams::fields::VWm>::type;
using Photon = UMSSM_cxx_diagrams::fields::VP;
using ZBoson = UMSSM_cxx_diagrams::fields::VZ;
using Gluon = UMSSM_cxx_diagrams::fields::VG;
using ChargedLepton = UMSSM_cxx_diagrams::fields::Fe;
using Neutrino = UMSSM_cxx_diagrams::fields::Fv;
using UpTypeQuark = UMSSM_cxx_diagrams::fields::Fu;
using DownTypeQuark = UMSSM_cxx_diagrams::fields::Fd;

}
