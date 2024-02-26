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
 * @file NUHMSSM_particle_aliases.hpp
 *
 * This file was generated with FlexibleSUSY @FlexibleSUSYVersion@ and SARAH @SARAHVersion@ .
 */

#include "cxx_qft/NUHMSSM_fields.hpp"

namespace flexiblesusy {

using Higgs = NUHMSSM_cxx_diagrams::fields::hh;
using PseudoscalarHiggs = NUHMSSM_cxx_diagrams::fields::Ah;
using Hm = NUHMSSM_cxx_diagrams::fields::Hpm;
using Hp = typename NUHMSSM_cxx_diagrams::fields::conj<NUHMSSM_cxx_diagrams::fields::Hpm>::type;
using WmBoson = NUHMSSM_cxx_diagrams::fields::VWm;
using WpBoson = typename NUHMSSM_cxx_diagrams::fields::conj<NUHMSSM_cxx_diagrams::fields::VWm>::type;
using Photon = NUHMSSM_cxx_diagrams::fields::VP;
using ZBoson = NUHMSSM_cxx_diagrams::fields::VZ;
using Gluon = NUHMSSM_cxx_diagrams::fields::VG;
using ChargedLepton = NUHMSSM_cxx_diagrams::fields::Fe;
using Neutrino = NUHMSSM_cxx_diagrams::fields::Fv;
using UpTypeQuark = NUHMSSM_cxx_diagrams::fields::Fu;
using DownTypeQuark = NUHMSSM_cxx_diagrams::fields::Fd;

}
