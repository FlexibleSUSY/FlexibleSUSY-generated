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
 * @file MSSM_particle_aliases.hpp
 *
 * This file was generated with FlexibleSUSY @FlexibleSUSYVersion@ and SARAH @SARAHVersion@ .
 */

#include "cxx_qft/MSSM_fields.hpp"

namespace flexiblesusy {

using Higgs = MSSM_cxx_diagrams::fields::hh;
using PseudoscalarHiggs = MSSM_cxx_diagrams::fields::Ah;
using Hm = MSSM_cxx_diagrams::fields::Hpm;
using Hp = typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Hpm>::type;
using WmBoson = MSSM_cxx_diagrams::fields::VWm;
using WpBoson = typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::VWm>::type;
using Photon = MSSM_cxx_diagrams::fields::VP;
using ZBoson = MSSM_cxx_diagrams::fields::VZ;
using Gluon = MSSM_cxx_diagrams::fields::VG;
using ChargedLepton = MSSM_cxx_diagrams::fields::Fe;
using Neutrino = MSSM_cxx_diagrams::fields::Fv;
using UpTypeQuark = MSSM_cxx_diagrams::fields::Fu;
using DownTypeQuark = MSSM_cxx_diagrams::fields::Fd;

}
