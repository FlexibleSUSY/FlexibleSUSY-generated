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
 * @file NMSSMEFTHiggs_particle_aliases.hpp
 *
 * This file was generated with FlexibleSUSY @FlexibleSUSYVersion@ and SARAH @SARAHVersion@ .
 */

#include "cxx_qft/NMSSMEFTHiggs_fields.hpp"

namespace flexiblesusy {

using Higgs = NMSSMEFTHiggs_cxx_diagrams::fields::hh;
using PseudoscalarHiggs = NMSSMEFTHiggs_cxx_diagrams::fields::Ah;
using Hm = NMSSMEFTHiggs_cxx_diagrams::fields::Hpm;
using Hp = typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type;
using WmBoson = NMSSMEFTHiggs_cxx_diagrams::fields::VWm;
using WpBoson = typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::VWm>::type;
using Photon = NMSSMEFTHiggs_cxx_diagrams::fields::VP;
using ZBoson = NMSSMEFTHiggs_cxx_diagrams::fields::VZ;
using Gluon = NMSSMEFTHiggs_cxx_diagrams::fields::VG;
using ChargedLepton = NMSSMEFTHiggs_cxx_diagrams::fields::Fe;
using Neutrino = NMSSMEFTHiggs_cxx_diagrams::fields::Fv;
using UpTypeQuark = NMSSMEFTHiggs_cxx_diagrams::fields::Fu;
using DownTypeQuark = NMSSMEFTHiggs_cxx_diagrams::fields::Fd;

}
