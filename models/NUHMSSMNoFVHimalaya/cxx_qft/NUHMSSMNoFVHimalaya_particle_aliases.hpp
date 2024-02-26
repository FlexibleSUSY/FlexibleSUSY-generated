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
 * @file NUHMSSMNoFVHimalaya_particle_aliases.hpp
 *
 * This file was generated with FlexibleSUSY @FlexibleSUSYVersion@ and SARAH @SARAHVersion@ .
 */

#include "cxx_qft/NUHMSSMNoFVHimalaya_fields.hpp"

namespace flexiblesusy {

using Higgs = NUHMSSMNoFVHimalaya_cxx_diagrams::fields::hh;
using PseudoscalarHiggs = NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ah;
using Hm = NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Hpm;
using Hp = typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Hpm>::type;
using WmBoson = NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VWm;
using WpBoson = typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VWm>::type;
using Photon = NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP;
using ZBoson = NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VZ;
using Gluon = NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VG;
using Electron = NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fe;
using Muon = NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm;
using Tauon = NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ftau;
using ElectronNeutrino = NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fve;
using MuonNeutrino = NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fvm;
using TauNeutrino = NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fvt;

}
