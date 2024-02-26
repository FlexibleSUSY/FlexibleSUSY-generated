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
 * @file MSSMNoFVHimalaya_muon_amm_wrapper.cpp
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#include "MSSMNoFVHimalaya_lepton_amm_wrapper.hpp"
#include "MSSMNoFVHimalaya_amm.hpp"
#include "cxx_qft/MSSMNoFVHimalaya_qft.hpp"

namespace flexiblesusy {
namespace MSSMNoFVHimalaya_lepton_amm_wrapper {

using namespace MSSMNoFVHimalaya_cxx_diagrams;

double calculate_Fm_amm(const MSSMNoFVHimalaya_mass_eigenstates& model, const softsusy::QedQcd& qedqcd, const Spectrum_generator_settings& settings) {
   return MSSMNoFVHimalaya_amm::calculate_amm<fields::Fm>(model, qedqcd, settings);
}

double calculate_Fm_amm_uncertainty(const MSSMNoFVHimalaya_mass_eigenstates& model, const softsusy::QedQcd& qedqcd, const Spectrum_generator_settings& settings) {
   return MSSMNoFVHimalaya_amm::calculate_amm_uncertainty<fields::Fm>(model, qedqcd, settings);
}

}
}
