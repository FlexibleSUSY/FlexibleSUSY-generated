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

// File generated at Tue 8 Mar 2016 16:18:11

#include "MRSSM_observables.hpp"
#include "MRSSM_mass_eigenstates.hpp"
#include "MRSSM_effective_couplings.hpp"
#include "gm2calc_interface.hpp"
#include "eigen_utils.hpp"
#include "numerics2.hpp"
#include "wrappers.hpp"
#include "lowe.h"
#include "physical_input.hpp"

#define MODEL model
#define AMUGM2CALC a_muon_gm2calc
#define AMUGM2CALCUNCERTAINTY a_muon_gm2calc_uncertainty
#define EFFCPHIGGSPHOTONPHOTON eff_cp_higgs_photon_photon
#define EFFCPHIGGSGLUONGLUON eff_cp_higgs_gluon_gluon
#define EFFCPPSEUDOSCALARPHOTONPHOTON eff_cp_pseudoscalar_photon_photon
#define EFFCPPSEUDOSCALARGLUONGLUON eff_cp_pseudoscalar_gluon_gluon

#define ALPHA_S_MZ qedqcd.displayAlpha(softsusy::ALPHAS)
#define MWPole qedqcd.displayPoleMW()
#define MZPole qedqcd.displayPoleMZ()
#define MTPole qedqcd.displayPoleMt()
#define MBMB qedqcd.displayMbMb()
#define MTauPole qedqcd.displayPoleMtau()
#define MMPole qedqcd.displayPoleMmuon()

namespace flexiblesusy {

const unsigned MRSSM_observables::NUMBER_OF_OBSERVABLES;

MRSSM_observables::MRSSM_observables()

{
}

Eigen::ArrayXd MRSSM_observables::get() const
{
   Eigen::ArrayXd vec(1);

   vec(0) = 0.;

   return vec;
}

std::vector<std::string> MRSSM_observables::get_names()
{
   std::vector<std::string> names(1);

   names[0] = "no observables defined";

   return names;
}

void MRSSM_observables::clear()
{

}

void MRSSM_observables::set(const Eigen::ArrayXd& vec)
{

}

MRSSM_observables calculate_observables(const MRSSM_mass_eigenstates& model,
                                              const softsusy::QedQcd& qedqcd,
                                              const Physical_input& physical_input)
{
   MRSSM_observables observables;

   


   return observables;
}

} // namespace flexiblesusy
