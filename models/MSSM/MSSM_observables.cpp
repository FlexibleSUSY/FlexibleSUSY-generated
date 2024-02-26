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


#include "MSSM_observables.hpp"
#include "MSSM_mass_eigenstates.hpp"
#include "MSSM_amm.hpp"
#include "MSSM_edm.hpp"
#include "MSSM_b_to_s_gamma.hpp"
#include "observables/l_to_l_conversion/settings.hpp"
#include "observables/MSSM_br_l_to_3l.hpp"
#include "observables/MSSM_br_l_to_l_gamma.hpp"
#include "observables/MSSM_l_to_l_conversion.hpp"
#include "cxx_qft/MSSM_qft.hpp"
#include "config.h"
#include "eigen_utils.hpp"
#include "numerics2.hpp"
#include "wrappers.hpp"
#include "lowe.h"
#include "physical_input.hpp"

#ifdef ENABLE_GM2CALC
#include "gm2calc_interface.hpp"
#endif

#define MODEL model
#define AMM0(p) amm_ ## p
#define AMM1(p,idx) amm_ ## p ## _ ## idx
#define AMMUNCERTAINTY0(p) amm_uncertainty_ ## p
#define AMMUNCERTAINTY1(p,idx) amm_uncertainty_ ## p ## _ ## idx
#define AMUGM2CALC a_muon_gm2calc
#define AMUGM2CALCUNCERTAINTY a_muon_gm2calc_uncertainty
#define DERIVEDPARAMETER(p) model.p()
#define EXTRAPARAMETER(p) model.get_##p()
#define INPUTPARAMETER(p) model.get_input().p
#define MODELPARAMETER(p) model.get_##p()
#define PHASE(p) model.get_##p()
#define LowEnergyConstant(p) Electroweak_constants::p
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model.get_physical().p
#define EDM0(p) edm_ ## p
#define EDM1(p,idx) edm_ ## p ## _ ## idx
#define BSGAMMA b_to_s_gamma

#define ALPHA_EM_MZ qedqcd.displayAlpha(softsusy::ALPHA)
#define ALPHA_EM_0 physical_input.get(Physical_input::alpha_em_0)
#define ALPHA_S_MZ qedqcd.displayAlpha(softsusy::ALPHAS)
#define MHPole physical_input.get(Physical_input::mh_pole)
#define MWPole qedqcd.displayPoleMW()
#define MZPole qedqcd.displayPoleMZ()
#define MU2GeV qedqcd.displayMu2GeV()
#define MS2GeV qedqcd.displayMs2GeV()
#define MTPole qedqcd.displayPoleMt()
#define MD2GeV qedqcd.displayMd2GeV()
#define MCMC qedqcd.displayMcMc()
#define MBMB qedqcd.displayMbMb()
#define Mv1Pole qedqcd.displayNeutrinoPoleMass(1)
#define Mv2Pole qedqcd.displayNeutrinoPoleMass(2)
#define Mv3Pole qedqcd.displayNeutrinoPoleMass(3)
#define MEPole qedqcd.displayPoleMel()
#define MMPole qedqcd.displayPoleMmuon()
#define MTauPole qedqcd.displayPoleMtau()
#define CKMInput qedqcd.get_complex_ckm()

namespace flexiblesusy {

const int MSSM_observables::NUMBER_OF_OBSERVABLES;

MSSM_observables::MSSM_observables()
   : amm_Fe_0(0)
   , amm_Fe_1(0)
   , amm_Fe_2(0)
   , edm_Fe_0(0)
   , edm_Fe_1(0)
   , edm_Fe_2(0)
   , Fe1_to_Fe0_VP(0)
   , Fe2to11bar1_All_1loop(Eigen::Array<std::complex<double>,13,1>::Zero())
   , Fe2Fe1inAl_All_1loop(Eigen::Array<std::complex<double>,13,1>::Zero())

{
}

Eigen::ArrayXd MSSM_observables::get() const
{
   Eigen::ArrayXd vec(MSSM_observables::NUMBER_OF_OBSERVABLES);

   vec(0) = amm_Fe_0;
   vec(1) = amm_Fe_1;
   vec(2) = amm_Fe_2;
   vec(3) = edm_Fe_0;
   vec(4) = edm_Fe_1;
   vec(5) = edm_Fe_2;
   vec(6) = Fe1_to_Fe0_VP;
   vec(7) = Re(Fe2to11bar1_All_1loop(0));
   vec(8) = Im(Fe2to11bar1_All_1loop(0));
   vec(9) = Re(Fe2to11bar1_All_1loop(1));
   vec(10) = Im(Fe2to11bar1_All_1loop(1));
   vec(11) = Re(Fe2to11bar1_All_1loop(2));
   vec(12) = Im(Fe2to11bar1_All_1loop(2));
   vec(13) = Re(Fe2to11bar1_All_1loop(3));
   vec(14) = Im(Fe2to11bar1_All_1loop(3));
   vec(15) = Re(Fe2to11bar1_All_1loop(4));
   vec(16) = Im(Fe2to11bar1_All_1loop(4));
   vec(17) = Re(Fe2to11bar1_All_1loop(5));
   vec(18) = Im(Fe2to11bar1_All_1loop(5));
   vec(19) = Re(Fe2to11bar1_All_1loop(6));
   vec(20) = Im(Fe2to11bar1_All_1loop(6));
   vec(21) = Re(Fe2to11bar1_All_1loop(7));
   vec(22) = Im(Fe2to11bar1_All_1loop(7));
   vec(23) = Re(Fe2to11bar1_All_1loop(8));
   vec(24) = Im(Fe2to11bar1_All_1loop(8));
   vec(25) = Re(Fe2to11bar1_All_1loop(9));
   vec(26) = Im(Fe2to11bar1_All_1loop(9));
   vec(27) = Re(Fe2to11bar1_All_1loop(10));
   vec(28) = Im(Fe2to11bar1_All_1loop(10));
   vec(29) = Re(Fe2to11bar1_All_1loop(11));
   vec(30) = Im(Fe2to11bar1_All_1loop(11));
   vec(31) = Re(Fe2to11bar1_All_1loop(12));
   vec(32) = Im(Fe2to11bar1_All_1loop(12));
   vec(33) = Re(Fe2Fe1inAl_All_1loop(0));
   vec(34) = Im(Fe2Fe1inAl_All_1loop(0));
   vec(35) = Re(Fe2Fe1inAl_All_1loop(1));
   vec(36) = Im(Fe2Fe1inAl_All_1loop(1));
   vec(37) = Re(Fe2Fe1inAl_All_1loop(2));
   vec(38) = Im(Fe2Fe1inAl_All_1loop(2));
   vec(39) = Re(Fe2Fe1inAl_All_1loop(3));
   vec(40) = Im(Fe2Fe1inAl_All_1loop(3));
   vec(41) = Re(Fe2Fe1inAl_All_1loop(4));
   vec(42) = Im(Fe2Fe1inAl_All_1loop(4));
   vec(43) = Re(Fe2Fe1inAl_All_1loop(5));
   vec(44) = Im(Fe2Fe1inAl_All_1loop(5));
   vec(45) = Re(Fe2Fe1inAl_All_1loop(6));
   vec(46) = Im(Fe2Fe1inAl_All_1loop(6));
   vec(47) = Re(Fe2Fe1inAl_All_1loop(7));
   vec(48) = Im(Fe2Fe1inAl_All_1loop(7));
   vec(49) = Re(Fe2Fe1inAl_All_1loop(8));
   vec(50) = Im(Fe2Fe1inAl_All_1loop(8));
   vec(51) = Re(Fe2Fe1inAl_All_1loop(9));
   vec(52) = Im(Fe2Fe1inAl_All_1loop(9));
   vec(53) = Re(Fe2Fe1inAl_All_1loop(10));
   vec(54) = Im(Fe2Fe1inAl_All_1loop(10));
   vec(55) = Re(Fe2Fe1inAl_All_1loop(11));
   vec(56) = Im(Fe2Fe1inAl_All_1loop(11));
   vec(57) = Re(Fe2Fe1inAl_All_1loop(12));
   vec(58) = Im(Fe2Fe1inAl_All_1loop(12));

   return vec;
}

std::vector<std::string> MSSM_observables::get_names()
{
   std::vector<std::string> names(MSSM_observables::NUMBER_OF_OBSERVABLES);

   names[0] = "amm_Fe_0";
   names[1] = "amm_Fe_1";
   names[2] = "amm_Fe_2";
   names[3] = "edm_Fe_0";
   names[4] = "edm_Fe_1";
   names[5] = "edm_Fe_2";
   names[6] = "Fe1_to_Fe0_VP";
   names[7] = "Re(Fe2to11bar1_All_1loop(0))";
   names[8] = "Im(Fe2to11bar1_All_1loop(0))";
   names[9] = "Re(Fe2to11bar1_All_1loop(1))";
   names[10] = "Im(Fe2to11bar1_All_1loop(1))";
   names[11] = "Re(Fe2to11bar1_All_1loop(2))";
   names[12] = "Im(Fe2to11bar1_All_1loop(2))";
   names[13] = "Re(Fe2to11bar1_All_1loop(3))";
   names[14] = "Im(Fe2to11bar1_All_1loop(3))";
   names[15] = "Re(Fe2to11bar1_All_1loop(4))";
   names[16] = "Im(Fe2to11bar1_All_1loop(4))";
   names[17] = "Re(Fe2to11bar1_All_1loop(5))";
   names[18] = "Im(Fe2to11bar1_All_1loop(5))";
   names[19] = "Re(Fe2to11bar1_All_1loop(6))";
   names[20] = "Im(Fe2to11bar1_All_1loop(6))";
   names[21] = "Re(Fe2to11bar1_All_1loop(7))";
   names[22] = "Im(Fe2to11bar1_All_1loop(7))";
   names[23] = "Re(Fe2to11bar1_All_1loop(8))";
   names[24] = "Im(Fe2to11bar1_All_1loop(8))";
   names[25] = "Re(Fe2to11bar1_All_1loop(9))";
   names[26] = "Im(Fe2to11bar1_All_1loop(9))";
   names[27] = "Re(Fe2to11bar1_All_1loop(10))";
   names[28] = "Im(Fe2to11bar1_All_1loop(10))";
   names[29] = "Re(Fe2to11bar1_All_1loop(11))";
   names[30] = "Im(Fe2to11bar1_All_1loop(11))";
   names[31] = "Re(Fe2to11bar1_All_1loop(12))";
   names[32] = "Im(Fe2to11bar1_All_1loop(12))";
   names[33] = "Re(Fe2Fe1inAl_All_1loop(0))";
   names[34] = "Im(Fe2Fe1inAl_All_1loop(0))";
   names[35] = "Re(Fe2Fe1inAl_All_1loop(1))";
   names[36] = "Im(Fe2Fe1inAl_All_1loop(1))";
   names[37] = "Re(Fe2Fe1inAl_All_1loop(2))";
   names[38] = "Im(Fe2Fe1inAl_All_1loop(2))";
   names[39] = "Re(Fe2Fe1inAl_All_1loop(3))";
   names[40] = "Im(Fe2Fe1inAl_All_1loop(3))";
   names[41] = "Re(Fe2Fe1inAl_All_1loop(4))";
   names[42] = "Im(Fe2Fe1inAl_All_1loop(4))";
   names[43] = "Re(Fe2Fe1inAl_All_1loop(5))";
   names[44] = "Im(Fe2Fe1inAl_All_1loop(5))";
   names[45] = "Re(Fe2Fe1inAl_All_1loop(6))";
   names[46] = "Im(Fe2Fe1inAl_All_1loop(6))";
   names[47] = "Re(Fe2Fe1inAl_All_1loop(7))";
   names[48] = "Im(Fe2Fe1inAl_All_1loop(7))";
   names[49] = "Re(Fe2Fe1inAl_All_1loop(8))";
   names[50] = "Im(Fe2Fe1inAl_All_1loop(8))";
   names[51] = "Re(Fe2Fe1inAl_All_1loop(9))";
   names[52] = "Im(Fe2Fe1inAl_All_1loop(9))";
   names[53] = "Re(Fe2Fe1inAl_All_1loop(10))";
   names[54] = "Im(Fe2Fe1inAl_All_1loop(10))";
   names[55] = "Re(Fe2Fe1inAl_All_1loop(11))";
   names[56] = "Im(Fe2Fe1inAl_All_1loop(11))";
   names[57] = "Re(Fe2Fe1inAl_All_1loop(12))";
   names[58] = "Im(Fe2Fe1inAl_All_1loop(12))";

   return names;
}

void MSSM_observables::clear()
{
   amm_Fe_0 = 0.;
   amm_Fe_1 = 0.;
   amm_Fe_2 = 0.;
   edm_Fe_0 = 0.;
   edm_Fe_1 = 0.;
   edm_Fe_2 = 0.;
   Fe1_to_Fe0_VP = 0.;
   Fe2to11bar1_All_1loop = Eigen::Array<std::complex<double>,13,1>::Zero();
   Fe2Fe1inAl_All_1loop = Eigen::Array<std::complex<double>,13,1>::Zero();

}

void MSSM_observables::set(const Eigen::ArrayXd& vec)
{
   assert(vec.rows() == MSSM_observables::NUMBER_OF_OBSERVABLES);

   amm_Fe_0 = vec(0);
   amm_Fe_1 = vec(1);
   amm_Fe_2 = vec(2);
   edm_Fe_0 = vec(3);
   edm_Fe_1 = vec(4);
   edm_Fe_2 = vec(5);
   Fe1_to_Fe0_VP = vec(6);
   Fe2to11bar1_All_1loop(0) = std::complex<double>(vec(7), vec(8));
   Fe2to11bar1_All_1loop(1) = std::complex<double>(vec(9), vec(10));
   Fe2to11bar1_All_1loop(2) = std::complex<double>(vec(11), vec(12));
   Fe2to11bar1_All_1loop(3) = std::complex<double>(vec(13), vec(14));
   Fe2to11bar1_All_1loop(4) = std::complex<double>(vec(15), vec(16));
   Fe2to11bar1_All_1loop(5) = std::complex<double>(vec(17), vec(18));
   Fe2to11bar1_All_1loop(6) = std::complex<double>(vec(19), vec(20));
   Fe2to11bar1_All_1loop(7) = std::complex<double>(vec(21), vec(22));
   Fe2to11bar1_All_1loop(8) = std::complex<double>(vec(23), vec(24));
   Fe2to11bar1_All_1loop(9) = std::complex<double>(vec(25), vec(26));
   Fe2to11bar1_All_1loop(10) = std::complex<double>(vec(27), vec(28));
   Fe2to11bar1_All_1loop(11) = std::complex<double>(vec(29), vec(30));
   Fe2to11bar1_All_1loop(12) = std::complex<double>(vec(31), vec(32));
   Fe2Fe1inAl_All_1loop(0) = std::complex<double>(vec(33), vec(34));
   Fe2Fe1inAl_All_1loop(1) = std::complex<double>(vec(35), vec(36));
   Fe2Fe1inAl_All_1loop(2) = std::complex<double>(vec(37), vec(38));
   Fe2Fe1inAl_All_1loop(3) = std::complex<double>(vec(39), vec(40));
   Fe2Fe1inAl_All_1loop(4) = std::complex<double>(vec(41), vec(42));
   Fe2Fe1inAl_All_1loop(5) = std::complex<double>(vec(43), vec(44));
   Fe2Fe1inAl_All_1loop(6) = std::complex<double>(vec(45), vec(46));
   Fe2Fe1inAl_All_1loop(7) = std::complex<double>(vec(47), vec(48));
   Fe2Fe1inAl_All_1loop(8) = std::complex<double>(vec(49), vec(50));
   Fe2Fe1inAl_All_1loop(9) = std::complex<double>(vec(51), vec(52));
   Fe2Fe1inAl_All_1loop(10) = std::complex<double>(vec(53), vec(54));
   Fe2Fe1inAl_All_1loop(11) = std::complex<double>(vec(55), vec(56));
   Fe2Fe1inAl_All_1loop(12) = std::complex<double>(vec(57), vec(58));

}

MSSM_observables calculate_observables(const MSSM_mass_eigenstates& model,
                                              const softsusy::QedQcd& qedqcd,
                                              const LToLConversion_settings& ltolconversion_settings,
                                              const Physical_input& physical_input,
                                              const Spectrum_generator_settings& settings,
                                              double scale)
{
   auto model_at_scale = model;

   if (scale > 0.) {
      try {
         model_at_scale.run_to(scale);
      } catch (const NonPerturbativeRunningError& e) {
         MSSM_observables observables;
         observables.problems.general.flag_non_perturbative_running(scale);
         return observables;
      } catch (const Error& e) {
         MSSM_observables observables;
         observables.problems.general.flag_thrown(e.what());
         return observables;
      } catch (const std::exception& e) {
         MSSM_observables observables;
         observables.problems.general.flag_thrown(e.what());
         return observables;
      }
   }

   return calculate_observables(model_at_scale,
                                qedqcd,
                                ltolconversion_settings,
                                physical_input,
                                settings);
}

MSSM_observables calculate_observables(const MSSM_mass_eigenstates& model,
                                              const softsusy::QedQcd& qedqcd,
                                              const LToLConversion_settings& ltolconversion_settings,
                                              const Physical_input& physical_input,
                                              const Spectrum_generator_settings& settings)
{
   MSSM_observables observables;

   try {
      
      observables.AMM1(Fe, 0) = MSSM_amm::calculate_amm<MSSM_cxx_diagrams::fields::Fe>(MODEL, qedqcd, settings,0);
      observables.AMM1(Fe, 1) = MSSM_amm::calculate_amm<MSSM_cxx_diagrams::fields::Fe>(MODEL, qedqcd, settings,1);
      observables.AMM1(Fe, 2) = MSSM_amm::calculate_amm<MSSM_cxx_diagrams::fields::Fe>(MODEL, qedqcd, settings,2);
      observables.EDM1(Fe, 0) = MSSM_edm::calculate_edm<MSSM_cxx_diagrams::fields::Fe>(MODEL, qedqcd, 0);
      observables.EDM1(Fe, 1) = MSSM_edm::calculate_edm<MSSM_cxx_diagrams::fields::Fe>(MODEL, qedqcd, 1);
      observables.EDM1(Fe, 2) = MSSM_edm::calculate_edm<MSSM_cxx_diagrams::fields::Fe>(MODEL, qedqcd, 2);
      observables.Fe1_to_Fe0_VP = MSSM_br_l_to_l_gamma::calculate_Fe_to_Fe_VP(1, 0, model, qedqcd, physical_input);
      observables.Fe2to11bar1_All_1loop = MSSM_br_l_to_3l::calculate_Fe_to_FeFebarFe_for_All_1loop(1, 0, 0, model, qedqcd);
      observables.Fe2Fe1inAl_All_1loop = MSSM_l_to_l_conversion::calculate_FeFe_forAll_1loop(1, 0, MSSM_l_to_l_conversion::Nucleus::Al, model, ltolconversion_settings, qedqcd);
   } catch (const NonPerturbativeRunningError& e) {
      observables.problems.general.flag_non_perturbative_running(e.get_scale());
   } catch (const Error& e) {
      observables.problems.general.flag_thrown(e.what());
   } catch (const std::exception& e) {
      observables.problems.general.flag_thrown(e.what());
   }

   return observables;
}

} // namespace flexiblesusy
