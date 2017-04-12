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

// File generated at Wed 12 Apr 2017 11:17:55

#include "SplitMSSM_effective_couplings.hpp"

#include "effective_couplings.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

using namespace effective_couplings;

#define INPUTPARAMETER(parameter) model.get_input().parameter
#define MODELPARAMETER(parameter) model.get_##parameter()
#define DERIVEDPARAMETER(parameter) model.##parameter()
#define PHASE(parameter) model.get_##parameter()
#define PHYSICAL(parameter) model.get_physical().parameter

SplitMSSM_effective_couplings::SplitMSSM_effective_couplings(
   const SplitMSSM_mass_eigenstates& model_,
   const softsusy::QedQcd& qedqcd_,
   const Physical_input& input_)
   : model(model_), qedqcd(qedqcd_), physical_input(input_)
   , rg_improve(true), include_qcd_corrections(true)
   , Vd(PHYSICAL(Vd)), Ud(PHYSICAL(Ud)), Vu(PHYSICAL(Vu)), Uu(PHYSICAL(Uu)), Ve
      (PHYSICAL(Ve)), Ue(PHYSICAL(Ue)), ZN(PHYSICAL(ZN)), UM(PHYSICAL(UM)), UP(
      PHYSICAL(UP)), ZZ(PHYSICAL(ZZ))

   , eff_CphhVPVP(0), eff_CphhVGVG(0), eff_CpAhVPVP(0), eff_CpAhVGVG(0)

{
}

void SplitMSSM_effective_couplings::calculate_effective_couplings()
{
   const standard_model::Standard_model sm(initialise_SM());

   const double scale = model.get_scale();
   const Eigen::ArrayXd saved_parameters(model.get());

   const double saved_mt = PHYSICAL(MFu(2));
   PHYSICAL(MFu(2)) = qedqcd.displayPoleMt();

   const auto Mhh = PHYSICAL(Mhh);
   run_SM_strong_coupling_to(sm, 0.5 * Mhh);
   calculate_eff_CphhVPVP();
   run_SM_strong_coupling_to(sm, Mhh);
   calculate_eff_CphhVGVG();

   PHYSICAL(MFu(2)) = saved_mt;
   model.set_scale(scale);
   model.set(saved_parameters);

}

void SplitMSSM_effective_couplings::set_model(const SplitMSSM_mass_eigenstates& model_)
{
   model = model_;
   copy_mixing_matrices_from_model();
}

void SplitMSSM_effective_couplings::copy_mixing_matrices_from_model()
{
   Vd = PHYSICAL(Vd);
   Ud = PHYSICAL(Ud);
   Vu = PHYSICAL(Vu);
   Uu = PHYSICAL(Uu);
   Ve = PHYSICAL(Ve);
   Ue = PHYSICAL(Ue);
   ZN = PHYSICAL(ZN);
   UM = PHYSICAL(UM);
   UP = PHYSICAL(UP);
   ZZ = PHYSICAL(ZZ);

}

standard_model::Standard_model SplitMSSM_effective_couplings::initialise_SM() const
{
   standard_model::Standard_model sm;

   sm.set_loops(2);
   sm.set_thresholds(2);
   sm.set_physical_input(physical_input);

   sm.initialise_from_input(qedqcd);

   return sm;
}

void SplitMSSM_effective_couplings::run_SM_strong_coupling_to(standard_model::Standard_model sm, double m)
{
   sm.run_to(m);

   model.set_g3(sm.get_g3());

}

std::complex<double> SplitMSSM_effective_couplings::scalar_scalar_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      if (m_loop > m_decay) {
         result = 1 + 0.06754745576155852*Sqr(g3);
      }

   }

   return result;
}

std::complex<double> SplitMSSM_effective_couplings::scalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) *
         scalar_diphoton_fermion_loop(m_decay, m_loop);

   }

   return result;
}

std::complex<double> SplitMSSM_effective_couplings::pseudoscalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) *
         pseudoscalar_diphoton_fermion_loop(m_decay, m_loop);

   }

   return result;
}

double SplitMSSM_effective_couplings::number_of_active_flavours(double m) const
{
   if (m < qedqcd.displayMbMb()) {
      return 4.0;
   } else if (m < qedqcd.displayPoleMt()) {
      return 5.0;
   } else {
      return 6.0;
   }
}

double SplitMSSM_effective_couplings::scalar_scaling_factor(double m) const
{
   const double Nf = number_of_active_flavours(m);
   const double mtpole = qedqcd.displayPoleMt();
   const double l = Log(Sqr(m) / Sqr(mtpole));

   const auto g3 = MODELPARAMETER(g3);

   const double nlo_qcd = 0.025330295910584444*(23.75 - 1.1666666666666667*Nf)*
      Sqr(g3);
   const double nnlo_qcd = 0.000641623890917771*Power(g3,4)*(370.1956513893174
      + 2.375*l + (-47.18640261449638 + 0.6666666666666666*l)*Nf +
      0.9017702481178881*Sqr(Nf));
   const double nnnlo_qcd = 0.000016252523020247696*Power(g3,6)*(467.683620788
      + 122.440972222*l + 10.9409722222*Sqr(l));

   return Sqrt(1.0 + nlo_qcd + nnlo_qcd + nnnlo_qcd);
}

double SplitMSSM_effective_couplings::pseudoscalar_scaling_factor(double m) const
{
   const double Nf = number_of_active_flavours(m);
   const double mtpole = qedqcd.displayPoleMt();
   const double l = Log(Sqr(m) / Sqr(mtpole));

   const auto g3 = MODELPARAMETER(g3);

   const double nlo_qcd = 0.025330295910584444*(24.25 - 1.1666666666666667*Nf)*
      Sqr(g3);
   const double nnlo_qcd = 0.000641623890917771*Power(g3,4)*(171.54400563089382
      + 5*l);
   const double nnnlo_qcd = 0;

   return Sqrt(1.0 + nlo_qcd + nnlo_qcd + nnnlo_qcd);
}

double SplitMSSM_effective_couplings::get_hhVPVP_partial_width() const
{
   const double mass = PHYSICAL(Mhh);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CphhVPVP);
}

double SplitMSSM_effective_couplings::get_hhVGVG_partial_width() const
{
   const double mass = PHYSICAL(Mhh);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CphhVGVG);
}

double SplitMSSM_effective_couplings::get_AhVPVP_partial_width() const
{
   const double mass = PHYSICAL(MAh);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CpAhVPVP);
}

double SplitMSSM_effective_couplings::get_AhVGVG_partial_width() const
{
   const double mass = PHYSICAL(MAh);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CpAhVGVG);
}

std::complex<double> SplitMSSM_effective_couplings::CpFdhhbarFdPL(unsigned gt1, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_617;
   std::complex<double> tmp_618;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_619;
      std::complex<double> tmp_620;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_620 += Conj(Ud(gt3,j1))*Yd(j1,j2);
      }
      tmp_619 += tmp_620;
      tmp_618 += (Conj(Vd(gt1,j2))) * tmp_619;
   }
   tmp_617 += tmp_618;
   result += (-0.7071067811865475) * tmp_617;

   return result;
}

std::complex<double> SplitMSSM_effective_couplings::CpFuhhbarFuPL(unsigned gt1, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_621;
   std::complex<double> tmp_622;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_623;
      std::complex<double> tmp_624;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_624 += Conj(Uu(gt3,j1))*Yu(j1,j2);
      }
      tmp_623 += tmp_624;
      tmp_622 += (Conj(Vu(gt1,j2))) * tmp_623;
   }
   tmp_621 += tmp_622;
   result += (-0.7071067811865475) * tmp_621;

   return result;
}

std::complex<double> SplitMSSM_effective_couplings::CpFehhbarFePL(unsigned gt1, unsigned gt3) const
{
   const auto Ye = MODELPARAMETER(Ye);

   std::complex<double> result;

   std::complex<double> tmp_625;
   std::complex<double> tmp_626;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_627;
      std::complex<double> tmp_628;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_628 += Conj(Ue(gt3,j1))*Ye(j1,j2);
      }
      tmp_627 += tmp_628;
      tmp_626 += (Conj(Ve(gt1,j2))) * tmp_627;
   }
   tmp_625 += tmp_626;
   result += (-0.7071067811865475) * tmp_625;

   return result;
}

std::complex<double> SplitMSSM_effective_couplings::CpChahhbarChaPL(unsigned gt1, unsigned gt3) const
{
   const auto g2d = MODELPARAMETER(g2d);
   const auto g2u = MODELPARAMETER(g2u);

   std::complex<double> result;

   result = -0.7071067811865475*(g2d*Conj(UM(gt1,1))*Conj(UP(gt3,0)) + g2u*Conj
      (UM(gt1,0))*Conj(UP(gt3,1)));

   return result;
}

double SplitMSSM_effective_couplings::CphhVWpconjVWp() const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);

   double result = 0.0;

   result = 0.7071067811865475*v*Sqr(g2);

   return result;
}

std::complex<double> SplitMSSM_effective_couplings::CpAhFdbarFdPL(unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_629;
   std::complex<double> tmp_630;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_631;
      std::complex<double> tmp_632;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_632 += Conj(Ud(gt3,j1))*Yd(j1,j2);
      }
      tmp_631 += tmp_632;
      tmp_630 += (Conj(Vd(gt2,j2))) * tmp_631;
   }
   tmp_629 += tmp_630;
   result += (std::complex<double>(0.,0.7071067811865475)) * tmp_629;

   return result;
}

std::complex<double> SplitMSSM_effective_couplings::CpAhFubarFuPL(unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_633;
   std::complex<double> tmp_634;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_635;
      std::complex<double> tmp_636;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_636 += Conj(Uu(gt3,j1))*Yu(j1,j2);
      }
      tmp_635 += tmp_636;
      tmp_634 += (Conj(Vu(gt2,j2))) * tmp_635;
   }
   tmp_633 += tmp_634;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_633;

   return result;
}

std::complex<double> SplitMSSM_effective_couplings::CpAhFebarFePL(unsigned gt2, unsigned gt3) const
{
   const auto Ye = MODELPARAMETER(Ye);

   std::complex<double> result;

   std::complex<double> tmp_637;
   std::complex<double> tmp_638;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_639;
      std::complex<double> tmp_640;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_640 += Conj(Ue(gt3,j1))*Ye(j1,j2);
      }
      tmp_639 += tmp_640;
      tmp_638 += (Conj(Ve(gt2,j2))) * tmp_639;
   }
   tmp_637 += tmp_638;
   result += (std::complex<double>(0.,0.7071067811865475)) * tmp_637;

   return result;
}

std::complex<double> SplitMSSM_effective_couplings::CpAhChabarChaPL(unsigned gt2, unsigned gt3) const
{
   const auto g2d = MODELPARAMETER(g2d);
   const auto g2u = MODELPARAMETER(g2u);

   std::complex<double> result;

   result = std::complex<double>(0.,-0.7071067811865475)*(g2d*Conj(UM(gt2,1))*
      Conj(UP(gt3,0)) - g2u*Conj(UM(gt2,0))*Conj(UP(gt3,1)));

   return result;
}

void SplitMSSM_effective_couplings::calculate_eff_CphhVPVP()
{
   const auto MFd = PHYSICAL(MFd);
   const auto MFu = PHYSICAL(MFu);
   const auto MFe = PHYSICAL(MFe);
   const auto MCha = PHYSICAL(MCha);
   const auto MVWp = PHYSICAL(MVWp);
   const auto decay_mass = PHYSICAL(Mhh);
   const auto decay_scale = 0.25 * Sqr(decay_mass);

   const auto vev = 1.0 / Sqrt(qedqcd.displayFermiConstant() * Sqrt(2.0));

   std::complex<double> result = 0;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += 0.3333333333333333 * scalar_fermion_qcd_factor(decay_mass,
         MFd(gI1)) * CpFdhhbarFdPL(gI1, gI1) * vev * AS12(decay_scale / Sqr(MFd(
         gI1))) / MFd(gI1);
   }
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += 1.3333333333333333 * scalar_fermion_qcd_factor(decay_mass,
         MFu(gI1)) * CpFuhhbarFuPL(gI1, gI1) * vev * AS12(decay_scale / Sqr(MFu(
         gI1))) / MFu(gI1);
   }
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += CpFehhbarFePL(gI1, gI1) * vev * AS12(decay_scale / Sqr(MFe(
         gI1))) / MFe(gI1);
   }
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      result += CpChahhbarChaPL(gI1, gI1) * vev * AS12(decay_scale / Sqr(
         MCha(gI1))) / MCha(gI1);
   }
   result += -0.5 * CphhVWpconjVWp() * vev * AS1(decay_scale / Sqr(MVWp)) / Sqr
      (MVWp);


   result *= 0.1892681907127351 * physical_input.get(Physical_input::alpha_em_0
      ) * Sqrt(qedqcd.displayFermiConstant());

   eff_CphhVPVP = result;

}

void SplitMSSM_effective_couplings::calculate_eff_CphhVGVG()
{
   const auto g3 = MODELPARAMETER(g3);
   const auto MFd = PHYSICAL(MFd);
   const auto MFu = PHYSICAL(MFu);
   const double alpha_s = 0.07957747154594767*Sqr(g3);
   const auto decay_mass = PHYSICAL(Mhh);
   const auto decay_scale = 0.25 * Sqr(decay_mass);

   const auto vev = 1.0 / Sqrt(qedqcd.displayFermiConstant() * Sqrt(2.0));

   std::complex<double> result = 0;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += CpFdhhbarFdPL(gI1, gI1) * vev * AS12(decay_scale / Sqr(MFd(
         gI1))) / MFd(gI1);
   }
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += CpFuhhbarFuPL(gI1, gI1) * vev * AS12(decay_scale / Sqr(MFu(
         gI1))) / MFu(gI1);
   }
   result *= std::complex<double>(0.75,0.);

   if (include_qcd_corrections) {
      result *= scalar_scaling_factor(decay_mass);
   }


   result *= 0.12617879380849006 * alpha_s * Sqrt(qedqcd.displayFermiConstant()
      );

   eff_CphhVGVG = result;

}

void SplitMSSM_effective_couplings::calculate_eff_CpAhVPVP()
{
   const auto MFd = PHYSICAL(MFd);
   const auto MFu = PHYSICAL(MFu);
   const auto MFe = PHYSICAL(MFe);
   const auto MCha = PHYSICAL(MCha);
   const auto decay_mass = PHYSICAL(MAh);
   const auto decay_scale = 0.25 * Sqr(decay_mass);

   const auto vev = 1.0 / Sqrt(qedqcd.displayFermiConstant() * Sqrt(2.0));

   std::complex<double> result = 0;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += 0.3333333333333333 * pseudoscalar_fermion_qcd_factor(
         decay_mass, MFd(gI1)) * CpAhFdbarFdPL(gI1, gI1) * vev * AP12(decay_scale
         / Sqr(MFd(gI1))) / MFd(gI1);
   }
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += 1.3333333333333333 * pseudoscalar_fermion_qcd_factor(
         decay_mass, MFu(gI1)) * CpAhFubarFuPL(gI1, gI1) * vev * AP12(decay_scale
         / Sqr(MFu(gI1))) / MFu(gI1);
   }
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += CpAhFebarFePL(gI1, gI1) * vev * AP12(decay_scale / Sqr(MFe(
         gI1))) / MFe(gI1);
   }
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      result += CpAhChabarChaPL(gI1, gI1) * vev * AP12(decay_scale / Sqr(
         MCha(gI1))) / MCha(gI1);
   }
   result *= std::complex<double>(2.0,0.);


   result *= 0.1892681907127351 * physical_input.get(Physical_input::alpha_em_0
      ) * Sqrt(qedqcd.displayFermiConstant());

   eff_CpAhVPVP = result;

}

void SplitMSSM_effective_couplings::calculate_eff_CpAhVGVG()
{
   const auto g3 = MODELPARAMETER(g3);
   const auto MFd = PHYSICAL(MFd);
   const auto MFu = PHYSICAL(MFu);
   const double alpha_s = 0.07957747154594767*Sqr(g3);
   const auto decay_mass = PHYSICAL(MAh);
   const auto decay_scale = 0.25 * Sqr(decay_mass);

   const auto vev = 1.0 / Sqrt(qedqcd.displayFermiConstant() * Sqrt(2.0));

   std::complex<double> result = 0;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += CpAhFdbarFdPL(gI1, gI1) * vev * AP12(decay_scale / Sqr(MFd(
         gI1))) / MFd(gI1);
   }
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += CpAhFubarFuPL(gI1, gI1) * vev * AP12(decay_scale / Sqr(MFu(
         gI1))) / MFu(gI1);
   }
   result *= std::complex<double>(1.5,0.);

   if (include_qcd_corrections) {
      result *= pseudoscalar_scaling_factor(decay_mass);
   }


   result *= 0.12617879380849006 * alpha_s * Sqrt(qedqcd.displayFermiConstant()
      );

   eff_CpAhVGVG = result;

}


} // namespace flexiblesusy
