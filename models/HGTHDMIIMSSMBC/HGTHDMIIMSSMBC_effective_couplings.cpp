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

// File generated at Sat 15 Oct 2016 15:28:46

#include "HGTHDMIIMSSMBC_effective_couplings.hpp"

#include "effective_couplings.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

using namespace effective_couplings;

#define INPUTPARAMETER(parameter) model.get_input().parameter
#define MODELPARAMETER(parameter) model.get_##parameter()
#define DERIVEDPARAMETER(parameter) model.##parameter()
#define PHASE(parameter) model.get_##parameter()
#define PHYSICAL(parameter) model.get_physical().parameter

HGTHDMIIMSSMBC_effective_couplings::HGTHDMIIMSSMBC_effective_couplings(
   const HGTHDMIIMSSMBC_mass_eigenstates& model_,
   const softsusy::QedQcd& qedqcd_,
   const Physical_input& input_)
   : model(model_), qedqcd(qedqcd_), physical_input(input_)
   , rg_improve(true), include_qcd_corrections(true)
   , ZH(PHYSICAL(ZH)), ZA(PHYSICAL(ZA)), ZP(PHYSICAL(ZP)), Vd(PHYSICAL(Vd)), Ud
      (PHYSICAL(Ud)), Vu(PHYSICAL(Vu)), Uu(PHYSICAL(Uu)), Ve(PHYSICAL(Ve)), Ue(
      PHYSICAL(Ue)), ZN(PHYSICAL(ZN)), UM(PHYSICAL(UM)), UP(PHYSICAL(UP)), ZZ(
      PHYSICAL(ZZ))

   , eff_CphhVPVP(Eigen::Array<std::complex<double>,2,1>::Zero()), eff_CphhVGVG
      (Eigen::Array<std::complex<double>,2,1>::Zero()), eff_CpAhVPVP(Eigen::Array<
      std::complex<double>,2,1>::Zero()), eff_CpAhVGVG(Eigen::Array<std::complex<
      double>,2,1>::Zero())

{
}

void HGTHDMIIMSSMBC_effective_couplings::calculate_effective_couplings()
{
   const standard_model::Standard_model sm(initialise_SM());

   const double scale = model.get_scale();
   const Eigen::ArrayXd saved_parameters(model.get());

   const double saved_mt = PHYSICAL(MFu(2));
   PHYSICAL(MFu(2)) = qedqcd.displayPoleMt();

   const auto Mhh = PHYSICAL(Mhh);
   for (unsigned gO1 = 0; gO1 < 2; ++gO1) {
      run_SM_strong_coupling_to(sm, 0.5 * Mhh(gO1));
      calculate_eff_CphhVPVP(gO1);
      run_SM_strong_coupling_to(sm, Mhh(gO1));
      calculate_eff_CphhVGVG(gO1);
   }

   const auto MAh = PHYSICAL(MAh);
   for (unsigned gO1 = 1; gO1 < 2; ++gO1) {
      run_SM_strong_coupling_to(sm, 0.5 * MAh(gO1));
      calculate_eff_CpAhVPVP(gO1);
      run_SM_strong_coupling_to(sm, MAh(gO1));
      calculate_eff_CpAhVGVG(gO1);
   }

   PHYSICAL(MFu(2)) = saved_mt;
   model.set_scale(scale);
   model.set(saved_parameters);

}

void HGTHDMIIMSSMBC_effective_couplings::set_model(const HGTHDMIIMSSMBC_mass_eigenstates& model_)
{
   model = model_;
   copy_mixing_matrices_from_model();
}

void HGTHDMIIMSSMBC_effective_couplings::copy_mixing_matrices_from_model()
{
   ZH = PHYSICAL(ZH);
   ZA = PHYSICAL(ZA);
   ZP = PHYSICAL(ZP);
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

standard_model::Standard_model HGTHDMIIMSSMBC_effective_couplings::initialise_SM() const
{
   standard_model::Standard_model sm;

   sm.set_loops(2);
   sm.set_thresholds(2);
   sm.set_low_energy_data(qedqcd);
   sm.set_physical_input(physical_input);

   sm.initialise_from_input();

   return sm;
}

void HGTHDMIIMSSMBC_effective_couplings::run_SM_strong_coupling_to(standard_model::Standard_model sm, double m)
{
   sm.run_to(m);

   model.set_g3(sm.get_g3());

}

std::complex<double> HGTHDMIIMSSMBC_effective_couplings::scalar_scalar_qcd_factor(double m_decay, double m_loop) const
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

std::complex<double> HGTHDMIIMSSMBC_effective_couplings::scalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) *
         scalar_diphoton_fermion_loop(m_decay, m_loop);

   }

   return result;
}

std::complex<double> HGTHDMIIMSSMBC_effective_couplings::pseudoscalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) *
         pseudoscalar_diphoton_fermion_loop(m_decay, m_loop);

   }

   return result;
}

double HGTHDMIIMSSMBC_effective_couplings::number_of_active_flavours(double m) const
{
   if (m < qedqcd.displayMbMb()) {
      return 4.0;
   } else if (m < qedqcd.displayPoleMt()) {
      return 5.0;
   } else {
      return 6.0;
   }
}

double HGTHDMIIMSSMBC_effective_couplings::scalar_scaling_factor(double m) const
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

double HGTHDMIIMSSMBC_effective_couplings::pseudoscalar_scaling_factor(double m) const
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

double HGTHDMIIMSSMBC_effective_couplings::get_hhVPVP_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(Mhh)(gO1);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CphhVPVP(gO1));
}

double HGTHDMIIMSSMBC_effective_couplings::get_hhVGVG_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(Mhh)(gO1);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CphhVGVG(gO1));
}

double HGTHDMIIMSSMBC_effective_couplings::get_AhVPVP_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(MAh)(gO1);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CpAhVPVP(gO1));
}

double HGTHDMIIMSSMBC_effective_couplings::get_AhVGVG_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(MAh)(gO1);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CpAhVGVG(gO1));
}

std::complex<double> HGTHDMIIMSSMBC_effective_couplings::CphhHmconjHm(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Lambda1 = MODELPARAMETER(Lambda1);
   const auto Lambda2 = MODELPARAMETER(Lambda2);
   const auto Lambda3 = MODELPARAMETER(Lambda3);
   const auto Lambda4 = MODELPARAMETER(Lambda4);
   const auto Lambda5 = MODELPARAMETER(Lambda5);
   const auto Lambda6 = MODELPARAMETER(Lambda6);
   const auto Lambda7 = MODELPARAMETER(Lambda7);
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);

   std::complex<double> result;

   result = 0.5*(-(ZH(gt1,1)*(ZP(gt2,0)*((Lambda6*v1 + 2*Lambda3*v2 + v1*Conj(
      Lambda6))*ZP(gt3,0) + (Lambda4*v1 + 2*Lambda7*v2 + v1*Conj(Lambda5))*ZP(gt3,
      1)) + ZP(gt2,1)*(((Lambda4 + Lambda5)*v1 + 2*v2*Conj(Lambda7))*ZP(gt3,0) + (
      Lambda7*v1 + 4*Lambda2*v2 + v1*Conj(Lambda7))*ZP(gt3,1)))) - ZH(gt1,0)*(ZP(
      gt2,0)*((4*Lambda1*v1 + Lambda6*v2 + v2*Conj(Lambda6))*ZP(gt3,0) + (2*
      Lambda6*v1 + Lambda4*v2 + v2*Conj(Lambda5))*ZP(gt3,1)) + ZP(gt2,1)*(((
      Lambda4 + Lambda5)*v2 + 2*v1*Conj(Lambda6))*ZP(gt3,0) + (2*Lambda3*v1 +
      Lambda7*v2 + v2*Conj(Lambda7))*ZP(gt3,1))));

   return result;
}

std::complex<double> HGTHDMIIMSSMBC_effective_couplings::CpFdhhbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_763;
   std::complex<double> tmp_764;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_765;
      std::complex<double> tmp_766;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_766 += Conj(Ud(gt3,j1))*Yd(j1,j2);
      }
      tmp_765 += tmp_766;
      tmp_764 += (Conj(Vd(gt1,j2))) * tmp_765;
   }
   tmp_763 += tmp_764;
   result += (-0.7071067811865475*ZH(gt2,0)) * tmp_763;

   return result;
}

std::complex<double> HGTHDMIIMSSMBC_effective_couplings::CpFuhhbarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_767;
   std::complex<double> tmp_768;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_769;
      std::complex<double> tmp_770;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_770 += Conj(Uu(gt3,j1))*Yu(j1,j2);
      }
      tmp_769 += tmp_770;
      tmp_768 += (Conj(Vu(gt1,j2))) * tmp_769;
   }
   tmp_767 += tmp_768;
   result += (0.7071067811865475*ZH(gt2,1)) * tmp_767;

   return result;
}

std::complex<double> HGTHDMIIMSSMBC_effective_couplings::CpFehhbarFePL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Ye = MODELPARAMETER(Ye);

   std::complex<double> result;

   std::complex<double> tmp_771;
   std::complex<double> tmp_772;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_773;
      std::complex<double> tmp_774;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_774 += Conj(Ue(gt3,j1))*Ye(j1,j2);
      }
      tmp_773 += tmp_774;
      tmp_772 += (Conj(Ve(gt1,j2))) * tmp_773;
   }
   tmp_771 += tmp_772;
   result += (-0.7071067811865475*ZH(gt2,0)) * tmp_771;

   return result;
}

std::complex<double> HGTHDMIIMSSMBC_effective_couplings::CpChahhbarChaPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g1d = MODELPARAMETER(g1d);
   const auto g2u = MODELPARAMETER(g2u);

   std::complex<double> result;

   result = -0.7071067811865475*(g1d*Conj(UM(gt1,1))*Conj(UP(gt3,0))*ZH(gt2,0)
      + g2u*Conj(UM(gt1,0))*Conj(UP(gt3,1))*ZH(gt2,1));

   return result;
}

std::complex<double> HGTHDMIIMSSMBC_effective_couplings::CphhVWmconjVWm(unsigned gt1) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);

   std::complex<double> result;

   result = 0.5*Sqr(g2)*(v1*ZH(gt1,0) + v2*ZH(gt1,1));

   return result;
}

std::complex<double> HGTHDMIIMSSMBC_effective_couplings::CpAhHmconjHm(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Lambda4 = MODELPARAMETER(Lambda4);
   const auto Lambda5 = MODELPARAMETER(Lambda5);
   const auto Lambda6 = MODELPARAMETER(Lambda6);
   const auto Lambda7 = MODELPARAMETER(Lambda7);
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);

   std::complex<double> result;

   result = std::complex<double>(0,0.5)*(v2*ZA(gt1,0) - v1*ZA(gt1,1))*(ZP(gt2,0
      )*((Lambda6 - Conj(Lambda6))*ZP(gt3,0) + (-Lambda4 + Conj(Lambda5))*ZP(gt3,1
      )) + ZP(gt2,1)*((Lambda4 - Lambda5)*ZP(gt3,0) + (Lambda7 - Conj(Lambda7))*ZP
      (gt3,1)));

   return result;
}

std::complex<double> HGTHDMIIMSSMBC_effective_couplings::CpAhFdbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_775;
   std::complex<double> tmp_776;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_777;
      std::complex<double> tmp_778;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_778 += Conj(Ud(gt3,j1))*Yd(j1,j2);
      }
      tmp_777 += tmp_778;
      tmp_776 += (Conj(Vd(gt2,j2))) * tmp_777;
   }
   tmp_775 += tmp_776;
   result += (std::complex<double>(0.,0.7071067811865475)*ZA(gt1,0)) * tmp_775;

   return result;
}

std::complex<double> HGTHDMIIMSSMBC_effective_couplings::CpAhFubarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_779;
   std::complex<double> tmp_780;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_781;
      std::complex<double> tmp_782;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_782 += Conj(Uu(gt3,j1))*Yu(j1,j2);
      }
      tmp_781 += tmp_782;
      tmp_780 += (Conj(Vu(gt2,j2))) * tmp_781;
   }
   tmp_779 += tmp_780;
   result += (std::complex<double>(0.,0.7071067811865475)*ZA(gt1,1)) * tmp_779;

   return result;
}

std::complex<double> HGTHDMIIMSSMBC_effective_couplings::CpAhFebarFePL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Ye = MODELPARAMETER(Ye);

   std::complex<double> result;

   std::complex<double> tmp_783;
   std::complex<double> tmp_784;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_785;
      std::complex<double> tmp_786;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_786 += Conj(Ue(gt3,j1))*Ye(j1,j2);
      }
      tmp_785 += tmp_786;
      tmp_784 += (Conj(Ve(gt2,j2))) * tmp_785;
   }
   tmp_783 += tmp_784;
   result += (std::complex<double>(0.,0.7071067811865475)*ZA(gt1,0)) * tmp_783;

   return result;
}

std::complex<double> HGTHDMIIMSSMBC_effective_couplings::CpAhChabarChaPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g1d = MODELPARAMETER(g1d);
   const auto g2u = MODELPARAMETER(g2u);

   std::complex<double> result;

   result = std::complex<double>(0.,-0.7071067811865475)*(g1d*Conj(UM(gt2,1))*
      Conj(UP(gt3,0))*ZA(gt1,0) - g2u*Conj(UM(gt2,0))*Conj(UP(gt3,1))*ZA(gt1,1));

   return result;
}

void HGTHDMIIMSSMBC_effective_couplings::calculate_eff_CphhVPVP(unsigned gO1)
{
   const auto MHm = PHYSICAL(MHm);
   const auto MFd = PHYSICAL(MFd);
   const auto MFu = PHYSICAL(MFu);
   const auto MFe = PHYSICAL(MFe);
   const auto MCha = PHYSICAL(MCha);
   const auto MVWm = PHYSICAL(MVWm);
   const auto decay_mass = PHYSICAL(Mhh)(gO1);
   const auto decay_scale = 0.25 * Sqr(decay_mass);
   const auto saved_ZH = ZH;
   ZH = PHYSICAL(ZH);

   const auto vev = 1.0 / Sqrt(qedqcd.displayFermiConstant() * Sqrt(2.0));

   std::complex<double> result = 0;
   for (unsigned gI1 = 1; gI1 < 2; ++gI1) {
      result += 0.5 * CphhHmconjHm(gO1, gI1, gI1) * vev * AS0(decay_scale /
         Sqr(MHm(gI1))) / Sqr(MHm(gI1));
   }
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += 0.3333333333333333 * scalar_fermion_qcd_factor(decay_mass,
         MFd(gI1)) * CpFdhhbarFdPL(gI1, gO1, gI1) * vev * AS12(decay_scale / Sqr(
         MFd(gI1))) / MFd(gI1);
   }
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += 1.3333333333333333 * scalar_fermion_qcd_factor(decay_mass,
         MFu(gI1)) * CpFuhhbarFuPL(gI1, gO1, gI1) * vev * AS12(decay_scale / Sqr(
         MFu(gI1))) / MFu(gI1);
   }
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += CpFehhbarFePL(gI1, gO1, gI1) * vev * AS12(decay_scale / Sqr(
         MFe(gI1))) / MFe(gI1);
   }
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      result += CpChahhbarChaPL(gI1, gO1, gI1) * vev * AS12(decay_scale /
         Sqr(MCha(gI1))) / MCha(gI1);
   }
   result += -0.5 * CphhVWmconjVWm(gO1) * vev * AS1(decay_scale / Sqr(MVWm)) /
      Sqr(MVWm);


   result *= 0.1892681907127351 * physical_input.get(Physical_input::alpha_em_0
      ) * Sqrt(qedqcd.displayFermiConstant());

   ZH = saved_ZH;
   eff_CphhVPVP(gO1) = result;

}

void HGTHDMIIMSSMBC_effective_couplings::calculate_eff_CphhVGVG(unsigned gO1)
{
   const auto g3 = MODELPARAMETER(g3);
   const auto MFd = PHYSICAL(MFd);
   const auto MFu = PHYSICAL(MFu);
   const double alpha_s = 0.07957747154594767*Sqr(g3);
   const auto decay_mass = PHYSICAL(Mhh)(gO1);
   const auto decay_scale = 0.25 * Sqr(decay_mass);
   const auto saved_ZH = ZH;
   ZH = PHYSICAL(ZH);

   const auto vev = 1.0 / Sqrt(qedqcd.displayFermiConstant() * Sqrt(2.0));

   std::complex<double> result = 0;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += CpFdhhbarFdPL(gI1, gO1, gI1) * vev * AS12(decay_scale / Sqr(
         MFd(gI1))) / MFd(gI1);
   }
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += CpFuhhbarFuPL(gI1, gO1, gI1) * vev * AS12(decay_scale / Sqr(
         MFu(gI1))) / MFu(gI1);
   }
   result *= std::complex<double>(0.75,0.);

   if (include_qcd_corrections) {
      result *= scalar_scaling_factor(decay_mass);
   }


   result *= 0.12617879380849006 * alpha_s * Sqrt(qedqcd.displayFermiConstant()
      );

   ZH = saved_ZH;
   eff_CphhVGVG(gO1) = result;

}

void HGTHDMIIMSSMBC_effective_couplings::calculate_eff_CpAhVPVP(unsigned gO1)
{
   const auto MFd = PHYSICAL(MFd);
   const auto MFu = PHYSICAL(MFu);
   const auto MFe = PHYSICAL(MFe);
   const auto MCha = PHYSICAL(MCha);
   const auto decay_mass = PHYSICAL(MAh)(gO1);
   const auto decay_scale = 0.25 * Sqr(decay_mass);
   const auto saved_ZA = ZA;
   ZA = PHYSICAL(ZA);

   const auto vev = 1.0 / Sqrt(qedqcd.displayFermiConstant() * Sqrt(2.0));

   std::complex<double> result = 0;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += 0.3333333333333333 * pseudoscalar_fermion_qcd_factor(
         decay_mass, MFd(gI1)) * CpAhFdbarFdPL(gO1, gI1, gI1) * vev * AP12(
         decay_scale / Sqr(MFd(gI1))) / MFd(gI1);
   }
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += 1.3333333333333333 * pseudoscalar_fermion_qcd_factor(
         decay_mass, MFu(gI1)) * CpAhFubarFuPL(gO1, gI1, gI1) * vev * AP12(
         decay_scale / Sqr(MFu(gI1))) / MFu(gI1);
   }
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += CpAhFebarFePL(gO1, gI1, gI1) * vev * AP12(decay_scale / Sqr(
         MFe(gI1))) / MFe(gI1);
   }
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      result += CpAhChabarChaPL(gO1, gI1, gI1) * vev * AP12(decay_scale /
         Sqr(MCha(gI1))) / MCha(gI1);
   }
   result *= std::complex<double>(2.0,0.);


   result *= 0.1892681907127351 * physical_input.get(Physical_input::alpha_em_0
      ) * Sqrt(qedqcd.displayFermiConstant());

   ZA = saved_ZA;
   eff_CpAhVPVP(gO1) = result;

}

void HGTHDMIIMSSMBC_effective_couplings::calculate_eff_CpAhVGVG(unsigned gO1)
{
   const auto g3 = MODELPARAMETER(g3);
   const auto MFd = PHYSICAL(MFd);
   const auto MFu = PHYSICAL(MFu);
   const double alpha_s = 0.07957747154594767*Sqr(g3);
   const auto decay_mass = PHYSICAL(MAh)(gO1);
   const auto decay_scale = 0.25 * Sqr(decay_mass);
   const auto saved_ZA = ZA;
   ZA = PHYSICAL(ZA);

   const auto vev = 1.0 / Sqrt(qedqcd.displayFermiConstant() * Sqrt(2.0));

   std::complex<double> result = 0;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += CpAhFdbarFdPL(gO1, gI1, gI1) * vev * AP12(decay_scale / Sqr(
         MFd(gI1))) / MFd(gI1);
   }
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += CpAhFubarFuPL(gO1, gI1, gI1) * vev * AP12(decay_scale / Sqr(
         MFu(gI1))) / MFu(gI1);
   }
   result *= std::complex<double>(1.5,0.);

   if (include_qcd_corrections) {
      result *= pseudoscalar_scaling_factor(decay_mass);
   }


   result *= 0.12617879380849006 * alpha_s * Sqrt(qedqcd.displayFermiConstant()
      );

   ZA = saved_ZA;
   eff_CpAhVGVG(gO1) = result;

}


} // namespace flexiblesusy
