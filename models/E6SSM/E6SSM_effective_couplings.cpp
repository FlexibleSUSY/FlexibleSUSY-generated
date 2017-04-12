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

// File generated at Wed 12 Apr 2017 13:02:56

#include "E6SSM_effective_couplings.hpp"

#include "effective_couplings.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

using namespace effective_couplings;

#define INPUTPARAMETER(parameter) model.get_input().parameter
#define MODELPARAMETER(parameter) model.get_##parameter()
#define DERIVEDPARAMETER(parameter) model.##parameter()
#define PHASE(parameter) model.get_##parameter()
#define PHYSICAL(parameter) model.get_physical().parameter

E6SSM_effective_couplings::E6SSM_effective_couplings(
   const E6SSM_mass_eigenstates& model_,
   const softsusy::QedQcd& qedqcd_,
   const Physical_input& input_)
   : model(model_), qedqcd(qedqcd_), physical_input(input_)
   , rg_improve(true), include_qcd_corrections(true)
   , ZD(MODELPARAMETER(ZD)), ZV(MODELPARAMETER(ZV)), ZU(MODELPARAMETER(ZU)), ZE
      (MODELPARAMETER(ZE)), ZDX(MODELPARAMETER(ZDX)), ZH(MODELPARAMETER(ZH)), ZA(
      MODELPARAMETER(ZA)), ZP(MODELPARAMETER(ZP)), ZN(MODELPARAMETER(ZN)), UM(
      MODELPARAMETER(UM)), UP(MODELPARAMETER(UP)), ZEL(MODELPARAMETER(ZEL)), ZER(
      MODELPARAMETER(ZER)), ZDL(MODELPARAMETER(ZDL)), ZDR(MODELPARAMETER(ZDR)),
      ZUL(MODELPARAMETER(ZUL)), ZUR(MODELPARAMETER(ZUR)), ZDXL(MODELPARAMETER(ZDXL
      )), ZDXR(MODELPARAMETER(ZDXR)), UHI0(MODELPARAMETER(UHI0)), UHIp(
      MODELPARAMETER(UHIp)), ZMI(MODELPARAMETER(ZMI)), ZPI(MODELPARAMETER(ZPI)),
      ZNI(MODELPARAMETER(ZNI)), ZSSI(MODELPARAMETER(ZSSI)), ZFSI(MODELPARAMETER(
      ZFSI)), UHp0(MODELPARAMETER(UHp0)), UHpp(MODELPARAMETER(UHpp)), ZNp(
      MODELPARAMETER(ZNp)), ZZ(MODELPARAMETER(ZZ))

   , eff_CphhVPVP(Eigen::Array<std::complex<double>,3,1>::Zero()), eff_CphhVGVG
      (Eigen::Array<std::complex<double>,3,1>::Zero()), eff_CpAhVPVP(Eigen::Array<
      std::complex<double>,3,1>::Zero()), eff_CpAhVGVG(Eigen::Array<std::complex<
      double>,3,1>::Zero())

{
}

void E6SSM_effective_couplings::calculate_effective_couplings()
{
   const standard_model::Standard_model sm(initialise_SM());

   const double scale = model.get_scale();
   const Eigen::ArrayXd saved_parameters(model.get());

   const double saved_mt = PHYSICAL(MFu(2));
   PHYSICAL(MFu(2)) = qedqcd.displayPoleMt();

   const auto Mhh = PHYSICAL(Mhh);
   for (unsigned gO1 = 0; gO1 < 3; ++gO1) {
      if (rg_improve && scale > Mhh(gO1)) {
         model.run_to(Mhh(gO1));
      }
      model.calculate_DRbar_masses();
      copy_mixing_matrices_from_model();
      run_SM_strong_coupling_to(sm, 0.5 * Mhh(gO1));
      calculate_eff_CphhVPVP(gO1);
      run_SM_strong_coupling_to(sm, Mhh(gO1));
      calculate_eff_CphhVGVG(gO1);
   }

   const auto MAh = PHYSICAL(MAh);
   for (unsigned gO1 = 2; gO1 < 3; ++gO1) {
      if (rg_improve && scale > MAh(gO1)) {
         model.run_to(MAh(gO1));
      }
      model.calculate_DRbar_masses();
      copy_mixing_matrices_from_model();
      run_SM_strong_coupling_to(sm, 0.5 * MAh(gO1));
      calculate_eff_CpAhVPVP(gO1);
      run_SM_strong_coupling_to(sm, MAh(gO1));
      calculate_eff_CpAhVGVG(gO1);
   }

   PHYSICAL(MFu(2)) = saved_mt;
   model.set_scale(scale);
   model.set(saved_parameters);

}

void E6SSM_effective_couplings::set_model(const E6SSM_mass_eigenstates& model_)
{
   model = model_;
   copy_mixing_matrices_from_model();
}

void E6SSM_effective_couplings::copy_mixing_matrices_from_model()
{
   ZD = MODELPARAMETER(ZD);
   ZV = MODELPARAMETER(ZV);
   ZU = MODELPARAMETER(ZU);
   ZE = MODELPARAMETER(ZE);
   ZDX = MODELPARAMETER(ZDX);
   ZH = MODELPARAMETER(ZH);
   ZA = MODELPARAMETER(ZA);
   ZP = MODELPARAMETER(ZP);
   ZN = MODELPARAMETER(ZN);
   UM = MODELPARAMETER(UM);
   UP = MODELPARAMETER(UP);
   ZEL = MODELPARAMETER(ZEL);
   ZER = MODELPARAMETER(ZER);
   ZDL = MODELPARAMETER(ZDL);
   ZDR = MODELPARAMETER(ZDR);
   ZUL = MODELPARAMETER(ZUL);
   ZUR = MODELPARAMETER(ZUR);
   ZDXL = MODELPARAMETER(ZDXL);
   ZDXR = MODELPARAMETER(ZDXR);
   UHI0 = MODELPARAMETER(UHI0);
   UHIp = MODELPARAMETER(UHIp);
   ZMI = MODELPARAMETER(ZMI);
   ZPI = MODELPARAMETER(ZPI);
   ZNI = MODELPARAMETER(ZNI);
   ZSSI = MODELPARAMETER(ZSSI);
   ZFSI = MODELPARAMETER(ZFSI);
   UHp0 = MODELPARAMETER(UHp0);
   UHpp = MODELPARAMETER(UHpp);
   ZNp = MODELPARAMETER(ZNp);
   ZZ = MODELPARAMETER(ZZ);

}

standard_model::Standard_model E6SSM_effective_couplings::initialise_SM() const
{
   standard_model::Standard_model sm;

   sm.set_loops(2);
   sm.set_thresholds(2);
   sm.set_physical_input(physical_input);

   sm.initialise_from_input(qedqcd);

   return sm;
}

void E6SSM_effective_couplings::run_SM_strong_coupling_to(standard_model::Standard_model sm, double m)
{
   sm.run_to(m);

   model.set_g3(sm.get_g3());

}

std::complex<double> E6SSM_effective_couplings::scalar_scalar_qcd_factor(double m_decay, double m_loop) const
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

std::complex<double> E6SSM_effective_couplings::scalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) *
         scalar_diphoton_fermion_loop(m_decay, m_loop);

   }

   return result;
}

std::complex<double> E6SSM_effective_couplings::pseudoscalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) *
         pseudoscalar_diphoton_fermion_loop(m_decay, m_loop);

   }

   return result;
}

double E6SSM_effective_couplings::number_of_active_flavours(double m) const
{
   if (m < qedqcd.displayMbMb()) {
      return 4.0;
   } else if (m < qedqcd.displayPoleMt()) {
      return 5.0;
   } else {
      return 6.0;
   }
}

double E6SSM_effective_couplings::scalar_scaling_factor(double m) const
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

double E6SSM_effective_couplings::pseudoscalar_scaling_factor(double m) const
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

double E6SSM_effective_couplings::get_hhVPVP_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(Mhh)(gO1);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CphhVPVP(gO1));
}

double E6SSM_effective_couplings::get_hhVGVG_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(Mhh)(gO1);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CphhVGVG(gO1));
}

double E6SSM_effective_couplings::get_AhVPVP_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(MAh)(gO1);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CpAhVPVP(gO1));
}

double E6SSM_effective_couplings::get_AhVGVG_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(MAh)(gO1);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CpAhVGVG(gO1));
}

std::complex<double> E6SSM_effective_couplings::CphhSdconjSd(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYd = MODELPARAMETER(TYd);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto vd = MODELPARAMETER(vd);
   const auto vs = MODELPARAMETER(vs);
   const auto vu = MODELPARAMETER(vu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   std::complex<double> result;

   std::complex<double> tmp_7690;
   std::complex<double> tmp_7691;
   std::complex<double> tmp_7692;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7692 += Conj(ZD(gt2,j1))*ZD(gt3,j1);
   }
   tmp_7691 += tmp_7692;
   tmp_7690 += (vd*(6*Sqr(g1) + 30*Sqr(g2) + 9*Sqr(gN))*ZH(gt1,0) - 2*vu*(3*Sqr
      (g1) + 15*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,1) - 15*vs*Sqr(gN)*ZH(gt1,2)) *
      tmp_7691;
   std::complex<double> tmp_7693;
   std::complex<double> tmp_7694;
   std::complex<double> tmp_7695;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7695 += Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j1))*ZD(gt3,j1);
   }
   tmp_7694 += tmp_7695;
   tmp_7693 += (-42.42640687119285*ZH(gt1,0)) * tmp_7694;
   std::complex<double> tmp_7696;
   std::complex<double> tmp_7697;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7697 += Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1);
   }
   tmp_7696 += tmp_7697;
   tmp_7693 += (vd*(6*Sqr(g1) + 9*Sqr(gN))*ZH(gt1,0) + (-6*vu*Sqr(g1) + 6*vu*
      Sqr(gN))*ZH(gt1,1) - 15*vs*Sqr(gN)*ZH(gt1,2)) * tmp_7696;
   std::complex<double> tmp_7698;
   std::complex<double> tmp_7699;
   std::complex<double> tmp_7700;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7700 += Conj(ZD(gt2,j1))*ZD(gt3,3 + j1)*TYd(j1,j1);
   }
   tmp_7699 += tmp_7700;
   tmp_7698 += (1.4142135623730951*ZH(gt1,0)) * tmp_7699;
   std::complex<double> tmp_7701;
   std::complex<double> tmp_7702;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_7702 += AbsSqr(Yd(j2,j2))*Conj(ZD(gt2,j2))*ZD(gt3,j2);
   }
   tmp_7701 += tmp_7702;
   tmp_7698 += (2*vd*ZH(gt1,0)) * tmp_7701;
   std::complex<double> tmp_7703;
   std::complex<double> tmp_7704;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_7704 += AbsSqr(Yd(j2,j2))*Conj(ZD(gt2,3 + j2))*ZD(gt3,3 + j2);
   }
   tmp_7703 += tmp_7704;
   tmp_7698 += (2*vd*ZH(gt1,0)) * tmp_7703;
   std::complex<double> tmp_7705;
   std::complex<double> tmp_7706;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7706 += Conj(Yd(j1,j1))*Conj(ZD(gt2,3 + j1))*ZD(gt3,j1);
   }
   tmp_7705 += tmp_7706;
   tmp_7698 += (-(vs*Lambdax*ZH(gt1,1))) * tmp_7705;
   std::complex<double> tmp_7707;
   std::complex<double> tmp_7708;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7708 += Conj(ZD(gt2,j1))*Yd(j1,j1)*ZD(gt3,3 + j1);
   }
   tmp_7707 += tmp_7708;
   tmp_7698 += (-(vs*Conj(Lambdax)*ZH(gt1,1))) * tmp_7707;
   std::complex<double> tmp_7709;
   std::complex<double> tmp_7710;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7710 += Conj(Yd(j1,j1))*Conj(ZD(gt2,3 + j1))*ZD(gt3,j1);
   }
   tmp_7709 += tmp_7710;
   tmp_7698 += (-(vu*Lambdax*ZH(gt1,2))) * tmp_7709;
   std::complex<double> tmp_7711;
   std::complex<double> tmp_7712;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7712 += Conj(ZD(gt2,j1))*Yd(j1,j1)*ZD(gt3,3 + j1);
   }
   tmp_7711 += tmp_7712;
   tmp_7698 += (-(vu*Conj(Lambdax)*ZH(gt1,2))) * tmp_7711;
   tmp_7693 += (-30) * tmp_7698;
   tmp_7690 += (2) * tmp_7693;
   result += (0.008333333333333333) * tmp_7690;

   return result;
}

std::complex<double> E6SSM_effective_couplings::CphhSuconjSu(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYu = MODELPARAMETER(TYu);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto vd = MODELPARAMETER(vd);
   const auto vs = MODELPARAMETER(vs);
   const auto vu = MODELPARAMETER(vu);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   std::complex<double> result;

   std::complex<double> tmp_7713;
   std::complex<double> tmp_7714;
   std::complex<double> tmp_7715;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7715 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_7714 += tmp_7715;
   tmp_7713 += (-24*vd*Sqr(g1)*ZH(gt1,0)) * tmp_7714;
   std::complex<double> tmp_7716;
   std::complex<double> tmp_7717;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7717 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_7716 += tmp_7717;
   tmp_7713 += (9*vd*Sqr(gN)*ZH(gt1,0)) * tmp_7716;
   std::complex<double> tmp_7718;
   std::complex<double> tmp_7719;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7719 += Conj(ZU(gt2,j1))*Yu(j1,j1)*ZU(gt3,3 + j1);
   }
   tmp_7718 += tmp_7719;
   tmp_7713 += (60*vs*Conj(Lambdax)*ZH(gt1,0)) * tmp_7718;
   std::complex<double> tmp_7720;
   std::complex<double> tmp_7721;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7721 += Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j1))*ZU(gt3,j1);
   }
   tmp_7720 += tmp_7721;
   tmp_7713 += (-84.8528137423857*ZH(gt1,1)) * tmp_7720;
   std::complex<double> tmp_7722;
   std::complex<double> tmp_7723;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7723 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_7722 += tmp_7723;
   tmp_7713 += (24*vu*Sqr(g1)*ZH(gt1,1)) * tmp_7722;
   std::complex<double> tmp_7724;
   std::complex<double> tmp_7725;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7725 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_7724 += tmp_7725;
   tmp_7713 += (6*vu*Sqr(gN)*ZH(gt1,1)) * tmp_7724;
   std::complex<double> tmp_7726;
   std::complex<double> tmp_7727;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7727 += Conj(ZU(gt2,j1))*ZU(gt3,3 + j1)*TYu(j1,j1);
   }
   tmp_7726 += tmp_7727;
   tmp_7713 += (-84.8528137423857*ZH(gt1,1)) * tmp_7726;
   std::complex<double> tmp_7728;
   std::complex<double> tmp_7729;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_7729 += AbsSqr(Yu(j2,j2))*Conj(ZU(gt2,j2))*ZU(gt3,j2);
   }
   tmp_7728 += tmp_7729;
   tmp_7713 += (-120*vu*ZH(gt1,1)) * tmp_7728;
   std::complex<double> tmp_7730;
   std::complex<double> tmp_7731;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_7731 += AbsSqr(Yu(j2,j2))*Conj(ZU(gt2,3 + j2))*ZU(gt3,3 + j2);
   }
   tmp_7730 += tmp_7731;
   tmp_7713 += (-120*vu*ZH(gt1,1)) * tmp_7730;
   std::complex<double> tmp_7732;
   std::complex<double> tmp_7733;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7733 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_7732 += tmp_7733;
   tmp_7713 += (-15*vs*Sqr(gN)*ZH(gt1,2)) * tmp_7732;
   std::complex<double> tmp_7734;
   std::complex<double> tmp_7735;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7735 += Conj(ZU(gt2,j1))*Yu(j1,j1)*ZU(gt3,3 + j1);
   }
   tmp_7734 += tmp_7735;
   tmp_7713 += (60*vd*Conj(Lambdax)*ZH(gt1,2)) * tmp_7734;
   std::complex<double> tmp_7736;
   std::complex<double> tmp_7737;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7737 += Conj(Yu(j1,j1))*Conj(ZU(gt2,3 + j1))*ZU(gt3,j1);
   }
   tmp_7736 += tmp_7737;
   tmp_7713 += (60*Lambdax*(vs*ZH(gt1,0) + vd*ZH(gt1,2))) * tmp_7736;
   std::complex<double> tmp_7738;
   std::complex<double> tmp_7739;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7739 += Conj(ZU(gt2,j1))*ZU(gt3,j1);
   }
   tmp_7738 += tmp_7739;
   tmp_7713 += (vd*(6*Sqr(g1) - 30*Sqr(g2) + 9*Sqr(gN))*ZH(gt1,0) + 2*vu*(-3*
      Sqr(g1) + 15*Sqr(g2) + 3*Sqr(gN))*ZH(gt1,1) - 15*vs*Sqr(gN)*ZH(gt1,2)) *
      tmp_7738;
   result += (0.008333333333333333) * tmp_7713;

   return result;
}

std::complex<double> E6SSM_effective_couplings::CphhSeconjSe(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYe = MODELPARAMETER(TYe);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto vd = MODELPARAMETER(vd);
   const auto vs = MODELPARAMETER(vs);
   const auto vu = MODELPARAMETER(vu);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   std::complex<double> result;

   std::complex<double> tmp_7740;
   std::complex<double> tmp_7741;
   std::complex<double> tmp_7742;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7742 += Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j1))*ZE(gt3,j1);
   }
   tmp_7741 += tmp_7742;
   tmp_7740 += (-28.284271247461902*ZH(gt1,0)) * tmp_7741;
   std::complex<double> tmp_7743;
   std::complex<double> tmp_7744;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7744 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_7743 += tmp_7744;
   tmp_7740 += (12*vd*Sqr(g1)*ZH(gt1,0)) * tmp_7743;
   std::complex<double> tmp_7745;
   std::complex<double> tmp_7746;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7746 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_7745 += tmp_7746;
   tmp_7740 += (3*vd*Sqr(gN)*ZH(gt1,0)) * tmp_7745;
   std::complex<double> tmp_7747;
   std::complex<double> tmp_7748;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7748 += Conj(ZE(gt2,j1))*ZE(gt3,3 + j1)*TYe(j1,j1);
   }
   tmp_7747 += tmp_7748;
   tmp_7740 += (-28.284271247461902*ZH(gt1,0)) * tmp_7747;
   std::complex<double> tmp_7749;
   std::complex<double> tmp_7750;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_7750 += AbsSqr(Ye(j2,j2))*Conj(ZE(gt2,j2))*ZE(gt3,j2);
   }
   tmp_7749 += tmp_7750;
   tmp_7740 += (-40*vd*ZH(gt1,0)) * tmp_7749;
   std::complex<double> tmp_7751;
   std::complex<double> tmp_7752;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_7752 += AbsSqr(Ye(j2,j2))*Conj(ZE(gt2,3 + j2))*ZE(gt3,3 + j2);
   }
   tmp_7751 += tmp_7752;
   tmp_7740 += (-40*vd*ZH(gt1,0)) * tmp_7751;
   std::complex<double> tmp_7753;
   std::complex<double> tmp_7754;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7754 += Conj(Ye(j1,j1))*Conj(ZE(gt2,3 + j1))*ZE(gt3,j1);
   }
   tmp_7753 += tmp_7754;
   tmp_7740 += (20*vs*Lambdax*ZH(gt1,1)) * tmp_7753;
   std::complex<double> tmp_7755;
   std::complex<double> tmp_7756;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7756 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_7755 += tmp_7756;
   tmp_7740 += (-12*vu*Sqr(g1)*ZH(gt1,1)) * tmp_7755;
   std::complex<double> tmp_7757;
   std::complex<double> tmp_7758;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7758 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_7757 += tmp_7758;
   tmp_7740 += (2*vu*Sqr(gN)*ZH(gt1,1)) * tmp_7757;
   std::complex<double> tmp_7759;
   std::complex<double> tmp_7760;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7760 += Conj(ZE(gt2,j1))*Ye(j1,j1)*ZE(gt3,3 + j1);
   }
   tmp_7759 += tmp_7760;
   tmp_7740 += (20*vs*Conj(Lambdax)*ZH(gt1,1)) * tmp_7759;
   std::complex<double> tmp_7761;
   std::complex<double> tmp_7762;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7762 += Conj(Ye(j1,j1))*Conj(ZE(gt2,3 + j1))*ZE(gt3,j1);
   }
   tmp_7761 += tmp_7762;
   tmp_7740 += (20*vu*Lambdax*ZH(gt1,2)) * tmp_7761;
   std::complex<double> tmp_7763;
   std::complex<double> tmp_7764;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7764 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_7763 += tmp_7764;
   tmp_7740 += (-5*vs*Sqr(gN)*ZH(gt1,2)) * tmp_7763;
   std::complex<double> tmp_7765;
   std::complex<double> tmp_7766;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7766 += Conj(ZE(gt2,j1))*Ye(j1,j1)*ZE(gt3,3 + j1);
   }
   tmp_7765 += tmp_7766;
   tmp_7740 += (20*vu*Conj(Lambdax)*ZH(gt1,2)) * tmp_7765;
   std::complex<double> tmp_7767;
   std::complex<double> tmp_7768;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7768 += Conj(ZE(gt2,j1))*ZE(gt3,j1);
   }
   tmp_7767 += tmp_7768;
   tmp_7740 += (-2*(vd*(3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,0) + vu*(-3*
      Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN))*ZH(gt1,1) + 5*vs*Sqr(gN)*ZH(gt1,2))) *
      tmp_7767;
   result += (0.025) * tmp_7740;

   return result;
}

std::complex<double> E6SSM_effective_couplings::CphhSDXconjSDX(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TKappa = MODELPARAMETER(TKappa);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto vd = MODELPARAMETER(vd);
   const auto vs = MODELPARAMETER(vs);
   const auto vu = MODELPARAMETER(vu);
   const auto Kappa = MODELPARAMETER(Kappa);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   std::complex<double> result;

   std::complex<double> tmp_7769;
   std::complex<double> tmp_7770;
   std::complex<double> tmp_7771;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7771 += Conj(ZDX(gt2,3 + j1))*ZDX(gt3,3 + j1);
   }
   tmp_7770 += tmp_7771;
   tmp_7769 += (12*vd*Sqr(g1)*ZH(gt1,0)) * tmp_7770;
   std::complex<double> tmp_7772;
   std::complex<double> tmp_7773;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7773 += Conj(ZDX(gt2,3 + j1))*ZDX(gt3,3 + j1);
   }
   tmp_7772 += tmp_7773;
   tmp_7769 += (-27*vd*Sqr(gN)*ZH(gt1,0)) * tmp_7772;
   std::complex<double> tmp_7774;
   std::complex<double> tmp_7775;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7775 += Conj(ZDX(gt2,j1))*ZDX(gt3,3 + j1)*Kappa(j1,j1);
   }
   tmp_7774 += tmp_7775;
   tmp_7769 += (60*vu*Conj(Lambdax)*ZH(gt1,0)) * tmp_7774;
   std::complex<double> tmp_7776;
   std::complex<double> tmp_7777;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7777 += Conj(ZDX(gt2,3 + j1))*ZDX(gt3,3 + j1);
   }
   tmp_7776 += tmp_7777;
   tmp_7769 += (-12*vu*Sqr(g1)*ZH(gt1,1)) * tmp_7776;
   std::complex<double> tmp_7778;
   std::complex<double> tmp_7779;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7779 += Conj(ZDX(gt2,3 + j1))*ZDX(gt3,3 + j1);
   }
   tmp_7778 += tmp_7779;
   tmp_7769 += (-18*vu*Sqr(gN)*ZH(gt1,1)) * tmp_7778;
   std::complex<double> tmp_7780;
   std::complex<double> tmp_7781;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7781 += Conj(ZDX(gt2,j1))*ZDX(gt3,3 + j1)*Kappa(j1,j1);
   }
   tmp_7780 += tmp_7781;
   tmp_7769 += (60*vd*Conj(Lambdax)*ZH(gt1,1)) * tmp_7780;
   std::complex<double> tmp_7782;
   std::complex<double> tmp_7783;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7783 += Conj(ZDX(gt2,3 + j1))*Conj(Kappa(j1,j1))*ZDX(gt3,j1);
   }
   tmp_7782 += tmp_7783;
   tmp_7769 += (60*Lambdax*(vu*ZH(gt1,0) + vd*ZH(gt1,1))) * tmp_7782;
   std::complex<double> tmp_7784;
   std::complex<double> tmp_7785;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7785 += Conj(ZDX(gt2,3 + j1))*Conj(TKappa(j1,j1))*ZDX(gt3,j1);
   }
   tmp_7784 += tmp_7785;
   tmp_7769 += (-84.8528137423857*ZH(gt1,2)) * tmp_7784;
   std::complex<double> tmp_7786;
   std::complex<double> tmp_7787;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7787 += Conj(ZDX(gt2,3 + j1))*ZDX(gt3,3 + j1);
   }
   tmp_7786 += tmp_7787;
   tmp_7769 += (45*vs*Sqr(gN)*ZH(gt1,2)) * tmp_7786;
   std::complex<double> tmp_7788;
   std::complex<double> tmp_7789;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7789 += Conj(ZDX(gt2,j1))*ZDX(gt3,3 + j1)*TKappa(j1,j1);
   }
   tmp_7788 += tmp_7789;
   tmp_7769 += (-84.8528137423857*ZH(gt1,2)) * tmp_7788;
   std::complex<double> tmp_7790;
   std::complex<double> tmp_7791;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_7791 += AbsSqr(Kappa(j2,j2))*Conj(ZDX(gt2,j2))*ZDX(gt3,j2);
   }
   tmp_7790 += tmp_7791;
   tmp_7769 += (-120*vs*ZH(gt1,2)) * tmp_7790;
   std::complex<double> tmp_7792;
   std::complex<double> tmp_7793;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_7793 += AbsSqr(Kappa(j2,j2))*Conj(ZDX(gt2,3 + j2))*ZDX(gt3,3 + j2)
         ;
   }
   tmp_7792 += tmp_7793;
   tmp_7769 += (-120*vs*ZH(gt1,2)) * tmp_7792;
   std::complex<double> tmp_7794;
   std::complex<double> tmp_7795;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7795 += Conj(ZDX(gt2,j1))*ZDX(gt3,j1);
   }
   tmp_7794 += tmp_7795;
   tmp_7769 += (-2*(vd*(6*Sqr(g1) + 9*Sqr(gN))*ZH(gt1,0) + (-6*vu*Sqr(g1) + 6*
      vu*Sqr(gN))*ZH(gt1,1) - 15*vs*Sqr(gN)*ZH(gt1,2))) * tmp_7794;
   result += (0.008333333333333333) * tmp_7769;

   return result;
}

std::complex<double> E6SSM_effective_couplings::CphhHpmconjHpm(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TLambdax = MODELPARAMETER(TLambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto vd = MODELPARAMETER(vd);
   const auto vs = MODELPARAMETER(vs);
   const auto vu = MODELPARAMETER(vu);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   std::complex<double> result;

   result = 0.025*(2*ZH(gt1,1)*(ZP(gt2,0)*(vu*(3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN
      ))*ZP(gt3,0) - 5*vd*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gt3,1)) - ZP(gt2,1)*(5
      *vd*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gt3,0) + vu*(3*Sqr(g1) + 5*Sqr(g2) + 2
      *Sqr(gN))*ZP(gt3,1))) - ZH(gt1,0)*(ZP(gt2,0)*(vd*(6*Sqr(g1) + 10*Sqr(g2) + 9
      *Sqr(gN))*ZP(gt3,0) + 10*vu*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gt3,1)) + 2*ZP
      (gt2,1)*(5*vu*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gt3,0) + vd*(-3*Sqr(g1) + 5*
      Sqr(g2) + 3*Sqr(gN))*ZP(gt3,1))) + 5*ZH(gt1,2)*(2*ZP(gt2,1)*(
      -2.8284271247461903*Conj(TLambdax)*ZP(gt3,0) + vs*(-4*AbsSqr(Lambdax) + Sqr(
      gN))*ZP(gt3,1)) + ZP(gt2,0)*(vs*(-8*AbsSqr(Lambdax) + 3*Sqr(gN))*ZP(gt3,0) -
      5.656854249492381*TLambdax*ZP(gt3,1))));

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpChahhbarChaPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   std::complex<double> result;

   result = -0.7071067811865475*(g2*Conj(UM(gt1,0))*Conj(UP(gt3,1))*ZH(gt2,1) +
      Conj(UM(gt1,1))*(g2*Conj(UP(gt3,0))*ZH(gt2,0) + Conj(UP(gt3,1))*Lambdax*ZH(
      gt2,2)));

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpFehhbarFePL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Ye = MODELPARAMETER(Ye);

   std::complex<double> result;

   std::complex<double> tmp_7796;
   std::complex<double> tmp_7797;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7797 += Conj(ZEL(gt1,j1))*Conj(ZER(gt3,j1))*Ye(j1,j1);
   }
   tmp_7796 += tmp_7797;
   result += (-0.7071067811865475*ZH(gt2,0)) * tmp_7796;

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpFdhhbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_7798;
   std::complex<double> tmp_7799;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7799 += Conj(ZDL(gt1,j1))*Conj(ZDR(gt3,j1))*Yd(j1,j1);
   }
   tmp_7798 += tmp_7799;
   result += (-0.7071067811865475*ZH(gt2,0)) * tmp_7798;

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpFuhhbarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_7800;
   std::complex<double> tmp_7801;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7801 += Conj(ZUL(gt1,j1))*Conj(ZUR(gt3,j1))*Yu(j1,j1);
   }
   tmp_7800 += tmp_7801;
   result += (-0.7071067811865475*ZH(gt2,1)) * tmp_7800;

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpFDXhhbarFDXPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Kappa = MODELPARAMETER(Kappa);

   std::complex<double> result;

   std::complex<double> tmp_7802;
   std::complex<double> tmp_7803;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7803 += Conj(ZDXL(gt1,j1))*Conj(ZDXR(gt3,j1))*Kappa(j1,j1);
   }
   tmp_7802 += tmp_7803;
   result += (-0.7071067811865475*ZH(gt2,2)) * tmp_7802;

   return result;
}

std::complex<double> E6SSM_effective_couplings::CphhSHIpconjSHIp(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TLambda12 = MODELPARAMETER(TLambda12);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto vd = MODELPARAMETER(vd);
   const auto vs = MODELPARAMETER(vs);
   const auto vu = MODELPARAMETER(vu);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Lambda12 = MODELPARAMETER(Lambda12);

   std::complex<double> result;

   std::complex<double> tmp_7804;
   std::complex<double> tmp_7805;
   std::complex<double> tmp_7806;
   for (unsigned j1 = 0; j1 < 2; ++j1) {
      tmp_7806 += Conj(UHIp(gt2,j1))*UHIp(gt3,j1);
   }
   tmp_7805 += tmp_7806;
   tmp_7804 += (vd*(-6*Sqr(g1) + 10*Sqr(g2) - 9*Sqr(gN))*ZH(gt1,0) + 2*vu*(3*
      Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,1) + 15*vs*Sqr(gN)*ZH(gt1,2)) *
      tmp_7805;
   std::complex<double> tmp_7807;
   std::complex<double> tmp_7808;
   std::complex<double> tmp_7809;
   for (unsigned j1 = 0; j1 < 2; ++j1) {
      tmp_7809 += Conj(UHIp(gt2,2 + j1))*Conj(Lambda12(j1,j1))*UHIp(gt3,j1);
   }
   tmp_7808 += tmp_7809;
   tmp_7807 += (10*Lambdax*(vu*ZH(gt1,0) + vd*ZH(gt1,1))) * tmp_7808;
   std::complex<double> tmp_7810;
   std::complex<double> tmp_7811;
   for (unsigned j1 = 0; j1 < 2; ++j1) {
      tmp_7811 += Conj(UHIp(gt2,2 + j1))*UHIp(gt3,2 + j1);
   }
   tmp_7810 += tmp_7811;
   tmp_7807 += (vd*(3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,0) + vu*(-3*Sqr(
      g1) + 5*Sqr(g2) - 2*Sqr(gN))*ZH(gt1,1) + 5*vs*Sqr(gN)*ZH(gt1,2)) * tmp_7810;
   std::complex<double> tmp_7812;
   std::complex<double> tmp_7813;
   std::complex<double> tmp_7814;
   for (unsigned j1 = 0; j1 < 2; ++j1) {
      tmp_7814 += Conj(UHIp(gt2,j1))*UHIp(gt3,2 + j1)*Lambda12(j1,j1);
   }
   tmp_7813 += tmp_7814;
   tmp_7812 += (Conj(Lambdax)*(vu*ZH(gt1,0) + vd*ZH(gt1,1))) * tmp_7813;
   std::complex<double> tmp_7815;
   std::complex<double> tmp_7816;
   std::complex<double> tmp_7817;
   for (unsigned j1 = 0; j1 < 2; ++j1) {
      tmp_7817 += Conj(UHIp(gt2,2 + j1))*Conj(TLambda12(j1,j1))*UHIp(gt3,j1)
         ;
   }
   tmp_7816 += tmp_7817;
   tmp_7815 += (1.4142135623730951) * tmp_7816;
   std::complex<double> tmp_7818;
   std::complex<double> tmp_7819;
   for (unsigned j1 = 0; j1 < 2; ++j1) {
      tmp_7819 += Conj(UHIp(gt2,j1))*UHIp(gt3,2 + j1)*TLambda12(j1,j1);
   }
   tmp_7818 += tmp_7819;
   tmp_7815 += (1.4142135623730951) * tmp_7818;
   std::complex<double> tmp_7820;
   std::complex<double> tmp_7821;
   for (unsigned j2 = 0; j2 < 2; ++j2) {
      tmp_7821 += AbsSqr(Lambda12(j2,j2))*Conj(UHIp(gt2,j2))*UHIp(gt3,j2);
   }
   tmp_7820 += tmp_7821;
   std::complex<double> tmp_7822;
   for (unsigned j2 = 0; j2 < 2; ++j2) {
      tmp_7822 += AbsSqr(Lambda12(j2,j2))*Conj(UHIp(gt2,2 + j2))*UHIp(gt3,2
         + j2);
   }
   tmp_7820 += tmp_7822;
   tmp_7815 += (2*vs) * tmp_7820;
   tmp_7812 += (-ZH(gt1,2)) * tmp_7815;
   tmp_7807 += (10) * tmp_7812;
   tmp_7804 += (2) * tmp_7807;
   result += (0.025) * tmp_7804;

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpChaIhhbarChaIPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Lambda12 = MODELPARAMETER(Lambda12);

   std::complex<double> result;

   std::complex<double> tmp_7823;
   std::complex<double> tmp_7824;
   for (unsigned j1 = 0; j1 < 2; ++j1) {
      tmp_7824 += Conj(ZMI(gt1,j1))*Conj(ZPI(gt3,j1))*Lambda12(j1,j1);
   }
   tmp_7823 += tmp_7824;
   result += (-0.7071067811865475*ZH(gt2,2)) * tmp_7823;

   return result;
}

std::complex<double> E6SSM_effective_couplings::CphhSHppconjSHpp(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto vd = MODELPARAMETER(vd);
   const auto vs = MODELPARAMETER(vs);
   const auto vu = MODELPARAMETER(vu);

   std::complex<double> result;

   result = -0.05*(Conj(UHpp(gt2,0))*UHpp(gt3,0) - Conj(UHpp(gt2,1))*UHpp(gt3,1
      ))*(vd*(3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,0) + vu*(-3*Sqr(g1) + 5*
      Sqr(g2) - 2*Sqr(gN))*ZH(gt1,1) + 5*vs*Sqr(gN)*ZH(gt1,2));

   return result;
}

std::complex<double> E6SSM_effective_couplings::CphhVWmconjVWm(unsigned gt1) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);

   std::complex<double> result;

   result = 0.5*Sqr(g2)*(vd*ZH(gt1,0) + vu*ZH(gt1,1));

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpAhSdconjSd(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYd = MODELPARAMETER(TYd);
   const auto vs = MODELPARAMETER(vs);
   const auto vu = MODELPARAMETER(vu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   std::complex<double> result;

   std::complex<double> tmp_7825;
   std::complex<double> tmp_7826;
   std::complex<double> tmp_7827;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7827 += Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j1))*ZD(gt3,j1);
   }
   tmp_7826 += tmp_7827;
   tmp_7825 += (-1.4142135623730951*ZA(gt1,0)) * tmp_7826;
   std::complex<double> tmp_7828;
   std::complex<double> tmp_7829;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7829 += Conj(ZD(gt2,j1))*ZD(gt3,3 + j1)*TYd(j1,j1);
   }
   tmp_7828 += tmp_7829;
   tmp_7825 += (1.4142135623730951*ZA(gt1,0)) * tmp_7828;
   std::complex<double> tmp_7830;
   std::complex<double> tmp_7831;
   std::complex<double> tmp_7832;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7832 += Conj(Yd(j1,j1))*Conj(ZD(gt2,3 + j1))*ZD(gt3,j1);
   }
   tmp_7831 += tmp_7832;
   tmp_7830 += (Lambdax) * tmp_7831;
   std::complex<double> tmp_7833;
   std::complex<double> tmp_7834;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7834 += Conj(ZD(gt2,j1))*Yd(j1,j1)*ZD(gt3,3 + j1);
   }
   tmp_7833 += tmp_7834;
   tmp_7830 += (-Conj(Lambdax)) * tmp_7833;
   tmp_7825 += (-(vs*ZA(gt1,1)) - vu*ZA(gt1,2)) * tmp_7830;
   result += (std::complex<double>(0,-0.5)) * tmp_7825;

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpAhSuconjSu(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYu = MODELPARAMETER(TYu);
   const auto vd = MODELPARAMETER(vd);
   const auto vs = MODELPARAMETER(vs);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   std::complex<double> result;

   std::complex<double> tmp_7835;
   std::complex<double> tmp_7836;
   std::complex<double> tmp_7837;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7837 += Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j1))*ZU(gt3,j1);
   }
   tmp_7836 += tmp_7837;
   std::complex<double> tmp_7838;
   std::complex<double> tmp_7839;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7839 += Conj(ZU(gt2,j1))*ZU(gt3,3 + j1)*TYu(j1,j1);
   }
   tmp_7838 += tmp_7839;
   tmp_7836 += (-1) * tmp_7838;
   tmp_7835 += (-1.4142135623730951*ZA(gt1,1)) * tmp_7836;
   std::complex<double> tmp_7840;
   std::complex<double> tmp_7841;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7841 += Conj(Yu(j1,j1))*Conj(ZU(gt2,3 + j1))*ZU(gt3,j1);
   }
   tmp_7840 += tmp_7841;
   tmp_7835 += (-(Lambdax*(vs*ZA(gt1,0) + vd*ZA(gt1,2)))) * tmp_7840;
   std::complex<double> tmp_7842;
   std::complex<double> tmp_7843;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7843 += Conj(ZU(gt2,j1))*Yu(j1,j1)*ZU(gt3,3 + j1);
   }
   tmp_7842 += tmp_7843;
   tmp_7835 += (Conj(Lambdax)*(vs*ZA(gt1,0) + vd*ZA(gt1,2))) * tmp_7842;
   result += (std::complex<double>(0,-0.5)) * tmp_7835;

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpAhSeconjSe(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYe = MODELPARAMETER(TYe);
   const auto vs = MODELPARAMETER(vs);
   const auto vu = MODELPARAMETER(vu);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   std::complex<double> result;

   std::complex<double> tmp_7844;
   std::complex<double> tmp_7845;
   std::complex<double> tmp_7846;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7846 += Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j1))*ZE(gt3,j1);
   }
   tmp_7845 += tmp_7846;
   tmp_7844 += (-1.4142135623730951*ZA(gt1,0)) * tmp_7845;
   std::complex<double> tmp_7847;
   std::complex<double> tmp_7848;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7848 += Conj(ZE(gt2,j1))*ZE(gt3,3 + j1)*TYe(j1,j1);
   }
   tmp_7847 += tmp_7848;
   tmp_7844 += (1.4142135623730951*ZA(gt1,0)) * tmp_7847;
   std::complex<double> tmp_7849;
   std::complex<double> tmp_7850;
   std::complex<double> tmp_7851;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7851 += Conj(Ye(j1,j1))*Conj(ZE(gt2,3 + j1))*ZE(gt3,j1);
   }
   tmp_7850 += tmp_7851;
   tmp_7849 += (Lambdax) * tmp_7850;
   std::complex<double> tmp_7852;
   std::complex<double> tmp_7853;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7853 += Conj(ZE(gt2,j1))*Ye(j1,j1)*ZE(gt3,3 + j1);
   }
   tmp_7852 += tmp_7853;
   tmp_7849 += (-Conj(Lambdax)) * tmp_7852;
   tmp_7844 += (-(vs*ZA(gt1,1)) - vu*ZA(gt1,2)) * tmp_7849;
   result += (std::complex<double>(0,-0.5)) * tmp_7844;

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpAhSDXconjSDX(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TKappa = MODELPARAMETER(TKappa);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Kappa = MODELPARAMETER(Kappa);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   std::complex<double> result;

   std::complex<double> tmp_7854;
   std::complex<double> tmp_7855;
   std::complex<double> tmp_7856;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7856 += Conj(ZDX(gt2,3 + j1))*Conj(Kappa(j1,j1))*ZDX(gt3,j1);
   }
   tmp_7855 += tmp_7856;
   tmp_7854 += (-(Lambdax*(vu*ZA(gt1,0) + vd*ZA(gt1,1)))) * tmp_7855;
   std::complex<double> tmp_7857;
   std::complex<double> tmp_7858;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7858 += Conj(ZDX(gt2,j1))*ZDX(gt3,3 + j1)*Kappa(j1,j1);
   }
   tmp_7857 += tmp_7858;
   tmp_7854 += (Conj(Lambdax)*(vu*ZA(gt1,0) + vd*ZA(gt1,1))) * tmp_7857;
   std::complex<double> tmp_7859;
   std::complex<double> tmp_7860;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7860 += Conj(ZDX(gt2,3 + j1))*Conj(TKappa(j1,j1))*ZDX(gt3,j1);
   }
   tmp_7859 += tmp_7860;
   std::complex<double> tmp_7861;
   std::complex<double> tmp_7862;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7862 += Conj(ZDX(gt2,j1))*ZDX(gt3,3 + j1)*TKappa(j1,j1);
   }
   tmp_7861 += tmp_7862;
   tmp_7859 += (-1) * tmp_7861;
   tmp_7854 += (-1.4142135623730951*ZA(gt1,2)) * tmp_7859;
   result += (std::complex<double>(0,-0.5)) * tmp_7854;

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpAhHpmconjHpm(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TLambdax = MODELPARAMETER(TLambdax);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   std::complex<double> result;

   result = std::complex<double>(0,-0.25)*(vu*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZA
      (gt1,0)*(ZP(gt2,1)*ZP(gt3,0) - ZP(gt2,0)*ZP(gt3,1)) + vd*(-2*AbsSqr(Lambdax)
      + Sqr(g2))*ZA(gt1,1)*(ZP(gt2,1)*ZP(gt3,0) - ZP(gt2,0)*ZP(gt3,1)) +
      2.8284271247461903*ZA(gt1,2)*(-(Conj(TLambdax)*ZP(gt2,1)*ZP(gt3,0)) +
      TLambdax*ZP(gt2,0)*ZP(gt3,1)));

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpAhChabarChaPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   std::complex<double> result;

   result = std::complex<double>(0.,0.7071067811865475)*(g2*Conj(UM(gt2,0))*
      Conj(UP(gt3,1))*ZA(gt1,1) + Conj(UM(gt2,1))*(g2*Conj(UP(gt3,0))*ZA(gt1,0) -
      Conj(UP(gt3,1))*Lambdax*ZA(gt1,2)));

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpAhFebarFePL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Ye = MODELPARAMETER(Ye);

   std::complex<double> result;

   std::complex<double> tmp_7863;
   std::complex<double> tmp_7864;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7864 += Conj(ZEL(gt2,j1))*Conj(ZER(gt3,j1))*Ye(j1,j1);
   }
   tmp_7863 += tmp_7864;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,0)) *
      tmp_7863;

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpAhFdbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_7865;
   std::complex<double> tmp_7866;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7866 += Conj(ZDL(gt2,j1))*Conj(ZDR(gt3,j1))*Yd(j1,j1);
   }
   tmp_7865 += tmp_7866;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,0)) *
      tmp_7865;

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpAhFubarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_7867;
   std::complex<double> tmp_7868;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7868 += Conj(ZUL(gt2,j1))*Conj(ZUR(gt3,j1))*Yu(j1,j1);
   }
   tmp_7867 += tmp_7868;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,1)) *
      tmp_7867;

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpAhFDXbarFDXPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Kappa = MODELPARAMETER(Kappa);

   std::complex<double> result;

   std::complex<double> tmp_7869;
   std::complex<double> tmp_7870;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7870 += Conj(ZDXL(gt2,j1))*Conj(ZDXR(gt3,j1))*Kappa(j1,j1);
   }
   tmp_7869 += tmp_7870;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,2)) *
      tmp_7869;

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpAhSHIpconjSHIp(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TLambda12 = MODELPARAMETER(TLambda12);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Lambda12 = MODELPARAMETER(Lambda12);

   std::complex<double> result;

   std::complex<double> tmp_7871;
   std::complex<double> tmp_7872;
   std::complex<double> tmp_7873;
   for (unsigned j1 = 0; j1 < 2; ++j1) {
      tmp_7873 += Conj(UHIp(gt2,2 + j1))*Conj(Lambda12(j1,j1))*UHIp(gt3,j1);
   }
   tmp_7872 += tmp_7873;
   tmp_7871 += (-(Lambdax*(vu*ZA(gt1,0) + vd*ZA(gt1,1)))) * tmp_7872;
   std::complex<double> tmp_7874;
   std::complex<double> tmp_7875;
   for (unsigned j1 = 0; j1 < 2; ++j1) {
      tmp_7875 += Conj(UHIp(gt2,j1))*UHIp(gt3,2 + j1)*Lambda12(j1,j1);
   }
   tmp_7874 += tmp_7875;
   tmp_7871 += (Conj(Lambdax)*(vu*ZA(gt1,0) + vd*ZA(gt1,1))) * tmp_7874;
   std::complex<double> tmp_7876;
   std::complex<double> tmp_7877;
   for (unsigned j1 = 0; j1 < 2; ++j1) {
      tmp_7877 += Conj(UHIp(gt2,2 + j1))*Conj(TLambda12(j1,j1))*UHIp(gt3,j1)
         ;
   }
   tmp_7876 += tmp_7877;
   std::complex<double> tmp_7878;
   std::complex<double> tmp_7879;
   for (unsigned j1 = 0; j1 < 2; ++j1) {
      tmp_7879 += Conj(UHIp(gt2,j1))*UHIp(gt3,2 + j1)*TLambda12(j1,j1);
   }
   tmp_7878 += tmp_7879;
   tmp_7876 += (-1) * tmp_7878;
   tmp_7871 += (-1.4142135623730951*ZA(gt1,2)) * tmp_7876;
   result += (std::complex<double>(0,-0.5)) * tmp_7871;

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpAhChaIbarChaIPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Lambda12 = MODELPARAMETER(Lambda12);

   std::complex<double> result;

   std::complex<double> tmp_7880;
   std::complex<double> tmp_7881;
   for (unsigned j1 = 0; j1 < 2; ++j1) {
      tmp_7881 += Conj(ZMI(gt2,j1))*Conj(ZPI(gt3,j1))*Lambda12(j1,j1);
   }
   tmp_7880 += tmp_7881;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,2)) *
      tmp_7880;

   return result;
}

void E6SSM_effective_couplings::calculate_eff_CphhVPVP(unsigned gO1)
{
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto MSd = MODELPARAMETER(MSd);
   const auto MSu = MODELPARAMETER(MSu);
   const auto MSe = MODELPARAMETER(MSe);
   const auto MSDX = MODELPARAMETER(MSDX);
   const auto MHpm = MODELPARAMETER(MHpm);
   const auto MCha = MODELPARAMETER(MCha);
   const auto MFe = MODELPARAMETER(MFe);
   const auto MFd = MODELPARAMETER(MFd);
   const auto MFu = MODELPARAMETER(MFu);
   const auto MFDX = MODELPARAMETER(MFDX);
   const auto MSHIp = MODELPARAMETER(MSHIp);
   const auto MChaI = MODELPARAMETER(MChaI);
   const auto MSHpp = MODELPARAMETER(MSHpp);
   const auto MVWm = MODELPARAMETER(MVWm);
   const auto decay_mass = PHYSICAL(Mhh)(gO1);
   const auto decay_scale = 0.25 * Sqr(decay_mass);
   const auto saved_ZH = ZH;
   ZH = PHYSICAL(ZH);

   const auto vev = Sqrt(Sqr(vd) + Sqr(vu));

   std::complex<double> result = 0;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      result += 0.16666666666666666 * scalar_scalar_qcd_factor(decay_mass,
         MSd(gI1)) * CphhSdconjSd(gO1, gI1, gI1) * vev * AS0(decay_scale / Sqr(MSd
         (gI1))) / Sqr(MSd(gI1));
   }
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      result += 0.6666666666666666 * scalar_scalar_qcd_factor(decay_mass,
         MSu(gI1)) * CphhSuconjSu(gO1, gI1, gI1) * vev * AS0(decay_scale / Sqr(MSu
         (gI1))) / Sqr(MSu(gI1));
   }
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      result += 0.5 * CphhSeconjSe(gO1, gI1, gI1) * vev * AS0(decay_scale /
         Sqr(MSe(gI1))) / Sqr(MSe(gI1));
   }
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      result += 0.16666666666666666 * scalar_scalar_qcd_factor(decay_mass,
         MSDX(gI1)) * CphhSDXconjSDX(gO1, gI1, gI1) * vev * AS0(decay_scale / Sqr(
         MSDX(gI1))) / Sqr(MSDX(gI1));
   }
   for (unsigned gI1 = 1; gI1 < 2; ++gI1) {
      result += 0.5 * CphhHpmconjHpm(gO1, gI1, gI1) * vev * AS0(decay_scale
         / Sqr(MHpm(gI1))) / Sqr(MHpm(gI1));
   }
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      result += CpChahhbarChaPL(gI1, gO1, gI1) * vev * AS12(decay_scale /
         Sqr(MCha(gI1))) / MCha(gI1);
   }
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += CpFehhbarFePL(gI1, gO1, gI1) * vev * AS12(decay_scale / Sqr(
         MFe(gI1))) / MFe(gI1);
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
      result += 0.3333333333333333 * scalar_fermion_qcd_factor(decay_mass,
         MFDX(gI1)) * CpFDXhhbarFDXPL(gI1, gO1, gI1) * vev * AS12(decay_scale /
         Sqr(MFDX(gI1))) / MFDX(gI1);
   }
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      result += 0.5 * CphhSHIpconjSHIp(gO1, gI1, gI1) * vev * AS0(
         decay_scale / Sqr(MSHIp(gI1))) / Sqr(MSHIp(gI1));
   }
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      result += CpChaIhhbarChaIPL(gI1, gO1, gI1) * vev * AS12(decay_scale /
         Sqr(MChaI(gI1))) / MChaI(gI1);
   }
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      result += 0.5 * CphhSHppconjSHpp(gO1, gI1, gI1) * vev * AS0(
         decay_scale / Sqr(MSHpp(gI1))) / Sqr(MSHpp(gI1));
   }
   result += -0.5 * CphhVWmconjVWm(gO1) * vev * AS1(decay_scale / Sqr(MVWm)) /
      Sqr(MVWm);

   result *= 0.1892681907127351 * physical_input.get(Physical_input::alpha_em_0
      ) * Sqrt(qedqcd.displayFermiConstant());

   ZH = saved_ZH;
   eff_CphhVPVP(gO1) = result;

}

void E6SSM_effective_couplings::calculate_eff_CphhVGVG(unsigned gO1)
{
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto g3 = MODELPARAMETER(g3);
   const auto MSd = MODELPARAMETER(MSd);
   const auto MSu = MODELPARAMETER(MSu);
   const auto MSDX = MODELPARAMETER(MSDX);
   const auto MFd = MODELPARAMETER(MFd);
   const auto MFu = MODELPARAMETER(MFu);
   const auto MFDX = MODELPARAMETER(MFDX);
   const double alpha_s = 0.07957747154594767*Sqr(g3);
   const auto decay_mass = PHYSICAL(Mhh)(gO1);
   const auto decay_scale = 0.25 * Sqr(decay_mass);
   const auto saved_ZH = ZH;
   ZH = PHYSICAL(ZH);

   const auto vev = Sqrt(Sqr(vd) + Sqr(vu));

   std::complex<double> result = 0;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      result += 0.5 * CphhSdconjSd(gO1, gI1, gI1) * vev * AS0(decay_scale /
         Sqr(MSd(gI1))) / Sqr(MSd(gI1));
   }
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      result += 0.5 * CphhSuconjSu(gO1, gI1, gI1) * vev * AS0(decay_scale /
         Sqr(MSu(gI1))) / Sqr(MSu(gI1));
   }
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      result += 0.5 * CphhSDXconjSDX(gO1, gI1, gI1) * vev * AS0(decay_scale
         / Sqr(MSDX(gI1))) / Sqr(MSDX(gI1));
   }
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += CpFdhhbarFdPL(gI1, gO1, gI1) * vev * AS12(decay_scale / Sqr(
         MFd(gI1))) / MFd(gI1);
   }
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += CpFuhhbarFuPL(gI1, gO1, gI1) * vev * AS12(decay_scale / Sqr(
         MFu(gI1))) / MFu(gI1);
   }
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += CpFDXhhbarFDXPL(gI1, gO1, gI1) * vev * AS12(decay_scale /
         Sqr(MFDX(gI1))) / MFDX(gI1);
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

void E6SSM_effective_couplings::calculate_eff_CpAhVPVP(unsigned gO1)
{
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto MCha = MODELPARAMETER(MCha);
   const auto MFe = MODELPARAMETER(MFe);
   const auto MFd = MODELPARAMETER(MFd);
   const auto MFu = MODELPARAMETER(MFu);
   const auto MFDX = MODELPARAMETER(MFDX);
   const auto MChaI = MODELPARAMETER(MChaI);
   const auto decay_mass = PHYSICAL(MAh)(gO1);
   const auto decay_scale = 0.25 * Sqr(decay_mass);
   const auto saved_ZA = ZA;
   ZA = PHYSICAL(ZA);

   const auto vev = Sqrt(Sqr(vd) + Sqr(vu));

   std::complex<double> result = 0;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      result += CpAhChabarChaPL(gO1, gI1, gI1) * vev * AP12(decay_scale /
         Sqr(MCha(gI1))) / MCha(gI1);
   }
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += CpAhFebarFePL(gO1, gI1, gI1) * vev * AP12(decay_scale / Sqr(
         MFe(gI1))) / MFe(gI1);
   }
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
      result += 0.3333333333333333 * pseudoscalar_fermion_qcd_factor(
         decay_mass, MFDX(gI1)) * CpAhFDXbarFDXPL(gO1, gI1, gI1) * vev * AP12(
         decay_scale / Sqr(MFDX(gI1))) / MFDX(gI1);
   }
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      result += CpAhChaIbarChaIPL(gO1, gI1, gI1) * vev * AP12(decay_scale /
         Sqr(MChaI(gI1))) / MChaI(gI1);
   }
   result *= std::complex<double>(2.0,0.);

   result *= 0.1892681907127351 * physical_input.get(Physical_input::alpha_em_0
      ) * Sqrt(qedqcd.displayFermiConstant());

   ZA = saved_ZA;
   eff_CpAhVPVP(gO1) = result;

}

void E6SSM_effective_couplings::calculate_eff_CpAhVGVG(unsigned gO1)
{
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto g3 = MODELPARAMETER(g3);
   const auto MFd = MODELPARAMETER(MFd);
   const auto MFu = MODELPARAMETER(MFu);
   const auto MFDX = MODELPARAMETER(MFDX);
   const double alpha_s = 0.07957747154594767*Sqr(g3);
   const auto decay_mass = PHYSICAL(MAh)(gO1);
   const auto decay_scale = 0.25 * Sqr(decay_mass);
   const auto saved_ZA = ZA;
   ZA = PHYSICAL(ZA);

   const auto vev = Sqrt(Sqr(vd) + Sqr(vu));

   std::complex<double> result = 0;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += CpAhFdbarFdPL(gO1, gI1, gI1) * vev * AP12(decay_scale / Sqr(
         MFd(gI1))) / MFd(gI1);
   }
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += CpAhFubarFuPL(gO1, gI1, gI1) * vev * AP12(decay_scale / Sqr(
         MFu(gI1))) / MFu(gI1);
   }
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      result += CpAhFDXbarFDXPL(gO1, gI1, gI1) * vev * AP12(decay_scale /
         Sqr(MFDX(gI1))) / MFDX(gI1);
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
