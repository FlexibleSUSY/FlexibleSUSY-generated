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

// File generated at Tue 8 Mar 2016 18:08:40

#include "E6SSM_effective_couplings.hpp"

#include "effective_couplings.hpp"
#include "standard_model.hpp"
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

E6SSM_effective_couplings::~E6SSM_effective_couplings()
{
}

void E6SSM_effective_couplings::calculate_effective_couplings()
{
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
      run_SM_strong_coupling_to(0.5 * Mhh(gO1));
      calculate_eff_CphhVPVP(gO1);
      run_SM_strong_coupling_to(Mhh(gO1));
      calculate_eff_CphhVGVG(gO1);
   }

   const auto MAh = PHYSICAL(MAh);
   for (unsigned gO1 = 2; gO1 < 3; ++gO1) {
      if (rg_improve && scale > MAh(gO1)) {
         model.run_to(MAh(gO1));
      }
      model.calculate_DRbar_masses();
      copy_mixing_matrices_from_model();
      run_SM_strong_coupling_to(0.5 * MAh(gO1));
      calculate_eff_CpAhVPVP(gO1);
      run_SM_strong_coupling_to(MAh(gO1));
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

void E6SSM_effective_couplings::run_SM_strong_coupling_to(double m)
{
   using namespace standard_model;

   Standard_model sm;

   sm.set_loops(2);
   sm.set_thresholds(2);
   sm.set_low_energy_data(qedqcd);
   sm.set_physical_input(physical_input);

   sm.initialise_from_input();
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

   std::complex<double> tmp_7247;
   std::complex<double> tmp_7248;
   std::complex<double> tmp_7249;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7249 += Conj(ZD(gt2,j1))*ZD(gt3,j1);
   }
   tmp_7248 += tmp_7249;
   tmp_7247 += (vd*(6*Sqr(g1) + 30*Sqr(g2) + 9*Sqr(gN))*ZH(gt1,0) - 2*vu*(3*Sqr
      (g1) + 15*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,1) - 15*vs*Sqr(gN)*ZH(gt1,2)) *
      tmp_7248;
   std::complex<double> tmp_7250;
   std::complex<double> tmp_7251;
   std::complex<double> tmp_7252;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7252 += Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j1))*ZD(gt3,j1);
   }
   tmp_7251 += tmp_7252;
   tmp_7250 += (-42.42640687119285*ZH(gt1,0)) * tmp_7251;
   std::complex<double> tmp_7253;
   std::complex<double> tmp_7254;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7254 += Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1);
   }
   tmp_7253 += tmp_7254;
   tmp_7250 += (vd*(6*Sqr(g1) + 9*Sqr(gN))*ZH(gt1,0) + (-6*vu*Sqr(g1) + 6*vu*
      Sqr(gN))*ZH(gt1,1) - 15*vs*Sqr(gN)*ZH(gt1,2)) * tmp_7253;
   std::complex<double> tmp_7255;
   std::complex<double> tmp_7256;
   std::complex<double> tmp_7257;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7257 += Conj(ZD(gt2,j1))*ZD(gt3,3 + j1)*TYd(j1,j1);
   }
   tmp_7256 += tmp_7257;
   tmp_7255 += (1.4142135623730951*ZH(gt1,0)) * tmp_7256;
   std::complex<double> tmp_7258;
   std::complex<double> tmp_7259;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_7259 += AbsSqr(Yd(j2,j2))*Conj(ZD(gt2,j2))*ZD(gt3,j2);
   }
   tmp_7258 += tmp_7259;
   tmp_7255 += (2*vd*ZH(gt1,0)) * tmp_7258;
   std::complex<double> tmp_7260;
   std::complex<double> tmp_7261;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_7261 += AbsSqr(Yd(j2,j2))*Conj(ZD(gt2,3 + j2))*ZD(gt3,3 + j2);
   }
   tmp_7260 += tmp_7261;
   tmp_7255 += (2*vd*ZH(gt1,0)) * tmp_7260;
   std::complex<double> tmp_7262;
   std::complex<double> tmp_7263;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7263 += Conj(Yd(j1,j1))*Conj(ZD(gt2,3 + j1))*ZD(gt3,j1);
   }
   tmp_7262 += tmp_7263;
   tmp_7255 += (-(vs*Lambdax*ZH(gt1,1))) * tmp_7262;
   std::complex<double> tmp_7264;
   std::complex<double> tmp_7265;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7265 += Conj(ZD(gt2,j1))*Yd(j1,j1)*ZD(gt3,3 + j1);
   }
   tmp_7264 += tmp_7265;
   tmp_7255 += (-(vs*Conj(Lambdax)*ZH(gt1,1))) * tmp_7264;
   std::complex<double> tmp_7266;
   std::complex<double> tmp_7267;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7267 += Conj(Yd(j1,j1))*Conj(ZD(gt2,3 + j1))*ZD(gt3,j1);
   }
   tmp_7266 += tmp_7267;
   tmp_7255 += (-(vu*Lambdax*ZH(gt1,2))) * tmp_7266;
   std::complex<double> tmp_7268;
   std::complex<double> tmp_7269;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7269 += Conj(ZD(gt2,j1))*Yd(j1,j1)*ZD(gt3,3 + j1);
   }
   tmp_7268 += tmp_7269;
   tmp_7255 += (-(vu*Conj(Lambdax)*ZH(gt1,2))) * tmp_7268;
   tmp_7250 += (-30) * tmp_7255;
   tmp_7247 += (2) * tmp_7250;
   result += (0.008333333333333333) * tmp_7247;

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

   std::complex<double> tmp_7270;
   std::complex<double> tmp_7271;
   std::complex<double> tmp_7272;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7272 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_7271 += tmp_7272;
   tmp_7270 += (-24*vd*Sqr(g1)*ZH(gt1,0)) * tmp_7271;
   std::complex<double> tmp_7273;
   std::complex<double> tmp_7274;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7274 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_7273 += tmp_7274;
   tmp_7270 += (9*vd*Sqr(gN)*ZH(gt1,0)) * tmp_7273;
   std::complex<double> tmp_7275;
   std::complex<double> tmp_7276;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7276 += Conj(ZU(gt2,j1))*Yu(j1,j1)*ZU(gt3,3 + j1);
   }
   tmp_7275 += tmp_7276;
   tmp_7270 += (60*vs*Conj(Lambdax)*ZH(gt1,0)) * tmp_7275;
   std::complex<double> tmp_7277;
   std::complex<double> tmp_7278;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7278 += Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j1))*ZU(gt3,j1);
   }
   tmp_7277 += tmp_7278;
   tmp_7270 += (-84.8528137423857*ZH(gt1,1)) * tmp_7277;
   std::complex<double> tmp_7279;
   std::complex<double> tmp_7280;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7280 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_7279 += tmp_7280;
   tmp_7270 += (24*vu*Sqr(g1)*ZH(gt1,1)) * tmp_7279;
   std::complex<double> tmp_7281;
   std::complex<double> tmp_7282;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7282 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_7281 += tmp_7282;
   tmp_7270 += (6*vu*Sqr(gN)*ZH(gt1,1)) * tmp_7281;
   std::complex<double> tmp_7283;
   std::complex<double> tmp_7284;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7284 += Conj(ZU(gt2,j1))*ZU(gt3,3 + j1)*TYu(j1,j1);
   }
   tmp_7283 += tmp_7284;
   tmp_7270 += (-84.8528137423857*ZH(gt1,1)) * tmp_7283;
   std::complex<double> tmp_7285;
   std::complex<double> tmp_7286;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_7286 += AbsSqr(Yu(j2,j2))*Conj(ZU(gt2,j2))*ZU(gt3,j2);
   }
   tmp_7285 += tmp_7286;
   tmp_7270 += (-120*vu*ZH(gt1,1)) * tmp_7285;
   std::complex<double> tmp_7287;
   std::complex<double> tmp_7288;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_7288 += AbsSqr(Yu(j2,j2))*Conj(ZU(gt2,3 + j2))*ZU(gt3,3 + j2);
   }
   tmp_7287 += tmp_7288;
   tmp_7270 += (-120*vu*ZH(gt1,1)) * tmp_7287;
   std::complex<double> tmp_7289;
   std::complex<double> tmp_7290;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7290 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_7289 += tmp_7290;
   tmp_7270 += (-15*vs*Sqr(gN)*ZH(gt1,2)) * tmp_7289;
   std::complex<double> tmp_7291;
   std::complex<double> tmp_7292;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7292 += Conj(ZU(gt2,j1))*Yu(j1,j1)*ZU(gt3,3 + j1);
   }
   tmp_7291 += tmp_7292;
   tmp_7270 += (60*vd*Conj(Lambdax)*ZH(gt1,2)) * tmp_7291;
   std::complex<double> tmp_7293;
   std::complex<double> tmp_7294;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7294 += Conj(Yu(j1,j1))*Conj(ZU(gt2,3 + j1))*ZU(gt3,j1);
   }
   tmp_7293 += tmp_7294;
   tmp_7270 += (60*Lambdax*(vs*ZH(gt1,0) + vd*ZH(gt1,2))) * tmp_7293;
   std::complex<double> tmp_7295;
   std::complex<double> tmp_7296;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7296 += Conj(ZU(gt2,j1))*ZU(gt3,j1);
   }
   tmp_7295 += tmp_7296;
   tmp_7270 += (vd*(6*Sqr(g1) - 30*Sqr(g2) + 9*Sqr(gN))*ZH(gt1,0) + 2*vu*(-3*
      Sqr(g1) + 15*Sqr(g2) + 3*Sqr(gN))*ZH(gt1,1) - 15*vs*Sqr(gN)*ZH(gt1,2)) *
      tmp_7295;
   result += (0.008333333333333333) * tmp_7270;

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

   std::complex<double> tmp_7297;
   std::complex<double> tmp_7298;
   std::complex<double> tmp_7299;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7299 += Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j1))*ZE(gt3,j1);
   }
   tmp_7298 += tmp_7299;
   tmp_7297 += (-28.284271247461902*ZH(gt1,0)) * tmp_7298;
   std::complex<double> tmp_7300;
   std::complex<double> tmp_7301;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7301 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_7300 += tmp_7301;
   tmp_7297 += (12*vd*Sqr(g1)*ZH(gt1,0)) * tmp_7300;
   std::complex<double> tmp_7302;
   std::complex<double> tmp_7303;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7303 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_7302 += tmp_7303;
   tmp_7297 += (3*vd*Sqr(gN)*ZH(gt1,0)) * tmp_7302;
   std::complex<double> tmp_7304;
   std::complex<double> tmp_7305;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7305 += Conj(ZE(gt2,j1))*ZE(gt3,3 + j1)*TYe(j1,j1);
   }
   tmp_7304 += tmp_7305;
   tmp_7297 += (-28.284271247461902*ZH(gt1,0)) * tmp_7304;
   std::complex<double> tmp_7306;
   std::complex<double> tmp_7307;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_7307 += AbsSqr(Ye(j2,j2))*Conj(ZE(gt2,j2))*ZE(gt3,j2);
   }
   tmp_7306 += tmp_7307;
   tmp_7297 += (-40*vd*ZH(gt1,0)) * tmp_7306;
   std::complex<double> tmp_7308;
   std::complex<double> tmp_7309;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_7309 += AbsSqr(Ye(j2,j2))*Conj(ZE(gt2,3 + j2))*ZE(gt3,3 + j2);
   }
   tmp_7308 += tmp_7309;
   tmp_7297 += (-40*vd*ZH(gt1,0)) * tmp_7308;
   std::complex<double> tmp_7310;
   std::complex<double> tmp_7311;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7311 += Conj(Ye(j1,j1))*Conj(ZE(gt2,3 + j1))*ZE(gt3,j1);
   }
   tmp_7310 += tmp_7311;
   tmp_7297 += (20*vs*Lambdax*ZH(gt1,1)) * tmp_7310;
   std::complex<double> tmp_7312;
   std::complex<double> tmp_7313;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7313 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_7312 += tmp_7313;
   tmp_7297 += (-12*vu*Sqr(g1)*ZH(gt1,1)) * tmp_7312;
   std::complex<double> tmp_7314;
   std::complex<double> tmp_7315;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7315 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_7314 += tmp_7315;
   tmp_7297 += (2*vu*Sqr(gN)*ZH(gt1,1)) * tmp_7314;
   std::complex<double> tmp_7316;
   std::complex<double> tmp_7317;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7317 += Conj(ZE(gt2,j1))*Ye(j1,j1)*ZE(gt3,3 + j1);
   }
   tmp_7316 += tmp_7317;
   tmp_7297 += (20*vs*Conj(Lambdax)*ZH(gt1,1)) * tmp_7316;
   std::complex<double> tmp_7318;
   std::complex<double> tmp_7319;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7319 += Conj(Ye(j1,j1))*Conj(ZE(gt2,3 + j1))*ZE(gt3,j1);
   }
   tmp_7318 += tmp_7319;
   tmp_7297 += (20*vu*Lambdax*ZH(gt1,2)) * tmp_7318;
   std::complex<double> tmp_7320;
   std::complex<double> tmp_7321;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7321 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_7320 += tmp_7321;
   tmp_7297 += (-5*vs*Sqr(gN)*ZH(gt1,2)) * tmp_7320;
   std::complex<double> tmp_7322;
   std::complex<double> tmp_7323;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7323 += Conj(ZE(gt2,j1))*Ye(j1,j1)*ZE(gt3,3 + j1);
   }
   tmp_7322 += tmp_7323;
   tmp_7297 += (20*vu*Conj(Lambdax)*ZH(gt1,2)) * tmp_7322;
   std::complex<double> tmp_7324;
   std::complex<double> tmp_7325;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7325 += Conj(ZE(gt2,j1))*ZE(gt3,j1);
   }
   tmp_7324 += tmp_7325;
   tmp_7297 += (-2*(vd*(3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,0) + vu*(-3*
      Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN))*ZH(gt1,1) + 5*vs*Sqr(gN)*ZH(gt1,2))) *
      tmp_7324;
   result += (0.025) * tmp_7297;

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

   std::complex<double> tmp_7326;
   std::complex<double> tmp_7327;
   std::complex<double> tmp_7328;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7328 += Conj(ZDX(gt2,3 + j1))*ZDX(gt3,3 + j1);
   }
   tmp_7327 += tmp_7328;
   tmp_7326 += (12*vd*Sqr(g1)*ZH(gt1,0)) * tmp_7327;
   std::complex<double> tmp_7329;
   std::complex<double> tmp_7330;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7330 += Conj(ZDX(gt2,3 + j1))*ZDX(gt3,3 + j1);
   }
   tmp_7329 += tmp_7330;
   tmp_7326 += (-27*vd*Sqr(gN)*ZH(gt1,0)) * tmp_7329;
   std::complex<double> tmp_7331;
   std::complex<double> tmp_7332;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7332 += Conj(ZDX(gt2,j1))*ZDX(gt3,3 + j1)*Kappa(j1,j1);
   }
   tmp_7331 += tmp_7332;
   tmp_7326 += (60*vu*Conj(Lambdax)*ZH(gt1,0)) * tmp_7331;
   std::complex<double> tmp_7333;
   std::complex<double> tmp_7334;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7334 += Conj(ZDX(gt2,3 + j1))*ZDX(gt3,3 + j1);
   }
   tmp_7333 += tmp_7334;
   tmp_7326 += (-12*vu*Sqr(g1)*ZH(gt1,1)) * tmp_7333;
   std::complex<double> tmp_7335;
   std::complex<double> tmp_7336;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7336 += Conj(ZDX(gt2,3 + j1))*ZDX(gt3,3 + j1);
   }
   tmp_7335 += tmp_7336;
   tmp_7326 += (-18*vu*Sqr(gN)*ZH(gt1,1)) * tmp_7335;
   std::complex<double> tmp_7337;
   std::complex<double> tmp_7338;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7338 += Conj(ZDX(gt2,j1))*ZDX(gt3,3 + j1)*Kappa(j1,j1);
   }
   tmp_7337 += tmp_7338;
   tmp_7326 += (60*vd*Conj(Lambdax)*ZH(gt1,1)) * tmp_7337;
   std::complex<double> tmp_7339;
   std::complex<double> tmp_7340;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7340 += Conj(ZDX(gt2,3 + j1))*Conj(Kappa(j1,j1))*ZDX(gt3,j1);
   }
   tmp_7339 += tmp_7340;
   tmp_7326 += (60*Lambdax*(vu*ZH(gt1,0) + vd*ZH(gt1,1))) * tmp_7339;
   std::complex<double> tmp_7341;
   std::complex<double> tmp_7342;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7342 += Conj(ZDX(gt2,3 + j1))*Conj(TKappa(j1,j1))*ZDX(gt3,j1);
   }
   tmp_7341 += tmp_7342;
   tmp_7326 += (-84.8528137423857*ZH(gt1,2)) * tmp_7341;
   std::complex<double> tmp_7343;
   std::complex<double> tmp_7344;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7344 += Conj(ZDX(gt2,3 + j1))*ZDX(gt3,3 + j1);
   }
   tmp_7343 += tmp_7344;
   tmp_7326 += (45*vs*Sqr(gN)*ZH(gt1,2)) * tmp_7343;
   std::complex<double> tmp_7345;
   std::complex<double> tmp_7346;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7346 += Conj(ZDX(gt2,j1))*ZDX(gt3,3 + j1)*TKappa(j1,j1);
   }
   tmp_7345 += tmp_7346;
   tmp_7326 += (-84.8528137423857*ZH(gt1,2)) * tmp_7345;
   std::complex<double> tmp_7347;
   std::complex<double> tmp_7348;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_7348 += AbsSqr(Kappa(j2,j2))*Conj(ZDX(gt2,j2))*ZDX(gt3,j2);
   }
   tmp_7347 += tmp_7348;
   tmp_7326 += (-120*vs*ZH(gt1,2)) * tmp_7347;
   std::complex<double> tmp_7349;
   std::complex<double> tmp_7350;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_7350 += AbsSqr(Kappa(j2,j2))*Conj(ZDX(gt2,3 + j2))*ZDX(gt3,3 + j2)
         ;
   }
   tmp_7349 += tmp_7350;
   tmp_7326 += (-120*vs*ZH(gt1,2)) * tmp_7349;
   std::complex<double> tmp_7351;
   std::complex<double> tmp_7352;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7352 += Conj(ZDX(gt2,j1))*ZDX(gt3,j1);
   }
   tmp_7351 += tmp_7352;
   tmp_7326 += (-2*(vd*(6*Sqr(g1) + 9*Sqr(gN))*ZH(gt1,0) + (-6*vu*Sqr(g1) + 6*
      vu*Sqr(gN))*ZH(gt1,1) - 15*vs*Sqr(gN)*ZH(gt1,2))) * tmp_7351;
   result += (0.008333333333333333) * tmp_7326;

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

   std::complex<double> tmp_7353;
   std::complex<double> tmp_7354;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7354 += Conj(ZEL(gt1,j1))*Conj(ZER(gt3,j1))*Ye(j1,j1);
   }
   tmp_7353 += tmp_7354;
   result += (-0.7071067811865475*ZH(gt2,0)) * tmp_7353;

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpFdhhbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_7355;
   std::complex<double> tmp_7356;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7356 += Conj(ZDL(gt1,j1))*Conj(ZDR(gt3,j1))*Yd(j1,j1);
   }
   tmp_7355 += tmp_7356;
   result += (-0.7071067811865475*ZH(gt2,0)) * tmp_7355;

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpFuhhbarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_7357;
   std::complex<double> tmp_7358;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7358 += Conj(ZUL(gt1,j1))*Conj(ZUR(gt3,j1))*Yu(j1,j1);
   }
   tmp_7357 += tmp_7358;
   result += (-0.7071067811865475*ZH(gt2,1)) * tmp_7357;

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpFDXhhbarFDXPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Kappa = MODELPARAMETER(Kappa);

   std::complex<double> result;

   std::complex<double> tmp_7359;
   std::complex<double> tmp_7360;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7360 += Conj(ZDXL(gt1,j1))*Conj(ZDXR(gt3,j1))*Kappa(j1,j1);
   }
   tmp_7359 += tmp_7360;
   result += (-0.7071067811865475*ZH(gt2,2)) * tmp_7359;

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

   std::complex<double> tmp_7361;
   std::complex<double> tmp_7362;
   std::complex<double> tmp_7363;
   for (unsigned j1 = 0; j1 < 2; ++j1) {
      tmp_7363 += Conj(UHIp(gt2,j1))*UHIp(gt3,j1);
   }
   tmp_7362 += tmp_7363;
   tmp_7361 += (vd*(-6*Sqr(g1) + 10*Sqr(g2) - 9*Sqr(gN))*ZH(gt1,0) + 2*vu*(3*
      Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,1) + 15*vs*Sqr(gN)*ZH(gt1,2)) *
      tmp_7362;
   std::complex<double> tmp_7364;
   std::complex<double> tmp_7365;
   std::complex<double> tmp_7366;
   for (unsigned j1 = 0; j1 < 2; ++j1) {
      tmp_7366 += Conj(UHIp(gt2,2 + j1))*Conj(Lambda12(j1,j1))*UHIp(gt3,j1);
   }
   tmp_7365 += tmp_7366;
   tmp_7364 += (10*Lambdax*(vu*ZH(gt1,0) + vd*ZH(gt1,1))) * tmp_7365;
   std::complex<double> tmp_7367;
   std::complex<double> tmp_7368;
   for (unsigned j1 = 0; j1 < 2; ++j1) {
      tmp_7368 += Conj(UHIp(gt2,2 + j1))*UHIp(gt3,2 + j1);
   }
   tmp_7367 += tmp_7368;
   tmp_7364 += (vd*(3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,0) + vu*(-3*Sqr(
      g1) + 5*Sqr(g2) - 2*Sqr(gN))*ZH(gt1,1) + 5*vs*Sqr(gN)*ZH(gt1,2)) * tmp_7367;
   std::complex<double> tmp_7369;
   std::complex<double> tmp_7370;
   std::complex<double> tmp_7371;
   for (unsigned j1 = 0; j1 < 2; ++j1) {
      tmp_7371 += Conj(UHIp(gt2,j1))*UHIp(gt3,2 + j1)*Lambda12(j1,j1);
   }
   tmp_7370 += tmp_7371;
   tmp_7369 += (Conj(Lambdax)*(vu*ZH(gt1,0) + vd*ZH(gt1,1))) * tmp_7370;
   std::complex<double> tmp_7372;
   std::complex<double> tmp_7373;
   std::complex<double> tmp_7374;
   for (unsigned j1 = 0; j1 < 2; ++j1) {
      tmp_7374 += Conj(UHIp(gt2,2 + j1))*Conj(TLambda12(j1,j1))*UHIp(gt3,j1)
         ;
   }
   tmp_7373 += tmp_7374;
   tmp_7372 += (1.4142135623730951) * tmp_7373;
   std::complex<double> tmp_7375;
   std::complex<double> tmp_7376;
   for (unsigned j1 = 0; j1 < 2; ++j1) {
      tmp_7376 += Conj(UHIp(gt2,j1))*UHIp(gt3,2 + j1)*TLambda12(j1,j1);
   }
   tmp_7375 += tmp_7376;
   tmp_7372 += (1.4142135623730951) * tmp_7375;
   std::complex<double> tmp_7377;
   std::complex<double> tmp_7378;
   for (unsigned j2 = 0; j2 < 2; ++j2) {
      tmp_7378 += AbsSqr(Lambda12(j2,j2))*Conj(UHIp(gt2,j2))*UHIp(gt3,j2);
   }
   tmp_7377 += tmp_7378;
   std::complex<double> tmp_7379;
   for (unsigned j2 = 0; j2 < 2; ++j2) {
      tmp_7379 += AbsSqr(Lambda12(j2,j2))*Conj(UHIp(gt2,2 + j2))*UHIp(gt3,2
         + j2);
   }
   tmp_7377 += tmp_7379;
   tmp_7372 += (2*vs) * tmp_7377;
   tmp_7369 += (-ZH(gt1,2)) * tmp_7372;
   tmp_7364 += (10) * tmp_7369;
   tmp_7361 += (2) * tmp_7364;
   result += (0.025) * tmp_7361;

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpChaIhhbarChaIPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Lambda12 = MODELPARAMETER(Lambda12);

   std::complex<double> result;

   std::complex<double> tmp_7380;
   std::complex<double> tmp_7381;
   for (unsigned j1 = 0; j1 < 2; ++j1) {
      tmp_7381 += Conj(ZMI(gt1,j1))*Conj(ZPI(gt3,j1))*Lambda12(j1,j1);
   }
   tmp_7380 += tmp_7381;
   result += (-0.7071067811865475*ZH(gt2,2)) * tmp_7380;

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

   std::complex<double> tmp_7382;
   std::complex<double> tmp_7383;
   std::complex<double> tmp_7384;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7384 += Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j1))*ZD(gt3,j1);
   }
   tmp_7383 += tmp_7384;
   tmp_7382 += (-1.4142135623730951*ZA(gt1,0)) * tmp_7383;
   std::complex<double> tmp_7385;
   std::complex<double> tmp_7386;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7386 += Conj(ZD(gt2,j1))*ZD(gt3,3 + j1)*TYd(j1,j1);
   }
   tmp_7385 += tmp_7386;
   tmp_7382 += (1.4142135623730951*ZA(gt1,0)) * tmp_7385;
   std::complex<double> tmp_7387;
   std::complex<double> tmp_7388;
   std::complex<double> tmp_7389;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7389 += Conj(Yd(j1,j1))*Conj(ZD(gt2,3 + j1))*ZD(gt3,j1);
   }
   tmp_7388 += tmp_7389;
   tmp_7387 += (Lambdax) * tmp_7388;
   std::complex<double> tmp_7390;
   std::complex<double> tmp_7391;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7391 += Conj(ZD(gt2,j1))*Yd(j1,j1)*ZD(gt3,3 + j1);
   }
   tmp_7390 += tmp_7391;
   tmp_7387 += (-Conj(Lambdax)) * tmp_7390;
   tmp_7382 += (-(vs*ZA(gt1,1)) - vu*ZA(gt1,2)) * tmp_7387;
   result += (std::complex<double>(0,-0.5)) * tmp_7382;

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

   std::complex<double> tmp_7392;
   std::complex<double> tmp_7393;
   std::complex<double> tmp_7394;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7394 += Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j1))*ZU(gt3,j1);
   }
   tmp_7393 += tmp_7394;
   std::complex<double> tmp_7395;
   std::complex<double> tmp_7396;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7396 += Conj(ZU(gt2,j1))*ZU(gt3,3 + j1)*TYu(j1,j1);
   }
   tmp_7395 += tmp_7396;
   tmp_7393 += (-1) * tmp_7395;
   tmp_7392 += (-1.4142135623730951*ZA(gt1,1)) * tmp_7393;
   std::complex<double> tmp_7397;
   std::complex<double> tmp_7398;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7398 += Conj(Yu(j1,j1))*Conj(ZU(gt2,3 + j1))*ZU(gt3,j1);
   }
   tmp_7397 += tmp_7398;
   tmp_7392 += (-(Lambdax*(vs*ZA(gt1,0) + vd*ZA(gt1,2)))) * tmp_7397;
   std::complex<double> tmp_7399;
   std::complex<double> tmp_7400;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7400 += Conj(ZU(gt2,j1))*Yu(j1,j1)*ZU(gt3,3 + j1);
   }
   tmp_7399 += tmp_7400;
   tmp_7392 += (Conj(Lambdax)*(vs*ZA(gt1,0) + vd*ZA(gt1,2))) * tmp_7399;
   result += (std::complex<double>(0,-0.5)) * tmp_7392;

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

   std::complex<double> tmp_7401;
   std::complex<double> tmp_7402;
   std::complex<double> tmp_7403;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7403 += Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j1))*ZE(gt3,j1);
   }
   tmp_7402 += tmp_7403;
   tmp_7401 += (-1.4142135623730951*ZA(gt1,0)) * tmp_7402;
   std::complex<double> tmp_7404;
   std::complex<double> tmp_7405;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7405 += Conj(ZE(gt2,j1))*ZE(gt3,3 + j1)*TYe(j1,j1);
   }
   tmp_7404 += tmp_7405;
   tmp_7401 += (1.4142135623730951*ZA(gt1,0)) * tmp_7404;
   std::complex<double> tmp_7406;
   std::complex<double> tmp_7407;
   std::complex<double> tmp_7408;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7408 += Conj(Ye(j1,j1))*Conj(ZE(gt2,3 + j1))*ZE(gt3,j1);
   }
   tmp_7407 += tmp_7408;
   tmp_7406 += (Lambdax) * tmp_7407;
   std::complex<double> tmp_7409;
   std::complex<double> tmp_7410;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7410 += Conj(ZE(gt2,j1))*Ye(j1,j1)*ZE(gt3,3 + j1);
   }
   tmp_7409 += tmp_7410;
   tmp_7406 += (-Conj(Lambdax)) * tmp_7409;
   tmp_7401 += (-(vs*ZA(gt1,1)) - vu*ZA(gt1,2)) * tmp_7406;
   result += (std::complex<double>(0,-0.5)) * tmp_7401;

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

   std::complex<double> tmp_7411;
   std::complex<double> tmp_7412;
   std::complex<double> tmp_7413;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7413 += Conj(ZDX(gt2,3 + j1))*Conj(Kappa(j1,j1))*ZDX(gt3,j1);
   }
   tmp_7412 += tmp_7413;
   tmp_7411 += (-(Lambdax*(vu*ZA(gt1,0) + vd*ZA(gt1,1)))) * tmp_7412;
   std::complex<double> tmp_7414;
   std::complex<double> tmp_7415;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7415 += Conj(ZDX(gt2,j1))*ZDX(gt3,3 + j1)*Kappa(j1,j1);
   }
   tmp_7414 += tmp_7415;
   tmp_7411 += (Conj(Lambdax)*(vu*ZA(gt1,0) + vd*ZA(gt1,1))) * tmp_7414;
   std::complex<double> tmp_7416;
   std::complex<double> tmp_7417;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7417 += Conj(ZDX(gt2,3 + j1))*Conj(TKappa(j1,j1))*ZDX(gt3,j1);
   }
   tmp_7416 += tmp_7417;
   std::complex<double> tmp_7418;
   std::complex<double> tmp_7419;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7419 += Conj(ZDX(gt2,j1))*ZDX(gt3,3 + j1)*TKappa(j1,j1);
   }
   tmp_7418 += tmp_7419;
   tmp_7416 += (-1) * tmp_7418;
   tmp_7411 += (-1.4142135623730951*ZA(gt1,2)) * tmp_7416;
   result += (std::complex<double>(0,-0.5)) * tmp_7411;

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

   std::complex<double> tmp_7420;
   std::complex<double> tmp_7421;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7421 += Conj(ZEL(gt2,j1))*Conj(ZER(gt3,j1))*Ye(j1,j1);
   }
   tmp_7420 += tmp_7421;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,0)) *
      tmp_7420;

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpAhFdbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_7422;
   std::complex<double> tmp_7423;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7423 += Conj(ZDL(gt2,j1))*Conj(ZDR(gt3,j1))*Yd(j1,j1);
   }
   tmp_7422 += tmp_7423;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,0)) *
      tmp_7422;

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpAhFubarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_7424;
   std::complex<double> tmp_7425;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7425 += Conj(ZUL(gt2,j1))*Conj(ZUR(gt3,j1))*Yu(j1,j1);
   }
   tmp_7424 += tmp_7425;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,1)) *
      tmp_7424;

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpAhFDXbarFDXPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Kappa = MODELPARAMETER(Kappa);

   std::complex<double> result;

   std::complex<double> tmp_7426;
   std::complex<double> tmp_7427;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7427 += Conj(ZDXL(gt2,j1))*Conj(ZDXR(gt3,j1))*Kappa(j1,j1);
   }
   tmp_7426 += tmp_7427;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,2)) *
      tmp_7426;

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

   std::complex<double> tmp_7428;
   std::complex<double> tmp_7429;
   std::complex<double> tmp_7430;
   for (unsigned j1 = 0; j1 < 2; ++j1) {
      tmp_7430 += Conj(UHIp(gt2,2 + j1))*Conj(Lambda12(j1,j1))*UHIp(gt3,j1);
   }
   tmp_7429 += tmp_7430;
   tmp_7428 += (-(Lambdax*(vu*ZA(gt1,0) + vd*ZA(gt1,1)))) * tmp_7429;
   std::complex<double> tmp_7431;
   std::complex<double> tmp_7432;
   for (unsigned j1 = 0; j1 < 2; ++j1) {
      tmp_7432 += Conj(UHIp(gt2,j1))*UHIp(gt3,2 + j1)*Lambda12(j1,j1);
   }
   tmp_7431 += tmp_7432;
   tmp_7428 += (Conj(Lambdax)*(vu*ZA(gt1,0) + vd*ZA(gt1,1))) * tmp_7431;
   std::complex<double> tmp_7433;
   std::complex<double> tmp_7434;
   for (unsigned j1 = 0; j1 < 2; ++j1) {
      tmp_7434 += Conj(UHIp(gt2,2 + j1))*Conj(TLambda12(j1,j1))*UHIp(gt3,j1)
         ;
   }
   tmp_7433 += tmp_7434;
   std::complex<double> tmp_7435;
   std::complex<double> tmp_7436;
   for (unsigned j1 = 0; j1 < 2; ++j1) {
      tmp_7436 += Conj(UHIp(gt2,j1))*UHIp(gt3,2 + j1)*TLambda12(j1,j1);
   }
   tmp_7435 += tmp_7436;
   tmp_7433 += (-1) * tmp_7435;
   tmp_7428 += (-1.4142135623730951*ZA(gt1,2)) * tmp_7433;
   result += (std::complex<double>(0,-0.5)) * tmp_7428;

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpAhChaIbarChaIPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Lambda12 = MODELPARAMETER(Lambda12);

   std::complex<double> result;

   std::complex<double> tmp_7437;
   std::complex<double> tmp_7438;
   for (unsigned j1 = 0; j1 < 2; ++j1) {
      tmp_7438 += Conj(ZMI(gt2,j1))*Conj(ZPI(gt3,j1))*Lambda12(j1,j1);
   }
   tmp_7437 += tmp_7438;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,2)) *
      tmp_7437;

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
