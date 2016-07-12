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

// File generated at Tue 12 Jul 2016 10:58:04

#include "TMSSM_effective_couplings.hpp"

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

TMSSM_effective_couplings::TMSSM_effective_couplings(
   const TMSSM_mass_eigenstates& model_,
   const softsusy::QedQcd& qedqcd_,
   const Physical_input& input_)
   : model(model_), qedqcd(qedqcd_), physical_input(input_)
   , rg_improve(true), include_qcd_corrections(true)
   , ZD(MODELPARAMETER(ZD)), ZV(MODELPARAMETER(ZV)), ZU(MODELPARAMETER(ZU)), ZE
      (MODELPARAMETER(ZE)), ZH(MODELPARAMETER(ZH)), ZA(MODELPARAMETER(ZA)), ZP(
      MODELPARAMETER(ZP)), ZN(MODELPARAMETER(ZN)), UM(MODELPARAMETER(UM)), UP(
      MODELPARAMETER(UP)), ZEL(MODELPARAMETER(ZEL)), ZER(MODELPARAMETER(ZER)), ZDL
      (MODELPARAMETER(ZDL)), ZDR(MODELPARAMETER(ZDR)), ZUL(MODELPARAMETER(ZUL)),
      ZUR(MODELPARAMETER(ZUR)), ZZ(MODELPARAMETER(ZZ))

   , eff_CphhVPVP(Eigen::Array<std::complex<double>,3,1>::Zero()), eff_CphhVGVG
      (Eigen::Array<std::complex<double>,3,1>::Zero()), eff_CpAhVPVP(Eigen::Array<
      std::complex<double>,3,1>::Zero()), eff_CpAhVGVG(Eigen::Array<std::complex<
      double>,3,1>::Zero())

{
}

TMSSM_effective_couplings::~TMSSM_effective_couplings()
{
}

void TMSSM_effective_couplings::calculate_effective_couplings()
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
   for (unsigned gO1 = 1; gO1 < 3; ++gO1) {
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

void TMSSM_effective_couplings::set_model(const TMSSM_mass_eigenstates& model_)
{
   model = model_;
   copy_mixing_matrices_from_model();
}

void TMSSM_effective_couplings::copy_mixing_matrices_from_model()
{
   ZD = MODELPARAMETER(ZD);
   ZV = MODELPARAMETER(ZV);
   ZU = MODELPARAMETER(ZU);
   ZE = MODELPARAMETER(ZE);
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
   ZZ = MODELPARAMETER(ZZ);

}

void TMSSM_effective_couplings::run_SM_strong_coupling_to(double m)
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

std::complex<double> TMSSM_effective_couplings::scalar_scalar_qcd_factor(double m_decay, double m_loop) const
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

std::complex<double> TMSSM_effective_couplings::scalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) *
         scalar_diphoton_fermion_loop(m_decay, m_loop);

   }

   return result;
}

std::complex<double> TMSSM_effective_couplings::pseudoscalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) *
         pseudoscalar_diphoton_fermion_loop(m_decay, m_loop);

   }

   return result;
}

double TMSSM_effective_couplings::number_of_active_flavours(double m) const
{
   if (m < qedqcd.displayMbMb()) {
      return 4.0;
   } else if (m < qedqcd.displayPoleMt()) {
      return 5.0;
   } else {
      return 6.0;
   }
}

double TMSSM_effective_couplings::scalar_scaling_factor(double m) const
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

double TMSSM_effective_couplings::pseudoscalar_scaling_factor(double m) const
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

double TMSSM_effective_couplings::get_hhVPVP_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(Mhh)(gO1);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CphhVPVP(gO1));
}

double TMSSM_effective_couplings::get_hhVGVG_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(Mhh)(gO1);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CphhVGVG(gO1));
}

double TMSSM_effective_couplings::get_AhVPVP_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(MAh)(gO1);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CpAhVPVP(gO1));
}

double TMSSM_effective_couplings::get_AhVGVG_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(MAh)(gO1);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CpAhVGVG(gO1));
}

std::complex<double> TMSSM_effective_couplings::CphhSdconjSd(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYd = MODELPARAMETER(TYd);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vT = MODELPARAMETER(vT);
   const auto vu = MODELPARAMETER(vu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3246;
   std::complex<double> tmp_3247;
   std::complex<double> tmp_3248;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3248 += Conj(ZD(gt2,j1))*ZD(gt3,j1);
   }
   tmp_3247 += tmp_3248;
   tmp_3246 += ((0.6*Sqr(g1) + 3*Sqr(g2))*(vd*ZH(gt1,0) - vu*ZH(gt1,1))) *
      tmp_3247;
   std::complex<double> tmp_3249;
   std::complex<double> tmp_3250;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3250 += Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1);
   }
   tmp_3249 += tmp_3250;
   tmp_3246 += (1.2*Sqr(g1)*(vd*ZH(gt1,0) - vu*ZH(gt1,1))) * tmp_3249;
   std::complex<double> tmp_3251;
   std::complex<double> tmp_3252;
   std::complex<double> tmp_3253;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3254;
      std::complex<double> tmp_3255;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3255 += ZD(gt3,3 + j1)*TYd(j1,j2);
      }
      tmp_3254 += tmp_3255;
      tmp_3253 += (Conj(ZD(gt2,j2))) * tmp_3254;
   }
   tmp_3252 += tmp_3253;
   tmp_3251 += (-2.8284271247461903*ZH(gt1,0)) * tmp_3252;
   std::complex<double> tmp_3256;
   std::complex<double> tmp_3257;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3258;
      std::complex<double> tmp_3259;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3259 += Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_3258 += tmp_3259;
      tmp_3257 += (ZD(gt3,j2)) * tmp_3258;
   }
   tmp_3256 += tmp_3257;
   tmp_3251 += (-2.8284271247461903*ZH(gt1,0)) * tmp_3256;
   std::complex<double> tmp_3260;
   std::complex<double> tmp_3261;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3262;
      std::complex<double> tmp_3263;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3264;
         std::complex<double> tmp_3265;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3265 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_3264 += tmp_3265;
         tmp_3263 += (ZD(gt3,3 + j2)) * tmp_3264;
      }
      tmp_3262 += tmp_3263;
      tmp_3261 += (Conj(ZD(gt2,3 + j3))) * tmp_3262;
   }
   tmp_3260 += tmp_3261;
   tmp_3251 += (-4*vd*ZH(gt1,0)) * tmp_3260;
   std::complex<double> tmp_3266;
   std::complex<double> tmp_3267;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3268;
      std::complex<double> tmp_3269;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3270;
         std::complex<double> tmp_3271;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3271 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_3270 += tmp_3271;
         tmp_3269 += (Conj(ZD(gt2,j2))) * tmp_3270;
      }
      tmp_3268 += tmp_3269;
      tmp_3267 += (ZD(gt3,j3)) * tmp_3268;
   }
   tmp_3266 += tmp_3267;
   tmp_3251 += (-4*vd*ZH(gt1,0)) * tmp_3266;
   std::complex<double> tmp_3272;
   std::complex<double> tmp_3273;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3274;
      std::complex<double> tmp_3275;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3275 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3274 += tmp_3275;
      tmp_3273 += (Conj(ZD(gt2,j2))) * tmp_3274;
   }
   tmp_3272 += tmp_3273;
   tmp_3251 += (1.4142135623730951*vT*Conj(Lambdax)*ZH(gt1,1)) * tmp_3272;
   std::complex<double> tmp_3276;
   std::complex<double> tmp_3277;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3278;
      std::complex<double> tmp_3279;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3279 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3278 += tmp_3279;
      tmp_3277 += (Conj(ZD(gt2,j2))) * tmp_3278;
   }
   tmp_3276 += tmp_3277;
   tmp_3251 += (2.8284271247461903*Conj(Mu)*ZH(gt1,1)) * tmp_3276;
   std::complex<double> tmp_3280;
   std::complex<double> tmp_3281;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3282;
      std::complex<double> tmp_3283;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3283 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3282 += tmp_3283;
      tmp_3281 += (ZD(gt3,j2)) * tmp_3282;
   }
   tmp_3280 += tmp_3281;
   tmp_3251 += (1.4142135623730951*vT*Lambdax*ZH(gt1,1)) * tmp_3280;
   std::complex<double> tmp_3284;
   std::complex<double> tmp_3285;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3286;
      std::complex<double> tmp_3287;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3287 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3286 += tmp_3287;
      tmp_3285 += (ZD(gt3,j2)) * tmp_3286;
   }
   tmp_3284 += tmp_3285;
   tmp_3251 += (2.8284271247461903*Mu*ZH(gt1,1)) * tmp_3284;
   std::complex<double> tmp_3288;
   std::complex<double> tmp_3289;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3290;
      std::complex<double> tmp_3291;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3291 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3290 += tmp_3291;
      tmp_3289 += (Conj(ZD(gt2,j2))) * tmp_3290;
   }
   tmp_3288 += tmp_3289;
   tmp_3251 += (1.4142135623730951*vu*Conj(Lambdax)*ZH(gt1,2)) * tmp_3288;
   std::complex<double> tmp_3292;
   std::complex<double> tmp_3293;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3294;
      std::complex<double> tmp_3295;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3295 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3294 += tmp_3295;
      tmp_3293 += (ZD(gt3,j2)) * tmp_3294;
   }
   tmp_3292 += tmp_3293;
   tmp_3251 += (1.4142135623730951*vu*Lambdax*ZH(gt1,2)) * tmp_3292;
   tmp_3246 += (3) * tmp_3251;
   result += (0.08333333333333333) * tmp_3246;

   return result;
}

std::complex<double> TMSSM_effective_couplings::CphhSuconjSu(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYu = MODELPARAMETER(TYu);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vT = MODELPARAMETER(vT);
   const auto vu = MODELPARAMETER(vu);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3296;
   std::complex<double> tmp_3297;
   std::complex<double> tmp_3298;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3298 += Conj(ZU(gt2,j1))*ZU(gt3,j1);
   }
   tmp_3297 += tmp_3298;
   tmp_3296 += ((0.6*Sqr(g1) - 3*Sqr(g2))*(vd*ZH(gt1,0) - vu*ZH(gt1,1))) *
      tmp_3297;
   std::complex<double> tmp_3299;
   std::complex<double> tmp_3300;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3300 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_3299 += tmp_3300;
   tmp_3296 += (2.4*Sqr(g1)*(-(vd*ZH(gt1,0)) + vu*ZH(gt1,1))) * tmp_3299;
   std::complex<double> tmp_3301;
   std::complex<double> tmp_3302;
   std::complex<double> tmp_3303;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3304;
      std::complex<double> tmp_3305;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3305 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_3304 += tmp_3305;
      tmp_3303 += (Conj(ZU(gt2,j2))) * tmp_3304;
   }
   tmp_3302 += tmp_3303;
   tmp_3301 += (2.8284271247461903*Conj(Mu)*ZH(gt1,0)) * tmp_3302;
   std::complex<double> tmp_3306;
   std::complex<double> tmp_3307;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3308;
      std::complex<double> tmp_3309;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3309 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3308 += tmp_3309;
      tmp_3307 += (ZU(gt3,j2)) * tmp_3308;
   }
   tmp_3306 += tmp_3307;
   tmp_3301 += (1.4142135623730951*vT*Lambdax*ZH(gt1,0)) * tmp_3306;
   std::complex<double> tmp_3310;
   std::complex<double> tmp_3311;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3312;
      std::complex<double> tmp_3313;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3313 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3312 += tmp_3313;
      tmp_3311 += (ZU(gt3,j2)) * tmp_3312;
   }
   tmp_3310 += tmp_3311;
   tmp_3301 += (2.8284271247461903*Mu*ZH(gt1,0)) * tmp_3310;
   std::complex<double> tmp_3314;
   std::complex<double> tmp_3315;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3316;
      std::complex<double> tmp_3317;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3317 += ZU(gt3,3 + j1)*TYu(j1,j2);
      }
      tmp_3316 += tmp_3317;
      tmp_3315 += (Conj(ZU(gt2,j2))) * tmp_3316;
   }
   tmp_3314 += tmp_3315;
   tmp_3301 += (-2.8284271247461903*ZH(gt1,1)) * tmp_3314;
   std::complex<double> tmp_3318;
   std::complex<double> tmp_3319;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3320;
      std::complex<double> tmp_3321;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3321 += Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_3320 += tmp_3321;
      tmp_3319 += (ZU(gt3,j2)) * tmp_3320;
   }
   tmp_3318 += tmp_3319;
   tmp_3301 += (-2.8284271247461903*ZH(gt1,1)) * tmp_3318;
   std::complex<double> tmp_3322;
   std::complex<double> tmp_3323;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3324;
      std::complex<double> tmp_3325;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3326;
         std::complex<double> tmp_3327;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3327 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_3326 += tmp_3327;
         tmp_3325 += (ZU(gt3,3 + j2)) * tmp_3326;
      }
      tmp_3324 += tmp_3325;
      tmp_3323 += (Conj(ZU(gt2,3 + j3))) * tmp_3324;
   }
   tmp_3322 += tmp_3323;
   tmp_3301 += (-4*vu*ZH(gt1,1)) * tmp_3322;
   std::complex<double> tmp_3328;
   std::complex<double> tmp_3329;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3330;
      std::complex<double> tmp_3331;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3332;
         std::complex<double> tmp_3333;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3333 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_3332 += tmp_3333;
         tmp_3331 += (Conj(ZU(gt2,j2))) * tmp_3332;
      }
      tmp_3330 += tmp_3331;
      tmp_3329 += (ZU(gt3,j3)) * tmp_3330;
   }
   tmp_3328 += tmp_3329;
   tmp_3301 += (-4*vu*ZH(gt1,1)) * tmp_3328;
   std::complex<double> tmp_3334;
   std::complex<double> tmp_3335;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3336;
      std::complex<double> tmp_3337;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3337 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3336 += tmp_3337;
      tmp_3335 += (ZU(gt3,j2)) * tmp_3336;
   }
   tmp_3334 += tmp_3335;
   tmp_3301 += (1.4142135623730951*vd*Lambdax*ZH(gt1,2)) * tmp_3334;
   std::complex<double> tmp_3338;
   std::complex<double> tmp_3339;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3340;
      std::complex<double> tmp_3341;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3341 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_3340 += tmp_3341;
      tmp_3339 += (Conj(ZU(gt2,j2))) * tmp_3340;
   }
   tmp_3338 += tmp_3339;
   tmp_3301 += (1.4142135623730951*Conj(Lambdax)*(vT*ZH(gt1,0) + vd*ZH(gt1,2)))
      * tmp_3338;
   tmp_3296 += (3) * tmp_3301;
   result += (0.08333333333333333) * tmp_3296;

   return result;
}

std::complex<double> TMSSM_effective_couplings::CphhSeconjSe(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYe = MODELPARAMETER(TYe);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vT = MODELPARAMETER(vT);
   const auto vu = MODELPARAMETER(vu);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3342;
   std::complex<double> tmp_3343;
   std::complex<double> tmp_3344;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3345;
      std::complex<double> tmp_3346;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3346 += ZE(gt3,3 + j1)*TYe(j1,j2);
      }
      tmp_3345 += tmp_3346;
      tmp_3344 += (Conj(ZE(gt2,j2))) * tmp_3345;
   }
   tmp_3343 += tmp_3344;
   tmp_3342 += (-2.8284271247461903*ZH(gt1,0)) * tmp_3343;
   std::complex<double> tmp_3347;
   std::complex<double> tmp_3348;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3349;
      std::complex<double> tmp_3350;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3350 += Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_3349 += tmp_3350;
      tmp_3348 += (ZE(gt3,j2)) * tmp_3349;
   }
   tmp_3347 += tmp_3348;
   tmp_3342 += (-2.8284271247461903*ZH(gt1,0)) * tmp_3347;
   std::complex<double> tmp_3351;
   std::complex<double> tmp_3352;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3353;
      std::complex<double> tmp_3354;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3355;
         std::complex<double> tmp_3356;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3356 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_3355 += tmp_3356;
         tmp_3354 += (ZE(gt3,3 + j2)) * tmp_3355;
      }
      tmp_3353 += tmp_3354;
      tmp_3352 += (Conj(ZE(gt2,3 + j3))) * tmp_3353;
   }
   tmp_3351 += tmp_3352;
   tmp_3342 += (-4*vd*ZH(gt1,0)) * tmp_3351;
   std::complex<double> tmp_3357;
   std::complex<double> tmp_3358;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3359;
      std::complex<double> tmp_3360;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3361;
         std::complex<double> tmp_3362;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3362 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_3361 += tmp_3362;
         tmp_3360 += (Conj(ZE(gt2,j2))) * tmp_3361;
      }
      tmp_3359 += tmp_3360;
      tmp_3358 += (ZE(gt3,j3)) * tmp_3359;
   }
   tmp_3357 += tmp_3358;
   tmp_3342 += (-4*vd*ZH(gt1,0)) * tmp_3357;
   std::complex<double> tmp_3363;
   std::complex<double> tmp_3364;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3365;
      std::complex<double> tmp_3366;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3366 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3365 += tmp_3366;
      tmp_3364 += (Conj(ZE(gt2,j2))) * tmp_3365;
   }
   tmp_3363 += tmp_3364;
   tmp_3342 += (1.4142135623730951*vT*Conj(Lambdax)*ZH(gt1,1)) * tmp_3363;
   std::complex<double> tmp_3367;
   std::complex<double> tmp_3368;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3369;
      std::complex<double> tmp_3370;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3370 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3369 += tmp_3370;
      tmp_3368 += (Conj(ZE(gt2,j2))) * tmp_3369;
   }
   tmp_3367 += tmp_3368;
   tmp_3342 += (2.8284271247461903*Conj(Mu)*ZH(gt1,1)) * tmp_3367;
   std::complex<double> tmp_3371;
   std::complex<double> tmp_3372;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3373;
      std::complex<double> tmp_3374;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3374 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3373 += tmp_3374;
      tmp_3372 += (ZE(gt3,j2)) * tmp_3373;
   }
   tmp_3371 += tmp_3372;
   tmp_3342 += (1.4142135623730951*vT*Lambdax*ZH(gt1,1)) * tmp_3371;
   std::complex<double> tmp_3375;
   std::complex<double> tmp_3376;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3377;
      std::complex<double> tmp_3378;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3378 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3377 += tmp_3378;
      tmp_3376 += (ZE(gt3,j2)) * tmp_3377;
   }
   tmp_3375 += tmp_3376;
   tmp_3342 += (2.8284271247461903*Mu*ZH(gt1,1)) * tmp_3375;
   std::complex<double> tmp_3379;
   std::complex<double> tmp_3380;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3380 += Conj(ZE(gt2,j1))*ZE(gt3,j1);
   }
   tmp_3379 += tmp_3380;
   tmp_3342 += (-((0.6*Sqr(g1) - Sqr(g2))*(vd*ZH(gt1,0) - vu*ZH(gt1,1)))) *
      tmp_3379;
   std::complex<double> tmp_3381;
   std::complex<double> tmp_3382;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3382 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_3381 += tmp_3382;
   tmp_3342 += (1.2*Sqr(g1)*(vd*ZH(gt1,0) - vu*ZH(gt1,1))) * tmp_3381;
   std::complex<double> tmp_3383;
   std::complex<double> tmp_3384;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3385;
      std::complex<double> tmp_3386;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3386 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3385 += tmp_3386;
      tmp_3384 += (Conj(ZE(gt2,j2))) * tmp_3385;
   }
   tmp_3383 += tmp_3384;
   tmp_3342 += (1.4142135623730951*vu*Conj(Lambdax)*ZH(gt1,2)) * tmp_3383;
   std::complex<double> tmp_3387;
   std::complex<double> tmp_3388;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3389;
      std::complex<double> tmp_3390;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3390 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3389 += tmp_3390;
      tmp_3388 += (ZE(gt3,j2)) * tmp_3389;
   }
   tmp_3387 += tmp_3388;
   tmp_3342 += (1.4142135623730951*vu*Lambdax*ZH(gt1,2)) * tmp_3387;
   result += (0.25) * tmp_3342;

   return result;
}

std::complex<double> TMSSM_effective_couplings::CphhHpmconjHpm(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TLambdax = MODELPARAMETER(TLambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto MT = MODELPARAMETER(MT);
   const auto vd = MODELPARAMETER(vd);
   const auto vT = MODELPARAMETER(vT);
   const auto vu = MODELPARAMETER(vu);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   result = 0.05*(-5*ZH(gt1,2)*(-2*Conj(Mu)*Lambdax*ZP(gt2,0)*ZP(gt3,0) - 2*
      Conj(TLambdax)*ZP(gt2,1)*ZP(gt3,0) + 1.4142135623730951*vd*Sqr(g2)*ZP(gt2,2)
      *ZP(gt3,0) - 1.4142135623730951*vd*Sqr(g2)*ZP(gt2,3)*ZP(gt3,0) - 4*Conj(MT)*
      Lambdax*ZP(gt2,0)*ZP(gt3,1) - 2*TLambdax*ZP(gt2,0)*ZP(gt3,1) - 2*Conj(Mu)*
      Lambdax*ZP(gt2,1)*ZP(gt3,1) + 1.4142135623730951*vu*Sqr(g2)*ZP(gt2,2)*ZP(gt3
      ,1) - 1.4142135623730951*vu*Sqr(g2)*ZP(gt2,3)*ZP(gt3,1) + 1.4142135623730951
      *vd*Sqr(g2)*ZP(gt2,0)*ZP(gt3,2) + 1.4142135623730951*vu*Sqr(g2)*ZP(gt2,1)*ZP
      (gt3,2) + 4*vT*Sqr(g2)*ZP(gt2,2)*ZP(gt3,2) - 4*vT*Sqr(g2)*ZP(gt2,3)*ZP(gt3,2
      ) - 1.4142135623730951*vd*Sqr(g2)*ZP(gt2,0)*ZP(gt3,3) - 1.4142135623730951*
      vu*Sqr(g2)*ZP(gt2,1)*ZP(gt3,3) - 4*vT*Sqr(g2)*ZP(gt2,2)*ZP(gt3,3) + 4*vT*Sqr
      (g2)*ZP(gt2,3)*ZP(gt3,3) + Conj(Lambdax)*(-1.4142135623730951*Lambdax*(ZP(
      gt2,2) - ZP(gt2,3))*(vd*ZP(gt3,0) + vu*ZP(gt3,1)) + ZP(gt2,0)*(2*(vT*Lambdax
      - Mu)*ZP(gt3,0) + 1.4142135623730951*vd*Lambdax*(-ZP(gt3,2) + ZP(gt3,3))) +
      ZP(gt2,1)*(-4*MT*ZP(gt3,0) + 2*(vT*Lambdax - Mu)*ZP(gt3,1) +
      1.4142135623730951*vu*Lambdax*(-ZP(gt3,2) + ZP(gt3,3))))) - ZH(gt1,0)*(ZP(
      gt2,1)*(5*vu*(AbsSqr(Lambdax) + Sqr(g2))*ZP(gt3,0) + vd*(20*AbsSqr(Lambdax)
      - 3*Sqr(g1) + 5*Sqr(g2))*ZP(gt3,1) + 14.142135623730951*(Conj(TLambdax)*ZP(
      gt3,2) + 2*MT*Conj(Lambdax)*ZP(gt3,3))) + 5*(ZP(gt2,2)*(1.4142135623730951*(
      -(vT*AbsSqr(Lambdax)) + 2*Conj(Mu)*Lambdax + vT*Sqr(g2))*ZP(gt3,0) +
      2.8284271247461903*TLambdax*ZP(gt3,1) - 2*vd*(-2*AbsSqr(Lambdax) + Sqr(g2))*
      ZP(gt3,2)) + ZP(gt2,3)*(1.4142135623730951*(Conj(Lambdax)*(vT*Lambdax + 2*Mu
      ) - vT*Sqr(g2))*ZP(gt3,0) + 5.656854249492381*Conj(MT)*Lambdax*ZP(gt3,1) + 2
      *vd*Sqr(g2)*ZP(gt3,3))) + ZP(gt2,0)*(vd*(3*Sqr(g1) + 5*Sqr(g2))*ZP(gt3,0) +
      5*(vu*(AbsSqr(Lambdax) + Sqr(g2))*ZP(gt3,1) + 1.4142135623730951*((-(vT*
      AbsSqr(Lambdax)) + 2*Conj(Lambdax)*Mu + vT*Sqr(g2))*ZP(gt3,2) + (vT*AbsSqr(
      Lambdax) + 2*Conj(Mu)*Lambdax - vT*Sqr(g2))*ZP(gt3,3))))) + ZH(gt1,1)*(-(ZP(
      gt2,1)*(5*vd*(AbsSqr(Lambdax) + Sqr(g2))*ZP(gt3,0) + vu*(3*Sqr(g1) + 5*Sqr(
      g2))*ZP(gt3,1) + 7.0710678118654755*((-(Conj(Lambdax)*(vT*Lambdax + 2*Mu)) +
      vT*Sqr(g2))*ZP(gt3,2) + (vT*AbsSqr(Lambdax) - 2*Conj(Mu)*Lambdax - vT*Sqr(
      g2))*ZP(gt3,3)))) + 5*(2.8284271247461903*Conj(TLambdax)*ZP(gt2,3)*ZP(gt3,0)
      + 2.8284271247461903*Conj(Mu)*Lambdax*ZP(gt2,2)*ZP(gt3,1) -
      1.4142135623730951*vT*Sqr(g2)*ZP(gt2,2)*ZP(gt3,1) + 1.4142135623730951*vT*
      Sqr(g2)*ZP(gt2,3)*ZP(gt3,1) - 2*vu*Sqr(g2)*ZP(gt2,2)*ZP(gt3,2) + 2*vu*Sqr(g2
      )*ZP(gt2,3)*ZP(gt3,3) + Conj(Lambdax)*(1.4142135623730951*ZP(gt2,2)*(4*MT*ZP
      (gt3,0) + vT*Lambdax*ZP(gt3,1)) - ZP(gt2,3)*(1.4142135623730951*(vT*Lambdax
      - 2*Mu)*ZP(gt3,1) + 4*vu*Lambdax*ZP(gt3,3)))) + ZP(gt2,0)*(vu*(-20*AbsSqr(
      Lambdax) + 3*Sqr(g1) - 5*Sqr(g2))*ZP(gt3,0) - 5*(vd*(AbsSqr(Lambdax) + Sqr(
      g2))*ZP(gt3,1) - 2.8284271247461903*(2*Conj(MT)*Lambdax*ZP(gt3,2) + TLambdax
      *ZP(gt3,3))))));

   return result;
}

std::complex<double> TMSSM_effective_couplings::CpChahhbarChaPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   std::complex<double> result;

   result = 0.5*(-(Conj(UM(gt1,2))*(1.4142135623730951*Conj(UP(gt3,1))*Lambdax*
      ZH(gt2,0) + 2*g2*Conj(UP(gt3,0))*ZH(gt2,2))) + g2*Conj(UM(gt1,0))*(
      -1.4142135623730951*Conj(UP(gt3,1))*ZH(gt2,1) + 2*Conj(UP(gt3,2))*ZH(gt2,2))
      + Conj(UM(gt1,1))*(-1.4142135623730951*g2*Conj(UP(gt3,0))*ZH(gt2,0) +
      1.4142135623730951*Conj(UP(gt3,2))*Lambdax*ZH(gt2,1) + Conj(UP(gt3,1))*
      Lambdax*ZH(gt2,2)));

   return result;
}

std::complex<double> TMSSM_effective_couplings::CpFehhbarFePL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Ye = MODELPARAMETER(Ye);

   std::complex<double> result;

   std::complex<double> tmp_3391;
   std::complex<double> tmp_3392;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3393;
      std::complex<double> tmp_3394;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3394 += Conj(ZER(gt3,j1))*Ye(j1,j2);
      }
      tmp_3393 += tmp_3394;
      tmp_3392 += (Conj(ZEL(gt1,j2))) * tmp_3393;
   }
   tmp_3391 += tmp_3392;
   result += (-0.7071067811865475*ZH(gt2,0)) * tmp_3391;

   return result;
}

std::complex<double> TMSSM_effective_couplings::CpFdhhbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_3395;
   std::complex<double> tmp_3396;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3397;
      std::complex<double> tmp_3398;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3398 += Conj(ZDR(gt3,j1))*Yd(j1,j2);
      }
      tmp_3397 += tmp_3398;
      tmp_3396 += (Conj(ZDL(gt1,j2))) * tmp_3397;
   }
   tmp_3395 += tmp_3396;
   result += (-0.7071067811865475*ZH(gt2,0)) * tmp_3395;

   return result;
}

std::complex<double> TMSSM_effective_couplings::CpFuhhbarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_3399;
   std::complex<double> tmp_3400;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3401;
      std::complex<double> tmp_3402;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3402 += Conj(ZUR(gt3,j1))*Yu(j1,j2);
      }
      tmp_3401 += tmp_3402;
      tmp_3400 += (Conj(ZUL(gt1,j2))) * tmp_3401;
   }
   tmp_3399 += tmp_3400;
   result += (-0.7071067811865475*ZH(gt2,1)) * tmp_3399;

   return result;
}

std::complex<double> TMSSM_effective_couplings::CphhVWmconjVWm(unsigned gt1) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vT = MODELPARAMETER(vT);
   const auto vu = MODELPARAMETER(vu);

   std::complex<double> result;

   result = 0.5*Sqr(g2)*(vd*ZH(gt1,0) + vu*ZH(gt1,1) + 4*vT*ZH(gt1,2));

   return result;
}

std::complex<double> TMSSM_effective_couplings::CpAhSdconjSd(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYd = MODELPARAMETER(TYd);
   const auto vT = MODELPARAMETER(vT);
   const auto vu = MODELPARAMETER(vu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3403;
   std::complex<double> tmp_3404;
   std::complex<double> tmp_3405;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3406;
      std::complex<double> tmp_3407;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3407 += ZD(gt3,3 + j1)*TYd(j1,j2);
      }
      tmp_3406 += tmp_3407;
      tmp_3405 += (Conj(ZD(gt2,j2))) * tmp_3406;
   }
   tmp_3404 += tmp_3405;
   tmp_3403 += (2*ZA(gt1,0)) * tmp_3404;
   std::complex<double> tmp_3408;
   std::complex<double> tmp_3409;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3410;
      std::complex<double> tmp_3411;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3411 += Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_3410 += tmp_3411;
      tmp_3409 += (ZD(gt3,j2)) * tmp_3410;
   }
   tmp_3408 += tmp_3409;
   tmp_3403 += (-2*ZA(gt1,0)) * tmp_3408;
   std::complex<double> tmp_3412;
   std::complex<double> tmp_3413;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3414;
      std::complex<double> tmp_3415;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3415 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3414 += tmp_3415;
      tmp_3413 += (Conj(ZD(gt2,j2))) * tmp_3414;
   }
   tmp_3412 += tmp_3413;
   tmp_3403 += (vT*Conj(Lambdax)*ZA(gt1,1)) * tmp_3412;
   std::complex<double> tmp_3416;
   std::complex<double> tmp_3417;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3418;
      std::complex<double> tmp_3419;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3419 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3418 += tmp_3419;
      tmp_3417 += (Conj(ZD(gt2,j2))) * tmp_3418;
   }
   tmp_3416 += tmp_3417;
   tmp_3403 += (2*Conj(Mu)*ZA(gt1,1)) * tmp_3416;
   std::complex<double> tmp_3420;
   std::complex<double> tmp_3421;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3422;
      std::complex<double> tmp_3423;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3423 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3422 += tmp_3423;
      tmp_3421 += (ZD(gt3,j2)) * tmp_3422;
   }
   tmp_3420 += tmp_3421;
   tmp_3403 += (-(vT*Lambdax*ZA(gt1,1))) * tmp_3420;
   std::complex<double> tmp_3424;
   std::complex<double> tmp_3425;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3426;
      std::complex<double> tmp_3427;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3427 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3426 += tmp_3427;
      tmp_3425 += (ZD(gt3,j2)) * tmp_3426;
   }
   tmp_3424 += tmp_3425;
   tmp_3403 += (-2*Mu*ZA(gt1,1)) * tmp_3424;
   std::complex<double> tmp_3428;
   std::complex<double> tmp_3429;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3430;
      std::complex<double> tmp_3431;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3431 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3430 += tmp_3431;
      tmp_3429 += (Conj(ZD(gt2,j2))) * tmp_3430;
   }
   tmp_3428 += tmp_3429;
   tmp_3403 += (vu*Conj(Lambdax)*ZA(gt1,2)) * tmp_3428;
   std::complex<double> tmp_3432;
   std::complex<double> tmp_3433;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3434;
      std::complex<double> tmp_3435;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3435 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3434 += tmp_3435;
      tmp_3433 += (ZD(gt3,j2)) * tmp_3434;
   }
   tmp_3432 += tmp_3433;
   tmp_3403 += (-(vu*Lambdax*ZA(gt1,2))) * tmp_3432;
   result += (std::complex<double>(0.,-0.35355339059327373)) * tmp_3403;

   return result;
}

std::complex<double> TMSSM_effective_couplings::CpAhSuconjSu(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYu = MODELPARAMETER(TYu);
   const auto vd = MODELPARAMETER(vd);
   const auto vT = MODELPARAMETER(vT);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3436;
   std::complex<double> tmp_3437;
   std::complex<double> tmp_3438;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3439;
      std::complex<double> tmp_3440;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3440 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_3439 += tmp_3440;
      tmp_3438 += (Conj(ZU(gt2,j2))) * tmp_3439;
   }
   tmp_3437 += tmp_3438;
   tmp_3436 += (2*Conj(Mu)*ZA(gt1,0)) * tmp_3437;
   std::complex<double> tmp_3441;
   std::complex<double> tmp_3442;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3443;
      std::complex<double> tmp_3444;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3444 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3443 += tmp_3444;
      tmp_3442 += (ZU(gt3,j2)) * tmp_3443;
   }
   tmp_3441 += tmp_3442;
   tmp_3436 += (-(vT*Lambdax*ZA(gt1,0))) * tmp_3441;
   std::complex<double> tmp_3445;
   std::complex<double> tmp_3446;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3447;
      std::complex<double> tmp_3448;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3448 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3447 += tmp_3448;
      tmp_3446 += (ZU(gt3,j2)) * tmp_3447;
   }
   tmp_3445 += tmp_3446;
   tmp_3436 += (-2*Mu*ZA(gt1,0)) * tmp_3445;
   std::complex<double> tmp_3449;
   std::complex<double> tmp_3450;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3451;
      std::complex<double> tmp_3452;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3452 += ZU(gt3,3 + j1)*TYu(j1,j2);
      }
      tmp_3451 += tmp_3452;
      tmp_3450 += (Conj(ZU(gt2,j2))) * tmp_3451;
   }
   tmp_3449 += tmp_3450;
   tmp_3436 += (2*ZA(gt1,1)) * tmp_3449;
   std::complex<double> tmp_3453;
   std::complex<double> tmp_3454;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3455;
      std::complex<double> tmp_3456;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3456 += Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_3455 += tmp_3456;
      tmp_3454 += (ZU(gt3,j2)) * tmp_3455;
   }
   tmp_3453 += tmp_3454;
   tmp_3436 += (-2*ZA(gt1,1)) * tmp_3453;
   std::complex<double> tmp_3457;
   std::complex<double> tmp_3458;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3459;
      std::complex<double> tmp_3460;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3460 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3459 += tmp_3460;
      tmp_3458 += (ZU(gt3,j2)) * tmp_3459;
   }
   tmp_3457 += tmp_3458;
   tmp_3436 += (-(vd*Lambdax*ZA(gt1,2))) * tmp_3457;
   std::complex<double> tmp_3461;
   std::complex<double> tmp_3462;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3463;
      std::complex<double> tmp_3464;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3464 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_3463 += tmp_3464;
      tmp_3462 += (Conj(ZU(gt2,j2))) * tmp_3463;
   }
   tmp_3461 += tmp_3462;
   tmp_3436 += (Conj(Lambdax)*(vT*ZA(gt1,0) + vd*ZA(gt1,2))) * tmp_3461;
   result += (std::complex<double>(0.,-0.35355339059327373)) * tmp_3436;

   return result;
}

std::complex<double> TMSSM_effective_couplings::CpAhSeconjSe(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYe = MODELPARAMETER(TYe);
   const auto vT = MODELPARAMETER(vT);
   const auto vu = MODELPARAMETER(vu);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3465;
   std::complex<double> tmp_3466;
   std::complex<double> tmp_3467;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3468;
      std::complex<double> tmp_3469;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3469 += ZE(gt3,3 + j1)*TYe(j1,j2);
      }
      tmp_3468 += tmp_3469;
      tmp_3467 += (Conj(ZE(gt2,j2))) * tmp_3468;
   }
   tmp_3466 += tmp_3467;
   tmp_3465 += (2*ZA(gt1,0)) * tmp_3466;
   std::complex<double> tmp_3470;
   std::complex<double> tmp_3471;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3472;
      std::complex<double> tmp_3473;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3473 += Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_3472 += tmp_3473;
      tmp_3471 += (ZE(gt3,j2)) * tmp_3472;
   }
   tmp_3470 += tmp_3471;
   tmp_3465 += (-2*ZA(gt1,0)) * tmp_3470;
   std::complex<double> tmp_3474;
   std::complex<double> tmp_3475;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3476;
      std::complex<double> tmp_3477;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3477 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3476 += tmp_3477;
      tmp_3475 += (Conj(ZE(gt2,j2))) * tmp_3476;
   }
   tmp_3474 += tmp_3475;
   tmp_3465 += (vT*Conj(Lambdax)*ZA(gt1,1)) * tmp_3474;
   std::complex<double> tmp_3478;
   std::complex<double> tmp_3479;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3480;
      std::complex<double> tmp_3481;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3481 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3480 += tmp_3481;
      tmp_3479 += (Conj(ZE(gt2,j2))) * tmp_3480;
   }
   tmp_3478 += tmp_3479;
   tmp_3465 += (2*Conj(Mu)*ZA(gt1,1)) * tmp_3478;
   std::complex<double> tmp_3482;
   std::complex<double> tmp_3483;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3484;
      std::complex<double> tmp_3485;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3485 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3484 += tmp_3485;
      tmp_3483 += (ZE(gt3,j2)) * tmp_3484;
   }
   tmp_3482 += tmp_3483;
   tmp_3465 += (-(vT*Lambdax*ZA(gt1,1))) * tmp_3482;
   std::complex<double> tmp_3486;
   std::complex<double> tmp_3487;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3488;
      std::complex<double> tmp_3489;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3489 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3488 += tmp_3489;
      tmp_3487 += (ZE(gt3,j2)) * tmp_3488;
   }
   tmp_3486 += tmp_3487;
   tmp_3465 += (-2*Mu*ZA(gt1,1)) * tmp_3486;
   std::complex<double> tmp_3490;
   std::complex<double> tmp_3491;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3492;
      std::complex<double> tmp_3493;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3493 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3492 += tmp_3493;
      tmp_3491 += (Conj(ZE(gt2,j2))) * tmp_3492;
   }
   tmp_3490 += tmp_3491;
   tmp_3465 += (vu*Conj(Lambdax)*ZA(gt1,2)) * tmp_3490;
   std::complex<double> tmp_3494;
   std::complex<double> tmp_3495;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3496;
      std::complex<double> tmp_3497;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3497 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3496 += tmp_3497;
      tmp_3495 += (ZE(gt3,j2)) * tmp_3496;
   }
   tmp_3494 += tmp_3495;
   tmp_3465 += (-(vu*Lambdax*ZA(gt1,2))) * tmp_3494;
   result += (std::complex<double>(0.,-0.35355339059327373)) * tmp_3465;

   return result;
}

std::complex<double> TMSSM_effective_couplings::CpAhHpmconjHpm(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TLambdax = MODELPARAMETER(TLambdax);
   const auto g2 = MODELPARAMETER(g2);
   const auto MT = MODELPARAMETER(MT);
   const auto vd = MODELPARAMETER(vd);
   const auto vT = MODELPARAMETER(vT);
   const auto vu = MODELPARAMETER(vu);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   result = std::complex<double>(0,-0.25)*(vu*Sqr(g2)*ZA(gt1,0)*ZP(gt2,1)*ZP(
      gt3,0) + vd*Sqr(g2)*ZA(gt1,1)*ZP(gt2,1)*ZP(gt3,0) + 2*Conj(TLambdax)*ZA(gt1,
      2)*ZP(gt2,1)*ZP(gt3,0) + 1.4142135623730951*vT*Sqr(g2)*ZA(gt1,0)*ZP(gt2,2)*
      ZP(gt3,0) - 1.4142135623730951*vd*Sqr(g2)*ZA(gt1,2)*ZP(gt2,2)*ZP(gt3,0) -
      1.4142135623730951*vT*Sqr(g2)*ZA(gt1,0)*ZP(gt2,3)*ZP(gt3,0) +
      2.8284271247461903*Conj(TLambdax)*ZA(gt1,1)*ZP(gt2,3)*ZP(gt3,0) -
      1.4142135623730951*vd*Sqr(g2)*ZA(gt1,2)*ZP(gt2,3)*ZP(gt3,0) - vu*Sqr(g2)*ZA(
      gt1,0)*ZP(gt2,0)*ZP(gt3,1) - vd*Sqr(g2)*ZA(gt1,1)*ZP(gt2,0)*ZP(gt3,1) + 4*
      Conj(MT)*Lambdax*ZA(gt1,2)*ZP(gt2,0)*ZP(gt3,1) - 2*TLambdax*ZA(gt1,2)*ZP(gt2
      ,0)*ZP(gt3,1) + 2.8284271247461903*TLambdax*ZA(gt1,0)*ZP(gt2,2)*ZP(gt3,1) -
      1.4142135623730951*vT*Sqr(g2)*ZA(gt1,1)*ZP(gt2,2)*ZP(gt3,1) -
      1.4142135623730951*vu*Sqr(g2)*ZA(gt1,2)*ZP(gt2,2)*ZP(gt3,1) +
      5.656854249492381*Conj(MT)*Lambdax*ZA(gt1,0)*ZP(gt2,3)*ZP(gt3,1) +
      1.4142135623730951*vT*Sqr(g2)*ZA(gt1,1)*ZP(gt2,3)*ZP(gt3,1) -
      1.4142135623730951*vu*Sqr(g2)*ZA(gt1,2)*ZP(gt2,3)*ZP(gt3,1) -
      1.4142135623730951*vT*Sqr(g2)*ZA(gt1,0)*ZP(gt2,0)*ZP(gt3,2) -
      5.656854249492381*Conj(MT)*Lambdax*ZA(gt1,1)*ZP(gt2,0)*ZP(gt3,2) +
      1.4142135623730951*vd*Sqr(g2)*ZA(gt1,2)*ZP(gt2,0)*ZP(gt3,2) -
      2.8284271247461903*Conj(TLambdax)*ZA(gt1,0)*ZP(gt2,1)*ZP(gt3,2) +
      1.4142135623730951*vT*Sqr(g2)*ZA(gt1,1)*ZP(gt2,1)*ZP(gt3,2) +
      1.4142135623730951*vu*Sqr(g2)*ZA(gt1,2)*ZP(gt2,1)*ZP(gt3,2) - 4*vT*Sqr(g2)*
      ZA(gt1,2)*ZP(gt2,3)*ZP(gt3,2) + 1.4142135623730951*vT*Sqr(g2)*ZA(gt1,0)*ZP(
      gt2,0)*ZP(gt3,3) - 2.8284271247461903*TLambdax*ZA(gt1,1)*ZP(gt2,0)*ZP(gt3,3)
      + 1.4142135623730951*vd*Sqr(g2)*ZA(gt1,2)*ZP(gt2,0)*ZP(gt3,3) -
      1.4142135623730951*vT*Sqr(g2)*ZA(gt1,1)*ZP(gt2,1)*ZP(gt3,3) +
      1.4142135623730951*vu*Sqr(g2)*ZA(gt1,2)*ZP(gt2,1)*ZP(gt3,3) + 4*vT*Sqr(g2)*
      ZA(gt1,2)*ZP(gt2,2)*ZP(gt3,3) - 2*Conj(Mu)*Lambdax*(ZA(gt1,2)*(ZP(gt2,0)*ZP(
      gt3,0) + ZP(gt2,1)*ZP(gt3,1)) + 1.4142135623730951*(ZA(gt1,0)*(-(ZP(gt2,2)*
      ZP(gt3,0)) + ZP(gt2,0)*ZP(gt3,3)) + ZA(gt1,1)*(-(ZP(gt2,2)*ZP(gt3,1)) + ZP(
      gt2,1)*ZP(gt3,3)))) + Conj(Lambdax)*(ZA(gt1,0)*(-1.4142135623730951*vT*
      Lambdax*ZP(gt2,2)*ZP(gt3,0) + 1.4142135623730951*vT*Lambdax*ZP(gt2,3)*ZP(gt3
      ,0) + 2.8284271247461903*Mu*ZP(gt2,3)*ZP(gt3,0) - vu*Lambdax*ZP(gt2,0)*ZP(
      gt3,1) + 1.4142135623730951*vT*Lambdax*ZP(gt2,0)*ZP(gt3,2) -
      2.8284271247461903*Mu*ZP(gt2,0)*ZP(gt3,2) - 1.4142135623730951*vT*Lambdax*ZP
      (gt2,0)*ZP(gt3,3) + ZP(gt2,1)*(vu*Lambdax*ZP(gt3,0) - 5.656854249492381*MT*
      ZP(gt3,3))) + ZA(gt1,2)*(1.4142135623730951*Lambdax*(ZP(gt2,2) + ZP(gt2,3))*
      (vd*ZP(gt3,0) + vu*ZP(gt3,1)) + ZP(gt2,0)*(2*Mu*ZP(gt3,0) -
      1.4142135623730951*vd*Lambdax*(ZP(gt3,2) + ZP(gt3,3))) - ZP(gt2,1)*(4*MT*ZP(
      gt3,0) - 2*Mu*ZP(gt3,1) + 1.4142135623730951*vu*Lambdax*(ZP(gt3,2) + ZP(gt3,
      3)))) + ZA(gt1,1)*(-((vd*Lambdax*ZP(gt2,0) + 1.4142135623730951*(vT*Lambdax
      - 2*Mu)*ZP(gt2,3))*ZP(gt3,1)) + 1.4142135623730951*ZP(gt2,2)*(4*MT*ZP(gt3,0)
      + vT*Lambdax*ZP(gt3,1)) + ZP(gt2,1)*(vd*Lambdax*ZP(gt3,0) +
      1.4142135623730951*(-((vT*Lambdax + 2*Mu)*ZP(gt3,2)) + vT*Lambdax*ZP(gt3,3))
      ))));

   return result;
}

std::complex<double> TMSSM_effective_couplings::CpAhChabarChaPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   std::complex<double> result;

   result = std::complex<double>(0,0.5)*(Conj(UM(gt2,2))*(-1.4142135623730951*
      Conj(UP(gt3,1))*Lambdax*ZA(gt1,0) + 2*g2*Conj(UP(gt3,0))*ZA(gt1,2)) + g2*
      Conj(UM(gt2,0))*(1.4142135623730951*Conj(UP(gt3,1))*ZA(gt1,1) - 2*Conj(UP(
      gt3,2))*ZA(gt1,2)) + Conj(UM(gt2,1))*(1.4142135623730951*g2*Conj(UP(gt3,0))*
      ZA(gt1,0) + 1.4142135623730951*Conj(UP(gt3,2))*Lambdax*ZA(gt1,1) + Conj(UP(
      gt3,1))*Lambdax*ZA(gt1,2)));

   return result;
}

std::complex<double> TMSSM_effective_couplings::CpAhFebarFePL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Ye = MODELPARAMETER(Ye);

   std::complex<double> result;

   std::complex<double> tmp_3498;
   std::complex<double> tmp_3499;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3500;
      std::complex<double> tmp_3501;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3501 += Conj(ZER(gt3,j1))*Ye(j1,j2);
      }
      tmp_3500 += tmp_3501;
      tmp_3499 += (Conj(ZEL(gt2,j2))) * tmp_3500;
   }
   tmp_3498 += tmp_3499;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,0)) *
      tmp_3498;

   return result;
}

std::complex<double> TMSSM_effective_couplings::CpAhFdbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_3502;
   std::complex<double> tmp_3503;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3504;
      std::complex<double> tmp_3505;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3505 += Conj(ZDR(gt3,j1))*Yd(j1,j2);
      }
      tmp_3504 += tmp_3505;
      tmp_3503 += (Conj(ZDL(gt2,j2))) * tmp_3504;
   }
   tmp_3502 += tmp_3503;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,0)) *
      tmp_3502;

   return result;
}

std::complex<double> TMSSM_effective_couplings::CpAhFubarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_3506;
   std::complex<double> tmp_3507;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3508;
      std::complex<double> tmp_3509;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3509 += Conj(ZUR(gt3,j1))*Yu(j1,j2);
      }
      tmp_3508 += tmp_3509;
      tmp_3507 += (Conj(ZUL(gt2,j2))) * tmp_3508;
   }
   tmp_3506 += tmp_3507;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,1)) *
      tmp_3506;

   return result;
}

void TMSSM_effective_couplings::calculate_eff_CphhVPVP(unsigned gO1)
{
   const auto vd = MODELPARAMETER(vd);
   const auto vT = MODELPARAMETER(vT);
   const auto vu = MODELPARAMETER(vu);
   const auto MSd = MODELPARAMETER(MSd);
   const auto MSu = MODELPARAMETER(MSu);
   const auto MSe = MODELPARAMETER(MSe);
   const auto MHpm = MODELPARAMETER(MHpm);
   const auto MCha = MODELPARAMETER(MCha);
   const auto MFe = MODELPARAMETER(MFe);
   const auto MFd = MODELPARAMETER(MFd);
   const auto MFu = MODELPARAMETER(MFu);
   const auto MVWm = MODELPARAMETER(MVWm);
   const auto decay_mass = PHYSICAL(Mhh)(gO1);
   const auto decay_scale = 0.25 * Sqr(decay_mass);
   const auto saved_ZH = ZH;
   ZH = PHYSICAL(ZH);

   const auto vev = Sqrt(Sqr(vd) + 4*Sqr(vT) + Sqr(vu));

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
   for (unsigned gI1 = 1; gI1 < 4; ++gI1) {
      result += 0.5 * CphhHpmconjHpm(gO1, gI1, gI1) * vev * AS0(decay_scale
         / Sqr(MHpm(gI1))) / Sqr(MHpm(gI1));
   }
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
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
   result += -0.5 * CphhVWmconjVWm(gO1) * vev * AS1(decay_scale / Sqr(MVWm)) /
      Sqr(MVWm);

   result *= 0.1892681907127351 * physical_input.get(Physical_input::alpha_em_0
      ) * Sqrt(qedqcd.displayFermiConstant());

   ZH = saved_ZH;
   eff_CphhVPVP(gO1) = result;

}

void TMSSM_effective_couplings::calculate_eff_CphhVGVG(unsigned gO1)
{
   const auto vd = MODELPARAMETER(vd);
   const auto vT = MODELPARAMETER(vT);
   const auto vu = MODELPARAMETER(vu);
   const auto g3 = MODELPARAMETER(g3);
   const auto MSd = MODELPARAMETER(MSd);
   const auto MSu = MODELPARAMETER(MSu);
   const auto MFd = MODELPARAMETER(MFd);
   const auto MFu = MODELPARAMETER(MFu);
   const double alpha_s = 0.07957747154594767*Sqr(g3);
   const auto decay_mass = PHYSICAL(Mhh)(gO1);
   const auto decay_scale = 0.25 * Sqr(decay_mass);
   const auto saved_ZH = ZH;
   ZH = PHYSICAL(ZH);

   const auto vev = Sqrt(Sqr(vd) + 4*Sqr(vT) + Sqr(vu));

   std::complex<double> result = 0;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      result += 0.5 * CphhSdconjSd(gO1, gI1, gI1) * vev * AS0(decay_scale /
         Sqr(MSd(gI1))) / Sqr(MSd(gI1));
   }
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      result += 0.5 * CphhSuconjSu(gO1, gI1, gI1) * vev * AS0(decay_scale /
         Sqr(MSu(gI1))) / Sqr(MSu(gI1));
   }
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

void TMSSM_effective_couplings::calculate_eff_CpAhVPVP(unsigned gO1)
{
   const auto vd = MODELPARAMETER(vd);
   const auto vT = MODELPARAMETER(vT);
   const auto vu = MODELPARAMETER(vu);
   const auto MCha = MODELPARAMETER(MCha);
   const auto MFe = MODELPARAMETER(MFe);
   const auto MFd = MODELPARAMETER(MFd);
   const auto MFu = MODELPARAMETER(MFu);
   const auto decay_mass = PHYSICAL(MAh)(gO1);
   const auto decay_scale = 0.25 * Sqr(decay_mass);
   const auto saved_ZA = ZA;
   ZA = PHYSICAL(ZA);

   const auto vev = Sqrt(Sqr(vd) + 4*Sqr(vT) + Sqr(vu));

   std::complex<double> result = 0;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
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
   result *= std::complex<double>(2.0,0.);

   result *= 0.1892681907127351 * physical_input.get(Physical_input::alpha_em_0
      ) * Sqrt(qedqcd.displayFermiConstant());

   ZA = saved_ZA;
   eff_CpAhVPVP(gO1) = result;

}

void TMSSM_effective_couplings::calculate_eff_CpAhVGVG(unsigned gO1)
{
   const auto vd = MODELPARAMETER(vd);
   const auto vT = MODELPARAMETER(vT);
   const auto vu = MODELPARAMETER(vu);
   const auto g3 = MODELPARAMETER(g3);
   const auto MFd = MODELPARAMETER(MFd);
   const auto MFu = MODELPARAMETER(MFu);
   const double alpha_s = 0.07957747154594767*Sqr(g3);
   const auto decay_mass = PHYSICAL(MAh)(gO1);
   const auto decay_scale = 0.25 * Sqr(decay_mass);
   const auto saved_ZA = ZA;
   ZA = PHYSICAL(ZA);

   const auto vev = Sqrt(Sqr(vd) + 4*Sqr(vT) + Sqr(vu));

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
