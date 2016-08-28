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

// File generated at Sun 28 Aug 2016 15:15:44

#include "NUTSMSSM_effective_couplings.hpp"

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

NUTSMSSM_effective_couplings::NUTSMSSM_effective_couplings(
   const NUTSMSSM_mass_eigenstates& model_,
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

void NUTSMSSM_effective_couplings::calculate_effective_couplings()
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

void NUTSMSSM_effective_couplings::set_model(const NUTSMSSM_mass_eigenstates& model_)
{
   model = model_;
   copy_mixing_matrices_from_model();
}

void NUTSMSSM_effective_couplings::copy_mixing_matrices_from_model()
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

void NUTSMSSM_effective_couplings::run_SM_strong_coupling_to(double m)
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

std::complex<double> NUTSMSSM_effective_couplings::scalar_scalar_qcd_factor(double m_decay, double m_loop) const
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

std::complex<double> NUTSMSSM_effective_couplings::scalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) *
         scalar_diphoton_fermion_loop(m_decay, m_loop);

   }

   return result;
}

std::complex<double> NUTSMSSM_effective_couplings::pseudoscalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) *
         pseudoscalar_diphoton_fermion_loop(m_decay, m_loop);

   }

   return result;
}

double NUTSMSSM_effective_couplings::number_of_active_flavours(double m) const
{
   if (m < qedqcd.displayMbMb()) {
      return 4.0;
   } else if (m < qedqcd.displayPoleMt()) {
      return 5.0;
   } else {
      return 6.0;
   }
}

double NUTSMSSM_effective_couplings::scalar_scaling_factor(double m) const
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

double NUTSMSSM_effective_couplings::pseudoscalar_scaling_factor(double m) const
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

double NUTSMSSM_effective_couplings::get_hhVPVP_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(Mhh)(gO1);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CphhVPVP(gO1));
}

double NUTSMSSM_effective_couplings::get_hhVGVG_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(Mhh)(gO1);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CphhVGVG(gO1));
}

double NUTSMSSM_effective_couplings::get_AhVPVP_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(MAh)(gO1);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CpAhVPVP(gO1));
}

double NUTSMSSM_effective_couplings::get_AhVGVG_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(MAh)(gO1);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CpAhVGVG(gO1));
}

std::complex<double> NUTSMSSM_effective_couplings::CphhSdconjSd(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYd = MODELPARAMETER(TYd);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vS = MODELPARAMETER(vS);
   const auto vu = MODELPARAMETER(vu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3159;
   std::complex<double> tmp_3160;
   std::complex<double> tmp_3161;
   std::complex<double> tmp_3162;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3163;
      std::complex<double> tmp_3164;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3164 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3163 += tmp_3164;
      tmp_3162 += (Conj(ZD(gt2,j2))) * tmp_3163;
   }
   tmp_3161 += tmp_3162;
   tmp_3160 += (Conj(Lambdax)) * tmp_3161;
   std::complex<double> tmp_3165;
   std::complex<double> tmp_3166;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3167;
      std::complex<double> tmp_3168;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3168 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3167 += tmp_3168;
      tmp_3166 += (ZD(gt3,j2)) * tmp_3167;
   }
   tmp_3165 += tmp_3166;
   tmp_3160 += (Lambdax) * tmp_3165;
   tmp_3159 += (6*vu*Conj(ZH(gt1,2))) * tmp_3160;
   std::complex<double> tmp_3169;
   std::complex<double> tmp_3170;
   std::complex<double> tmp_3171;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3171 += Conj(ZD(gt2,j1))*ZD(gt3,j1);
   }
   tmp_3170 += tmp_3171;
   tmp_3169 += (-(vu*(0.6*Sqr(g1) + 3*Sqr(g2)))) * tmp_3170;
   std::complex<double> tmp_3172;
   std::complex<double> tmp_3173;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3173 += Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1);
   }
   tmp_3172 += tmp_3173;
   tmp_3169 += (-1.2*vu*Sqr(g1)) * tmp_3172;
   std::complex<double> tmp_3174;
   std::complex<double> tmp_3175;
   std::complex<double> tmp_3176;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3177;
      std::complex<double> tmp_3178;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3178 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3177 += tmp_3178;
      tmp_3176 += (Conj(ZD(gt2,j2))) * tmp_3177;
   }
   tmp_3175 += tmp_3176;
   tmp_3174 += (vS*Conj(Lambdax)) * tmp_3175;
   std::complex<double> tmp_3179;
   std::complex<double> tmp_3180;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3181;
      std::complex<double> tmp_3182;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3182 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3181 += tmp_3182;
      tmp_3180 += (Conj(ZD(gt2,j2))) * tmp_3181;
   }
   tmp_3179 += tmp_3180;
   tmp_3174 += (1.4142135623730951*Conj(Mu)) * tmp_3179;
   std::complex<double> tmp_3183;
   std::complex<double> tmp_3184;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3185;
      std::complex<double> tmp_3186;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3186 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3185 += tmp_3186;
      tmp_3184 += (ZD(gt3,j2)) * tmp_3185;
   }
   tmp_3183 += tmp_3184;
   tmp_3174 += (vS*Lambdax + 1.4142135623730951*Mu) * tmp_3183;
   tmp_3169 += (6) * tmp_3174;
   tmp_3159 += (Conj(ZH(gt1,1))) * tmp_3169;
   std::complex<double> tmp_3187;
   std::complex<double> tmp_3188;
   std::complex<double> tmp_3189;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3189 += Conj(ZD(gt2,j1))*ZD(gt3,j1);
   }
   tmp_3188 += tmp_3189;
   tmp_3187 += (vd*(0.6*Sqr(g1) + 3*Sqr(g2))) * tmp_3188;
   std::complex<double> tmp_3190;
   std::complex<double> tmp_3191;
   std::complex<double> tmp_3192;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3192 += Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1);
   }
   tmp_3191 += tmp_3192;
   tmp_3190 += (0.6*vd*Sqr(g1)) * tmp_3191;
   std::complex<double> tmp_3193;
   std::complex<double> tmp_3194;
   std::complex<double> tmp_3195;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3196;
      std::complex<double> tmp_3197;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3197 += ZD(gt3,3 + j1)*TYd(j1,j2);
      }
      tmp_3196 += tmp_3197;
      tmp_3195 += (Conj(ZD(gt2,j2))) * tmp_3196;
   }
   tmp_3194 += tmp_3195;
   tmp_3193 += (1.4142135623730951) * tmp_3194;
   std::complex<double> tmp_3198;
   std::complex<double> tmp_3199;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3200;
      std::complex<double> tmp_3201;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3201 += Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_3200 += tmp_3201;
      tmp_3199 += (ZD(gt3,j2)) * tmp_3200;
   }
   tmp_3198 += tmp_3199;
   tmp_3193 += (1.4142135623730951) * tmp_3198;
   std::complex<double> tmp_3202;
   std::complex<double> tmp_3203;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3204;
      std::complex<double> tmp_3205;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3206;
         std::complex<double> tmp_3207;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3207 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_3206 += tmp_3207;
         tmp_3205 += (ZD(gt3,3 + j2)) * tmp_3206;
      }
      tmp_3204 += tmp_3205;
      tmp_3203 += (Conj(ZD(gt2,3 + j3))) * tmp_3204;
   }
   tmp_3202 += tmp_3203;
   std::complex<double> tmp_3208;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3209;
      std::complex<double> tmp_3210;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3211;
         std::complex<double> tmp_3212;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3212 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_3211 += tmp_3212;
         tmp_3210 += (Conj(ZD(gt2,j2))) * tmp_3211;
      }
      tmp_3209 += tmp_3210;
      tmp_3208 += (ZD(gt3,j3)) * tmp_3209;
   }
   tmp_3202 += tmp_3208;
   tmp_3193 += (2*vd) * tmp_3202;
   tmp_3190 += (-3) * tmp_3193;
   tmp_3187 += (2) * tmp_3190;
   tmp_3159 += (Conj(ZH(gt1,0))) * tmp_3187;
   result += (0.08333333333333333) * tmp_3159;

   return result;
}

std::complex<double> NUTSMSSM_effective_couplings::CphhSuconjSu(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYu = MODELPARAMETER(TYu);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vS = MODELPARAMETER(vS);
   const auto vu = MODELPARAMETER(vu);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3213;
   std::complex<double> tmp_3214;
   std::complex<double> tmp_3215;
   std::complex<double> tmp_3216;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3217;
      std::complex<double> tmp_3218;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3218 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_3217 += tmp_3218;
      tmp_3216 += (Conj(ZU(gt2,j2))) * tmp_3217;
   }
   tmp_3215 += tmp_3216;
   tmp_3214 += (Conj(Lambdax)) * tmp_3215;
   std::complex<double> tmp_3219;
   std::complex<double> tmp_3220;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3221;
      std::complex<double> tmp_3222;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3222 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3221 += tmp_3222;
      tmp_3220 += (ZU(gt3,j2)) * tmp_3221;
   }
   tmp_3219 += tmp_3220;
   tmp_3214 += (Lambdax) * tmp_3219;
   tmp_3213 += (6*vd*Conj(ZH(gt1,2))) * tmp_3214;
   std::complex<double> tmp_3223;
   std::complex<double> tmp_3224;
   std::complex<double> tmp_3225;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3225 += Conj(ZU(gt2,j1))*ZU(gt3,j1);
   }
   tmp_3224 += tmp_3225;
   tmp_3223 += (vd*(0.6*Sqr(g1) - 3*Sqr(g2))) * tmp_3224;
   std::complex<double> tmp_3226;
   std::complex<double> tmp_3227;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3227 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_3226 += tmp_3227;
   tmp_3223 += (-2.4*vd*Sqr(g1)) * tmp_3226;
   std::complex<double> tmp_3228;
   std::complex<double> tmp_3229;
   std::complex<double> tmp_3230;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3231;
      std::complex<double> tmp_3232;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3232 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_3231 += tmp_3232;
      tmp_3230 += (Conj(ZU(gt2,j2))) * tmp_3231;
   }
   tmp_3229 += tmp_3230;
   tmp_3228 += (vS*Conj(Lambdax)) * tmp_3229;
   std::complex<double> tmp_3233;
   std::complex<double> tmp_3234;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3235;
      std::complex<double> tmp_3236;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3236 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_3235 += tmp_3236;
      tmp_3234 += (Conj(ZU(gt2,j2))) * tmp_3235;
   }
   tmp_3233 += tmp_3234;
   tmp_3228 += (1.4142135623730951*Conj(Mu)) * tmp_3233;
   std::complex<double> tmp_3237;
   std::complex<double> tmp_3238;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3239;
      std::complex<double> tmp_3240;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3240 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3239 += tmp_3240;
      tmp_3238 += (ZU(gt3,j2)) * tmp_3239;
   }
   tmp_3237 += tmp_3238;
   tmp_3228 += (vS*Lambdax + 1.4142135623730951*Mu) * tmp_3237;
   tmp_3223 += (6) * tmp_3228;
   tmp_3213 += (Conj(ZH(gt1,0))) * tmp_3223;
   std::complex<double> tmp_3241;
   std::complex<double> tmp_3242;
   std::complex<double> tmp_3243;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3243 += Conj(ZU(gt2,j1))*ZU(gt3,j1);
   }
   tmp_3242 += tmp_3243;
   tmp_3241 += (vu*(0.6*Sqr(g1) - 3*Sqr(g2))) * tmp_3242;
   std::complex<double> tmp_3244;
   std::complex<double> tmp_3245;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3245 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_3244 += tmp_3245;
   tmp_3241 += (-2.4*vu*Sqr(g1)) * tmp_3244;
   std::complex<double> tmp_3246;
   std::complex<double> tmp_3247;
   std::complex<double> tmp_3248;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3249;
      std::complex<double> tmp_3250;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3250 += ZU(gt3,3 + j1)*TYu(j1,j2);
      }
      tmp_3249 += tmp_3250;
      tmp_3248 += (Conj(ZU(gt2,j2))) * tmp_3249;
   }
   tmp_3247 += tmp_3248;
   tmp_3246 += (1.4142135623730951) * tmp_3247;
   std::complex<double> tmp_3251;
   std::complex<double> tmp_3252;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3253;
      std::complex<double> tmp_3254;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3254 += Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_3253 += tmp_3254;
      tmp_3252 += (ZU(gt3,j2)) * tmp_3253;
   }
   tmp_3251 += tmp_3252;
   tmp_3246 += (1.4142135623730951) * tmp_3251;
   std::complex<double> tmp_3255;
   std::complex<double> tmp_3256;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3257;
      std::complex<double> tmp_3258;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3259;
         std::complex<double> tmp_3260;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3260 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_3259 += tmp_3260;
         tmp_3258 += (ZU(gt3,3 + j2)) * tmp_3259;
      }
      tmp_3257 += tmp_3258;
      tmp_3256 += (Conj(ZU(gt2,3 + j3))) * tmp_3257;
   }
   tmp_3255 += tmp_3256;
   std::complex<double> tmp_3261;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3262;
      std::complex<double> tmp_3263;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3264;
         std::complex<double> tmp_3265;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3265 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_3264 += tmp_3265;
         tmp_3263 += (Conj(ZU(gt2,j2))) * tmp_3264;
      }
      tmp_3262 += tmp_3263;
      tmp_3261 += (ZU(gt3,j3)) * tmp_3262;
   }
   tmp_3255 += tmp_3261;
   tmp_3246 += (2*vu) * tmp_3255;
   tmp_3241 += (6) * tmp_3246;
   tmp_3213 += (-Conj(ZH(gt1,1))) * tmp_3241;
   result += (0.08333333333333333) * tmp_3213;

   return result;
}

std::complex<double> NUTSMSSM_effective_couplings::CphhSeconjSe(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYe = MODELPARAMETER(TYe);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vS = MODELPARAMETER(vS);
   const auto vu = MODELPARAMETER(vu);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3266;
   std::complex<double> tmp_3267;
   std::complex<double> tmp_3268;
   std::complex<double> tmp_3269;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3270;
      std::complex<double> tmp_3271;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3271 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3270 += tmp_3271;
      tmp_3269 += (Conj(ZE(gt2,j2))) * tmp_3270;
   }
   tmp_3268 += tmp_3269;
   tmp_3267 += (Conj(Lambdax)) * tmp_3268;
   std::complex<double> tmp_3272;
   std::complex<double> tmp_3273;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3274;
      std::complex<double> tmp_3275;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3275 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3274 += tmp_3275;
      tmp_3273 += (ZE(gt3,j2)) * tmp_3274;
   }
   tmp_3272 += tmp_3273;
   tmp_3267 += (Lambdax) * tmp_3272;
   tmp_3266 += (2*vu*Conj(ZH(gt1,2))) * tmp_3267;
   std::complex<double> tmp_3276;
   std::complex<double> tmp_3277;
   std::complex<double> tmp_3278;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3278 += Conj(ZE(gt2,j1))*ZE(gt3,j1);
   }
   tmp_3277 += tmp_3278;
   tmp_3276 += (vu*(0.6*Sqr(g1) - Sqr(g2))) * tmp_3277;
   std::complex<double> tmp_3279;
   std::complex<double> tmp_3280;
   std::complex<double> tmp_3281;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3281 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_3280 += tmp_3281;
   tmp_3279 += (-0.6*vu*Sqr(g1)) * tmp_3280;
   std::complex<double> tmp_3282;
   std::complex<double> tmp_3283;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3284;
      std::complex<double> tmp_3285;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3285 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3284 += tmp_3285;
      tmp_3283 += (Conj(ZE(gt2,j2))) * tmp_3284;
   }
   tmp_3282 += tmp_3283;
   tmp_3279 += (vS*Conj(Lambdax)) * tmp_3282;
   std::complex<double> tmp_3286;
   std::complex<double> tmp_3287;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3288;
      std::complex<double> tmp_3289;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3289 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3288 += tmp_3289;
      tmp_3287 += (Conj(ZE(gt2,j2))) * tmp_3288;
   }
   tmp_3286 += tmp_3287;
   tmp_3279 += (1.4142135623730951*Conj(Mu)) * tmp_3286;
   std::complex<double> tmp_3290;
   std::complex<double> tmp_3291;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3292;
      std::complex<double> tmp_3293;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3293 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3292 += tmp_3293;
      tmp_3291 += (ZE(gt3,j2)) * tmp_3292;
   }
   tmp_3290 += tmp_3291;
   tmp_3279 += (vS*Lambdax) * tmp_3290;
   std::complex<double> tmp_3294;
   std::complex<double> tmp_3295;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3296;
      std::complex<double> tmp_3297;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3297 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3296 += tmp_3297;
      tmp_3295 += (ZE(gt3,j2)) * tmp_3296;
   }
   tmp_3294 += tmp_3295;
   tmp_3279 += (1.4142135623730951*Mu) * tmp_3294;
   tmp_3276 += (2) * tmp_3279;
   tmp_3266 += (Conj(ZH(gt1,1))) * tmp_3276;
   std::complex<double> tmp_3298;
   std::complex<double> tmp_3299;
   std::complex<double> tmp_3300;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3300 += Conj(ZE(gt2,j1))*ZE(gt3,j1);
   }
   tmp_3299 += tmp_3300;
   tmp_3298 += (vd*(0.6*Sqr(g1) - Sqr(g2))) * tmp_3299;
   std::complex<double> tmp_3301;
   std::complex<double> tmp_3302;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3302 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_3301 += tmp_3302;
   tmp_3298 += (-1.2*vd*Sqr(g1)) * tmp_3301;
   std::complex<double> tmp_3303;
   std::complex<double> tmp_3304;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3305;
      std::complex<double> tmp_3306;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3306 += ZE(gt3,3 + j1)*TYe(j1,j2);
      }
      tmp_3305 += tmp_3306;
      tmp_3304 += (Conj(ZE(gt2,j2))) * tmp_3305;
   }
   tmp_3303 += tmp_3304;
   tmp_3298 += (2.8284271247461903) * tmp_3303;
   std::complex<double> tmp_3307;
   std::complex<double> tmp_3308;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3309;
      std::complex<double> tmp_3310;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3310 += Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_3309 += tmp_3310;
      tmp_3308 += (ZE(gt3,j2)) * tmp_3309;
   }
   tmp_3307 += tmp_3308;
   tmp_3298 += (2.8284271247461903) * tmp_3307;
   std::complex<double> tmp_3311;
   std::complex<double> tmp_3312;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3313;
      std::complex<double> tmp_3314;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3315;
         std::complex<double> tmp_3316;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3316 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_3315 += tmp_3316;
         tmp_3314 += (ZE(gt3,3 + j2)) * tmp_3315;
      }
      tmp_3313 += tmp_3314;
      tmp_3312 += (Conj(ZE(gt2,3 + j3))) * tmp_3313;
   }
   tmp_3311 += tmp_3312;
   tmp_3298 += (4*vd) * tmp_3311;
   std::complex<double> tmp_3317;
   std::complex<double> tmp_3318;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3319;
      std::complex<double> tmp_3320;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3321;
         std::complex<double> tmp_3322;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3322 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_3321 += tmp_3322;
         tmp_3320 += (Conj(ZE(gt2,j2))) * tmp_3321;
      }
      tmp_3319 += tmp_3320;
      tmp_3318 += (ZE(gt3,j3)) * tmp_3319;
   }
   tmp_3317 += tmp_3318;
   tmp_3298 += (4*vd) * tmp_3317;
   tmp_3266 += (-Conj(ZH(gt1,0))) * tmp_3298;
   result += (0.25) * tmp_3266;

   return result;
}

std::complex<double> NUTSMSSM_effective_couplings::CphhHpmconjHpm(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TLambdax = MODELPARAMETER(TLambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto MS = MODELPARAMETER(MS);
   const auto vd = MODELPARAMETER(vd);
   const auto vS = MODELPARAMETER(vS);
   const auto vu = MODELPARAMETER(vu);
   const auto Kappa = MODELPARAMETER(Kappa);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   result = 0.25*(-(Conj(ZH(gt1,0))*(ZP(gt2,0)*(vd*(0.6*Sqr(g1) + Sqr(g2))*ZP(
      gt3,0) + vu*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gt3,1)) + ZP(gt2,1)*(vu*(-2*
      AbsSqr(Lambdax) + Sqr(g2))*ZP(gt3,0) + vd*(-0.6*Sqr(g1) + Sqr(g2))*ZP(gt3,1)
      ))) + Conj(ZH(gt1,1))*(ZP(gt2,0)*(vu*(0.6*Sqr(g1) - Sqr(g2))*ZP(gt3,0) - vd*
      (-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gt3,1)) - ZP(gt2,1)*(vd*(-2*AbsSqr(Lambdax
      ) + Sqr(g2))*ZP(gt3,0) + vu*(0.6*Sqr(g1) + Sqr(g2))*ZP(gt3,1))) - 2*Conj(ZH(
      gt1,2))*(1.4142135623730951*Conj(TLambdax)*ZP(gt2,1)*ZP(gt3,0) +
      1.4142135623730951*Conj(MS)*Lambdax*ZP(gt2,0)*ZP(gt3,1) + 2*vS*Conj(Kappa)*
      Lambdax*ZP(gt2,0)*ZP(gt3,1) + 1.4142135623730951*TLambdax*ZP(gt2,0)*ZP(gt3,1
      ) + 1.4142135623730951*Conj(Mu)*Lambdax*(ZP(gt2,0)*ZP(gt3,0) + ZP(gt2,1)*ZP(
      gt3,1)) + Conj(Lambdax)*((2*vS*Lambdax + 1.4142135623730951*Mu)*ZP(gt2,0)*ZP
      (gt3,0) + ZP(gt2,1)*((1.4142135623730951*MS + 2*vS*Kappa)*ZP(gt3,0) + (2*vS*
      Lambdax + 1.4142135623730951*Mu)*ZP(gt3,1)))));

   return result;
}

std::complex<double> NUTSMSSM_effective_couplings::CpChahhbarChaPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   std::complex<double> result;

   result = -0.7071067811865475*(g2*Conj(UM(gt1,0))*Conj(UP(gt3,1))*Conj(ZH(gt2
      ,1)) + Conj(UM(gt1,1))*(g2*Conj(UP(gt3,0))*Conj(ZH(gt2,0)) + Conj(UP(gt3,1))
      *Conj(ZH(gt2,2))*Lambdax));

   return result;
}

std::complex<double> NUTSMSSM_effective_couplings::CpFehhbarFePL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Ye = MODELPARAMETER(Ye);

   std::complex<double> result;

   std::complex<double> tmp_3323;
   std::complex<double> tmp_3324;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3325;
      std::complex<double> tmp_3326;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3326 += Conj(ZER(gt3,j1))*Ye(j1,j2);
      }
      tmp_3325 += tmp_3326;
      tmp_3324 += (Conj(ZEL(gt1,j2))) * tmp_3325;
   }
   tmp_3323 += tmp_3324;
   result += (-0.7071067811865475*Conj(ZH(gt2,0))) * tmp_3323;

   return result;
}

std::complex<double> NUTSMSSM_effective_couplings::CpFdhhbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_3327;
   std::complex<double> tmp_3328;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3329;
      std::complex<double> tmp_3330;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3330 += Conj(ZDR(gt3,j1))*Yd(j1,j2);
      }
      tmp_3329 += tmp_3330;
      tmp_3328 += (Conj(ZDL(gt1,j2))) * tmp_3329;
   }
   tmp_3327 += tmp_3328;
   result += (-0.7071067811865475*Conj(ZH(gt2,0))) * tmp_3327;

   return result;
}

std::complex<double> NUTSMSSM_effective_couplings::CpFuhhbarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_3331;
   std::complex<double> tmp_3332;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3333;
      std::complex<double> tmp_3334;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3334 += Conj(ZUR(gt3,j1))*Yu(j1,j2);
      }
      tmp_3333 += tmp_3334;
      tmp_3332 += (Conj(ZUL(gt1,j2))) * tmp_3333;
   }
   tmp_3331 += tmp_3332;
   result += (-0.7071067811865475*Conj(ZH(gt2,1))) * tmp_3331;

   return result;
}

std::complex<double> NUTSMSSM_effective_couplings::CphhVWmconjVWm(unsigned gt1) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);

   std::complex<double> result;

   result = 0.5*(vd*Conj(ZH(gt1,0)) + vu*Conj(ZH(gt1,1)))*Sqr(g2);

   return result;
}

std::complex<double> NUTSMSSM_effective_couplings::CpAhSdconjSd(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYd = MODELPARAMETER(TYd);
   const auto vS = MODELPARAMETER(vS);
   const auto vu = MODELPARAMETER(vu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3335;
   std::complex<double> tmp_3336;
   std::complex<double> tmp_3337;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3338;
      std::complex<double> tmp_3339;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3339 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3338 += tmp_3339;
      tmp_3337 += (Conj(ZD(gt2,j2))) * tmp_3338;
   }
   tmp_3336 += tmp_3337;
   tmp_3335 += (1.4142135623730951*Conj(Mu)*Conj(ZA(gt1,1))) * tmp_3336;
   std::complex<double> tmp_3340;
   std::complex<double> tmp_3341;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3342;
      std::complex<double> tmp_3343;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3343 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3342 += tmp_3343;
      tmp_3341 += (Conj(ZD(gt2,j2))) * tmp_3342;
   }
   tmp_3340 += tmp_3341;
   tmp_3335 += (Conj(Lambdax)*(vS*Conj(ZA(gt1,1)) + vu*Conj(ZA(gt1,2)))) *
      tmp_3340;
   std::complex<double> tmp_3344;
   std::complex<double> tmp_3345;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3346;
      std::complex<double> tmp_3347;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3347 += ZD(gt3,3 + j1)*TYd(j1,j2);
      }
      tmp_3346 += tmp_3347;
      tmp_3345 += (Conj(ZD(gt2,j2))) * tmp_3346;
   }
   tmp_3344 += tmp_3345;
   tmp_3335 += (1.4142135623730951*Conj(ZA(gt1,0))) * tmp_3344;
   std::complex<double> tmp_3348;
   std::complex<double> tmp_3349;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3350;
      std::complex<double> tmp_3351;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3351 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3350 += tmp_3351;
      tmp_3349 += (ZD(gt3,j2)) * tmp_3350;
   }
   tmp_3348 += tmp_3349;
   tmp_3335 += (-(vS*Conj(ZA(gt1,1))*Lambdax)) * tmp_3348;
   std::complex<double> tmp_3352;
   std::complex<double> tmp_3353;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3354;
      std::complex<double> tmp_3355;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3355 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3354 += tmp_3355;
      tmp_3353 += (ZD(gt3,j2)) * tmp_3354;
   }
   tmp_3352 += tmp_3353;
   tmp_3335 += (-1.4142135623730951*Conj(ZA(gt1,1))*Mu) * tmp_3352;
   std::complex<double> tmp_3356;
   std::complex<double> tmp_3357;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3358;
      std::complex<double> tmp_3359;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3359 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3358 += tmp_3359;
      tmp_3357 += (ZD(gt3,j2)) * tmp_3358;
   }
   tmp_3356 += tmp_3357;
   tmp_3335 += (-(vu*Conj(ZA(gt1,2))*Lambdax)) * tmp_3356;
   std::complex<double> tmp_3360;
   std::complex<double> tmp_3361;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3362;
      std::complex<double> tmp_3363;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3363 += Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_3362 += tmp_3363;
      tmp_3361 += (ZD(gt3,j2)) * tmp_3362;
   }
   tmp_3360 += tmp_3361;
   tmp_3335 += (-1.4142135623730951*Conj(ZA(gt1,0))) * tmp_3360;
   result += (std::complex<double>(0,-0.5)) * tmp_3335;

   return result;
}

std::complex<double> NUTSMSSM_effective_couplings::CpAhSuconjSu(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYu = MODELPARAMETER(TYu);
   const auto vd = MODELPARAMETER(vd);
   const auto vS = MODELPARAMETER(vS);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3364;
   std::complex<double> tmp_3365;
   std::complex<double> tmp_3366;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3367;
      std::complex<double> tmp_3368;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3368 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_3367 += tmp_3368;
      tmp_3366 += (Conj(ZU(gt2,j2))) * tmp_3367;
   }
   tmp_3365 += tmp_3366;
   tmp_3364 += (1.4142135623730951*Conj(Mu)*Conj(ZA(gt1,0))) * tmp_3365;
   std::complex<double> tmp_3369;
   std::complex<double> tmp_3370;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3371;
      std::complex<double> tmp_3372;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3372 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_3371 += tmp_3372;
      tmp_3370 += (Conj(ZU(gt2,j2))) * tmp_3371;
   }
   tmp_3369 += tmp_3370;
   tmp_3364 += (Conj(Lambdax)*(vS*Conj(ZA(gt1,0)) + vd*Conj(ZA(gt1,2)))) *
      tmp_3369;
   std::complex<double> tmp_3373;
   std::complex<double> tmp_3374;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3375;
      std::complex<double> tmp_3376;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3376 += ZU(gt3,3 + j1)*TYu(j1,j2);
      }
      tmp_3375 += tmp_3376;
      tmp_3374 += (Conj(ZU(gt2,j2))) * tmp_3375;
   }
   tmp_3373 += tmp_3374;
   tmp_3364 += (1.4142135623730951*Conj(ZA(gt1,1))) * tmp_3373;
   std::complex<double> tmp_3377;
   std::complex<double> tmp_3378;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3379;
      std::complex<double> tmp_3380;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3380 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3379 += tmp_3380;
      tmp_3378 += (ZU(gt3,j2)) * tmp_3379;
   }
   tmp_3377 += tmp_3378;
   tmp_3364 += (-(vS*Conj(ZA(gt1,0))*Lambdax)) * tmp_3377;
   std::complex<double> tmp_3381;
   std::complex<double> tmp_3382;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3383;
      std::complex<double> tmp_3384;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3384 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3383 += tmp_3384;
      tmp_3382 += (ZU(gt3,j2)) * tmp_3383;
   }
   tmp_3381 += tmp_3382;
   tmp_3364 += (-1.4142135623730951*Conj(ZA(gt1,0))*Mu) * tmp_3381;
   std::complex<double> tmp_3385;
   std::complex<double> tmp_3386;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3387;
      std::complex<double> tmp_3388;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3388 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3387 += tmp_3388;
      tmp_3386 += (ZU(gt3,j2)) * tmp_3387;
   }
   tmp_3385 += tmp_3386;
   tmp_3364 += (-(vd*Conj(ZA(gt1,2))*Lambdax)) * tmp_3385;
   std::complex<double> tmp_3389;
   std::complex<double> tmp_3390;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3391;
      std::complex<double> tmp_3392;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3392 += Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_3391 += tmp_3392;
      tmp_3390 += (ZU(gt3,j2)) * tmp_3391;
   }
   tmp_3389 += tmp_3390;
   tmp_3364 += (-1.4142135623730951*Conj(ZA(gt1,1))) * tmp_3389;
   result += (std::complex<double>(0,-0.5)) * tmp_3364;

   return result;
}

std::complex<double> NUTSMSSM_effective_couplings::CpAhSeconjSe(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYe = MODELPARAMETER(TYe);
   const auto vS = MODELPARAMETER(vS);
   const auto vu = MODELPARAMETER(vu);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3393;
   std::complex<double> tmp_3394;
   std::complex<double> tmp_3395;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3396;
      std::complex<double> tmp_3397;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3397 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3396 += tmp_3397;
      tmp_3395 += (Conj(ZE(gt2,j2))) * tmp_3396;
   }
   tmp_3394 += tmp_3395;
   tmp_3393 += (1.4142135623730951*Conj(Mu)*Conj(ZA(gt1,1))) * tmp_3394;
   std::complex<double> tmp_3398;
   std::complex<double> tmp_3399;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3400;
      std::complex<double> tmp_3401;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3401 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3400 += tmp_3401;
      tmp_3399 += (Conj(ZE(gt2,j2))) * tmp_3400;
   }
   tmp_3398 += tmp_3399;
   tmp_3393 += (Conj(Lambdax)*(vS*Conj(ZA(gt1,1)) + vu*Conj(ZA(gt1,2)))) *
      tmp_3398;
   std::complex<double> tmp_3402;
   std::complex<double> tmp_3403;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3404;
      std::complex<double> tmp_3405;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3405 += ZE(gt3,3 + j1)*TYe(j1,j2);
      }
      tmp_3404 += tmp_3405;
      tmp_3403 += (Conj(ZE(gt2,j2))) * tmp_3404;
   }
   tmp_3402 += tmp_3403;
   tmp_3393 += (1.4142135623730951*Conj(ZA(gt1,0))) * tmp_3402;
   std::complex<double> tmp_3406;
   std::complex<double> tmp_3407;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3408;
      std::complex<double> tmp_3409;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3409 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3408 += tmp_3409;
      tmp_3407 += (ZE(gt3,j2)) * tmp_3408;
   }
   tmp_3406 += tmp_3407;
   tmp_3393 += (-(vS*Conj(ZA(gt1,1))*Lambdax)) * tmp_3406;
   std::complex<double> tmp_3410;
   std::complex<double> tmp_3411;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3412;
      std::complex<double> tmp_3413;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3413 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3412 += tmp_3413;
      tmp_3411 += (ZE(gt3,j2)) * tmp_3412;
   }
   tmp_3410 += tmp_3411;
   tmp_3393 += (-1.4142135623730951*Conj(ZA(gt1,1))*Mu) * tmp_3410;
   std::complex<double> tmp_3414;
   std::complex<double> tmp_3415;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3416;
      std::complex<double> tmp_3417;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3417 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3416 += tmp_3417;
      tmp_3415 += (ZE(gt3,j2)) * tmp_3416;
   }
   tmp_3414 += tmp_3415;
   tmp_3393 += (-(vu*Conj(ZA(gt1,2))*Lambdax)) * tmp_3414;
   std::complex<double> tmp_3418;
   std::complex<double> tmp_3419;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3420;
      std::complex<double> tmp_3421;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3421 += Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_3420 += tmp_3421;
      tmp_3419 += (ZE(gt3,j2)) * tmp_3420;
   }
   tmp_3418 += tmp_3419;
   tmp_3393 += (-1.4142135623730951*Conj(ZA(gt1,0))) * tmp_3418;
   result += (std::complex<double>(0,-0.5)) * tmp_3393;

   return result;
}

std::complex<double> NUTSMSSM_effective_couplings::CpAhHpmconjHpm(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TLambdax = MODELPARAMETER(TLambdax);
   const auto g2 = MODELPARAMETER(g2);
   const auto MS = MODELPARAMETER(MS);
   const auto vd = MODELPARAMETER(vd);
   const auto vS = MODELPARAMETER(vS);
   const auto vu = MODELPARAMETER(vu);
   const auto Kappa = MODELPARAMETER(Kappa);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   result = std::complex<double>(0,0.25)*(2.8284271247461903*Conj(TLambdax)*
      Conj(ZA(gt1,2))*ZP(gt2,1)*ZP(gt3,0) - vu*Conj(ZA(gt1,0))*Sqr(g2)*ZP(gt2,1)*
      ZP(gt3,0) - vd*Conj(ZA(gt1,1))*Sqr(g2)*ZP(gt2,1)*ZP(gt3,0) +
      2.8284271247461903*Conj(MS)*Conj(ZA(gt1,2))*Lambdax*ZP(gt2,0)*ZP(gt3,1) + 4*
      vS*Conj(Kappa)*Conj(ZA(gt1,2))*Lambdax*ZP(gt2,0)*ZP(gt3,1) + vu*Conj(ZA(gt1,
      0))*Sqr(g2)*ZP(gt2,0)*ZP(gt3,1) + vd*Conj(ZA(gt1,1))*Sqr(g2)*ZP(gt2,0)*ZP(
      gt3,1) - 2.8284271247461903*Conj(ZA(gt1,2))*TLambdax*ZP(gt2,0)*ZP(gt3,1) -
      2.8284271247461903*Conj(Mu)*Conj(ZA(gt1,2))*Lambdax*(ZP(gt2,0)*ZP(gt3,0) +
      ZP(gt2,1)*ZP(gt3,1)) + 2*Conj(Lambdax)*((vu*Conj(ZA(gt1,0)) + vd*Conj(ZA(gt1
      ,1)))*Lambdax*(ZP(gt2,1)*ZP(gt3,0) - ZP(gt2,0)*ZP(gt3,1)) + Conj(ZA(gt1,2))*
      (1.4142135623730951*Mu*ZP(gt2,0)*ZP(gt3,0) + ZP(gt2,1)*(-((
      1.4142135623730951*MS + 2*vS*Kappa)*ZP(gt3,0)) + 1.4142135623730951*Mu*ZP(
      gt3,1)))));

   return result;
}

std::complex<double> NUTSMSSM_effective_couplings::CpAhChabarChaPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   std::complex<double> result;

   result = std::complex<double>(0.,0.7071067811865475)*(g2*Conj(UM(gt2,0))*
      Conj(UP(gt3,1))*Conj(ZA(gt1,1)) + Conj(UM(gt2,1))*(g2*Conj(UP(gt3,0))*Conj(
      ZA(gt1,0)) - Conj(UP(gt3,1))*Conj(ZA(gt1,2))*Lambdax));

   return result;
}

std::complex<double> NUTSMSSM_effective_couplings::CpAhFebarFePL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Ye = MODELPARAMETER(Ye);

   std::complex<double> result;

   std::complex<double> tmp_3422;
   std::complex<double> tmp_3423;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3424;
      std::complex<double> tmp_3425;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3425 += Conj(ZER(gt3,j1))*Ye(j1,j2);
      }
      tmp_3424 += tmp_3425;
      tmp_3423 += (Conj(ZEL(gt2,j2))) * tmp_3424;
   }
   tmp_3422 += tmp_3423;
   result += (std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gt1,0))) *
      tmp_3422;

   return result;
}

std::complex<double> NUTSMSSM_effective_couplings::CpAhFdbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_3426;
   std::complex<double> tmp_3427;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3428;
      std::complex<double> tmp_3429;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3429 += Conj(ZDR(gt3,j1))*Yd(j1,j2);
      }
      tmp_3428 += tmp_3429;
      tmp_3427 += (Conj(ZDL(gt2,j2))) * tmp_3428;
   }
   tmp_3426 += tmp_3427;
   result += (std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gt1,0))) *
      tmp_3426;

   return result;
}

std::complex<double> NUTSMSSM_effective_couplings::CpAhFubarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_3430;
   std::complex<double> tmp_3431;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3432;
      std::complex<double> tmp_3433;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3433 += Conj(ZUR(gt3,j1))*Yu(j1,j2);
      }
      tmp_3432 += tmp_3433;
      tmp_3431 += (Conj(ZUL(gt2,j2))) * tmp_3432;
   }
   tmp_3430 += tmp_3431;
   result += (std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gt1,1))) *
      tmp_3430;

   return result;
}

void NUTSMSSM_effective_couplings::calculate_eff_CphhVPVP(unsigned gO1)
{
   const auto vd = MODELPARAMETER(vd);
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
   result += -0.5 * CphhVWmconjVWm(gO1) * vev * AS1(decay_scale / Sqr(MVWm)) /
      Sqr(MVWm);

   result *= 0.1892681907127351 * physical_input.get(Physical_input::alpha_em_0
      ) * Sqrt(qedqcd.displayFermiConstant());

   ZH = saved_ZH;
   eff_CphhVPVP(gO1) = result;

}

void NUTSMSSM_effective_couplings::calculate_eff_CphhVGVG(unsigned gO1)
{
   const auto vd = MODELPARAMETER(vd);
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

void NUTSMSSM_effective_couplings::calculate_eff_CpAhVPVP(unsigned gO1)
{
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto MCha = MODELPARAMETER(MCha);
   const auto MFe = MODELPARAMETER(MFe);
   const auto MFd = MODELPARAMETER(MFd);
   const auto MFu = MODELPARAMETER(MFu);
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
   result *= std::complex<double>(2.0,0.);

   result *= 0.1892681907127351 * physical_input.get(Physical_input::alpha_em_0
      ) * Sqrt(qedqcd.displayFermiConstant());

   ZA = saved_ZA;
   eff_CpAhVPVP(gO1) = result;

}

void NUTSMSSM_effective_couplings::calculate_eff_CpAhVGVG(unsigned gO1)
{
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto g3 = MODELPARAMETER(g3);
   const auto MFd = MODELPARAMETER(MFd);
   const auto MFu = MODELPARAMETER(MFu);
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
