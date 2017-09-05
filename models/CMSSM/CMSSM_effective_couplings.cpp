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

// File generated at Tue 5 Sep 2017 13:08:49

#include "CMSSM_effective_couplings.hpp"

#include "effective_couplings.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

using namespace effective_couplings;

#define INPUTPARAMETER(parameter) model.get_input().parameter
#define MODELPARAMETER(parameter) model.get_##parameter()
#define DERIVEDPARAMETER(parameter) model.##parameter()
#define PHASE(parameter) model.get_##parameter()
#define PHYSICAL(parameter) model.get_physical().parameter

CMSSM_effective_couplings::CMSSM_effective_couplings(
   const CMSSM_mass_eigenstates& model_,
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

   , eff_CphhVPVP(Eigen::Array<std::complex<double>,2,1>::Zero()), eff_CphhVGVG
      (Eigen::Array<std::complex<double>,2,1>::Zero()), eff_CpAhVPVP(Eigen::Array<
      std::complex<double>,2,1>::Zero()), eff_CpAhVGVG(Eigen::Array<std::complex<
      double>,2,1>::Zero())

{
}

void CMSSM_effective_couplings::calculate_effective_couplings()
{
   const standard_model::Standard_model sm(initialise_SM());

   const double scale = model.get_scale();
   const Eigen::ArrayXd saved_parameters(model.get());

   const double saved_mt = PHYSICAL(MFu(2));
   PHYSICAL(MFu(2)) = qedqcd.displayPoleMt();

   const auto Mhh = PHYSICAL(Mhh);
   for (unsigned gO1 = 0; gO1 < 2; ++gO1) {
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
   for (unsigned gO1 = 1; gO1 < 2; ++gO1) {
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

void CMSSM_effective_couplings::set_model(const CMSSM_mass_eigenstates& model_)
{
   model = model_;
   copy_mixing_matrices_from_model();
}

void CMSSM_effective_couplings::copy_mixing_matrices_from_model()
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

standard_model::Standard_model CMSSM_effective_couplings::initialise_SM() const
{
   standard_model::Standard_model sm;

   sm.set_loops(2);
   sm.set_thresholds(2);
   sm.set_physical_input(physical_input);

   sm.initialise_from_input(qedqcd);

   return sm;
}

void CMSSM_effective_couplings::run_SM_strong_coupling_to(standard_model::Standard_model sm, double m)
{
   sm.run_to(m);

   model.set_g3(sm.get_g3());

}

std::complex<double> CMSSM_effective_couplings::scalar_scalar_qcd_factor(double m_decay, double m_loop) const
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

std::complex<double> CMSSM_effective_couplings::scalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) *
         scalar_diphoton_fermion_loop(m_decay, m_loop);

   }

   return result;
}

std::complex<double> CMSSM_effective_couplings::pseudoscalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) *
         pseudoscalar_diphoton_fermion_loop(m_decay, m_loop);

   }

   return result;
}

double CMSSM_effective_couplings::number_of_active_flavours(double m) const
{
   if (m < qedqcd.displayMbMb()) {
      return 4.0;
   } else if (m < qedqcd.displayPoleMt()) {
      return 5.0;
   } else {
      return 6.0;
   }
}

double CMSSM_effective_couplings::scalar_scaling_factor(double m) const
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

double CMSSM_effective_couplings::pseudoscalar_scaling_factor(double m) const
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

double CMSSM_effective_couplings::get_hhVPVP_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(Mhh)(gO1);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CphhVPVP(gO1));
}

double CMSSM_effective_couplings::get_hhVGVG_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(Mhh)(gO1);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CphhVGVG(gO1));
}

double CMSSM_effective_couplings::get_AhVPVP_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(MAh)(gO1);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CpAhVPVP(gO1));
}

double CMSSM_effective_couplings::get_AhVGVG_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(MAh)(gO1);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CpAhVGVG(gO1));
}

std::complex<double> CMSSM_effective_couplings::CphhSdconjSd(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYd = MODELPARAMETER(TYd);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3065;
   std::complex<double> tmp_3066;
   std::complex<double> tmp_3067;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3067 += Conj(ZD(gt2,j1))*ZD(gt3,j1);
   }
   tmp_3066 += tmp_3067;
   tmp_3065 += ((0.6*Sqr(g1) + 3*Sqr(g2))*(vd*ZH(gt1,0) - vu*ZH(gt1,1))) *
      tmp_3066;
   std::complex<double> tmp_3068;
   std::complex<double> tmp_3069;
   std::complex<double> tmp_3070;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3070 += Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1);
   }
   tmp_3069 += tmp_3070;
   tmp_3068 += (0.6*Sqr(g1)*(vd*ZH(gt1,0) - vu*ZH(gt1,1))) * tmp_3069;
   std::complex<double> tmp_3071;
   std::complex<double> tmp_3072;
   std::complex<double> tmp_3073;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3074;
      std::complex<double> tmp_3075;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3075 += ZD(gt3,3 + j1)*TYd(j1,j2);
      }
      tmp_3074 += tmp_3075;
      tmp_3073 += (Conj(ZD(gt2,j2))) * tmp_3074;
   }
   tmp_3072 += tmp_3073;
   tmp_3071 += (1.4142135623730951*ZH(gt1,0)) * tmp_3072;
   std::complex<double> tmp_3076;
   std::complex<double> tmp_3077;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3078;
      std::complex<double> tmp_3079;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3079 += Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_3078 += tmp_3079;
      tmp_3077 += (ZD(gt3,j2)) * tmp_3078;
   }
   tmp_3076 += tmp_3077;
   tmp_3071 += (1.4142135623730951*ZH(gt1,0)) * tmp_3076;
   std::complex<double> tmp_3080;
   std::complex<double> tmp_3081;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3082;
      std::complex<double> tmp_3083;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3084;
         std::complex<double> tmp_3085;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3085 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_3084 += tmp_3085;
         tmp_3083 += (ZD(gt3,3 + j2)) * tmp_3084;
      }
      tmp_3082 += tmp_3083;
      tmp_3081 += (Conj(ZD(gt2,3 + j3))) * tmp_3082;
   }
   tmp_3080 += tmp_3081;
   tmp_3071 += (2*vd*ZH(gt1,0)) * tmp_3080;
   std::complex<double> tmp_3086;
   std::complex<double> tmp_3087;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3088;
      std::complex<double> tmp_3089;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3090;
         std::complex<double> tmp_3091;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3091 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_3090 += tmp_3091;
         tmp_3089 += (Conj(ZD(gt2,j2))) * tmp_3090;
      }
      tmp_3088 += tmp_3089;
      tmp_3087 += (ZD(gt3,j3)) * tmp_3088;
   }
   tmp_3086 += tmp_3087;
   tmp_3071 += (2*vd*ZH(gt1,0)) * tmp_3086;
   std::complex<double> tmp_3092;
   std::complex<double> tmp_3093;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3094;
      std::complex<double> tmp_3095;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3095 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3094 += tmp_3095;
      tmp_3093 += (Conj(ZD(gt2,j2))) * tmp_3094;
   }
   tmp_3092 += tmp_3093;
   tmp_3071 += (-1.4142135623730951*Conj(Mu)*ZH(gt1,1)) * tmp_3092;
   std::complex<double> tmp_3096;
   std::complex<double> tmp_3097;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3098;
      std::complex<double> tmp_3099;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3099 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3098 += tmp_3099;
      tmp_3097 += (ZD(gt3,j2)) * tmp_3098;
   }
   tmp_3096 += tmp_3097;
   tmp_3071 += (-1.4142135623730951*Mu*ZH(gt1,1)) * tmp_3096;
   tmp_3068 += (-3) * tmp_3071;
   tmp_3065 += (2) * tmp_3068;
   result += (0.08333333333333333) * tmp_3065;

   return result;
}

std::complex<double> CMSSM_effective_couplings::CphhSuconjSu(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYu = MODELPARAMETER(TYu);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3100;
   std::complex<double> tmp_3101;
   std::complex<double> tmp_3102;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3102 += Conj(ZU(gt2,j1))*ZU(gt3,j1);
   }
   tmp_3101 += tmp_3102;
   tmp_3100 += ((0.6*Sqr(g1) - 3*Sqr(g2))*(vd*ZH(gt1,0) - vu*ZH(gt1,1))) *
      tmp_3101;
   std::complex<double> tmp_3103;
   std::complex<double> tmp_3104;
   std::complex<double> tmp_3105;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3106;
      std::complex<double> tmp_3107;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3107 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_3106 += tmp_3107;
      tmp_3105 += (Conj(ZU(gt2,j2))) * tmp_3106;
   }
   tmp_3104 += tmp_3105;
   tmp_3103 += (-4.242640687119286*Conj(Mu)*ZH(gt1,0)) * tmp_3104;
   std::complex<double> tmp_3108;
   std::complex<double> tmp_3109;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3110;
      std::complex<double> tmp_3111;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3111 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3110 += tmp_3111;
      tmp_3109 += (ZU(gt3,j2)) * tmp_3110;
   }
   tmp_3108 += tmp_3109;
   tmp_3103 += (-4.242640687119286*Mu*ZH(gt1,0)) * tmp_3108;
   std::complex<double> tmp_3112;
   std::complex<double> tmp_3113;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3114;
      std::complex<double> tmp_3115;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3115 += ZU(gt3,3 + j1)*TYu(j1,j2);
      }
      tmp_3114 += tmp_3115;
      tmp_3113 += (Conj(ZU(gt2,j2))) * tmp_3114;
   }
   tmp_3112 += tmp_3113;
   tmp_3103 += (4.242640687119286*ZH(gt1,1)) * tmp_3112;
   std::complex<double> tmp_3116;
   std::complex<double> tmp_3117;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3118;
      std::complex<double> tmp_3119;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3119 += Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_3118 += tmp_3119;
      tmp_3117 += (ZU(gt3,j2)) * tmp_3118;
   }
   tmp_3116 += tmp_3117;
   tmp_3103 += (4.242640687119286*ZH(gt1,1)) * tmp_3116;
   std::complex<double> tmp_3120;
   std::complex<double> tmp_3121;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3122;
      std::complex<double> tmp_3123;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3124;
         std::complex<double> tmp_3125;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3125 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_3124 += tmp_3125;
         tmp_3123 += (ZU(gt3,3 + j2)) * tmp_3124;
      }
      tmp_3122 += tmp_3123;
      tmp_3121 += (Conj(ZU(gt2,3 + j3))) * tmp_3122;
   }
   tmp_3120 += tmp_3121;
   tmp_3103 += (6*vu*ZH(gt1,1)) * tmp_3120;
   std::complex<double> tmp_3126;
   std::complex<double> tmp_3127;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3128;
      std::complex<double> tmp_3129;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3130;
         std::complex<double> tmp_3131;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3131 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_3130 += tmp_3131;
         tmp_3129 += (Conj(ZU(gt2,j2))) * tmp_3130;
      }
      tmp_3128 += tmp_3129;
      tmp_3127 += (ZU(gt3,j3)) * tmp_3128;
   }
   tmp_3126 += tmp_3127;
   tmp_3103 += (6*vu*ZH(gt1,1)) * tmp_3126;
   std::complex<double> tmp_3132;
   std::complex<double> tmp_3133;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3133 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_3132 += tmp_3133;
   tmp_3103 += (1.2*Sqr(g1)*(vd*ZH(gt1,0) - vu*ZH(gt1,1))) * tmp_3132;
   tmp_3100 += (-2) * tmp_3103;
   result += (0.08333333333333333) * tmp_3100;

   return result;
}

std::complex<double> CMSSM_effective_couplings::CphhSeconjSe(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYe = MODELPARAMETER(TYe);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3134;
   std::complex<double> tmp_3135;
   std::complex<double> tmp_3136;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3136 += Conj(ZE(gt2,j1))*ZE(gt3,j1);
   }
   tmp_3135 += tmp_3136;
   tmp_3134 += (-((0.6*Sqr(g1) - Sqr(g2))*(vd*ZH(gt1,0) - vu*ZH(gt1,1)))) *
      tmp_3135;
   std::complex<double> tmp_3137;
   std::complex<double> tmp_3138;
   std::complex<double> tmp_3139;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3140;
      std::complex<double> tmp_3141;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3141 += ZE(gt3,3 + j1)*TYe(j1,j2);
      }
      tmp_3140 += tmp_3141;
      tmp_3139 += (Conj(ZE(gt2,j2))) * tmp_3140;
   }
   tmp_3138 += tmp_3139;
   tmp_3137 += (-1.4142135623730951*ZH(gt1,0)) * tmp_3138;
   std::complex<double> tmp_3142;
   std::complex<double> tmp_3143;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3144;
      std::complex<double> tmp_3145;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3145 += Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_3144 += tmp_3145;
      tmp_3143 += (ZE(gt3,j2)) * tmp_3144;
   }
   tmp_3142 += tmp_3143;
   tmp_3137 += (-1.4142135623730951*ZH(gt1,0)) * tmp_3142;
   std::complex<double> tmp_3146;
   std::complex<double> tmp_3147;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3148;
      std::complex<double> tmp_3149;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3150;
         std::complex<double> tmp_3151;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3151 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_3150 += tmp_3151;
         tmp_3149 += (ZE(gt3,3 + j2)) * tmp_3150;
      }
      tmp_3148 += tmp_3149;
      tmp_3147 += (Conj(ZE(gt2,3 + j3))) * tmp_3148;
   }
   tmp_3146 += tmp_3147;
   tmp_3137 += (-2*vd*ZH(gt1,0)) * tmp_3146;
   std::complex<double> tmp_3152;
   std::complex<double> tmp_3153;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3154;
      std::complex<double> tmp_3155;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3156;
         std::complex<double> tmp_3157;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3157 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_3156 += tmp_3157;
         tmp_3155 += (Conj(ZE(gt2,j2))) * tmp_3156;
      }
      tmp_3154 += tmp_3155;
      tmp_3153 += (ZE(gt3,j3)) * tmp_3154;
   }
   tmp_3152 += tmp_3153;
   tmp_3137 += (-2*vd*ZH(gt1,0)) * tmp_3152;
   std::complex<double> tmp_3158;
   std::complex<double> tmp_3159;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3160;
      std::complex<double> tmp_3161;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3161 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3160 += tmp_3161;
      tmp_3159 += (Conj(ZE(gt2,j2))) * tmp_3160;
   }
   tmp_3158 += tmp_3159;
   tmp_3137 += (1.4142135623730951*Conj(Mu)*ZH(gt1,1)) * tmp_3158;
   std::complex<double> tmp_3162;
   std::complex<double> tmp_3163;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3164;
      std::complex<double> tmp_3165;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3165 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3164 += tmp_3165;
      tmp_3163 += (ZE(gt3,j2)) * tmp_3164;
   }
   tmp_3162 += tmp_3163;
   tmp_3137 += (1.4142135623730951*Mu*ZH(gt1,1)) * tmp_3162;
   std::complex<double> tmp_3166;
   std::complex<double> tmp_3167;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3167 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_3166 += tmp_3167;
   tmp_3137 += (0.6*Sqr(g1)*(vd*ZH(gt1,0) - vu*ZH(gt1,1))) * tmp_3166;
   tmp_3134 += (2) * tmp_3137;
   result += (0.25) * tmp_3134;

   return result;
}

std::complex<double> CMSSM_effective_couplings::CphhHpmconjHpm(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);

   std::complex<double> result;

   result = 0.25*(-(ZH(gt1,0)*(ZP(gt2,0)*(vd*(0.6*Sqr(g1) + Sqr(g2))*ZP(gt3,0)
      + vu*Sqr(g2)*ZP(gt3,1)) + ZP(gt2,1)*(vu*Sqr(g2)*ZP(gt3,0) + vd*(-0.6*Sqr(g1)
      + Sqr(g2))*ZP(gt3,1)))) + ZH(gt1,1)*(ZP(gt2,0)*(vu*(0.6*Sqr(g1) - Sqr(g2))*
      ZP(gt3,0) - vd*Sqr(g2)*ZP(gt3,1)) - ZP(gt2,1)*(vd*Sqr(g2)*ZP(gt3,0) + vu*(
      0.6*Sqr(g1) + Sqr(g2))*ZP(gt3,1))));

   return result;
}

std::complex<double> CMSSM_effective_couplings::CpChahhbarChaPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g2 = MODELPARAMETER(g2);

   std::complex<double> result;

   result = -0.7071067811865475*g2*(Conj(UM(gt1,1))*Conj(UP(gt3,0))*ZH(gt2,0) +
      Conj(UM(gt1,0))*Conj(UP(gt3,1))*ZH(gt2,1));

   return result;
}

std::complex<double> CMSSM_effective_couplings::CpFehhbarFePL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Ye = MODELPARAMETER(Ye);

   std::complex<double> result;

   std::complex<double> tmp_3168;
   std::complex<double> tmp_3169;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3170;
      std::complex<double> tmp_3171;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3171 += Conj(ZER(gt3,j1))*Ye(j1,j2);
      }
      tmp_3170 += tmp_3171;
      tmp_3169 += (Conj(ZEL(gt1,j2))) * tmp_3170;
   }
   tmp_3168 += tmp_3169;
   result += (-0.7071067811865475*ZH(gt2,0)) * tmp_3168;

   return result;
}

std::complex<double> CMSSM_effective_couplings::CpFdhhbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_3172;
   std::complex<double> tmp_3173;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3174;
      std::complex<double> tmp_3175;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3175 += Conj(ZDR(gt3,j1))*Yd(j1,j2);
      }
      tmp_3174 += tmp_3175;
      tmp_3173 += (Conj(ZDL(gt1,j2))) * tmp_3174;
   }
   tmp_3172 += tmp_3173;
   result += (-0.7071067811865475*ZH(gt2,0)) * tmp_3172;

   return result;
}

std::complex<double> CMSSM_effective_couplings::CpFuhhbarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_3176;
   std::complex<double> tmp_3177;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3178;
      std::complex<double> tmp_3179;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3179 += Conj(ZUR(gt3,j1))*Yu(j1,j2);
      }
      tmp_3178 += tmp_3179;
      tmp_3177 += (Conj(ZUL(gt1,j2))) * tmp_3178;
   }
   tmp_3176 += tmp_3177;
   result += (-0.7071067811865475*ZH(gt2,1)) * tmp_3176;

   return result;
}

std::complex<double> CMSSM_effective_couplings::CphhVWmconjVWm(unsigned gt1) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);

   std::complex<double> result;

   result = 0.5*Sqr(g2)*(vd*ZH(gt1,0) + vu*ZH(gt1,1));

   return result;
}

std::complex<double> CMSSM_effective_couplings::CpAhSdconjSd(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYd = MODELPARAMETER(TYd);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3180;
   std::complex<double> tmp_3181;
   std::complex<double> tmp_3182;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3183;
      std::complex<double> tmp_3184;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3184 += ZD(gt3,3 + j1)*TYd(j1,j2);
      }
      tmp_3183 += tmp_3184;
      tmp_3182 += (Conj(ZD(gt2,j2))) * tmp_3183;
   }
   tmp_3181 += tmp_3182;
   tmp_3180 += (ZA(gt1,0)) * tmp_3181;
   std::complex<double> tmp_3185;
   std::complex<double> tmp_3186;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3187;
      std::complex<double> tmp_3188;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3188 += Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_3187 += tmp_3188;
      tmp_3186 += (ZD(gt3,j2)) * tmp_3187;
   }
   tmp_3185 += tmp_3186;
   tmp_3180 += (-ZA(gt1,0)) * tmp_3185;
   std::complex<double> tmp_3189;
   std::complex<double> tmp_3190;
   std::complex<double> tmp_3191;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3192;
      std::complex<double> tmp_3193;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3193 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3192 += tmp_3193;
      tmp_3191 += (Conj(ZD(gt2,j2))) * tmp_3192;
   }
   tmp_3190 += tmp_3191;
   tmp_3189 += (Conj(Mu)) * tmp_3190;
   std::complex<double> tmp_3194;
   std::complex<double> tmp_3195;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3196;
      std::complex<double> tmp_3197;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3197 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3196 += tmp_3197;
      tmp_3195 += (ZD(gt3,j2)) * tmp_3196;
   }
   tmp_3194 += tmp_3195;
   tmp_3189 += (-Mu) * tmp_3194;
   tmp_3180 += (ZA(gt1,1)) * tmp_3189;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_3180;

   return result;
}

std::complex<double> CMSSM_effective_couplings::CpAhSuconjSu(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYu = MODELPARAMETER(TYu);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3198;
   std::complex<double> tmp_3199;
   std::complex<double> tmp_3200;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3201;
      std::complex<double> tmp_3202;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3202 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_3201 += tmp_3202;
      tmp_3200 += (Conj(ZU(gt2,j2))) * tmp_3201;
   }
   tmp_3199 += tmp_3200;
   tmp_3198 += (Conj(Mu)*ZA(gt1,0)) * tmp_3199;
   std::complex<double> tmp_3203;
   std::complex<double> tmp_3204;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3205;
      std::complex<double> tmp_3206;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3206 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3205 += tmp_3206;
      tmp_3204 += (ZU(gt3,j2)) * tmp_3205;
   }
   tmp_3203 += tmp_3204;
   tmp_3198 += (-(Mu*ZA(gt1,0))) * tmp_3203;
   std::complex<double> tmp_3207;
   std::complex<double> tmp_3208;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3209;
      std::complex<double> tmp_3210;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3210 += ZU(gt3,3 + j1)*TYu(j1,j2);
      }
      tmp_3209 += tmp_3210;
      tmp_3208 += (Conj(ZU(gt2,j2))) * tmp_3209;
   }
   tmp_3207 += tmp_3208;
   std::complex<double> tmp_3211;
   std::complex<double> tmp_3212;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3213;
      std::complex<double> tmp_3214;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3214 += Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_3213 += tmp_3214;
      tmp_3212 += (ZU(gt3,j2)) * tmp_3213;
   }
   tmp_3211 += tmp_3212;
   tmp_3207 += (-1) * tmp_3211;
   tmp_3198 += (ZA(gt1,1)) * tmp_3207;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_3198;

   return result;
}

std::complex<double> CMSSM_effective_couplings::CpAhSeconjSe(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYe = MODELPARAMETER(TYe);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3215;
   std::complex<double> tmp_3216;
   std::complex<double> tmp_3217;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3218;
      std::complex<double> tmp_3219;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3219 += ZE(gt3,3 + j1)*TYe(j1,j2);
      }
      tmp_3218 += tmp_3219;
      tmp_3217 += (Conj(ZE(gt2,j2))) * tmp_3218;
   }
   tmp_3216 += tmp_3217;
   tmp_3215 += (ZA(gt1,0)) * tmp_3216;
   std::complex<double> tmp_3220;
   std::complex<double> tmp_3221;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3222;
      std::complex<double> tmp_3223;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3223 += Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_3222 += tmp_3223;
      tmp_3221 += (ZE(gt3,j2)) * tmp_3222;
   }
   tmp_3220 += tmp_3221;
   tmp_3215 += (-ZA(gt1,0)) * tmp_3220;
   std::complex<double> tmp_3224;
   std::complex<double> tmp_3225;
   std::complex<double> tmp_3226;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3227;
      std::complex<double> tmp_3228;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3228 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3227 += tmp_3228;
      tmp_3226 += (Conj(ZE(gt2,j2))) * tmp_3227;
   }
   tmp_3225 += tmp_3226;
   tmp_3224 += (Conj(Mu)) * tmp_3225;
   std::complex<double> tmp_3229;
   std::complex<double> tmp_3230;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3231;
      std::complex<double> tmp_3232;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3232 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3231 += tmp_3232;
      tmp_3230 += (ZE(gt3,j2)) * tmp_3231;
   }
   tmp_3229 += tmp_3230;
   tmp_3224 += (-Mu) * tmp_3229;
   tmp_3215 += (ZA(gt1,1)) * tmp_3224;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_3215;

   return result;
}

std::complex<double> CMSSM_effective_couplings::CpAhHpmconjHpm(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);

   std::complex<double> result;

   result = std::complex<double>(0,-0.25)*Sqr(g2)*(vu*ZA(gt1,0) + vd*ZA(gt1,1))
      *(ZP(gt2,1)*ZP(gt3,0) - ZP(gt2,0)*ZP(gt3,1));

   return result;
}

std::complex<double> CMSSM_effective_couplings::CpAhChabarChaPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g2 = MODELPARAMETER(g2);

   std::complex<double> result;

   result = std::complex<double>(0.,0.7071067811865475)*g2*(Conj(UM(gt2,1))*
      Conj(UP(gt3,0))*ZA(gt1,0) + Conj(UM(gt2,0))*Conj(UP(gt3,1))*ZA(gt1,1));

   return result;
}

std::complex<double> CMSSM_effective_couplings::CpAhFebarFePL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Ye = MODELPARAMETER(Ye);

   std::complex<double> result;

   std::complex<double> tmp_3233;
   std::complex<double> tmp_3234;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3235;
      std::complex<double> tmp_3236;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3236 += Conj(ZER(gt3,j1))*Ye(j1,j2);
      }
      tmp_3235 += tmp_3236;
      tmp_3234 += (Conj(ZEL(gt2,j2))) * tmp_3235;
   }
   tmp_3233 += tmp_3234;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,0)) *
      tmp_3233;

   return result;
}

std::complex<double> CMSSM_effective_couplings::CpAhFdbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_3237;
   std::complex<double> tmp_3238;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3239;
      std::complex<double> tmp_3240;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3240 += Conj(ZDR(gt3,j1))*Yd(j1,j2);
      }
      tmp_3239 += tmp_3240;
      tmp_3238 += (Conj(ZDL(gt2,j2))) * tmp_3239;
   }
   tmp_3237 += tmp_3238;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,0)) *
      tmp_3237;

   return result;
}

std::complex<double> CMSSM_effective_couplings::CpAhFubarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_3241;
   std::complex<double> tmp_3242;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3243;
      std::complex<double> tmp_3244;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3244 += Conj(ZUR(gt3,j1))*Yu(j1,j2);
      }
      tmp_3243 += tmp_3244;
      tmp_3242 += (Conj(ZUL(gt2,j2))) * tmp_3243;
   }
   tmp_3241 += tmp_3242;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,1)) *
      tmp_3241;

   return result;
}

void CMSSM_effective_couplings::calculate_eff_CphhVPVP(unsigned gO1)
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

void CMSSM_effective_couplings::calculate_eff_CphhVGVG(unsigned gO1)
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

void CMSSM_effective_couplings::calculate_eff_CpAhVPVP(unsigned gO1)
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

void CMSSM_effective_couplings::calculate_eff_CpAhVGVG(unsigned gO1)
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
