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

// File generated at Wed 12 Apr 2017 13:01:39

#include "UMSSM_effective_couplings.hpp"

#include "effective_couplings.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

using namespace effective_couplings;

#define INPUTPARAMETER(parameter) model.get_input().parameter
#define MODELPARAMETER(parameter) model.get_##parameter()
#define DERIVEDPARAMETER(parameter) model.##parameter()
#define PHASE(parameter) model.get_##parameter()
#define PHYSICAL(parameter) model.get_physical().parameter

UMSSM_effective_couplings::UMSSM_effective_couplings(
   const UMSSM_mass_eigenstates& model_,
   const softsusy::QedQcd& qedqcd_,
   const Physical_input& input_)
   : model(model_), qedqcd(qedqcd_), physical_input(input_)
   , rg_improve(true), include_qcd_corrections(true)
   , ZD(MODELPARAMETER(ZD)), ZV(MODELPARAMETER(ZV)), ZU(MODELPARAMETER(ZU)), ZE
      (MODELPARAMETER(ZE)), ZH(MODELPARAMETER(ZH)), ZA(MODELPARAMETER(ZA)), ZP(
      MODELPARAMETER(ZP)), ZN(MODELPARAMETER(ZN)), ZVL(MODELPARAMETER(ZVL)), ZVR(
      MODELPARAMETER(ZVR)), UM(MODELPARAMETER(UM)), UP(MODELPARAMETER(UP)), ZEL(
      MODELPARAMETER(ZEL)), ZER(MODELPARAMETER(ZER)), ZDL(MODELPARAMETER(ZDL)),
      ZDR(MODELPARAMETER(ZDR)), ZUL(MODELPARAMETER(ZUL)), ZUR(MODELPARAMETER(ZUR))
      , ZZ(MODELPARAMETER(ZZ))

   , eff_CphhVPVP(Eigen::Array<std::complex<double>,3,1>::Zero()), eff_CphhVGVG
      (Eigen::Array<std::complex<double>,3,1>::Zero()), eff_CpAhVPVP(Eigen::Array<
      std::complex<double>,3,1>::Zero()), eff_CpAhVGVG(Eigen::Array<std::complex<
      double>,3,1>::Zero())

{
}

void UMSSM_effective_couplings::calculate_effective_couplings()
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

void UMSSM_effective_couplings::set_model(const UMSSM_mass_eigenstates& model_)
{
   model = model_;
   copy_mixing_matrices_from_model();
}

void UMSSM_effective_couplings::copy_mixing_matrices_from_model()
{
   ZD = MODELPARAMETER(ZD);
   ZV = MODELPARAMETER(ZV);
   ZU = MODELPARAMETER(ZU);
   ZE = MODELPARAMETER(ZE);
   ZH = MODELPARAMETER(ZH);
   ZA = MODELPARAMETER(ZA);
   ZP = MODELPARAMETER(ZP);
   ZN = MODELPARAMETER(ZN);
   ZVL = MODELPARAMETER(ZVL);
   ZVR = MODELPARAMETER(ZVR);
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

standard_model::Standard_model UMSSM_effective_couplings::initialise_SM() const
{
   standard_model::Standard_model sm;

   sm.set_loops(2);
   sm.set_thresholds(2);
   sm.set_physical_input(physical_input);

   sm.initialise_from_input(qedqcd);

   return sm;
}

void UMSSM_effective_couplings::run_SM_strong_coupling_to(standard_model::Standard_model sm, double m)
{
   sm.run_to(m);

   model.set_g3(sm.get_g3());

}

std::complex<double> UMSSM_effective_couplings::scalar_scalar_qcd_factor(double m_decay, double m_loop) const
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

std::complex<double> UMSSM_effective_couplings::scalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) *
         scalar_diphoton_fermion_loop(m_decay, m_loop);

   }

   return result;
}

std::complex<double> UMSSM_effective_couplings::pseudoscalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) *
         pseudoscalar_diphoton_fermion_loop(m_decay, m_loop);

   }

   return result;
}

double UMSSM_effective_couplings::number_of_active_flavours(double m) const
{
   if (m < qedqcd.displayMbMb()) {
      return 4.0;
   } else if (m < qedqcd.displayPoleMt()) {
      return 5.0;
   } else {
      return 6.0;
   }
}

double UMSSM_effective_couplings::scalar_scaling_factor(double m) const
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

double UMSSM_effective_couplings::pseudoscalar_scaling_factor(double m) const
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

double UMSSM_effective_couplings::get_hhVPVP_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(Mhh)(gO1);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CphhVPVP(gO1));
}

double UMSSM_effective_couplings::get_hhVGVG_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(Mhh)(gO1);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CphhVGVG(gO1));
}

double UMSSM_effective_couplings::get_AhVPVP_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(MAh)(gO1);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CpAhVPVP(gO1));
}

double UMSSM_effective_couplings::get_AhVGVG_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(MAh)(gO1);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CpAhVGVG(gO1));
}

std::complex<double> UMSSM_effective_couplings::CphhSdconjSd(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Qd = INPUTPARAMETER(Qd);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qs = INPUTPARAMETER(Qs);
   const auto TYd = MODELPARAMETER(TYd);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto vd = MODELPARAMETER(vd);
   const auto vS = MODELPARAMETER(vS);
   const auto vu = MODELPARAMETER(vu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   std::complex<double> result;

   std::complex<double> tmp_4973;
   std::complex<double> tmp_4974;
   std::complex<double> tmp_4975;
   std::complex<double> tmp_4976;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4976 += Conj(ZD(gt2,j1))*ZD(gt3,j1);
   }
   tmp_4975 += tmp_4976;
   tmp_4974 += (-2*Qq*Qs*vS*Sqr(gp)) * tmp_4975;
   std::complex<double> tmp_4977;
   std::complex<double> tmp_4978;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4978 += Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1);
   }
   tmp_4977 += tmp_4978;
   tmp_4974 += (-2*Qd*Qs*vS*Sqr(gp)) * tmp_4977;
   std::complex<double> tmp_4979;
   std::complex<double> tmp_4980;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4981;
      std::complex<double> tmp_4982;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4982 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_4981 += tmp_4982;
      tmp_4980 += (Conj(ZD(gt2,j2))) * tmp_4981;
   }
   tmp_4979 += tmp_4980;
   tmp_4974 += (vu*Conj(Lambdax)) * tmp_4979;
   std::complex<double> tmp_4983;
   std::complex<double> tmp_4984;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4985;
      std::complex<double> tmp_4986;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4986 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_4985 += tmp_4986;
      tmp_4984 += (ZD(gt3,j2)) * tmp_4985;
   }
   tmp_4983 += tmp_4984;
   tmp_4974 += (vu*Lambdax) * tmp_4983;
   tmp_4973 += (6*Conj(ZH(gt1,2))) * tmp_4974;
   std::complex<double> tmp_4987;
   std::complex<double> tmp_4988;
   std::complex<double> tmp_4989;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4989 += Conj(ZD(gt2,j1))*ZD(gt3,j1);
   }
   tmp_4988 += tmp_4989;
   tmp_4987 += (vu*(0.6*Sqr(g1) + 3*(Sqr(g2) + 4*QHu*Qq*Sqr(gp)))) * tmp_4988;
   std::complex<double> tmp_4990;
   std::complex<double> tmp_4991;
   std::complex<double> tmp_4992;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4992 += Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1);
   }
   tmp_4991 += tmp_4992;
   tmp_4990 += (vu*(0.6*Sqr(g1) + 6*Qd*QHu*Sqr(gp))) * tmp_4991;
   std::complex<double> tmp_4993;
   std::complex<double> tmp_4994;
   std::complex<double> tmp_4995;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4996;
      std::complex<double> tmp_4997;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4997 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_4996 += tmp_4997;
      tmp_4995 += (Conj(ZD(gt2,j2))) * tmp_4996;
   }
   tmp_4994 += tmp_4995;
   tmp_4993 += (Conj(Lambdax)) * tmp_4994;
   std::complex<double> tmp_4998;
   std::complex<double> tmp_4999;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5000;
      std::complex<double> tmp_5001;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5001 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_5000 += tmp_5001;
      tmp_4999 += (ZD(gt3,j2)) * tmp_5000;
   }
   tmp_4998 += tmp_4999;
   tmp_4993 += (Lambdax) * tmp_4998;
   tmp_4990 += (-3*vS) * tmp_4993;
   tmp_4987 += (2) * tmp_4990;
   tmp_4973 += (-Conj(ZH(gt1,1))) * tmp_4987;
   std::complex<double> tmp_5002;
   std::complex<double> tmp_5003;
   std::complex<double> tmp_5004;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_5004 += Conj(ZD(gt2,j1))*ZD(gt3,j1);
   }
   tmp_5003 += tmp_5004;
   tmp_5002 += (vd*(0.6*Sqr(g1) + 3*(Sqr(g2) - 4*QHd*Qq*Sqr(gp)))) * tmp_5003;
   std::complex<double> tmp_5005;
   std::complex<double> tmp_5006;
   std::complex<double> tmp_5007;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_5007 += Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1);
   }
   tmp_5006 += tmp_5007;
   tmp_5005 += (vd*(0.6*Sqr(g1) - 6*Qd*QHd*Sqr(gp))) * tmp_5006;
   std::complex<double> tmp_5008;
   std::complex<double> tmp_5009;
   std::complex<double> tmp_5010;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5011;
      std::complex<double> tmp_5012;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5012 += ZD(gt3,3 + j1)*TYd(j1,j2);
      }
      tmp_5011 += tmp_5012;
      tmp_5010 += (Conj(ZD(gt2,j2))) * tmp_5011;
   }
   tmp_5009 += tmp_5010;
   tmp_5008 += (1.4142135623730951) * tmp_5009;
   std::complex<double> tmp_5013;
   std::complex<double> tmp_5014;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5015;
      std::complex<double> tmp_5016;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5016 += Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_5015 += tmp_5016;
      tmp_5014 += (ZD(gt3,j2)) * tmp_5015;
   }
   tmp_5013 += tmp_5014;
   tmp_5008 += (1.4142135623730951) * tmp_5013;
   std::complex<double> tmp_5017;
   std::complex<double> tmp_5018;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_5019;
      std::complex<double> tmp_5020;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_5021;
         std::complex<double> tmp_5022;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_5022 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_5021 += tmp_5022;
         tmp_5020 += (ZD(gt3,3 + j2)) * tmp_5021;
      }
      tmp_5019 += tmp_5020;
      tmp_5018 += (Conj(ZD(gt2,3 + j3))) * tmp_5019;
   }
   tmp_5017 += tmp_5018;
   std::complex<double> tmp_5023;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_5024;
      std::complex<double> tmp_5025;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_5026;
         std::complex<double> tmp_5027;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_5027 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_5026 += tmp_5027;
         tmp_5025 += (Conj(ZD(gt2,j2))) * tmp_5026;
      }
      tmp_5024 += tmp_5025;
      tmp_5023 += (ZD(gt3,j3)) * tmp_5024;
   }
   tmp_5017 += tmp_5023;
   tmp_5008 += (2*vd) * tmp_5017;
   tmp_5005 += (-3) * tmp_5008;
   tmp_5002 += (2) * tmp_5005;
   tmp_4973 += (Conj(ZH(gt1,0))) * tmp_5002;
   result += (0.08333333333333333) * tmp_4973;

   return result;
}

std::complex<double> UMSSM_effective_couplings::CphhSuconjSu(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qs = INPUTPARAMETER(Qs);
   const auto Qu = INPUTPARAMETER(Qu);
   const auto TYu = MODELPARAMETER(TYu);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto vd = MODELPARAMETER(vd);
   const auto vS = MODELPARAMETER(vS);
   const auto vu = MODELPARAMETER(vu);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   std::complex<double> result;

   std::complex<double> tmp_5028;
   std::complex<double> tmp_5029;
   std::complex<double> tmp_5030;
   std::complex<double> tmp_5031;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_5031 += Conj(ZU(gt2,j1))*ZU(gt3,j1);
   }
   tmp_5030 += tmp_5031;
   tmp_5029 += (-2*Qq*Qs*vS*Sqr(gp)) * tmp_5030;
   std::complex<double> tmp_5032;
   std::complex<double> tmp_5033;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_5033 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_5032 += tmp_5033;
   tmp_5029 += (-2*Qs*Qu*vS*Sqr(gp)) * tmp_5032;
   std::complex<double> tmp_5034;
   std::complex<double> tmp_5035;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5036;
      std::complex<double> tmp_5037;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5037 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_5036 += tmp_5037;
      tmp_5035 += (Conj(ZU(gt2,j2))) * tmp_5036;
   }
   tmp_5034 += tmp_5035;
   tmp_5029 += (vd*Conj(Lambdax)) * tmp_5034;
   std::complex<double> tmp_5038;
   std::complex<double> tmp_5039;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5040;
      std::complex<double> tmp_5041;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5041 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_5040 += tmp_5041;
      tmp_5039 += (ZU(gt3,j2)) * tmp_5040;
   }
   tmp_5038 += tmp_5039;
   tmp_5029 += (vd*Lambdax) * tmp_5038;
   tmp_5028 += (6*Conj(ZH(gt1,2))) * tmp_5029;
   std::complex<double> tmp_5042;
   std::complex<double> tmp_5043;
   std::complex<double> tmp_5044;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_5044 += Conj(ZU(gt2,j1))*ZU(gt3,j1);
   }
   tmp_5043 += tmp_5044;
   tmp_5042 += (vd*(0.6*Sqr(g1) - 3*(Sqr(g2) + 4*QHd*Qq*Sqr(gp)))) * tmp_5043;
   std::complex<double> tmp_5045;
   std::complex<double> tmp_5046;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_5046 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_5045 += tmp_5046;
   tmp_5042 += (-4*vd*(0.6*Sqr(g1) + 3*QHd*Qu*Sqr(gp))) * tmp_5045;
   std::complex<double> tmp_5047;
   std::complex<double> tmp_5048;
   std::complex<double> tmp_5049;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5050;
      std::complex<double> tmp_5051;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5051 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_5050 += tmp_5051;
      tmp_5049 += (Conj(ZU(gt2,j2))) * tmp_5050;
   }
   tmp_5048 += tmp_5049;
   tmp_5047 += (Conj(Lambdax)) * tmp_5048;
   std::complex<double> tmp_5052;
   std::complex<double> tmp_5053;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5054;
      std::complex<double> tmp_5055;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5055 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_5054 += tmp_5055;
      tmp_5053 += (ZU(gt3,j2)) * tmp_5054;
   }
   tmp_5052 += tmp_5053;
   tmp_5047 += (Lambdax) * tmp_5052;
   tmp_5042 += (6*vS) * tmp_5047;
   tmp_5028 += (Conj(ZH(gt1,0))) * tmp_5042;
   std::complex<double> tmp_5056;
   std::complex<double> tmp_5057;
   std::complex<double> tmp_5058;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_5058 += Conj(ZU(gt2,j1))*ZU(gt3,j1);
   }
   tmp_5057 += tmp_5058;
   tmp_5056 += (vu*(0.6*Sqr(g1) - 3*Sqr(g2) + 12*QHu*Qq*Sqr(gp))) * tmp_5057;
   std::complex<double> tmp_5059;
   std::complex<double> tmp_5060;
   std::complex<double> tmp_5061;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_5061 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_5060 += tmp_5061;
   tmp_5059 += (-2*vu*(0.6*Sqr(g1) - 3*QHu*Qu*Sqr(gp))) * tmp_5060;
   std::complex<double> tmp_5062;
   std::complex<double> tmp_5063;
   std::complex<double> tmp_5064;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5065;
      std::complex<double> tmp_5066;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5066 += ZU(gt3,3 + j1)*TYu(j1,j2);
      }
      tmp_5065 += tmp_5066;
      tmp_5064 += (Conj(ZU(gt2,j2))) * tmp_5065;
   }
   tmp_5063 += tmp_5064;
   tmp_5062 += (1.4142135623730951) * tmp_5063;
   std::complex<double> tmp_5067;
   std::complex<double> tmp_5068;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5069;
      std::complex<double> tmp_5070;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5070 += Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_5069 += tmp_5070;
      tmp_5068 += (ZU(gt3,j2)) * tmp_5069;
   }
   tmp_5067 += tmp_5068;
   tmp_5062 += (1.4142135623730951) * tmp_5067;
   std::complex<double> tmp_5071;
   std::complex<double> tmp_5072;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_5073;
      std::complex<double> tmp_5074;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_5075;
         std::complex<double> tmp_5076;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_5076 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_5075 += tmp_5076;
         tmp_5074 += (ZU(gt3,3 + j2)) * tmp_5075;
      }
      tmp_5073 += tmp_5074;
      tmp_5072 += (Conj(ZU(gt2,3 + j3))) * tmp_5073;
   }
   tmp_5071 += tmp_5072;
   std::complex<double> tmp_5077;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_5078;
      std::complex<double> tmp_5079;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_5080;
         std::complex<double> tmp_5081;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_5081 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_5080 += tmp_5081;
         tmp_5079 += (Conj(ZU(gt2,j2))) * tmp_5080;
      }
      tmp_5078 += tmp_5079;
      tmp_5077 += (ZU(gt3,j3)) * tmp_5078;
   }
   tmp_5071 += tmp_5077;
   tmp_5062 += (2*vu) * tmp_5071;
   tmp_5059 += (3) * tmp_5062;
   tmp_5056 += (2) * tmp_5059;
   tmp_5028 += (-Conj(ZH(gt1,1))) * tmp_5056;
   result += (0.08333333333333333) * tmp_5028;

   return result;
}

std::complex<double> UMSSM_effective_couplings::CphhSeconjSe(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Qe = INPUTPARAMETER(Qe);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qs = INPUTPARAMETER(Qs);
   const auto TYe = MODELPARAMETER(TYe);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto vd = MODELPARAMETER(vd);
   const auto vS = MODELPARAMETER(vS);
   const auto vu = MODELPARAMETER(vu);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   std::complex<double> result;

   std::complex<double> tmp_5082;
   std::complex<double> tmp_5083;
   std::complex<double> tmp_5084;
   std::complex<double> tmp_5085;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_5085 += Conj(ZE(gt2,j1))*ZE(gt3,j1);
   }
   tmp_5084 += tmp_5085;
   tmp_5083 += (-2*Ql*Qs*vS*Sqr(gp)) * tmp_5084;
   std::complex<double> tmp_5086;
   std::complex<double> tmp_5087;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_5087 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_5086 += tmp_5087;
   tmp_5083 += (-2*Qe*Qs*vS*Sqr(gp)) * tmp_5086;
   std::complex<double> tmp_5088;
   std::complex<double> tmp_5089;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5090;
      std::complex<double> tmp_5091;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5091 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_5090 += tmp_5091;
      tmp_5089 += (Conj(ZE(gt2,j2))) * tmp_5090;
   }
   tmp_5088 += tmp_5089;
   tmp_5083 += (vu*Conj(Lambdax)) * tmp_5088;
   std::complex<double> tmp_5092;
   std::complex<double> tmp_5093;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5094;
      std::complex<double> tmp_5095;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5095 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_5094 += tmp_5095;
      tmp_5093 += (ZE(gt3,j2)) * tmp_5094;
   }
   tmp_5092 += tmp_5093;
   tmp_5083 += (vu*Lambdax) * tmp_5092;
   tmp_5082 += (2*Conj(ZH(gt1,2))) * tmp_5083;
   std::complex<double> tmp_5096;
   std::complex<double> tmp_5097;
   std::complex<double> tmp_5098;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_5098 += Conj(ZE(gt2,j1))*ZE(gt3,j1);
   }
   tmp_5097 += tmp_5098;
   tmp_5096 += (vu*(0.6*Sqr(g1) - Sqr(g2) - 4*QHu*Ql*Sqr(gp))) * tmp_5097;
   std::complex<double> tmp_5099;
   std::complex<double> tmp_5100;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_5100 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_5099 += tmp_5100;
   tmp_5096 += (-2*vu*(0.6*Sqr(g1) + 2*Qe*QHu*Sqr(gp))) * tmp_5099;
   std::complex<double> tmp_5101;
   std::complex<double> tmp_5102;
   std::complex<double> tmp_5103;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5104;
      std::complex<double> tmp_5105;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5105 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_5104 += tmp_5105;
      tmp_5103 += (Conj(ZE(gt2,j2))) * tmp_5104;
   }
   tmp_5102 += tmp_5103;
   tmp_5101 += (Conj(Lambdax)) * tmp_5102;
   std::complex<double> tmp_5106;
   std::complex<double> tmp_5107;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5108;
      std::complex<double> tmp_5109;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5109 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_5108 += tmp_5109;
      tmp_5107 += (ZE(gt3,j2)) * tmp_5108;
   }
   tmp_5106 += tmp_5107;
   tmp_5101 += (Lambdax) * tmp_5106;
   tmp_5096 += (2*vS) * tmp_5101;
   tmp_5082 += (Conj(ZH(gt1,1))) * tmp_5096;
   std::complex<double> tmp_5110;
   std::complex<double> tmp_5111;
   std::complex<double> tmp_5112;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_5112 += Conj(ZE(gt2,j1))*ZE(gt3,j1);
   }
   tmp_5111 += tmp_5112;
   tmp_5110 += (vd*(0.6*Sqr(g1) - Sqr(g2) + 4*QHd*Ql*Sqr(gp))) * tmp_5111;
   std::complex<double> tmp_5113;
   std::complex<double> tmp_5114;
   std::complex<double> tmp_5115;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_5115 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_5114 += tmp_5115;
   tmp_5113 += (-(vd*(0.6*Sqr(g1) - 2*Qe*QHd*Sqr(gp)))) * tmp_5114;
   std::complex<double> tmp_5116;
   std::complex<double> tmp_5117;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5118;
      std::complex<double> tmp_5119;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5119 += ZE(gt3,3 + j1)*TYe(j1,j2);
      }
      tmp_5118 += tmp_5119;
      tmp_5117 += (Conj(ZE(gt2,j2))) * tmp_5118;
   }
   tmp_5116 += tmp_5117;
   tmp_5113 += (1.4142135623730951) * tmp_5116;
   std::complex<double> tmp_5120;
   std::complex<double> tmp_5121;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5122;
      std::complex<double> tmp_5123;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5123 += Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_5122 += tmp_5123;
      tmp_5121 += (ZE(gt3,j2)) * tmp_5122;
   }
   tmp_5120 += tmp_5121;
   tmp_5113 += (1.4142135623730951) * tmp_5120;
   std::complex<double> tmp_5124;
   std::complex<double> tmp_5125;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_5126;
      std::complex<double> tmp_5127;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_5128;
         std::complex<double> tmp_5129;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_5129 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_5128 += tmp_5129;
         tmp_5127 += (ZE(gt3,3 + j2)) * tmp_5128;
      }
      tmp_5126 += tmp_5127;
      tmp_5125 += (Conj(ZE(gt2,3 + j3))) * tmp_5126;
   }
   tmp_5124 += tmp_5125;
   tmp_5113 += (2*vd) * tmp_5124;
   std::complex<double> tmp_5130;
   std::complex<double> tmp_5131;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_5132;
      std::complex<double> tmp_5133;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_5134;
         std::complex<double> tmp_5135;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_5135 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_5134 += tmp_5135;
         tmp_5133 += (Conj(ZE(gt2,j2))) * tmp_5134;
      }
      tmp_5132 += tmp_5133;
      tmp_5131 += (ZE(gt3,j3)) * tmp_5132;
   }
   tmp_5130 += tmp_5131;
   tmp_5113 += (2*vd) * tmp_5130;
   tmp_5110 += (2) * tmp_5113;
   tmp_5082 += (-Conj(ZH(gt1,0))) * tmp_5110;
   result += (0.25) * tmp_5082;

   return result;
}

std::complex<double> UMSSM_effective_couplings::CphhHpmconjHpm(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto Qs = INPUTPARAMETER(Qs);
   const auto TLambdax = MODELPARAMETER(TLambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto vd = MODELPARAMETER(vd);
   const auto vS = MODELPARAMETER(vS);
   const auto vu = MODELPARAMETER(vu);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   std::complex<double> result;

   result = 0.25*(-(Conj(ZH(gt1,0))*(ZP(gt2,0)*(vd*(0.6*Sqr(g1) + Sqr(g2) + 4*
      Sqr(gp)*Sqr(QHd))*ZP(gt3,0) + vu*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gt3,1)) +
      ZP(gt2,1)*(vu*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gt3,0) + vd*(-0.6*Sqr(g1) +
      Sqr(g2) + 4*QHd*QHu*Sqr(gp))*ZP(gt3,1)))) + Conj(ZH(gt1,1))*(ZP(gt2,0)*(vu*
      (0.6*Sqr(g1) - Sqr(g2) - 4*QHd*QHu*Sqr(gp))*ZP(gt3,0) - vd*(-2*AbsSqr(
      Lambdax) + Sqr(g2))*ZP(gt3,1)) - ZP(gt2,1)*(vd*(-2*AbsSqr(Lambdax) + Sqr(g2)
      )*ZP(gt3,0) + vu*(0.6*Sqr(g1) + Sqr(g2) + 4*Sqr(gp)*Sqr(QHu))*ZP(gt3,1))) -
      2*Conj(ZH(gt1,2))*(ZP(gt2,1)*(1.4142135623730951*Conj(TLambdax)*ZP(gt3,0) +
      2*vS*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp))*ZP(gt3,1)) + ZP(gt2,0)*(2*vS*(AbsSqr
      (Lambdax) + QHd*Qs*Sqr(gp))*ZP(gt3,0) + 1.4142135623730951*TLambdax*ZP(gt3,1
      ))));

   return result;
}

std::complex<double> UMSSM_effective_couplings::CpChahhbarChaPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   std::complex<double> result;

   result = -0.7071067811865475*(g2*Conj(UM(gt1,0))*Conj(UP(gt3,1))*Conj(ZH(gt2
      ,1)) + Conj(UM(gt1,1))*(g2*Conj(UP(gt3,0))*Conj(ZH(gt2,0)) + Conj(UP(gt3,1))
      *Conj(ZH(gt2,2))*Lambdax));

   return result;
}

std::complex<double> UMSSM_effective_couplings::CpFehhbarFePL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Ye = MODELPARAMETER(Ye);

   std::complex<double> result;

   std::complex<double> tmp_5136;
   std::complex<double> tmp_5137;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5138;
      std::complex<double> tmp_5139;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5139 += Conj(ZER(gt3,j1))*Ye(j1,j2);
      }
      tmp_5138 += tmp_5139;
      tmp_5137 += (Conj(ZEL(gt1,j2))) * tmp_5138;
   }
   tmp_5136 += tmp_5137;
   result += (-0.7071067811865475*Conj(ZH(gt2,0))) * tmp_5136;

   return result;
}

std::complex<double> UMSSM_effective_couplings::CpFdhhbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_5140;
   std::complex<double> tmp_5141;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5142;
      std::complex<double> tmp_5143;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5143 += Conj(ZDR(gt3,j1))*Yd(j1,j2);
      }
      tmp_5142 += tmp_5143;
      tmp_5141 += (Conj(ZDL(gt1,j2))) * tmp_5142;
   }
   tmp_5140 += tmp_5141;
   result += (-0.7071067811865475*Conj(ZH(gt2,0))) * tmp_5140;

   return result;
}

std::complex<double> UMSSM_effective_couplings::CpFuhhbarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_5144;
   std::complex<double> tmp_5145;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5146;
      std::complex<double> tmp_5147;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5147 += Conj(ZUR(gt3,j1))*Yu(j1,j2);
      }
      tmp_5146 += tmp_5147;
      tmp_5145 += (Conj(ZUL(gt1,j2))) * tmp_5146;
   }
   tmp_5144 += tmp_5145;
   result += (-0.7071067811865475*Conj(ZH(gt2,1))) * tmp_5144;

   return result;
}

std::complex<double> UMSSM_effective_couplings::CphhVWmconjVWm(unsigned gt1) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);

   std::complex<double> result;

   result = 0.5*(vd*Conj(ZH(gt1,0)) + vu*Conj(ZH(gt1,1)))*Sqr(g2);

   return result;
}

std::complex<double> UMSSM_effective_couplings::CpAhSdconjSd(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYd = MODELPARAMETER(TYd);
   const auto vS = MODELPARAMETER(vS);
   const auto vu = MODELPARAMETER(vu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   std::complex<double> result;

   std::complex<double> tmp_5148;
   std::complex<double> tmp_5149;
   std::complex<double> tmp_5150;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5151;
      std::complex<double> tmp_5152;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5152 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_5151 += tmp_5152;
      tmp_5150 += (Conj(ZD(gt2,j2))) * tmp_5151;
   }
   tmp_5149 += tmp_5150;
   tmp_5148 += (Conj(Lambdax)*(vS*Conj(ZA(gt1,1)) + vu*Conj(ZA(gt1,2)))) *
      tmp_5149;
   std::complex<double> tmp_5153;
   std::complex<double> tmp_5154;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5155;
      std::complex<double> tmp_5156;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5156 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_5155 += tmp_5156;
      tmp_5154 += (ZD(gt3,j2)) * tmp_5155;
   }
   tmp_5153 += tmp_5154;
   tmp_5148 += (-((vS*Conj(ZA(gt1,1)) + vu*Conj(ZA(gt1,2)))*Lambdax)) *
      tmp_5153;
   std::complex<double> tmp_5157;
   std::complex<double> tmp_5158;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5159;
      std::complex<double> tmp_5160;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5160 += ZD(gt3,3 + j1)*TYd(j1,j2);
      }
      tmp_5159 += tmp_5160;
      tmp_5158 += (Conj(ZD(gt2,j2))) * tmp_5159;
   }
   tmp_5157 += tmp_5158;
   std::complex<double> tmp_5161;
   std::complex<double> tmp_5162;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5163;
      std::complex<double> tmp_5164;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5164 += Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_5163 += tmp_5164;
      tmp_5162 += (ZD(gt3,j2)) * tmp_5163;
   }
   tmp_5161 += tmp_5162;
   tmp_5157 += (-1) * tmp_5161;
   tmp_5148 += (1.4142135623730951*Conj(ZA(gt1,0))) * tmp_5157;
   result += (std::complex<double>(0,-0.5)) * tmp_5148;

   return result;
}

std::complex<double> UMSSM_effective_couplings::CpAhSuconjSu(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYu = MODELPARAMETER(TYu);
   const auto vd = MODELPARAMETER(vd);
   const auto vS = MODELPARAMETER(vS);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   std::complex<double> result;

   std::complex<double> tmp_5165;
   std::complex<double> tmp_5166;
   std::complex<double> tmp_5167;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5168;
      std::complex<double> tmp_5169;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5169 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_5168 += tmp_5169;
      tmp_5167 += (Conj(ZU(gt2,j2))) * tmp_5168;
   }
   tmp_5166 += tmp_5167;
   tmp_5165 += (Conj(Lambdax)*(vS*Conj(ZA(gt1,0)) + vd*Conj(ZA(gt1,2)))) *
      tmp_5166;
   std::complex<double> tmp_5170;
   std::complex<double> tmp_5171;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5172;
      std::complex<double> tmp_5173;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5173 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_5172 += tmp_5173;
      tmp_5171 += (ZU(gt3,j2)) * tmp_5172;
   }
   tmp_5170 += tmp_5171;
   tmp_5165 += (-((vS*Conj(ZA(gt1,0)) + vd*Conj(ZA(gt1,2)))*Lambdax)) *
      tmp_5170;
   std::complex<double> tmp_5174;
   std::complex<double> tmp_5175;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5176;
      std::complex<double> tmp_5177;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5177 += ZU(gt3,3 + j1)*TYu(j1,j2);
      }
      tmp_5176 += tmp_5177;
      tmp_5175 += (Conj(ZU(gt2,j2))) * tmp_5176;
   }
   tmp_5174 += tmp_5175;
   std::complex<double> tmp_5178;
   std::complex<double> tmp_5179;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5180;
      std::complex<double> tmp_5181;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5181 += Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_5180 += tmp_5181;
      tmp_5179 += (ZU(gt3,j2)) * tmp_5180;
   }
   tmp_5178 += tmp_5179;
   tmp_5174 += (-1) * tmp_5178;
   tmp_5165 += (1.4142135623730951*Conj(ZA(gt1,1))) * tmp_5174;
   result += (std::complex<double>(0,-0.5)) * tmp_5165;

   return result;
}

std::complex<double> UMSSM_effective_couplings::CpAhSeconjSe(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYe = MODELPARAMETER(TYe);
   const auto vS = MODELPARAMETER(vS);
   const auto vu = MODELPARAMETER(vu);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   std::complex<double> result;

   std::complex<double> tmp_5182;
   std::complex<double> tmp_5183;
   std::complex<double> tmp_5184;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5185;
      std::complex<double> tmp_5186;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5186 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_5185 += tmp_5186;
      tmp_5184 += (Conj(ZE(gt2,j2))) * tmp_5185;
   }
   tmp_5183 += tmp_5184;
   tmp_5182 += (Conj(Lambdax)*(vS*Conj(ZA(gt1,1)) + vu*Conj(ZA(gt1,2)))) *
      tmp_5183;
   std::complex<double> tmp_5187;
   std::complex<double> tmp_5188;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5189;
      std::complex<double> tmp_5190;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5190 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_5189 += tmp_5190;
      tmp_5188 += (ZE(gt3,j2)) * tmp_5189;
   }
   tmp_5187 += tmp_5188;
   tmp_5182 += (-((vS*Conj(ZA(gt1,1)) + vu*Conj(ZA(gt1,2)))*Lambdax)) *
      tmp_5187;
   std::complex<double> tmp_5191;
   std::complex<double> tmp_5192;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5193;
      std::complex<double> tmp_5194;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5194 += ZE(gt3,3 + j1)*TYe(j1,j2);
      }
      tmp_5193 += tmp_5194;
      tmp_5192 += (Conj(ZE(gt2,j2))) * tmp_5193;
   }
   tmp_5191 += tmp_5192;
   std::complex<double> tmp_5195;
   std::complex<double> tmp_5196;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5197;
      std::complex<double> tmp_5198;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5198 += Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_5197 += tmp_5198;
      tmp_5196 += (ZE(gt3,j2)) * tmp_5197;
   }
   tmp_5195 += tmp_5196;
   tmp_5191 += (-1) * tmp_5195;
   tmp_5182 += (1.4142135623730951*Conj(ZA(gt1,0))) * tmp_5191;
   result += (std::complex<double>(0,-0.5)) * tmp_5182;

   return result;
}

std::complex<double> UMSSM_effective_couplings::CpAhHpmconjHpm(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TLambdax = MODELPARAMETER(TLambdax);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   std::complex<double> result;

   result = std::complex<double>(0,-0.25)*(vu*Conj(ZA(gt1,0))*(-2*AbsSqr(
      Lambdax) + Sqr(g2))*(ZP(gt2,1)*ZP(gt3,0) - ZP(gt2,0)*ZP(gt3,1)) + vd*Conj(ZA
      (gt1,1))*(-2*AbsSqr(Lambdax) + Sqr(g2))*(ZP(gt2,1)*ZP(gt3,0) - ZP(gt2,0)*ZP(
      gt3,1)) + 2.8284271247461903*Conj(ZA(gt1,2))*(-(Conj(TLambdax)*ZP(gt2,1)*ZP(
      gt3,0)) + TLambdax*ZP(gt2,0)*ZP(gt3,1)));

   return result;
}

std::complex<double> UMSSM_effective_couplings::CpAhChabarChaPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   std::complex<double> result;

   result = std::complex<double>(0.,0.7071067811865475)*(g2*Conj(UM(gt2,0))*
      Conj(UP(gt3,1))*Conj(ZA(gt1,1)) + Conj(UM(gt2,1))*(g2*Conj(UP(gt3,0))*Conj(
      ZA(gt1,0)) - Conj(UP(gt3,1))*Conj(ZA(gt1,2))*Lambdax));

   return result;
}

std::complex<double> UMSSM_effective_couplings::CpAhFebarFePL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Ye = MODELPARAMETER(Ye);

   std::complex<double> result;

   std::complex<double> tmp_5199;
   std::complex<double> tmp_5200;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5201;
      std::complex<double> tmp_5202;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5202 += Conj(ZER(gt3,j1))*Ye(j1,j2);
      }
      tmp_5201 += tmp_5202;
      tmp_5200 += (Conj(ZEL(gt2,j2))) * tmp_5201;
   }
   tmp_5199 += tmp_5200;
   result += (std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gt1,0))) *
      tmp_5199;

   return result;
}

std::complex<double> UMSSM_effective_couplings::CpAhFdbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_5203;
   std::complex<double> tmp_5204;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5205;
      std::complex<double> tmp_5206;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5206 += Conj(ZDR(gt3,j1))*Yd(j1,j2);
      }
      tmp_5205 += tmp_5206;
      tmp_5204 += (Conj(ZDL(gt2,j2))) * tmp_5205;
   }
   tmp_5203 += tmp_5204;
   result += (std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gt1,0))) *
      tmp_5203;

   return result;
}

std::complex<double> UMSSM_effective_couplings::CpAhFubarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_5207;
   std::complex<double> tmp_5208;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_5209;
      std::complex<double> tmp_5210;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5210 += Conj(ZUR(gt3,j1))*Yu(j1,j2);
      }
      tmp_5209 += tmp_5210;
      tmp_5208 += (Conj(ZUL(gt2,j2))) * tmp_5209;
   }
   tmp_5207 += tmp_5208;
   result += (std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gt1,1))) *
      tmp_5207;

   return result;
}

void UMSSM_effective_couplings::calculate_eff_CphhVPVP(unsigned gO1)
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

void UMSSM_effective_couplings::calculate_eff_CphhVGVG(unsigned gO1)
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

void UMSSM_effective_couplings::calculate_eff_CpAhVPVP(unsigned gO1)
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

void UMSSM_effective_couplings::calculate_eff_CpAhVGVG(unsigned gO1)
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
