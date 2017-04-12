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

// File generated at Wed 12 Apr 2017 12:19:04

#include "MRSSM_effective_couplings.hpp"

#include "effective_couplings.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

using namespace effective_couplings;

#define INPUTPARAMETER(parameter) model.get_input().parameter
#define MODELPARAMETER(parameter) model.get_##parameter()
#define DERIVEDPARAMETER(parameter) model.##parameter()
#define PHASE(parameter) model.get_##parameter()
#define PHYSICAL(parameter) model.get_physical().parameter

MRSSM_effective_couplings::MRSSM_effective_couplings(
   const MRSSM_mass_eigenstates& model_,
   const softsusy::QedQcd& qedqcd_,
   const Physical_input& input_)
   : model(model_), qedqcd(qedqcd_), physical_input(input_)
   , rg_improve(true), include_qcd_corrections(true)
   , ZD(MODELPARAMETER(ZD)), ZV(MODELPARAMETER(ZV)), ZU(MODELPARAMETER(ZU)), ZE
      (MODELPARAMETER(ZE)), ZH(MODELPARAMETER(ZH)), ZA(MODELPARAMETER(ZA)), ZHR(
      MODELPARAMETER(ZHR)), ZP(MODELPARAMETER(ZP)), ZN1(MODELPARAMETER(ZN1)), ZN2(
      MODELPARAMETER(ZN2)), UM1(MODELPARAMETER(UM1)), UP1(MODELPARAMETER(UP1)),
      UM2(MODELPARAMETER(UM2)), UP2(MODELPARAMETER(UP2)), ZEL(MODELPARAMETER(ZEL))
      , ZER(MODELPARAMETER(ZER)), ZDL(MODELPARAMETER(ZDL)), ZDR(MODELPARAMETER(ZDR
      )), ZUL(MODELPARAMETER(ZUL)), ZUR(MODELPARAMETER(ZUR)), ZZ(MODELPARAMETER(ZZ
      ))

   , eff_CphhVPVP(Eigen::Array<std::complex<double>,4,1>::Zero()), eff_CphhVGVG
      (Eigen::Array<std::complex<double>,4,1>::Zero()), eff_CpAhVPVP(Eigen::Array<
      std::complex<double>,4,1>::Zero()), eff_CpAhVGVG(Eigen::Array<std::complex<
      double>,4,1>::Zero())

{
}

void MRSSM_effective_couplings::calculate_effective_couplings()
{
   const standard_model::Standard_model sm(initialise_SM());

   const double scale = model.get_scale();
   const Eigen::ArrayXd saved_parameters(model.get());

   const double saved_mt = PHYSICAL(MFu(2));
   PHYSICAL(MFu(2)) = qedqcd.displayPoleMt();

   const auto Mhh = PHYSICAL(Mhh);
   for (unsigned gO1 = 0; gO1 < 4; ++gO1) {
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
   for (unsigned gO1 = 1; gO1 < 4; ++gO1) {
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

void MRSSM_effective_couplings::set_model(const MRSSM_mass_eigenstates& model_)
{
   model = model_;
   copy_mixing_matrices_from_model();
}

void MRSSM_effective_couplings::copy_mixing_matrices_from_model()
{
   ZD = MODELPARAMETER(ZD);
   ZV = MODELPARAMETER(ZV);
   ZU = MODELPARAMETER(ZU);
   ZE = MODELPARAMETER(ZE);
   ZH = MODELPARAMETER(ZH);
   ZA = MODELPARAMETER(ZA);
   ZHR = MODELPARAMETER(ZHR);
   ZP = MODELPARAMETER(ZP);
   ZN1 = MODELPARAMETER(ZN1);
   ZN2 = MODELPARAMETER(ZN2);
   UM1 = MODELPARAMETER(UM1);
   UP1 = MODELPARAMETER(UP1);
   UM2 = MODELPARAMETER(UM2);
   UP2 = MODELPARAMETER(UP2);
   ZEL = MODELPARAMETER(ZEL);
   ZER = MODELPARAMETER(ZER);
   ZDL = MODELPARAMETER(ZDL);
   ZDR = MODELPARAMETER(ZDR);
   ZUL = MODELPARAMETER(ZUL);
   ZUR = MODELPARAMETER(ZUR);
   ZZ = MODELPARAMETER(ZZ);

}

standard_model::Standard_model MRSSM_effective_couplings::initialise_SM() const
{
   standard_model::Standard_model sm;

   sm.set_loops(2);
   sm.set_thresholds(2);
   sm.set_physical_input(physical_input);

   sm.initialise_from_input(qedqcd);

   return sm;
}

void MRSSM_effective_couplings::run_SM_strong_coupling_to(standard_model::Standard_model sm, double m)
{
   sm.run_to(m);

   model.set_g3(sm.get_g3());

}

std::complex<double> MRSSM_effective_couplings::scalar_scalar_qcd_factor(double m_decay, double m_loop) const
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

std::complex<double> MRSSM_effective_couplings::scalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) *
         scalar_diphoton_fermion_loop(m_decay, m_loop);

   }

   return result;
}

std::complex<double> MRSSM_effective_couplings::pseudoscalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) *
         pseudoscalar_diphoton_fermion_loop(m_decay, m_loop);

   }

   return result;
}

double MRSSM_effective_couplings::number_of_active_flavours(double m) const
{
   if (m < qedqcd.displayMbMb()) {
      return 4.0;
   } else if (m < qedqcd.displayPoleMt()) {
      return 5.0;
   } else {
      return 6.0;
   }
}

double MRSSM_effective_couplings::scalar_scaling_factor(double m) const
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

double MRSSM_effective_couplings::pseudoscalar_scaling_factor(double m) const
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

double MRSSM_effective_couplings::get_hhVPVP_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(Mhh)(gO1);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CphhVPVP(gO1));
}

double MRSSM_effective_couplings::get_hhVGVG_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(Mhh)(gO1);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CphhVGVG(gO1));
}

double MRSSM_effective_couplings::get_AhVPVP_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(MAh)(gO1);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CpAhVPVP(gO1));
}

double MRSSM_effective_couplings::get_AhVGVG_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(MAh)(gO1);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CpAhVGVG(gO1));
}

std::complex<double> MRSSM_effective_couplings::CphhSRdpconjSRdp(unsigned gt1) const
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto MuD = MODELPARAMETER(MuD);
   const auto vd = MODELPARAMETER(vd);
   const auto vS = MODELPARAMETER(vS);
   const auto vT = MODELPARAMETER(vT);
   const auto vu = MODELPARAMETER(vu);

   std::complex<double> result;

   result = 0.25*(vd*(-4*AbsSqr(LamTD) + 0.6*Sqr(g1) - Sqr(g2))*ZH(gt1,0) + vu*
      (-0.6*Sqr(g1) + Sqr(g2))*ZH(gt1,1) - 1.5491933384829668*g1*MDBS*ZH(gt1,2) -
      4*vS*AbsSqr(LamSD)*ZH(gt1,2) - 2.8284271247461903*MuD*Conj(LamSD)*ZH(gt1,2)
      + 1.4142135623730951*LamTD*vT*Conj(LamSD)*ZH(gt1,2) + 1.4142135623730951*
      LamSD*vT*Conj(LamTD)*ZH(gt1,2) - 1.5491933384829668*g1*Conj(MDBS)*ZH(gt1,2)
      - 2.8284271247461903*LamSD*Conj(MuD)*ZH(gt1,2) - 2*g2*MDWBT*ZH(gt1,3) - 2*vT
      *AbsSqr(LamTD)*ZH(gt1,3) + 1.4142135623730951*LamTD*vS*Conj(LamSD)*ZH(gt1,3)
      + 2*MuD*Conj(LamTD)*ZH(gt1,3) + 1.4142135623730951*LamSD*vS*Conj(LamTD)*ZH(
      gt1,3) - 2*g2*Conj(MDWBT)*ZH(gt1,3) + 2*LamTD*Conj(MuD)*ZH(gt1,3));

   return result;
}

std::complex<double> MRSSM_effective_couplings::CphhSRumconjSRum(unsigned gt1) const
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto MuU = MODELPARAMETER(MuU);
   const auto vd = MODELPARAMETER(vd);
   const auto vS = MODELPARAMETER(vS);
   const auto vT = MODELPARAMETER(vT);
   const auto vu = MODELPARAMETER(vu);

   std::complex<double> result;

   result = 0.25*(vd*(-0.6*Sqr(g1) + Sqr(g2))*ZH(gt1,0) + vu*(-4*AbsSqr(LamTU)
      + 0.6*Sqr(g1) - Sqr(g2))*ZH(gt1,1) + 1.5491933384829668*g1*MDBS*ZH(gt1,2) -
      4*vS*AbsSqr(LamSU)*ZH(gt1,2) - 2.8284271247461903*MuU*Conj(LamSU)*ZH(gt1,2)
      - 1.4142135623730951*LamTU*vT*Conj(LamSU)*ZH(gt1,2) - 1.4142135623730951*
      LamSU*vT*Conj(LamTU)*ZH(gt1,2) + 1.5491933384829668*g1*Conj(MDBS)*ZH(gt1,2)
      - 2.8284271247461903*LamSU*Conj(MuU)*ZH(gt1,2) + 2*g2*MDWBT*ZH(gt1,3) - 2*vT
      *AbsSqr(LamTU)*ZH(gt1,3) - 1.4142135623730951*LamTU*vS*Conj(LamSU)*ZH(gt1,3)
      - 2*MuU*Conj(LamTU)*ZH(gt1,3) - 1.4142135623730951*LamSU*vS*Conj(LamTU)*ZH(
      gt1,3) + 2*g2*Conj(MDWBT)*ZH(gt1,3) - 2*LamTU*Conj(MuU)*ZH(gt1,3));

   return result;
}

std::complex<double> MRSSM_effective_couplings::CphhSdconjSd(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3916;
   std::complex<double> tmp_3917;
   std::complex<double> tmp_3918;
   std::complex<double> tmp_3919;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3920;
      std::complex<double> tmp_3921;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3922;
         std::complex<double> tmp_3923;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3923 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_3922 += tmp_3923;
         tmp_3921 += (ZD(gt3,3 + j2)) * tmp_3922;
      }
      tmp_3920 += tmp_3921;
      tmp_3919 += (Conj(ZD(gt2,3 + j3))) * tmp_3920;
   }
   tmp_3918 += tmp_3919;
   tmp_3917 += (-6*vd*ZH(gt1,0)) * tmp_3918;
   std::complex<double> tmp_3924;
   std::complex<double> tmp_3925;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3926;
      std::complex<double> tmp_3927;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3928;
         std::complex<double> tmp_3929;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3929 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_3928 += tmp_3929;
         tmp_3927 += (Conj(ZD(gt2,j2))) * tmp_3928;
      }
      tmp_3926 += tmp_3927;
      tmp_3925 += (ZD(gt3,j3)) * tmp_3926;
   }
   tmp_3924 += tmp_3925;
   tmp_3917 += (-6*vd*ZH(gt1,0)) * tmp_3924;
   std::complex<double> tmp_3930;
   std::complex<double> tmp_3931;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3932;
      std::complex<double> tmp_3933;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3933 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3932 += tmp_3933;
      tmp_3931 += (Conj(ZD(gt2,j2))) * tmp_3932;
   }
   tmp_3930 += tmp_3931;
   tmp_3917 += (4.242640687119286*Conj(Mu)*ZH(gt1,1)) * tmp_3930;
   std::complex<double> tmp_3934;
   std::complex<double> tmp_3935;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3936;
      std::complex<double> tmp_3937;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3937 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3936 += tmp_3937;
      tmp_3935 += (ZD(gt3,j2)) * tmp_3936;
   }
   tmp_3934 += tmp_3935;
   tmp_3917 += (4.242640687119286*Mu*ZH(gt1,1)) * tmp_3934;
   std::complex<double> tmp_3938;
   std::complex<double> tmp_3939;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3939 += Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1);
   }
   tmp_3938 += tmp_3939;
   tmp_3917 += (0.7745966692414834*g1*(0.7745966692414834*g1*vd*ZH(gt1,0) -
      0.7745966692414834*g1*vu*ZH(gt1,1) - 2*(MDBS + Conj(MDBS))*ZH(gt1,2))) *
      tmp_3938;
   tmp_3916 += (2) * tmp_3917;
   std::complex<double> tmp_3940;
   std::complex<double> tmp_3941;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3941 += Conj(ZD(gt2,j1))*ZD(gt3,j1);
   }
   tmp_3940 += tmp_3941;
   tmp_3916 += (vd*(0.6*Sqr(g1) + 3*Sqr(g2))*ZH(gt1,0) - vu*(0.6*Sqr(g1) + 3*
      Sqr(g2))*ZH(gt1,1) - 1.5491933384829668*g1*MDBS*ZH(gt1,2) -
      1.5491933384829668*g1*Conj(MDBS)*ZH(gt1,2) + 6*g2*MDWBT*ZH(gt1,3) + 6*g2*
      Conj(MDWBT)*ZH(gt1,3)) * tmp_3940;
   result += (0.08333333333333333) * tmp_3916;

   return result;
}

std::complex<double> MRSSM_effective_couplings::CphhSuconjSu(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3942;
   std::complex<double> tmp_3943;
   std::complex<double> tmp_3944;
   std::complex<double> tmp_3945;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3946;
      std::complex<double> tmp_3947;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3947 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_3946 += tmp_3947;
      tmp_3945 += (Conj(ZU(gt2,j2))) * tmp_3946;
   }
   tmp_3944 += tmp_3945;
   tmp_3943 += (1.4142135623730951*Conj(Mu)*ZH(gt1,0)) * tmp_3944;
   std::complex<double> tmp_3948;
   std::complex<double> tmp_3949;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3950;
      std::complex<double> tmp_3951;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3951 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3950 += tmp_3951;
      tmp_3949 += (ZU(gt3,j2)) * tmp_3950;
   }
   tmp_3948 += tmp_3949;
   tmp_3943 += (1.4142135623730951*Mu*ZH(gt1,0)) * tmp_3948;
   std::complex<double> tmp_3952;
   std::complex<double> tmp_3953;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3954;
      std::complex<double> tmp_3955;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3956;
         std::complex<double> tmp_3957;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3957 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_3956 += tmp_3957;
         tmp_3955 += (ZU(gt3,3 + j2)) * tmp_3956;
      }
      tmp_3954 += tmp_3955;
      tmp_3953 += (Conj(ZU(gt2,3 + j3))) * tmp_3954;
   }
   tmp_3952 += tmp_3953;
   std::complex<double> tmp_3958;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3959;
      std::complex<double> tmp_3960;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3961;
         std::complex<double> tmp_3962;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3962 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_3961 += tmp_3962;
         tmp_3960 += (Conj(ZU(gt2,j2))) * tmp_3961;
      }
      tmp_3959 += tmp_3960;
      tmp_3958 += (ZU(gt3,j3)) * tmp_3959;
   }
   tmp_3952 += tmp_3958;
   tmp_3943 += (-2*vu*ZH(gt1,1)) * tmp_3952;
   tmp_3942 += (6) * tmp_3943;
   std::complex<double> tmp_3963;
   std::complex<double> tmp_3964;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3964 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_3963 += tmp_3964;
   tmp_3942 += (3.0983866769659336*g1*(-0.7745966692414834*g1*vd*ZH(gt1,0) +
      0.7745966692414834*g1*vu*ZH(gt1,1) + 2*(MDBS + Conj(MDBS))*ZH(gt1,2))) *
      tmp_3963;
   std::complex<double> tmp_3965;
   std::complex<double> tmp_3966;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3966 += Conj(ZU(gt2,j1))*ZU(gt3,j1);
   }
   tmp_3965 += tmp_3966;
   tmp_3942 += (vd*(0.6*Sqr(g1) - 3*Sqr(g2))*ZH(gt1,0) - vu*(0.6*Sqr(g1) - 3*
      Sqr(g2))*ZH(gt1,1) - 2*(0.7745966692414834*g1*(MDBS + Conj(MDBS))*ZH(gt1,2)
      + 3*g2*(MDWBT + Conj(MDWBT))*ZH(gt1,3))) * tmp_3965;
   result += (0.08333333333333333) * tmp_3942;

   return result;
}

std::complex<double> MRSSM_effective_couplings::CphhSeconjSe(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3967;
   std::complex<double> tmp_3968;
   std::complex<double> tmp_3969;
   std::complex<double> tmp_3970;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3971;
      std::complex<double> tmp_3972;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3973;
         std::complex<double> tmp_3974;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3974 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_3973 += tmp_3974;
         tmp_3972 += (ZE(gt3,3 + j2)) * tmp_3973;
      }
      tmp_3971 += tmp_3972;
      tmp_3970 += (Conj(ZE(gt2,3 + j3))) * tmp_3971;
   }
   tmp_3969 += tmp_3970;
   tmp_3968 += (-2*vd*ZH(gt1,0)) * tmp_3969;
   std::complex<double> tmp_3975;
   std::complex<double> tmp_3976;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3977;
      std::complex<double> tmp_3978;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3979;
         std::complex<double> tmp_3980;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3980 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_3979 += tmp_3980;
         tmp_3978 += (Conj(ZE(gt2,j2))) * tmp_3979;
      }
      tmp_3977 += tmp_3978;
      tmp_3976 += (ZE(gt3,j3)) * tmp_3977;
   }
   tmp_3975 += tmp_3976;
   tmp_3968 += (-2*vd*ZH(gt1,0)) * tmp_3975;
   std::complex<double> tmp_3981;
   std::complex<double> tmp_3982;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3983;
      std::complex<double> tmp_3984;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3984 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3983 += tmp_3984;
      tmp_3982 += (Conj(ZE(gt2,j2))) * tmp_3983;
   }
   tmp_3981 += tmp_3982;
   tmp_3968 += (1.4142135623730951*Conj(Mu)*ZH(gt1,1)) * tmp_3981;
   std::complex<double> tmp_3985;
   std::complex<double> tmp_3986;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3987;
      std::complex<double> tmp_3988;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3988 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3987 += tmp_3988;
      tmp_3986 += (ZE(gt3,j2)) * tmp_3987;
   }
   tmp_3985 += tmp_3986;
   tmp_3968 += (1.4142135623730951*Mu*ZH(gt1,1)) * tmp_3985;
   std::complex<double> tmp_3989;
   std::complex<double> tmp_3990;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3990 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_3989 += tmp_3990;
   tmp_3968 += (0.7745966692414834*g1*(0.7745966692414834*g1*vd*ZH(gt1,0) -
      0.7745966692414834*g1*vu*ZH(gt1,1) - 2*(MDBS + Conj(MDBS))*ZH(gt1,2))) *
      tmp_3989;
   tmp_3967 += (2) * tmp_3968;
   std::complex<double> tmp_3991;
   std::complex<double> tmp_3992;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3992 += Conj(ZE(gt2,j1))*ZE(gt3,j1);
   }
   tmp_3991 += tmp_3992;
   tmp_3967 += (vd*(-0.6*Sqr(g1) + Sqr(g2))*ZH(gt1,0) + vu*(0.6*Sqr(g1) - Sqr(
      g2))*ZH(gt1,1) + 2*(0.7745966692414834*g1*(MDBS + Conj(MDBS))*ZH(gt1,2) + g2
      *(MDWBT + Conj(MDWBT))*ZH(gt1,3))) * tmp_3991;
   result += (0.25) * tmp_3967;

   return result;
}

std::complex<double> MRSSM_effective_couplings::CphhHpmconjHpm(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto MuD = MODELPARAMETER(MuD);
   const auto MuU = MODELPARAMETER(MuU);
   const auto vd = MODELPARAMETER(vd);
   const auto vS = MODELPARAMETER(vS);
   const auto vT = MODELPARAMETER(vT);
   const auto vu = MODELPARAMETER(vu);

   std::complex<double> result;

   result = 0.05*(7.745966692414834*g1*MDBS*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) - 20*
      vS*AbsSqr(LamSD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) - 14.142135623730951*MuD*Conj
      (LamSD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) + 7.0710678118654755*LamTD*vT*Conj(
      LamSD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) + 7.0710678118654755*LamSD*vT*Conj(
      LamTD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) + 7.745966692414834*g1*Conj(MDBS)*ZH(
      gt1,2)*ZP(gt2,0)*ZP(gt3,0) - 14.142135623730951*LamSD*Conj(MuD)*ZH(gt1,2)*ZP
      (gt2,0)*ZP(gt3,0) + 10*g2*MDWBT*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) - 10*vT*AbsSqr
      (LamTD)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) + 7.0710678118654755*LamTD*vS*Conj(
      LamSD)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) + 10*MuD*Conj(LamTD)*ZH(gt1,3)*ZP(gt2,0
      )*ZP(gt3,0) + 7.0710678118654755*LamSD*vS*Conj(LamTD)*ZH(gt1,3)*ZP(gt2,0)*ZP
      (gt3,0) + 10*g2*Conj(MDWBT)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) + 10*LamTD*Conj(
      MuD)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) - 10*LamTD*vd*Conj(LamSD)*ZH(gt1,2)*ZP(
      gt2,2)*ZP(gt3,0) + 7.0710678118654755*vd*AbsSqr(LamTD)*ZH(gt1,3)*ZP(gt2,2)*
      ZP(gt3,0) - 7.0710678118654755*vd*Sqr(g2)*ZH(gt1,3)*ZP(gt2,2)*ZP(gt3,0) - 10
      *LamSD*vd*Conj(LamTD)*ZH(gt1,2)*ZP(gt2,3)*ZP(gt3,0) - 7.0710678118654755*vd*
      AbsSqr(LamTD)*ZH(gt1,3)*ZP(gt2,3)*ZP(gt3,0) + 7.0710678118654755*vd*Sqr(g2)*
      ZH(gt1,3)*ZP(gt2,3)*ZP(gt3,0) - 7.745966692414834*g1*MDBS*ZH(gt1,2)*ZP(gt2,1
      )*ZP(gt3,1) - 20*vS*AbsSqr(LamSU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) -
      14.142135623730951*MuU*Conj(LamSU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) -
      7.0710678118654755*LamTU*vT*Conj(LamSU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) -
      7.0710678118654755*LamSU*vT*Conj(LamTU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) -
      7.745966692414834*g1*Conj(MDBS)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) -
      14.142135623730951*LamSU*Conj(MuU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) - 10*g2*
      MDWBT*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,1) - 10*vT*AbsSqr(LamTU)*ZH(gt1,3)*ZP(gt2,1
      )*ZP(gt3,1) - 7.0710678118654755*LamTU*vS*Conj(LamSU)*ZH(gt1,3)*ZP(gt2,1)*ZP
      (gt3,1) - 10*MuU*Conj(LamTU)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,1) -
      7.0710678118654755*LamSU*vS*Conj(LamTU)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,1) - 10*
      g2*Conj(MDWBT)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,1) - 10*LamTU*Conj(MuU)*ZH(gt1,3)*
      ZP(gt2,1)*ZP(gt3,1) - 10*LamTU*vu*Conj(LamSU)*ZH(gt1,2)*ZP(gt2,2)*ZP(gt3,1)
      + 7.0710678118654755*vu*AbsSqr(LamTU)*ZH(gt1,3)*ZP(gt2,2)*ZP(gt3,1) -
      7.0710678118654755*vu*Sqr(g2)*ZH(gt1,3)*ZP(gt2,2)*ZP(gt3,1) - 10*LamSU*vu*
      Conj(LamTU)*ZH(gt1,2)*ZP(gt2,3)*ZP(gt3,1) - 7.0710678118654755*vu*AbsSqr(
      LamTU)*ZH(gt1,3)*ZP(gt2,3)*ZP(gt3,1) + 7.0710678118654755*vu*Sqr(g2)*ZH(gt1,
      3)*ZP(gt2,3)*ZP(gt3,1) - 10*LamSD*vd*Conj(LamTD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,
      2) + 7.0710678118654755*vd*AbsSqr(LamTD)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,2) -
      7.0710678118654755*vd*Sqr(g2)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,2) - 10*LamSU*vu*
      Conj(LamTU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,2) + 7.0710678118654755*vu*AbsSqr(
      LamTU)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,2) - 7.0710678118654755*vu*Sqr(g2)*ZH(gt1,
      3)*ZP(gt2,1)*ZP(gt3,2) - 20*vT*Sqr(g2)*ZH(gt1,3)*ZP(gt2,2)*ZP(gt3,2) + 20*vT
      *Sqr(g2)*ZH(gt1,3)*ZP(gt2,3)*ZP(gt3,2) - 10*LamTD*vd*Conj(LamSD)*ZH(gt1,2)*
      ZP(gt2,0)*ZP(gt3,3) - 7.0710678118654755*vd*AbsSqr(LamTD)*ZH(gt1,3)*ZP(gt2,0
      )*ZP(gt3,3) + 7.0710678118654755*vd*Sqr(g2)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,3) -
      10*LamTU*vu*Conj(LamSU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,3) - 7.0710678118654755*
      vu*AbsSqr(LamTU)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,3) + 7.0710678118654755*vu*Sqr(
      g2)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,3) + 20*vT*Sqr(g2)*ZH(gt1,3)*ZP(gt2,2)*ZP(gt3
      ,3) - 20*vT*Sqr(g2)*ZH(gt1,3)*ZP(gt2,3)*ZP(gt3,3) - ZH(gt1,0)*(ZP(gt2,1)*(5*
      vu*Sqr(g2)*ZP(gt3,0) + vd*(-3*Sqr(g1) + 5*Sqr(g2))*ZP(gt3,1)) + 5*(ZP(gt2,2)
      *((2*LamTD*vS*Conj(LamSD) + 1.4142135623730951*(2*g2*MDWBT - vT*AbsSqr(LamTD
      ) + 2*LamTD*Conj(MuD) + vT*Sqr(g2)))*ZP(gt3,0) - 2*vd*(-2*AbsSqr(LamTD) +
      Sqr(g2))*ZP(gt3,2)) + ZP(gt2,3)*(((2.8284271247461903*MuD + 2*LamSD*vS +
      1.4142135623730951*LamTD*vT)*Conj(LamTD) + 1.4142135623730951*g2*(-(g2*vT) +
      2*Conj(MDWBT)))*ZP(gt3,0) + 2*vd*Sqr(g2)*ZP(gt3,3))) + ZP(gt2,0)*(vd*(3*Sqr
      (g1) + 5*Sqr(g2))*ZP(gt3,0) + 5*(vu*Sqr(g2)*ZP(gt3,1) + ((2.8284271247461903
      *MuD + 2*LamSD*vS - 1.4142135623730951*LamTD*vT)*Conj(LamTD) +
      1.4142135623730951*g2*(g2*vT + 2*Conj(MDWBT)))*ZP(gt3,2) + (2*LamTD*vS*Conj(
      LamSD) + 1.4142135623730951*(2*g2*MDWBT + vT*AbsSqr(LamTD) + 2*LamTD*Conj(
      MuD) - vT*Sqr(g2)))*ZP(gt3,3)))) - ZH(gt1,1)*(ZP(gt2,0)*((-3*vu*Sqr(g1) + 5*
      vu*Sqr(g2))*ZP(gt3,0) + 5*vd*Sqr(g2)*ZP(gt3,1)) + 5*(ZP(gt2,2)*((2*LamTU*vS*
      Conj(LamSU) + 1.4142135623730951*(2*g2*MDWBT - vT*AbsSqr(LamTU) + 2*LamTU*
      Conj(MuU) + vT*Sqr(g2)))*ZP(gt3,1) + 2*vu*Sqr(g2)*ZP(gt3,2)) + ZP(gt2,3)*(((
      2.8284271247461903*MuU + 2*LamSU*vS + 1.4142135623730951*LamTU*vT)*Conj(
      LamTU) + 1.4142135623730951*g2*(-(g2*vT) + 2*Conj(MDWBT)))*ZP(gt3,1) - 2*vu*
      (-2*AbsSqr(LamTU) + Sqr(g2))*ZP(gt3,3))) + ZP(gt2,1)*(5*vd*Sqr(g2)*ZP(gt3,0)
      + vu*(3*Sqr(g1) + 5*Sqr(g2))*ZP(gt3,1) + 5*(((2.8284271247461903*MuU + 2*
      LamSU*vS - 1.4142135623730951*LamTU*vT)*Conj(LamTU) + 1.4142135623730951*g2*
      (g2*vT + 2*Conj(MDWBT)))*ZP(gt3,2) + (2*LamTU*vS*Conj(LamSU) +
      1.4142135623730951*(2*g2*MDWBT + vT*AbsSqr(LamTU) + 2*LamTU*Conj(MuU) - vT*
      Sqr(g2)))*ZP(gt3,3)))));

   return result;
}

std::complex<double> MRSSM_effective_couplings::CpCha1hhbarCha1PL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);

   std::complex<double> result;

   result = 0.5*(-(Conj(UM1(gt3,0))*(1.4142135623730951*LamTD*Conj(UP1(gt1,1))*
      ZH(gt2,0) + 2*g2*Conj(UP1(gt1,0))*ZH(gt2,3))) - Conj(UM1(gt3,1))*(
      1.4142135623730951*g2*Conj(UP1(gt1,0))*ZH(gt2,0) + Conj(UP1(gt1,1))*(
      1.4142135623730951*LamSD*ZH(gt2,2) - LamTD*ZH(gt2,3))));

   return result;
}

std::complex<double> MRSSM_effective_couplings::CpCha2hhbarCha2PL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);

   std::complex<double> result;

   result = 0.5*(g2*Conj(UM2(gt1,0))*(-1.4142135623730951*Conj(UP2(gt3,1))*ZH(
      gt2,1) + 2*Conj(UP2(gt3,0))*ZH(gt2,3)) + Conj(UM2(gt1,1))*(
      1.4142135623730951*LamTU*Conj(UP2(gt3,0))*ZH(gt2,1) + Conj(UP2(gt3,1))*(
      1.4142135623730951*LamSU*ZH(gt2,2) + LamTU*ZH(gt2,3))));

   return result;
}

std::complex<double> MRSSM_effective_couplings::CpFehhbarFePL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Ye = MODELPARAMETER(Ye);

   std::complex<double> result;

   std::complex<double> tmp_3993;
   std::complex<double> tmp_3994;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3995;
      std::complex<double> tmp_3996;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3996 += Conj(ZER(gt3,j1))*Ye(j1,j2);
      }
      tmp_3995 += tmp_3996;
      tmp_3994 += (Conj(ZEL(gt1,j2))) * tmp_3995;
   }
   tmp_3993 += tmp_3994;
   result += (-0.7071067811865475*ZH(gt2,0)) * tmp_3993;

   return result;
}

std::complex<double> MRSSM_effective_couplings::CpFdhhbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_3997;
   std::complex<double> tmp_3998;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3999;
      std::complex<double> tmp_4000;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4000 += Conj(ZDR(gt3,j1))*Yd(j1,j2);
      }
      tmp_3999 += tmp_4000;
      tmp_3998 += (Conj(ZDL(gt1,j2))) * tmp_3999;
   }
   tmp_3997 += tmp_3998;
   result += (-0.7071067811865475*ZH(gt2,0)) * tmp_3997;

   return result;
}

std::complex<double> MRSSM_effective_couplings::CpFuhhbarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_4001;
   std::complex<double> tmp_4002;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4003;
      std::complex<double> tmp_4004;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4004 += Conj(ZUR(gt3,j1))*Yu(j1,j2);
      }
      tmp_4003 += tmp_4004;
      tmp_4002 += (Conj(ZUL(gt1,j2))) * tmp_4003;
   }
   tmp_4001 += tmp_4002;
   result += (-0.7071067811865475*ZH(gt2,1)) * tmp_4001;

   return result;
}

std::complex<double> MRSSM_effective_couplings::CphhVWmconjVWm(unsigned gt1) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vT = MODELPARAMETER(vT);
   const auto vu = MODELPARAMETER(vu);

   std::complex<double> result;

   result = 0.5*Sqr(g2)*(vd*ZH(gt1,0) + vu*ZH(gt1,1) + 4*vT*ZH(gt1,3));

   return result;
}

std::complex<double> MRSSM_effective_couplings::CpAhSRdpconjSRdp(unsigned gt1) const
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto MuD = MODELPARAMETER(MuD);
   const auto vS = MODELPARAMETER(vS);
   const auto vT = MODELPARAMETER(vT);

   std::complex<double> result;

   result = std::complex<double>(0,-0.25)*((1.5491933384829668*g1*MDBS +
      1.4142135623730951*(-2*MuD + LamTD*vT)*Conj(LamSD) - 1.4142135623730951*
      LamSD*vT*Conj(LamTD) - 1.5491933384829668*g1*Conj(MDBS) + 2.8284271247461903
      *LamSD*Conj(MuD))*ZA(gt1,2) + (2*g2*MDWBT - 1.4142135623730951*LamTD*vS*Conj
      (LamSD) + (2*MuD + 1.4142135623730951*LamSD*vS)*Conj(LamTD) - 2*g2*Conj(
      MDWBT) - 2*LamTD*Conj(MuD))*ZA(gt1,3));

   return result;
}

std::complex<double> MRSSM_effective_couplings::CpAhSRumconjSRum(unsigned gt1) const
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto MuU = MODELPARAMETER(MuU);
   const auto vS = MODELPARAMETER(vS);
   const auto vT = MODELPARAMETER(vT);

   std::complex<double> result;

   result = std::complex<double>(0,0.05)*((7.745966692414834*g1*MDBS +
      7.0710678118654755*(2*MuU + LamTU*vT)*Conj(LamSU) - 7.0710678118654755*LamSU
      *vT*Conj(LamTU) - 7.745966692414834*g1*Conj(MDBS) - 14.142135623730951*LamSU
      *Conj(MuU))*ZA(gt1,2) + 5*(2*g2*MDWBT - 1.4142135623730951*LamTU*vS*Conj(
      LamSU) + 2*MuU*Conj(LamTU) + 1.4142135623730951*LamSU*vS*Conj(LamTU) - 2*g2*
      Conj(MDWBT) - 2*LamTU*Conj(MuU))*ZA(gt1,3));

   return result;
}

std::complex<double> MRSSM_effective_couplings::CpAhSdconjSd(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_4005;
   std::complex<double> tmp_4006;
   std::complex<double> tmp_4007;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4008;
      std::complex<double> tmp_4009;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4009 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_4008 += tmp_4009;
      tmp_4007 += (Conj(ZD(gt2,j2))) * tmp_4008;
   }
   tmp_4006 += tmp_4007;
   tmp_4005 += (4.242640687119286*Conj(Mu)*ZA(gt1,1)) * tmp_4006;
   std::complex<double> tmp_4010;
   std::complex<double> tmp_4011;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4012;
      std::complex<double> tmp_4013;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4013 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_4012 += tmp_4013;
      tmp_4011 += (ZD(gt3,j2)) * tmp_4012;
   }
   tmp_4010 += tmp_4011;
   tmp_4005 += (-4.242640687119286*Mu*ZA(gt1,1)) * tmp_4010;
   std::complex<double> tmp_4014;
   std::complex<double> tmp_4015;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4015 += Conj(ZD(gt2,j1))*ZD(gt3,j1);
   }
   tmp_4014 += tmp_4015;
   tmp_4005 += (0.7745966692414834*g1*MDBS*ZA(gt1,2)) * tmp_4014;
   std::complex<double> tmp_4016;
   std::complex<double> tmp_4017;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4017 += Conj(ZD(gt2,j1))*ZD(gt3,j1);
   }
   tmp_4016 += tmp_4017;
   tmp_4005 += (-0.7745966692414834*g1*Conj(MDBS)*ZA(gt1,2)) * tmp_4016;
   std::complex<double> tmp_4018;
   std::complex<double> tmp_4019;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4019 += Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1);
   }
   tmp_4018 += tmp_4019;
   tmp_4005 += (1.5491933384829668*g1*MDBS*ZA(gt1,2)) * tmp_4018;
   std::complex<double> tmp_4020;
   std::complex<double> tmp_4021;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4021 += Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1);
   }
   tmp_4020 += tmp_4021;
   tmp_4005 += (-1.5491933384829668*g1*Conj(MDBS)*ZA(gt1,2)) * tmp_4020;
   std::complex<double> tmp_4022;
   std::complex<double> tmp_4023;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4023 += Conj(ZD(gt2,j1))*ZD(gt3,j1);
   }
   tmp_4022 += tmp_4023;
   tmp_4005 += (-3*g2*MDWBT*ZA(gt1,3)) * tmp_4022;
   std::complex<double> tmp_4024;
   std::complex<double> tmp_4025;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4025 += Conj(ZD(gt2,j1))*ZD(gt3,j1);
   }
   tmp_4024 += tmp_4025;
   tmp_4005 += (3*g2*Conj(MDWBT)*ZA(gt1,3)) * tmp_4024;
   result += (std::complex<double>(0,-0.16666666666666666)) * tmp_4005;

   return result;
}

std::complex<double> MRSSM_effective_couplings::CpAhSuconjSu(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_4026;
   std::complex<double> tmp_4027;
   std::complex<double> tmp_4028;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4029;
      std::complex<double> tmp_4030;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4030 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_4029 += tmp_4030;
      tmp_4028 += (Conj(ZU(gt2,j2))) * tmp_4029;
   }
   tmp_4027 += tmp_4028;
   tmp_4026 += (4.242640687119286*Conj(Mu)*ZA(gt1,0)) * tmp_4027;
   std::complex<double> tmp_4031;
   std::complex<double> tmp_4032;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4033;
      std::complex<double> tmp_4034;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4034 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_4033 += tmp_4034;
      tmp_4032 += (ZU(gt3,j2)) * tmp_4033;
   }
   tmp_4031 += tmp_4032;
   tmp_4026 += (-4.242640687119286*Mu*ZA(gt1,0)) * tmp_4031;
   std::complex<double> tmp_4035;
   std::complex<double> tmp_4036;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4036 += Conj(ZU(gt2,j1))*ZU(gt3,j1);
   }
   tmp_4035 += tmp_4036;
   tmp_4026 += (0.7745966692414834*g1*MDBS*ZA(gt1,2)) * tmp_4035;
   std::complex<double> tmp_4037;
   std::complex<double> tmp_4038;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4038 += Conj(ZU(gt2,j1))*ZU(gt3,j1);
   }
   tmp_4037 += tmp_4038;
   tmp_4026 += (-0.7745966692414834*g1*Conj(MDBS)*ZA(gt1,2)) * tmp_4037;
   std::complex<double> tmp_4039;
   std::complex<double> tmp_4040;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4040 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_4039 += tmp_4040;
   tmp_4026 += (-3.0983866769659336*g1*MDBS*ZA(gt1,2)) * tmp_4039;
   std::complex<double> tmp_4041;
   std::complex<double> tmp_4042;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4042 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_4041 += tmp_4042;
   tmp_4026 += (3.0983866769659336*g1*Conj(MDBS)*ZA(gt1,2)) * tmp_4041;
   std::complex<double> tmp_4043;
   std::complex<double> tmp_4044;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4044 += Conj(ZU(gt2,j1))*ZU(gt3,j1);
   }
   tmp_4043 += tmp_4044;
   tmp_4026 += (3*g2*MDWBT*ZA(gt1,3)) * tmp_4043;
   std::complex<double> tmp_4045;
   std::complex<double> tmp_4046;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4046 += Conj(ZU(gt2,j1))*ZU(gt3,j1);
   }
   tmp_4045 += tmp_4046;
   tmp_4026 += (-3*g2*Conj(MDWBT)*ZA(gt1,3)) * tmp_4045;
   result += (std::complex<double>(0,-0.16666666666666666)) * tmp_4026;

   return result;
}

std::complex<double> MRSSM_effective_couplings::CpAhSeconjSe(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_4047;
   std::complex<double> tmp_4048;
   std::complex<double> tmp_4049;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4050;
      std::complex<double> tmp_4051;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4051 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_4050 += tmp_4051;
      tmp_4049 += (Conj(ZE(gt2,j2))) * tmp_4050;
   }
   tmp_4048 += tmp_4049;
   tmp_4047 += (1.4142135623730951*Conj(Mu)*ZA(gt1,1)) * tmp_4048;
   std::complex<double> tmp_4052;
   std::complex<double> tmp_4053;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4054;
      std::complex<double> tmp_4055;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4055 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_4054 += tmp_4055;
      tmp_4053 += (ZE(gt3,j2)) * tmp_4054;
   }
   tmp_4052 += tmp_4053;
   tmp_4047 += (-1.4142135623730951*Mu*ZA(gt1,1)) * tmp_4052;
   std::complex<double> tmp_4056;
   std::complex<double> tmp_4057;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4057 += Conj(ZE(gt2,j1))*ZE(gt3,j1);
   }
   tmp_4056 += tmp_4057;
   tmp_4047 += (-0.7745966692414834*g1*MDBS*ZA(gt1,2)) * tmp_4056;
   std::complex<double> tmp_4058;
   std::complex<double> tmp_4059;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4059 += Conj(ZE(gt2,j1))*ZE(gt3,j1);
   }
   tmp_4058 += tmp_4059;
   tmp_4047 += (0.7745966692414834*g1*Conj(MDBS)*ZA(gt1,2)) * tmp_4058;
   std::complex<double> tmp_4060;
   std::complex<double> tmp_4061;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4061 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_4060 += tmp_4061;
   tmp_4047 += (1.5491933384829668*g1*MDBS*ZA(gt1,2)) * tmp_4060;
   std::complex<double> tmp_4062;
   std::complex<double> tmp_4063;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4063 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_4062 += tmp_4063;
   tmp_4047 += (-1.5491933384829668*g1*Conj(MDBS)*ZA(gt1,2)) * tmp_4062;
   std::complex<double> tmp_4064;
   std::complex<double> tmp_4065;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4065 += Conj(ZE(gt2,j1))*ZE(gt3,j1);
   }
   tmp_4064 += tmp_4065;
   tmp_4047 += (-(g2*MDWBT*ZA(gt1,3))) * tmp_4064;
   std::complex<double> tmp_4066;
   std::complex<double> tmp_4067;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4067 += Conj(ZE(gt2,j1))*ZE(gt3,j1);
   }
   tmp_4066 += tmp_4067;
   tmp_4047 += (g2*Conj(MDWBT)*ZA(gt1,3)) * tmp_4066;
   result += (std::complex<double>(0,-0.5)) * tmp_4047;

   return result;
}

std::complex<double> MRSSM_effective_couplings::CpAhHpmconjHpm(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto MuD = MODELPARAMETER(MuD);
   const auto MuU = MODELPARAMETER(MuU);
   const auto vd = MODELPARAMETER(vd);
   const auto vS = MODELPARAMETER(vS);
   const auto vT = MODELPARAMETER(vT);
   const auto vu = MODELPARAMETER(vu);

   std::complex<double> result;

   result = std::complex<double>(0,0.05)*(ZA(gt1,2)*(10*LamTD*vd*Conj(LamSD)*ZP
      (gt2,2)*ZP(gt3,0) - 10*LamSD*vd*Conj(LamTD)*ZP(gt2,3)*ZP(gt3,0) -
      7.745966692414834*g1*MDBS*ZP(gt2,1)*ZP(gt3,1) + 14.142135623730951*MuU*Conj(
      LamSU)*ZP(gt2,1)*ZP(gt3,1) + 7.0710678118654755*LamTU*vT*Conj(LamSU)*ZP(gt2,
      1)*ZP(gt3,1) - 7.0710678118654755*LamSU*vT*Conj(LamTU)*ZP(gt2,1)*ZP(gt3,1) +
      7.745966692414834*g1*Conj(MDBS)*ZP(gt2,1)*ZP(gt3,1) - 14.142135623730951*
      LamSU*Conj(MuU)*ZP(gt2,1)*ZP(gt3,1) + 10*LamTU*vu*Conj(LamSU)*ZP(gt2,2)*ZP(
      gt3,1) - 10*LamSU*vu*Conj(LamTU)*ZP(gt2,3)*ZP(gt3,1) - 10*LamSU*vu*Conj(
      LamTU)*ZP(gt2,1)*ZP(gt3,2) + 10*LamTU*vu*Conj(LamSU)*ZP(gt2,1)*ZP(gt3,3) +
      ZP(gt2,0)*((7.745966692414834*g1*MDBS + 7.0710678118654755*(2*MuD - LamTD*vT
      )*Conj(LamSD) + 7.0710678118654755*LamSD*vT*Conj(LamTD) - 7.745966692414834*
      g1*Conj(MDBS) - 14.142135623730951*LamSD*Conj(MuD))*ZP(gt3,0) + 10*vd*(-(
      LamSD*Conj(LamTD)*ZP(gt3,2)) + LamTD*Conj(LamSD)*ZP(gt3,3)))) - 5*(ZA(gt1,0)
      *(vu*Sqr(g2)*ZP(gt2,1)*ZP(gt3,0) + (2*LamTD*vS*Conj(LamSD) +
      1.4142135623730951*(2*g2*MDWBT - vT*AbsSqr(LamTD) + 2*LamTD*Conj(MuD) + vT*
      Sqr(g2)))*ZP(gt2,2)*ZP(gt3,0) + 1.4142135623730951*vT*AbsSqr(LamTD)*ZP(gt2,3
      )*ZP(gt3,0) + 2.8284271247461903*MuD*Conj(LamTD)*ZP(gt2,3)*ZP(gt3,0) + 2*
      LamSD*vS*Conj(LamTD)*ZP(gt2,3)*ZP(gt3,0) + 2.8284271247461903*g2*Conj(MDWBT)
      *ZP(gt2,3)*ZP(gt3,0) - 1.4142135623730951*vT*Sqr(g2)*ZP(gt2,3)*ZP(gt3,0) -
      vu*Sqr(g2)*ZP(gt2,0)*ZP(gt3,1) + 1.4142135623730951*vT*AbsSqr(LamTD)*ZP(gt2,
      0)*ZP(gt3,2) - 2.8284271247461903*MuD*Conj(LamTD)*ZP(gt2,0)*ZP(gt3,2) - 2*
      LamSD*vS*Conj(LamTD)*ZP(gt2,0)*ZP(gt3,2) - 2.8284271247461903*g2*Conj(MDWBT)
      *ZP(gt2,0)*ZP(gt3,2) - 1.4142135623730951*vT*Sqr(g2)*ZP(gt2,0)*ZP(gt3,2) -
      2.8284271247461903*g2*MDWBT*ZP(gt2,0)*ZP(gt3,3) - 1.4142135623730951*vT*
      AbsSqr(LamTD)*ZP(gt2,0)*ZP(gt3,3) - 2*LamTD*vS*Conj(LamSD)*ZP(gt2,0)*ZP(gt3,
      3) - 2.8284271247461903*LamTD*Conj(MuD)*ZP(gt2,0)*ZP(gt3,3) +
      1.4142135623730951*vT*Sqr(g2)*ZP(gt2,0)*ZP(gt3,3)) + ZA(gt1,1)*(-((vd*Sqr(g2
      )*ZP(gt2,0) + (2*LamTU*vS*Conj(LamSU) + 1.4142135623730951*(2*g2*MDWBT - vT*
      AbsSqr(LamTU) + 2*LamTU*Conj(MuU) + vT*Sqr(g2)))*ZP(gt2,2) + ((
      2.8284271247461903*MuU + 2*LamSU*vS + 1.4142135623730951*LamTU*vT)*Conj(
      LamTU) + 1.4142135623730951*g2*(-(g2*vT) + 2*Conj(MDWBT)))*ZP(gt2,3))*ZP(gt3
      ,1)) + ZP(gt2,1)*(vd*Sqr(g2)*ZP(gt3,0) + ((2.8284271247461903*MuU + 2*LamSU*
      vS - 1.4142135623730951*LamTU*vT)*Conj(LamTU) + 1.4142135623730951*g2*(g2*vT
      + 2*Conj(MDWBT)))*ZP(gt3,2) + (2*LamTU*vS*Conj(LamSU) + 1.4142135623730951*
      (2*g2*MDWBT + vT*AbsSqr(LamTU) + 2*LamTU*Conj(MuU) - vT*Sqr(g2)))*ZP(gt3,3))
      ) + ZA(gt1,3)*(1.4142135623730951*vd*AbsSqr(LamTD)*ZP(gt2,3)*ZP(gt3,0) -
      1.4142135623730951*vd*Sqr(g2)*ZP(gt2,3)*ZP(gt3,0) + 2*g2*MDWBT*ZP(gt2,1)*ZP(
      gt3,1) + 1.4142135623730951*LamTU*vS*Conj(LamSU)*ZP(gt2,1)*ZP(gt3,1) - 2*MuU
      *Conj(LamTU)*ZP(gt2,1)*ZP(gt3,1) - 1.4142135623730951*LamSU*vS*Conj(LamTU)*
      ZP(gt2,1)*ZP(gt3,1) - 2*g2*Conj(MDWBT)*ZP(gt2,1)*ZP(gt3,1) + 2*LamTU*Conj(
      MuU)*ZP(gt2,1)*ZP(gt3,1) + 1.4142135623730951*vu*AbsSqr(LamTU)*ZP(gt2,3)*ZP(
      gt3,1) - 1.4142135623730951*vu*Sqr(g2)*ZP(gt2,3)*ZP(gt3,1) -
      1.4142135623730951*vu*AbsSqr(LamTU)*ZP(gt2,1)*ZP(gt3,2) + 1.4142135623730951
      *vu*Sqr(g2)*ZP(gt2,1)*ZP(gt3,2) - 4*vT*Sqr(g2)*ZP(gt2,3)*ZP(gt3,2) -
      1.4142135623730951*vu*AbsSqr(LamTU)*ZP(gt2,1)*ZP(gt3,3) + 1.4142135623730951
      *vu*Sqr(g2)*ZP(gt2,1)*ZP(gt3,3) + ZP(gt2,2)*(1.4142135623730951*vd*(AbsSqr(
      LamTD) - Sqr(g2))*ZP(gt3,0) + 1.4142135623730951*vu*(AbsSqr(LamTU) - Sqr(g2)
      )*ZP(gt3,1) + 4*vT*Sqr(g2)*ZP(gt3,3)) + ZP(gt2,0)*((-2*g2*MDWBT -
      1.4142135623730951*LamTD*vS*Conj(LamSD) + 2*MuD*Conj(LamTD) +
      1.4142135623730951*LamSD*vS*Conj(LamTD) + 2*g2*Conj(MDWBT) - 2*LamTD*Conj(
      MuD))*ZP(gt3,0) + 1.4142135623730951*vd*(-AbsSqr(LamTD) + Sqr(g2))*(ZP(gt3,2
      ) + ZP(gt3,3))))));

   return result;
}

std::complex<double> MRSSM_effective_couplings::CpAhCha1barCha1PL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);

   std::complex<double> result;

   result = std::complex<double>(0,0.5)*(Conj(UM1(gt3,0))*(-1.4142135623730951*
      LamTD*Conj(UP1(gt2,1))*ZA(gt1,0) + 2*g2*Conj(UP1(gt2,0))*ZA(gt1,3)) + Conj(
      UM1(gt3,1))*(1.4142135623730951*g2*Conj(UP1(gt2,0))*ZA(gt1,0) + Conj(UP1(gt2
      ,1))*(-1.4142135623730951*LamSD*ZA(gt1,2) + LamTD*ZA(gt1,3))));

   return result;
}

std::complex<double> MRSSM_effective_couplings::CpAhCha2barCha2PL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);

   std::complex<double> result;

   result = std::complex<double>(0,0.5)*(g2*Conj(UM2(gt2,0))*(
      1.4142135623730951*Conj(UP2(gt3,1))*ZA(gt1,1) - 2*Conj(UP2(gt3,0))*ZA(gt1,3)
      ) + Conj(UM2(gt2,1))*(1.4142135623730951*LamTU*Conj(UP2(gt3,0))*ZA(gt1,1) +
      Conj(UP2(gt3,1))*(1.4142135623730951*LamSU*ZA(gt1,2) + LamTU*ZA(gt1,3))));

   return result;
}

std::complex<double> MRSSM_effective_couplings::CpAhFebarFePL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Ye = MODELPARAMETER(Ye);

   std::complex<double> result;

   std::complex<double> tmp_4068;
   std::complex<double> tmp_4069;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4070;
      std::complex<double> tmp_4071;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4071 += Conj(ZER(gt3,j1))*Ye(j1,j2);
      }
      tmp_4070 += tmp_4071;
      tmp_4069 += (Conj(ZEL(gt2,j2))) * tmp_4070;
   }
   tmp_4068 += tmp_4069;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,0)) *
      tmp_4068;

   return result;
}

std::complex<double> MRSSM_effective_couplings::CpAhFdbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_4072;
   std::complex<double> tmp_4073;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4074;
      std::complex<double> tmp_4075;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4075 += Conj(ZDR(gt3,j1))*Yd(j1,j2);
      }
      tmp_4074 += tmp_4075;
      tmp_4073 += (Conj(ZDL(gt2,j2))) * tmp_4074;
   }
   tmp_4072 += tmp_4073;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,0)) *
      tmp_4072;

   return result;
}

std::complex<double> MRSSM_effective_couplings::CpAhFubarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_4076;
   std::complex<double> tmp_4077;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4078;
      std::complex<double> tmp_4079;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4079 += Conj(ZUR(gt3,j1))*Yu(j1,j2);
      }
      tmp_4078 += tmp_4079;
      tmp_4077 += (Conj(ZUL(gt2,j2))) * tmp_4078;
   }
   tmp_4076 += tmp_4077;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,1)) *
      tmp_4076;

   return result;
}

void MRSSM_effective_couplings::calculate_eff_CphhVPVP(unsigned gO1)
{
   const auto vd = MODELPARAMETER(vd);
   const auto vT = MODELPARAMETER(vT);
   const auto vu = MODELPARAMETER(vu);
   const auto MSRdp = MODELPARAMETER(MSRdp);
   const auto MSRum = MODELPARAMETER(MSRum);
   const auto MSd = MODELPARAMETER(MSd);
   const auto MSu = MODELPARAMETER(MSu);
   const auto MSe = MODELPARAMETER(MSe);
   const auto MHpm = MODELPARAMETER(MHpm);
   const auto MCha1 = MODELPARAMETER(MCha1);
   const auto MCha2 = MODELPARAMETER(MCha2);
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
   result += 0.5 * CphhSRdpconjSRdp(gO1) * vev * AS0(decay_scale / Sqr(MSRdp))
      / Sqr(MSRdp);
   result += 0.5 * CphhSRumconjSRum(gO1) * vev * AS0(decay_scale / Sqr(MSRum))
      / Sqr(MSRum);
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
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      result += CpCha1hhbarCha1PL(gI1, gO1, gI1) * vev * AS12(decay_scale /
         Sqr(MCha1(gI1))) / MCha1(gI1);
   }
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      result += CpCha2hhbarCha2PL(gI1, gO1, gI1) * vev * AS12(decay_scale /
         Sqr(MCha2(gI1))) / MCha2(gI1);
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

void MRSSM_effective_couplings::calculate_eff_CphhVGVG(unsigned gO1)
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

void MRSSM_effective_couplings::calculate_eff_CpAhVPVP(unsigned gO1)
{
   const auto vd = MODELPARAMETER(vd);
   const auto vT = MODELPARAMETER(vT);
   const auto vu = MODELPARAMETER(vu);
   const auto MCha1 = MODELPARAMETER(MCha1);
   const auto MCha2 = MODELPARAMETER(MCha2);
   const auto MFe = MODELPARAMETER(MFe);
   const auto MFd = MODELPARAMETER(MFd);
   const auto MFu = MODELPARAMETER(MFu);
   const auto decay_mass = PHYSICAL(MAh)(gO1);
   const auto decay_scale = 0.25 * Sqr(decay_mass);
   const auto saved_ZA = ZA;
   ZA = PHYSICAL(ZA);

   const auto vev = Sqrt(Sqr(vd) + 4*Sqr(vT) + Sqr(vu));

   std::complex<double> result = 0;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      result += CpAhCha1barCha1PL(gO1, gI1, gI1) * vev * AP12(decay_scale /
         Sqr(MCha1(gI1))) / MCha1(gI1);
   }
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      result += CpAhCha2barCha2PL(gO1, gI1, gI1) * vev * AP12(decay_scale /
         Sqr(MCha2(gI1))) / MCha2(gI1);
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

void MRSSM_effective_couplings::calculate_eff_CpAhVGVG(unsigned gO1)
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
