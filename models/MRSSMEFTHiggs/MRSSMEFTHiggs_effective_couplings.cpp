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

// File generated at Sun 4 Aug 2019 18:15:17

#include "MRSSMEFTHiggs_effective_couplings.hpp"

#include "effective_couplings.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

using namespace effective_couplings;

#define INPUTPARAMETER(parameter) model.get_input().parameter
#define MODELPARAMETER(parameter) model.get_##parameter()
#define DERIVEDPARAMETER(parameter) model.##parameter()
#define PHASE(parameter) model.get_##parameter()
#define PHYSICAL(parameter) model.get_physical().parameter

MRSSMEFTHiggs_effective_couplings::MRSSMEFTHiggs_effective_couplings(
   const MRSSMEFTHiggs_mass_eigenstates& model_,
   const softsusy::QedQcd& qedqcd_,
   const Physical_input& input_)
   : model(model_), qedqcd(qedqcd_), physical_input(input_)
   , rg_improve(true), include_qcd_corrections(true)
   , ZD(MODELPARAMETER(ZD)), ZV(MODELPARAMETER(ZV)), ZU(MODELPARAMETER(ZU)), ZE(
      MODELPARAMETER(ZE)), ZH(MODELPARAMETER(ZH)), ZA(MODELPARAMETER(ZA)), ZHR(
      MODELPARAMETER(ZHR)), ZP(MODELPARAMETER(ZP)), ZN1(MODELPARAMETER(ZN1)), ZN2(
      MODELPARAMETER(ZN2)), UM1(MODELPARAMETER(UM1)), UP1(MODELPARAMETER(UP1)),
      UM2(MODELPARAMETER(UM2)), UP2(MODELPARAMETER(UP2)), ZEL(MODELPARAMETER(ZEL))
      , ZER(MODELPARAMETER(ZER)), ZDL(MODELPARAMETER(ZDL)), ZDR(MODELPARAMETER(ZDR
      )), ZUL(MODELPARAMETER(ZUL)), ZUR(MODELPARAMETER(ZUR)), ZZ(MODELPARAMETER(ZZ
      ))
   , eff_CphhVPVP(Eigen::Array<std::complex<double>,4,1>::Zero()), eff_CphhVGVG(
      Eigen::Array<std::complex<double>,4,1>::Zero()), eff_CpAhVPVP(Eigen::Array<
      std::complex<double>,4,1>::Zero()), eff_CpAhVGVG(Eigen::Array<std::complex<
      double>,4,1>::Zero())
{
}

void MRSSMEFTHiggs_effective_couplings::calculate_effective_couplings()
{
   const standard_model::Standard_model sm(initialise_SM());

   const double scale = model.get_scale();
   const Eigen::ArrayXd saved_parameters(model.get());

   const double saved_mt = PHYSICAL(MFu(2));
   PHYSICAL(MFu(2)) = qedqcd.displayPoleMt();

   const auto Mhh = PHYSICAL(Mhh);
   for (int gO1 = 0; gO1 < 4; ++gO1) {
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
   for (int gO1 = 1; gO1 < 4; ++gO1) {
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

void MRSSMEFTHiggs_effective_couplings::set_model(const MRSSMEFTHiggs_mass_eigenstates& model_)
{
   model = model_;
   copy_mixing_matrices_from_model();
}

void MRSSMEFTHiggs_effective_couplings::copy_mixing_matrices_from_model()
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

standard_model::Standard_model MRSSMEFTHiggs_effective_couplings::initialise_SM() const
{
   standard_model::Standard_model sm;

   sm.set_loops(2);
   sm.set_thresholds(2);
   sm.set_physical_input(physical_input);

   sm.initialise_from_input(qedqcd);

   return sm;
}

void MRSSMEFTHiggs_effective_couplings::run_SM_strong_coupling_to(standard_model::Standard_model sm, double m)
{
   sm.run_to(m);

   model.set_g3(sm.get_g3());

}

std::complex<double> MRSSMEFTHiggs_effective_couplings::scalar_scalar_qcd_factor(double m_decay, double m_loop) const
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

std::complex<double> MRSSMEFTHiggs_effective_couplings::scalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) * scalar_diphoton_fermion_loop(
         m_decay, m_loop);

   }

   return result;
}

std::complex<double> MRSSMEFTHiggs_effective_couplings::pseudoscalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) *
         pseudoscalar_diphoton_fermion_loop(m_decay, m_loop);

   }

   return result;
}

double MRSSMEFTHiggs_effective_couplings::number_of_active_flavours(double m) const
{
   if (m < qedqcd.displayMbMb()) {
      return 4.0;
   } else if (m < qedqcd.displayPoleMt()) {
      return 5.0;
   } else {
      return 6.0;
   }
}

double MRSSMEFTHiggs_effective_couplings::scalar_scaling_factor(double m) const
{
   const double Nf = number_of_active_flavours(m);
   const double mtpole = qedqcd.displayPoleMt();
   const double l = Log(Sqr(m) / Sqr(mtpole));

   const auto g3 = MODELPARAMETER(g3);

   const double nlo_qcd = 0.025330295910584444*(23.75 - 1.1666666666666667*Nf)*Sqr
      (g3);
   const double nnlo_qcd = 0.000641623890917771*Quad(g3)*(370.1956513893174 +
      2.375*l + (-47.18640261449638 + 0.6666666666666666*l)*Nf +
      0.9017702481178881*Sqr(Nf));
   const double nnnlo_qcd = 0.000016252523020247696*Power6(g3)*(467.683620788 +
      122.440972222*l + 10.9409722222*Sqr(l));

   return Sqrt(1.0 + nlo_qcd + nnlo_qcd + nnnlo_qcd);
}

double MRSSMEFTHiggs_effective_couplings::pseudoscalar_scaling_factor(double m) const
{
   const double Nf = number_of_active_flavours(m);
   const double mtpole = qedqcd.displayPoleMt();
   const double l = Log(Sqr(m) / Sqr(mtpole));

   const auto g3 = MODELPARAMETER(g3);

   const double nlo_qcd = 0.025330295910584444*(24.25 - 1.1666666666666667*Nf)*Sqr
      (g3);
   const double nnlo_qcd = 0.000641623890917771*(171.54400563089382 + 5*l)*Quad(g3
      );
   const double nnnlo_qcd = 0;

   return Sqrt(1.0 + nlo_qcd + nnlo_qcd + nnnlo_qcd);
}

double MRSSMEFTHiggs_effective_couplings::get_hhVPVP_partial_width(int gO1) const
{
   const double mass = PHYSICAL(Mhh)(gO1);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CphhVPVP(gO1));
}

double MRSSMEFTHiggs_effective_couplings::get_hhVGVG_partial_width(int gO1) const
{
   const double mass = PHYSICAL(Mhh)(gO1);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CphhVGVG(gO1));
}

double MRSSMEFTHiggs_effective_couplings::get_AhVPVP_partial_width(int gO1) const
{
   const double mass = PHYSICAL(MAh)(gO1);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CpAhVPVP(gO1));
}

double MRSSMEFTHiggs_effective_couplings::get_AhVGVG_partial_width(int gO1) const
{
   const double mass = PHYSICAL(MAh)(gO1);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CpAhVGVG(gO1));
}

std::complex<double> MRSSMEFTHiggs_effective_couplings::CpSRdpAhconjSRdp(int gI2) const
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

   const std::complex<double> result = std::complex<double>(0,-0.25)*((
      1.5491933384829668*g1*MDBS + 1.4142135623730951*(-2*MuD + LamTD*vT)*Conj(
      LamSD) - 1.4142135623730951*LamSD*vT*Conj(LamTD) - 1.5491933384829668*g1*
      Conj(MDBS) + 2.8284271247461903*LamSD*Conj(MuD))*ZA(gI2,2) + (2*g2*MDWBT -
      1.4142135623730951*LamTD*vS*Conj(LamSD) + (2*MuD + 1.4142135623730951*LamSD*
      vS)*Conj(LamTD) - 2*g2*Conj(MDWBT) - 2*LamTD*Conj(MuD))*ZA(gI2,3));

   return result;
}

std::complex<double> MRSSMEFTHiggs_effective_couplings::CpSRdphhconjSRdp(int gI2) const
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

   const std::complex<double> result = 0.25*(vd*(-4*AbsSqr(LamTD) + 0.6*Sqr(g1) -
      Sqr(g2))*ZH(gI2,0) + vu*(-0.6*Sqr(g1) + Sqr(g2))*ZH(gI2,1) -
      1.5491933384829668*g1*MDBS*ZH(gI2,2) - 4*vS*AbsSqr(LamSD)*ZH(gI2,2) -
      2.8284271247461903*MuD*Conj(LamSD)*ZH(gI2,2) + 1.4142135623730951*LamTD*vT*
      Conj(LamSD)*ZH(gI2,2) + 1.4142135623730951*LamSD*vT*Conj(LamTD)*ZH(gI2,2) -
      1.5491933384829668*g1*Conj(MDBS)*ZH(gI2,2) - 2.8284271247461903*LamSD*Conj(
      MuD)*ZH(gI2,2) - 2*g2*MDWBT*ZH(gI2,3) - 2*vT*AbsSqr(LamTD)*ZH(gI2,3) +
      1.4142135623730951*LamTD*vS*Conj(LamSD)*ZH(gI2,3) + 2*MuD*Conj(LamTD)*ZH(gI2
      ,3) + 1.4142135623730951*LamSD*vS*Conj(LamTD)*ZH(gI2,3) - 2*g2*Conj(MDWBT)*
      ZH(gI2,3) + 2*LamTD*Conj(MuD)*ZH(gI2,3));

   return result;
}

std::complex<double> MRSSMEFTHiggs_effective_couplings::CpSRumAhconjSRum(int gI2) const
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

   const std::complex<double> result = std::complex<double>(0,0.05)*((
      7.745966692414834*g1*MDBS + 7.0710678118654755*(2*MuU + LamTU*vT)*Conj(LamSU
      ) - 7.0710678118654755*LamSU*vT*Conj(LamTU) - 7.745966692414834*g1*Conj(MDBS
      ) - 14.142135623730951*LamSU*Conj(MuU))*ZA(gI2,2) + 5*(2*g2*MDWBT -
      1.4142135623730951*LamTU*vS*Conj(LamSU) + 2*MuU*Conj(LamTU) +
      1.4142135623730951*LamSU*vS*Conj(LamTU) - 2*g2*Conj(MDWBT) - 2*LamTU*Conj(
      MuU))*ZA(gI2,3));

   return result;
}

std::complex<double> MRSSMEFTHiggs_effective_couplings::CpSRumhhconjSRum(int gI2) const
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

   const std::complex<double> result = 0.25*(vd*(-0.6*Sqr(g1) + Sqr(g2))*ZH(gI2,0)
      + vu*(-4*AbsSqr(LamTU) + 0.6*Sqr(g1) - Sqr(g2))*ZH(gI2,1) +
      1.5491933384829668*g1*MDBS*ZH(gI2,2) - 4*vS*AbsSqr(LamSU)*ZH(gI2,2) -
      2.8284271247461903*MuU*Conj(LamSU)*ZH(gI2,2) - 1.4142135623730951*LamTU*vT*
      Conj(LamSU)*ZH(gI2,2) - 1.4142135623730951*LamSU*vT*Conj(LamTU)*ZH(gI2,2) +
      1.5491933384829668*g1*Conj(MDBS)*ZH(gI2,2) - 2.8284271247461903*LamSU*Conj(
      MuU)*ZH(gI2,2) + 2*g2*MDWBT*ZH(gI2,3) - 2*vT*AbsSqr(LamTU)*ZH(gI2,3) -
      1.4142135623730951*LamTU*vS*Conj(LamSU)*ZH(gI2,3) - 2*MuU*Conj(LamTU)*ZH(gI2
      ,3) - 1.4142135623730951*LamSU*vS*Conj(LamTU)*ZH(gI2,3) + 2*g2*Conj(MDWBT)*
      ZH(gI2,3) - 2*LamTU*Conj(MuU)*ZH(gI2,3));

   return result;
}

std::complex<double> MRSSMEFTHiggs_effective_couplings::CphhconjVWmVWm(int gI2) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vT = MODELPARAMETER(vT);
   const auto vu = MODELPARAMETER(vu);

   const std::complex<double> result = 0.5*Sqr(g2)*(vd*ZH(gI2,0) + vu*ZH(gI2,1) +
      4*vT*ZH(gI2,3));

   return result;
}

std::complex<double> MRSSMEFTHiggs_effective_couplings::CpbarFeFeAhPL(int gO2, int gI1, int gI2) const
{
   const auto Ye = MODELPARAMETER(Ye);

   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*SUM(j2,0,2,Conj(ZEL(gI1,j2))*SUM(j1,0,2,Conj(ZER(gO2,j1))*Ye(j1,j2)))*ZA(
      gI2,0);

   return result;
}

std::complex<double> MRSSMEFTHiggs_effective_couplings::CpbarFeFehhPL(int gO2, int gI2, int gI1) const
{
   const auto Ye = MODELPARAMETER(Ye);

   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,Conj(ZEL(gI2
      ,j2))*SUM(j1,0,2,Conj(ZER(gO2,j1))*Ye(j1,j2)))*ZH(gI1,0);

   return result;
}

std::complex<double> MRSSMEFTHiggs_effective_couplings::CpbarFdFdAhPL(int gO2, int gI1, int gI2) const
{
   const auto Yd = MODELPARAMETER(Yd);

   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*SUM(j2,0,2,Conj(ZDL(gI1,j2))*SUM(j1,0,2,Conj(ZDR(gO2,j1))*Yd(j1,j2)))*ZA(
      gI2,0);

   return result;
}

std::complex<double> MRSSMEFTHiggs_effective_couplings::CpbarFdFdhhPL(int gO2, int gI2, int gI1) const
{
   const auto Yd = MODELPARAMETER(Yd);

   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,Conj(ZDL(gI2
      ,j2))*SUM(j1,0,2,Conj(ZDR(gO2,j1))*Yd(j1,j2)))*ZH(gI1,0);

   return result;
}

std::complex<double> MRSSMEFTHiggs_effective_couplings::CpbarFuFuAhPL(int gO2, int gI1, int gI2) const
{
   const auto Yu = MODELPARAMETER(Yu);

   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*SUM(j2,0,2,Conj(ZUL(gI1,j2))*SUM(j1,0,2,Conj(ZUR(gO2,j1))*Yu(j1,j2)))*ZA(
      gI2,1);

   return result;
}

std::complex<double> MRSSMEFTHiggs_effective_couplings::CpbarFuFuhhPL(int gO2, int gI2, int gI1) const
{
   const auto Yu = MODELPARAMETER(Yu);

   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,Conj(ZUL(gI2
      ,j2))*SUM(j1,0,2,Conj(ZUR(gO2,j1))*Yu(j1,j2)))*ZH(gI1,1);

   return result;
}

std::complex<double> MRSSMEFTHiggs_effective_couplings::CphhSdconjSd(int gt1, int gt2, int gt3) const
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Mu = MODELPARAMETER(Mu);

   const std::complex<double> result = 0.016666666666666666*(30*(-2*vd*SUM(j3,0,2,
      Conj(ZD(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,j1))*ZD(gt3
      ,3 + j2)))*ZH(gt1,0) - 2*vd*SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,
      2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gt3,j3))*ZH(gt1,0) + 1.4142135623730951*(
      Conj(Mu)*SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1))) +
      Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1)))*ZD(gt3,j2)))*
      ZH(gt1,1)) + 2*g1*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*(3*g1*vd*
      ZH(gt1,0) - 3*g1*vu*ZH(gt1,1) - 7.745966692414834*(MDBS + Conj(MDBS))*ZH(gt1
      ,2)) + SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*(3*vd*(Sqr(g1) + 5*Sqr(g2))*
      ZH(gt1,0) - 3*vu*(Sqr(g1) + 5*Sqr(g2))*ZH(gt1,1) - 7.745966692414834*g1*MDBS
      *ZH(gt1,2) - 7.745966692414834*g1*Conj(MDBS)*ZH(gt1,2) + 30*g2*MDWBT*ZH(gt1,
      3) + 30*g2*Conj(MDWBT)*ZH(gt1,3)));

   return result;
}

std::complex<double> MRSSMEFTHiggs_effective_couplings::CphhSuconjSu(int gt1, int gt2, int gt3) const
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Mu = MODELPARAMETER(Mu);

   const std::complex<double> result = 0.016666666666666666*(30*(
      1.4142135623730951*Conj(Mu)*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Yu(j1,j2)
      *ZU(gt3,3 + j1)))*ZH(gt1,0) + 1.4142135623730951*Mu*SUM(j2,0,2,SUM(j1,0,2,
      Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1)))*ZU(gt3,j2))*ZH(gt1,0) - 2*vu*(SUM(j3,0
      ,2,Conj(ZU(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))*ZU(
      gt3,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Conj(Yu(j1
      ,j3))*Yu(j1,j2)))*ZU(gt3,j3)))*ZH(gt1,1)) + 4*g1*SUM(j1,0,2,Conj(ZU(gt2,3 +
      j1))*ZU(gt3,3 + j1))*(-3*g1*vd*ZH(gt1,0) + 3*g1*vu*ZH(gt1,1) +
      7.745966692414834*(MDBS + Conj(MDBS))*ZH(gt1,2)) + SUM(j1,0,2,Conj(ZU(gt2,j1
      ))*ZU(gt3,j1))*(3*vd*(Sqr(g1) - 5*Sqr(g2))*ZH(gt1,0) - 3*vu*(Sqr(g1) - 5*Sqr
      (g2))*ZH(gt1,1) - 2*(3.872983346207417*g1*(MDBS + Conj(MDBS))*ZH(gt1,2) + 15
      *g2*(MDWBT + Conj(MDWBT))*ZH(gt1,3))));

   return result;
}

std::complex<double> MRSSMEFTHiggs_effective_couplings::CphhSeconjSe(int gt1, int gt2, int gt3) const
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Mu = MODELPARAMETER(Mu);

   const std::complex<double> result = 0.05*(2*(5*(-2*vd*SUM(j3,0,2,Conj(ZE(gt2,3
      + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Ye(j2,j1))*ZE(gt3,3 + j2)))*ZH(
      gt1,0) - 2*vd*SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Conj(Ye(j1,
      j3))*Ye(j1,j2)))*ZE(gt3,j3))*ZH(gt1,0) + 1.4142135623730951*(Conj(Mu)*SUM(j2
      ,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1))) + Mu*SUM(j2,0,2,
      SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1)))*ZE(gt3,j2)))*ZH(gt1,1)) +
      g1*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*(3*g1*vd*ZH(gt1,0) - 3*g1
      *vu*ZH(gt1,1) - 7.745966692414834*(MDBS + Conj(MDBS))*ZH(gt1,2))) + SUM(j1,0
      ,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*((-3*vd*Sqr(g1) + 5*vd*Sqr(g2))*ZH(gt1,0) +
      vu*(3*Sqr(g1) - 5*Sqr(g2))*ZH(gt1,1) + 2*(3.872983346207417*g1*(MDBS + Conj(
      MDBS))*ZH(gt1,2) + 5*g2*(MDWBT + Conj(MDWBT))*ZH(gt1,3))));

   return result;
}

std::complex<double> MRSSMEFTHiggs_effective_couplings::CphhHpmconjHpm(int gt1, int gt2, int gt3) const
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

   const std::complex<double> result = 0.05*(7.745966692414834*g1*MDBS*ZH(gt1,2)*
      ZP(gt2,0)*ZP(gt3,0) - 20*vS*AbsSqr(LamSD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) -
      14.142135623730951*MuD*Conj(LamSD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) +
      7.0710678118654755*LamTD*vT*Conj(LamSD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) +
      7.0710678118654755*LamSD*vT*Conj(LamTD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) +
      7.745966692414834*g1*Conj(MDBS)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) -
      14.142135623730951*LamSD*Conj(MuD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) + 10*g2*
      MDWBT*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) - 10*vT*AbsSqr(LamTD)*ZH(gt1,3)*ZP(gt2,0
      )*ZP(gt3,0) + 7.0710678118654755*LamTD*vS*Conj(LamSD)*ZH(gt1,3)*ZP(gt2,0)*ZP
      (gt3,0) + 10*MuD*Conj(LamTD)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) +
      7.0710678118654755*LamSD*vS*Conj(LamTD)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) + 10*
      g2*Conj(MDWBT)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) + 10*LamTD*Conj(MuD)*ZH(gt1,3)*
      ZP(gt2,0)*ZP(gt3,0) - 10*LamTD*vd*Conj(LamSD)*ZH(gt1,2)*ZP(gt2,2)*ZP(gt3,0)
      + 7.0710678118654755*vd*AbsSqr(LamTD)*ZH(gt1,3)*ZP(gt2,2)*ZP(gt3,0) -
      7.0710678118654755*vd*Sqr(g2)*ZH(gt1,3)*ZP(gt2,2)*ZP(gt3,0) - 10*LamSD*vd*
      Conj(LamTD)*ZH(gt1,2)*ZP(gt2,3)*ZP(gt3,0) - 7.0710678118654755*vd*AbsSqr(
      LamTD)*ZH(gt1,3)*ZP(gt2,3)*ZP(gt3,0) + 7.0710678118654755*vd*Sqr(g2)*ZH(gt1,
      3)*ZP(gt2,3)*ZP(gt3,0) - 7.745966692414834*g1*MDBS*ZH(gt1,2)*ZP(gt2,1)*ZP(
      gt3,1) - 20*vS*AbsSqr(LamSU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) -
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
      2*Conj(MDWBT)))*ZP(gt3,0) + 2*vd*Sqr(g2)*ZP(gt3,3))) + ZP(gt2,0)*(vd*(3*Sqr(
      g1) + 5*Sqr(g2))*ZP(gt3,0) + 5*(vu*Sqr(g2)*ZP(gt3,1) + ((2.8284271247461903*
      MuD + 2*LamSD*vS - 1.4142135623730951*LamTD*vT)*Conj(LamTD) +
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

std::complex<double> MRSSMEFTHiggs_effective_couplings::CpbarCha1Cha1hhPL(int gt3, int gt1, int gt2) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);

   const std::complex<double> result = 0.5*(-(Conj(UM1(gt3,0))*(1.4142135623730951
      *LamTD*Conj(UP1(gt1,1))*ZH(gt2,0) + 2*g2*Conj(UP1(gt1,0))*ZH(gt2,3))) - Conj
      (UM1(gt3,1))*(1.4142135623730951*g2*Conj(UP1(gt1,0))*ZH(gt2,0) + Conj(UP1(
      gt1,1))*(1.4142135623730951*LamSD*ZH(gt2,2) - LamTD*ZH(gt2,3))));

   return result;
}

std::complex<double> MRSSMEFTHiggs_effective_couplings::CpbarCha2Cha2hhPL(int gt3, int gt1, int gt2) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);

   const std::complex<double> result = 0.5*(g2*Conj(UM2(gt1,0))*(-
      1.4142135623730951*Conj(UP2(gt3,1))*ZH(gt2,1) + 2*Conj(UP2(gt3,0))*ZH(gt2,3)
      ) + Conj(UM2(gt1,1))*(1.4142135623730951*LamTU*Conj(UP2(gt3,0))*ZH(gt2,1) +
      Conj(UP2(gt3,1))*(1.4142135623730951*LamSU*ZH(gt2,2) + LamTU*ZH(gt2,3))));

   return result;
}

std::complex<double> MRSSMEFTHiggs_effective_couplings::CpAhSdconjSd(int gt1, int gt2, int gt3) const
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Mu = MODELPARAMETER(Mu);

   const std::complex<double> result = std::complex<double>(0,-0.16666666666666666
      )*(4.242640687119286*Conj(Mu)*SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Yd(j1,
      j2)*ZD(gt3,3 + j1)))*ZA(gt1,1) - 4.242640687119286*Mu*SUM(j2,0,2,SUM(j1,0,2,
      Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1)))*ZD(gt3,j2))*ZA(gt1,1) +
      0.7745966692414834*g1*MDBS*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*ZA(gt1,2)
      - 0.7745966692414834*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*
      ZA(gt1,2) + 1.5491933384829668*g1*MDBS*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(
      gt3,3 + j1))*ZA(gt1,2) - 1.5491933384829668*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZD
      (gt2,3 + j1))*ZD(gt3,3 + j1))*ZA(gt1,2) - 3*g2*MDWBT*SUM(j1,0,2,Conj(ZD(gt2,
      j1))*ZD(gt3,j1))*ZA(gt1,3) + 3*g2*Conj(MDWBT)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD
      (gt3,j1))*ZA(gt1,3));

   return result;
}

std::complex<double> MRSSMEFTHiggs_effective_couplings::CpAhSuconjSu(int gt1, int gt2, int gt3) const
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Mu = MODELPARAMETER(Mu);

   const std::complex<double> result = std::complex<double>(0,-0.16666666666666666
      )*(4.242640687119286*Conj(Mu)*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Yu(j1,
      j2)*ZU(gt3,3 + j1)))*ZA(gt1,0) - 4.242640687119286*Mu*SUM(j2,0,2,SUM(j1,0,2,
      Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1)))*ZU(gt3,j2))*ZA(gt1,0) +
      0.7745966692414834*g1*MDBS*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*ZA(gt1,2)
      - 0.7745966692414834*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*
      ZA(gt1,2) - 3.0983866769659336*g1*MDBS*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(
      gt3,3 + j1))*ZA(gt1,2) + 3.0983866769659336*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZU
      (gt2,3 + j1))*ZU(gt3,3 + j1))*ZA(gt1,2) + 3*g2*MDWBT*SUM(j1,0,2,Conj(ZU(gt2,
      j1))*ZU(gt3,j1))*ZA(gt1,3) - 3*g2*Conj(MDWBT)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU
      (gt3,j1))*ZA(gt1,3));

   return result;
}

std::complex<double> MRSSMEFTHiggs_effective_couplings::CpAhSeconjSe(int gt1, int gt2, int gt3) const
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Mu = MODELPARAMETER(Mu);

   const std::complex<double> result = std::complex<double>(0,-0.1)*(
      7.0710678118654755*Conj(Mu)*SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)
      *ZE(gt3,3 + j1)))*ZA(gt1,1) - 7.0710678118654755*Mu*SUM(j2,0,2,SUM(j1,0,2,
      Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1)))*ZE(gt3,j2))*ZA(gt1,1) -
      3.872983346207417*g1*MDBS*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*ZA(gt1,2)
      + 3.872983346207417*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*ZA
      (gt1,2) + 7.745966692414834*g1*MDBS*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3
       + j1))*ZA(gt1,2) - 7.745966692414834*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZE(gt2,3
       + j1))*ZE(gt3,3 + j1))*ZA(gt1,2) - 5*g2*MDWBT*SUM(j1,0,2,Conj(ZE(gt2,j1))*
      ZE(gt3,j1))*ZA(gt1,3) + 5*g2*Conj(MDWBT)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,
      j1))*ZA(gt1,3));

   return result;
}

std::complex<double> MRSSMEFTHiggs_effective_couplings::CpAhHpmconjHpm(int gt1, int gt2, int gt3) const
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

   const std::complex<double> result = std::complex<double>(0,0.05)*(ZA(gt1,2)*(10
      *LamTD*vd*Conj(LamSD)*ZP(gt2,2)*ZP(gt3,0) - 10*LamSD*vd*Conj(LamTD)*ZP(gt2,3
      )*ZP(gt3,0) - 7.745966692414834*g1*MDBS*ZP(gt2,1)*ZP(gt3,1) +
      14.142135623730951*MuU*Conj(LamSU)*ZP(gt2,1)*ZP(gt3,1) + 7.0710678118654755*
      LamTU*vT*Conj(LamSU)*ZP(gt2,1)*ZP(gt3,1) - 7.0710678118654755*LamSU*vT*Conj(
      LamTU)*ZP(gt2,1)*ZP(gt3,1) + 7.745966692414834*g1*Conj(MDBS)*ZP(gt2,1)*ZP(
      gt3,1) - 14.142135623730951*LamSU*Conj(MuU)*ZP(gt2,1)*ZP(gt3,1) + 10*LamTU*
      vu*Conj(LamSU)*ZP(gt2,2)*ZP(gt3,1) - 10*LamSU*vu*Conj(LamTU)*ZP(gt2,3)*ZP(
      gt3,1) - 10*LamSU*vu*Conj(LamTU)*ZP(gt2,1)*ZP(gt3,2) + 10*LamTU*vu*Conj(
      LamSU)*ZP(gt2,1)*ZP(gt3,3) + ZP(gt2,0)*((7.745966692414834*g1*MDBS +
      7.0710678118654755*(2*MuD - LamTD*vT)*Conj(LamSD) + 7.0710678118654755*LamSD
      *vT*Conj(LamTD) - 7.745966692414834*g1*Conj(MDBS) - 14.142135623730951*LamSD
      *Conj(MuD))*ZP(gt3,0) + 10*vd*(-(LamSD*Conj(LamTD)*ZP(gt3,2)) + LamTD*Conj(
      LamSD)*ZP(gt3,3)))) - 5*(ZA(gt1,0)*(vu*Sqr(g2)*ZP(gt2,1)*ZP(gt3,0) + (2*
      LamTD*vS*Conj(LamSD) + 1.4142135623730951*(2*g2*MDWBT - vT*AbsSqr(LamTD) + 2
      *LamTD*Conj(MuD) + vT*Sqr(g2)))*ZP(gt2,2)*ZP(gt3,0) + 1.4142135623730951*vT*
      AbsSqr(LamTD)*ZP(gt2,3)*ZP(gt3,0) + 2.8284271247461903*MuD*Conj(LamTD)*ZP(
      gt2,3)*ZP(gt3,0) + 2*LamSD*vS*Conj(LamTD)*ZP(gt2,3)*ZP(gt3,0) +
      2.8284271247461903*g2*Conj(MDWBT)*ZP(gt2,3)*ZP(gt3,0) - 1.4142135623730951*
      vT*Sqr(g2)*ZP(gt2,3)*ZP(gt3,0) - vu*Sqr(g2)*ZP(gt2,0)*ZP(gt3,1) +
      1.4142135623730951*vT*AbsSqr(LamTD)*ZP(gt2,0)*ZP(gt3,2) - 2.8284271247461903
      *MuD*Conj(LamTD)*ZP(gt2,0)*ZP(gt3,2) - 2*LamSD*vS*Conj(LamTD)*ZP(gt2,0)*ZP(
      gt3,2) - 2.8284271247461903*g2*Conj(MDWBT)*ZP(gt2,0)*ZP(gt3,2) -
      1.4142135623730951*vT*Sqr(g2)*ZP(gt2,0)*ZP(gt3,2) - 2.8284271247461903*g2*
      MDWBT*ZP(gt2,0)*ZP(gt3,3) - 1.4142135623730951*vT*AbsSqr(LamTD)*ZP(gt2,0)*ZP
      (gt3,3) - 2*LamTD*vS*Conj(LamSD)*ZP(gt2,0)*ZP(gt3,3) - 2.8284271247461903*
      LamTD*Conj(MuD)*ZP(gt2,0)*ZP(gt3,3) + 1.4142135623730951*vT*Sqr(g2)*ZP(gt2,0
      )*ZP(gt3,3)) + ZA(gt1,1)*(-((vd*Sqr(g2)*ZP(gt2,0) + (2*LamTU*vS*Conj(LamSU)
      + 1.4142135623730951*(2*g2*MDWBT - vT*AbsSqr(LamTU) + 2*LamTU*Conj(MuU) + vT
      *Sqr(g2)))*ZP(gt2,2) + ((2.8284271247461903*MuU + 2*LamSU*vS +
      1.4142135623730951*LamTU*vT)*Conj(LamTU) + 1.4142135623730951*g2*(-(g2*vT) +
      2*Conj(MDWBT)))*ZP(gt2,3))*ZP(gt3,1)) + ZP(gt2,1)*(vd*Sqr(g2)*ZP(gt3,0) + ((
      2.8284271247461903*MuU + 2*LamSU*vS - 1.4142135623730951*LamTU*vT)*Conj(
      LamTU) + 1.4142135623730951*g2*(g2*vT + 2*Conj(MDWBT)))*ZP(gt3,2) + (2*LamTU
      *vS*Conj(LamSU) + 1.4142135623730951*(2*g2*MDWBT + vT*AbsSqr(LamTU) + 2*
      LamTU*Conj(MuU) - vT*Sqr(g2)))*ZP(gt3,3))) + ZA(gt1,3)*(1.4142135623730951*
      vd*AbsSqr(LamTD)*ZP(gt2,3)*ZP(gt3,0) - 1.4142135623730951*vd*Sqr(g2)*ZP(gt2,
      3)*ZP(gt3,0) + 2*g2*MDWBT*ZP(gt2,1)*ZP(gt3,1) + 1.4142135623730951*LamTU*vS*
      Conj(LamSU)*ZP(gt2,1)*ZP(gt3,1) - 2*MuU*Conj(LamTU)*ZP(gt2,1)*ZP(gt3,1) -
      1.4142135623730951*LamSU*vS*Conj(LamTU)*ZP(gt2,1)*ZP(gt3,1) - 2*g2*Conj(
      MDWBT)*ZP(gt2,1)*ZP(gt3,1) + 2*LamTU*Conj(MuU)*ZP(gt2,1)*ZP(gt3,1) +
      1.4142135623730951*vu*AbsSqr(LamTU)*ZP(gt2,3)*ZP(gt3,1) - 1.4142135623730951
      *vu*Sqr(g2)*ZP(gt2,3)*ZP(gt3,1) - 1.4142135623730951*vu*AbsSqr(LamTU)*ZP(gt2
      ,1)*ZP(gt3,2) + 1.4142135623730951*vu*Sqr(g2)*ZP(gt2,1)*ZP(gt3,2) - 4*vT*Sqr
      (g2)*ZP(gt2,3)*ZP(gt3,2) - 1.4142135623730951*vu*AbsSqr(LamTU)*ZP(gt2,1)*ZP(
      gt3,3) + 1.4142135623730951*vu*Sqr(g2)*ZP(gt2,1)*ZP(gt3,3) + ZP(gt2,2)*(
      1.4142135623730951*vd*(AbsSqr(LamTD) - Sqr(g2))*ZP(gt3,0) +
      1.4142135623730951*vu*(AbsSqr(LamTU) - Sqr(g2))*ZP(gt3,1) + 4*vT*Sqr(g2)*ZP(
      gt3,3)) + ZP(gt2,0)*((-2*g2*MDWBT - 1.4142135623730951*LamTD*vS*Conj(LamSD)
      + 2*MuD*Conj(LamTD) + 1.4142135623730951*LamSD*vS*Conj(LamTD) + 2*g2*Conj(
      MDWBT) - 2*LamTD*Conj(MuD))*ZP(gt3,0) + 1.4142135623730951*vd*(-AbsSqr(LamTD
      ) + Sqr(g2))*(ZP(gt3,2) + ZP(gt3,3))))));

   return result;
}

std::complex<double> MRSSMEFTHiggs_effective_couplings::CpbarCha1Cha1AhPL(int gt3, int gt2, int gt1) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);

   const std::complex<double> result = std::complex<double>(0,0.5)*(Conj(UM1(gt3,0
      ))*(-1.4142135623730951*LamTD*Conj(UP1(gt2,1))*ZA(gt1,0) + 2*g2*Conj(UP1(gt2
      ,0))*ZA(gt1,3)) + Conj(UM1(gt3,1))*(1.4142135623730951*g2*Conj(UP1(gt2,0))*
      ZA(gt1,0) + Conj(UP1(gt2,1))*(-1.4142135623730951*LamSD*ZA(gt1,2) + LamTD*ZA
      (gt1,3))));

   return result;
}

std::complex<double> MRSSMEFTHiggs_effective_couplings::CpbarCha2Cha2AhPL(int gt3, int gt2, int gt1) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);

   const std::complex<double> result = std::complex<double>(0,0.5)*(g2*Conj(UM2(
      gt2,0))*(1.4142135623730951*Conj(UP2(gt3,1))*ZA(gt1,1) - 2*Conj(UP2(gt3,0))*
      ZA(gt1,3)) + Conj(UM2(gt2,1))*(1.4142135623730951*LamTU*Conj(UP2(gt3,0))*ZA(
      gt1,1) + Conj(UP2(gt3,1))*(1.4142135623730951*LamSU*ZA(gt1,2) + LamTU*ZA(gt1
      ,3))));

   return result;
}

void MRSSMEFTHiggs_effective_couplings::calculate_eff_CphhVPVP(int gO1)
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
   result += 0.5 * CpSRdphhconjSRdp(gO1) * vev * AS0(decay_scale / Sqr(MSRdp)) /
      Sqr(MSRdp);
   result += 0.5 * CpSRumhhconjSRum(gO1) * vev * AS0(decay_scale / Sqr(MSRum)) /
      Sqr(MSRum);
   for (int gI1 = 0; gI1 < 6; ++gI1) {
      result += 0.16666666666666666 * scalar_scalar_qcd_factor(decay_mass, MSd(gI1
         )) * CphhSdconjSd(gO1, gI1, gI1) * vev * AS0(decay_scale / Sqr(MSd(gI1)))
         / Sqr(MSd(gI1));
   }
   for (int gI1 = 0; gI1 < 6; ++gI1) {
      result += 0.6666666666666666 * scalar_scalar_qcd_factor(decay_mass, MSu(gI1)
         ) * CphhSuconjSu(gO1, gI1, gI1) * vev * AS0(decay_scale / Sqr(MSu(gI1)))
         / Sqr(MSu(gI1));
   }
   for (int gI1 = 0; gI1 < 6; ++gI1) {
      result += 0.5 * CphhSeconjSe(gO1, gI1, gI1) * vev * AS0(decay_scale / Sqr(
         MSe(gI1))) / Sqr(MSe(gI1));
   }
   for (int gI1 = 1; gI1 < 4; ++gI1) {
      result += 0.5 * CphhHpmconjHpm(gO1, gI1, gI1) * vev * AS0(decay_scale / Sqr(
         MHpm(gI1))) / Sqr(MHpm(gI1));
   }
   for (int gI1 = 0; gI1 < 2; ++gI1) {
      result += CpbarCha1Cha1hhPL(gI1, gI1, gO1) * vev * AS12(decay_scale / Sqr(
         MCha1(gI1))) / MCha1(gI1);
   }
   for (int gI1 = 0; gI1 < 2; ++gI1) {
      result += CpbarCha2Cha2hhPL(gI1, gI1, gO1) * vev * AS12(decay_scale / Sqr(
         MCha2(gI1))) / MCha2(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpbarFeFehhPL(gI1, gI1, gO1) * vev * AS12(decay_scale / Sqr(MFe(
         gI1))) / MFe(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += 0.3333333333333333 * scalar_fermion_qcd_factor(decay_mass, MFd(gI1
         )) * CpbarFdFdhhPL(gI1, gI1, gO1) * vev * AS12(decay_scale / Sqr(MFd(gI1)
         )) / MFd(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += 1.3333333333333333 * scalar_fermion_qcd_factor(decay_mass, MFu(gI1
         )) * CpbarFuFuhhPL(gI1, gI1, gO1) * vev * AS12(decay_scale / Sqr(MFu(gI1)
         )) / MFu(gI1);
   }
   result += -0.5 * CphhconjVWmVWm(gO1) * vev * AS1(decay_scale / Sqr(MVWm)) / Sqr
      (MVWm);

   result *= 0.1892681907127351 * physical_input.get(Physical_input::alpha_em_0) *
      Sqrt(qedqcd.displayFermiConstant());

   ZH = saved_ZH;
   eff_CphhVPVP(gO1) = result;

}

void MRSSMEFTHiggs_effective_couplings::calculate_eff_CphhVGVG(int gO1)
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
   for (int gI1 = 0; gI1 < 6; ++gI1) {
      result += 0.5 * CphhSdconjSd(gO1, gI1, gI1) * vev * AS0(decay_scale / Sqr(
         MSd(gI1))) / Sqr(MSd(gI1));
   }
   for (int gI1 = 0; gI1 < 6; ++gI1) {
      result += 0.5 * CphhSuconjSu(gO1, gI1, gI1) * vev * AS0(decay_scale / Sqr(
         MSu(gI1))) / Sqr(MSu(gI1));
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpbarFdFdhhPL(gI1, gI1, gO1) * vev * AS12(decay_scale / Sqr(MFd(
         gI1))) / MFd(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpbarFuFuhhPL(gI1, gI1, gO1) * vev * AS12(decay_scale / Sqr(MFu(
         gI1))) / MFu(gI1);
   }
   result *= 0.75;

   if (include_qcd_corrections) {
      result *= scalar_scaling_factor(decay_mass);
   }

   result *= 0.12617879380849006 * alpha_s * Sqrt(qedqcd.displayFermiConstant());

   ZH = saved_ZH;
   eff_CphhVGVG(gO1) = result;

}

void MRSSMEFTHiggs_effective_couplings::calculate_eff_CpAhVPVP(int gO1)
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
   for (int gI1 = 0; gI1 < 2; ++gI1) {
      result += CpbarCha1Cha1AhPL(gI1, gI1, gO1) * vev * AP12(decay_scale / Sqr(
         MCha1(gI1))) / MCha1(gI1);
   }
   for (int gI1 = 0; gI1 < 2; ++gI1) {
      result += CpbarCha2Cha2AhPL(gI1, gI1, gO1) * vev * AP12(decay_scale / Sqr(
         MCha2(gI1))) / MCha2(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpbarFeFeAhPL(gI1, gI1, gO1) * vev * AP12(decay_scale / Sqr(MFe(
         gI1))) / MFe(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += 0.3333333333333333 * pseudoscalar_fermion_qcd_factor(decay_mass,
         MFd(gI1)) * CpbarFdFdAhPL(gI1, gI1, gO1) * vev * AP12(decay_scale / Sqr(
         MFd(gI1))) / MFd(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += 1.3333333333333333 * pseudoscalar_fermion_qcd_factor(decay_mass,
         MFu(gI1)) * CpbarFuFuAhPL(gI1, gI1, gO1) * vev * AP12(decay_scale / Sqr(
         MFu(gI1))) / MFu(gI1);
   }
   result *= 2.0;

   result *= 0.1892681907127351 * physical_input.get(Physical_input::alpha_em_0) *
      Sqrt(qedqcd.displayFermiConstant());

   ZA = saved_ZA;
   eff_CpAhVPVP(gO1) = result;

}

void MRSSMEFTHiggs_effective_couplings::calculate_eff_CpAhVGVG(int gO1)
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
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpbarFdFdAhPL(gI1, gI1, gO1) * vev * AP12(decay_scale / Sqr(MFd(
         gI1))) / MFd(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpbarFuFuAhPL(gI1, gI1, gO1) * vev * AP12(decay_scale / Sqr(MFu(
         gI1))) / MFu(gI1);
   }
   result *= 1.5;

   if (include_qcd_corrections) {
      result *= pseudoscalar_scaling_factor(decay_mass);
   }

   result *= 0.12617879380849006 * alpha_s * Sqrt(qedqcd.displayFermiConstant());

   ZA = saved_ZA;
   eff_CpAhVGVG(gO1) = result;

}


} // namespace flexiblesusy
