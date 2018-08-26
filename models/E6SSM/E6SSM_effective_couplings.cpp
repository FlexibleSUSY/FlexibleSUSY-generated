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

// File generated at Sun 26 Aug 2018 14:28:27

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
   , ZD(MODELPARAMETER(ZD)), ZV(MODELPARAMETER(ZV)), ZU(MODELPARAMETER(ZU)), ZE(
      MODELPARAMETER(ZE)), ZDX(MODELPARAMETER(ZDX)), ZH(MODELPARAMETER(ZH)), ZA(
      MODELPARAMETER(ZA)), ZP(MODELPARAMETER(ZP)), ZN(MODELPARAMETER(ZN)), UM(
      MODELPARAMETER(UM)), UP(MODELPARAMETER(UP)), ZEL(MODELPARAMETER(ZEL)), ZER(
      MODELPARAMETER(ZER)), ZDL(MODELPARAMETER(ZDL)), ZDR(MODELPARAMETER(ZDR)),
      ZUL(MODELPARAMETER(ZUL)), ZUR(MODELPARAMETER(ZUR)), ZDXL(MODELPARAMETER(ZDXL
      )), ZDXR(MODELPARAMETER(ZDXR)), UHI0(MODELPARAMETER(UHI0)), UHIp(
      MODELPARAMETER(UHIp)), ZMI(MODELPARAMETER(ZMI)), ZPI(MODELPARAMETER(ZPI)),
      ZNI(MODELPARAMETER(ZNI)), ZSSI(MODELPARAMETER(ZSSI)), ZFSI(MODELPARAMETER(
      ZFSI)), UHp0(MODELPARAMETER(UHp0)), UHpp(MODELPARAMETER(UHpp)), ZNp(
      MODELPARAMETER(ZNp)), ZZ(MODELPARAMETER(ZZ))
   , eff_CphhVPVP(Eigen::Array<std::complex<double>,3,1>::Zero()), eff_CphhVGVG(
      Eigen::Array<std::complex<double>,3,1>::Zero()), eff_CpAhVPVP(Eigen::Array<
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
   for (int gO1 = 0; gO1 < 3; ++gO1) {
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
   for (int gO1 = 2; gO1 < 3; ++gO1) {
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
      result = 1.0 + 0.025330295910584444*Sqr(g3) * scalar_diphoton_fermion_loop(
         m_decay, m_loop);

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

   const double nlo_qcd = 0.025330295910584444*(23.75 - 1.1666666666666667*Nf)*Sqr
      (g3);
   const double nnlo_qcd = 0.000641623890917771*Quad(g3)*(370.1956513893174 +
      2.375*l + (-47.18640261449638 + 0.6666666666666666*l)*Nf +
      0.9017702481178881*Sqr(Nf));
   const double nnnlo_qcd = 0.000016252523020247696*Power6(g3)*(467.683620788 +
      122.440972222*l + 10.9409722222*Sqr(l));

   return Sqrt(1.0 + nlo_qcd + nnlo_qcd + nnnlo_qcd);
}

double E6SSM_effective_couplings::pseudoscalar_scaling_factor(double m) const
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

double E6SSM_effective_couplings::get_hhVPVP_partial_width(int gO1) const
{
   const double mass = PHYSICAL(Mhh)(gO1);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CphhVPVP(gO1));
}

double E6SSM_effective_couplings::get_hhVGVG_partial_width(int gO1) const
{
   const double mass = PHYSICAL(Mhh)(gO1);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CphhVGVG(gO1));
}

double E6SSM_effective_couplings::get_AhVPVP_partial_width(int gO1) const
{
   const double mass = PHYSICAL(MAh)(gO1);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CpAhVPVP(gO1));
}

double E6SSM_effective_couplings::get_AhVGVG_partial_width(int gO1) const
{
   const double mass = PHYSICAL(MAh)(gO1);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CpAhVGVG(gO1));
}

std::complex<double> E6SSM_effective_couplings::CphhconjVWmVWm(int gI2) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);

   const std::complex<double> result = 0.5*Sqr(g2)*(vd*ZH(gI2,0) + vu*ZH(gI2,1));

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpbarFeFeAhPL(int gO2, int gI1, int gI2) const
{
   const auto Ye = MODELPARAMETER(Ye);

   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*SUM(j1,0,2,Conj(ZEL(gI1,j1))*Conj(ZER(gO2,j1))*Ye(j1,j1))*ZA(gI2,0);

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpbarFeFehhPL(int gO2, int gI2, int gI1) const
{
   const auto Ye = MODELPARAMETER(Ye);

   const std::complex<double> result = -0.7071067811865475*SUM(j1,0,2,Conj(ZEL(gI2
      ,j1))*Conj(ZER(gO2,j1))*Ye(j1,j1))*ZH(gI1,0);

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpbarFdFdAhPL(int gO2, int gI1, int gI2) const
{
   const auto Yd = MODELPARAMETER(Yd);

   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*SUM(j1,0,2,Conj(ZDL(gI1,j1))*Conj(ZDR(gO2,j1))*Yd(j1,j1))*ZA(gI2,0);

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpbarFdFdhhPL(int gO2, int gI2, int gI1) const
{
   const auto Yd = MODELPARAMETER(Yd);

   const std::complex<double> result = -0.7071067811865475*SUM(j1,0,2,Conj(ZDL(gI2
      ,j1))*Conj(ZDR(gO2,j1))*Yd(j1,j1))*ZH(gI1,0);

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpbarFuFuAhPL(int gO2, int gI1, int gI2) const
{
   const auto Yu = MODELPARAMETER(Yu);

   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*SUM(j1,0,2,Conj(ZUL(gI1,j1))*Conj(ZUR(gO2,j1))*Yu(j1,j1))*ZA(gI2,1);

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpbarFuFuhhPL(int gO2, int gI2, int gI1) const
{
   const auto Yu = MODELPARAMETER(Yu);

   const std::complex<double> result = -0.7071067811865475*SUM(j1,0,2,Conj(ZUL(gI2
      ,j1))*Conj(ZUR(gO2,j1))*Yu(j1,j1))*ZH(gI1,1);

   return result;
}

std::complex<double> E6SSM_effective_couplings::CphhSdconjSd(int gt1, int gt2, int gt3) const
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

   const std::complex<double> result = 0.008333333333333333*(SUM(j1,0,2,Conj(ZD(
      gt2,j1))*ZD(gt3,j1))*(vd*(6*Sqr(g1) + 30*Sqr(g2) + 9*Sqr(gN))*ZH(gt1,0) - 6*
      vu*(Sqr(g1) + 5*Sqr(g2) - Sqr(gN))*ZH(gt1,1) - 15*vs*Sqr(gN)*ZH(gt1,2)) + 2*
      (-42.42640687119285*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j1))*ZD(gt3,
      j1))*ZH(gt1,0) + SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*(vd*(6*Sqr(
      g1) + 9*Sqr(gN))*ZH(gt1,0) + 6*vu*(-Sqr(g1) + Sqr(gN))*ZH(gt1,1) - 15*vs*Sqr
      (gN)*ZH(gt1,2)) - 30*(1.4142135623730951*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,
      3 + j1)*TYd(j1,j1))*ZH(gt1,0) + 2*vd*SUM(j2,0,2,AbsSqr(Yd(j2,j2))*Conj(ZD(
      gt2,j2))*ZD(gt3,j2))*ZH(gt1,0) + 2*vd*SUM(j2,0,2,AbsSqr(Yd(j2,j2))*Conj(ZD(
      gt2,3 + j2))*ZD(gt3,3 + j2))*ZH(gt1,0) - vs*Lambdax*SUM(j1,0,2,Conj(Yd(j1,j1
      ))*Conj(ZD(gt2,3 + j1))*ZD(gt3,j1))*ZH(gt1,1) - vs*Conj(Lambdax)*SUM(j1,0,2,
      Conj(ZD(gt2,j1))*Yd(j1,j1)*ZD(gt3,3 + j1))*ZH(gt1,1) - vu*Lambdax*SUM(j1,0,2
      ,Conj(Yd(j1,j1))*Conj(ZD(gt2,3 + j1))*ZD(gt3,j1))*ZH(gt1,2) - vu*Conj(
      Lambdax)*SUM(j1,0,2,Conj(ZD(gt2,j1))*Yd(j1,j1)*ZD(gt3,3 + j1))*ZH(gt1,2))));

   return result;
}

std::complex<double> E6SSM_effective_couplings::CphhSuconjSu(int gt1, int gt2, int gt3) const
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

   const std::complex<double> result = 0.025*(-8*vd*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt2
      ,3 + j1))*ZU(gt3,3 + j1))*ZH(gt1,0) + 3*vd*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gt2,3
      + j1))*ZU(gt3,3 + j1))*ZH(gt1,0) + 20*vs*Conj(Lambdax)*SUM(j1,0,2,Conj(ZU(
      gt2,j1))*Yu(j1,j1)*ZU(gt3,3 + j1))*ZH(gt1,0) - 28.284271247461902*SUM(j1,0,2
      ,Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j1))*ZU(gt3,j1))*ZH(gt1,1) + 8*vu*Sqr(g1)*
      SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*ZH(gt1,1) + 2*vu*Sqr(gN)*SUM
      (j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*ZH(gt1,1) - 28.284271247461902*
      SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,3 + j1)*TYu(j1,j1))*ZH(gt1,1) - 40*vu*SUM
      (j2,0,2,AbsSqr(Yu(j2,j2))*Conj(ZU(gt2,j2))*ZU(gt3,j2))*ZH(gt1,1) - 40*vu*SUM
      (j2,0,2,AbsSqr(Yu(j2,j2))*Conj(ZU(gt2,3 + j2))*ZU(gt3,3 + j2))*ZH(gt1,1) - 5
      *vs*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*ZH(gt1,2) + 20*
      vd*Conj(Lambdax)*SUM(j1,0,2,Conj(ZU(gt2,j1))*Yu(j1,j1)*ZU(gt3,3 + j1))*ZH(
      gt1,2) + 20*Lambdax*SUM(j1,0,2,Conj(Yu(j1,j1))*Conj(ZU(gt2,3 + j1))*ZU(gt3,
      j1))*(vs*ZH(gt1,0) + vd*ZH(gt1,2)) + SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))
      *(vd*(2*Sqr(g1) - 10*Sqr(g2) + 3*Sqr(gN))*ZH(gt1,0) + 2*vu*(-Sqr(g1) + 5*Sqr
      (g2) + Sqr(gN))*ZH(gt1,1) - 5*vs*Sqr(gN)*ZH(gt1,2)));

   return result;
}

std::complex<double> E6SSM_effective_couplings::CphhSeconjSe(int gt1, int gt2, int gt3) const
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

   const std::complex<double> result = 0.025*(-28.284271247461902*SUM(j1,0,2,Conj(
      ZE(gt2,3 + j1))*Conj(TYe(j1,j1))*ZE(gt3,j1))*ZH(gt1,0) + 12*vd*Sqr(g1)*SUM(
      j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*ZH(gt1,0) + 3*vd*Sqr(gN)*SUM(j1,
      0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*ZH(gt1,0) - 28.284271247461902*SUM(
      j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,3 + j1)*TYe(j1,j1))*ZH(gt1,0) - 40*vd*SUM(j2,
      0,2,AbsSqr(Ye(j2,j2))*Conj(ZE(gt2,j2))*ZE(gt3,j2))*ZH(gt1,0) - 40*vd*SUM(j2,
      0,2,AbsSqr(Ye(j2,j2))*Conj(ZE(gt2,3 + j2))*ZE(gt3,3 + j2))*ZH(gt1,0) + 20*vs
      *Lambdax*SUM(j1,0,2,Conj(Ye(j1,j1))*Conj(ZE(gt2,3 + j1))*ZE(gt3,j1))*ZH(gt1,
      1) - 12*vu*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*ZH(gt1,1)
      + 2*vu*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*ZH(gt1,1) +
      20*vs*Conj(Lambdax)*SUM(j1,0,2,Conj(ZE(gt2,j1))*Ye(j1,j1)*ZE(gt3,3 + j1))*ZH
      (gt1,1) + 20*vu*Lambdax*SUM(j1,0,2,Conj(Ye(j1,j1))*Conj(ZE(gt2,3 + j1))*ZE(
      gt3,j1))*ZH(gt1,2) - 5*vs*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 +
      j1))*ZH(gt1,2) + 20*vu*Conj(Lambdax)*SUM(j1,0,2,Conj(ZE(gt2,j1))*Ye(j1,j1)*
      ZE(gt3,3 + j1))*ZH(gt1,2) - 2*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*(vd*(3
      *Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,0) + vu*(-3*Sqr(g1) + 5*Sqr(g2) - 2
      *Sqr(gN))*ZH(gt1,1) + 5*vs*Sqr(gN)*ZH(gt1,2)));

   return result;
}

std::complex<double> E6SSM_effective_couplings::CphhSDXconjSDX(int gt1, int gt2, int gt3) const
{
   const auto TKappa = MODELPARAMETER(TKappa);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto vd = MODELPARAMETER(vd);
   const auto vs = MODELPARAMETER(vs);
   const auto vu = MODELPARAMETER(vu);
   const auto Kappa = MODELPARAMETER(Kappa);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   const std::complex<double> result = 0.025*(4*vd*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gt2
      ,3 + j1))*ZDX(gt3,3 + j1))*ZH(gt1,0) - 9*vd*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gt2,
      3 + j1))*ZDX(gt3,3 + j1))*ZH(gt1,0) + 20*vu*Conj(Lambdax)*SUM(j1,0,2,Conj(
      ZDX(gt2,j1))*ZDX(gt3,3 + j1)*Kappa(j1,j1))*ZH(gt1,0) - 4*vu*Sqr(g1)*SUM(j1,0
      ,2,Conj(ZDX(gt2,3 + j1))*ZDX(gt3,3 + j1))*ZH(gt1,1) - 6*vu*Sqr(gN)*SUM(j1,0,
      2,Conj(ZDX(gt2,3 + j1))*ZDX(gt3,3 + j1))*ZH(gt1,1) + 20*vd*Conj(Lambdax)*SUM
      (j1,0,2,Conj(ZDX(gt2,j1))*ZDX(gt3,3 + j1)*Kappa(j1,j1))*ZH(gt1,1) + 20*
      Lambdax*SUM(j1,0,2,Conj(ZDX(gt2,3 + j1))*Conj(Kappa(j1,j1))*ZDX(gt3,j1))*(vu
      *ZH(gt1,0) + vd*ZH(gt1,1)) - 28.284271247461902*SUM(j1,0,2,Conj(ZDX(gt2,3 +
      j1))*Conj(TKappa(j1,j1))*ZDX(gt3,j1))*ZH(gt1,2) + 15*vs*Sqr(gN)*SUM(j1,0,2,
      Conj(ZDX(gt2,3 + j1))*ZDX(gt3,3 + j1))*ZH(gt1,2) - 28.284271247461902*SUM(j1
      ,0,2,Conj(ZDX(gt2,j1))*ZDX(gt3,3 + j1)*TKappa(j1,j1))*ZH(gt1,2) - 40*vs*SUM(
      j2,0,2,AbsSqr(Kappa(j2,j2))*Conj(ZDX(gt2,j2))*ZDX(gt3,j2))*ZH(gt1,2) - 40*vs
      *SUM(j2,0,2,AbsSqr(Kappa(j2,j2))*Conj(ZDX(gt2,3 + j2))*ZDX(gt3,3 + j2))*ZH(
      gt1,2) - 2*SUM(j1,0,2,Conj(ZDX(gt2,j1))*ZDX(gt3,j1))*(vd*(2*Sqr(g1) + 3*Sqr(
      gN))*ZH(gt1,0) + 2*vu*(-Sqr(g1) + Sqr(gN))*ZH(gt1,1) - 5*vs*Sqr(gN)*ZH(gt1,2
      )));

   return result;
}

std::complex<double> E6SSM_effective_couplings::CphhHpmconjHpm(int gt1, int gt2, int gt3) const
{
   const auto TLambdax = MODELPARAMETER(TLambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto vd = MODELPARAMETER(vd);
   const auto vs = MODELPARAMETER(vs);
   const auto vu = MODELPARAMETER(vu);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   const std::complex<double> result = 0.025*(2*ZH(gt1,1)*(ZP(gt2,0)*(vu*(3*Sqr(g1
      ) - 5*Sqr(g2) - 3*Sqr(gN))*ZP(gt3,0) - 5*vd*(-2*AbsSqr(Lambdax) + Sqr(g2))*
      ZP(gt3,1)) - ZP(gt2,1)*(5*vd*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gt3,0) + vu*(
      3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZP(gt3,1))) - ZH(gt1,0)*(ZP(gt2,0)*(vd*(6
      *Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*ZP(gt3,0) + 10*vu*(-2*AbsSqr(Lambdax) +
      Sqr(g2))*ZP(gt3,1)) + 2*ZP(gt2,1)*(5*vu*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(
      gt3,0) + vd*(-3*Sqr(g1) + 5*Sqr(g2) + 3*Sqr(gN))*ZP(gt3,1))) + 5*ZH(gt1,2)*(
      2*ZP(gt2,1)*(-2.8284271247461903*Conj(TLambdax)*ZP(gt3,0) + vs*(-4*AbsSqr(
      Lambdax) + Sqr(gN))*ZP(gt3,1)) + ZP(gt2,0)*(vs*(-8*AbsSqr(Lambdax) + 3*Sqr(
      gN))*ZP(gt3,0) - 5.656854249492381*TLambdax*ZP(gt3,1))));

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpbarChaChahhPL(int gt3, int gt1, int gt2) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   const std::complex<double> result = -0.7071067811865475*(g2*Conj(UM(gt1,0))*
      Conj(UP(gt3,1))*ZH(gt2,1) + Conj(UM(gt1,1))*(g2*Conj(UP(gt3,0))*ZH(gt2,0) +
      Conj(UP(gt3,1))*Lambdax*ZH(gt2,2)));

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpbarFDXFDXhhPL(int gt3, int gt1, int gt2) const
{
   const auto Kappa = MODELPARAMETER(Kappa);

   const std::complex<double> result = -0.7071067811865475*SUM(j1,0,2,Conj(ZDXL(
      gt1,j1))*Conj(ZDXR(gt3,j1))*Kappa(j1,j1))*ZH(gt2,2);

   return result;
}

std::complex<double> E6SSM_effective_couplings::CphhSHIpconjSHIp(int gt1, int gt2, int gt3) const
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

   const std::complex<double> result = 0.025*(SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(
      gt3,j1))*(vd*(-6*Sqr(g1) + 10*Sqr(g2) - 9*Sqr(gN))*ZH(gt1,0) + 2*vu*(3*Sqr(
      g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,1) + 15*vs*Sqr(gN)*ZH(gt1,2)) + 2*(10*
      Lambdax*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*Conj(Lambda12(j1,j1))*UHIp(gt3,j1)
      )*(vu*ZH(gt1,0) + vd*ZH(gt1,1)) + SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*UHIp(gt3
      ,2 + j1))*(vd*(3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,0) + vu*(-3*Sqr(g1)
      + 5*Sqr(g2) - 2*Sqr(gN))*ZH(gt1,1) + 5*vs*Sqr(gN)*ZH(gt1,2)) + 10*(Conj(
      Lambdax)*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt3,2 + j1)*Lambda12(j1,j1))*(vu
      *ZH(gt1,0) + vd*ZH(gt1,1)) - (1.4142135623730951*SUM(j1,0,1,Conj(UHIp(gt2,2
      + j1))*Conj(TLambda12(j1,j1))*UHIp(gt3,j1)) + 1.4142135623730951*SUM(j1,0,1,
      Conj(UHIp(gt2,j1))*UHIp(gt3,2 + j1)*TLambda12(j1,j1)) + 2*vs*(SUM(j2,0,1,
      AbsSqr(Lambda12(j2,j2))*Conj(UHIp(gt2,j2))*UHIp(gt3,j2)) + SUM(j2,0,1,AbsSqr
      (Lambda12(j2,j2))*Conj(UHIp(gt2,2 + j2))*UHIp(gt3,2 + j2))))*ZH(gt1,2))));

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpbarChaIChaIhhPL(int gt3, int gt1, int gt2) const
{
   const auto Lambda12 = MODELPARAMETER(Lambda12);

   const std::complex<double> result = -0.7071067811865475*SUM(j1,0,1,Conj(ZMI(gt1
      ,j1))*Conj(ZPI(gt3,j1))*Lambda12(j1,j1))*ZH(gt2,2);

   return result;
}

std::complex<double> E6SSM_effective_couplings::CphhSHppconjSHpp(int gt1, int gt2, int gt3) const
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto vd = MODELPARAMETER(vd);
   const auto vs = MODELPARAMETER(vs);
   const auto vu = MODELPARAMETER(vu);

   const std::complex<double> result = -0.05*(Conj(UHpp(gt2,0))*UHpp(gt3,0) - Conj
      (UHpp(gt2,1))*UHpp(gt3,1))*(vd*(3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,0)
      + vu*(-3*Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN))*ZH(gt1,1) + 5*vs*Sqr(gN)*ZH(gt1,2)
      );

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpAhSdconjSd(int gt1, int gt2, int gt3) const
{
   const auto TYd = MODELPARAMETER(TYd);
   const auto vs = MODELPARAMETER(vs);
   const auto vu = MODELPARAMETER(vu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   const std::complex<double> result = std::complex<double>(0,0.5)*(
      1.4142135623730951*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j1))*ZD(gt3,
      j1))*ZA(gt1,0) - 1.4142135623730951*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,3 +
      j1)*TYd(j1,j1))*ZA(gt1,0) + (Lambdax*SUM(j1,0,2,Conj(Yd(j1,j1))*Conj(ZD(gt2,
      3 + j1))*ZD(gt3,j1)) - Conj(Lambdax)*SUM(j1,0,2,Conj(ZD(gt2,j1))*Yd(j1,j1)*
      ZD(gt3,3 + j1)))*(vs*ZA(gt1,1) + vu*ZA(gt1,2)));

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpAhSuconjSu(int gt1, int gt2, int gt3) const
{
   const auto TYu = MODELPARAMETER(TYu);
   const auto vd = MODELPARAMETER(vd);
   const auto vs = MODELPARAMETER(vs);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   const std::complex<double> result = std::complex<double>(0,0.5)*(
      1.4142135623730951*(SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j1))*ZU(gt3,
      j1)) - SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,3 + j1)*TYu(j1,j1)))*ZA(gt1,1) +
      Lambdax*SUM(j1,0,2,Conj(Yu(j1,j1))*Conj(ZU(gt2,3 + j1))*ZU(gt3,j1))*(vs*ZA(
      gt1,0) + vd*ZA(gt1,2)) - Conj(Lambdax)*SUM(j1,0,2,Conj(ZU(gt2,j1))*Yu(j1,j1)
      *ZU(gt3,3 + j1))*(vs*ZA(gt1,0) + vd*ZA(gt1,2)));

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpAhSeconjSe(int gt1, int gt2, int gt3) const
{
   const auto TYe = MODELPARAMETER(TYe);
   const auto vs = MODELPARAMETER(vs);
   const auto vu = MODELPARAMETER(vu);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   const std::complex<double> result = std::complex<double>(0,0.5)*(
      1.4142135623730951*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j1))*ZE(gt3,
      j1))*ZA(gt1,0) - 1.4142135623730951*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,3 +
      j1)*TYe(j1,j1))*ZA(gt1,0) + (Lambdax*SUM(j1,0,2,Conj(Ye(j1,j1))*Conj(ZE(gt2,
      3 + j1))*ZE(gt3,j1)) - Conj(Lambdax)*SUM(j1,0,2,Conj(ZE(gt2,j1))*Ye(j1,j1)*
      ZE(gt3,3 + j1)))*(vs*ZA(gt1,1) + vu*ZA(gt1,2)));

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpAhSDXconjSDX(int gt1, int gt2, int gt3) const
{
   const auto TKappa = MODELPARAMETER(TKappa);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Kappa = MODELPARAMETER(Kappa);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   const std::complex<double> result = std::complex<double>(0,0.5)*(Lambdax*SUM(j1
      ,0,2,Conj(ZDX(gt2,3 + j1))*Conj(Kappa(j1,j1))*ZDX(gt3,j1))*(vu*ZA(gt1,0) +
      vd*ZA(gt1,1)) - Conj(Lambdax)*SUM(j1,0,2,Conj(ZDX(gt2,j1))*ZDX(gt3,3 + j1)*
      Kappa(j1,j1))*(vu*ZA(gt1,0) + vd*ZA(gt1,1)) + 1.4142135623730951*(SUM(j1,0,2
      ,Conj(ZDX(gt2,3 + j1))*Conj(TKappa(j1,j1))*ZDX(gt3,j1)) - SUM(j1,0,2,Conj(
      ZDX(gt2,j1))*ZDX(gt3,3 + j1)*TKappa(j1,j1)))*ZA(gt1,2));

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpAhHpmconjHpm(int gt1, int gt2, int gt3) const
{
   const auto TLambdax = MODELPARAMETER(TLambdax);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   const std::complex<double> result = std::complex<double>(0,-0.25)*(vu*(-2*
      AbsSqr(Lambdax) + Sqr(g2))*ZA(gt1,0)*(ZP(gt2,1)*ZP(gt3,0) - ZP(gt2,0)*ZP(gt3
      ,1)) + vd*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZA(gt1,1)*(ZP(gt2,1)*ZP(gt3,0) - ZP
      (gt2,0)*ZP(gt3,1)) + 2.8284271247461903*ZA(gt1,2)*(-(Conj(TLambdax)*ZP(gt2,1
      )*ZP(gt3,0)) + TLambdax*ZP(gt2,0)*ZP(gt3,1)));

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpbarChaChaAhPL(int gt3, int gt2, int gt1) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *(g2*Conj(UM(gt2,0))*Conj(UP(gt3,1))*ZA(gt1,1) + Conj(UM(gt2,1))*(g2*Conj(UP
      (gt3,0))*ZA(gt1,0) - Conj(UP(gt3,1))*Lambdax*ZA(gt1,2)));

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpbarFDXFDXAhPL(int gt3, int gt2, int gt1) const
{
   const auto Kappa = MODELPARAMETER(Kappa);

   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*SUM(j1,0,2,Conj(ZDXL(gt2,j1))*Conj(ZDXR(gt3,j1))*Kappa(j1,j1))*ZA(gt1,2);

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpAhSHIpconjSHIp(int gt1, int gt2, int gt3) const
{
   const auto TLambda12 = MODELPARAMETER(TLambda12);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Lambda12 = MODELPARAMETER(Lambda12);

   const std::complex<double> result = std::complex<double>(0,0.5)*(Lambdax*SUM(j1
      ,0,1,Conj(UHIp(gt2,2 + j1))*Conj(Lambda12(j1,j1))*UHIp(gt3,j1))*(vu*ZA(gt1,0
      ) + vd*ZA(gt1,1)) - Conj(Lambdax)*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt3,2 +
      j1)*Lambda12(j1,j1))*(vu*ZA(gt1,0) + vd*ZA(gt1,1)) + 1.4142135623730951*(SUM
      (j1,0,1,Conj(UHIp(gt2,2 + j1))*Conj(TLambda12(j1,j1))*UHIp(gt3,j1)) - SUM(j1
      ,0,1,Conj(UHIp(gt2,j1))*UHIp(gt3,2 + j1)*TLambda12(j1,j1)))*ZA(gt1,2));

   return result;
}

std::complex<double> E6SSM_effective_couplings::CpbarChaIChaIAhPL(int gt3, int gt2, int gt1) const
{
   const auto Lambda12 = MODELPARAMETER(Lambda12);

   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*SUM(j1,0,1,Conj(ZMI(gt2,j1))*Conj(ZPI(gt3,j1))*Lambda12(j1,j1))*ZA(gt1,2);

   return result;
}

void E6SSM_effective_couplings::calculate_eff_CphhVPVP(int gO1)
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
   for (int gI1 = 0; gI1 < 6; ++gI1) {
      result += 0.16666666666666666 * scalar_scalar_qcd_factor(decay_mass, MSDX(
         gI1)) * CphhSDXconjSDX(gO1, gI1, gI1) * vev * AS0(decay_scale / Sqr(MSDX(
         gI1))) / Sqr(MSDX(gI1));
   }
   for (int gI1 = 1; gI1 < 2; ++gI1) {
      result += 0.5 * CphhHpmconjHpm(gO1, gI1, gI1) * vev * AS0(decay_scale / Sqr(
         MHpm(gI1))) / Sqr(MHpm(gI1));
   }
   for (int gI1 = 0; gI1 < 2; ++gI1) {
      result += CpbarChaChahhPL(gI1, gI1, gO1) * vev * AS12(decay_scale / Sqr(MCha
         (gI1))) / MCha(gI1);
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
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += 0.3333333333333333 * scalar_fermion_qcd_factor(decay_mass, MFDX(
         gI1)) * CpbarFDXFDXhhPL(gI1, gI1, gO1) * vev * AS12(decay_scale / Sqr(
         MFDX(gI1))) / MFDX(gI1);
   }
   for (int gI1 = 0; gI1 < 4; ++gI1) {
      result += 0.5 * CphhSHIpconjSHIp(gO1, gI1, gI1) * vev * AS0(decay_scale /
         Sqr(MSHIp(gI1))) / Sqr(MSHIp(gI1));
   }
   for (int gI1 = 0; gI1 < 2; ++gI1) {
      result += CpbarChaIChaIhhPL(gI1, gI1, gO1) * vev * AS12(decay_scale / Sqr(
         MChaI(gI1))) / MChaI(gI1);
   }
   for (int gI1 = 0; gI1 < 2; ++gI1) {
      result += 0.5 * CphhSHppconjSHpp(gO1, gI1, gI1) * vev * AS0(decay_scale /
         Sqr(MSHpp(gI1))) / Sqr(MSHpp(gI1));
   }
   result += -0.5 * CphhconjVWmVWm(gO1) * vev * AS1(decay_scale / Sqr(MVWm)) / Sqr
      (MVWm);

   result *= 0.1892681907127351 * physical_input.get(Physical_input::alpha_em_0) *
      Sqrt(qedqcd.displayFermiConstant());

   ZH = saved_ZH;
   eff_CphhVPVP(gO1) = result;

}

void E6SSM_effective_couplings::calculate_eff_CphhVGVG(int gO1)
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
   for (int gI1 = 0; gI1 < 6; ++gI1) {
      result += 0.5 * CphhSdconjSd(gO1, gI1, gI1) * vev * AS0(decay_scale / Sqr(
         MSd(gI1))) / Sqr(MSd(gI1));
   }
   for (int gI1 = 0; gI1 < 6; ++gI1) {
      result += 0.5 * CphhSuconjSu(gO1, gI1, gI1) * vev * AS0(decay_scale / Sqr(
         MSu(gI1))) / Sqr(MSu(gI1));
   }
   for (int gI1 = 0; gI1 < 6; ++gI1) {
      result += 0.5 * CphhSDXconjSDX(gO1, gI1, gI1) * vev * AS0(decay_scale / Sqr(
         MSDX(gI1))) / Sqr(MSDX(gI1));
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpbarFdFdhhPL(gI1, gI1, gO1) * vev * AS12(decay_scale / Sqr(MFd(
         gI1))) / MFd(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpbarFuFuhhPL(gI1, gI1, gO1) * vev * AS12(decay_scale / Sqr(MFu(
         gI1))) / MFu(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpbarFDXFDXhhPL(gI1, gI1, gO1) * vev * AS12(decay_scale / Sqr(MFDX
         (gI1))) / MFDX(gI1);
   }
   result *= 0.75;

   if (include_qcd_corrections) {
      result *= scalar_scaling_factor(decay_mass);
   }

   result *= 0.12617879380849006 * alpha_s * Sqrt(qedqcd.displayFermiConstant());

   ZH = saved_ZH;
   eff_CphhVGVG(gO1) = result;

}

void E6SSM_effective_couplings::calculate_eff_CpAhVPVP(int gO1)
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
   for (int gI1 = 0; gI1 < 2; ++gI1) {
      result += CpbarChaChaAhPL(gI1, gI1, gO1) * vev * AP12(decay_scale / Sqr(MCha
         (gI1))) / MCha(gI1);
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
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += 0.3333333333333333 * pseudoscalar_fermion_qcd_factor(decay_mass,
         MFDX(gI1)) * CpbarFDXFDXAhPL(gI1, gI1, gO1) * vev * AP12(decay_scale /
         Sqr(MFDX(gI1))) / MFDX(gI1);
   }
   for (int gI1 = 0; gI1 < 2; ++gI1) {
      result += CpbarChaIChaIAhPL(gI1, gI1, gO1) * vev * AP12(decay_scale / Sqr(
         MChaI(gI1))) / MChaI(gI1);
   }
   result *= 2.0;

   result *= 0.1892681907127351 * physical_input.get(Physical_input::alpha_em_0) *
      Sqrt(qedqcd.displayFermiConstant());

   ZA = saved_ZA;
   eff_CpAhVPVP(gO1) = result;

}

void E6SSM_effective_couplings::calculate_eff_CpAhVGVG(int gO1)
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
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpbarFdFdAhPL(gI1, gI1, gO1) * vev * AP12(decay_scale / Sqr(MFd(
         gI1))) / MFd(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpbarFuFuAhPL(gI1, gI1, gO1) * vev * AP12(decay_scale / Sqr(MFu(
         gI1))) / MFu(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpbarFDXFDXAhPL(gI1, gI1, gO1) * vev * AP12(decay_scale / Sqr(MFDX
         (gI1))) / MFDX(gI1);
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
