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

// File generated at Wed 12 Apr 2017 11:55:10

#include "TMSSM_effective_couplings.hpp"

#include "effective_couplings.hpp"
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

void TMSSM_effective_couplings::calculate_effective_couplings()
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
   for (unsigned gO1 = 1; gO1 < 3; ++gO1) {
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

standard_model::Standard_model TMSSM_effective_couplings::initialise_SM() const
{
   standard_model::Standard_model sm;

   sm.set_loops(2);
   sm.set_thresholds(2);
   sm.set_physical_input(physical_input);

   sm.initialise_from_input(qedqcd);

   return sm;
}

void TMSSM_effective_couplings::run_SM_strong_coupling_to(standard_model::Standard_model sm, double m)
{
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

   std::complex<double> tmp_3503;
   std::complex<double> tmp_3504;
   std::complex<double> tmp_3505;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3505 += Conj(ZD(gt2,j1))*ZD(gt3,j1);
   }
   tmp_3504 += tmp_3505;
   tmp_3503 += ((0.6*Sqr(g1) + 3*Sqr(g2))*(vd*ZH(gt1,0) - vu*ZH(gt1,1))) *
      tmp_3504;
   std::complex<double> tmp_3506;
   std::complex<double> tmp_3507;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3507 += Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1);
   }
   tmp_3506 += tmp_3507;
   tmp_3503 += (1.2*Sqr(g1)*(vd*ZH(gt1,0) - vu*ZH(gt1,1))) * tmp_3506;
   std::complex<double> tmp_3508;
   std::complex<double> tmp_3509;
   std::complex<double> tmp_3510;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3511;
      std::complex<double> tmp_3512;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3512 += ZD(gt3,3 + j1)*TYd(j1,j2);
      }
      tmp_3511 += tmp_3512;
      tmp_3510 += (Conj(ZD(gt2,j2))) * tmp_3511;
   }
   tmp_3509 += tmp_3510;
   tmp_3508 += (-2.8284271247461903*ZH(gt1,0)) * tmp_3509;
   std::complex<double> tmp_3513;
   std::complex<double> tmp_3514;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3515;
      std::complex<double> tmp_3516;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3516 += Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_3515 += tmp_3516;
      tmp_3514 += (ZD(gt3,j2)) * tmp_3515;
   }
   tmp_3513 += tmp_3514;
   tmp_3508 += (-2.8284271247461903*ZH(gt1,0)) * tmp_3513;
   std::complex<double> tmp_3517;
   std::complex<double> tmp_3518;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3519;
      std::complex<double> tmp_3520;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3521;
         std::complex<double> tmp_3522;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3522 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_3521 += tmp_3522;
         tmp_3520 += (ZD(gt3,3 + j2)) * tmp_3521;
      }
      tmp_3519 += tmp_3520;
      tmp_3518 += (Conj(ZD(gt2,3 + j3))) * tmp_3519;
   }
   tmp_3517 += tmp_3518;
   tmp_3508 += (-4*vd*ZH(gt1,0)) * tmp_3517;
   std::complex<double> tmp_3523;
   std::complex<double> tmp_3524;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3525;
      std::complex<double> tmp_3526;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3527;
         std::complex<double> tmp_3528;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3528 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_3527 += tmp_3528;
         tmp_3526 += (Conj(ZD(gt2,j2))) * tmp_3527;
      }
      tmp_3525 += tmp_3526;
      tmp_3524 += (ZD(gt3,j3)) * tmp_3525;
   }
   tmp_3523 += tmp_3524;
   tmp_3508 += (-4*vd*ZH(gt1,0)) * tmp_3523;
   std::complex<double> tmp_3529;
   std::complex<double> tmp_3530;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3531;
      std::complex<double> tmp_3532;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3532 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3531 += tmp_3532;
      tmp_3530 += (Conj(ZD(gt2,j2))) * tmp_3531;
   }
   tmp_3529 += tmp_3530;
   tmp_3508 += (1.4142135623730951*vT*Conj(Lambdax)*ZH(gt1,1)) * tmp_3529;
   std::complex<double> tmp_3533;
   std::complex<double> tmp_3534;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3535;
      std::complex<double> tmp_3536;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3536 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3535 += tmp_3536;
      tmp_3534 += (Conj(ZD(gt2,j2))) * tmp_3535;
   }
   tmp_3533 += tmp_3534;
   tmp_3508 += (2.8284271247461903*Conj(Mu)*ZH(gt1,1)) * tmp_3533;
   std::complex<double> tmp_3537;
   std::complex<double> tmp_3538;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3539;
      std::complex<double> tmp_3540;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3540 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3539 += tmp_3540;
      tmp_3538 += (ZD(gt3,j2)) * tmp_3539;
   }
   tmp_3537 += tmp_3538;
   tmp_3508 += (1.4142135623730951*vT*Lambdax*ZH(gt1,1)) * tmp_3537;
   std::complex<double> tmp_3541;
   std::complex<double> tmp_3542;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3543;
      std::complex<double> tmp_3544;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3544 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3543 += tmp_3544;
      tmp_3542 += (ZD(gt3,j2)) * tmp_3543;
   }
   tmp_3541 += tmp_3542;
   tmp_3508 += (2.8284271247461903*Mu*ZH(gt1,1)) * tmp_3541;
   std::complex<double> tmp_3545;
   std::complex<double> tmp_3546;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3547;
      std::complex<double> tmp_3548;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3548 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3547 += tmp_3548;
      tmp_3546 += (Conj(ZD(gt2,j2))) * tmp_3547;
   }
   tmp_3545 += tmp_3546;
   tmp_3508 += (1.4142135623730951*vu*Conj(Lambdax)*ZH(gt1,2)) * tmp_3545;
   std::complex<double> tmp_3549;
   std::complex<double> tmp_3550;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3551;
      std::complex<double> tmp_3552;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3552 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3551 += tmp_3552;
      tmp_3550 += (ZD(gt3,j2)) * tmp_3551;
   }
   tmp_3549 += tmp_3550;
   tmp_3508 += (1.4142135623730951*vu*Lambdax*ZH(gt1,2)) * tmp_3549;
   tmp_3503 += (3) * tmp_3508;
   result += (0.08333333333333333) * tmp_3503;

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

   std::complex<double> tmp_3553;
   std::complex<double> tmp_3554;
   std::complex<double> tmp_3555;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3555 += Conj(ZU(gt2,j1))*ZU(gt3,j1);
   }
   tmp_3554 += tmp_3555;
   tmp_3553 += ((0.6*Sqr(g1) - 3*Sqr(g2))*(vd*ZH(gt1,0) - vu*ZH(gt1,1))) *
      tmp_3554;
   std::complex<double> tmp_3556;
   std::complex<double> tmp_3557;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3557 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_3556 += tmp_3557;
   tmp_3553 += (2.4*Sqr(g1)*(-(vd*ZH(gt1,0)) + vu*ZH(gt1,1))) * tmp_3556;
   std::complex<double> tmp_3558;
   std::complex<double> tmp_3559;
   std::complex<double> tmp_3560;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3561;
      std::complex<double> tmp_3562;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3562 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_3561 += tmp_3562;
      tmp_3560 += (Conj(ZU(gt2,j2))) * tmp_3561;
   }
   tmp_3559 += tmp_3560;
   tmp_3558 += (2.8284271247461903*Conj(Mu)*ZH(gt1,0)) * tmp_3559;
   std::complex<double> tmp_3563;
   std::complex<double> tmp_3564;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3565;
      std::complex<double> tmp_3566;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3566 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3565 += tmp_3566;
      tmp_3564 += (ZU(gt3,j2)) * tmp_3565;
   }
   tmp_3563 += tmp_3564;
   tmp_3558 += (1.4142135623730951*vT*Lambdax*ZH(gt1,0)) * tmp_3563;
   std::complex<double> tmp_3567;
   std::complex<double> tmp_3568;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3569;
      std::complex<double> tmp_3570;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3570 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3569 += tmp_3570;
      tmp_3568 += (ZU(gt3,j2)) * tmp_3569;
   }
   tmp_3567 += tmp_3568;
   tmp_3558 += (2.8284271247461903*Mu*ZH(gt1,0)) * tmp_3567;
   std::complex<double> tmp_3571;
   std::complex<double> tmp_3572;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3573;
      std::complex<double> tmp_3574;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3574 += ZU(gt3,3 + j1)*TYu(j1,j2);
      }
      tmp_3573 += tmp_3574;
      tmp_3572 += (Conj(ZU(gt2,j2))) * tmp_3573;
   }
   tmp_3571 += tmp_3572;
   tmp_3558 += (-2.8284271247461903*ZH(gt1,1)) * tmp_3571;
   std::complex<double> tmp_3575;
   std::complex<double> tmp_3576;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3577;
      std::complex<double> tmp_3578;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3578 += Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_3577 += tmp_3578;
      tmp_3576 += (ZU(gt3,j2)) * tmp_3577;
   }
   tmp_3575 += tmp_3576;
   tmp_3558 += (-2.8284271247461903*ZH(gt1,1)) * tmp_3575;
   std::complex<double> tmp_3579;
   std::complex<double> tmp_3580;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3581;
      std::complex<double> tmp_3582;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3583;
         std::complex<double> tmp_3584;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3584 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_3583 += tmp_3584;
         tmp_3582 += (ZU(gt3,3 + j2)) * tmp_3583;
      }
      tmp_3581 += tmp_3582;
      tmp_3580 += (Conj(ZU(gt2,3 + j3))) * tmp_3581;
   }
   tmp_3579 += tmp_3580;
   tmp_3558 += (-4*vu*ZH(gt1,1)) * tmp_3579;
   std::complex<double> tmp_3585;
   std::complex<double> tmp_3586;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3587;
      std::complex<double> tmp_3588;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3589;
         std::complex<double> tmp_3590;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3590 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_3589 += tmp_3590;
         tmp_3588 += (Conj(ZU(gt2,j2))) * tmp_3589;
      }
      tmp_3587 += tmp_3588;
      tmp_3586 += (ZU(gt3,j3)) * tmp_3587;
   }
   tmp_3585 += tmp_3586;
   tmp_3558 += (-4*vu*ZH(gt1,1)) * tmp_3585;
   std::complex<double> tmp_3591;
   std::complex<double> tmp_3592;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3593;
      std::complex<double> tmp_3594;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3594 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3593 += tmp_3594;
      tmp_3592 += (ZU(gt3,j2)) * tmp_3593;
   }
   tmp_3591 += tmp_3592;
   tmp_3558 += (1.4142135623730951*vd*Lambdax*ZH(gt1,2)) * tmp_3591;
   std::complex<double> tmp_3595;
   std::complex<double> tmp_3596;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3597;
      std::complex<double> tmp_3598;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3598 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_3597 += tmp_3598;
      tmp_3596 += (Conj(ZU(gt2,j2))) * tmp_3597;
   }
   tmp_3595 += tmp_3596;
   tmp_3558 += (1.4142135623730951*Conj(Lambdax)*(vT*ZH(gt1,0) + vd*ZH(gt1,2)))
      * tmp_3595;
   tmp_3553 += (3) * tmp_3558;
   result += (0.08333333333333333) * tmp_3553;

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

   std::complex<double> tmp_3599;
   std::complex<double> tmp_3600;
   std::complex<double> tmp_3601;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3602;
      std::complex<double> tmp_3603;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3603 += ZE(gt3,3 + j1)*TYe(j1,j2);
      }
      tmp_3602 += tmp_3603;
      tmp_3601 += (Conj(ZE(gt2,j2))) * tmp_3602;
   }
   tmp_3600 += tmp_3601;
   tmp_3599 += (-2.8284271247461903*ZH(gt1,0)) * tmp_3600;
   std::complex<double> tmp_3604;
   std::complex<double> tmp_3605;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3606;
      std::complex<double> tmp_3607;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3607 += Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_3606 += tmp_3607;
      tmp_3605 += (ZE(gt3,j2)) * tmp_3606;
   }
   tmp_3604 += tmp_3605;
   tmp_3599 += (-2.8284271247461903*ZH(gt1,0)) * tmp_3604;
   std::complex<double> tmp_3608;
   std::complex<double> tmp_3609;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3610;
      std::complex<double> tmp_3611;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3612;
         std::complex<double> tmp_3613;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3613 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_3612 += tmp_3613;
         tmp_3611 += (ZE(gt3,3 + j2)) * tmp_3612;
      }
      tmp_3610 += tmp_3611;
      tmp_3609 += (Conj(ZE(gt2,3 + j3))) * tmp_3610;
   }
   tmp_3608 += tmp_3609;
   tmp_3599 += (-4*vd*ZH(gt1,0)) * tmp_3608;
   std::complex<double> tmp_3614;
   std::complex<double> tmp_3615;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3616;
      std::complex<double> tmp_3617;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3618;
         std::complex<double> tmp_3619;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3619 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_3618 += tmp_3619;
         tmp_3617 += (Conj(ZE(gt2,j2))) * tmp_3618;
      }
      tmp_3616 += tmp_3617;
      tmp_3615 += (ZE(gt3,j3)) * tmp_3616;
   }
   tmp_3614 += tmp_3615;
   tmp_3599 += (-4*vd*ZH(gt1,0)) * tmp_3614;
   std::complex<double> tmp_3620;
   std::complex<double> tmp_3621;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3622;
      std::complex<double> tmp_3623;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3623 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3622 += tmp_3623;
      tmp_3621 += (Conj(ZE(gt2,j2))) * tmp_3622;
   }
   tmp_3620 += tmp_3621;
   tmp_3599 += (1.4142135623730951*vT*Conj(Lambdax)*ZH(gt1,1)) * tmp_3620;
   std::complex<double> tmp_3624;
   std::complex<double> tmp_3625;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3626;
      std::complex<double> tmp_3627;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3627 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3626 += tmp_3627;
      tmp_3625 += (Conj(ZE(gt2,j2))) * tmp_3626;
   }
   tmp_3624 += tmp_3625;
   tmp_3599 += (2.8284271247461903*Conj(Mu)*ZH(gt1,1)) * tmp_3624;
   std::complex<double> tmp_3628;
   std::complex<double> tmp_3629;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3630;
      std::complex<double> tmp_3631;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3631 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3630 += tmp_3631;
      tmp_3629 += (ZE(gt3,j2)) * tmp_3630;
   }
   tmp_3628 += tmp_3629;
   tmp_3599 += (1.4142135623730951*vT*Lambdax*ZH(gt1,1)) * tmp_3628;
   std::complex<double> tmp_3632;
   std::complex<double> tmp_3633;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3634;
      std::complex<double> tmp_3635;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3635 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3634 += tmp_3635;
      tmp_3633 += (ZE(gt3,j2)) * tmp_3634;
   }
   tmp_3632 += tmp_3633;
   tmp_3599 += (2.8284271247461903*Mu*ZH(gt1,1)) * tmp_3632;
   std::complex<double> tmp_3636;
   std::complex<double> tmp_3637;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3637 += Conj(ZE(gt2,j1))*ZE(gt3,j1);
   }
   tmp_3636 += tmp_3637;
   tmp_3599 += (-((0.6*Sqr(g1) - Sqr(g2))*(vd*ZH(gt1,0) - vu*ZH(gt1,1)))) *
      tmp_3636;
   std::complex<double> tmp_3638;
   std::complex<double> tmp_3639;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3639 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_3638 += tmp_3639;
   tmp_3599 += (1.2*Sqr(g1)*(vd*ZH(gt1,0) - vu*ZH(gt1,1))) * tmp_3638;
   std::complex<double> tmp_3640;
   std::complex<double> tmp_3641;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3642;
      std::complex<double> tmp_3643;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3643 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3642 += tmp_3643;
      tmp_3641 += (Conj(ZE(gt2,j2))) * tmp_3642;
   }
   tmp_3640 += tmp_3641;
   tmp_3599 += (1.4142135623730951*vu*Conj(Lambdax)*ZH(gt1,2)) * tmp_3640;
   std::complex<double> tmp_3644;
   std::complex<double> tmp_3645;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3646;
      std::complex<double> tmp_3647;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3647 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3646 += tmp_3647;
      tmp_3645 += (ZE(gt3,j2)) * tmp_3646;
   }
   tmp_3644 += tmp_3645;
   tmp_3599 += (1.4142135623730951*vu*Lambdax*ZH(gt1,2)) * tmp_3644;
   result += (0.25) * tmp_3599;

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

   std::complex<double> tmp_3648;
   std::complex<double> tmp_3649;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3650;
      std::complex<double> tmp_3651;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3651 += Conj(ZER(gt3,j1))*Ye(j1,j2);
      }
      tmp_3650 += tmp_3651;
      tmp_3649 += (Conj(ZEL(gt1,j2))) * tmp_3650;
   }
   tmp_3648 += tmp_3649;
   result += (-0.7071067811865475*ZH(gt2,0)) * tmp_3648;

   return result;
}

std::complex<double> TMSSM_effective_couplings::CpFdhhbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_3652;
   std::complex<double> tmp_3653;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3654;
      std::complex<double> tmp_3655;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3655 += Conj(ZDR(gt3,j1))*Yd(j1,j2);
      }
      tmp_3654 += tmp_3655;
      tmp_3653 += (Conj(ZDL(gt1,j2))) * tmp_3654;
   }
   tmp_3652 += tmp_3653;
   result += (-0.7071067811865475*ZH(gt2,0)) * tmp_3652;

   return result;
}

std::complex<double> TMSSM_effective_couplings::CpFuhhbarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_3656;
   std::complex<double> tmp_3657;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3658;
      std::complex<double> tmp_3659;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3659 += Conj(ZUR(gt3,j1))*Yu(j1,j2);
      }
      tmp_3658 += tmp_3659;
      tmp_3657 += (Conj(ZUL(gt1,j2))) * tmp_3658;
   }
   tmp_3656 += tmp_3657;
   result += (-0.7071067811865475*ZH(gt2,1)) * tmp_3656;

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

   std::complex<double> tmp_3660;
   std::complex<double> tmp_3661;
   std::complex<double> tmp_3662;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3663;
      std::complex<double> tmp_3664;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3664 += ZD(gt3,3 + j1)*TYd(j1,j2);
      }
      tmp_3663 += tmp_3664;
      tmp_3662 += (Conj(ZD(gt2,j2))) * tmp_3663;
   }
   tmp_3661 += tmp_3662;
   tmp_3660 += (2*ZA(gt1,0)) * tmp_3661;
   std::complex<double> tmp_3665;
   std::complex<double> tmp_3666;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3667;
      std::complex<double> tmp_3668;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3668 += Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_3667 += tmp_3668;
      tmp_3666 += (ZD(gt3,j2)) * tmp_3667;
   }
   tmp_3665 += tmp_3666;
   tmp_3660 += (-2*ZA(gt1,0)) * tmp_3665;
   std::complex<double> tmp_3669;
   std::complex<double> tmp_3670;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3671;
      std::complex<double> tmp_3672;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3672 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3671 += tmp_3672;
      tmp_3670 += (Conj(ZD(gt2,j2))) * tmp_3671;
   }
   tmp_3669 += tmp_3670;
   tmp_3660 += (vT*Conj(Lambdax)*ZA(gt1,1)) * tmp_3669;
   std::complex<double> tmp_3673;
   std::complex<double> tmp_3674;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3675;
      std::complex<double> tmp_3676;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3676 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3675 += tmp_3676;
      tmp_3674 += (Conj(ZD(gt2,j2))) * tmp_3675;
   }
   tmp_3673 += tmp_3674;
   tmp_3660 += (2*Conj(Mu)*ZA(gt1,1)) * tmp_3673;
   std::complex<double> tmp_3677;
   std::complex<double> tmp_3678;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3679;
      std::complex<double> tmp_3680;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3680 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3679 += tmp_3680;
      tmp_3678 += (ZD(gt3,j2)) * tmp_3679;
   }
   tmp_3677 += tmp_3678;
   tmp_3660 += (-(vT*Lambdax*ZA(gt1,1))) * tmp_3677;
   std::complex<double> tmp_3681;
   std::complex<double> tmp_3682;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3683;
      std::complex<double> tmp_3684;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3684 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3683 += tmp_3684;
      tmp_3682 += (ZD(gt3,j2)) * tmp_3683;
   }
   tmp_3681 += tmp_3682;
   tmp_3660 += (-2*Mu*ZA(gt1,1)) * tmp_3681;
   std::complex<double> tmp_3685;
   std::complex<double> tmp_3686;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3687;
      std::complex<double> tmp_3688;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3688 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3687 += tmp_3688;
      tmp_3686 += (Conj(ZD(gt2,j2))) * tmp_3687;
   }
   tmp_3685 += tmp_3686;
   tmp_3660 += (vu*Conj(Lambdax)*ZA(gt1,2)) * tmp_3685;
   std::complex<double> tmp_3689;
   std::complex<double> tmp_3690;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3691;
      std::complex<double> tmp_3692;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3692 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3691 += tmp_3692;
      tmp_3690 += (ZD(gt3,j2)) * tmp_3691;
   }
   tmp_3689 += tmp_3690;
   tmp_3660 += (-(vu*Lambdax*ZA(gt1,2))) * tmp_3689;
   result += (std::complex<double>(0.,-0.35355339059327373)) * tmp_3660;

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

   std::complex<double> tmp_3693;
   std::complex<double> tmp_3694;
   std::complex<double> tmp_3695;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3696;
      std::complex<double> tmp_3697;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3697 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_3696 += tmp_3697;
      tmp_3695 += (Conj(ZU(gt2,j2))) * tmp_3696;
   }
   tmp_3694 += tmp_3695;
   tmp_3693 += (2*Conj(Mu)*ZA(gt1,0)) * tmp_3694;
   std::complex<double> tmp_3698;
   std::complex<double> tmp_3699;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3700;
      std::complex<double> tmp_3701;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3701 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3700 += tmp_3701;
      tmp_3699 += (ZU(gt3,j2)) * tmp_3700;
   }
   tmp_3698 += tmp_3699;
   tmp_3693 += (-(vT*Lambdax*ZA(gt1,0))) * tmp_3698;
   std::complex<double> tmp_3702;
   std::complex<double> tmp_3703;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3704;
      std::complex<double> tmp_3705;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3705 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3704 += tmp_3705;
      tmp_3703 += (ZU(gt3,j2)) * tmp_3704;
   }
   tmp_3702 += tmp_3703;
   tmp_3693 += (-2*Mu*ZA(gt1,0)) * tmp_3702;
   std::complex<double> tmp_3706;
   std::complex<double> tmp_3707;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3708;
      std::complex<double> tmp_3709;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3709 += ZU(gt3,3 + j1)*TYu(j1,j2);
      }
      tmp_3708 += tmp_3709;
      tmp_3707 += (Conj(ZU(gt2,j2))) * tmp_3708;
   }
   tmp_3706 += tmp_3707;
   tmp_3693 += (2*ZA(gt1,1)) * tmp_3706;
   std::complex<double> tmp_3710;
   std::complex<double> tmp_3711;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3712;
      std::complex<double> tmp_3713;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3713 += Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_3712 += tmp_3713;
      tmp_3711 += (ZU(gt3,j2)) * tmp_3712;
   }
   tmp_3710 += tmp_3711;
   tmp_3693 += (-2*ZA(gt1,1)) * tmp_3710;
   std::complex<double> tmp_3714;
   std::complex<double> tmp_3715;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3716;
      std::complex<double> tmp_3717;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3717 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3716 += tmp_3717;
      tmp_3715 += (ZU(gt3,j2)) * tmp_3716;
   }
   tmp_3714 += tmp_3715;
   tmp_3693 += (-(vd*Lambdax*ZA(gt1,2))) * tmp_3714;
   std::complex<double> tmp_3718;
   std::complex<double> tmp_3719;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3720;
      std::complex<double> tmp_3721;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3721 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_3720 += tmp_3721;
      tmp_3719 += (Conj(ZU(gt2,j2))) * tmp_3720;
   }
   tmp_3718 += tmp_3719;
   tmp_3693 += (Conj(Lambdax)*(vT*ZA(gt1,0) + vd*ZA(gt1,2))) * tmp_3718;
   result += (std::complex<double>(0.,-0.35355339059327373)) * tmp_3693;

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

   std::complex<double> tmp_3722;
   std::complex<double> tmp_3723;
   std::complex<double> tmp_3724;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3725;
      std::complex<double> tmp_3726;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3726 += ZE(gt3,3 + j1)*TYe(j1,j2);
      }
      tmp_3725 += tmp_3726;
      tmp_3724 += (Conj(ZE(gt2,j2))) * tmp_3725;
   }
   tmp_3723 += tmp_3724;
   tmp_3722 += (2*ZA(gt1,0)) * tmp_3723;
   std::complex<double> tmp_3727;
   std::complex<double> tmp_3728;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3729;
      std::complex<double> tmp_3730;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3730 += Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_3729 += tmp_3730;
      tmp_3728 += (ZE(gt3,j2)) * tmp_3729;
   }
   tmp_3727 += tmp_3728;
   tmp_3722 += (-2*ZA(gt1,0)) * tmp_3727;
   std::complex<double> tmp_3731;
   std::complex<double> tmp_3732;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3733;
      std::complex<double> tmp_3734;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3734 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3733 += tmp_3734;
      tmp_3732 += (Conj(ZE(gt2,j2))) * tmp_3733;
   }
   tmp_3731 += tmp_3732;
   tmp_3722 += (vT*Conj(Lambdax)*ZA(gt1,1)) * tmp_3731;
   std::complex<double> tmp_3735;
   std::complex<double> tmp_3736;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3737;
      std::complex<double> tmp_3738;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3738 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3737 += tmp_3738;
      tmp_3736 += (Conj(ZE(gt2,j2))) * tmp_3737;
   }
   tmp_3735 += tmp_3736;
   tmp_3722 += (2*Conj(Mu)*ZA(gt1,1)) * tmp_3735;
   std::complex<double> tmp_3739;
   std::complex<double> tmp_3740;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3741;
      std::complex<double> tmp_3742;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3742 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3741 += tmp_3742;
      tmp_3740 += (ZE(gt3,j2)) * tmp_3741;
   }
   tmp_3739 += tmp_3740;
   tmp_3722 += (-(vT*Lambdax*ZA(gt1,1))) * tmp_3739;
   std::complex<double> tmp_3743;
   std::complex<double> tmp_3744;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3745;
      std::complex<double> tmp_3746;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3746 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3745 += tmp_3746;
      tmp_3744 += (ZE(gt3,j2)) * tmp_3745;
   }
   tmp_3743 += tmp_3744;
   tmp_3722 += (-2*Mu*ZA(gt1,1)) * tmp_3743;
   std::complex<double> tmp_3747;
   std::complex<double> tmp_3748;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3749;
      std::complex<double> tmp_3750;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3750 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3749 += tmp_3750;
      tmp_3748 += (Conj(ZE(gt2,j2))) * tmp_3749;
   }
   tmp_3747 += tmp_3748;
   tmp_3722 += (vu*Conj(Lambdax)*ZA(gt1,2)) * tmp_3747;
   std::complex<double> tmp_3751;
   std::complex<double> tmp_3752;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3753;
      std::complex<double> tmp_3754;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3754 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3753 += tmp_3754;
      tmp_3752 += (ZE(gt3,j2)) * tmp_3753;
   }
   tmp_3751 += tmp_3752;
   tmp_3722 += (-(vu*Lambdax*ZA(gt1,2))) * tmp_3751;
   result += (std::complex<double>(0.,-0.35355339059327373)) * tmp_3722;

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

   std::complex<double> tmp_3755;
   std::complex<double> tmp_3756;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3757;
      std::complex<double> tmp_3758;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3758 += Conj(ZER(gt3,j1))*Ye(j1,j2);
      }
      tmp_3757 += tmp_3758;
      tmp_3756 += (Conj(ZEL(gt2,j2))) * tmp_3757;
   }
   tmp_3755 += tmp_3756;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,0)) *
      tmp_3755;

   return result;
}

std::complex<double> TMSSM_effective_couplings::CpAhFdbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_3759;
   std::complex<double> tmp_3760;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3761;
      std::complex<double> tmp_3762;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3762 += Conj(ZDR(gt3,j1))*Yd(j1,j2);
      }
      tmp_3761 += tmp_3762;
      tmp_3760 += (Conj(ZDL(gt2,j2))) * tmp_3761;
   }
   tmp_3759 += tmp_3760;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,0)) *
      tmp_3759;

   return result;
}

std::complex<double> TMSSM_effective_couplings::CpAhFubarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_3763;
   std::complex<double> tmp_3764;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3765;
      std::complex<double> tmp_3766;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3766 += Conj(ZUR(gt3,j1))*Yu(j1,j2);
      }
      tmp_3765 += tmp_3766;
      tmp_3764 += (Conj(ZUL(gt2,j2))) * tmp_3765;
   }
   tmp_3763 += tmp_3764;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,1)) *
      tmp_3763;

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
