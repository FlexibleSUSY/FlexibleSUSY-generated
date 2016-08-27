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

// File generated at Sat 27 Aug 2016 13:02:24

#include "MSSMRHN_effective_couplings.hpp"

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

MSSMRHN_effective_couplings::MSSMRHN_effective_couplings(
   const MSSMRHN_mass_eigenstates& model_,
   const softsusy::QedQcd& qedqcd_,
   const Physical_input& input_)
   : model(model_), qedqcd(qedqcd_), physical_input(input_)
   , rg_improve(true), include_qcd_corrections(true)
   , ZD(MODELPARAMETER(ZD)), ZU(MODELPARAMETER(ZU)), ZE(MODELPARAMETER(ZE)), ZV
      (MODELPARAMETER(ZV)), ZH(MODELPARAMETER(ZH)), ZA(MODELPARAMETER(ZA)), ZP(
      MODELPARAMETER(ZP)), ZN(MODELPARAMETER(ZN)), UV(MODELPARAMETER(UV)), UM(
      MODELPARAMETER(UM)), UP(MODELPARAMETER(UP)), ZEL(MODELPARAMETER(ZEL)), ZER(
      MODELPARAMETER(ZER)), ZDL(MODELPARAMETER(ZDL)), ZDR(MODELPARAMETER(ZDR)),
      ZUL(MODELPARAMETER(ZUL)), ZUR(MODELPARAMETER(ZUR)), ZZ(MODELPARAMETER(ZZ))

   , eff_CphhVPVP(Eigen::Array<std::complex<double>,2,1>::Zero()), eff_CphhVGVG
      (Eigen::Array<std::complex<double>,2,1>::Zero()), eff_CpAhVPVP(Eigen::Array<
      std::complex<double>,2,1>::Zero()), eff_CpAhVGVG(Eigen::Array<std::complex<
      double>,2,1>::Zero())

{
}

void MSSMRHN_effective_couplings::calculate_effective_couplings()
{
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
      run_SM_strong_coupling_to(0.5 * Mhh(gO1));
      calculate_eff_CphhVPVP(gO1);
      run_SM_strong_coupling_to(Mhh(gO1));
      calculate_eff_CphhVGVG(gO1);
   }

   const auto MAh = PHYSICAL(MAh);
   for (unsigned gO1 = 1; gO1 < 2; ++gO1) {
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

void MSSMRHN_effective_couplings::set_model(const MSSMRHN_mass_eigenstates& model_)
{
   model = model_;
   copy_mixing_matrices_from_model();
}

void MSSMRHN_effective_couplings::copy_mixing_matrices_from_model()
{
   ZD = MODELPARAMETER(ZD);
   ZU = MODELPARAMETER(ZU);
   ZE = MODELPARAMETER(ZE);
   ZV = MODELPARAMETER(ZV);
   ZH = MODELPARAMETER(ZH);
   ZA = MODELPARAMETER(ZA);
   ZP = MODELPARAMETER(ZP);
   ZN = MODELPARAMETER(ZN);
   UV = MODELPARAMETER(UV);
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

void MSSMRHN_effective_couplings::run_SM_strong_coupling_to(double m)
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

std::complex<double> MSSMRHN_effective_couplings::scalar_scalar_qcd_factor(double m_decay, double m_loop) const
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

std::complex<double> MSSMRHN_effective_couplings::scalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) *
         scalar_diphoton_fermion_loop(m_decay, m_loop);

   }

   return result;
}

std::complex<double> MSSMRHN_effective_couplings::pseudoscalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) *
         pseudoscalar_diphoton_fermion_loop(m_decay, m_loop);

   }

   return result;
}

double MSSMRHN_effective_couplings::number_of_active_flavours(double m) const
{
   if (m < qedqcd.displayMbMb()) {
      return 4.0;
   } else if (m < qedqcd.displayPoleMt()) {
      return 5.0;
   } else {
      return 6.0;
   }
}

double MSSMRHN_effective_couplings::scalar_scaling_factor(double m) const
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

double MSSMRHN_effective_couplings::pseudoscalar_scaling_factor(double m) const
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

double MSSMRHN_effective_couplings::get_hhVPVP_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(Mhh)(gO1);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CphhVPVP(gO1));
}

double MSSMRHN_effective_couplings::get_hhVGVG_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(Mhh)(gO1);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CphhVGVG(gO1));
}

double MSSMRHN_effective_couplings::get_AhVPVP_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(MAh)(gO1);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CpAhVPVP(gO1));
}

double MSSMRHN_effective_couplings::get_AhVGVG_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(MAh)(gO1);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CpAhVGVG(gO1));
}

std::complex<double> MSSMRHN_effective_couplings::CphhSdconjSd(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYd = MODELPARAMETER(TYd);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3532;
   std::complex<double> tmp_3533;
   std::complex<double> tmp_3534;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3534 += Conj(ZD(gt2,j1))*ZD(gt3,j1);
   }
   tmp_3533 += tmp_3534;
   tmp_3532 += ((0.6*Sqr(g1) + 3*Sqr(g2))*(vd*ZH(gt1,0) - vu*ZH(gt1,1))) *
      tmp_3533;
   std::complex<double> tmp_3535;
   std::complex<double> tmp_3536;
   std::complex<double> tmp_3537;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3537 += Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1);
   }
   tmp_3536 += tmp_3537;
   tmp_3535 += (0.6*Sqr(g1)*(vd*ZH(gt1,0) - vu*ZH(gt1,1))) * tmp_3536;
   std::complex<double> tmp_3538;
   std::complex<double> tmp_3539;
   std::complex<double> tmp_3540;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3541;
      std::complex<double> tmp_3542;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3542 += ZD(gt3,3 + j1)*TYd(j1,j2);
      }
      tmp_3541 += tmp_3542;
      tmp_3540 += (Conj(ZD(gt2,j2))) * tmp_3541;
   }
   tmp_3539 += tmp_3540;
   tmp_3538 += (1.4142135623730951*ZH(gt1,0)) * tmp_3539;
   std::complex<double> tmp_3543;
   std::complex<double> tmp_3544;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3545;
      std::complex<double> tmp_3546;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3546 += Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_3545 += tmp_3546;
      tmp_3544 += (ZD(gt3,j2)) * tmp_3545;
   }
   tmp_3543 += tmp_3544;
   tmp_3538 += (1.4142135623730951*ZH(gt1,0)) * tmp_3543;
   std::complex<double> tmp_3547;
   std::complex<double> tmp_3548;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3549;
      std::complex<double> tmp_3550;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3551;
         std::complex<double> tmp_3552;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3552 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_3551 += tmp_3552;
         tmp_3550 += (ZD(gt3,3 + j2)) * tmp_3551;
      }
      tmp_3549 += tmp_3550;
      tmp_3548 += (Conj(ZD(gt2,3 + j3))) * tmp_3549;
   }
   tmp_3547 += tmp_3548;
   tmp_3538 += (2*vd*ZH(gt1,0)) * tmp_3547;
   std::complex<double> tmp_3553;
   std::complex<double> tmp_3554;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3555;
      std::complex<double> tmp_3556;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3557;
         std::complex<double> tmp_3558;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3558 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_3557 += tmp_3558;
         tmp_3556 += (Conj(ZD(gt2,j2))) * tmp_3557;
      }
      tmp_3555 += tmp_3556;
      tmp_3554 += (ZD(gt3,j3)) * tmp_3555;
   }
   tmp_3553 += tmp_3554;
   tmp_3538 += (2*vd*ZH(gt1,0)) * tmp_3553;
   std::complex<double> tmp_3559;
   std::complex<double> tmp_3560;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3561;
      std::complex<double> tmp_3562;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3562 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3561 += tmp_3562;
      tmp_3560 += (Conj(ZD(gt2,j2))) * tmp_3561;
   }
   tmp_3559 += tmp_3560;
   tmp_3538 += (-1.4142135623730951*Conj(Mu)*ZH(gt1,1)) * tmp_3559;
   std::complex<double> tmp_3563;
   std::complex<double> tmp_3564;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3565;
      std::complex<double> tmp_3566;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3566 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3565 += tmp_3566;
      tmp_3564 += (ZD(gt3,j2)) * tmp_3565;
   }
   tmp_3563 += tmp_3564;
   tmp_3538 += (-1.4142135623730951*Mu*ZH(gt1,1)) * tmp_3563;
   tmp_3535 += (-3) * tmp_3538;
   tmp_3532 += (2) * tmp_3535;
   result += (0.08333333333333333) * tmp_3532;

   return result;
}

std::complex<double> MSSMRHN_effective_couplings::CphhSuconjSu(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYu = MODELPARAMETER(TYu);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3567;
   std::complex<double> tmp_3568;
   std::complex<double> tmp_3569;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3569 += Conj(ZU(gt2,j1))*ZU(gt3,j1);
   }
   tmp_3568 += tmp_3569;
   tmp_3567 += ((0.6*Sqr(g1) - 3*Sqr(g2))*(vd*ZH(gt1,0) - vu*ZH(gt1,1))) *
      tmp_3568;
   std::complex<double> tmp_3570;
   std::complex<double> tmp_3571;
   std::complex<double> tmp_3572;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3573;
      std::complex<double> tmp_3574;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3574 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_3573 += tmp_3574;
      tmp_3572 += (Conj(ZU(gt2,j2))) * tmp_3573;
   }
   tmp_3571 += tmp_3572;
   tmp_3570 += (-4.242640687119286*Conj(Mu)*ZH(gt1,0)) * tmp_3571;
   std::complex<double> tmp_3575;
   std::complex<double> tmp_3576;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3577;
      std::complex<double> tmp_3578;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3578 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3577 += tmp_3578;
      tmp_3576 += (ZU(gt3,j2)) * tmp_3577;
   }
   tmp_3575 += tmp_3576;
   tmp_3570 += (-4.242640687119286*Mu*ZH(gt1,0)) * tmp_3575;
   std::complex<double> tmp_3579;
   std::complex<double> tmp_3580;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3581;
      std::complex<double> tmp_3582;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3582 += ZU(gt3,3 + j1)*TYu(j1,j2);
      }
      tmp_3581 += tmp_3582;
      tmp_3580 += (Conj(ZU(gt2,j2))) * tmp_3581;
   }
   tmp_3579 += tmp_3580;
   tmp_3570 += (4.242640687119286*ZH(gt1,1)) * tmp_3579;
   std::complex<double> tmp_3583;
   std::complex<double> tmp_3584;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3585;
      std::complex<double> tmp_3586;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3586 += Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_3585 += tmp_3586;
      tmp_3584 += (ZU(gt3,j2)) * tmp_3585;
   }
   tmp_3583 += tmp_3584;
   tmp_3570 += (4.242640687119286*ZH(gt1,1)) * tmp_3583;
   std::complex<double> tmp_3587;
   std::complex<double> tmp_3588;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3589;
      std::complex<double> tmp_3590;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3591;
         std::complex<double> tmp_3592;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3592 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_3591 += tmp_3592;
         tmp_3590 += (ZU(gt3,3 + j2)) * tmp_3591;
      }
      tmp_3589 += tmp_3590;
      tmp_3588 += (Conj(ZU(gt2,3 + j3))) * tmp_3589;
   }
   tmp_3587 += tmp_3588;
   tmp_3570 += (6*vu*ZH(gt1,1)) * tmp_3587;
   std::complex<double> tmp_3593;
   std::complex<double> tmp_3594;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3595;
      std::complex<double> tmp_3596;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3597;
         std::complex<double> tmp_3598;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3598 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_3597 += tmp_3598;
         tmp_3596 += (Conj(ZU(gt2,j2))) * tmp_3597;
      }
      tmp_3595 += tmp_3596;
      tmp_3594 += (ZU(gt3,j3)) * tmp_3595;
   }
   tmp_3593 += tmp_3594;
   tmp_3570 += (6*vu*ZH(gt1,1)) * tmp_3593;
   std::complex<double> tmp_3599;
   std::complex<double> tmp_3600;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3600 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_3599 += tmp_3600;
   tmp_3570 += (1.2*Sqr(g1)*(vd*ZH(gt1,0) - vu*ZH(gt1,1))) * tmp_3599;
   tmp_3567 += (-2) * tmp_3570;
   result += (0.08333333333333333) * tmp_3567;

   return result;
}

std::complex<double> MSSMRHN_effective_couplings::CphhSeconjSe(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYe = MODELPARAMETER(TYe);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3601;
   std::complex<double> tmp_3602;
   std::complex<double> tmp_3603;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3603 += Conj(ZE(gt2,j1))*ZE(gt3,j1);
   }
   tmp_3602 += tmp_3603;
   tmp_3601 += (-((0.6*Sqr(g1) - Sqr(g2))*(vd*ZH(gt1,0) - vu*ZH(gt1,1)))) *
      tmp_3602;
   std::complex<double> tmp_3604;
   std::complex<double> tmp_3605;
   std::complex<double> tmp_3606;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3607;
      std::complex<double> tmp_3608;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3608 += ZE(gt3,3 + j1)*TYe(j1,j2);
      }
      tmp_3607 += tmp_3608;
      tmp_3606 += (Conj(ZE(gt2,j2))) * tmp_3607;
   }
   tmp_3605 += tmp_3606;
   tmp_3604 += (-1.4142135623730951*ZH(gt1,0)) * tmp_3605;
   std::complex<double> tmp_3609;
   std::complex<double> tmp_3610;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3611;
      std::complex<double> tmp_3612;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3612 += Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_3611 += tmp_3612;
      tmp_3610 += (ZE(gt3,j2)) * tmp_3611;
   }
   tmp_3609 += tmp_3610;
   tmp_3604 += (-1.4142135623730951*ZH(gt1,0)) * tmp_3609;
   std::complex<double> tmp_3613;
   std::complex<double> tmp_3614;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3615;
      std::complex<double> tmp_3616;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3617;
         std::complex<double> tmp_3618;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3618 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_3617 += tmp_3618;
         tmp_3616 += (ZE(gt3,3 + j2)) * tmp_3617;
      }
      tmp_3615 += tmp_3616;
      tmp_3614 += (Conj(ZE(gt2,3 + j3))) * tmp_3615;
   }
   tmp_3613 += tmp_3614;
   tmp_3604 += (-2*vd*ZH(gt1,0)) * tmp_3613;
   std::complex<double> tmp_3619;
   std::complex<double> tmp_3620;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3621;
      std::complex<double> tmp_3622;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3623;
         std::complex<double> tmp_3624;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3624 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_3623 += tmp_3624;
         tmp_3622 += (Conj(ZE(gt2,j2))) * tmp_3623;
      }
      tmp_3621 += tmp_3622;
      tmp_3620 += (ZE(gt3,j3)) * tmp_3621;
   }
   tmp_3619 += tmp_3620;
   tmp_3604 += (-2*vd*ZH(gt1,0)) * tmp_3619;
   std::complex<double> tmp_3625;
   std::complex<double> tmp_3626;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3627;
      std::complex<double> tmp_3628;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3628 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3627 += tmp_3628;
      tmp_3626 += (Conj(ZE(gt2,j2))) * tmp_3627;
   }
   tmp_3625 += tmp_3626;
   tmp_3604 += (1.4142135623730951*Conj(Mu)*ZH(gt1,1)) * tmp_3625;
   std::complex<double> tmp_3629;
   std::complex<double> tmp_3630;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3631;
      std::complex<double> tmp_3632;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3632 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3631 += tmp_3632;
      tmp_3630 += (ZE(gt3,j2)) * tmp_3631;
   }
   tmp_3629 += tmp_3630;
   tmp_3604 += (1.4142135623730951*Mu*ZH(gt1,1)) * tmp_3629;
   std::complex<double> tmp_3633;
   std::complex<double> tmp_3634;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3634 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_3633 += tmp_3634;
   tmp_3604 += (0.6*Sqr(g1)*(vd*ZH(gt1,0) - vu*ZH(gt1,1))) * tmp_3633;
   tmp_3601 += (2) * tmp_3604;
   result += (0.25) * tmp_3601;

   return result;
}

std::complex<double> MSSMRHN_effective_couplings::CphhHpmconjHpm(unsigned gt1, unsigned gt2, unsigned gt3) const
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

std::complex<double> MSSMRHN_effective_couplings::CpChahhbarChaPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g2 = MODELPARAMETER(g2);

   std::complex<double> result;

   result = -0.7071067811865475*g2*(Conj(UM(gt1,1))*Conj(UP(gt3,0))*ZH(gt2,0) +
      Conj(UM(gt1,0))*Conj(UP(gt3,1))*ZH(gt2,1));

   return result;
}

std::complex<double> MSSMRHN_effective_couplings::CpFehhbarFePL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Ye = MODELPARAMETER(Ye);

   std::complex<double> result;

   std::complex<double> tmp_3635;
   std::complex<double> tmp_3636;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3637;
      std::complex<double> tmp_3638;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3638 += Conj(ZER(gt3,j1))*Ye(j1,j2);
      }
      tmp_3637 += tmp_3638;
      tmp_3636 += (Conj(ZEL(gt1,j2))) * tmp_3637;
   }
   tmp_3635 += tmp_3636;
   result += (-0.7071067811865475*ZH(gt2,0)) * tmp_3635;

   return result;
}

std::complex<double> MSSMRHN_effective_couplings::CpFdhhbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_3639;
   std::complex<double> tmp_3640;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3641;
      std::complex<double> tmp_3642;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3642 += Conj(ZDR(gt3,j1))*Yd(j1,j2);
      }
      tmp_3641 += tmp_3642;
      tmp_3640 += (Conj(ZDL(gt1,j2))) * tmp_3641;
   }
   tmp_3639 += tmp_3640;
   result += (-0.7071067811865475*ZH(gt2,0)) * tmp_3639;

   return result;
}

std::complex<double> MSSMRHN_effective_couplings::CpFuhhbarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_3643;
   std::complex<double> tmp_3644;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3645;
      std::complex<double> tmp_3646;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3646 += Conj(ZUR(gt3,j1))*Yu(j1,j2);
      }
      tmp_3645 += tmp_3646;
      tmp_3644 += (Conj(ZUL(gt1,j2))) * tmp_3645;
   }
   tmp_3643 += tmp_3644;
   result += (-0.7071067811865475*ZH(gt2,1)) * tmp_3643;

   return result;
}

std::complex<double> MSSMRHN_effective_couplings::CphhVWmconjVWm(unsigned gt1) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);

   std::complex<double> result;

   result = 0.5*Sqr(g2)*(vd*ZH(gt1,0) + vu*ZH(gt1,1));

   return result;
}

std::complex<double> MSSMRHN_effective_couplings::CpAhSdconjSd(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYd = MODELPARAMETER(TYd);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3647;
   std::complex<double> tmp_3648;
   std::complex<double> tmp_3649;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3650;
      std::complex<double> tmp_3651;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3651 += ZD(gt3,3 + j1)*TYd(j1,j2);
      }
      tmp_3650 += tmp_3651;
      tmp_3649 += (Conj(ZD(gt2,j2))) * tmp_3650;
   }
   tmp_3648 += tmp_3649;
   tmp_3647 += (ZA(gt1,0)) * tmp_3648;
   std::complex<double> tmp_3652;
   std::complex<double> tmp_3653;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3654;
      std::complex<double> tmp_3655;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3655 += Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_3654 += tmp_3655;
      tmp_3653 += (ZD(gt3,j2)) * tmp_3654;
   }
   tmp_3652 += tmp_3653;
   tmp_3647 += (-ZA(gt1,0)) * tmp_3652;
   std::complex<double> tmp_3656;
   std::complex<double> tmp_3657;
   std::complex<double> tmp_3658;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3659;
      std::complex<double> tmp_3660;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3660 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3659 += tmp_3660;
      tmp_3658 += (Conj(ZD(gt2,j2))) * tmp_3659;
   }
   tmp_3657 += tmp_3658;
   tmp_3656 += (Conj(Mu)) * tmp_3657;
   std::complex<double> tmp_3661;
   std::complex<double> tmp_3662;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3663;
      std::complex<double> tmp_3664;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3664 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3663 += tmp_3664;
      tmp_3662 += (ZD(gt3,j2)) * tmp_3663;
   }
   tmp_3661 += tmp_3662;
   tmp_3656 += (-Mu) * tmp_3661;
   tmp_3647 += (ZA(gt1,1)) * tmp_3656;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_3647;

   return result;
}

std::complex<double> MSSMRHN_effective_couplings::CpAhSuconjSu(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYu = MODELPARAMETER(TYu);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3665;
   std::complex<double> tmp_3666;
   std::complex<double> tmp_3667;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3668;
      std::complex<double> tmp_3669;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3669 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_3668 += tmp_3669;
      tmp_3667 += (Conj(ZU(gt2,j2))) * tmp_3668;
   }
   tmp_3666 += tmp_3667;
   tmp_3665 += (Conj(Mu)*ZA(gt1,0)) * tmp_3666;
   std::complex<double> tmp_3670;
   std::complex<double> tmp_3671;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3672;
      std::complex<double> tmp_3673;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3673 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3672 += tmp_3673;
      tmp_3671 += (ZU(gt3,j2)) * tmp_3672;
   }
   tmp_3670 += tmp_3671;
   tmp_3665 += (-(Mu*ZA(gt1,0))) * tmp_3670;
   std::complex<double> tmp_3674;
   std::complex<double> tmp_3675;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3676;
      std::complex<double> tmp_3677;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3677 += ZU(gt3,3 + j1)*TYu(j1,j2);
      }
      tmp_3676 += tmp_3677;
      tmp_3675 += (Conj(ZU(gt2,j2))) * tmp_3676;
   }
   tmp_3674 += tmp_3675;
   std::complex<double> tmp_3678;
   std::complex<double> tmp_3679;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3680;
      std::complex<double> tmp_3681;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3681 += Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_3680 += tmp_3681;
      tmp_3679 += (ZU(gt3,j2)) * tmp_3680;
   }
   tmp_3678 += tmp_3679;
   tmp_3674 += (-1) * tmp_3678;
   tmp_3665 += (ZA(gt1,1)) * tmp_3674;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_3665;

   return result;
}

std::complex<double> MSSMRHN_effective_couplings::CpAhSeconjSe(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYe = MODELPARAMETER(TYe);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3682;
   std::complex<double> tmp_3683;
   std::complex<double> tmp_3684;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3685;
      std::complex<double> tmp_3686;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3686 += ZE(gt3,3 + j1)*TYe(j1,j2);
      }
      tmp_3685 += tmp_3686;
      tmp_3684 += (Conj(ZE(gt2,j2))) * tmp_3685;
   }
   tmp_3683 += tmp_3684;
   tmp_3682 += (ZA(gt1,0)) * tmp_3683;
   std::complex<double> tmp_3687;
   std::complex<double> tmp_3688;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3689;
      std::complex<double> tmp_3690;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3690 += Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_3689 += tmp_3690;
      tmp_3688 += (ZE(gt3,j2)) * tmp_3689;
   }
   tmp_3687 += tmp_3688;
   tmp_3682 += (-ZA(gt1,0)) * tmp_3687;
   std::complex<double> tmp_3691;
   std::complex<double> tmp_3692;
   std::complex<double> tmp_3693;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3694;
      std::complex<double> tmp_3695;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3695 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3694 += tmp_3695;
      tmp_3693 += (Conj(ZE(gt2,j2))) * tmp_3694;
   }
   tmp_3692 += tmp_3693;
   tmp_3691 += (Conj(Mu)) * tmp_3692;
   std::complex<double> tmp_3696;
   std::complex<double> tmp_3697;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3698;
      std::complex<double> tmp_3699;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3699 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3698 += tmp_3699;
      tmp_3697 += (ZE(gt3,j2)) * tmp_3698;
   }
   tmp_3696 += tmp_3697;
   tmp_3691 += (-Mu) * tmp_3696;
   tmp_3682 += (ZA(gt1,1)) * tmp_3691;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_3682;

   return result;
}

std::complex<double> MSSMRHN_effective_couplings::CpAhHpmconjHpm(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);

   std::complex<double> result;

   result = std::complex<double>(0,-0.25)*Sqr(g2)*(vu*ZA(gt1,0) + vd*ZA(gt1,1))
      *(ZP(gt2,1)*ZP(gt3,0) - ZP(gt2,0)*ZP(gt3,1));

   return result;
}

std::complex<double> MSSMRHN_effective_couplings::CpAhChabarChaPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g2 = MODELPARAMETER(g2);

   std::complex<double> result;

   result = std::complex<double>(0.,0.7071067811865475)*g2*(Conj(UM(gt2,1))*
      Conj(UP(gt3,0))*ZA(gt1,0) + Conj(UM(gt2,0))*Conj(UP(gt3,1))*ZA(gt1,1));

   return result;
}

std::complex<double> MSSMRHN_effective_couplings::CpAhFebarFePL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Ye = MODELPARAMETER(Ye);

   std::complex<double> result;

   std::complex<double> tmp_3700;
   std::complex<double> tmp_3701;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3702;
      std::complex<double> tmp_3703;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3703 += Conj(ZER(gt3,j1))*Ye(j1,j2);
      }
      tmp_3702 += tmp_3703;
      tmp_3701 += (Conj(ZEL(gt2,j2))) * tmp_3702;
   }
   tmp_3700 += tmp_3701;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,0)) *
      tmp_3700;

   return result;
}

std::complex<double> MSSMRHN_effective_couplings::CpAhFdbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_3704;
   std::complex<double> tmp_3705;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3706;
      std::complex<double> tmp_3707;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3707 += Conj(ZDR(gt3,j1))*Yd(j1,j2);
      }
      tmp_3706 += tmp_3707;
      tmp_3705 += (Conj(ZDL(gt2,j2))) * tmp_3706;
   }
   tmp_3704 += tmp_3705;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,0)) *
      tmp_3704;

   return result;
}

std::complex<double> MSSMRHN_effective_couplings::CpAhFubarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_3708;
   std::complex<double> tmp_3709;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3710;
      std::complex<double> tmp_3711;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3711 += Conj(ZUR(gt3,j1))*Yu(j1,j2);
      }
      tmp_3710 += tmp_3711;
      tmp_3709 += (Conj(ZUL(gt2,j2))) * tmp_3710;
   }
   tmp_3708 += tmp_3709;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,1)) *
      tmp_3708;

   return result;
}

void MSSMRHN_effective_couplings::calculate_eff_CphhVPVP(unsigned gO1)
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

void MSSMRHN_effective_couplings::calculate_eff_CphhVGVG(unsigned gO1)
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

void MSSMRHN_effective_couplings::calculate_eff_CpAhVPVP(unsigned gO1)
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

void MSSMRHN_effective_couplings::calculate_eff_CpAhVGVG(unsigned gO1)
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
