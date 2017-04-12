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

// File generated at Wed 12 Apr 2017 13:01:22

#include "NUTSMSSM_effective_couplings.hpp"

#include "effective_couplings.hpp"
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

standard_model::Standard_model NUTSMSSM_effective_couplings::initialise_SM() const
{
   standard_model::Standard_model sm;

   sm.set_loops(2);
   sm.set_thresholds(2);
   sm.set_physical_input(physical_input);

   sm.initialise_from_input(qedqcd);

   return sm;
}

void NUTSMSSM_effective_couplings::run_SM_strong_coupling_to(standard_model::Standard_model sm, double m)
{
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

   std::complex<double> tmp_3416;
   std::complex<double> tmp_3417;
   std::complex<double> tmp_3418;
   std::complex<double> tmp_3419;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3420;
      std::complex<double> tmp_3421;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3421 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3420 += tmp_3421;
      tmp_3419 += (Conj(ZD(gt2,j2))) * tmp_3420;
   }
   tmp_3418 += tmp_3419;
   tmp_3417 += (Conj(Lambdax)) * tmp_3418;
   std::complex<double> tmp_3422;
   std::complex<double> tmp_3423;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3424;
      std::complex<double> tmp_3425;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3425 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3424 += tmp_3425;
      tmp_3423 += (ZD(gt3,j2)) * tmp_3424;
   }
   tmp_3422 += tmp_3423;
   tmp_3417 += (Lambdax) * tmp_3422;
   tmp_3416 += (6*vu*Conj(ZH(gt1,2))) * tmp_3417;
   std::complex<double> tmp_3426;
   std::complex<double> tmp_3427;
   std::complex<double> tmp_3428;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3428 += Conj(ZD(gt2,j1))*ZD(gt3,j1);
   }
   tmp_3427 += tmp_3428;
   tmp_3426 += (-(vu*(0.6*Sqr(g1) + 3*Sqr(g2)))) * tmp_3427;
   std::complex<double> tmp_3429;
   std::complex<double> tmp_3430;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3430 += Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1);
   }
   tmp_3429 += tmp_3430;
   tmp_3426 += (-1.2*vu*Sqr(g1)) * tmp_3429;
   std::complex<double> tmp_3431;
   std::complex<double> tmp_3432;
   std::complex<double> tmp_3433;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3434;
      std::complex<double> tmp_3435;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3435 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3434 += tmp_3435;
      tmp_3433 += (Conj(ZD(gt2,j2))) * tmp_3434;
   }
   tmp_3432 += tmp_3433;
   tmp_3431 += (vS*Conj(Lambdax)) * tmp_3432;
   std::complex<double> tmp_3436;
   std::complex<double> tmp_3437;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3438;
      std::complex<double> tmp_3439;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3439 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3438 += tmp_3439;
      tmp_3437 += (Conj(ZD(gt2,j2))) * tmp_3438;
   }
   tmp_3436 += tmp_3437;
   tmp_3431 += (1.4142135623730951*Conj(Mu)) * tmp_3436;
   std::complex<double> tmp_3440;
   std::complex<double> tmp_3441;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3442;
      std::complex<double> tmp_3443;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3443 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3442 += tmp_3443;
      tmp_3441 += (ZD(gt3,j2)) * tmp_3442;
   }
   tmp_3440 += tmp_3441;
   tmp_3431 += (vS*Lambdax + 1.4142135623730951*Mu) * tmp_3440;
   tmp_3426 += (6) * tmp_3431;
   tmp_3416 += (Conj(ZH(gt1,1))) * tmp_3426;
   std::complex<double> tmp_3444;
   std::complex<double> tmp_3445;
   std::complex<double> tmp_3446;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3446 += Conj(ZD(gt2,j1))*ZD(gt3,j1);
   }
   tmp_3445 += tmp_3446;
   tmp_3444 += (vd*(0.6*Sqr(g1) + 3*Sqr(g2))) * tmp_3445;
   std::complex<double> tmp_3447;
   std::complex<double> tmp_3448;
   std::complex<double> tmp_3449;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3449 += Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1);
   }
   tmp_3448 += tmp_3449;
   tmp_3447 += (0.6*vd*Sqr(g1)) * tmp_3448;
   std::complex<double> tmp_3450;
   std::complex<double> tmp_3451;
   std::complex<double> tmp_3452;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3453;
      std::complex<double> tmp_3454;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3454 += ZD(gt3,3 + j1)*TYd(j1,j2);
      }
      tmp_3453 += tmp_3454;
      tmp_3452 += (Conj(ZD(gt2,j2))) * tmp_3453;
   }
   tmp_3451 += tmp_3452;
   tmp_3450 += (1.4142135623730951) * tmp_3451;
   std::complex<double> tmp_3455;
   std::complex<double> tmp_3456;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3457;
      std::complex<double> tmp_3458;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3458 += Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_3457 += tmp_3458;
      tmp_3456 += (ZD(gt3,j2)) * tmp_3457;
   }
   tmp_3455 += tmp_3456;
   tmp_3450 += (1.4142135623730951) * tmp_3455;
   std::complex<double> tmp_3459;
   std::complex<double> tmp_3460;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3461;
      std::complex<double> tmp_3462;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3463;
         std::complex<double> tmp_3464;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3464 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_3463 += tmp_3464;
         tmp_3462 += (ZD(gt3,3 + j2)) * tmp_3463;
      }
      tmp_3461 += tmp_3462;
      tmp_3460 += (Conj(ZD(gt2,3 + j3))) * tmp_3461;
   }
   tmp_3459 += tmp_3460;
   std::complex<double> tmp_3465;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3466;
      std::complex<double> tmp_3467;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3468;
         std::complex<double> tmp_3469;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3469 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_3468 += tmp_3469;
         tmp_3467 += (Conj(ZD(gt2,j2))) * tmp_3468;
      }
      tmp_3466 += tmp_3467;
      tmp_3465 += (ZD(gt3,j3)) * tmp_3466;
   }
   tmp_3459 += tmp_3465;
   tmp_3450 += (2*vd) * tmp_3459;
   tmp_3447 += (-3) * tmp_3450;
   tmp_3444 += (2) * tmp_3447;
   tmp_3416 += (Conj(ZH(gt1,0))) * tmp_3444;
   result += (0.08333333333333333) * tmp_3416;

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

   std::complex<double> tmp_3470;
   std::complex<double> tmp_3471;
   std::complex<double> tmp_3472;
   std::complex<double> tmp_3473;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3474;
      std::complex<double> tmp_3475;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3475 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_3474 += tmp_3475;
      tmp_3473 += (Conj(ZU(gt2,j2))) * tmp_3474;
   }
   tmp_3472 += tmp_3473;
   tmp_3471 += (Conj(Lambdax)) * tmp_3472;
   std::complex<double> tmp_3476;
   std::complex<double> tmp_3477;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3478;
      std::complex<double> tmp_3479;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3479 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3478 += tmp_3479;
      tmp_3477 += (ZU(gt3,j2)) * tmp_3478;
   }
   tmp_3476 += tmp_3477;
   tmp_3471 += (Lambdax) * tmp_3476;
   tmp_3470 += (6*vd*Conj(ZH(gt1,2))) * tmp_3471;
   std::complex<double> tmp_3480;
   std::complex<double> tmp_3481;
   std::complex<double> tmp_3482;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3482 += Conj(ZU(gt2,j1))*ZU(gt3,j1);
   }
   tmp_3481 += tmp_3482;
   tmp_3480 += (vd*(0.6*Sqr(g1) - 3*Sqr(g2))) * tmp_3481;
   std::complex<double> tmp_3483;
   std::complex<double> tmp_3484;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3484 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_3483 += tmp_3484;
   tmp_3480 += (-2.4*vd*Sqr(g1)) * tmp_3483;
   std::complex<double> tmp_3485;
   std::complex<double> tmp_3486;
   std::complex<double> tmp_3487;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3488;
      std::complex<double> tmp_3489;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3489 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_3488 += tmp_3489;
      tmp_3487 += (Conj(ZU(gt2,j2))) * tmp_3488;
   }
   tmp_3486 += tmp_3487;
   tmp_3485 += (vS*Conj(Lambdax)) * tmp_3486;
   std::complex<double> tmp_3490;
   std::complex<double> tmp_3491;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3492;
      std::complex<double> tmp_3493;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3493 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_3492 += tmp_3493;
      tmp_3491 += (Conj(ZU(gt2,j2))) * tmp_3492;
   }
   tmp_3490 += tmp_3491;
   tmp_3485 += (1.4142135623730951*Conj(Mu)) * tmp_3490;
   std::complex<double> tmp_3494;
   std::complex<double> tmp_3495;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3496;
      std::complex<double> tmp_3497;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3497 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3496 += tmp_3497;
      tmp_3495 += (ZU(gt3,j2)) * tmp_3496;
   }
   tmp_3494 += tmp_3495;
   tmp_3485 += (vS*Lambdax + 1.4142135623730951*Mu) * tmp_3494;
   tmp_3480 += (6) * tmp_3485;
   tmp_3470 += (Conj(ZH(gt1,0))) * tmp_3480;
   std::complex<double> tmp_3498;
   std::complex<double> tmp_3499;
   std::complex<double> tmp_3500;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3500 += Conj(ZU(gt2,j1))*ZU(gt3,j1);
   }
   tmp_3499 += tmp_3500;
   tmp_3498 += (vu*(0.6*Sqr(g1) - 3*Sqr(g2))) * tmp_3499;
   std::complex<double> tmp_3501;
   std::complex<double> tmp_3502;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3502 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_3501 += tmp_3502;
   tmp_3498 += (-2.4*vu*Sqr(g1)) * tmp_3501;
   std::complex<double> tmp_3503;
   std::complex<double> tmp_3504;
   std::complex<double> tmp_3505;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3506;
      std::complex<double> tmp_3507;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3507 += ZU(gt3,3 + j1)*TYu(j1,j2);
      }
      tmp_3506 += tmp_3507;
      tmp_3505 += (Conj(ZU(gt2,j2))) * tmp_3506;
   }
   tmp_3504 += tmp_3505;
   tmp_3503 += (1.4142135623730951) * tmp_3504;
   std::complex<double> tmp_3508;
   std::complex<double> tmp_3509;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3510;
      std::complex<double> tmp_3511;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3511 += Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_3510 += tmp_3511;
      tmp_3509 += (ZU(gt3,j2)) * tmp_3510;
   }
   tmp_3508 += tmp_3509;
   tmp_3503 += (1.4142135623730951) * tmp_3508;
   std::complex<double> tmp_3512;
   std::complex<double> tmp_3513;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3514;
      std::complex<double> tmp_3515;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3516;
         std::complex<double> tmp_3517;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3517 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_3516 += tmp_3517;
         tmp_3515 += (ZU(gt3,3 + j2)) * tmp_3516;
      }
      tmp_3514 += tmp_3515;
      tmp_3513 += (Conj(ZU(gt2,3 + j3))) * tmp_3514;
   }
   tmp_3512 += tmp_3513;
   std::complex<double> tmp_3518;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3519;
      std::complex<double> tmp_3520;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3521;
         std::complex<double> tmp_3522;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3522 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_3521 += tmp_3522;
         tmp_3520 += (Conj(ZU(gt2,j2))) * tmp_3521;
      }
      tmp_3519 += tmp_3520;
      tmp_3518 += (ZU(gt3,j3)) * tmp_3519;
   }
   tmp_3512 += tmp_3518;
   tmp_3503 += (2*vu) * tmp_3512;
   tmp_3498 += (6) * tmp_3503;
   tmp_3470 += (-Conj(ZH(gt1,1))) * tmp_3498;
   result += (0.08333333333333333) * tmp_3470;

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

   std::complex<double> tmp_3523;
   std::complex<double> tmp_3524;
   std::complex<double> tmp_3525;
   std::complex<double> tmp_3526;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3527;
      std::complex<double> tmp_3528;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3528 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3527 += tmp_3528;
      tmp_3526 += (Conj(ZE(gt2,j2))) * tmp_3527;
   }
   tmp_3525 += tmp_3526;
   tmp_3524 += (Conj(Lambdax)) * tmp_3525;
   std::complex<double> tmp_3529;
   std::complex<double> tmp_3530;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3531;
      std::complex<double> tmp_3532;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3532 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3531 += tmp_3532;
      tmp_3530 += (ZE(gt3,j2)) * tmp_3531;
   }
   tmp_3529 += tmp_3530;
   tmp_3524 += (Lambdax) * tmp_3529;
   tmp_3523 += (2*vu*Conj(ZH(gt1,2))) * tmp_3524;
   std::complex<double> tmp_3533;
   std::complex<double> tmp_3534;
   std::complex<double> tmp_3535;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3535 += Conj(ZE(gt2,j1))*ZE(gt3,j1);
   }
   tmp_3534 += tmp_3535;
   tmp_3533 += (vu*(0.6*Sqr(g1) - Sqr(g2))) * tmp_3534;
   std::complex<double> tmp_3536;
   std::complex<double> tmp_3537;
   std::complex<double> tmp_3538;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3538 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_3537 += tmp_3538;
   tmp_3536 += (-0.6*vu*Sqr(g1)) * tmp_3537;
   std::complex<double> tmp_3539;
   std::complex<double> tmp_3540;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3541;
      std::complex<double> tmp_3542;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3542 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3541 += tmp_3542;
      tmp_3540 += (Conj(ZE(gt2,j2))) * tmp_3541;
   }
   tmp_3539 += tmp_3540;
   tmp_3536 += (vS*Conj(Lambdax)) * tmp_3539;
   std::complex<double> tmp_3543;
   std::complex<double> tmp_3544;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3545;
      std::complex<double> tmp_3546;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3546 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3545 += tmp_3546;
      tmp_3544 += (Conj(ZE(gt2,j2))) * tmp_3545;
   }
   tmp_3543 += tmp_3544;
   tmp_3536 += (1.4142135623730951*Conj(Mu)) * tmp_3543;
   std::complex<double> tmp_3547;
   std::complex<double> tmp_3548;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3549;
      std::complex<double> tmp_3550;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3550 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3549 += tmp_3550;
      tmp_3548 += (ZE(gt3,j2)) * tmp_3549;
   }
   tmp_3547 += tmp_3548;
   tmp_3536 += (vS*Lambdax) * tmp_3547;
   std::complex<double> tmp_3551;
   std::complex<double> tmp_3552;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3553;
      std::complex<double> tmp_3554;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3554 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3553 += tmp_3554;
      tmp_3552 += (ZE(gt3,j2)) * tmp_3553;
   }
   tmp_3551 += tmp_3552;
   tmp_3536 += (1.4142135623730951*Mu) * tmp_3551;
   tmp_3533 += (2) * tmp_3536;
   tmp_3523 += (Conj(ZH(gt1,1))) * tmp_3533;
   std::complex<double> tmp_3555;
   std::complex<double> tmp_3556;
   std::complex<double> tmp_3557;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3557 += Conj(ZE(gt2,j1))*ZE(gt3,j1);
   }
   tmp_3556 += tmp_3557;
   tmp_3555 += (vd*(0.6*Sqr(g1) - Sqr(g2))) * tmp_3556;
   std::complex<double> tmp_3558;
   std::complex<double> tmp_3559;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3559 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_3558 += tmp_3559;
   tmp_3555 += (-1.2*vd*Sqr(g1)) * tmp_3558;
   std::complex<double> tmp_3560;
   std::complex<double> tmp_3561;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3562;
      std::complex<double> tmp_3563;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3563 += ZE(gt3,3 + j1)*TYe(j1,j2);
      }
      tmp_3562 += tmp_3563;
      tmp_3561 += (Conj(ZE(gt2,j2))) * tmp_3562;
   }
   tmp_3560 += tmp_3561;
   tmp_3555 += (2.8284271247461903) * tmp_3560;
   std::complex<double> tmp_3564;
   std::complex<double> tmp_3565;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3566;
      std::complex<double> tmp_3567;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3567 += Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_3566 += tmp_3567;
      tmp_3565 += (ZE(gt3,j2)) * tmp_3566;
   }
   tmp_3564 += tmp_3565;
   tmp_3555 += (2.8284271247461903) * tmp_3564;
   std::complex<double> tmp_3568;
   std::complex<double> tmp_3569;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3570;
      std::complex<double> tmp_3571;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3572;
         std::complex<double> tmp_3573;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3573 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_3572 += tmp_3573;
         tmp_3571 += (ZE(gt3,3 + j2)) * tmp_3572;
      }
      tmp_3570 += tmp_3571;
      tmp_3569 += (Conj(ZE(gt2,3 + j3))) * tmp_3570;
   }
   tmp_3568 += tmp_3569;
   tmp_3555 += (4*vd) * tmp_3568;
   std::complex<double> tmp_3574;
   std::complex<double> tmp_3575;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3576;
      std::complex<double> tmp_3577;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3578;
         std::complex<double> tmp_3579;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3579 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_3578 += tmp_3579;
         tmp_3577 += (Conj(ZE(gt2,j2))) * tmp_3578;
      }
      tmp_3576 += tmp_3577;
      tmp_3575 += (ZE(gt3,j3)) * tmp_3576;
   }
   tmp_3574 += tmp_3575;
   tmp_3555 += (4*vd) * tmp_3574;
   tmp_3523 += (-Conj(ZH(gt1,0))) * tmp_3555;
   result += (0.25) * tmp_3523;

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

   std::complex<double> tmp_3580;
   std::complex<double> tmp_3581;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3582;
      std::complex<double> tmp_3583;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3583 += Conj(ZER(gt3,j1))*Ye(j1,j2);
      }
      tmp_3582 += tmp_3583;
      tmp_3581 += (Conj(ZEL(gt1,j2))) * tmp_3582;
   }
   tmp_3580 += tmp_3581;
   result += (-0.7071067811865475*Conj(ZH(gt2,0))) * tmp_3580;

   return result;
}

std::complex<double> NUTSMSSM_effective_couplings::CpFdhhbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_3584;
   std::complex<double> tmp_3585;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3586;
      std::complex<double> tmp_3587;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3587 += Conj(ZDR(gt3,j1))*Yd(j1,j2);
      }
      tmp_3586 += tmp_3587;
      tmp_3585 += (Conj(ZDL(gt1,j2))) * tmp_3586;
   }
   tmp_3584 += tmp_3585;
   result += (-0.7071067811865475*Conj(ZH(gt2,0))) * tmp_3584;

   return result;
}

std::complex<double> NUTSMSSM_effective_couplings::CpFuhhbarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_3588;
   std::complex<double> tmp_3589;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3590;
      std::complex<double> tmp_3591;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3591 += Conj(ZUR(gt3,j1))*Yu(j1,j2);
      }
      tmp_3590 += tmp_3591;
      tmp_3589 += (Conj(ZUL(gt1,j2))) * tmp_3590;
   }
   tmp_3588 += tmp_3589;
   result += (-0.7071067811865475*Conj(ZH(gt2,1))) * tmp_3588;

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

   std::complex<double> tmp_3592;
   std::complex<double> tmp_3593;
   std::complex<double> tmp_3594;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3595;
      std::complex<double> tmp_3596;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3596 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3595 += tmp_3596;
      tmp_3594 += (Conj(ZD(gt2,j2))) * tmp_3595;
   }
   tmp_3593 += tmp_3594;
   tmp_3592 += (1.4142135623730951*Conj(Mu)*Conj(ZA(gt1,1))) * tmp_3593;
   std::complex<double> tmp_3597;
   std::complex<double> tmp_3598;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3599;
      std::complex<double> tmp_3600;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3600 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3599 += tmp_3600;
      tmp_3598 += (Conj(ZD(gt2,j2))) * tmp_3599;
   }
   tmp_3597 += tmp_3598;
   tmp_3592 += (Conj(Lambdax)*(vS*Conj(ZA(gt1,1)) + vu*Conj(ZA(gt1,2)))) *
      tmp_3597;
   std::complex<double> tmp_3601;
   std::complex<double> tmp_3602;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3603;
      std::complex<double> tmp_3604;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3604 += ZD(gt3,3 + j1)*TYd(j1,j2);
      }
      tmp_3603 += tmp_3604;
      tmp_3602 += (Conj(ZD(gt2,j2))) * tmp_3603;
   }
   tmp_3601 += tmp_3602;
   tmp_3592 += (1.4142135623730951*Conj(ZA(gt1,0))) * tmp_3601;
   std::complex<double> tmp_3605;
   std::complex<double> tmp_3606;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3607;
      std::complex<double> tmp_3608;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3608 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3607 += tmp_3608;
      tmp_3606 += (ZD(gt3,j2)) * tmp_3607;
   }
   tmp_3605 += tmp_3606;
   tmp_3592 += (-(vS*Conj(ZA(gt1,1))*Lambdax)) * tmp_3605;
   std::complex<double> tmp_3609;
   std::complex<double> tmp_3610;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3611;
      std::complex<double> tmp_3612;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3612 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3611 += tmp_3612;
      tmp_3610 += (ZD(gt3,j2)) * tmp_3611;
   }
   tmp_3609 += tmp_3610;
   tmp_3592 += (-1.4142135623730951*Conj(ZA(gt1,1))*Mu) * tmp_3609;
   std::complex<double> tmp_3613;
   std::complex<double> tmp_3614;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3615;
      std::complex<double> tmp_3616;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3616 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3615 += tmp_3616;
      tmp_3614 += (ZD(gt3,j2)) * tmp_3615;
   }
   tmp_3613 += tmp_3614;
   tmp_3592 += (-(vu*Conj(ZA(gt1,2))*Lambdax)) * tmp_3613;
   std::complex<double> tmp_3617;
   std::complex<double> tmp_3618;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3619;
      std::complex<double> tmp_3620;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3620 += Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_3619 += tmp_3620;
      tmp_3618 += (ZD(gt3,j2)) * tmp_3619;
   }
   tmp_3617 += tmp_3618;
   tmp_3592 += (-1.4142135623730951*Conj(ZA(gt1,0))) * tmp_3617;
   result += (std::complex<double>(0,-0.5)) * tmp_3592;

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

   std::complex<double> tmp_3621;
   std::complex<double> tmp_3622;
   std::complex<double> tmp_3623;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3624;
      std::complex<double> tmp_3625;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3625 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_3624 += tmp_3625;
      tmp_3623 += (Conj(ZU(gt2,j2))) * tmp_3624;
   }
   tmp_3622 += tmp_3623;
   tmp_3621 += (1.4142135623730951*Conj(Mu)*Conj(ZA(gt1,0))) * tmp_3622;
   std::complex<double> tmp_3626;
   std::complex<double> tmp_3627;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3628;
      std::complex<double> tmp_3629;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3629 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_3628 += tmp_3629;
      tmp_3627 += (Conj(ZU(gt2,j2))) * tmp_3628;
   }
   tmp_3626 += tmp_3627;
   tmp_3621 += (Conj(Lambdax)*(vS*Conj(ZA(gt1,0)) + vd*Conj(ZA(gt1,2)))) *
      tmp_3626;
   std::complex<double> tmp_3630;
   std::complex<double> tmp_3631;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3632;
      std::complex<double> tmp_3633;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3633 += ZU(gt3,3 + j1)*TYu(j1,j2);
      }
      tmp_3632 += tmp_3633;
      tmp_3631 += (Conj(ZU(gt2,j2))) * tmp_3632;
   }
   tmp_3630 += tmp_3631;
   tmp_3621 += (1.4142135623730951*Conj(ZA(gt1,1))) * tmp_3630;
   std::complex<double> tmp_3634;
   std::complex<double> tmp_3635;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3636;
      std::complex<double> tmp_3637;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3637 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3636 += tmp_3637;
      tmp_3635 += (ZU(gt3,j2)) * tmp_3636;
   }
   tmp_3634 += tmp_3635;
   tmp_3621 += (-(vS*Conj(ZA(gt1,0))*Lambdax)) * tmp_3634;
   std::complex<double> tmp_3638;
   std::complex<double> tmp_3639;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3640;
      std::complex<double> tmp_3641;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3641 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3640 += tmp_3641;
      tmp_3639 += (ZU(gt3,j2)) * tmp_3640;
   }
   tmp_3638 += tmp_3639;
   tmp_3621 += (-1.4142135623730951*Conj(ZA(gt1,0))*Mu) * tmp_3638;
   std::complex<double> tmp_3642;
   std::complex<double> tmp_3643;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3644;
      std::complex<double> tmp_3645;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3645 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3644 += tmp_3645;
      tmp_3643 += (ZU(gt3,j2)) * tmp_3644;
   }
   tmp_3642 += tmp_3643;
   tmp_3621 += (-(vd*Conj(ZA(gt1,2))*Lambdax)) * tmp_3642;
   std::complex<double> tmp_3646;
   std::complex<double> tmp_3647;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3648;
      std::complex<double> tmp_3649;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3649 += Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_3648 += tmp_3649;
      tmp_3647 += (ZU(gt3,j2)) * tmp_3648;
   }
   tmp_3646 += tmp_3647;
   tmp_3621 += (-1.4142135623730951*Conj(ZA(gt1,1))) * tmp_3646;
   result += (std::complex<double>(0,-0.5)) * tmp_3621;

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

   std::complex<double> tmp_3650;
   std::complex<double> tmp_3651;
   std::complex<double> tmp_3652;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3653;
      std::complex<double> tmp_3654;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3654 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3653 += tmp_3654;
      tmp_3652 += (Conj(ZE(gt2,j2))) * tmp_3653;
   }
   tmp_3651 += tmp_3652;
   tmp_3650 += (1.4142135623730951*Conj(Mu)*Conj(ZA(gt1,1))) * tmp_3651;
   std::complex<double> tmp_3655;
   std::complex<double> tmp_3656;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3657;
      std::complex<double> tmp_3658;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3658 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3657 += tmp_3658;
      tmp_3656 += (Conj(ZE(gt2,j2))) * tmp_3657;
   }
   tmp_3655 += tmp_3656;
   tmp_3650 += (Conj(Lambdax)*(vS*Conj(ZA(gt1,1)) + vu*Conj(ZA(gt1,2)))) *
      tmp_3655;
   std::complex<double> tmp_3659;
   std::complex<double> tmp_3660;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3661;
      std::complex<double> tmp_3662;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3662 += ZE(gt3,3 + j1)*TYe(j1,j2);
      }
      tmp_3661 += tmp_3662;
      tmp_3660 += (Conj(ZE(gt2,j2))) * tmp_3661;
   }
   tmp_3659 += tmp_3660;
   tmp_3650 += (1.4142135623730951*Conj(ZA(gt1,0))) * tmp_3659;
   std::complex<double> tmp_3663;
   std::complex<double> tmp_3664;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3665;
      std::complex<double> tmp_3666;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3666 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3665 += tmp_3666;
      tmp_3664 += (ZE(gt3,j2)) * tmp_3665;
   }
   tmp_3663 += tmp_3664;
   tmp_3650 += (-(vS*Conj(ZA(gt1,1))*Lambdax)) * tmp_3663;
   std::complex<double> tmp_3667;
   std::complex<double> tmp_3668;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3669;
      std::complex<double> tmp_3670;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3670 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3669 += tmp_3670;
      tmp_3668 += (ZE(gt3,j2)) * tmp_3669;
   }
   tmp_3667 += tmp_3668;
   tmp_3650 += (-1.4142135623730951*Conj(ZA(gt1,1))*Mu) * tmp_3667;
   std::complex<double> tmp_3671;
   std::complex<double> tmp_3672;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3673;
      std::complex<double> tmp_3674;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3674 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3673 += tmp_3674;
      tmp_3672 += (ZE(gt3,j2)) * tmp_3673;
   }
   tmp_3671 += tmp_3672;
   tmp_3650 += (-(vu*Conj(ZA(gt1,2))*Lambdax)) * tmp_3671;
   std::complex<double> tmp_3675;
   std::complex<double> tmp_3676;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3677;
      std::complex<double> tmp_3678;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3678 += Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_3677 += tmp_3678;
      tmp_3676 += (ZE(gt3,j2)) * tmp_3677;
   }
   tmp_3675 += tmp_3676;
   tmp_3650 += (-1.4142135623730951*Conj(ZA(gt1,0))) * tmp_3675;
   result += (std::complex<double>(0,-0.5)) * tmp_3650;

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

   std::complex<double> tmp_3679;
   std::complex<double> tmp_3680;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3681;
      std::complex<double> tmp_3682;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3682 += Conj(ZER(gt3,j1))*Ye(j1,j2);
      }
      tmp_3681 += tmp_3682;
      tmp_3680 += (Conj(ZEL(gt2,j2))) * tmp_3681;
   }
   tmp_3679 += tmp_3680;
   result += (std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gt1,0))) *
      tmp_3679;

   return result;
}

std::complex<double> NUTSMSSM_effective_couplings::CpAhFdbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_3683;
   std::complex<double> tmp_3684;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3685;
      std::complex<double> tmp_3686;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3686 += Conj(ZDR(gt3,j1))*Yd(j1,j2);
      }
      tmp_3685 += tmp_3686;
      tmp_3684 += (Conj(ZDL(gt2,j2))) * tmp_3685;
   }
   tmp_3683 += tmp_3684;
   result += (std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gt1,0))) *
      tmp_3683;

   return result;
}

std::complex<double> NUTSMSSM_effective_couplings::CpAhFubarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_3687;
   std::complex<double> tmp_3688;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3689;
      std::complex<double> tmp_3690;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3690 += Conj(ZUR(gt3,j1))*Yu(j1,j2);
      }
      tmp_3689 += tmp_3690;
      tmp_3688 += (Conj(ZUL(gt2,j2))) * tmp_3689;
   }
   tmp_3687 += tmp_3688;
   result += (std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gt1,1))) *
      tmp_3687;

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
