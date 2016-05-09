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

// File generated at Mon 9 May 2016 13:00:21

#include "UMSSM_effective_couplings.hpp"

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

UMSSM_effective_couplings::~UMSSM_effective_couplings()
{
}

void UMSSM_effective_couplings::calculate_effective_couplings()
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
   for (unsigned gO1 = 2; gO1 < 3; ++gO1) {
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

void UMSSM_effective_couplings::run_SM_strong_coupling_to(double m)
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

   std::complex<double> tmp_4727;
   std::complex<double> tmp_4728;
   std::complex<double> tmp_4729;
   std::complex<double> tmp_4730;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4730 += Conj(ZD(gt2,j1))*ZD(gt3,j1);
   }
   tmp_4729 += tmp_4730;
   tmp_4728 += (-2*Qq*Qs*vS*Sqr(gp)) * tmp_4729;
   std::complex<double> tmp_4731;
   std::complex<double> tmp_4732;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4732 += Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1);
   }
   tmp_4731 += tmp_4732;
   tmp_4728 += (-2*Qd*Qs*vS*Sqr(gp)) * tmp_4731;
   std::complex<double> tmp_4733;
   std::complex<double> tmp_4734;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4735;
      std::complex<double> tmp_4736;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4736 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_4735 += tmp_4736;
      tmp_4734 += (Conj(ZD(gt2,j2))) * tmp_4735;
   }
   tmp_4733 += tmp_4734;
   tmp_4728 += (vu*Conj(Lambdax)) * tmp_4733;
   std::complex<double> tmp_4737;
   std::complex<double> tmp_4738;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4739;
      std::complex<double> tmp_4740;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4740 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_4739 += tmp_4740;
      tmp_4738 += (ZD(gt3,j2)) * tmp_4739;
   }
   tmp_4737 += tmp_4738;
   tmp_4728 += (vu*Lambdax) * tmp_4737;
   tmp_4727 += (6*Conj(ZH(gt1,2))) * tmp_4728;
   std::complex<double> tmp_4741;
   std::complex<double> tmp_4742;
   std::complex<double> tmp_4743;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4743 += Conj(ZD(gt2,j1))*ZD(gt3,j1);
   }
   tmp_4742 += tmp_4743;
   tmp_4741 += (vu*(0.6*Sqr(g1) + 3*(Sqr(g2) + 4*QHu*Qq*Sqr(gp)))) * tmp_4742;
   std::complex<double> tmp_4744;
   std::complex<double> tmp_4745;
   std::complex<double> tmp_4746;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4746 += Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1);
   }
   tmp_4745 += tmp_4746;
   tmp_4744 += (vu*(0.6*Sqr(g1) + 6*Qd*QHu*Sqr(gp))) * tmp_4745;
   std::complex<double> tmp_4747;
   std::complex<double> tmp_4748;
   std::complex<double> tmp_4749;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4750;
      std::complex<double> tmp_4751;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4751 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_4750 += tmp_4751;
      tmp_4749 += (Conj(ZD(gt2,j2))) * tmp_4750;
   }
   tmp_4748 += tmp_4749;
   tmp_4747 += (Conj(Lambdax)) * tmp_4748;
   std::complex<double> tmp_4752;
   std::complex<double> tmp_4753;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4754;
      std::complex<double> tmp_4755;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4755 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_4754 += tmp_4755;
      tmp_4753 += (ZD(gt3,j2)) * tmp_4754;
   }
   tmp_4752 += tmp_4753;
   tmp_4747 += (Lambdax) * tmp_4752;
   tmp_4744 += (-3*vS) * tmp_4747;
   tmp_4741 += (2) * tmp_4744;
   tmp_4727 += (-Conj(ZH(gt1,1))) * tmp_4741;
   std::complex<double> tmp_4756;
   std::complex<double> tmp_4757;
   std::complex<double> tmp_4758;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4758 += Conj(ZD(gt2,j1))*ZD(gt3,j1);
   }
   tmp_4757 += tmp_4758;
   tmp_4756 += (vd*(0.6*Sqr(g1) + 3*(Sqr(g2) - 4*QHd*Qq*Sqr(gp)))) * tmp_4757;
   std::complex<double> tmp_4759;
   std::complex<double> tmp_4760;
   std::complex<double> tmp_4761;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4761 += Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1);
   }
   tmp_4760 += tmp_4761;
   tmp_4759 += (vd*(0.6*Sqr(g1) - 6*Qd*QHd*Sqr(gp))) * tmp_4760;
   std::complex<double> tmp_4762;
   std::complex<double> tmp_4763;
   std::complex<double> tmp_4764;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4765;
      std::complex<double> tmp_4766;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4766 += ZD(gt3,3 + j1)*TYd(j1,j2);
      }
      tmp_4765 += tmp_4766;
      tmp_4764 += (Conj(ZD(gt2,j2))) * tmp_4765;
   }
   tmp_4763 += tmp_4764;
   tmp_4762 += (1.4142135623730951) * tmp_4763;
   std::complex<double> tmp_4767;
   std::complex<double> tmp_4768;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4769;
      std::complex<double> tmp_4770;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4770 += Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_4769 += tmp_4770;
      tmp_4768 += (ZD(gt3,j2)) * tmp_4769;
   }
   tmp_4767 += tmp_4768;
   tmp_4762 += (1.4142135623730951) * tmp_4767;
   std::complex<double> tmp_4771;
   std::complex<double> tmp_4772;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_4773;
      std::complex<double> tmp_4774;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_4775;
         std::complex<double> tmp_4776;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_4776 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_4775 += tmp_4776;
         tmp_4774 += (ZD(gt3,3 + j2)) * tmp_4775;
      }
      tmp_4773 += tmp_4774;
      tmp_4772 += (Conj(ZD(gt2,3 + j3))) * tmp_4773;
   }
   tmp_4771 += tmp_4772;
   std::complex<double> tmp_4777;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_4778;
      std::complex<double> tmp_4779;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_4780;
         std::complex<double> tmp_4781;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_4781 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_4780 += tmp_4781;
         tmp_4779 += (Conj(ZD(gt2,j2))) * tmp_4780;
      }
      tmp_4778 += tmp_4779;
      tmp_4777 += (ZD(gt3,j3)) * tmp_4778;
   }
   tmp_4771 += tmp_4777;
   tmp_4762 += (2*vd) * tmp_4771;
   tmp_4759 += (-3) * tmp_4762;
   tmp_4756 += (2) * tmp_4759;
   tmp_4727 += (Conj(ZH(gt1,0))) * tmp_4756;
   result += (0.08333333333333333) * tmp_4727;

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

   std::complex<double> tmp_4782;
   std::complex<double> tmp_4783;
   std::complex<double> tmp_4784;
   std::complex<double> tmp_4785;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4785 += Conj(ZU(gt2,j1))*ZU(gt3,j1);
   }
   tmp_4784 += tmp_4785;
   tmp_4783 += (-2*Qq*Qs*vS*Sqr(gp)) * tmp_4784;
   std::complex<double> tmp_4786;
   std::complex<double> tmp_4787;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4787 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_4786 += tmp_4787;
   tmp_4783 += (-2*Qs*Qu*vS*Sqr(gp)) * tmp_4786;
   std::complex<double> tmp_4788;
   std::complex<double> tmp_4789;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4790;
      std::complex<double> tmp_4791;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4791 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_4790 += tmp_4791;
      tmp_4789 += (Conj(ZU(gt2,j2))) * tmp_4790;
   }
   tmp_4788 += tmp_4789;
   tmp_4783 += (vd*Conj(Lambdax)) * tmp_4788;
   std::complex<double> tmp_4792;
   std::complex<double> tmp_4793;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4794;
      std::complex<double> tmp_4795;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4795 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_4794 += tmp_4795;
      tmp_4793 += (ZU(gt3,j2)) * tmp_4794;
   }
   tmp_4792 += tmp_4793;
   tmp_4783 += (vd*Lambdax) * tmp_4792;
   tmp_4782 += (6*Conj(ZH(gt1,2))) * tmp_4783;
   std::complex<double> tmp_4796;
   std::complex<double> tmp_4797;
   std::complex<double> tmp_4798;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4798 += Conj(ZU(gt2,j1))*ZU(gt3,j1);
   }
   tmp_4797 += tmp_4798;
   tmp_4796 += (vd*(0.6*Sqr(g1) - 3*(Sqr(g2) + 4*QHd*Qq*Sqr(gp)))) * tmp_4797;
   std::complex<double> tmp_4799;
   std::complex<double> tmp_4800;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4800 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_4799 += tmp_4800;
   tmp_4796 += (-4*vd*(0.6*Sqr(g1) + 3*QHd*Qu*Sqr(gp))) * tmp_4799;
   std::complex<double> tmp_4801;
   std::complex<double> tmp_4802;
   std::complex<double> tmp_4803;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4804;
      std::complex<double> tmp_4805;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4805 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_4804 += tmp_4805;
      tmp_4803 += (Conj(ZU(gt2,j2))) * tmp_4804;
   }
   tmp_4802 += tmp_4803;
   tmp_4801 += (Conj(Lambdax)) * tmp_4802;
   std::complex<double> tmp_4806;
   std::complex<double> tmp_4807;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4808;
      std::complex<double> tmp_4809;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4809 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_4808 += tmp_4809;
      tmp_4807 += (ZU(gt3,j2)) * tmp_4808;
   }
   tmp_4806 += tmp_4807;
   tmp_4801 += (Lambdax) * tmp_4806;
   tmp_4796 += (6*vS) * tmp_4801;
   tmp_4782 += (Conj(ZH(gt1,0))) * tmp_4796;
   std::complex<double> tmp_4810;
   std::complex<double> tmp_4811;
   std::complex<double> tmp_4812;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4812 += Conj(ZU(gt2,j1))*ZU(gt3,j1);
   }
   tmp_4811 += tmp_4812;
   tmp_4810 += (vu*(0.6*Sqr(g1) - 3*Sqr(g2) + 12*QHu*Qq*Sqr(gp))) * tmp_4811;
   std::complex<double> tmp_4813;
   std::complex<double> tmp_4814;
   std::complex<double> tmp_4815;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4815 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_4814 += tmp_4815;
   tmp_4813 += (-2*vu*(0.6*Sqr(g1) - 3*QHu*Qu*Sqr(gp))) * tmp_4814;
   std::complex<double> tmp_4816;
   std::complex<double> tmp_4817;
   std::complex<double> tmp_4818;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4819;
      std::complex<double> tmp_4820;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4820 += ZU(gt3,3 + j1)*TYu(j1,j2);
      }
      tmp_4819 += tmp_4820;
      tmp_4818 += (Conj(ZU(gt2,j2))) * tmp_4819;
   }
   tmp_4817 += tmp_4818;
   tmp_4816 += (1.4142135623730951) * tmp_4817;
   std::complex<double> tmp_4821;
   std::complex<double> tmp_4822;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4823;
      std::complex<double> tmp_4824;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4824 += Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_4823 += tmp_4824;
      tmp_4822 += (ZU(gt3,j2)) * tmp_4823;
   }
   tmp_4821 += tmp_4822;
   tmp_4816 += (1.4142135623730951) * tmp_4821;
   std::complex<double> tmp_4825;
   std::complex<double> tmp_4826;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_4827;
      std::complex<double> tmp_4828;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_4829;
         std::complex<double> tmp_4830;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_4830 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_4829 += tmp_4830;
         tmp_4828 += (ZU(gt3,3 + j2)) * tmp_4829;
      }
      tmp_4827 += tmp_4828;
      tmp_4826 += (Conj(ZU(gt2,3 + j3))) * tmp_4827;
   }
   tmp_4825 += tmp_4826;
   std::complex<double> tmp_4831;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_4832;
      std::complex<double> tmp_4833;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_4834;
         std::complex<double> tmp_4835;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_4835 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_4834 += tmp_4835;
         tmp_4833 += (Conj(ZU(gt2,j2))) * tmp_4834;
      }
      tmp_4832 += tmp_4833;
      tmp_4831 += (ZU(gt3,j3)) * tmp_4832;
   }
   tmp_4825 += tmp_4831;
   tmp_4816 += (2*vu) * tmp_4825;
   tmp_4813 += (3) * tmp_4816;
   tmp_4810 += (2) * tmp_4813;
   tmp_4782 += (-Conj(ZH(gt1,1))) * tmp_4810;
   result += (0.08333333333333333) * tmp_4782;

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

   std::complex<double> tmp_4836;
   std::complex<double> tmp_4837;
   std::complex<double> tmp_4838;
   std::complex<double> tmp_4839;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4839 += Conj(ZE(gt2,j1))*ZE(gt3,j1);
   }
   tmp_4838 += tmp_4839;
   tmp_4837 += (-2*Ql*Qs*vS*Sqr(gp)) * tmp_4838;
   std::complex<double> tmp_4840;
   std::complex<double> tmp_4841;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4841 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_4840 += tmp_4841;
   tmp_4837 += (-2*Qe*Qs*vS*Sqr(gp)) * tmp_4840;
   std::complex<double> tmp_4842;
   std::complex<double> tmp_4843;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4844;
      std::complex<double> tmp_4845;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4845 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_4844 += tmp_4845;
      tmp_4843 += (Conj(ZE(gt2,j2))) * tmp_4844;
   }
   tmp_4842 += tmp_4843;
   tmp_4837 += (vu*Conj(Lambdax)) * tmp_4842;
   std::complex<double> tmp_4846;
   std::complex<double> tmp_4847;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4848;
      std::complex<double> tmp_4849;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4849 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_4848 += tmp_4849;
      tmp_4847 += (ZE(gt3,j2)) * tmp_4848;
   }
   tmp_4846 += tmp_4847;
   tmp_4837 += (vu*Lambdax) * tmp_4846;
   tmp_4836 += (2*Conj(ZH(gt1,2))) * tmp_4837;
   std::complex<double> tmp_4850;
   std::complex<double> tmp_4851;
   std::complex<double> tmp_4852;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4852 += Conj(ZE(gt2,j1))*ZE(gt3,j1);
   }
   tmp_4851 += tmp_4852;
   tmp_4850 += (vu*(0.6*Sqr(g1) - Sqr(g2) - 4*QHu*Ql*Sqr(gp))) * tmp_4851;
   std::complex<double> tmp_4853;
   std::complex<double> tmp_4854;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4854 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_4853 += tmp_4854;
   tmp_4850 += (-2*vu*(0.6*Sqr(g1) + 2*Qe*QHu*Sqr(gp))) * tmp_4853;
   std::complex<double> tmp_4855;
   std::complex<double> tmp_4856;
   std::complex<double> tmp_4857;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4858;
      std::complex<double> tmp_4859;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4859 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_4858 += tmp_4859;
      tmp_4857 += (Conj(ZE(gt2,j2))) * tmp_4858;
   }
   tmp_4856 += tmp_4857;
   tmp_4855 += (Conj(Lambdax)) * tmp_4856;
   std::complex<double> tmp_4860;
   std::complex<double> tmp_4861;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4862;
      std::complex<double> tmp_4863;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4863 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_4862 += tmp_4863;
      tmp_4861 += (ZE(gt3,j2)) * tmp_4862;
   }
   tmp_4860 += tmp_4861;
   tmp_4855 += (Lambdax) * tmp_4860;
   tmp_4850 += (2*vS) * tmp_4855;
   tmp_4836 += (Conj(ZH(gt1,1))) * tmp_4850;
   std::complex<double> tmp_4864;
   std::complex<double> tmp_4865;
   std::complex<double> tmp_4866;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4866 += Conj(ZE(gt2,j1))*ZE(gt3,j1);
   }
   tmp_4865 += tmp_4866;
   tmp_4864 += (vd*(0.6*Sqr(g1) - Sqr(g2) + 4*QHd*Ql*Sqr(gp))) * tmp_4865;
   std::complex<double> tmp_4867;
   std::complex<double> tmp_4868;
   std::complex<double> tmp_4869;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_4869 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_4868 += tmp_4869;
   tmp_4867 += (-(vd*(0.6*Sqr(g1) - 2*Qe*QHd*Sqr(gp)))) * tmp_4868;
   std::complex<double> tmp_4870;
   std::complex<double> tmp_4871;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4872;
      std::complex<double> tmp_4873;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4873 += ZE(gt3,3 + j1)*TYe(j1,j2);
      }
      tmp_4872 += tmp_4873;
      tmp_4871 += (Conj(ZE(gt2,j2))) * tmp_4872;
   }
   tmp_4870 += tmp_4871;
   tmp_4867 += (1.4142135623730951) * tmp_4870;
   std::complex<double> tmp_4874;
   std::complex<double> tmp_4875;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4876;
      std::complex<double> tmp_4877;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4877 += Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_4876 += tmp_4877;
      tmp_4875 += (ZE(gt3,j2)) * tmp_4876;
   }
   tmp_4874 += tmp_4875;
   tmp_4867 += (1.4142135623730951) * tmp_4874;
   std::complex<double> tmp_4878;
   std::complex<double> tmp_4879;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_4880;
      std::complex<double> tmp_4881;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_4882;
         std::complex<double> tmp_4883;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_4883 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_4882 += tmp_4883;
         tmp_4881 += (ZE(gt3,3 + j2)) * tmp_4882;
      }
      tmp_4880 += tmp_4881;
      tmp_4879 += (Conj(ZE(gt2,3 + j3))) * tmp_4880;
   }
   tmp_4878 += tmp_4879;
   tmp_4867 += (2*vd) * tmp_4878;
   std::complex<double> tmp_4884;
   std::complex<double> tmp_4885;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_4886;
      std::complex<double> tmp_4887;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_4888;
         std::complex<double> tmp_4889;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_4889 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_4888 += tmp_4889;
         tmp_4887 += (Conj(ZE(gt2,j2))) * tmp_4888;
      }
      tmp_4886 += tmp_4887;
      tmp_4885 += (ZE(gt3,j3)) * tmp_4886;
   }
   tmp_4884 += tmp_4885;
   tmp_4867 += (2*vd) * tmp_4884;
   tmp_4864 += (2) * tmp_4867;
   tmp_4836 += (-Conj(ZH(gt1,0))) * tmp_4864;
   result += (0.25) * tmp_4836;

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

   std::complex<double> tmp_4890;
   std::complex<double> tmp_4891;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4892;
      std::complex<double> tmp_4893;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4893 += Conj(ZER(gt3,j1))*Ye(j1,j2);
      }
      tmp_4892 += tmp_4893;
      tmp_4891 += (Conj(ZEL(gt1,j2))) * tmp_4892;
   }
   tmp_4890 += tmp_4891;
   result += (-0.7071067811865475*Conj(ZH(gt2,0))) * tmp_4890;

   return result;
}

std::complex<double> UMSSM_effective_couplings::CpFdhhbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_4894;
   std::complex<double> tmp_4895;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4896;
      std::complex<double> tmp_4897;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4897 += Conj(ZDR(gt3,j1))*Yd(j1,j2);
      }
      tmp_4896 += tmp_4897;
      tmp_4895 += (Conj(ZDL(gt1,j2))) * tmp_4896;
   }
   tmp_4894 += tmp_4895;
   result += (-0.7071067811865475*Conj(ZH(gt2,0))) * tmp_4894;

   return result;
}

std::complex<double> UMSSM_effective_couplings::CpFuhhbarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_4898;
   std::complex<double> tmp_4899;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4900;
      std::complex<double> tmp_4901;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4901 += Conj(ZUR(gt3,j1))*Yu(j1,j2);
      }
      tmp_4900 += tmp_4901;
      tmp_4899 += (Conj(ZUL(gt1,j2))) * tmp_4900;
   }
   tmp_4898 += tmp_4899;
   result += (-0.7071067811865475*Conj(ZH(gt2,1))) * tmp_4898;

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

   std::complex<double> tmp_4902;
   std::complex<double> tmp_4903;
   std::complex<double> tmp_4904;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4905;
      std::complex<double> tmp_4906;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4906 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_4905 += tmp_4906;
      tmp_4904 += (Conj(ZD(gt2,j2))) * tmp_4905;
   }
   tmp_4903 += tmp_4904;
   tmp_4902 += (Conj(Lambdax)*(vS*Conj(ZA(gt1,1)) + vu*Conj(ZA(gt1,2)))) *
      tmp_4903;
   std::complex<double> tmp_4907;
   std::complex<double> tmp_4908;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4909;
      std::complex<double> tmp_4910;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4910 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_4909 += tmp_4910;
      tmp_4908 += (ZD(gt3,j2)) * tmp_4909;
   }
   tmp_4907 += tmp_4908;
   tmp_4902 += (-((vS*Conj(ZA(gt1,1)) + vu*Conj(ZA(gt1,2)))*Lambdax)) *
      tmp_4907;
   std::complex<double> tmp_4911;
   std::complex<double> tmp_4912;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4913;
      std::complex<double> tmp_4914;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4914 += ZD(gt3,3 + j1)*TYd(j1,j2);
      }
      tmp_4913 += tmp_4914;
      tmp_4912 += (Conj(ZD(gt2,j2))) * tmp_4913;
   }
   tmp_4911 += tmp_4912;
   std::complex<double> tmp_4915;
   std::complex<double> tmp_4916;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4917;
      std::complex<double> tmp_4918;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4918 += Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_4917 += tmp_4918;
      tmp_4916 += (ZD(gt3,j2)) * tmp_4917;
   }
   tmp_4915 += tmp_4916;
   tmp_4911 += (-1) * tmp_4915;
   tmp_4902 += (1.4142135623730951*Conj(ZA(gt1,0))) * tmp_4911;
   result += (std::complex<double>(0,-0.5)) * tmp_4902;

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

   std::complex<double> tmp_4919;
   std::complex<double> tmp_4920;
   std::complex<double> tmp_4921;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4922;
      std::complex<double> tmp_4923;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4923 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_4922 += tmp_4923;
      tmp_4921 += (Conj(ZU(gt2,j2))) * tmp_4922;
   }
   tmp_4920 += tmp_4921;
   tmp_4919 += (Conj(Lambdax)*(vS*Conj(ZA(gt1,0)) + vd*Conj(ZA(gt1,2)))) *
      tmp_4920;
   std::complex<double> tmp_4924;
   std::complex<double> tmp_4925;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4926;
      std::complex<double> tmp_4927;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4927 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_4926 += tmp_4927;
      tmp_4925 += (ZU(gt3,j2)) * tmp_4926;
   }
   tmp_4924 += tmp_4925;
   tmp_4919 += (-((vS*Conj(ZA(gt1,0)) + vd*Conj(ZA(gt1,2)))*Lambdax)) *
      tmp_4924;
   std::complex<double> tmp_4928;
   std::complex<double> tmp_4929;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4930;
      std::complex<double> tmp_4931;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4931 += ZU(gt3,3 + j1)*TYu(j1,j2);
      }
      tmp_4930 += tmp_4931;
      tmp_4929 += (Conj(ZU(gt2,j2))) * tmp_4930;
   }
   tmp_4928 += tmp_4929;
   std::complex<double> tmp_4932;
   std::complex<double> tmp_4933;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4934;
      std::complex<double> tmp_4935;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4935 += Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_4934 += tmp_4935;
      tmp_4933 += (ZU(gt3,j2)) * tmp_4934;
   }
   tmp_4932 += tmp_4933;
   tmp_4928 += (-1) * tmp_4932;
   tmp_4919 += (1.4142135623730951*Conj(ZA(gt1,1))) * tmp_4928;
   result += (std::complex<double>(0,-0.5)) * tmp_4919;

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

   std::complex<double> tmp_4936;
   std::complex<double> tmp_4937;
   std::complex<double> tmp_4938;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4939;
      std::complex<double> tmp_4940;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4940 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_4939 += tmp_4940;
      tmp_4938 += (Conj(ZE(gt2,j2))) * tmp_4939;
   }
   tmp_4937 += tmp_4938;
   tmp_4936 += (Conj(Lambdax)*(vS*Conj(ZA(gt1,1)) + vu*Conj(ZA(gt1,2)))) *
      tmp_4937;
   std::complex<double> tmp_4941;
   std::complex<double> tmp_4942;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4943;
      std::complex<double> tmp_4944;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4944 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_4943 += tmp_4944;
      tmp_4942 += (ZE(gt3,j2)) * tmp_4943;
   }
   tmp_4941 += tmp_4942;
   tmp_4936 += (-((vS*Conj(ZA(gt1,1)) + vu*Conj(ZA(gt1,2)))*Lambdax)) *
      tmp_4941;
   std::complex<double> tmp_4945;
   std::complex<double> tmp_4946;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4947;
      std::complex<double> tmp_4948;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4948 += ZE(gt3,3 + j1)*TYe(j1,j2);
      }
      tmp_4947 += tmp_4948;
      tmp_4946 += (Conj(ZE(gt2,j2))) * tmp_4947;
   }
   tmp_4945 += tmp_4946;
   std::complex<double> tmp_4949;
   std::complex<double> tmp_4950;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4951;
      std::complex<double> tmp_4952;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4952 += Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_4951 += tmp_4952;
      tmp_4950 += (ZE(gt3,j2)) * tmp_4951;
   }
   tmp_4949 += tmp_4950;
   tmp_4945 += (-1) * tmp_4949;
   tmp_4936 += (1.4142135623730951*Conj(ZA(gt1,0))) * tmp_4945;
   result += (std::complex<double>(0,-0.5)) * tmp_4936;

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

   std::complex<double> tmp_4953;
   std::complex<double> tmp_4954;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4955;
      std::complex<double> tmp_4956;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4956 += Conj(ZER(gt3,j1))*Ye(j1,j2);
      }
      tmp_4955 += tmp_4956;
      tmp_4954 += (Conj(ZEL(gt2,j2))) * tmp_4955;
   }
   tmp_4953 += tmp_4954;
   result += (std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gt1,0))) *
      tmp_4953;

   return result;
}

std::complex<double> UMSSM_effective_couplings::CpAhFdbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_4957;
   std::complex<double> tmp_4958;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4959;
      std::complex<double> tmp_4960;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4960 += Conj(ZDR(gt3,j1))*Yd(j1,j2);
      }
      tmp_4959 += tmp_4960;
      tmp_4958 += (Conj(ZDL(gt2,j2))) * tmp_4959;
   }
   tmp_4957 += tmp_4958;
   result += (std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gt1,0))) *
      tmp_4957;

   return result;
}

std::complex<double> UMSSM_effective_couplings::CpAhFubarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_4961;
   std::complex<double> tmp_4962;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4963;
      std::complex<double> tmp_4964;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_4964 += Conj(ZUR(gt3,j1))*Yu(j1,j2);
      }
      tmp_4963 += tmp_4964;
      tmp_4962 += (Conj(ZUL(gt2,j2))) * tmp_4963;
   }
   tmp_4961 += tmp_4962;
   result += (std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gt1,1))) *
      tmp_4961;

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
