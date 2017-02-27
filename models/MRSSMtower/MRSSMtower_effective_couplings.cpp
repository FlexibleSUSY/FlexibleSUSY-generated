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

// File generated at Mon 27 Feb 2017 13:27:20

#include "MRSSMtower_effective_couplings.hpp"

#include "effective_couplings.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

using namespace effective_couplings;

#define INPUTPARAMETER(parameter) model.get_input().parameter
#define MODELPARAMETER(parameter) model.get_##parameter()
#define DERIVEDPARAMETER(parameter) model.##parameter()
#define PHASE(parameter) model.get_##parameter()
#define PHYSICAL(parameter) model.get_physical().parameter

MRSSMtower_effective_couplings::MRSSMtower_effective_couplings(
   const MRSSMtower_mass_eigenstates& model_,
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

void MRSSMtower_effective_couplings::calculate_effective_couplings()
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

void MRSSMtower_effective_couplings::set_model(const MRSSMtower_mass_eigenstates& model_)
{
   model = model_;
   copy_mixing_matrices_from_model();
}

void MRSSMtower_effective_couplings::copy_mixing_matrices_from_model()
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

standard_model::Standard_model MRSSMtower_effective_couplings::initialise_SM() const
{
   standard_model::Standard_model sm;

   sm.set_loops(2);
   sm.set_thresholds(2);
   sm.set_physical_input(physical_input);

   sm.initialise_from_input(qedqcd);

   return sm;
}

void MRSSMtower_effective_couplings::run_SM_strong_coupling_to(standard_model::Standard_model sm, double m)
{
   sm.run_to(m);

   model.set_g3(sm.get_g3());

}

std::complex<double> MRSSMtower_effective_couplings::scalar_scalar_qcd_factor(double m_decay, double m_loop) const
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

std::complex<double> MRSSMtower_effective_couplings::scalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) *
         scalar_diphoton_fermion_loop(m_decay, m_loop);

   }

   return result;
}

std::complex<double> MRSSMtower_effective_couplings::pseudoscalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) *
         pseudoscalar_diphoton_fermion_loop(m_decay, m_loop);

   }

   return result;
}

double MRSSMtower_effective_couplings::number_of_active_flavours(double m) const
{
   if (m < qedqcd.displayMbMb()) {
      return 4.0;
   } else if (m < qedqcd.displayPoleMt()) {
      return 5.0;
   } else {
      return 6.0;
   }
}

double MRSSMtower_effective_couplings::scalar_scaling_factor(double m) const
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

double MRSSMtower_effective_couplings::pseudoscalar_scaling_factor(double m) const
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

double MRSSMtower_effective_couplings::get_hhVPVP_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(Mhh)(gO1);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CphhVPVP(gO1));
}

double MRSSMtower_effective_couplings::get_hhVGVG_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(Mhh)(gO1);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CphhVGVG(gO1));
}

double MRSSMtower_effective_couplings::get_AhVPVP_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(MAh)(gO1);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CpAhVPVP(gO1));
}

double MRSSMtower_effective_couplings::get_AhVGVG_partial_width(unsigned gO1) const
{
   const double mass = PHYSICAL(MAh)(gO1);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CpAhVGVG(gO1));
}

std::complex<double> MRSSMtower_effective_couplings::CphhSRdpconjSRdp(unsigned gt1) const
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

std::complex<double> MRSSMtower_effective_couplings::CphhSRumconjSRum(unsigned gt1) const
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

std::complex<double> MRSSMtower_effective_couplings::CphhSdconjSd(unsigned gt1, unsigned gt2, unsigned gt3) const
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

   std::complex<double> tmp_3776;
   std::complex<double> tmp_3777;
   std::complex<double> tmp_3778;
   std::complex<double> tmp_3779;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3780;
      std::complex<double> tmp_3781;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3782;
         std::complex<double> tmp_3783;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3783 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_3782 += tmp_3783;
         tmp_3781 += (ZD(gt3,3 + j2)) * tmp_3782;
      }
      tmp_3780 += tmp_3781;
      tmp_3779 += (Conj(ZD(gt2,3 + j3))) * tmp_3780;
   }
   tmp_3778 += tmp_3779;
   tmp_3777 += (-6*vd*ZH(gt1,0)) * tmp_3778;
   std::complex<double> tmp_3784;
   std::complex<double> tmp_3785;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3786;
      std::complex<double> tmp_3787;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3788;
         std::complex<double> tmp_3789;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3789 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_3788 += tmp_3789;
         tmp_3787 += (Conj(ZD(gt2,j2))) * tmp_3788;
      }
      tmp_3786 += tmp_3787;
      tmp_3785 += (ZD(gt3,j3)) * tmp_3786;
   }
   tmp_3784 += tmp_3785;
   tmp_3777 += (-6*vd*ZH(gt1,0)) * tmp_3784;
   std::complex<double> tmp_3790;
   std::complex<double> tmp_3791;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3792;
      std::complex<double> tmp_3793;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3793 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3792 += tmp_3793;
      tmp_3791 += (Conj(ZD(gt2,j2))) * tmp_3792;
   }
   tmp_3790 += tmp_3791;
   tmp_3777 += (4.242640687119286*Conj(Mu)*ZH(gt1,1)) * tmp_3790;
   std::complex<double> tmp_3794;
   std::complex<double> tmp_3795;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3796;
      std::complex<double> tmp_3797;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3797 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3796 += tmp_3797;
      tmp_3795 += (ZD(gt3,j2)) * tmp_3796;
   }
   tmp_3794 += tmp_3795;
   tmp_3777 += (4.242640687119286*Mu*ZH(gt1,1)) * tmp_3794;
   std::complex<double> tmp_3798;
   std::complex<double> tmp_3799;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3799 += Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1);
   }
   tmp_3798 += tmp_3799;
   tmp_3777 += (0.7745966692414834*g1*(0.7745966692414834*g1*vd*ZH(gt1,0) -
      0.7745966692414834*g1*vu*ZH(gt1,1) - 2*(MDBS + Conj(MDBS))*ZH(gt1,2))) *
      tmp_3798;
   tmp_3776 += (2) * tmp_3777;
   std::complex<double> tmp_3800;
   std::complex<double> tmp_3801;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3801 += Conj(ZD(gt2,j1))*ZD(gt3,j1);
   }
   tmp_3800 += tmp_3801;
   tmp_3776 += (vd*(0.6*Sqr(g1) + 3*Sqr(g2))*ZH(gt1,0) - vu*(0.6*Sqr(g1) + 3*
      Sqr(g2))*ZH(gt1,1) - 1.5491933384829668*g1*MDBS*ZH(gt1,2) -
      1.5491933384829668*g1*Conj(MDBS)*ZH(gt1,2) + 6*g2*MDWBT*ZH(gt1,3) + 6*g2*
      Conj(MDWBT)*ZH(gt1,3)) * tmp_3800;
   result += (0.08333333333333333) * tmp_3776;

   return result;
}

std::complex<double> MRSSMtower_effective_couplings::CphhSuconjSu(unsigned gt1, unsigned gt2, unsigned gt3) const
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

   std::complex<double> tmp_3802;
   std::complex<double> tmp_3803;
   std::complex<double> tmp_3804;
   std::complex<double> tmp_3805;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3806;
      std::complex<double> tmp_3807;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3807 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_3806 += tmp_3807;
      tmp_3805 += (Conj(ZU(gt2,j2))) * tmp_3806;
   }
   tmp_3804 += tmp_3805;
   tmp_3803 += (1.4142135623730951*Conj(Mu)*ZH(gt1,0)) * tmp_3804;
   std::complex<double> tmp_3808;
   std::complex<double> tmp_3809;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3810;
      std::complex<double> tmp_3811;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3811 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3810 += tmp_3811;
      tmp_3809 += (ZU(gt3,j2)) * tmp_3810;
   }
   tmp_3808 += tmp_3809;
   tmp_3803 += (1.4142135623730951*Mu*ZH(gt1,0)) * tmp_3808;
   std::complex<double> tmp_3812;
   std::complex<double> tmp_3813;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3814;
      std::complex<double> tmp_3815;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3816;
         std::complex<double> tmp_3817;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3817 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_3816 += tmp_3817;
         tmp_3815 += (ZU(gt3,3 + j2)) * tmp_3816;
      }
      tmp_3814 += tmp_3815;
      tmp_3813 += (Conj(ZU(gt2,3 + j3))) * tmp_3814;
   }
   tmp_3812 += tmp_3813;
   std::complex<double> tmp_3818;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3819;
      std::complex<double> tmp_3820;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3821;
         std::complex<double> tmp_3822;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3822 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_3821 += tmp_3822;
         tmp_3820 += (Conj(ZU(gt2,j2))) * tmp_3821;
      }
      tmp_3819 += tmp_3820;
      tmp_3818 += (ZU(gt3,j3)) * tmp_3819;
   }
   tmp_3812 += tmp_3818;
   tmp_3803 += (-2*vu*ZH(gt1,1)) * tmp_3812;
   tmp_3802 += (6) * tmp_3803;
   std::complex<double> tmp_3823;
   std::complex<double> tmp_3824;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3824 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_3823 += tmp_3824;
   tmp_3802 += (3.0983866769659336*g1*(-0.7745966692414834*g1*vd*ZH(gt1,0) +
      0.7745966692414834*g1*vu*ZH(gt1,1) + 2*(MDBS + Conj(MDBS))*ZH(gt1,2))) *
      tmp_3823;
   std::complex<double> tmp_3825;
   std::complex<double> tmp_3826;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3826 += Conj(ZU(gt2,j1))*ZU(gt3,j1);
   }
   tmp_3825 += tmp_3826;
   tmp_3802 += (vd*(0.6*Sqr(g1) - 3*Sqr(g2))*ZH(gt1,0) - vu*(0.6*Sqr(g1) - 3*
      Sqr(g2))*ZH(gt1,1) - 2*(0.7745966692414834*g1*(MDBS + Conj(MDBS))*ZH(gt1,2)
      + 3*g2*(MDWBT + Conj(MDWBT))*ZH(gt1,3))) * tmp_3825;
   result += (0.08333333333333333) * tmp_3802;

   return result;
}

std::complex<double> MRSSMtower_effective_couplings::CphhSeconjSe(unsigned gt1, unsigned gt2, unsigned gt3) const
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

   std::complex<double> tmp_3827;
   std::complex<double> tmp_3828;
   std::complex<double> tmp_3829;
   std::complex<double> tmp_3830;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3831;
      std::complex<double> tmp_3832;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3833;
         std::complex<double> tmp_3834;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3834 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_3833 += tmp_3834;
         tmp_3832 += (ZE(gt3,3 + j2)) * tmp_3833;
      }
      tmp_3831 += tmp_3832;
      tmp_3830 += (Conj(ZE(gt2,3 + j3))) * tmp_3831;
   }
   tmp_3829 += tmp_3830;
   tmp_3828 += (-2*vd*ZH(gt1,0)) * tmp_3829;
   std::complex<double> tmp_3835;
   std::complex<double> tmp_3836;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3837;
      std::complex<double> tmp_3838;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3839;
         std::complex<double> tmp_3840;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3840 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_3839 += tmp_3840;
         tmp_3838 += (Conj(ZE(gt2,j2))) * tmp_3839;
      }
      tmp_3837 += tmp_3838;
      tmp_3836 += (ZE(gt3,j3)) * tmp_3837;
   }
   tmp_3835 += tmp_3836;
   tmp_3828 += (-2*vd*ZH(gt1,0)) * tmp_3835;
   std::complex<double> tmp_3841;
   std::complex<double> tmp_3842;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3843;
      std::complex<double> tmp_3844;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3844 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3843 += tmp_3844;
      tmp_3842 += (Conj(ZE(gt2,j2))) * tmp_3843;
   }
   tmp_3841 += tmp_3842;
   tmp_3828 += (1.4142135623730951*Conj(Mu)*ZH(gt1,1)) * tmp_3841;
   std::complex<double> tmp_3845;
   std::complex<double> tmp_3846;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3847;
      std::complex<double> tmp_3848;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3848 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3847 += tmp_3848;
      tmp_3846 += (ZE(gt3,j2)) * tmp_3847;
   }
   tmp_3845 += tmp_3846;
   tmp_3828 += (1.4142135623730951*Mu*ZH(gt1,1)) * tmp_3845;
   std::complex<double> tmp_3849;
   std::complex<double> tmp_3850;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3850 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_3849 += tmp_3850;
   tmp_3828 += (0.7745966692414834*g1*(0.7745966692414834*g1*vd*ZH(gt1,0) -
      0.7745966692414834*g1*vu*ZH(gt1,1) - 2*(MDBS + Conj(MDBS))*ZH(gt1,2))) *
      tmp_3849;
   tmp_3827 += (2) * tmp_3828;
   std::complex<double> tmp_3851;
   std::complex<double> tmp_3852;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3852 += Conj(ZE(gt2,j1))*ZE(gt3,j1);
   }
   tmp_3851 += tmp_3852;
   tmp_3827 += (vd*(-0.6*Sqr(g1) + Sqr(g2))*ZH(gt1,0) + vu*(0.6*Sqr(g1) - Sqr(
      g2))*ZH(gt1,1) + 2*(0.7745966692414834*g1*(MDBS + Conj(MDBS))*ZH(gt1,2) + g2
      *(MDWBT + Conj(MDWBT))*ZH(gt1,3))) * tmp_3851;
   result += (0.25) * tmp_3827;

   return result;
}

std::complex<double> MRSSMtower_effective_couplings::CphhHpmconjHpm(unsigned gt1, unsigned gt2, unsigned gt3) const
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

std::complex<double> MRSSMtower_effective_couplings::CpCha1hhbarCha1PL(unsigned gt1, unsigned gt2, unsigned gt3) const
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

std::complex<double> MRSSMtower_effective_couplings::CpCha2hhbarCha2PL(unsigned gt1, unsigned gt2, unsigned gt3) const
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

std::complex<double> MRSSMtower_effective_couplings::CpFehhbarFePL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Ye = MODELPARAMETER(Ye);

   std::complex<double> result;

   std::complex<double> tmp_3853;
   std::complex<double> tmp_3854;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3855;
      std::complex<double> tmp_3856;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3856 += Conj(ZER(gt3,j1))*Ye(j1,j2);
      }
      tmp_3855 += tmp_3856;
      tmp_3854 += (Conj(ZEL(gt1,j2))) * tmp_3855;
   }
   tmp_3853 += tmp_3854;
   result += (-0.7071067811865475*ZH(gt2,0)) * tmp_3853;

   return result;
}

std::complex<double> MRSSMtower_effective_couplings::CpFdhhbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_3857;
   std::complex<double> tmp_3858;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3859;
      std::complex<double> tmp_3860;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3860 += Conj(ZDR(gt3,j1))*Yd(j1,j2);
      }
      tmp_3859 += tmp_3860;
      tmp_3858 += (Conj(ZDL(gt1,j2))) * tmp_3859;
   }
   tmp_3857 += tmp_3858;
   result += (-0.7071067811865475*ZH(gt2,0)) * tmp_3857;

   return result;
}

std::complex<double> MRSSMtower_effective_couplings::CpFuhhbarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_3861;
   std::complex<double> tmp_3862;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3863;
      std::complex<double> tmp_3864;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3864 += Conj(ZUR(gt3,j1))*Yu(j1,j2);
      }
      tmp_3863 += tmp_3864;
      tmp_3862 += (Conj(ZUL(gt1,j2))) * tmp_3863;
   }
   tmp_3861 += tmp_3862;
   result += (-0.7071067811865475*ZH(gt2,1)) * tmp_3861;

   return result;
}

std::complex<double> MRSSMtower_effective_couplings::CphhVWmconjVWm(unsigned gt1) const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vT = MODELPARAMETER(vT);
   const auto vu = MODELPARAMETER(vu);

   std::complex<double> result;

   result = 0.5*Sqr(g2)*(vd*ZH(gt1,0) + vu*ZH(gt1,1) + 4*vT*ZH(gt1,3));

   return result;
}

std::complex<double> MRSSMtower_effective_couplings::CpAhSRdpconjSRdp(unsigned gt1) const
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

std::complex<double> MRSSMtower_effective_couplings::CpAhSRumconjSRum(unsigned gt1) const
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

std::complex<double> MRSSMtower_effective_couplings::CpAhSdconjSd(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3865;
   std::complex<double> tmp_3866;
   std::complex<double> tmp_3867;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3868;
      std::complex<double> tmp_3869;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3869 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3868 += tmp_3869;
      tmp_3867 += (Conj(ZD(gt2,j2))) * tmp_3868;
   }
   tmp_3866 += tmp_3867;
   tmp_3865 += (4.242640687119286*Conj(Mu)*ZA(gt1,1)) * tmp_3866;
   std::complex<double> tmp_3870;
   std::complex<double> tmp_3871;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3872;
      std::complex<double> tmp_3873;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3873 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3872 += tmp_3873;
      tmp_3871 += (ZD(gt3,j2)) * tmp_3872;
   }
   tmp_3870 += tmp_3871;
   tmp_3865 += (-4.242640687119286*Mu*ZA(gt1,1)) * tmp_3870;
   std::complex<double> tmp_3874;
   std::complex<double> tmp_3875;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3875 += Conj(ZD(gt2,j1))*ZD(gt3,j1);
   }
   tmp_3874 += tmp_3875;
   tmp_3865 += (0.7745966692414834*g1*MDBS*ZA(gt1,2)) * tmp_3874;
   std::complex<double> tmp_3876;
   std::complex<double> tmp_3877;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3877 += Conj(ZD(gt2,j1))*ZD(gt3,j1);
   }
   tmp_3876 += tmp_3877;
   tmp_3865 += (-0.7745966692414834*g1*Conj(MDBS)*ZA(gt1,2)) * tmp_3876;
   std::complex<double> tmp_3878;
   std::complex<double> tmp_3879;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3879 += Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1);
   }
   tmp_3878 += tmp_3879;
   tmp_3865 += (1.5491933384829668*g1*MDBS*ZA(gt1,2)) * tmp_3878;
   std::complex<double> tmp_3880;
   std::complex<double> tmp_3881;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3881 += Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1);
   }
   tmp_3880 += tmp_3881;
   tmp_3865 += (-1.5491933384829668*g1*Conj(MDBS)*ZA(gt1,2)) * tmp_3880;
   std::complex<double> tmp_3882;
   std::complex<double> tmp_3883;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3883 += Conj(ZD(gt2,j1))*ZD(gt3,j1);
   }
   tmp_3882 += tmp_3883;
   tmp_3865 += (-3*g2*MDWBT*ZA(gt1,3)) * tmp_3882;
   std::complex<double> tmp_3884;
   std::complex<double> tmp_3885;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3885 += Conj(ZD(gt2,j1))*ZD(gt3,j1);
   }
   tmp_3884 += tmp_3885;
   tmp_3865 += (3*g2*Conj(MDWBT)*ZA(gt1,3)) * tmp_3884;
   result += (std::complex<double>(0,-0.16666666666666666)) * tmp_3865;

   return result;
}

std::complex<double> MRSSMtower_effective_couplings::CpAhSuconjSu(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3886;
   std::complex<double> tmp_3887;
   std::complex<double> tmp_3888;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3889;
      std::complex<double> tmp_3890;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3890 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_3889 += tmp_3890;
      tmp_3888 += (Conj(ZU(gt2,j2))) * tmp_3889;
   }
   tmp_3887 += tmp_3888;
   tmp_3886 += (4.242640687119286*Conj(Mu)*ZA(gt1,0)) * tmp_3887;
   std::complex<double> tmp_3891;
   std::complex<double> tmp_3892;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3893;
      std::complex<double> tmp_3894;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3894 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3893 += tmp_3894;
      tmp_3892 += (ZU(gt3,j2)) * tmp_3893;
   }
   tmp_3891 += tmp_3892;
   tmp_3886 += (-4.242640687119286*Mu*ZA(gt1,0)) * tmp_3891;
   std::complex<double> tmp_3895;
   std::complex<double> tmp_3896;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3896 += Conj(ZU(gt2,j1))*ZU(gt3,j1);
   }
   tmp_3895 += tmp_3896;
   tmp_3886 += (0.7745966692414834*g1*MDBS*ZA(gt1,2)) * tmp_3895;
   std::complex<double> tmp_3897;
   std::complex<double> tmp_3898;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3898 += Conj(ZU(gt2,j1))*ZU(gt3,j1);
   }
   tmp_3897 += tmp_3898;
   tmp_3886 += (-0.7745966692414834*g1*Conj(MDBS)*ZA(gt1,2)) * tmp_3897;
   std::complex<double> tmp_3899;
   std::complex<double> tmp_3900;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3900 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_3899 += tmp_3900;
   tmp_3886 += (-3.0983866769659336*g1*MDBS*ZA(gt1,2)) * tmp_3899;
   std::complex<double> tmp_3901;
   std::complex<double> tmp_3902;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3902 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_3901 += tmp_3902;
   tmp_3886 += (3.0983866769659336*g1*Conj(MDBS)*ZA(gt1,2)) * tmp_3901;
   std::complex<double> tmp_3903;
   std::complex<double> tmp_3904;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3904 += Conj(ZU(gt2,j1))*ZU(gt3,j1);
   }
   tmp_3903 += tmp_3904;
   tmp_3886 += (3*g2*MDWBT*ZA(gt1,3)) * tmp_3903;
   std::complex<double> tmp_3905;
   std::complex<double> tmp_3906;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3906 += Conj(ZU(gt2,j1))*ZU(gt3,j1);
   }
   tmp_3905 += tmp_3906;
   tmp_3886 += (-3*g2*Conj(MDWBT)*ZA(gt1,3)) * tmp_3905;
   result += (std::complex<double>(0,-0.16666666666666666)) * tmp_3886;

   return result;
}

std::complex<double> MRSSMtower_effective_couplings::CpAhSeconjSe(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3907;
   std::complex<double> tmp_3908;
   std::complex<double> tmp_3909;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3910;
      std::complex<double> tmp_3911;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3911 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3910 += tmp_3911;
      tmp_3909 += (Conj(ZE(gt2,j2))) * tmp_3910;
   }
   tmp_3908 += tmp_3909;
   tmp_3907 += (1.4142135623730951*Conj(Mu)*ZA(gt1,1)) * tmp_3908;
   std::complex<double> tmp_3912;
   std::complex<double> tmp_3913;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3914;
      std::complex<double> tmp_3915;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3915 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3914 += tmp_3915;
      tmp_3913 += (ZE(gt3,j2)) * tmp_3914;
   }
   tmp_3912 += tmp_3913;
   tmp_3907 += (-1.4142135623730951*Mu*ZA(gt1,1)) * tmp_3912;
   std::complex<double> tmp_3916;
   std::complex<double> tmp_3917;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3917 += Conj(ZE(gt2,j1))*ZE(gt3,j1);
   }
   tmp_3916 += tmp_3917;
   tmp_3907 += (-0.7745966692414834*g1*MDBS*ZA(gt1,2)) * tmp_3916;
   std::complex<double> tmp_3918;
   std::complex<double> tmp_3919;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3919 += Conj(ZE(gt2,j1))*ZE(gt3,j1);
   }
   tmp_3918 += tmp_3919;
   tmp_3907 += (0.7745966692414834*g1*Conj(MDBS)*ZA(gt1,2)) * tmp_3918;
   std::complex<double> tmp_3920;
   std::complex<double> tmp_3921;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3921 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_3920 += tmp_3921;
   tmp_3907 += (1.5491933384829668*g1*MDBS*ZA(gt1,2)) * tmp_3920;
   std::complex<double> tmp_3922;
   std::complex<double> tmp_3923;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3923 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_3922 += tmp_3923;
   tmp_3907 += (-1.5491933384829668*g1*Conj(MDBS)*ZA(gt1,2)) * tmp_3922;
   std::complex<double> tmp_3924;
   std::complex<double> tmp_3925;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3925 += Conj(ZE(gt2,j1))*ZE(gt3,j1);
   }
   tmp_3924 += tmp_3925;
   tmp_3907 += (-(g2*MDWBT*ZA(gt1,3))) * tmp_3924;
   std::complex<double> tmp_3926;
   std::complex<double> tmp_3927;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3927 += Conj(ZE(gt2,j1))*ZE(gt3,j1);
   }
   tmp_3926 += tmp_3927;
   tmp_3907 += (g2*Conj(MDWBT)*ZA(gt1,3)) * tmp_3926;
   result += (std::complex<double>(0,-0.5)) * tmp_3907;

   return result;
}

std::complex<double> MRSSMtower_effective_couplings::CpAhHpmconjHpm(unsigned gt1, unsigned gt2, unsigned gt3) const
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

std::complex<double> MRSSMtower_effective_couplings::CpAhCha1barCha1PL(unsigned gt1, unsigned gt2, unsigned gt3) const
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

std::complex<double> MRSSMtower_effective_couplings::CpAhCha2barCha2PL(unsigned gt1, unsigned gt2, unsigned gt3) const
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

std::complex<double> MRSSMtower_effective_couplings::CpAhFebarFePL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Ye = MODELPARAMETER(Ye);

   std::complex<double> result;

   std::complex<double> tmp_3928;
   std::complex<double> tmp_3929;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3930;
      std::complex<double> tmp_3931;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3931 += Conj(ZER(gt3,j1))*Ye(j1,j2);
      }
      tmp_3930 += tmp_3931;
      tmp_3929 += (Conj(ZEL(gt2,j2))) * tmp_3930;
   }
   tmp_3928 += tmp_3929;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,0)) *
      tmp_3928;

   return result;
}

std::complex<double> MRSSMtower_effective_couplings::CpAhFdbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_3932;
   std::complex<double> tmp_3933;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3934;
      std::complex<double> tmp_3935;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3935 += Conj(ZDR(gt3,j1))*Yd(j1,j2);
      }
      tmp_3934 += tmp_3935;
      tmp_3933 += (Conj(ZDL(gt2,j2))) * tmp_3934;
   }
   tmp_3932 += tmp_3933;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,0)) *
      tmp_3932;

   return result;
}

std::complex<double> MRSSMtower_effective_couplings::CpAhFubarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_3936;
   std::complex<double> tmp_3937;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3938;
      std::complex<double> tmp_3939;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3939 += Conj(ZUR(gt3,j1))*Yu(j1,j2);
      }
      tmp_3938 += tmp_3939;
      tmp_3937 += (Conj(ZUL(gt2,j2))) * tmp_3938;
   }
   tmp_3936 += tmp_3937;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,1)) *
      tmp_3936;

   return result;
}

void MRSSMtower_effective_couplings::calculate_eff_CphhVPVP(unsigned gO1)
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

void MRSSMtower_effective_couplings::calculate_eff_CphhVGVG(unsigned gO1)
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

void MRSSMtower_effective_couplings::calculate_eff_CpAhVPVP(unsigned gO1)
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

void MRSSMtower_effective_couplings::calculate_eff_CpAhVGVG(unsigned gO1)
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
