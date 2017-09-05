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

// File generated at Tue 5 Sep 2017 12:50:23

#include "MSSMRHN_effective_couplings.hpp"

#include "effective_couplings.hpp"
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

standard_model::Standard_model MSSMRHN_effective_couplings::initialise_SM() const
{
   standard_model::Standard_model sm;

   sm.set_loops(2);
   sm.set_thresholds(2);
   sm.set_physical_input(physical_input);

   sm.initialise_from_input(qedqcd);

   return sm;
}

void MSSMRHN_effective_couplings::run_SM_strong_coupling_to(standard_model::Standard_model sm, double m)
{
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

   std::complex<double> tmp_3807;
   std::complex<double> tmp_3808;
   std::complex<double> tmp_3809;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3809 += Conj(ZD(gt2,j1))*ZD(gt3,j1);
   }
   tmp_3808 += tmp_3809;
   tmp_3807 += ((0.6*Sqr(g1) + 3*Sqr(g2))*(vd*ZH(gt1,0) - vu*ZH(gt1,1))) *
      tmp_3808;
   std::complex<double> tmp_3810;
   std::complex<double> tmp_3811;
   std::complex<double> tmp_3812;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3812 += Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1);
   }
   tmp_3811 += tmp_3812;
   tmp_3810 += (0.6*Sqr(g1)*(vd*ZH(gt1,0) - vu*ZH(gt1,1))) * tmp_3811;
   std::complex<double> tmp_3813;
   std::complex<double> tmp_3814;
   std::complex<double> tmp_3815;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3816;
      std::complex<double> tmp_3817;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3817 += ZD(gt3,3 + j1)*TYd(j1,j2);
      }
      tmp_3816 += tmp_3817;
      tmp_3815 += (Conj(ZD(gt2,j2))) * tmp_3816;
   }
   tmp_3814 += tmp_3815;
   tmp_3813 += (1.4142135623730951*ZH(gt1,0)) * tmp_3814;
   std::complex<double> tmp_3818;
   std::complex<double> tmp_3819;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3820;
      std::complex<double> tmp_3821;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3821 += Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_3820 += tmp_3821;
      tmp_3819 += (ZD(gt3,j2)) * tmp_3820;
   }
   tmp_3818 += tmp_3819;
   tmp_3813 += (1.4142135623730951*ZH(gt1,0)) * tmp_3818;
   std::complex<double> tmp_3822;
   std::complex<double> tmp_3823;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3824;
      std::complex<double> tmp_3825;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3826;
         std::complex<double> tmp_3827;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3827 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_3826 += tmp_3827;
         tmp_3825 += (ZD(gt3,3 + j2)) * tmp_3826;
      }
      tmp_3824 += tmp_3825;
      tmp_3823 += (Conj(ZD(gt2,3 + j3))) * tmp_3824;
   }
   tmp_3822 += tmp_3823;
   tmp_3813 += (2*vd*ZH(gt1,0)) * tmp_3822;
   std::complex<double> tmp_3828;
   std::complex<double> tmp_3829;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3830;
      std::complex<double> tmp_3831;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3832;
         std::complex<double> tmp_3833;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3833 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_3832 += tmp_3833;
         tmp_3831 += (Conj(ZD(gt2,j2))) * tmp_3832;
      }
      tmp_3830 += tmp_3831;
      tmp_3829 += (ZD(gt3,j3)) * tmp_3830;
   }
   tmp_3828 += tmp_3829;
   tmp_3813 += (2*vd*ZH(gt1,0)) * tmp_3828;
   std::complex<double> tmp_3834;
   std::complex<double> tmp_3835;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3836;
      std::complex<double> tmp_3837;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3837 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3836 += tmp_3837;
      tmp_3835 += (Conj(ZD(gt2,j2))) * tmp_3836;
   }
   tmp_3834 += tmp_3835;
   tmp_3813 += (-1.4142135623730951*Conj(Mu)*ZH(gt1,1)) * tmp_3834;
   std::complex<double> tmp_3838;
   std::complex<double> tmp_3839;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3840;
      std::complex<double> tmp_3841;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3841 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3840 += tmp_3841;
      tmp_3839 += (ZD(gt3,j2)) * tmp_3840;
   }
   tmp_3838 += tmp_3839;
   tmp_3813 += (-1.4142135623730951*Mu*ZH(gt1,1)) * tmp_3838;
   tmp_3810 += (-3) * tmp_3813;
   tmp_3807 += (2) * tmp_3810;
   result += (0.08333333333333333) * tmp_3807;

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

   std::complex<double> tmp_3842;
   std::complex<double> tmp_3843;
   std::complex<double> tmp_3844;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3844 += Conj(ZU(gt2,j1))*ZU(gt3,j1);
   }
   tmp_3843 += tmp_3844;
   tmp_3842 += ((0.6*Sqr(g1) - 3*Sqr(g2))*(vd*ZH(gt1,0) - vu*ZH(gt1,1))) *
      tmp_3843;
   std::complex<double> tmp_3845;
   std::complex<double> tmp_3846;
   std::complex<double> tmp_3847;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3848;
      std::complex<double> tmp_3849;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3849 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_3848 += tmp_3849;
      tmp_3847 += (Conj(ZU(gt2,j2))) * tmp_3848;
   }
   tmp_3846 += tmp_3847;
   tmp_3845 += (-4.242640687119286*Conj(Mu)*ZH(gt1,0)) * tmp_3846;
   std::complex<double> tmp_3850;
   std::complex<double> tmp_3851;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3852;
      std::complex<double> tmp_3853;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3853 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3852 += tmp_3853;
      tmp_3851 += (ZU(gt3,j2)) * tmp_3852;
   }
   tmp_3850 += tmp_3851;
   tmp_3845 += (-4.242640687119286*Mu*ZH(gt1,0)) * tmp_3850;
   std::complex<double> tmp_3854;
   std::complex<double> tmp_3855;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3856;
      std::complex<double> tmp_3857;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3857 += ZU(gt3,3 + j1)*TYu(j1,j2);
      }
      tmp_3856 += tmp_3857;
      tmp_3855 += (Conj(ZU(gt2,j2))) * tmp_3856;
   }
   tmp_3854 += tmp_3855;
   tmp_3845 += (4.242640687119286*ZH(gt1,1)) * tmp_3854;
   std::complex<double> tmp_3858;
   std::complex<double> tmp_3859;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3860;
      std::complex<double> tmp_3861;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3861 += Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_3860 += tmp_3861;
      tmp_3859 += (ZU(gt3,j2)) * tmp_3860;
   }
   tmp_3858 += tmp_3859;
   tmp_3845 += (4.242640687119286*ZH(gt1,1)) * tmp_3858;
   std::complex<double> tmp_3862;
   std::complex<double> tmp_3863;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3864;
      std::complex<double> tmp_3865;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3866;
         std::complex<double> tmp_3867;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3867 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_3866 += tmp_3867;
         tmp_3865 += (ZU(gt3,3 + j2)) * tmp_3866;
      }
      tmp_3864 += tmp_3865;
      tmp_3863 += (Conj(ZU(gt2,3 + j3))) * tmp_3864;
   }
   tmp_3862 += tmp_3863;
   tmp_3845 += (6*vu*ZH(gt1,1)) * tmp_3862;
   std::complex<double> tmp_3868;
   std::complex<double> tmp_3869;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3870;
      std::complex<double> tmp_3871;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3872;
         std::complex<double> tmp_3873;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3873 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_3872 += tmp_3873;
         tmp_3871 += (Conj(ZU(gt2,j2))) * tmp_3872;
      }
      tmp_3870 += tmp_3871;
      tmp_3869 += (ZU(gt3,j3)) * tmp_3870;
   }
   tmp_3868 += tmp_3869;
   tmp_3845 += (6*vu*ZH(gt1,1)) * tmp_3868;
   std::complex<double> tmp_3874;
   std::complex<double> tmp_3875;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3875 += Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1);
   }
   tmp_3874 += tmp_3875;
   tmp_3845 += (1.2*Sqr(g1)*(vd*ZH(gt1,0) - vu*ZH(gt1,1))) * tmp_3874;
   tmp_3842 += (-2) * tmp_3845;
   result += (0.08333333333333333) * tmp_3842;

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

   std::complex<double> tmp_3876;
   std::complex<double> tmp_3877;
   std::complex<double> tmp_3878;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3878 += Conj(ZE(gt2,j1))*ZE(gt3,j1);
   }
   tmp_3877 += tmp_3878;
   tmp_3876 += (-((0.6*Sqr(g1) - Sqr(g2))*(vd*ZH(gt1,0) - vu*ZH(gt1,1)))) *
      tmp_3877;
   std::complex<double> tmp_3879;
   std::complex<double> tmp_3880;
   std::complex<double> tmp_3881;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3882;
      std::complex<double> tmp_3883;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3883 += ZE(gt3,3 + j1)*TYe(j1,j2);
      }
      tmp_3882 += tmp_3883;
      tmp_3881 += (Conj(ZE(gt2,j2))) * tmp_3882;
   }
   tmp_3880 += tmp_3881;
   tmp_3879 += (-1.4142135623730951*ZH(gt1,0)) * tmp_3880;
   std::complex<double> tmp_3884;
   std::complex<double> tmp_3885;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3886;
      std::complex<double> tmp_3887;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3887 += Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_3886 += tmp_3887;
      tmp_3885 += (ZE(gt3,j2)) * tmp_3886;
   }
   tmp_3884 += tmp_3885;
   tmp_3879 += (-1.4142135623730951*ZH(gt1,0)) * tmp_3884;
   std::complex<double> tmp_3888;
   std::complex<double> tmp_3889;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3890;
      std::complex<double> tmp_3891;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3892;
         std::complex<double> tmp_3893;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3893 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_3892 += tmp_3893;
         tmp_3891 += (ZE(gt3,3 + j2)) * tmp_3892;
      }
      tmp_3890 += tmp_3891;
      tmp_3889 += (Conj(ZE(gt2,3 + j3))) * tmp_3890;
   }
   tmp_3888 += tmp_3889;
   tmp_3879 += (-2*vd*ZH(gt1,0)) * tmp_3888;
   std::complex<double> tmp_3894;
   std::complex<double> tmp_3895;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3896;
      std::complex<double> tmp_3897;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3898;
         std::complex<double> tmp_3899;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3899 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_3898 += tmp_3899;
         tmp_3897 += (Conj(ZE(gt2,j2))) * tmp_3898;
      }
      tmp_3896 += tmp_3897;
      tmp_3895 += (ZE(gt3,j3)) * tmp_3896;
   }
   tmp_3894 += tmp_3895;
   tmp_3879 += (-2*vd*ZH(gt1,0)) * tmp_3894;
   std::complex<double> tmp_3900;
   std::complex<double> tmp_3901;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3902;
      std::complex<double> tmp_3903;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3903 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3902 += tmp_3903;
      tmp_3901 += (Conj(ZE(gt2,j2))) * tmp_3902;
   }
   tmp_3900 += tmp_3901;
   tmp_3879 += (1.4142135623730951*Conj(Mu)*ZH(gt1,1)) * tmp_3900;
   std::complex<double> tmp_3904;
   std::complex<double> tmp_3905;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3906;
      std::complex<double> tmp_3907;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3907 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3906 += tmp_3907;
      tmp_3905 += (ZE(gt3,j2)) * tmp_3906;
   }
   tmp_3904 += tmp_3905;
   tmp_3879 += (1.4142135623730951*Mu*ZH(gt1,1)) * tmp_3904;
   std::complex<double> tmp_3908;
   std::complex<double> tmp_3909;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3909 += Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1);
   }
   tmp_3908 += tmp_3909;
   tmp_3879 += (0.6*Sqr(g1)*(vd*ZH(gt1,0) - vu*ZH(gt1,1))) * tmp_3908;
   tmp_3876 += (2) * tmp_3879;
   result += (0.25) * tmp_3876;

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

   std::complex<double> tmp_3910;
   std::complex<double> tmp_3911;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3912;
      std::complex<double> tmp_3913;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3913 += Conj(ZER(gt3,j1))*Ye(j1,j2);
      }
      tmp_3912 += tmp_3913;
      tmp_3911 += (Conj(ZEL(gt1,j2))) * tmp_3912;
   }
   tmp_3910 += tmp_3911;
   result += (-0.7071067811865475*ZH(gt2,0)) * tmp_3910;

   return result;
}

std::complex<double> MSSMRHN_effective_couplings::CpFdhhbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_3914;
   std::complex<double> tmp_3915;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3916;
      std::complex<double> tmp_3917;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3917 += Conj(ZDR(gt3,j1))*Yd(j1,j2);
      }
      tmp_3916 += tmp_3917;
      tmp_3915 += (Conj(ZDL(gt1,j2))) * tmp_3916;
   }
   tmp_3914 += tmp_3915;
   result += (-0.7071067811865475*ZH(gt2,0)) * tmp_3914;

   return result;
}

std::complex<double> MSSMRHN_effective_couplings::CpFuhhbarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_3918;
   std::complex<double> tmp_3919;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3920;
      std::complex<double> tmp_3921;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3921 += Conj(ZUR(gt3,j1))*Yu(j1,j2);
      }
      tmp_3920 += tmp_3921;
      tmp_3919 += (Conj(ZUL(gt1,j2))) * tmp_3920;
   }
   tmp_3918 += tmp_3919;
   result += (-0.7071067811865475*ZH(gt2,1)) * tmp_3918;

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

   std::complex<double> tmp_3922;
   std::complex<double> tmp_3923;
   std::complex<double> tmp_3924;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3925;
      std::complex<double> tmp_3926;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3926 += ZD(gt3,3 + j1)*TYd(j1,j2);
      }
      tmp_3925 += tmp_3926;
      tmp_3924 += (Conj(ZD(gt2,j2))) * tmp_3925;
   }
   tmp_3923 += tmp_3924;
   tmp_3922 += (ZA(gt1,0)) * tmp_3923;
   std::complex<double> tmp_3927;
   std::complex<double> tmp_3928;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3929;
      std::complex<double> tmp_3930;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3930 += Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_3929 += tmp_3930;
      tmp_3928 += (ZD(gt3,j2)) * tmp_3929;
   }
   tmp_3927 += tmp_3928;
   tmp_3922 += (-ZA(gt1,0)) * tmp_3927;
   std::complex<double> tmp_3931;
   std::complex<double> tmp_3932;
   std::complex<double> tmp_3933;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3934;
      std::complex<double> tmp_3935;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3935 += Yd(j1,j2)*ZD(gt3,3 + j1);
      }
      tmp_3934 += tmp_3935;
      tmp_3933 += (Conj(ZD(gt2,j2))) * tmp_3934;
   }
   tmp_3932 += tmp_3933;
   tmp_3931 += (Conj(Mu)) * tmp_3932;
   std::complex<double> tmp_3936;
   std::complex<double> tmp_3937;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3938;
      std::complex<double> tmp_3939;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3939 += Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1));
      }
      tmp_3938 += tmp_3939;
      tmp_3937 += (ZD(gt3,j2)) * tmp_3938;
   }
   tmp_3936 += tmp_3937;
   tmp_3931 += (-Mu) * tmp_3936;
   tmp_3922 += (ZA(gt1,1)) * tmp_3931;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_3922;

   return result;
}

std::complex<double> MSSMRHN_effective_couplings::CpAhSuconjSu(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYu = MODELPARAMETER(TYu);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3940;
   std::complex<double> tmp_3941;
   std::complex<double> tmp_3942;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3943;
      std::complex<double> tmp_3944;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3944 += Yu(j1,j2)*ZU(gt3,3 + j1);
      }
      tmp_3943 += tmp_3944;
      tmp_3942 += (Conj(ZU(gt2,j2))) * tmp_3943;
   }
   tmp_3941 += tmp_3942;
   tmp_3940 += (Conj(Mu)*ZA(gt1,0)) * tmp_3941;
   std::complex<double> tmp_3945;
   std::complex<double> tmp_3946;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3947;
      std::complex<double> tmp_3948;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3948 += Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1));
      }
      tmp_3947 += tmp_3948;
      tmp_3946 += (ZU(gt3,j2)) * tmp_3947;
   }
   tmp_3945 += tmp_3946;
   tmp_3940 += (-(Mu*ZA(gt1,0))) * tmp_3945;
   std::complex<double> tmp_3949;
   std::complex<double> tmp_3950;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3951;
      std::complex<double> tmp_3952;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3952 += ZU(gt3,3 + j1)*TYu(j1,j2);
      }
      tmp_3951 += tmp_3952;
      tmp_3950 += (Conj(ZU(gt2,j2))) * tmp_3951;
   }
   tmp_3949 += tmp_3950;
   std::complex<double> tmp_3953;
   std::complex<double> tmp_3954;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3955;
      std::complex<double> tmp_3956;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3956 += Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_3955 += tmp_3956;
      tmp_3954 += (ZU(gt3,j2)) * tmp_3955;
   }
   tmp_3953 += tmp_3954;
   tmp_3949 += (-1) * tmp_3953;
   tmp_3940 += (ZA(gt1,1)) * tmp_3949;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_3940;

   return result;
}

std::complex<double> MSSMRHN_effective_couplings::CpAhSeconjSe(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto TYe = MODELPARAMETER(TYe);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Mu = MODELPARAMETER(Mu);

   std::complex<double> result;

   std::complex<double> tmp_3957;
   std::complex<double> tmp_3958;
   std::complex<double> tmp_3959;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3960;
      std::complex<double> tmp_3961;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3961 += ZE(gt3,3 + j1)*TYe(j1,j2);
      }
      tmp_3960 += tmp_3961;
      tmp_3959 += (Conj(ZE(gt2,j2))) * tmp_3960;
   }
   tmp_3958 += tmp_3959;
   tmp_3957 += (ZA(gt1,0)) * tmp_3958;
   std::complex<double> tmp_3962;
   std::complex<double> tmp_3963;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3964;
      std::complex<double> tmp_3965;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3965 += Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_3964 += tmp_3965;
      tmp_3963 += (ZE(gt3,j2)) * tmp_3964;
   }
   tmp_3962 += tmp_3963;
   tmp_3957 += (-ZA(gt1,0)) * tmp_3962;
   std::complex<double> tmp_3966;
   std::complex<double> tmp_3967;
   std::complex<double> tmp_3968;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3969;
      std::complex<double> tmp_3970;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3970 += Ye(j1,j2)*ZE(gt3,3 + j1);
      }
      tmp_3969 += tmp_3970;
      tmp_3968 += (Conj(ZE(gt2,j2))) * tmp_3969;
   }
   tmp_3967 += tmp_3968;
   tmp_3966 += (Conj(Mu)) * tmp_3967;
   std::complex<double> tmp_3971;
   std::complex<double> tmp_3972;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3973;
      std::complex<double> tmp_3974;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3974 += Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1));
      }
      tmp_3973 += tmp_3974;
      tmp_3972 += (ZE(gt3,j2)) * tmp_3973;
   }
   tmp_3971 += tmp_3972;
   tmp_3966 += (-Mu) * tmp_3971;
   tmp_3957 += (ZA(gt1,1)) * tmp_3966;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_3957;

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

   std::complex<double> tmp_3975;
   std::complex<double> tmp_3976;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3977;
      std::complex<double> tmp_3978;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3978 += Conj(ZER(gt3,j1))*Ye(j1,j2);
      }
      tmp_3977 += tmp_3978;
      tmp_3976 += (Conj(ZEL(gt2,j2))) * tmp_3977;
   }
   tmp_3975 += tmp_3976;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,0)) *
      tmp_3975;

   return result;
}

std::complex<double> MSSMRHN_effective_couplings::CpAhFdbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_3979;
   std::complex<double> tmp_3980;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3981;
      std::complex<double> tmp_3982;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3982 += Conj(ZDR(gt3,j1))*Yd(j1,j2);
      }
      tmp_3981 += tmp_3982;
      tmp_3980 += (Conj(ZDL(gt2,j2))) * tmp_3981;
   }
   tmp_3979 += tmp_3980;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,0)) *
      tmp_3979;

   return result;
}

std::complex<double> MSSMRHN_effective_couplings::CpAhFubarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_3983;
   std::complex<double> tmp_3984;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3985;
      std::complex<double> tmp_3986;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3986 += Conj(ZUR(gt3,j1))*Yu(j1,j2);
      }
      tmp_3985 += tmp_3986;
      tmp_3984 += (Conj(ZUL(gt2,j2))) * tmp_3985;
   }
   tmp_3983 += tmp_3984;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gt1,1)) *
      tmp_3983;

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
