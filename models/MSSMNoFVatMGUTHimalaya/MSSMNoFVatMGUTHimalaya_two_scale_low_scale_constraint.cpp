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


#include "MSSMNoFVatMGUTHimalaya_two_scale_low_scale_constraint.hpp"
#include "MSSMNoFVatMGUTHimalaya_two_scale_model.hpp"
#include "MSSMNoFVatMGUTHimalaya_info.hpp"
#include "MSSMNoFVatMGUTHimalaya_weinberg_angle.hpp"
#include "config.h"
#include "wrappers.hpp"
#include "logger.hpp"
#include "error.hpp"
#include "ew_input.hpp"
#include "minimizer.hpp"
#include "raii.hpp"
#include "root_finder.hpp"
#include "threshold_loop_functions.hpp"
#include "numerics2.hpp"
#include "mssm_twoloop_mb.hpp"
#include "mssm_twoloop_mt.hpp"
#include "mssm_twoloop_mtau.hpp"
#include "mssm_twoloop_as.hpp"


#include <algorithm>
#include <cmath>
#include <limits>

namespace flexiblesusy {

#define DERIVEDPARAMETER(p) model->p()
#define EXTRAPARAMETER(p) model->get_##p()
#define INPUTPARAMETER(p) model->get_input().p
#define MODELPARAMETER(p) model->get_##p()
#define PHASE(p) model->get_##p()
#define BETAPARAMETER(p) beta_functions.get_##p()
#define BETAPARAMETER1(l,p) beta_functions_##l##L.get_##p()
#define BETA(p) beta_##p
#define BETA1(l,p) beta_##l##L_##p
#define LowEnergyGaugeCoupling(i) new_g##i
#define LowEnergyConstant(p) Electroweak_constants::p
#define MZPole qedqcd.displayPoleMZ()
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define SCALE model->get_scale()
#define MODEL model
#define MODELCLASSNAME MSSMNoFVatMGUTHimalaya<Two_scale>
#define MWMSbar mW_run
#define MWDRbar mW_run
#define MZMSbar mZ_run
#define MZDRbar mZ_run
#define EDRbar e_run
#define EMSbar e_run
#define CKM ckm
#define PMNS pmns
#define THETAW theta_w
#define THRESHOLD static_cast<int>(model->get_thresholds())
#define ALPHA_EM_DRBAR alpha_em_drbar
#define CALCULATE_DRBAR_MASSES() model->calculate_DRbar_masses()

MSSMNoFVatMGUTHimalaya_low_scale_constraint<Two_scale>::MSSMNoFVatMGUTHimalaya_low_scale_constraint(
   MSSMNoFVatMGUTHimalaya<Two_scale>* model_, const softsusy::QedQcd& qedqcd_)
   : model(model_)
   , qedqcd(qedqcd_)
{
   initialize();
}

void MSSMNoFVatMGUTHimalaya_low_scale_constraint<Two_scale>::apply()
{
   check_model_ptr();

   


   model->calculate_DRbar_masses();
   update_scale();
   qedqcd.run_to(scale, 1.0e-5);
   calculate_DRbar_gauge_couplings();
   calculate_running_SM_masses();

   const auto TanBeta = INPUTPARAMETER(TanBeta);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);

   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();
   MODEL->set_vd(Re((2*MZDRbar)/(Sqrt(0.6*Sqr(g1) + Sqr(g2))*Sqrt(1 + Sqr(TanBeta)
      ))));
   MODEL->set_vu(Re((2*MZDRbar*TanBeta)/(Sqrt(0.6*Sqr(g1) + Sqr(g2))*Sqrt(1 + Sqr(
      TanBeta)))));
   MODEL->set_g1(new_g1);
   MODEL->set_g2(new_g2);
   MODEL->set_g3(new_g3);

}

const Eigen::Matrix<std::complex<double>,3,3>& MSSMNoFVatMGUTHimalaya_low_scale_constraint<Two_scale>::get_ckm()
{
   return ckm;
}

const Eigen::Matrix<std::complex<double>,3,3>& MSSMNoFVatMGUTHimalaya_low_scale_constraint<Two_scale>::get_pmns()
{
   return pmns;
}

double MSSMNoFVatMGUTHimalaya_low_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double MSSMNoFVatMGUTHimalaya_low_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

void MSSMNoFVatMGUTHimalaya_low_scale_constraint<Two_scale>::set_model(Model* model_)
{
   model = cast_model<MSSMNoFVatMGUTHimalaya<Two_scale>*>(model_);
}

void MSSMNoFVatMGUTHimalaya_low_scale_constraint<Two_scale>::set_sm_parameters(
   const softsusy::QedQcd& qedqcd_)
{
   qedqcd = qedqcd_;
}

const softsusy::QedQcd& MSSMNoFVatMGUTHimalaya_low_scale_constraint<Two_scale>::get_sm_parameters() const
{
   return qedqcd;
}

void MSSMNoFVatMGUTHimalaya_low_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = nullptr;
   qedqcd = softsusy::QedQcd();
   ckm.setIdentity();
   pmns.setIdentity();
   upQuarksDRbar.setZero();
   downQuarksDRbar.setZero();
   downLeptonsDRbar.setZero();
   neutrinoDRbar.setZero();
   mW_run = 0.;
   mZ_run = 0.;
   AlphaS = 0.;
   e_run = 0.;
   ThetaWDRbar = 0.;
   new_g1 = 0.;
   new_g2 = 0.;
   new_g3 = 0.;
}

void MSSMNoFVatMGUTHimalaya_low_scale_constraint<Two_scale>::initialize()
{
   check_model_ptr();

   const auto Mlow = INPUTPARAMETER(Mlow);

   initial_scale_guess = IF(Mlow != 0, Mlow, LowEnergyConstant(MZ));

   scale = initial_scale_guess;

   ckm = qedqcd.get_complex_ckm();
   pmns = qedqcd.get_complex_pmns();
   upQuarksDRbar.setZero();
   downQuarksDRbar.setZero();
   downLeptonsDRbar.setZero();
   neutrinoDRbar.setZero();
   mW_run = 0.;
   mZ_run = 0.;
   AlphaS = 0.;
   e_run = 0.;
   ThetaWDRbar = 0.;
   new_g1 = 0.;
   new_g2 = 0.;
   new_g3 = 0.;
}

void MSSMNoFVatMGUTHimalaya_low_scale_constraint<Two_scale>::update_scale()
{
   check_model_ptr();

   const auto Mlow = INPUTPARAMETER(Mlow);

   scale = IF(Mlow != 0, Mlow, LowEnergyConstant(MZ));


}

void MSSMNoFVatMGUTHimalaya_low_scale_constraint<Two_scale>::calculate_threshold_corrections()
{
   check_model_ptr();

   if (!is_zero(qedqcd.get_scale() - get_scale()))
      throw SetupError("Error: low-energy data"
          " set is not defined at the same scale as the low-energy"
          " constraint.  You need to run the low-energy data set to this"
          " scale!");

   const double alpha_em = qedqcd.displayAlpha(softsusy::ALPHA);
   const double alpha_s  = qedqcd.displayAlpha(softsusy::ALPHAS);
   const double mw_pole  = qedqcd.displayPoleMW();
   const double mz_pole  = qedqcd.displayPoleMZ();

   double delta_alpha_em = 0.;
   double delta_alpha_s  = 0.;

   if (model->get_thresholds() && model->get_threshold_corrections().alpha_em > 0)
      delta_alpha_em = calculate_delta_alpha_em(alpha_em);

   if (model->get_thresholds() && model->get_threshold_corrections().alpha_s > 0)
      delta_alpha_s  = calculate_delta_alpha_s(alpha_s);

   const double alpha_em_drbar = alpha_em / (1.0 - delta_alpha_em);
   const double alpha_s_drbar  = alpha_s  / (1.0 - delta_alpha_s);
   const double e_drbar        = Sqrt(4.0 * Pi * alpha_em_drbar);

   // interface variables
   mZ_run = mz_pole;
   mW_run = mw_pole;

   if (model->get_thresholds() && model->get_threshold_corrections().mz > 0)
      mZ_run = model->calculate_MVZ_DRbar(mz_pole);

   if (model->get_thresholds() && model->get_threshold_corrections().mw > 0)
      mW_run = model->calculate_MVWm_DRbar(mw_pole);

   AlphaS = alpha_s_drbar;
   e_run = e_drbar;
   ThetaWDRbar = calculate_theta_w();
}

double MSSMNoFVatMGUTHimalaya_low_scale_constraint<Two_scale>::calculate_theta_w()
{
   check_model_ptr();

   double theta_w = std::asin(Electroweak_constants::sinThetaW);

   const auto get_mh_pole = [&] () {
      double mh_pole = MODEL->get_physical().Mhh(higgs_idx);
      if (mh_pole == 0) {
         mh_pole = Electroweak_constants::MH;
      }
      return mh_pole;
   };

   MSSMNoFVatMGUTHimalaya_weinberg_angle::Sm_parameters sm_pars;
   sm_pars.fermi_constant = qedqcd.displayFermiConstant();
   sm_pars.mw_pole = qedqcd.displayPoleMW();
   sm_pars.mz_pole = qedqcd.displayPoleMZ();
   sm_pars.mt_pole = qedqcd.displayPoleMt();
   sm_pars.mh_pole = get_mh_pole();
   sm_pars.alpha_s = calculate_alpha_s_SM5_at(qedqcd, qedqcd.displayPoleMt());
   sm_pars.alpha_s_mz = qedqcd.displayAlphaSInput();
   sm_pars.dalpha_s_5_had = Electroweak_constants::delta_alpha_s_5_had;
   sm_pars.higgs_index = higgs_idx;

   const int number_of_iterations =
       std::max(20, static_cast<int>(std::abs(-log10(MODEL->get_precision()) * 10)
          ));

   MSSMNoFVatMGUTHimalaya_weinberg_angle weinberg(MODEL, sm_pars);
   weinberg.set_number_of_loops(MODEL->get_threshold_corrections().sin_theta_w);
   weinberg.set_number_of_iterations(number_of_iterations);

   try {
      const auto result = weinberg.calculate();
      THETAW = ArcSin(result.first);

      if (MODEL->get_thresholds() && MODEL->get_threshold_corrections().
         sin_theta_w > 0)
         Pole(MVWm) = result.second;

      MODEL->get_problems().unflag_no_sinThetaW_convergence();
   } catch (const Error& e) {
      VERBOSE_MSG(e.what_detailed());
      MODEL->get_problems().flag_no_sinThetaW_convergence();
   }

   return theta_w;
}

void MSSMNoFVatMGUTHimalaya_low_scale_constraint<Two_scale>::calculate_DRbar_gauge_couplings()
{
   check_model_ptr();
   calculate_threshold_corrections();

   new_g1 = 1.2909944487358056*EDRbar*Sec(ThetaWDRbar);
   new_g2 = EDRbar*Csc(ThetaWDRbar);
   new_g3 = 3.5449077018110318*Sqrt(AlphaS);

   if (IsFinite(new_g1)) {
      model->get_problems().unflag_non_perturbative_parameter(
         MSSMNoFVatMGUTHimalaya_info::g1);
   } else {
      model->get_problems().flag_non_perturbative_parameter(
         MSSMNoFVatMGUTHimalaya_info::g1, new_g1, get_scale());
      new_g1 = Electroweak_constants::g1;
   }

   if (IsFinite(new_g2)) {
      model->get_problems().unflag_non_perturbative_parameter(
         MSSMNoFVatMGUTHimalaya_info::g2);
   } else {
      model->get_problems().flag_non_perturbative_parameter(
         MSSMNoFVatMGUTHimalaya_info::g2, new_g2, get_scale());
      new_g2 = Electroweak_constants::g2;
   }
}

double MSSMNoFVatMGUTHimalaya_low_scale_constraint<Two_scale>::calculate_delta_alpha_em(double alphaEm) const
{
   check_model_ptr();

   const double currentScale = model->get_scale();
   const auto MCha = MODELPARAMETER(MCha);
   const auto MHpm = MODELPARAMETER(MHpm);
   const auto MSb = MODELPARAMETER(MSb);
   const auto MSc = MODELPARAMETER(MSc);
   const auto MSd = MODELPARAMETER(MSd);
   const auto MSe = MODELPARAMETER(MSe);
   const auto MSm = MODELPARAMETER(MSm);
   const auto MSs = MODELPARAMETER(MSs);
   const auto MSt = MODELPARAMETER(MSt);
   const auto MStau = MODELPARAMETER(MStau);
   const auto MSu = MODELPARAMETER(MSu);
   const auto MFt = MODELPARAMETER(MFt);

   const double delta_alpha_em_SM = -0.28294212105225836*alphaEm*FiniteLog(Abs(MFt
      /currentScale));

   const double delta_alpha_em = 0.15915494309189535*alphaEm*(0.3333333333333333 -
      1.3333333333333333*FiniteLog(Abs(MCha(0)/currentScale)) - 1.3333333333333333
      *FiniteLog(Abs(MCha(1)/currentScale)) - 0.3333333333333333*FiniteLog(Abs(
      MHpm(1)/currentScale)) - 0.1111111111111111*FiniteLog(Abs(MSb(0)/
      currentScale)) - 0.1111111111111111*FiniteLog(Abs(MSb(1)/currentScale)) -
      0.4444444444444444*FiniteLog(Abs(MSc(0)/currentScale)) - 0.4444444444444444*
      FiniteLog(Abs(MSc(1)/currentScale)) - 0.1111111111111111*FiniteLog(Abs(MSd(0
      )/currentScale)) - 0.1111111111111111*FiniteLog(Abs(MSd(1)/currentScale)) -
      0.3333333333333333*FiniteLog(Abs(MSe(0)/currentScale)) - 0.3333333333333333*
      FiniteLog(Abs(MSe(1)/currentScale)) - 0.3333333333333333*FiniteLog(Abs(MSm(0
      )/currentScale)) - 0.3333333333333333*FiniteLog(Abs(MSm(1)/currentScale)) -
      0.1111111111111111*FiniteLog(Abs(MSs(0)/currentScale)) - 0.1111111111111111*
      FiniteLog(Abs(MSs(1)/currentScale)) - 0.4444444444444444*FiniteLog(Abs(MSt(0
      )/currentScale)) - 0.4444444444444444*FiniteLog(Abs(MSt(1)/currentScale)) -
      0.3333333333333333*FiniteLog(Abs(MStau(0)/currentScale)) -
      0.3333333333333333*FiniteLog(Abs(MStau(1)/currentScale)) -
      0.4444444444444444*FiniteLog(Abs(MSu(0)/currentScale)) - 0.4444444444444444*
      FiniteLog(Abs(MSu(1)/currentScale)));

   return delta_alpha_em + delta_alpha_em_SM;

}

double MSSMNoFVatMGUTHimalaya_low_scale_constraint<Two_scale>::calculate_delta_alpha_s(double alphaS) const
{
   check_model_ptr();

   const double currentScale = model->get_scale();
   const auto MSb = MODELPARAMETER(MSb);
   const auto MSc = MODELPARAMETER(MSc);
   const auto MSd = MODELPARAMETER(MSd);
   const auto MSs = MODELPARAMETER(MSs);
   const auto MSt = MODELPARAMETER(MSt);
   const auto MSu = MODELPARAMETER(MSu);
   const auto MGlu = MODELPARAMETER(MGlu);
   const auto MFt = MODELPARAMETER(MFt);

   const double delta_alpha_s_SM = -0.1061032953945969*alphaS*FiniteLog(Abs(MFt/
      currentScale));

   const double delta_alpha_s = 0.15915494309189535*alphaS*(0.5 - 2*FiniteLog(Abs(
      MGlu/currentScale)) - 0.16666666666666666*FiniteLog(Abs(MSb(0)/currentScale)
      ) - 0.16666666666666666*FiniteLog(Abs(MSb(1)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSc(0)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSc(1)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSd(0)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSd(1)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSs(0)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSs(1)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSt(0)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSt(1)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSu(0)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSu(1)/currentScale)));

   const double delta_alpha_s_1loop = delta_alpha_s + delta_alpha_s_SM;
   double delta_alpha_s_2loop = 0.;
   double delta_alpha_s_3loop = 0.;
   double delta_alpha_s_4loop = 0.;

   if (model->get_thresholds() > 1 && model->get_threshold_corrections().alpha_s >
      1) {
      const auto Mu = MODELPARAMETER(Mu);


      double mst_1, mst_2, theta_t;
      double msb_1, msb_2, theta_b;
      double msd_1, msd_2, theta_d;

      model->calculate_MTopSquark_3rd_generation(mst_1, mst_2, theta_t);
      model->calculate_MBottomSquark_3rd_generation(msb_1, msb_2, theta_b);
      model->calculate_MBottomSquark_2nd_generation(msd_1, msd_2, theta_d);

      mssm_twoloop_as::Parameters pars;
      pars.g3   = model->get_g3();
      pars.yt   = model->get_Yu(2,2);
      pars.yb   = model->get_Yd(2,2);
      pars.mt   = model->get_MFt();
      pars.mb   = model->get_MFb();
      pars.mg   = model->get_MGlu();
      pars.mst1 = mst_1;
      pars.mst2 = mst_2;
      pars.msb1 = msb_1;
      pars.msb2 = msb_2;
      pars.msd1 = msd_1;
      pars.msd2 = msd_2;
      pars.xt   = Sin(2*theta_t) * (Sqr(mst_1) - Sqr(mst_2)) / (2. * pars.mt);
      pars.xb   = Sin(2*theta_b) * (Sqr(msb_1) - Sqr(msb_2)) / (2. * pars.mb);
      pars.mw   = model->get_MVWm();
      pars.mz   = model->get_MVZ();
      pars.mh   = model->get_Mhh(0);
      pars.mH   = model->get_Mhh(1);
      pars.mC   = model->get_MHpm(1);
      pars.mA   = model->get_MAh(1);
      pars.mu   = Mu;
      pars.tb   = model->get_vu() / model->get_vd();
      pars.Q    = model->get_scale();

      delta_alpha_s_2loop =
         - Sqr(delta_alpha_s_1loop)/4.
         - 2.*(
            + mssm_twoloop_as::delta_alpha_s_2loop_as_as(pars)
            + mssm_twoloop_as::delta_alpha_s_2loop_at_as(pars)
            + mssm_twoloop_as::delta_alpha_s_2loop_ab_as(pars)
           );
   }

   return delta_alpha_s_1loop + delta_alpha_s_2loop + delta_alpha_s_3loop +
      delta_alpha_s_4loop;

}

double MSSMNoFVatMGUTHimalaya_low_scale_constraint<Two_scale>::calculate_alpha_s_SM5_at(
   softsusy::QedQcd qedqcd_tmp, double scale) const
{
   qedqcd_tmp.run_to(scale); // running in SM(5)
   return qedqcd_tmp.displayAlpha(softsusy::ALPHAS);
}

void MSSMNoFVatMGUTHimalaya_low_scale_constraint<Two_scale>::calculate_running_SM_masses()
{
   check_model_ptr();

   upQuarksDRbar.setZero();
   upQuarksDRbar(0,0) = qedqcd.displayMass(softsusy::mUp);
   upQuarksDRbar(1,1) = qedqcd.displayMass(softsusy::mCharm);
   upQuarksDRbar(2,2) = qedqcd.displayPoleMt();

   downQuarksDRbar.setZero();
   downQuarksDRbar(0,0) = qedqcd.displayMass(softsusy::mDown);
   downQuarksDRbar(1,1) = qedqcd.displayMass(softsusy::mStrange);
   downQuarksDRbar(2,2) = qedqcd.displayMass(softsusy::mBottom);

   downLeptonsDRbar.setZero();
   downLeptonsDRbar(0,0) = qedqcd.displayPoleMel();
   downLeptonsDRbar(1,1) = qedqcd.displayPoleMmuon();
   downLeptonsDRbar(2,2) = qedqcd.displayPoleMtau();

   neutrinoDRbar.setZero();
   neutrinoDRbar(0,0) = qedqcd.displayNeutrinoPoleMass(1);
   neutrinoDRbar(1,1) = qedqcd.displayNeutrinoPoleMass(2);
   neutrinoDRbar(2,2) = qedqcd.displayNeutrinoPoleMass(3);

   if (model->get_thresholds() && model->get_threshold_corrections().mt > 0) {
      upQuarksDRbar(2,2) = MODEL->calculate_MFt_DRbar(qedqcd.displayPoleMt());
   }

   if (model->get_thresholds() && model->get_threshold_corrections().mb > 0) {
      downQuarksDRbar(2,2) = MODEL->calculate_MFb_DRbar(qedqcd.displayMass(softsusy::mBottom));
   }

   if (model->get_thresholds()) {
      downLeptonsDRbar(0,0) = MODEL->calculate_MFe_DRbar(qedqcd.displayMass(softsusy::mElectron));
      downLeptonsDRbar(1,1) = MODEL->calculate_MFm_DRbar(qedqcd.displayMass(softsusy::mMuon));
   }

   if (model->get_thresholds() && model->get_threshold_corrections().mtau > 0) {
      downLeptonsDRbar(2,2) = MODEL->calculate_MFtau_DRbar(qedqcd.displayMass(softsusy::mTau));
   }
}

void MSSMNoFVatMGUTHimalaya_low_scale_constraint<Two_scale>::calculate_DRbar_yukawa_couplings()
{
   calculate_running_SM_masses();
   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();
}

void MSSMNoFVatMGUTHimalaya_low_scale_constraint<Two_scale>::calculate_Yu_DRbar()
{
   check_model_ptr();

   const auto vu = MODELPARAMETER(vu);
   MODEL->set_Yu(((1.4142135623730951*upQuarksDRbar)/vu).real());

}

void MSSMNoFVatMGUTHimalaya_low_scale_constraint<Two_scale>::calculate_Yd_DRbar()
{
   check_model_ptr();

   const auto vd = MODELPARAMETER(vd);
   MODEL->set_Yd(((1.4142135623730951*downQuarksDRbar)/vd).real());

}

void MSSMNoFVatMGUTHimalaya_low_scale_constraint<Two_scale>::calculate_Ye_DRbar()
{
   check_model_ptr();

   const auto vd = MODELPARAMETER(vd);
   MODEL->set_Ye(((1.4142135623730951*downLeptonsDRbar)/vd).real());

}

void MSSMNoFVatMGUTHimalaya_low_scale_constraint<Two_scale>::check_model_ptr() const
{
   if (!model)
      throw SetupError("MSSMNoFVatMGUTHimalaya_low_scale_constraint<Two_scale>: "
                       "model pointer is zero!");
}

} // namespace flexiblesusy
