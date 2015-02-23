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

// File generated at Mon 23 Feb 2015 12:25:20

#include "SM_two_scale_low_scale_constraint.hpp"
#include "SM_two_scale_model.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "root_finder.hpp"
#include "weinberg_angle.hpp"

#include <cassert>
#include <cmath>
#include <limits>

namespace flexiblesusy {

#define INPUTPARAMETER(p) inputPars.p
#define MODELPARAMETER(p) model->get_##p()
#define BETAPARAMETER(p) beta_functions.get_##p()
#define BETA(p) beta_##p
#define SM(p) Electroweak_constants::p
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define MODEL model
#define MODELCLASSNAME SM<Two_scale>
#define THETAW theta_w
#define ALPHA_EM_DRBAR alpha_em_drbar
#define CALCULATE_DRBAR_MASSES() model->calculate_DRbar_masses()

SM_low_scale_constraint<Two_scale>::SM_low_scale_constraint()
   : Constraint<Two_scale>()
   , scale(0.)
   , initial_scale_guess(0.)
   , model(0)
   , inputPars()
   , oneset()
   , MWDRbar(0.)
   , MZDRbar(0.)
   , AlphaS(0.)
   , EDRbar(0.)
   , ThetaWDRbar(0.)
   , new_g1(0.)
   , new_g2(0.)
   , new_g3(0.)
   , self_energy_w_at_mw(0.)
   , threshold_corrections_loop_order(1)
{
}

SM_low_scale_constraint<Two_scale>::SM_low_scale_constraint(
   SM<Two_scale>* model_,
   const SM_input_parameters& inputPars_, const QedQcd& oneset_)
   : Constraint<Two_scale>()
   , model(model_)
   , inputPars(inputPars_)
   , oneset(oneset_)
   , new_g1(0.)
   , new_g2(0.)
   , new_g3(0.)
   , self_energy_w_at_mw(0.)
{
   initialize();
}

SM_low_scale_constraint<Two_scale>::~SM_low_scale_constraint()
{
}

void SM_low_scale_constraint<Two_scale>::apply()
{
   assert(model && "Error: SM_low_scale_constraint::apply():"
          " model pointer must not be zero");

   model->calculate_DRbar_masses();
   update_scale();
   calculate_DRbar_gauge_couplings();

   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);

   MODEL->set_v((2*MZDRbar)/Sqrt(0.6*Sqr(g1) + Sqr(g2)));
   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();


   model->set_g1(new_g1);
   model->set_g2(new_g2);
   model->set_g3(new_g3);

   recalculate_mw_pole();
}

double SM_low_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double SM_low_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

void SM_low_scale_constraint<Two_scale>::set_model(Two_scale_model* model_)
{
   model = cast_model<SM<Two_scale>*>(model_);
}

void SM_low_scale_constraint<Two_scale>::set_input_parameters(const SM_input_parameters& inputPars_)
{
   inputPars = inputPars_;
}

void SM_low_scale_constraint<Two_scale>::set_sm_parameters(const QedQcd& oneset_)
{
   oneset = oneset_;
}

const QedQcd& SM_low_scale_constraint<Two_scale>::get_sm_parameters() const
{
   return oneset;
}

void SM_low_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = NULL;
   oneset = QedQcd();
   MWDRbar = 0.;
   MZDRbar = 0.;
   AlphaS = 0.;
   EDRbar = 0.;
   ThetaWDRbar = 0.;
   new_g1 = 0.;
   new_g2 = 0.;
   new_g3 = 0.;
   self_energy_w_at_mw = 0.;
}

void SM_low_scale_constraint<Two_scale>::initialize()
{
   assert(model && "SM_low_scale_constraint<Two_scale>::"
          "initialize(): model pointer is zero.");

   initial_scale_guess = SM(MZ);

   scale = initial_scale_guess;

   MWDRbar = 0.;
   MZDRbar = 0.;
   AlphaS = 0.;
   EDRbar = 0.;
   ThetaWDRbar = 0.;
   new_g1 = 0.;
   new_g2 = 0.;
   new_g3 = 0.;
   self_energy_w_at_mw = 0.;
}

void SM_low_scale_constraint<Two_scale>::update_scale()
{
   assert(model && "SM_low_scale_constraint<Two_scale>::"
          "update_scale(): model pointer is zero.");

   scale = SM(MZ);


}

void SM_low_scale_constraint<Two_scale>::calculate_threshold_corrections()
{
   assert(oneset.displayMu() == get_scale() && "Error: low-energy data"
          " set is not defined at the same scale as the low-energy"
          " constraint.  You need to run the low-energy data set to this"
          " scale!");

   const double alpha_em = oneset.displayAlpha(ALPHA);
   const double alpha_s  = oneset.displayAlpha(ALPHAS);
   const double mw_pole  = oneset.displayPoleMW();
   const double mz_pole  = oneset.displayPoleMZ();

   double delta_alpha_em = 0.;
   double delta_alpha_s  = 0.;

   if (model->get_thresholds()) {
      delta_alpha_em = calculate_delta_alpha_em(alpha_em);
      delta_alpha_s  = calculate_delta_alpha_s(alpha_s);
   }

   const double alpha_em_drbar = alpha_em / (1.0 - delta_alpha_em);
   const double alpha_s_drbar  = alpha_s  / (1.0 - delta_alpha_s);
   const double e_drbar        = Sqrt(4.0 * Pi * alpha_em_drbar);

   // interface variables
   MZDRbar = mz_pole;
   MWDRbar = mw_pole;

   if (model->get_thresholds()) {
      MZDRbar = model->calculate_MVZ_DRbar(mz_pole);
      MWDRbar = model->calculate_MVWp_DRbar(mw_pole);
   }

   AlphaS = alpha_s_drbar;
   EDRbar = e_drbar;
   ThetaWDRbar = calculate_theta_w(alpha_em_drbar);
}

double SM_low_scale_constraint<Two_scale>::calculate_theta_w(double alpha_em_drbar)
{
   double theta_w = 0.;

   using namespace weinberg_angle;

   const double scale         = MODEL->get_scale();
   const double mw_pole       = oneset.displayPoleMW();
   const double mz_pole       = oneset.displayPoleMZ();
   const double mt_pole       = oneset.displayPoleMt();
   const double mt_drbar      = MODEL->get_MFu(2);
   const double mb_drbar      = MODEL->get_MFd(2);
   const double mh_drbar      = MODEL->get_Mhh();
   const double gY            = MODEL->get_g1() * 0.7745966692414834;
   const double g2            = MODEL->get_g2();
   const double g3            = MODEL->get_g3();
   const double ymu           = MODEL->get_Ye(1,1);
   const double pizztMZ       = Re(MODEL->self_energy_VZ(mz_pole));
   const double piwwt0        = Re(MODEL->self_energy_VWp(0.));
   self_energy_w_at_mw        = Re(MODEL->self_energy_VWp(mw_pole));

   Weinberg_angle::Self_energy_data se_data;
   se_data.scale    = scale;
   se_data.mt_pole  = mt_pole;
   se_data.mt_drbar = mt_drbar;
   se_data.mb_drbar = mb_drbar;
   se_data.gY       = gY;
   se_data.g2       = g2;

   const double pizztMZ_corrected =
      Weinberg_angle::replace_mtop_in_self_energy_z(pizztMZ, mz_pole,
         se_data);

   const double piwwtMW_corrected =
      Weinberg_angle::replace_mtop_in_self_energy_w(self_energy_w_at_mw,
         mw_pole, se_data);

   const double piwwt0_corrected =
      Weinberg_angle::replace_mtop_in_self_energy_w(piwwt0, 0., se_data);

   Weinberg_angle::Data data;
   data.scale               = scale;
   data.alpha_em_drbar      = ALPHA_EM_DRBAR;
   data.fermi_contant       = oneset.displayFermiConstant();
   data.self_energy_z_at_mz = pizztMZ_corrected;
   data.self_energy_w_at_mw = piwwtMW_corrected;
   data.self_energy_w_at_0  = piwwt0_corrected;
   data.mw_pole             = mw_pole;
   data.mz_pole             = mz_pole;
   data.mt_pole             = mt_pole;
   data.mh_drbar            = mh_drbar;
   data.gY                  = gY;
   data.g2                  = g2;
   data.g3                  = g3;
   data.ymu                 = ymu;

   Weinberg_angle weinberg;
   weinberg.disable_susy_contributions();
   weinberg.set_number_of_loops(MODEL->get_thresholds());
   weinberg.set_data(data);

   const int error = weinberg.calculate();

   THETAW = ArcSin(weinberg.get_sin_theta());

   if (error)
      MODEL->get_problems().flag_no_rho_convergence();
   else
      MODEL->get_problems().unflag_no_rho_convergence();


   return theta_w;
}

void SM_low_scale_constraint<Two_scale>::calculate_DRbar_gauge_couplings()
{
   calculate_threshold_corrections();

   new_g1 = 1.2909944487358056*EDRbar*Sec(ThetaWDRbar);
   new_g2 = EDRbar*Csc(ThetaWDRbar);
   new_g3 = 3.5449077018110318*Sqrt(AlphaS);

}

double SM_low_scale_constraint<Two_scale>::calculate_delta_alpha_em(double alphaEm) const
{
   const double currentScale = model->get_scale();
   const auto MFu = MODELPARAMETER(MFu);

   const double delta_alpha_em_SM = 0.15915494309189535*alphaEm*(
      0.3333333333333333 - 1.7777777777777777*FiniteLog(Abs(MFu(2)/currentScale)))
      ;

   const double delta_alpha_em = 0;

   return delta_alpha_em + delta_alpha_em_SM;

}

double SM_low_scale_constraint<Two_scale>::calculate_delta_alpha_s(double alphaS) const
{
   const double currentScale = model->get_scale();
   const auto MFu = MODELPARAMETER(MFu);

   const double delta_alpha_s_SM = -0.1061032953945969*alphaS*FiniteLog(Abs(MFu
      (2)/currentScale));

   const double delta_alpha_s = 0;

   return delta_alpha_s + delta_alpha_s_SM;

}

void SM_low_scale_constraint<Two_scale>::calculate_DRbar_yukawa_couplings()
{
   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();
}

void SM_low_scale_constraint<Two_scale>::calculate_Yu_DRbar()
{
   Eigen::Matrix<double,3,3> topDRbar(Eigen::Matrix<double,3,3>::Zero());
   topDRbar(0,0)      = oneset.displayMass(mUp);
   topDRbar(1,1)      = oneset.displayMass(mCharm);
   topDRbar(2,2)      = oneset.displayMass(mTop);

   if (model->get_thresholds())
      topDRbar(2,2) = model->calculate_MFu_DRbar(oneset.displayPoleMt(), 2);

   const auto v = MODELPARAMETER(v);
   MODEL->set_Yu(-((1.4142135623730951*topDRbar)/v).transpose());

}

void SM_low_scale_constraint<Two_scale>::calculate_Yd_DRbar()
{
   Eigen::Matrix<double,3,3> bottomDRbar(Eigen::Matrix<double,3,3>::Zero());
   bottomDRbar(0,0)   = oneset.displayMass(mDown);
   bottomDRbar(1,1)   = oneset.displayMass(mStrange);
   bottomDRbar(2,2)   = oneset.displayMass(mBottom);

   if (model->get_thresholds())
      bottomDRbar(2,2) = model->calculate_MFd_DRbar(oneset.displayMass(mBottom), 2);

   const auto v = MODELPARAMETER(v);
   MODEL->set_Yd(((1.4142135623730951*bottomDRbar)/v).transpose());

}

void SM_low_scale_constraint<Two_scale>::calculate_Ye_DRbar()
{
   Eigen::Matrix<double,3,3> electronDRbar(Eigen::Matrix<double,3,3>::Zero());
   electronDRbar(0,0) = oneset.displayMass(mElectron);
   electronDRbar(1,1) = oneset.displayMass(mMuon);
   electronDRbar(2,2) = oneset.displayMass(mTau);

   if (model->get_thresholds())
      electronDRbar(2,2) = model->calculate_MFe_DRbar(oneset.displayMass(mTau), 2);

   const auto v = MODELPARAMETER(v);
   MODEL->set_Ye(((1.4142135623730951*electronDRbar)/v).transpose());

}

/**
 * Recalculates the W boson pole mass using the new gauge couplings.
 */
void SM_low_scale_constraint<Two_scale>::recalculate_mw_pole()
{
   if (!model->get_thresholds())
      return;

   MODEL->calculate_MVWp();

   const double mw_drbar    = MODEL->get_MVWp();
   const double mw_pole_sqr = Sqr(mw_drbar) - self_energy_w_at_mw;

   if (mw_pole_sqr < 0.)
      MODEL->get_problems().flag_tachyon(SM_info::VWp);

   const double mw_pole = AbsSqrt(mw_pole_sqr);

   oneset.setPoleMW(mw_pole);

}

} // namespace flexiblesusy
