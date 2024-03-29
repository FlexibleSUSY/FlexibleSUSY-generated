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


#include "CNMSSM_semi_analytic_initial_guesser.hpp"
#include "CNMSSM_semi_analytic_model.hpp"
#include "CNMSSM_semi_analytic_susy_convergence_tester.hpp"
#include "lowe.h"
#include "error.hpp"
#include "ew_input.hpp"
#include "raii.hpp"
#include "wrappers.hpp"
#include "two_scale_running_precision.hpp"
#include "two_scale_solver.hpp"

#include <Eigen/Core>

namespace flexiblesusy {

#define DERIVEDPARAMETER(p) model->p()
#define INPUTPARAMETER(p) model->get_input().p
#define MODELPARAMETER(p) model->get_##p()
#define PHASE(p) model->get_##p()
#define LowEnergyConstant(p) Electroweak_constants::p
#define MODEL model

CNMSSM_initial_guesser<Semi_analytic>::CNMSSM_initial_guesser(
   CNMSSM<Semi_analytic>* model_,
   const softsusy::QedQcd& qedqcd_,
   CNMSSM_low_scale_constraint<Semi_analytic>& low_constraint_,
   CNMSSM_susy_scale_constraint<Semi_analytic>& susy_constraint_,
   CNMSSM_high_scale_constraint<Semi_analytic>& high_constraint_
)
   : model(model_)
   , qedqcd(qedqcd_)
   , low_constraint(low_constraint_)
   , susy_constraint(susy_constraint_)
   , high_constraint(high_constraint_)
{
   if (!model)
      throw SetupError("CNMSSM_initial_guesser: Error: pointer to model"
                       " CNMSSM<Semi_analytic> must not be zero");
}

/**
 * Guesses the DR-bar model parameters by calling
 * guess_susy_parameters() and guess_soft_parameters() .
 */
void CNMSSM_initial_guesser<Semi_analytic>::guess()
{
   guess_susy_parameters();
   guess_soft_parameters();
}

/**
 * Guesses the SUSY parameters (gauge, Yukawa couplings) at
 * \f$m_\text{top}^\text{pole}\f$ from the Standard Model gauge
 * couplings and fermion masses.  Threshold corrections are ignored.
 * The user-defined initial guess at the low-scale and high-scale
 * are applied, before a two-scale iteration is used to obtain a full
 * initial guess for the SUSY parameters.
 */
void CNMSSM_initial_guesser<Semi_analytic>::guess_susy_parameters()
{
   initial_guess_low_scale_parameters();
   initial_guess_high_scale_parameters();
   solve_susy_parameters();
}

void CNMSSM_initial_guesser<Semi_analytic>::initial_guess_low_scale_parameters()
{
   softsusy::QedQcd leAtMt(qedqcd);
   const double mtpole = leAtMt.displayPoleMt();

   mu_guess = leAtMt.displayMass(softsusy::mUp);
   mc_guess = leAtMt.displayMass(softsusy::mCharm);
   mt_guess = model->get_thresholds() > 0 && model->get_threshold_corrections().mt > 0 ?
      leAtMt.displayMass(softsusy::mTop) - 30.0 :
      leAtMt.displayPoleMt();
   md_guess = leAtMt.displayMass(softsusy::mDown);
   ms_guess = leAtMt.displayMass(softsusy::mStrange);
   mb_guess = leAtMt.displayMass(softsusy::mBottom);
   me_guess = model->get_thresholds() > 0 ?
      leAtMt.displayMass(softsusy::mElectron) :
      leAtMt.displayPoleMel();
   mm_guess = model->get_thresholds() > 0 ?
      leAtMt.displayMass(softsusy::mMuon) :
      leAtMt.displayPoleMmuon();
   mtau_guess = leAtMt.displayMass(softsusy::mTau);

   calculate_running_SM_masses();

   // guess gauge couplings at mt
   const auto alpha_sm(leAtMt.guess_alpha_SM5(mtpole));

   MODEL->set_g1(Sqrt(4. * Pi * alpha_sm(0)));
   MODEL->set_g2(Sqrt(4. * Pi * alpha_sm(1)));
   MODEL->set_g3(Sqrt(4. * Pi * alpha_sm(2)));


   model->set_scale(mtpole);

   // apply user-defined initial guess at the low scale
   const auto TanBeta = INPUTPARAMETER(TanBeta);
   const auto LambdaInput = INPUTPARAMETER(LambdaInput);

   MODEL->set_vd(Re(LowEnergyConstant(vev)/Sqrt(1 + Sqr(TanBeta))));
   MODEL->set_vu(Re((TanBeta*LowEnergyConstant(vev))/Sqrt(1 + Sqr(TanBeta))));
   MODEL->set_Lambdax(Re(LambdaInput));
   MODEL->set_Kappa(Re(0.1));
   MODEL->set_vS(Re(1000));
   MODEL->set_m0Sq(Re(Sqr(LowEnergyConstant(MZ))));
   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();

}

void CNMSSM_initial_guesser<Semi_analytic>::calculate_DRbar_yukawa_couplings()
{
   calculate_running_SM_masses();
   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();
}

void CNMSSM_initial_guesser<Semi_analytic>::calculate_running_SM_masses()
{
   upQuarksDRbar.setZero();
   upQuarksDRbar(0,0) = mu_guess;
   upQuarksDRbar(1,1) = mc_guess;
   upQuarksDRbar(2,2) = mt_guess;

   downQuarksDRbar.setZero();
   downQuarksDRbar(0,0) = md_guess;
   downQuarksDRbar(1,1) = ms_guess;
   downQuarksDRbar(2,2) = mb_guess;

   downLeptonsDRbar.setZero();
   downLeptonsDRbar(0,0) = me_guess;
   downLeptonsDRbar(1,1) = mm_guess;
   downLeptonsDRbar(2,2) = mtau_guess;
}

/**
 * Calculates the Yukawa couplings Yu of the up-type quarks
 * from the Standard Model up-type quark masses (ignoring threshold
 * corrections).
 */
void CNMSSM_initial_guesser<Semi_analytic>::calculate_Yu_DRbar()
{
   const auto vu = MODELPARAMETER(vu);
   MODEL->set_Yu((((1.4142135623730951*upQuarksDRbar)/vu).transpose()).real());

}

/**
 * Calculates the Yukawa couplings Yd of the down-type
 * quarks from the Standard Model down-type quark masses (ignoring
 * threshold corrections).
 */
void CNMSSM_initial_guesser<Semi_analytic>::calculate_Yd_DRbar()
{
   const auto vd = MODELPARAMETER(vd);
   MODEL->set_Yd((((1.4142135623730951*downQuarksDRbar)/vd).transpose()).real());

}

/**
 * Calculates the Yukawa couplings Ye of the leptons
 * from the Standard Model down-type lepton masses (ignoring threshold
 * corrections).
 */
void CNMSSM_initial_guesser<Semi_analytic>::calculate_Ye_DRbar()
{
   const auto vd = MODELPARAMETER(vd);
   MODEL->set_Ye((((1.4142135623730951*downLeptonsDRbar)/vd).transpose()).real());

}

void CNMSSM_initial_guesser<Semi_analytic>::initial_guess_high_scale_parameters()
{
   const double high_scale_guess = high_constraint.get_initial_scale_guess();

   model->run_to(high_scale_guess, running_precision);

   // apply high-scale constraint
   high_constraint.apply();

   // apply user-defined initial guess at the high scale
   

}

/**
 * Performs an initial iteration for the SUSY parameters, ignoring
 * threshold corrections.
 */
void CNMSSM_initial_guesser<Semi_analytic>::solve_susy_parameters()
{
   const double low_scale_guess = low_constraint.get_initial_scale_guess();
   model->run_to(low_scale_guess, running_precision);

   const auto old_thresholds = model->get_thresholds();
   const auto old_threshold_corrections = model->get_threshold_corrections();
   const auto old_loop_order = model->get_ewsb_loop_order();
   const softsusy::QedQcd old_low_sm_parameters
      = low_constraint.get_sm_parameters();
   const softsusy::QedQcd old_susy_sm_parameters
      = susy_constraint.get_sm_parameters();

   const auto save_old_settings = make_raii_guard(
      [this, old_thresholds, old_threshold_corrections, old_loop_order,
       &old_low_sm_parameters, &old_susy_sm_parameters] () {
         this->model->set_thresholds(old_thresholds);
         this->model->set_threshold_corrections(old_threshold_corrections);
         this->model->set_ewsb_loop_order(old_loop_order);
         this->low_constraint.set_sm_parameters(old_low_sm_parameters);
         this->susy_constraint.set_sm_parameters(old_susy_sm_parameters);
      });

   Threshold_corrections temp_thresholds;
   temp_thresholds.set(0);
   model->set_thresholds(0);
   model->set_threshold_corrections(temp_thresholds);
   model->set_ewsb_loop_order(0);

   softsusy::QedQcd tmp_sm_parameters(old_low_sm_parameters);

   tmp_sm_parameters.setMass(softsusy::mUp, mu_guess);
   tmp_sm_parameters.setMass(softsusy::mCharm, mc_guess);
   tmp_sm_parameters.setPoleMt(mt_guess);
   tmp_sm_parameters.setMass(softsusy::mDown, md_guess);
   tmp_sm_parameters.setMass(softsusy::mStrange, ms_guess);
   tmp_sm_parameters.setMass(softsusy::mBottom, mb_guess);
   tmp_sm_parameters.setPoleMel(me_guess);
   tmp_sm_parameters.setPoleMmuon(mm_guess);
   tmp_sm_parameters.setPoleMtau(mtau_guess);

   low_constraint.set_sm_parameters(tmp_sm_parameters);
   susy_constraint.set_sm_parameters(tmp_sm_parameters);

   low_constraint.set_is_initial_guess(true);

   CNMSSM_susy_convergence_tester<Semi_analytic> convergence_tester(
      model, running_precision);
   Two_scale_increasing_precision precision(10.0, running_precision);

   RGFlow<Two_scale> initial_solver;
   initial_solver.set_convergence_tester(&convergence_tester);
   initial_solver.set_running_precision(&precision);
   initial_solver.add(&low_constraint, model);
   initial_solver.add(&high_constraint, model);
   initial_solver.solve();

   low_constraint.set_is_initial_guess(false);
}

void CNMSSM_initial_guesser<Semi_analytic>::guess_soft_parameters()
{
   const double low_scale_guess = low_constraint.get_scale();
   const double input_scale = high_constraint.get_scale();

   model->run_to(low_scale_guess, running_precision);

   // apply EWSB constraint
   model->calculate_semi_analytic_solutions(input_scale);
   model->solve_ewsb_tree_level();

   // calculate tree-level spectrum
   model->calculate_DRbar_masses();
}

} // namespace flexiblesusy
