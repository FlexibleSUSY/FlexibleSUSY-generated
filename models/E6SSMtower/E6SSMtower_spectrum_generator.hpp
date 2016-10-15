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

// File generated at Sat 15 Oct 2016 16:33:12

#ifndef E6SSMtower_STANDARD_MODEL_SPECTRUM_GENERATOR_H
#define E6SSMtower_STANDARD_MODEL_SPECTRUM_GENERATOR_H

#include "E6SSMtower_spectrum_generator_interface.hpp"
#include "E6SSMtower_two_scale_susy_scale_constraint.hpp"
#include "E6SSMtower_two_scale_convergence_tester.hpp"
#include "E6SSMtower_two_scale_initial_guesser.hpp"
#include "E6SSMtower_standard_model_matching.hpp"
#include "E6SSMtower_input_parameters.hpp"
#include "E6SSMtower_standard_model_two_scale_matching.hpp"
#include "standard_model_two_scale_model.hpp"
#include "standard_model_two_scale_low_scale_constraint.hpp"

#include "lowe.h"
#include "error.hpp"
#include "numerics2.hpp"
#include "two_scale_running_precision.hpp"
#include "two_scale_solver.hpp"

#include <limits>

namespace flexiblesusy {

template <class T>
class E6SSMtower_spectrum_generator
   : public E6SSMtower_spectrum_generator_interface<T> {
public:
   E6SSMtower_spectrum_generator()
      : E6SSMtower_spectrum_generator_interface<T>()
      , solver()
      , susy_scale_constraint()
      , low_scale_constraint()
      , matching()
      , eft()
      , susy_scale(0.)
      , low_scale(0.)
   {}
   virtual ~E6SSMtower_spectrum_generator() {}

   const standard_model::StandardModel<T>& get_eft() const { return eft; }
   double get_high_scale() const { return 0.; }
   double get_susy_scale() const { return susy_scale; }
   double get_low_scale()  const { return low_scale;  }

   virtual void run(const softsusy::QedQcd&, const E6SSMtower_input_parameters&);
   void write_running_couplings(const std::string& filename = "E6SSMtower_rgflow.dat") const;

private:
   RGFlow<T> solver;
   E6SSMtower_susy_scale_constraint<T> susy_scale_constraint;
   standard_model::Standard_model_low_scale_constraint<T>  low_scale_constraint;
   E6SSMtower_standard_model_Matching<T> matching;
   standard_model::StandardModel<T> eft;
   double susy_scale, low_scale;
};

/**
 * @brief Run's the RG solver with the given input parameters
 *
 * This function sets up the RG solver using a susy-scale
 * and low-scale constraint.  Afterwards the solver is run until
 * convergence is reached or an error occours.  Finally the particle
 * spectrum (pole masses) is calculated.
 *
 * @param oneset Standard Model input parameters
 * @param input model input parameters
 */
template <class T>
void E6SSMtower_spectrum_generator<T>::run(const softsusy::QedQcd& oneset,
                                const E6SSMtower_input_parameters& input)
{
   E6SSMtower<T>& model = this->model;
   model.clear();
   model.set_input_parameters(input);
   model.do_calculate_sm_pole_masses(this->settings.get(Spectrum_generator_settings::calculate_sm_masses));
   model.do_force_output(this->settings.get(Spectrum_generator_settings::force_output));
   model.set_loops(this->settings.get(Spectrum_generator_settings::beta_loop_order));
   model.set_thresholds(this->settings.get(Spectrum_generator_settings::threshold_corrections_loop_order));
   model.set_zero_threshold(this->settings.get(Spectrum_generator_settings::beta_zero_threshold));

   eft.clear();
   eft.do_force_output(this->settings.get(Spectrum_generator_settings::force_output));
   eft.set_loops(this->settings.get(Spectrum_generator_settings::beta_loop_order));
   eft.set_thresholds(this->settings.get(Spectrum_generator_settings::threshold_corrections_loop_order));
   eft.set_zero_threshold(this->settings.get(Spectrum_generator_settings::beta_zero_threshold));
   eft.set_pole_mass_loop_order(this->model.get_pole_mass_loop_order());
   eft.set_ewsb_loop_order(this->model.get_ewsb_loop_order());
   eft.set_ewsb_iteration_precision(this->model.get_ewsb_iteration_precision());
   eft.set_number_of_ewsb_iterations(this->model.get_number_of_ewsb_iterations());
   eft.set_number_of_mass_iterations(this->model.get_number_of_mass_iterations());
   eft.set_two_loop_corrections(this->model.get_two_loop_corrections());

   susy_scale_constraint.clear();
   low_scale_constraint .clear();

   // needed for constraint::initialize()
   susy_scale_constraint.set_model(&model);
   low_scale_constraint .set_model(&eft);

   low_scale_constraint .set_sm_parameters(oneset);

   const unsigned index = this->settings.get(Spectrum_generator_settings::eft_higgs_index);

   matching.set_models(&eft, &model);
   matching.set_constraint(&susy_scale_constraint);
   matching.set_scale(this->settings.get(Spectrum_generator_settings::eft_matching_scale));
   matching.set_loop_order_up(this->settings.get(Spectrum_generator_settings::eft_matching_loop_order_up));
   matching.set_loop_order_down(this->settings.get(Spectrum_generator_settings::eft_matching_loop_order_down));
   matching.set_higgs_index(index);

   susy_scale_constraint.initialize();
   low_scale_constraint .initialize();


   std::vector<Constraint<T>*> model_constraints(1);
   model_constraints[0] = &susy_scale_constraint;

   std::vector<Constraint<T>*> eft_constraints(1);
   eft_constraints[0] = &low_scale_constraint;

   E6SSMtower_convergence_tester<T> convergence_tester(
      &model, this->settings.get(Spectrum_generator_settings::precision));
   if (this->settings.get(Spectrum_generator_settings::max_iterations) > 0)
      convergence_tester.set_max_iterations(
         this->settings.get(Spectrum_generator_settings::max_iterations));

   E6SSMtower_standard_model_initial_guesser<T> initial_guesser(&model, &eft, oneset,
                                                  low_scale_constraint,
                                                  susy_scale_constraint);
   Two_scale_increasing_precision precision(
      10.0, this->settings.get(Spectrum_generator_settings::precision));

   solver.reset();
   solver.set_convergence_tester(&convergence_tester);
   solver.set_running_precision(&precision);
   solver.set_initial_guesser(&initial_guesser);
   solver.add_model(&eft, &matching, eft_constraints);
   solver.add_model(&model, model_constraints);

   susy_scale = low_scale = 0.;
   this->reached_precision = std::numeric_limits<double>::infinity();

   try {
      solver.solve();
      susy_scale = susy_scale_constraint.get_scale();
      low_scale  = low_scale_constraint.get_scale();
      this->reached_precision = convergence_tester.get_current_accuracy();

      const double mass_scale =
         this->settings.get(Spectrum_generator_settings::pole_mass_scale) != 0. ?
         this->settings.get(Spectrum_generator_settings::pole_mass_scale) : susy_scale;

      model.run_to(mass_scale);
      model.solve_ewsb();
      model.calculate_spectrum();

      // start logarithmic resummation
      {
         double Q_higgs = this->settings.get(Spectrum_generator_settings::eft_pole_mass_scale);

         if (Q_higgs == 0.) {
            const double Mt = low_scale_constraint.get_sm_parameters().displayPoleMt();
            Q_higgs = std::min(susy_scale, Mt);
         }

         eft.run_to(Q_higgs);

         // computation of pole mass spectrum in the SM
         const unsigned eft_pole_loops = eft.get_pole_mass_loop_order();
         const unsigned eft_ewsb_loops = eft.get_ewsb_loop_order();

         if (eft_pole_loops > 1) {
            WARNING("Pole mass loop order " << eft_pole_loops
                    << " in the EFT is currently not supported!  I'm using 1-loop.");
            eft.set_pole_mass_loop_order(1);
         }
         if (eft_ewsb_loops > 1) {
            WARNING("EWSB loop order " << eft_ewsb_loops
                    << " in the EFT is currently not supported!  I'm using 1-loop.");
            eft.set_ewsb_loop_order(1);
         }

         eft.calculate_DRbar_masses();
         eft.solve_ewsb();
         eft.calculate_spectrum();

         eft.set_pole_mass_loop_order(eft_pole_loops);
         eft.set_ewsb_loop_order(eft_ewsb_loops);

         this->model.get_physical().Mhh(index) = eft.get_physical().Mhh;
         this->model.get_physical().MVZ = eft.get_physical().MVZ;
         this->model.get_physical().MVWm = eft.get_physical().MVWp;
         this->model.get_physical().MFu(0) = eft.get_physical().MFu(0);
         this->model.get_physical().MFu(1) = eft.get_physical().MFu(1);
         this->model.get_physical().MFu(2) = eft.get_physical().MFu(2);
         this->model.get_physical().MFd(0) = eft.get_physical().MFd(0);
         this->model.get_physical().MFd(1) = eft.get_physical().MFd(1);
         this->model.get_physical().MFd(2) = eft.get_physical().MFd(2);
         this->model.get_physical().MFe(0) = eft.get_physical().MFe(0);
         this->model.get_physical().MFe(1) = eft.get_physical().MFe(1);
         this->model.get_physical().MFe(2) = eft.get_physical().MFe(2);

      }

      // copy calculated W pole mass
      model.get_physical().MVWm
         = low_scale_constraint.get_sm_parameters().displayPoleMW();

      // run to output scale (if scale > 0)
      if (!is_zero(this->parameter_output_scale)) {
         model.run_to(this->parameter_output_scale);
         eft.run_to(this->parameter_output_scale);
      }
   } catch (...) {
      this->translate_exception_to_problem(model);
   }
}

/**
 * Create a text file which contains the values of all model
 * parameters at all scales between the low-scale and the high-scale.
 *
 * @param filename name of output file
 */
template <class T>
void E6SSMtower_spectrum_generator<T>::write_running_couplings(
   const std::string& filename) const
{
   E6SSMtower_spectrum_generator_interface<T>::write_running_couplings(filename, low_scale, susy_scale);
}

} // namespace flexiblesusy

#endif
