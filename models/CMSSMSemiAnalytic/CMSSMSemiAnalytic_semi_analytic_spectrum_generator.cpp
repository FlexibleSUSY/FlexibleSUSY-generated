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

// File generated at Wed 16 Oct 2019 22:56:08

#include "CMSSMSemiAnalytic_semi_analytic_spectrum_generator.hpp"
#include "CMSSMSemiAnalytic_input_parameters.hpp"
#include "CMSSMSemiAnalytic_semi_analytic_convergence_tester.hpp"
#include "CMSSMSemiAnalytic_semi_analytic_ewsb_solver.hpp"
#include "CMSSMSemiAnalytic_semi_analytic_high_scale_constraint.hpp"
#include "CMSSMSemiAnalytic_semi_analytic_initial_guesser.hpp"
#include "CMSSMSemiAnalytic_semi_analytic_low_scale_constraint.hpp"
#include "CMSSMSemiAnalytic_semi_analytic_soft_parameters_constraint.hpp"
#include "CMSSMSemiAnalytic_semi_analytic_susy_convergence_tester.hpp"
#include "CMSSMSemiAnalytic_semi_analytic_susy_scale_constraint.hpp"

#include "error.hpp"
#include "lowe.h"
#include "numerics2.hpp"
#include "two_scale_running_precision.hpp"
#include "semi_analytic_solver.hpp"

#include <limits>

namespace flexiblesusy {

double CMSSMSemiAnalytic_spectrum_generator<Semi_analytic>::get_pole_mass_scale() const
{
   return settings.get(Spectrum_generator_settings::pole_mass_scale) != 0. ?
      settings.get(Spectrum_generator_settings::pole_mass_scale) :
      get_susy_scale();
}

/**
 * @brief Run's the RG solver with the given input parameters
 *
 * This function sets up the RG solver using a high-scale, susy-scale,
 * low-scale and soft parameters constraint.  Afterwards the solver
 * is run until convergence is reached or an error occours.  Finally
 * the particle spectrum (pole masses) is calculated.
 *
 * @param qedqcd Standard Model input parameters
 * @param input model input parameters
 */
void CMSSMSemiAnalytic_spectrum_generator<Semi_analytic>::run_except(const softsusy::QedQcd& qedqcd,
                                const CMSSMSemiAnalytic_input_parameters& input)
{
   VERBOSE_MSG("Solving BVP using semi-analytic solver");

   problems.set_bvp_solver_problems({ BVP_solver_problems("SemiAnalyticSolver") });

   auto& model = this->model;
   model.clear();
   model.set_input_parameters(input);
   model.do_calculate_sm_pole_masses(
      settings.get(Spectrum_generator_settings::calculate_sm_masses));
   model.do_calculate_bsm_pole_masses(
      settings.get(Spectrum_generator_settings::calculate_bsm_masses));
   model.do_force_output(
      settings.get(Spectrum_generator_settings::force_output));
   model.set_loops(
      settings.get(Spectrum_generator_settings::beta_loop_order));
   model.set_thresholds(
      settings.get(
         Spectrum_generator_settings::threshold_corrections_loop_order));
   model.set_zero_threshold(
      settings.get(Spectrum_generator_settings::beta_zero_threshold));

   CMSSMSemiAnalytic_semi_analytic_solutions& solutions(
      model.get_semi_analytic_solutions());

   CMSSMSemiAnalytic_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_semi_analytic_solutions(&solutions);
   model.set_ewsb_solver(
      std::make_shared<CMSSMSemiAnalytic_ewsb_solver<Semi_analytic> >(ewsb_solver));

   CMSSMSemiAnalytic_high_scale_constraint<Semi_analytic> high_scale_constraint(&model);
   CMSSMSemiAnalytic_susy_scale_constraint<Semi_analytic> susy_scale_constraint(&model, qedqcd);
   CMSSMSemiAnalytic_low_scale_constraint<Semi_analytic>  low_scale_constraint(&model, qedqcd);
   CMSSMSemiAnalytic_soft_parameters_constraint<Semi_analytic> soft_constraint(&model, qedqcd);

   soft_constraint.set_boundary_scale(
      [this,&high_scale_constraint] () {
         return high_scale_constraint.get_scale(); });
   susy_scale_constraint.set_scale(
      [this,&soft_constraint] () {
         return soft_constraint.get_scale(); });

   high_scale_constraint.initialize();
   susy_scale_constraint.initialize();
   low_scale_constraint .initialize();
   soft_constraint      .initialize();

   CMSSMSemiAnalytic_susy_convergence_tester<Semi_analytic> inner_ct(
      &model, settings.get(Spectrum_generator_settings::precision));
   CMSSMSemiAnalytic_convergence_tester<Semi_analytic> outer_ct(
      &model, settings.get(Spectrum_generator_settings::precision));

   if (settings.get(Spectrum_generator_settings::pole_mass_scale) != 0.) {
      outer_ct.set_scale_getter(
         [this](){
            return settings.get(Spectrum_generator_settings::pole_mass_scale);
         });
   }

   if (settings.get(Spectrum_generator_settings::max_iterations) > 0) {
      inner_ct.set_max_iterations(
         settings.get(Spectrum_generator_settings::max_iterations));
      outer_ct.set_max_iterations(
         settings.get(Spectrum_generator_settings::max_iterations));
   }

   CMSSMSemiAnalytic_initial_guesser<Semi_analytic> initial_guesser(&model, qedqcd,
                                                              low_scale_constraint,
                                                              susy_scale_constraint,
                                                              high_scale_constraint);

   Two_scale_increasing_precision precision(
      10.0, settings.get(Spectrum_generator_settings::precision));

   RGFlow<Semi_analytic> solver;
   solver.set_inner_convergence_tester(&inner_ct);
   solver.set_outer_convergence_tester(&outer_ct);
   solver.set_running_precision(&precision);
   solver.set_initial_guesser(&initial_guesser);
   solver.add_inner(&low_scale_constraint, &model);
   solver.add_inner(&high_scale_constraint, &model);
   solver.add_inner(&susy_scale_constraint, &model);
   solver.add_outer(&soft_constraint, &model);

   high_scale = susy_scale = low_scale = 0.;
   reached_precision = std::numeric_limits<double>::infinity();

   solver.solve();

   // impose low-scale constraint one last time
   model.run_to(low_scale_constraint.get_scale());
   low_scale_constraint.apply();

   high_scale = high_scale_constraint.get_scale();
   susy_scale = susy_scale_constraint.get_scale();
   low_scale  = low_scale_constraint.get_scale();
   reached_precision = outer_ct.get_current_accuracy();

   calculate_spectrum();

   // copy calculated W pole mass
   model.get_physical().MVWm
      = low_scale_constraint.get_sm_parameters().displayPoleMW();

   // run to output scale (if scale > 0)
   if (!is_zero(parameter_output_scale))
      model.run_to(parameter_output_scale);
}

/**
 * Create a text file which contains the values of all model
 * parameters at all scales between the low-scale and the high-scale.
 *
 * @param filename name of output file
 */
void CMSSMSemiAnalytic_spectrum_generator<Semi_analytic>::write_running_couplings(
   const std::string& filename) const
{
   CMSSMSemiAnalytic_spectrum_generator_interface<Semi_analytic>::write_running_couplings(
      filename, get_low_scale(), get_high_scale());
}

void CMSSMSemiAnalytic_spectrum_generator<Semi_analytic>::calculate_spectrum()
{
   model.run_to(get_pole_mass_scale());
   model.calculate_semi_analytic_solutions(get_high_scale());
   model.solve_ewsb();
   model.calculate_spectrum();
}

} // namespace flexiblesusy
