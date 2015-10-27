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

// File generated at Tue 27 Oct 2015 15:07:32

#ifndef SplitMSSM_SPECTRUM_GENERATOR_INTERFACE_H
#define SplitMSSM_SPECTRUM_GENERATOR_INTERFACE_H

#include "SplitMSSM_two_scale_model.hpp"
#include "SplitMSSM_utilities.hpp"

#include "spectrum_generator_settings.hpp"
#include "coupling_monitor.hpp"
#include "two_loop_corrections.hpp"
#include "error.hpp"
#include "logger.hpp"

namespace softsusy {
   class QedQcd;
}

namespace flexiblesusy {

struct SplitMSSM_input_parameters;

template <class T>
class SplitMSSM_spectrum_generator_interface {
public:
   SplitMSSM_spectrum_generator_interface()
      : model()
      , parameter_output_scale(0.)
      , precision_goal(1.0e-4)
      , reached_precision(std::numeric_limits<double>::infinity())
      , beta_zero_threshold(1.0e-11)
      , max_iterations(0)
      , beta_loop_order(2)
      , threshold_corrections_loop_order(2)
      , calculate_sm_masses(false)
      , force_output(false)
   {}
   virtual ~SplitMSSM_spectrum_generator_interface() {}

   const SplitMSSM<T>& get_model() const { return model; }
   const Problems<SplitMSSM_info::NUMBER_OF_PARTICLES>& get_problems() const {
      return model.get_problems();
   }
   int get_exit_code() const { return get_problems().have_problem(); }
   double get_reached_precision() const { return reached_precision; }
   void set_parameter_output_scale(double s) { parameter_output_scale = s; }
   void set_precision_goal(double precision_goal_) { precision_goal = precision_goal_; }
   void set_pole_mass_loop_order(unsigned l) { model.set_pole_mass_loop_order(l); }
   void set_ewsb_loop_order(unsigned l) { model.set_ewsb_loop_order(l); }
   void set_beta_loop_order(unsigned l) { beta_loop_order = l; }
   void set_beta_zero_threshold(double t) { beta_zero_threshold = t; }
   void set_max_iterations(unsigned n) { max_iterations = n; }
   void set_calculate_sm_masses(bool flag) { calculate_sm_masses = flag; }
   void set_force_output(bool flag) { force_output = flag; }
   void set_threshold_corrections_loop_order(unsigned t) { threshold_corrections_loop_order = t; }
   void set_two_loop_corrections(const Two_loop_corrections& c) { model.set_two_loop_corrections(c); }
   void set_settings(const Spectrum_generator_settings&);

   virtual void run(const softsusy::QedQcd&, const SplitMSSM_input_parameters&) = 0;
   void write_running_couplings(const std::string& filename, double, double) const;
   void write_spectrum(const std::string& filename = "SplitMSSM_spectrum.dat") const;

protected:
   SplitMSSM<T> model;
   double parameter_output_scale; ///< output scale for running parameters
   double precision_goal; ///< precision goal
   double reached_precision; ///< the precision that was reached
   double beta_zero_threshold; ///< beta function zero threshold
   unsigned max_iterations; ///< maximum number of iterations
   unsigned beta_loop_order; ///< beta-function loop order
   unsigned threshold_corrections_loop_order; ///< threshold corrections loop order
   bool calculate_sm_masses; ///< calculate SM pole masses
   bool force_output; ///< force output
};

/**
 * Setup spectrum generator from a Spectrum_generator_settings object.
 *
 * @param settings spectrum generator settings
 */
template <class T>
void SplitMSSM_spectrum_generator_interface<T>::set_settings(
   const Spectrum_generator_settings& settings)
{
   set_precision_goal(
      settings.get(Spectrum_generator_settings::precision));
   set_max_iterations(
      settings.get(Spectrum_generator_settings::max_iterations));
   set_calculate_sm_masses(
      settings.get(Spectrum_generator_settings::calculate_sm_masses) >= 1.0);
   set_force_output(
      settings.get(Spectrum_generator_settings::force_output) >= 1.0);
   set_pole_mass_loop_order(
      settings.get(Spectrum_generator_settings::pole_mass_loop_order));
   set_ewsb_loop_order(
      settings.get(Spectrum_generator_settings::ewsb_loop_order));
   set_beta_loop_order(
      settings.get(Spectrum_generator_settings::beta_loop_order));
   set_beta_zero_threshold(
      settings.get(Spectrum_generator_settings::beta_zero_threshold));
   set_threshold_corrections_loop_order(
      settings.get(Spectrum_generator_settings::threshold_corrections_loop_order));
   set_two_loop_corrections(
      settings.get_two_loop_corrections());
}

/**
 * Create a text file which contains the values of all model
 * parameters at all scales between the low-scale and the high-scale.
 *
 * @param filename name of output file
 * @param start lowest scale
 * @param stop highest scale
 */
template <class T>
void SplitMSSM_spectrum_generator_interface<T>::write_running_couplings(
   const std::string& filename,
   double start, double stop) const
{
   SplitMSSM_mass_eigenstates tmp_model(model);
   try {
      tmp_model.run_to(start);
   } catch (const Error& error) {
      ERROR("write_running_couplings: running to scale "
            << start << " failed: " << error.what());
      return;
   }

   SplitMSSM_parameter_getter parameter_getter;
   Coupling_monitor<SplitMSSM_mass_eigenstates, SplitMSSM_parameter_getter>
      coupling_monitor(tmp_model, parameter_getter);

   coupling_monitor.run(start, stop, 100, true);
   coupling_monitor.write_to_file(filename);
}

/**
 * Write spectrum (pole masses) to a text file
 *
 * @param filename output file name
 */
template <class T>
void SplitMSSM_spectrum_generator_interface<T>::write_spectrum(const std::string& filename) const
{
   SplitMSSM_spectrum_plotter plotter;
   plotter.extract_spectrum(model);
   plotter.write_to_file(filename);
}

} // namespace flexiblesusy

#endif
