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

// File generated at Fri 8 Jan 2016 12:13:56

#include "MRSSM_input_parameters.hpp"
#include "MRSSM_observables.hpp"
#include "MRSSM_slha_io.hpp"
#include "MRSSM_spectrum_generator.hpp"
#include "MRSSM_utilities.hpp"

#include "spectrum_generator_settings.hpp"
#include "lowe.h"
#include "command_line_options.hpp"

#include <iostream>
#include <cstdlib>

int main(int argc, const char* argv[])
{
   using namespace flexiblesusy;
   typedef Two_scale algorithm_type;

   Command_line_options options(argc, argv);
   if (options.must_print_model_info())
      MRSSM_info::print(std::cout);
   if (options.must_exit())
      return options.status();

   const std::string database_output_file(options.get_database_output_file());
   const std::string rgflow_file(options.get_rgflow_file());
   const std::string slha_input_source(options.get_slha_input_file());
   const std::string slha_output_file(options.get_slha_output_file());
   const std::string spectrum_file(options.get_spectrum_file());
   MRSSM_slha_io slha_io;
   Spectrum_generator_settings spectrum_generator_settings;
   softsusy::QedQcd qedqcd;
   MRSSM_input_parameters input;

   if (slha_input_source.empty()) {
      ERROR("No SLHA input source given!\n"
            "   Please provide one via the option --slha-input-file=");
      return EXIT_FAILURE;
   }

   try {
      slha_io.read_from_source(slha_input_source);
      slha_io.fill(qedqcd);
      slha_io.fill(input);
      slha_io.fill(spectrum_generator_settings);
   } catch (const Error& error) {
      ERROR(error.what());
      return EXIT_FAILURE;
   }

   qedqcd.toMz(); // run SM fermion masses to MZ

   MRSSM_spectrum_generator<algorithm_type> spectrum_generator;
   spectrum_generator.set_settings(spectrum_generator_settings);
   spectrum_generator.set_parameter_output_scale(
      slha_io.get_parameter_output_scale());

   spectrum_generator.run(qedqcd, input);

   const MRSSM_slha<algorithm_type> model(spectrum_generator.get_model());
   const Problems<MRSSM_info::NUMBER_OF_PARTICLES>& problems
      = spectrum_generator.get_problems();

   MRSSM_scales scales;
   scales.HighScale = spectrum_generator.get_high_scale();
   scales.SUSYScale = spectrum_generator.get_susy_scale();
   scales.LowScale  = spectrum_generator.get_low_scale();

   Observables observables;
   if (spectrum_generator_settings.get(Spectrum_generator_settings::calculate_observables))
      observables = calculate_observables(model, qedqcd);

   // SLHA output
   if (!slha_output_file.empty()) {
      slha_io.set_spinfo(problems);
      slha_io.set_minpar(input);
      slha_io.set_extpar(input);
      if (!problems.have_problem() ||
          spectrum_generator_settings.get(Spectrum_generator_settings::force_output)) {
         slha_io.set_spectrum(model);
         slha_io.set_extra(model, scales, observables);
      }

      if (slha_output_file == "-") {
         slha_io.write_to_stream(std::cout);
      } else {
         slha_io.write_to_file(slha_output_file);
      }
   }

   if (!database_output_file.empty() &&
       (!problems.have_problem() ||
        spectrum_generator_settings.get(Spectrum_generator_settings::force_output))) {
      MRSSM_database::to_database(database_output_file, model, &qedqcd, &observables);
   }

   if (!spectrum_file.empty())
      spectrum_generator.write_spectrum(spectrum_file);

   if (!rgflow_file.empty())
      spectrum_generator.write_running_couplings(rgflow_file);

   const int exit_code = spectrum_generator.get_exit_code();

   return exit_code;
}
