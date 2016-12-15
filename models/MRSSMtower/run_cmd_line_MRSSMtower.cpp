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

// File generated at Thu 15 Dec 2016 12:44:45

#include "MRSSMtower_input_parameters.hpp"
#include "MRSSMtower_spectrum_generator.hpp"
#include "MRSSMtower_observables.hpp"
#include "MRSSMtower_slha_io.hpp"

#include "command_line_options.hpp"
#include "lowe.h"
#include "logger.hpp"
#include "physical_input.hpp"

#include <iostream>
#include <cstring>

namespace flexiblesusy {

void print_usage()
{
   std::cout <<
      "Usage: run_cmd_line_MRSSMtower.x [options]\n"
      "Options:\n"
      "  --TanBeta=<value>\n"
      "  --MS=<value>\n"
      "  --LamTDInput=<value>\n"
      "  --LamTUInput=<value>\n"
      "  --LamSDInput=<value>\n"
      "  --LamSUInput=<value>\n"
      "  --MuDInput=<value>\n"
      "  --MuUInput=<value>\n"
      "  --BMuInput=<value>\n"
      "  --mS2Input=<value>\n"
      "  --mT2Input=<value>\n"
      "  --moc2Input=<value>\n"
      "  --mRd2Input=<value>\n"
      "  --mRu2Input=<value>\n"
      "  --MDBSInput=<value>\n"
      "  --MDWBTInput=<value>\n"
      "  --MDGocInput=<value>\n"

      "  --help,-h                         print this help message"
             << std::endl;
}

void set_command_line_parameters(int argc, char* argv[],
                                 MRSSMtower_input_parameters& input)
{
   for (int i = 1; i < argc; ++i) {
      const char* option = argv[i];

      if(Command_line_options::get_parameter_value(option, "--TanBeta=", input.TanBeta))
         continue;

      if(Command_line_options::get_parameter_value(option, "--MS=", input.MS))
         continue;

      if(Command_line_options::get_parameter_value(option, "--LamTDInput=", input.LamTDInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--LamTUInput=", input.LamTUInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--LamSDInput=", input.LamSDInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--LamSUInput=", input.LamSUInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--MuDInput=", input.MuDInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--MuUInput=", input.MuUInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--BMuInput=", input.BMuInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--mS2Input=", input.mS2Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--mT2Input=", input.mT2Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--moc2Input=", input.moc2Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--mRd2Input=", input.mRd2Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--mRu2Input=", input.mRu2Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--MDBSInput=", input.MDBSInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--MDWBTInput=", input.MDWBTInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--MDGocInput=", input.MDGocInput))
         continue;

      
      if (strcmp(option,"--help") == 0 || strcmp(option,"-h") == 0) {
         print_usage();
         exit(EXIT_SUCCESS);
      }

      ERROR("Unrecognized command line option: " << option);
      exit(EXIT_FAILURE);
   }
}

} // namespace flexiblesusy


int main(int argc, char* argv[])
{
   using namespace flexiblesusy;
   typedef Two_scale algorithm_type;

   MRSSMtower_input_parameters input;
   set_command_line_parameters(argc, argv, input);

   Physical_input physical_input;
   softsusy::QedQcd qedqcd;

   try {
      qedqcd.to(qedqcd.displayPoleMZ()); // run SM fermion masses to MZ
   } catch (const Error& e) {
      ERROR(e.what());
      return EXIT_FAILURE;
   }

   MRSSMtower_spectrum_generator<algorithm_type> spectrum_generator;
   spectrum_generator.set_precision_goal(1.0e-4);
   spectrum_generator.set_beta_zero_threshold(1e-11);
   spectrum_generator.set_max_iterations(0);         // 0 == automatic
   spectrum_generator.set_calculate_sm_masses(0);    // 0 == no
   spectrum_generator.set_parameter_output_scale(0); // 0 == susy scale
   spectrum_generator.set_pole_mass_loop_order(2);   // 2-loop
   spectrum_generator.set_ewsb_loop_order(2);        // 2-loop
   spectrum_generator.set_beta_loop_order(2);        // 2-loop
   spectrum_generator.set_threshold_corrections_loop_order(1); // 1-loop

   spectrum_generator.run(qedqcd, input);

   const int exit_code = spectrum_generator.get_exit_code();
   const MRSSMtower_slha<algorithm_type> model(spectrum_generator.get_model());

   MRSSMtower_scales scales;
   scales.HighScale = spectrum_generator.get_high_scale();
   scales.SUSYScale = spectrum_generator.get_susy_scale();
   scales.LowScale  = spectrum_generator.get_low_scale();

   const MRSSMtower_observables observables(calculate_observables(model, qedqcd, physical_input));

   // SLHA output
   SLHAea::Coll slhaea(MRSSMtower_slha_io::fill_slhaea(model, qedqcd, scales, observables));

   std::cout << slhaea;

   return exit_code;
}
