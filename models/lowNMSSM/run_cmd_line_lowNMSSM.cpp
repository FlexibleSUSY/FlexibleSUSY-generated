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

// File generated at Wed 12 Apr 2017 12:41:02

#include "lowNMSSM_input_parameters.hpp"
#include "lowNMSSM_spectrum_generator.hpp"
#include "lowNMSSM_observables.hpp"
#include "lowNMSSM_slha_io.hpp"

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
      "Usage: run_cmd_line_lowNMSSM.x [options]\n"
      "Options:\n"
      "  --Qin=<value>\n"
      "  --M1Input=<value>\n"
      "  --M2Input=<value>\n"
      "  --M3Input=<value>\n"
      "  --AtInput=<value>\n"
      "  --AbInput=<value>\n"
      "  --ATauInput=<value>\n"
      "  --TanBeta=<value>\n"
      "  --ml1Input=<value>\n"
      "  --ml2Input=<value>\n"
      "  --ml3Input=<value>\n"
      "  --me1Input=<value>\n"
      "  --me2Input=<value>\n"
      "  --me3Input=<value>\n"
      "  --mq1Input=<value>\n"
      "  --mq2Input=<value>\n"
      "  --mq3Input=<value>\n"
      "  --md1Input=<value>\n"
      "  --md2Input=<value>\n"
      "  --md3Input=<value>\n"
      "  --mu1Input=<value>\n"
      "  --mu2Input=<value>\n"
      "  --mu3Input=<value>\n"
      "  --LambdaInput=<value>\n"
      "  --KappaInput=<value>\n"
      "  --ALambdaInput=<value>\n"
      "  --AKappaInput=<value>\n"
      "  --MuEffInput=<value>\n"

      "  --help,-h                         print this help message"
             << std::endl;
}

void set_command_line_parameters(int argc, char* argv[],
                                 lowNMSSM_input_parameters& input)
{
   for (int i = 1; i < argc; ++i) {
      const char* option = argv[i];

      if(Command_line_options::get_parameter_value(option, "--Qin=", input.Qin))
         continue;

      if(Command_line_options::get_parameter_value(option, "--M1Input=", input.M1Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--M2Input=", input.M2Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--M3Input=", input.M3Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--AtInput=", input.AtInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--AbInput=", input.AbInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--ATauInput=", input.ATauInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--TanBeta=", input.TanBeta))
         continue;

      if(Command_line_options::get_parameter_value(option, "--ml1Input=", input.ml1Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--ml2Input=", input.ml2Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--ml3Input=", input.ml3Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--me1Input=", input.me1Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--me2Input=", input.me2Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--me3Input=", input.me3Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--mq1Input=", input.mq1Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--mq2Input=", input.mq2Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--mq3Input=", input.mq3Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--md1Input=", input.md1Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--md2Input=", input.md2Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--md3Input=", input.md3Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--mu1Input=", input.mu1Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--mu2Input=", input.mu2Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--mu3Input=", input.mu3Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--LambdaInput=", input.LambdaInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--KappaInput=", input.KappaInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--ALambdaInput=", input.ALambdaInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--AKappaInput=", input.AKappaInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--MuEffInput=", input.MuEffInput))
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

   lowNMSSM_input_parameters input;
   set_command_line_parameters(argc, argv, input);

   Physical_input physical_input;
   softsusy::QedQcd qedqcd;

   try {
      qedqcd.to(qedqcd.displayPoleMZ()); // run SM fermion masses to MZ
   } catch (const Error& e) {
      ERROR(e.what());
      return EXIT_FAILURE;
   }

   lowNMSSM_spectrum_generator<algorithm_type> spectrum_generator;
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
   const lowNMSSM_slha<algorithm_type> model(spectrum_generator.get_model());

   lowNMSSM_scales scales;
   scales.HighScale = spectrum_generator.get_high_scale();
   scales.SUSYScale = spectrum_generator.get_susy_scale();
   scales.LowScale  = spectrum_generator.get_low_scale();

   const lowNMSSM_observables observables(calculate_observables(model, qedqcd, physical_input));

   // SLHA output
   SLHAea::Coll slhaea(lowNMSSM_slha_io::fill_slhaea(model, qedqcd, scales, observables));

   std::cout << slhaea;

   return exit_code;
}
