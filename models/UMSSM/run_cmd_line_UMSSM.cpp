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

// File generated at Fri 8 Jan 2016 15:18:35

#include "UMSSM_input_parameters.hpp"
#include "UMSSM_observables.hpp"
#include "UMSSM_spectrum_generator.hpp"
#include "UMSSM_slha_io.hpp"

#include "command_line_options.hpp"
#include "lowe.h"
#include "logger.hpp"

#include <iostream>
#include <cstring>

namespace flexiblesusy {

void print_usage()
{
   std::cout <<
      "Usage: run_cmd_line_UMSSM.x [options]\n"
      "Options:\n"
      "  --m0=<value>\n"
      "  --m12=<value>\n"
      "  --TanBeta=<value>\n"
      "  --Azero=<value>\n"
      "  --LambdaInput=<value>\n"
      "  --ALambdaInput=<value>\n"
      "  --vSInput=<value>\n"
      "  --Qq=<value>\n"
      "  --Ql=<value>\n"
      "  --QHd=<value>\n"
      "  --QHu=<value>\n"
      "  --Qd=<value>\n"
      "  --Qu=<value>\n"
      "  --Qe=<value>\n"
      "  --Qs=<value>\n"
      "  --Qv=<value>\n"

      "  --help,-h                         print this help message"
             << std::endl;
}

void set_command_line_parameters(int argc, char* argv[],
                                 UMSSM_input_parameters& input)
{
   for (int i = 1; i < argc; ++i) {
      const char* option = argv[i];

      if(Command_line_options::get_parameter_value(option, "--m0=", input.m0))
         continue;

      if(Command_line_options::get_parameter_value(option, "--m12=", input.m12))
         continue;

      if(Command_line_options::get_parameter_value(option, "--TanBeta=", input.TanBeta))
         continue;

      if(Command_line_options::get_parameter_value(option, "--Azero=", input.Azero))
         continue;

      if(Command_line_options::get_parameter_value(option, "--LambdaInput=", input.LambdaInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--ALambdaInput=", input.ALambdaInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--vSInput=", input.vSInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--Qq=", input.Qq))
         continue;

      if(Command_line_options::get_parameter_value(option, "--Ql=", input.Ql))
         continue;

      if(Command_line_options::get_parameter_value(option, "--QHd=", input.QHd))
         continue;

      if(Command_line_options::get_parameter_value(option, "--QHu=", input.QHu))
         continue;

      if(Command_line_options::get_parameter_value(option, "--Qd=", input.Qd))
         continue;

      if(Command_line_options::get_parameter_value(option, "--Qu=", input.Qu))
         continue;

      if(Command_line_options::get_parameter_value(option, "--Qe=", input.Qe))
         continue;

      if(Command_line_options::get_parameter_value(option, "--Qs=", input.Qs))
         continue;

      if(Command_line_options::get_parameter_value(option, "--Qv=", input.Qv))
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

   UMSSM_input_parameters input;
   set_command_line_parameters(argc, argv, input);

   softsusy::QedQcd qedqcd;
   qedqcd.toMz();

   UMSSM_spectrum_generator<algorithm_type> spectrum_generator;
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
   const UMSSM_slha<algorithm_type> model(spectrum_generator.get_model());

   UMSSM_scales scales;
   scales.HighScale = spectrum_generator.get_high_scale();
   scales.SUSYScale = spectrum_generator.get_susy_scale();
   scales.LowScale  = spectrum_generator.get_low_scale();

   const Observables observables(calculate_observables(model, qedqcd));

   // SLHA output
   SLHAea::Coll slhaea(UMSSM_slha_io::fill_slhaea(model, qedqcd, scales, observables));

   std::cout << slhaea;

   return exit_code;
}
