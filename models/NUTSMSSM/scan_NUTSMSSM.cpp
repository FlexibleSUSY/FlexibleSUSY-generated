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

// File generated at Fri 8 Jan 2016 15:17:31

#include "NUTSMSSM_input_parameters.hpp"
#include "NUTSMSSM_spectrum_generator.hpp"
#include "NUTSMSSM_two_scale_model_slha.hpp"

#include "command_line_options.hpp"
#include "scan.hpp"
#include "lowe.h"
#include "logger.hpp"

#include <iostream>
#include <cstring>

namespace flexiblesusy {

void print_usage()
{
   std::cout <<
      "Usage: scan_NUTSMSSM.x [options]\n"
      "Options:\n"
      "  --m0=<value>\n"
      "  --m12=<value>\n"
      "  --TanBeta=<value>\n"
      "  --Azero=<value>\n"
      "  --LambdaInput=<value>\n"
      "  --KappaInput=<value>\n"
      "  --LambdaSInput=<value>\n"
      "  --L1Input=<value>\n"
      "  --MSInput=<value>\n"
      "  --BInput=<value>\n"
      "  --MuInput=<value>\n"
      "  --LInput=<value>\n"

      "  --help,-h                         print this help message"
             << std::endl;
}

void set_command_line_parameters(int argc, char* argv[],
                                 NUTSMSSM_input_parameters& input)
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

      if(Command_line_options::get_parameter_value(option, "--KappaInput=", input.KappaInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--LambdaSInput=", input.LambdaSInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--L1Input=", input.L1Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--MSInput=", input.MSInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--BInput=", input.BInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--MuInput=", input.MuInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--LInput=", input.LInput))
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

   NUTSMSSM_input_parameters input;
   set_command_line_parameters(argc, argv, input);

   softsusy::QedQcd qedqcd;
   qedqcd.toMz();

   NUTSMSSM_spectrum_generator<algorithm_type> spectrum_generator;
   spectrum_generator.set_precision_goal(1.0e-4);
   spectrum_generator.set_max_iterations(0);         // 0 == automatic
   spectrum_generator.set_calculate_sm_masses(0);    // 0 == no
   spectrum_generator.set_parameter_output_scale(0); // 0 == susy scale

   const std::vector<double> range(float_range(0., 100., 10));

   cout << "# "
        << std::setw(12) << std::left << "m0" << ' '
        << std::setw(12) << std::left << "Mhh(0)/GeV" << ' '
        << std::setw(12) << std::left << "error"
        << '\n';

   for (std::vector<double>::const_iterator it = range.begin(),
           end = range.end(); it != end; ++it) {
      input.m0 = *it;

      spectrum_generator.run(qedqcd, input);

      const NUTSMSSM_slha<algorithm_type> model(spectrum_generator.get_model());
      const NUTSMSSM_physical& pole_masses = model.get_physical_slha();
      const Problems<NUTSMSSM_info::NUMBER_OF_PARTICLES>& problems
         = spectrum_generator.get_problems();
      const double higgs = pole_masses.Mhh(0);
      const bool error = problems.have_problem();

      cout << "  "
           << std::setw(12) << std::left << input.m0 << ' '
           << std::setw(12) << std::left << higgs << ' '
           << std::setw(12) << std::left << error;
      if (error) {
         cout << "\t# " << problems;
      }
      cout << '\n';
   }

   return 0;
}