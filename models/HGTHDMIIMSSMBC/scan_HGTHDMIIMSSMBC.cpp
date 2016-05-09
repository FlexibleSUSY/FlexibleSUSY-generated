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

// File generated at Mon 9 May 2016 12:07:08

#include "HGTHDMIIMSSMBC_input_parameters.hpp"
#include "HGTHDMIIMSSMBC_spectrum_generator.hpp"
#include "HGTHDMIIMSSMBC_two_scale_model_slha.hpp"

#include "command_line_options.hpp"
#include "scan.hpp"
#include "lowe.h"
#include "logger.hpp"

#include <iostream>
#include <cstring>

#define INPUTPARAMETER(p) input.p

namespace flexiblesusy {

void print_usage()
{
   std::cout <<
      "Usage: scan_HGTHDMIIMSSMBC.x [options]\n"
      "Options:\n"
      "  --TanBeta=<value>\n"
      "  --MSUSY=<value>\n"
      "  --MEWSB=<value>\n"
      "  --MuInput=<value>\n"
      "  --M1Input=<value>\n"
      "  --M2Input=<value>\n"
      "  --M3Input=<value>\n"
      "  --MAInput=<value>\n"
      "  --AtInput=<value>\n"
      "  --AbInput=<value>\n"
      "  --AtauInput=<value>\n"
      "  --LambdaLoopOrder=<value>\n"

      "  --help,-h                         print this help message"
             << std::endl;
}

void set_command_line_parameters(int argc, char* argv[],
                                 HGTHDMIIMSSMBC_input_parameters& input)
{
   for (int i = 1; i < argc; ++i) {
      const char* option = argv[i];

      if(Command_line_options::get_parameter_value(option, "--TanBeta=", input.TanBeta))
         continue;

      if(Command_line_options::get_parameter_value(option, "--MSUSY=", input.MSUSY))
         continue;

      if(Command_line_options::get_parameter_value(option, "--MEWSB=", input.MEWSB))
         continue;

      if(Command_line_options::get_parameter_value(option, "--MuInput=", input.MuInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--M1Input=", input.M1Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--M2Input=", input.M2Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--M3Input=", input.M3Input))
         continue;

      if(Command_line_options::get_parameter_value(option, "--MAInput=", input.MAInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--AtInput=", input.AtInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--AbInput=", input.AbInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--AtauInput=", input.AtauInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--LambdaLoopOrder=", input.LambdaLoopOrder))
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

   HGTHDMIIMSSMBC_input_parameters input;
   set_command_line_parameters(argc, argv, input);

   softsusy::QedQcd qedqcd;

   try {
      qedqcd.to(qedqcd.displayPoleMZ()); // run SM fermion masses to MZ
   } catch (const std::string& s) {
      ERROR(s);
      return EXIT_FAILURE;
   }

   HGTHDMIIMSSMBC_spectrum_generator<algorithm_type> spectrum_generator;
   spectrum_generator.set_precision_goal(1.0e-4);
   spectrum_generator.set_max_iterations(0);         // 0 == automatic
   spectrum_generator.set_calculate_sm_masses(0);    // 0 == no
   spectrum_generator.set_parameter_output_scale(0); // 0 == susy scale

   const std::vector<double> range(float_range(0., 100., 10));

   cout << "# "
        << std::setw(12) << std::left << "TanBeta" << ' '
        << std::setw(12) << std::left << "Mhh(0)/GeV" << ' '
        << std::setw(12) << std::left << "error"
        << '\n';

   for (std::vector<double>::const_iterator it = range.begin(),
           end = range.end(); it != end; ++it) {
      INPUTPARAMETER(TanBeta) = *it;


      spectrum_generator.run(qedqcd, input);

      const HGTHDMIIMSSMBC_slha<algorithm_type> model(spectrum_generator.get_model());
      const HGTHDMIIMSSMBC_physical& pole_masses = model.get_physical_slha();
      const Problems<HGTHDMIIMSSMBC_info::NUMBER_OF_PARTICLES>& problems
         = spectrum_generator.get_problems();
      const double higgs = pole_masses.Mhh(0);
      const bool error = problems.have_problem();

      cout << "  "
           << std::setw(12) << std::left << *it << ' '
           << std::setw(12) << std::left << higgs << ' '
           << std::setw(12) << std::left << error;
      if (error) {
         cout << "\t# " << problems;
      }
      cout << '\n';
   }

   return 0;
}
