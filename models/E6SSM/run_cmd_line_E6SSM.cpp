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

// File generated at Sun 26 Aug 2018 14:29:15

#include "config.h"

#include "E6SSM_input_parameters.hpp"
#include "E6SSM_observables.hpp"
#include "E6SSM_slha_io.hpp"
#include "E6SSM_spectrum_generator.hpp"

#ifdef ENABLE_TWO_SCALE_SOLVER
#include "E6SSM_two_scale_spectrum_generator.hpp"
#endif

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
      "Usage: run_cmd_line_E6SSM.x [options]\n"
      "Options:\n"
      "  --m0=<value>\n"
      "  --m12=<value>\n"
      "  --TanBeta=<value>\n"
      "  --Azero=<value>\n"
      "  --LambdaInput=<value>\n"
      "  --KappaInput=<value>\n"
      "  --muPrimeInput=<value>\n"
      "  --BmuPrimeInput=<value>\n"
      "  --vSInput=<value>\n"
      "  --Lambda12Input=<value>\n"

      "  --solver-type=<value>             an integer corresponding\n"
      "                                    to the solver type to use\n"
      "  --help,-h                         print this help message"
             << std::endl;
}

void set_command_line_parameters(const Dynamic_array_view<char*>& args,
                                 E6SSM_input_parameters& input,
                                 int& solver_type)
{
   for (int i = 1; i < args.size(); ++i) {
      const auto option = args[i];

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

      if(Command_line_options::get_parameter_value(option, "--muPrimeInput=", input.muPrimeInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--BmuPrimeInput=", input.BmuPrimeInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--vSInput=", input.vSInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--Lambda12Input=", input.Lambda12Input))
         continue;

      
      if (Command_line_options::get_parameter_value(
             option, "--solver-type=", solver_type))
         continue;

      if (strcmp(option,"--help") == 0 || strcmp(option,"-h") == 0) {
         print_usage();
         exit(EXIT_SUCCESS);
      }

      ERROR("Unrecognized command line option: " << option);
      exit(EXIT_FAILURE);
   }
}

template<class solver_type>
int run_solver(const E6SSM_input_parameters& input)
{
   Physical_input physical_input;
   softsusy::QedQcd qedqcd;

   Spectrum_generator_settings settings;
   settings.set(Spectrum_generator_settings::precision, 1.0e-4);

   E6SSM_spectrum_generator<solver_type> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);

   auto model = std::get<0>(spectrum_generator.get_models_slha());

   E6SSM_scales scales;
   scales.HighScale = spectrum_generator.get_high_scale();
   scales.SUSYScale = spectrum_generator.get_susy_scale();
   scales.LowScale  = spectrum_generator.get_low_scale();
   scales.pole_mass_scale = spectrum_generator.get_pole_mass_scale();

   const auto observables = calculate_observables(
      model, qedqcd, physical_input, scales.pole_mass_scale);

   // SLHA output
   SLHAea::Coll slhaea(E6SSM_slha_io::fill_slhaea(
                          model, qedqcd, scales, observables));

   std::cout << slhaea;

   return spectrum_generator.get_exit_code();
}

int run(int solver_type, const E6SSM_input_parameters& input)
{
   int exit_code = 0;

   switch (solver_type) {
   case 0:
#ifdef ENABLE_TWO_SCALE_SOLVER
   case 1:
      exit_code = run_solver<Two_scale>(input);
      if (!exit_code || solver_type != 0) break;
#endif

   default:
      if (solver_type != 0) {
         ERROR("unknown solver type: " << solver_type);
         exit_code = -1;
      }
      break;
   }

   return exit_code;
}

} // namespace flexiblesusy


int main(int argc, char* argv[])
{
   using namespace flexiblesusy;

   E6SSM_input_parameters input;
   int solver_type = 0;
   set_command_line_parameters(make_dynamic_array_view(&argv[0], argc), input,
                               solver_type);

   const int exit_code = run(solver_type, input);

   return exit_code;
}
