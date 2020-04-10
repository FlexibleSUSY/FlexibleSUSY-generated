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

// File generated at Fri 10 Apr 2020 20:21:32

#include "config.h"

#include "lowNMSSMTanBetaAtMZ_input_parameters.hpp"
#include "lowNMSSMTanBetaAtMZ_observables.hpp"
#include "lowNMSSMTanBetaAtMZ_slha_io.hpp"
#include "lowNMSSMTanBetaAtMZ_spectrum_generator.hpp"

#ifdef ENABLE_TWO_SCALE_SOLVER
#include "lowNMSSMTanBetaAtMZ_two_scale_spectrum_generator.hpp"
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
      "Usage: run_cmd_line_lowNMSSMTanBetaAtMZ.x [options]\n"
      "Options:\n"
      "  --TanBeta=<value>\n"
      "  --Qin=<value>\n"
      "  --M1Input=<value>\n"
      "  --M2Input=<value>\n"
      "  --M3Input=<value>\n"
      "  --AtInput=<value>\n"
      "  --AbInput=<value>\n"
      "  --ATauInput=<value>\n"
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

      "  --solver-type=<value>             an integer corresponding\n"
      "                                    to the solver type to use\n"
      "  --help,-h                         print this help message"
             << std::endl;
}

void set_command_line_parameters(const Dynamic_array_view<char*>& args,
                                 lowNMSSMTanBetaAtMZ_input_parameters& input,
                                 int& solver_type)
{
   for (int i = 1; i < args.size(); ++i) {
      const auto option = args[i];

      if(Command_line_options::get_parameter_value(option, "--TanBeta=", input.TanBeta))
         continue;

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
int run_solver(const lowNMSSMTanBetaAtMZ_input_parameters& input)
{
   Physical_input physical_input;
   softsusy::QedQcd qedqcd;

   Spectrum_generator_settings settings;
   settings.set(Spectrum_generator_settings::precision, 1.0e-4);

   lowNMSSMTanBetaAtMZ_spectrum_generator<solver_type> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);

   auto model = std::get<0>(spectrum_generator.get_models_slha());

   lowNMSSMTanBetaAtMZ_scales scales;
   scales.HighScale = spectrum_generator.get_high_scale();
   scales.SUSYScale = spectrum_generator.get_susy_scale();
   scales.LowScale  = spectrum_generator.get_low_scale();
   scales.pole_mass_scale = spectrum_generator.get_pole_mass_scale();

   const auto observables = calculate_observables(
      model, qedqcd, physical_input, scales.pole_mass_scale);

   // SLHA output
   SLHAea::Coll slhaea(lowNMSSMTanBetaAtMZ_slha_io::fill_slhaea(
                          model, qedqcd, scales, observables));

   std::cout << slhaea;

   return spectrum_generator.get_exit_code();
}

int run(int solver_type, const lowNMSSMTanBetaAtMZ_input_parameters& input)
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

   lowNMSSMTanBetaAtMZ_input_parameters input;
   int solver_type = 0;
   set_command_line_parameters(make_dynamic_array_view(&argv[0], argc), input,
                               solver_type);

   const int exit_code = run(solver_type, input);

   return exit_code;
}
