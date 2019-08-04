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

// File generated at Sun 4 Aug 2019 19:30:52

#include "config.h"

#include "MRSSM_input_parameters.hpp"
#include "MRSSM_model_slha.hpp"
#include "MRSSM_spectrum_generator.hpp"

#ifdef ENABLE_TWO_SCALE_SOLVER
#include "MRSSM_two_scale_spectrum_generator.hpp"
#endif

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
      "Usage: scan_MRSSM.x [options]\n"
      "Options:\n"
      "  --TanBeta=<value>\n"
      "  --LamTDInput=<value>\n"
      "  --LamTUInput=<value>\n"
      "  --LamSDInput=<value>\n"
      "  --LamSUInput=<value>\n"
      "  --MuInput=<value>\n"
      "  --MuDInput=<value>\n"
      "  --MuUInput=<value>\n"
      "  --vTInput=<value>\n"
      "  --vSInput=<value>\n"
      "  --BMuInput=<value>\n"
      "  --BMuDInput=<value>\n"
      "  --BMuUInput=<value>\n"
      "  --moc2Input=<value>\n"
      "  --mRd2Input=<value>\n"
      "  --mRu2Input=<value>\n"
      "  --MDBSInput=<value>\n"
      "  --MDWBTInput=<value>\n"
      "  --MDGocInput=<value>\n"

      "  --solver-type=<value>             an integer corresponding\n"
      "                                    to the solver type to use\n"
      "  --help,-h                         print this help message"
             << std::endl;
}

void set_command_line_parameters(const Dynamic_array_view<char*>& args,
                                 MRSSM_input_parameters& input,
                                 int& solver_type)
{
   for (int i = 1; i < args.size(); ++i) {
      const auto option = args[i];

      if(Command_line_options::get_parameter_value(option, "--TanBeta=", input.TanBeta))
         continue;

      if(Command_line_options::get_parameter_value(option, "--LamTDInput=", input.LamTDInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--LamTUInput=", input.LamTUInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--LamSDInput=", input.LamSDInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--LamSUInput=", input.LamSUInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--MuInput=", input.MuInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--MuDInput=", input.MuDInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--MuUInput=", input.MuUInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--vTInput=", input.vTInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--vSInput=", input.vSInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--BMuInput=", input.BMuInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--BMuDInput=", input.BMuDInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--BMuUInput=", input.BMuUInput))
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

struct MRSSM_scan_result {
   Spectrum_generator_problems problems;
   double higgs{0.};
};

template <class solver_type>
MRSSM_scan_result run_parameter_point(const softsusy::QedQcd& qedqcd,
   MRSSM_input_parameters& input)
{
   Spectrum_generator_settings settings;
   settings.set(Spectrum_generator_settings::precision, 1.0e-4);

   MRSSM_spectrum_generator<solver_type> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);

   const auto model = std::get<0>(spectrum_generator.get_models_slha());
   const auto& pole_masses = model.get_physical_slha();

   MRSSM_scan_result result;
   result.problems = spectrum_generator.get_problems();
   result.higgs = pole_masses.Mhh(0);

   return result;
}

void scan(int solver_type, MRSSM_input_parameters& input,
          const std::vector<double>& range)
{
   softsusy::QedQcd qedqcd;

   for (const auto p: range) {
      INPUTPARAMETER(TanBeta) = p;

      MRSSM_scan_result result;
      switch (solver_type) {
      case 0:
#ifdef ENABLE_TWO_SCALE_SOLVER
      case 1:
         result = run_parameter_point<Two_scale>(qedqcd, input);
         if (!result.problems.have_problem() || solver_type != 0) break;
#endif

      default:
         if (solver_type != 0) {
            ERROR("unknown solver type: " << solver_type);
            exit(EXIT_FAILURE);
         }
      }

      const int error = result.problems.have_problem();
      std::cout << "  "
                << std::setw(12) << std::left << p << ' '
                << std::setw(12) << std::left << result.higgs << ' '
                << std::setw(12) << std::left << error;
      if (error) {
         std::cout << "\t# " << result.problems;
      }
      std::cout << '\n';
   }
}

} // namespace flexiblesusy


int main(int argc, char* argv[])
{
   using namespace flexiblesusy;

   MRSSM_input_parameters input;
   int solver_type = 1;
   set_command_line_parameters(make_dynamic_array_view(&argv[0], argc), input,
                               solver_type);

   std::cout << "# "
             << std::setw(12) << std::left << "TanBeta" << ' '
             << std::setw(12) << std::left << "Mhh(0)/GeV" << ' '
             << std::setw(12) << std::left << "error"
             << '\n';

   const std::vector<double> range(float_range(0., 100., 10));

   scan(solver_type, input, range);

   return 0;
}
