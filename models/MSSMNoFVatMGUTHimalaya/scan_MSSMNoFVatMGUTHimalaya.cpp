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


#include "config.h"

#include "MSSMNoFVatMGUTHimalaya_input_parameters.hpp"
#include "MSSMNoFVatMGUTHimalaya_model_slha.hpp"
#include "MSSMNoFVatMGUTHimalaya_spectrum_generator.hpp"

#ifdef ENABLE_TWO_SCALE_SOLVER
#include "MSSMNoFVatMGUTHimalaya_two_scale_spectrum_generator.hpp"
#endif

#include "command_line_options.hpp"
#include "array_view.hpp"
#include "scan.hpp"
#include "lowe.h"
#include "logger.hpp"

#include <iostream>
#include <iomanip>
#include <string>

#define INPUTPARAMETER(p) input.p

namespace flexiblesusy {

void print_usage()
{
   std::cout <<
      "Usage: scan_MSSMNoFVatMGUTHimalaya.x [options]\n"
      "Options:\n"
      "  --TanBeta=<value>\n"
      "  --SignMu=<value>\n"
      "  --M1=<value>\n"
      "  --M2=<value>\n"
      "  --M3=<value>\n"
      "  --AtIN=<value>\n"
      "  --AbIN=<value>\n"
      "  --AtauIN=<value>\n"
      "  --AcIN=<value>\n"
      "  --AsIN=<value>\n"
      "  --AmuonIN=<value>\n"
      "  --AuIN=<value>\n"
      "  --AdIN=<value>\n"
      "  --AeIN=<value>\n"
      "  --mHd2IN=<value>\n"
      "  --mHu2IN=<value>\n"
      "  --ml11IN=<value>\n"
      "  --ml22IN=<value>\n"
      "  --ml33IN=<value>\n"
      "  --me11IN=<value>\n"
      "  --me22IN=<value>\n"
      "  --me33IN=<value>\n"
      "  --mq11IN=<value>\n"
      "  --mq22IN=<value>\n"
      "  --mq33IN=<value>\n"
      "  --mu11IN=<value>\n"
      "  --mu22IN=<value>\n"
      "  --mu33IN=<value>\n"
      "  --md11IN=<value>\n"
      "  --md22IN=<value>\n"
      "  --md33IN=<value>\n"
      "  --Mlow=<value>\n"

      "  --solver-type=<value>             an integer corresponding\n"
      "                                    to the solver type to use\n"
      "  --loop-library=<value>            an integer corresponding to the used\n"
      "                                    realization of loop library\n"
      "  --help,-h                         print this help message"
             << std::endl;
}

void set_command_line_parameters(const Dynamic_array_view<char*>& args,
                                 MSSMNoFVatMGUTHimalaya_input_parameters& input,
                                 int& solver_type,
                                 int& loop_library)
{
   for (int i = 1; i < args.size(); ++i) {
      const std::string option = args[i];

      if(Command_line_options::get_parameter_value(option, "--TanBeta=", input.TanBeta))
         continue;

      if(Command_line_options::get_parameter_value(option, "--SignMu=", input.SignMu))
         continue;

      if(Command_line_options::get_parameter_value(option, "--M1=", input.M1))
         continue;

      if(Command_line_options::get_parameter_value(option, "--M2=", input.M2))
         continue;

      if(Command_line_options::get_parameter_value(option, "--M3=", input.M3))
         continue;

      if(Command_line_options::get_parameter_value(option, "--AtIN=", input.AtIN))
         continue;

      if(Command_line_options::get_parameter_value(option, "--AbIN=", input.AbIN))
         continue;

      if(Command_line_options::get_parameter_value(option, "--AtauIN=", input.AtauIN))
         continue;

      if(Command_line_options::get_parameter_value(option, "--AcIN=", input.AcIN))
         continue;

      if(Command_line_options::get_parameter_value(option, "--AsIN=", input.AsIN))
         continue;

      if(Command_line_options::get_parameter_value(option, "--AmuonIN=", input.AmuonIN))
         continue;

      if(Command_line_options::get_parameter_value(option, "--AuIN=", input.AuIN))
         continue;

      if(Command_line_options::get_parameter_value(option, "--AdIN=", input.AdIN))
         continue;

      if(Command_line_options::get_parameter_value(option, "--AeIN=", input.AeIN))
         continue;

      if(Command_line_options::get_parameter_value(option, "--mHd2IN=", input.mHd2IN))
         continue;

      if(Command_line_options::get_parameter_value(option, "--mHu2IN=", input.mHu2IN))
         continue;

      if(Command_line_options::get_parameter_value(option, "--ml11IN=", input.ml11IN))
         continue;

      if(Command_line_options::get_parameter_value(option, "--ml22IN=", input.ml22IN))
         continue;

      if(Command_line_options::get_parameter_value(option, "--ml33IN=", input.ml33IN))
         continue;

      if(Command_line_options::get_parameter_value(option, "--me11IN=", input.me11IN))
         continue;

      if(Command_line_options::get_parameter_value(option, "--me22IN=", input.me22IN))
         continue;

      if(Command_line_options::get_parameter_value(option, "--me33IN=", input.me33IN))
         continue;

      if(Command_line_options::get_parameter_value(option, "--mq11IN=", input.mq11IN))
         continue;

      if(Command_line_options::get_parameter_value(option, "--mq22IN=", input.mq22IN))
         continue;

      if(Command_line_options::get_parameter_value(option, "--mq33IN=", input.mq33IN))
         continue;

      if(Command_line_options::get_parameter_value(option, "--mu11IN=", input.mu11IN))
         continue;

      if(Command_line_options::get_parameter_value(option, "--mu22IN=", input.mu22IN))
         continue;

      if(Command_line_options::get_parameter_value(option, "--mu33IN=", input.mu33IN))
         continue;

      if(Command_line_options::get_parameter_value(option, "--md11IN=", input.md11IN))
         continue;

      if(Command_line_options::get_parameter_value(option, "--md22IN=", input.md22IN))
         continue;

      if(Command_line_options::get_parameter_value(option, "--md33IN=", input.md33IN))
         continue;

      if(Command_line_options::get_parameter_value(option, "--Mlow=", input.Mlow))
         continue;

      
      if (Command_line_options::get_parameter_value(
             option, "--solver-type=", solver_type))
         continue;

      if (Command_line_options::get_parameter_value(
             option, "--loop-library=", loop_library))
         continue;

      if (option == "--help" || option == "-h") {
         print_usage();
         exit(EXIT_SUCCESS);
      }

      ERROR("Unrecognized command line option: " << option);
      exit(EXIT_FAILURE);
   }
}

struct MSSMNoFVatMGUTHimalaya_scan_result {
   Spectrum_generator_problems problems;
   double higgs{0.};
};

template <class solver_type>
MSSMNoFVatMGUTHimalaya_scan_result run_parameter_point(int loop_library, const softsusy::QedQcd& qedqcd,
   MSSMNoFVatMGUTHimalaya_input_parameters& input)
{
   Spectrum_generator_settings settings;
   settings.set(Spectrum_generator_settings::precision, 1.0e-4);
   settings.set(Spectrum_generator_settings::loop_library, loop_library);

   MSSMNoFVatMGUTHimalaya_spectrum_generator<solver_type> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);

   const auto model = std::get<0>(spectrum_generator.get_models_slha());
   const auto& pole_masses = model.get_physical_slha();

   MSSMNoFVatMGUTHimalaya_scan_result result;
   result.problems = spectrum_generator.get_problems();
   result.higgs = pole_masses.Mhh(0);

   return result;
}

void scan(int solver_type, int loop_library, MSSMNoFVatMGUTHimalaya_input_parameters& input,
          const std::vector<double>& range)
{
   softsusy::QedQcd qedqcd;

   for (const auto p: range) {
      INPUTPARAMETER(TanBeta) = p;

      MSSMNoFVatMGUTHimalaya_scan_result result;
      switch (solver_type) {
      case 0:
#ifdef ENABLE_TWO_SCALE_SOLVER
      case 1:
         result = run_parameter_point<Two_scale>(loop_library, qedqcd, input);
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

   MSSMNoFVatMGUTHimalaya_input_parameters input;
   int solver_type = 1;
   int loop_library = 0;
   set_command_line_parameters(make_dynamic_array_view(&argv[0], argc), input,
                               solver_type, loop_library);

   std::cout << "# "
             << std::setw(12) << std::left << "TanBeta" << ' '
             << std::setw(12) << std::left << "Mhh(0)/GeV" << ' '
             << std::setw(12) << std::left << "error"
             << '\n';

   const std::vector<double> range(float_range(0., 100., 10));

   scan(solver_type, loop_library, input, range);

   return 0;
}
