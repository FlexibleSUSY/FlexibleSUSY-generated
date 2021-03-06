#include <fstream>
#include <string>
#include <iostream>

#include "config.h"

#include "decays/CMSSMSemiAnalytic_decays.hpp"
#include "CMSSMSemiAnalytic_slha_io.hpp"
#ifdef ENABLE_SEMI_ANALYTIC_SOLVER
#include "CMSSMSemiAnalytic_semi_analytic_spectrum_generator.hpp"
#endif

#include "physical_input.hpp"

int main(int argc, char* argv[])
{
   if (argc < 2) {
      ERROR("Need at least one argument - the name of the file containing a spectrum.");
      return EXIT_FAILURE;
   }
   else if (argc > 3) {
      ERROR("The program expects at most 2 arguments.");
      return EXIT_FAILURE;
   }

   std::ifstream ifs(argv[1]);
   std::string slha_input {
      std::istreambuf_iterator<char>(ifs),
      std::istreambuf_iterator<char>()
   };

   std::stringstream istr(slha_input);

   using namespace flexiblesusy;
   CMSSMSemiAnalytic_slha_io slha_io;
   slha_io.read_from_stream(istr);

   softsusy::QedQcd qedqcd;
   CMSSMSemiAnalytic_input_parameters input;
   Spectrum_generator_settings settings;
   FlexibleDecay_settings flexibledecay_settings;

   // extract the input parameters from spectrum string
   slha_io.fill(settings);
   slha_io.fill(flexibledecay_settings);
   slha_io.fill(qedqcd);
   slha_io.fill(input);

   CMSSMSemiAnalytic_spectrum_generator<Semi_analytic> spectrum_generator;

   settings.set(Spectrum_generator_settings::calculate_sm_masses, 1.0);
   settings.set(Spectrum_generator_settings::calculate_bsm_masses, 1.0);
   settings.set(Spectrum_generator_settings::loop_library, 1.0);
   spectrum_generator.set_settings(settings);

   spectrum_generator.run(qedqcd, input);

   CMSSMSemiAnalytic_slha m = std::get<0>(spectrum_generator.get_models_slha());

   Loop_library::set(1);
   Physical_input physical_input;
   CMSSMSemiAnalytic_decays decays(m, qedqcd, physical_input, flexibledecay_settings);
   decays.calculate_decays();

   const bool show_decays = !decays.get_problems().have_problem() ||
       settings.get(Spectrum_generator_settings::force_output);

   if (show_decays) {
      slha_io.set_dcinfo(decays.get_problems());
      slha_io.set_decays(decays.get_decay_table(), flexibledecay_settings);
   }

   slha_io.write_to(argc > 2 ? argv[2] : "-");

   return 0;
}
