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

#include "THDMII_info.hpp"
#include "THDMII_input_parameters.hpp"
#include "THDMII_observables.hpp"
#include "THDMII_physical.hpp"
#include "THDMII_slha_io.hpp"
#include "THDMII_model_slha.hpp"

#ifdef ENABLE_TWO_SCALE_SOLVER
#include "THDMII_two_scale_model.hpp"
#include "THDMII_two_scale_spectrum_generator.hpp"
#endif

#include "array_view.hpp"
#include "bvp_solver_problems_format_mathlink.hpp"
#include "error.hpp"
#include "for_each.hpp"
#include "observable_problems_format_mathlink.hpp"
#include "physical_input.hpp"
#include "problems_format_mathlink.hpp"
#include "slha_io.hpp"
#include "spectrum_generator_settings.hpp"
#include "decays/flexibledecay_settings.hpp"
#include "standard_model_two_scale_model.hpp"
#include "lowe.h"
#include "THDMII_unitarity.hpp"


#include <mathlink.h>
#include "mathlink_utils.hpp"
#include <WolframLibrary.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>

#define INPUTPARAMETER(p) input.p
#define MODELPARAMETER(p) model.get_##p()
#define PHYSICALPARAMETER(p) model.get_physical().p
#define OBSERVABLE(o) observables.o

namespace flexiblesusy {
namespace THDMII_librarylink {

using Handle = long;

class Redirect_output {
public:
   explicit Redirect_output(MLINK link_)
      : link(link_)
      , buffer()
      , old_cout(std::cout.rdbuf(buffer.rdbuf()))
      , old_cerr(std::cerr.rdbuf(buffer.rdbuf()))
      {}

   ~Redirect_output() {
      std::cout.rdbuf(old_cout);
      std::cerr.rdbuf(old_cerr);
      flush();
   }

private:
   MLINK link;               ///< redirect to this link
   std::stringstream buffer; ///< buffer caching stdout
   std::streambuf* old_cout; ///< original stdout buffer
   std::streambuf* old_cerr; ///< original stderr buffer

   void flush() {
      std::string line;
      while (std::getline(buffer, line)) {
         MLPutFunction(link, "CompoundExpression", 2);
         MLPutFunction(link, "FSTHDMIIMessage", 1);
         MLPutString(link, line.c_str());
      }
   }
};

class EUnknownHandle : public Error {
public:
   explicit EUnknownHandle(Handle hid_) : Error("Unknown handle"), hid(hid_) {}
   virtual ~EUnknownHandle() = default;
   std::string what_detailed() const override {
      return std::string(what()) + ": " + ToString(hid);
   }
   Handle hid;
};

class ENotEnoughFreeHandles : public Error {
public:
   explicit ENotEnoughFreeHandles(std::size_t max_handles_)
      : Error("Maximum number of open handles reached")
      , max_handles(max_handles_) {}
   virtual ~ENotEnoughFreeHandles() = default;
   std::string what_detailed() const override {
      return std::string(what()) + ": "
         + ToString(max_handles) + ".  Please close some handles!";
   }
   std::size_t max_handles;
};

class EWrongNumberOfParameters : public Error {
public:
   EWrongNumberOfParameters(mint pars_, mint expected_)
      : Error("Wrong number of arguments")
      , pars(pars_), expected(expected_) {}
   virtual ~EWrongNumberOfParameters() = default;
   std::string what_detailed() const override {
      return std::string(what()) + ": " + ToString(pars)
         + ".  Expected: " + ToString(expected);
   }
   mint pars, expected;
};

class EInvalidSpectrum : public Error {
public:
   EInvalidSpectrum() : Error("Invalid spectrum") {}
   virtual ~EInvalidSpectrum() = default;
};

class THDMII_spectrum {
public:
   virtual ~THDMII_spectrum() = default;

   virtual void put_model_spectra(MLINK link) const = 0;

   virtual const Spectrum_generator_problems& get_problems() const = 0;
   virtual void fill_slha_io(THDMII_slha_io&, const Spectrum_generator_settings&, const FlexibleDecay_settings&) const = 0;
   virtual double get_model_scale() const = 0;
   virtual const THDMII_observables& get_observables() const = 0;

   virtual const UnitarityInfiniteS& get_unitarity() const = 0;
   virtual void calculate_spectrum(
      const Spectrum_generator_settings&, const SLHA_io::Modsel&,
      const softsusy::QedQcd&, const THDMII_input_parameters&) = 0;
   virtual void calculate_model_observables(
      const softsusy::QedQcd&,
      
      const Physical_input&,
      const Spectrum_generator_settings&) = 0;

   virtual void calculate_unitarity() = 0;
};

template <typename Solver_type>
class THDMII_spectrum_impl : public THDMII_spectrum
{
public:
   virtual ~THDMII_spectrum_impl() = default;

   virtual void put_model_spectra(MLINK link) const override;

   virtual const Spectrum_generator_problems& get_problems() const override { return problems; }
   virtual void fill_slha_io(THDMII_slha_io&, const Spectrum_generator_settings&, const FlexibleDecay_settings&) const override;
   virtual double get_model_scale() const override { return std::get<0>(models).get_scale(); }
   virtual const THDMII_observables& get_observables() const override { return observables; }

   virtual const UnitarityInfiniteS& get_unitarity() const { return unitarityData; };
   virtual void calculate_spectrum(
      const Spectrum_generator_settings&, const SLHA_io::Modsel&,
      const softsusy::QedQcd&, const THDMII_input_parameters&) override;
   virtual void calculate_model_observables(
      const softsusy::QedQcd&,
      
      const Physical_input&,
      const Spectrum_generator_settings&) override;

   virtual void calculate_unitarity() override;
private:
   std::tuple<THDMII<Solver_type>> models{};        ///< running parameters and pole masses
   Spectrum_generator_problems problems{};   ///< spectrum generator problems
   THDMII_scales scales{};              ///< scale information
   THDMII_observables observables{};    ///< observables
   UnitarityInfiniteS unitarityData = {};    ///< unitarity constraints

};

class Model_data {
public:
   Model_data() = default;
   Model_data(const Model_data&) = delete;
   Model_data(Model_data&&) = default;
   ~Model_data() = default;
   Model_data& operator=(const Model_data&) = delete;
   Model_data& operator=(Model_data&&) = default;

   void set_input_parameters(const THDMII_input_parameters& input_) { input = input_; }
   void set_physical_input(const Physical_input& p) { physical_input = p; }
   void set_sm_input_parameters(const softsusy::QedQcd& qedqcd_) { qedqcd = qedqcd_; }
   void set_settings(const Spectrum_generator_settings& s) { settings = s; }
   
   void set_fd_settings(const FlexibleDecay_settings& s) { flexibledecay_settings = s; }
   void set_modsel(const SLHA_io::Modsel& m) { modsel = m; }

   const Spectrum_generator_settings& get_settings() const { return settings; }
   const FlexibleDecay_settings& get_fd_settings() const { return flexibledecay_settings; }

   void put_settings(MLINK link) const;
   void put_sm_input_parameters(MLINK link) const;
   void put_input_parameters(MLINK link) const;
   void put_observables(MLINK link) const;
   void put_slha(MLINK link) const;

   void put_problems(MLINK link) const;
   void put_warnings(MLINK link) const;
   void put_model_spectra(MLINK link) const;
   void put_unitarity(MLINK link) const;
   void calculate_spectrum();
   void check_spectrum(MLINK link) const;
   void calculate_model_observables();

   void calculate_unitarity();
   double get_model_scale() const {
      check_spectrum_pointer();
      return spectrum->get_model_scale();
   }
private:
   THDMII_input_parameters input{};                    ///< model input parameters
   Physical_input physical_input{};                         ///< extra non-SLHA physical input
   softsusy::QedQcd qedqcd{};                               ///< SLHA physical input
   Spectrum_generator_settings settings{};                  ///< spectrum generator settings
   FlexibleDecay_settings flexibledecay_settings {};        ///< FlexibleDecay settings
   SLHA_io::Modsel modsel{};                                ///< MODSEL input
   
   std::unique_ptr<THDMII_spectrum> spectrum{nullptr}; ///< spectrum information

   THDMII_slha_io get_slha_io() const;

   void check_spectrum_pointer() const;
};

/// current handles
using Handle_map = std::map<Handle, Model_data>;
Handle_map handles;

/******************************************************************/

Handle get_new_handle()
{
   static const std::size_t max_handles =
      static_cast<std::size_t>(std::numeric_limits<Handle>::max());

   if (handles.size() >= max_handles)
      throw ENotEnoughFreeHandles(handles.size());

   Handle hid = 0;

   while (handles.find(hid) != handles.end())
      hid++;

   return hid;
}

/******************************************************************/

auto find_handle(Handle hid) -> decltype(handles.find(hid))
{
   const auto handle = handles.find(hid);

   if (handle == handles.end())
      throw EUnknownHandle(hid);

   return handle;
}

/******************************************************************/

Model_data& find_data(Handle hid)
{
   return find_handle(hid)->second;
}

/******************************************************************/

Handle get_handle_from(MLINK link)
{
   Handle hid;
   MLGet(link, &hid);

   return hid;
}

/******************************************************************/

Handle get_handle_from(MArgument arg)
{
   return MArgument_getInteger(arg);
}

/******************************************************************/

long number_of_args(MLINK link, const std::string& head)
{
   long argc;

   if (!MLCheckFunction(link, head.c_str(), &argc))
      std::cerr << "Error: argument is not a " << head << std::endl;

   return argc;
}

/******************************************************************/

bool check_number_of_args(MLINK link, long number_of_arguments,
                          const std::string& function_name)
{
   const auto n_given = number_of_args(link, "List");
   const bool ok = n_given == number_of_arguments;

   if (!ok) {
      std::cerr << "Error: " << function_name << " expects "
                << ToString(number_of_arguments) << " argument ("
                << n_given << " given)." << std::endl;
   }

   return ok;
}

/******************************************************************/

void put_error_output(MLINK link)
{
   MLPutSymbol(link, "$Failed");
}

/******************************************************************/

void put_message(MLINK link,
                 const std::string& function_name,
                 const std::string& message_tag,
                 const std::string& message_str)
{
   MLPutFunction(link, "CompoundExpression", 2);
   MLPutFunction(link, "Message", 2);
   MLPutFunction(link, "MessageName", 2);
   MLPutSymbol(link, function_name.c_str());
   MLPutString(link, message_tag.c_str());
   MLPutString(link, message_str.c_str());
}

/******************************************************************/

void Model_data::check_spectrum_pointer() const
{
   if (!spectrum) {
      throw SetupError("No spectrum generator set! "
                       "Did you run FSTHDMIICalculateSpectrum[]?");
   }
}

/******************************************************************/

void Model_data::put_settings(MLINK link) const
{
   MLPutFunction(link, "List", Spectrum_generator_settings::NUMBER_OF_OPTIONS);

   MLPutRuleTo(link, settings.get(Spectrum_generator_settings::precision), "precisionGoal");
   MLPutRuleTo(link, static_cast<int>(settings.get(Spectrum_generator_settings::max_iterations)), "maxIterations");
   MLPutRuleTo(link, static_cast<int>(settings.get(Spectrum_generator_settings::solver)), "solver");
   MLPutRuleTo(link, static_cast<int>(settings.get(Spectrum_generator_settings::calculate_sm_masses)), "calculateStandardModelMasses");
   MLPutRuleTo(link, static_cast<int>(settings.get(Spectrum_generator_settings::pole_mass_loop_order)), "poleMassLoopOrder");
   MLPutRuleTo(link, static_cast<int>(settings.get(Spectrum_generator_settings::ewsb_loop_order)), "ewsbLoopOrder");
   MLPutRuleTo(link, static_cast<int>(settings.get(Spectrum_generator_settings::beta_loop_order)), "betaFunctionLoopOrder");
   MLPutRuleTo(link, static_cast<int>(settings.get(Spectrum_generator_settings::threshold_corrections_loop_order)), "thresholdCorrectionsLoopOrder");
   MLPutRuleTo(link, static_cast<int>(settings.get(Spectrum_generator_settings::higgs_2loop_correction_at_as)), "higgs2loopCorrectionAtAs");
   MLPutRuleTo(link, static_cast<int>(settings.get(Spectrum_generator_settings::higgs_2loop_correction_ab_as)), "higgs2loopCorrectionAbAs");
   MLPutRuleTo(link, static_cast<int>(settings.get(Spectrum_generator_settings::higgs_2loop_correction_at_at)), "higgs2loopCorrectionAtAt");
   MLPutRuleTo(link, static_cast<int>(settings.get(Spectrum_generator_settings::higgs_2loop_correction_atau_atau)), "higgs2loopCorrectionAtauAtau");
   MLPutRuleTo(link, static_cast<int>(settings.get(Spectrum_generator_settings::force_output)), "forceOutput");
   MLPutRuleTo(link, static_cast<int>(settings.get(Spectrum_generator_settings::top_pole_qcd_corrections)), "topPoleQCDCorrections");
   MLPutRuleTo(link, settings.get(Spectrum_generator_settings::beta_zero_threshold), "betaZeroThreshold");
   MLPutRuleTo(link, static_cast<int>(settings.get(Spectrum_generator_settings::force_positive_masses)), "forcePositiveMasses");
   MLPutRuleTo(link, settings.get(Spectrum_generator_settings::pole_mass_scale), "poleMassScale");
   MLPutRuleTo(link, settings.get(Spectrum_generator_settings::eft_pole_mass_scale), "eftPoleMassScale");
   MLPutRuleTo(link, settings.get(Spectrum_generator_settings::eft_matching_scale), "eftMatchingScale");
   MLPutRuleTo(link, static_cast<int>(settings.get(Spectrum_generator_settings::eft_matching_loop_order_up)), "eftMatchingLoopOrderUp");
   MLPutRuleTo(link, static_cast<int>(settings.get(Spectrum_generator_settings::eft_matching_loop_order_down)), "eftMatchingLoopOrderDown");
   MLPutRuleTo(link, static_cast<int>(settings.get(Spectrum_generator_settings::eft_higgs_index)), "eftHiggsIndex");
   MLPutRuleTo(link, static_cast<int>(settings.get(Spectrum_generator_settings::calculate_bsm_masses)), "calculateBSMMasses");
   MLPutRuleTo(link, static_cast<int>(settings.get(Spectrum_generator_settings::threshold_corrections)), "thresholdCorrections");
   MLPutRuleTo(link, static_cast<int>(settings.get(Spectrum_generator_settings::higgs_3loop_ren_scheme_atb_as2)), "higgs3loopCorrectionRenScheme");
   MLPutRuleTo(link, static_cast<int>(settings.get(Spectrum_generator_settings::higgs_3loop_correction_at_as2)), "higgs3loopCorrectionAtAsAs");
   MLPutRuleTo(link, static_cast<int>(settings.get(Spectrum_generator_settings::higgs_3loop_correction_ab_as2)), "higgs3loopCorrectionAbAsAs");
   MLPutRuleTo(link, static_cast<int>(settings.get(Spectrum_generator_settings::higgs_3loop_correction_at2_as)), "higgs3loopCorrectionAtAtAs");
   MLPutRuleTo(link, static_cast<int>(settings.get(Spectrum_generator_settings::higgs_3loop_correction_at3)), "higgs3loopCorrectionAtAtAt");
   MLPutRuleTo(link, static_cast<int>(settings.get(Spectrum_generator_settings::higgs_4loop_correction_at_as3)), "higgs4loopCorrectionAtAsAsAs");
   MLPutRuleTo(link, static_cast<int>(settings.get(Spectrum_generator_settings::loop_library)), "loopLibrary");
   MLPutRuleTo(link, settings.get(Spectrum_generator_settings::calculate_amm), "calculateAMM");
   MLPutRuleTo(link, modsel.parameter_output_scale, "parameterOutputScale");

   MLEndPacket(link);
}

/******************************************************************/

void Model_data::put_sm_input_parameters(MLINK link) const
{
   MLPutFunction(link, "List",softsusy::NUMBER_OF_LOW_ENERGY_INPUT_PARAMETERS
                              + Physical_input::NUMBER_OF_INPUT_PARAMETERS);

   MLPutRuleTo(link, qedqcd.displayAlphaEmInput(), "alphaEmMZ");
   MLPutRuleTo(link, qedqcd.displayFermiConstant(), "GF");
   MLPutRuleTo(link, qedqcd.displayAlphaSInput(), "alphaSMZ");
   MLPutRuleTo(link, qedqcd.displayPoleMZ(), "MZ");
   MLPutRuleTo(link, qedqcd.displayMbMb(), "mbmb");
   MLPutRuleTo(link, qedqcd.displayPoleMt(), "Mt");
   MLPutRuleTo(link, qedqcd.displayPoleMtau(), "Mtau");
   MLPutRuleTo(link, qedqcd.displayNeutrinoPoleMass(3), "Mv3");
   MLPutRuleTo(link, qedqcd.displayPoleMW(), "MW");
   MLPutRuleTo(link, qedqcd.displayPoleMel(), "Me");
   MLPutRuleTo(link, qedqcd.displayNeutrinoPoleMass(1), "Mv1");
   MLPutRuleTo(link, qedqcd.displayPoleMmuon(), "Mm");
   MLPutRuleTo(link, qedqcd.displayNeutrinoPoleMass(2), "Mv2");
   MLPutRuleTo(link, qedqcd.displayMd2GeV(), "md2GeV");
   MLPutRuleTo(link, qedqcd.displayMu2GeV(), "mu2GeV");
   MLPutRuleTo(link, qedqcd.displayMs2GeV(), "ms2GeV");
   MLPutRuleTo(link, qedqcd.displayMcMc(), "mcmc");

   const flexiblesusy::CKM_parameters ckm(qedqcd.displayCKM());
   MLPutRuleTo(link, ckm.theta_12, "CKMTheta12");
   MLPutRuleTo(link, ckm.theta_13, "CKMTheta13");
   MLPutRuleTo(link, ckm.theta_23, "CKMTheta23");
   MLPutRuleTo(link, ckm.delta   , "CKMDelta");

   const flexiblesusy::PMNS_parameters pmns(qedqcd.displayPMNS());
   MLPutRuleTo(link, pmns.theta_12, "PMNSTheta12");
   MLPutRuleTo(link, pmns.theta_13, "PMNSTheta13");
   MLPutRuleTo(link, pmns.theta_23, "PMNSTheta23");
   MLPutRuleTo(link, pmns.delta   , "PMNSDelta");
   MLPutRuleTo(link, pmns.alpha_1 , "PMNSAlpha1");
   MLPutRuleTo(link, pmns.alpha_2 , "PMNSAlpha2");

   MLPutRuleTo(link, physical_input.get(Physical_input::alpha_em_0), "alphaEm0");
   MLPutRuleTo(link, physical_input.get(Physical_input::mh_pole), "Mh");

   MLEndPacket(link);
}

/******************************************************************/

void Model_data::put_input_parameters(MLINK link) const
{
   MLPutFunction(link, "List", 10);

   MLPutRuleTo(link, INPUTPARAMETER(Lambda1IN), "Lambda1IN");
   MLPutRuleTo(link, INPUTPARAMETER(Lambda2IN), "Lambda2IN");
   MLPutRuleTo(link, INPUTPARAMETER(Lambda3IN), "Lambda3IN");
   MLPutRuleTo(link, INPUTPARAMETER(Lambda4IN), "Lambda4IN");
   MLPutRuleTo(link, INPUTPARAMETER(Lambda5IN), "Lambda5IN");
   MLPutRuleTo(link, INPUTPARAMETER(Lambda6IN), "Lambda6IN");
   MLPutRuleTo(link, INPUTPARAMETER(Lambda7IN), "Lambda7IN");
   MLPutRuleTo(link, INPUTPARAMETER(M122IN), "M122IN");
   MLPutRuleTo(link, INPUTPARAMETER(TanBeta), "TanBeta");
   MLPutRuleTo(link, INPUTPARAMETER(Qin), "Qin");


   MLEndPacket(link);
}

/******************************************************************/

void Model_data::put_problems(MLINK link) const
{
   check_spectrum_pointer();
   const auto problems = spectrum->get_problems();
   const auto models = problems.get_model_problems();
   const auto n_models = models.size();
   const auto solvers = problems.get_bvp_solver_problems();
   const auto n_solvers = solvers.size();
   const auto observables = spectrum->get_observables().problems;

   MLPutFunction(link, "List", n_models + n_solvers + 1);

   for (const auto& m: models) {
      MLPutRule(link, m.get_model_name());
      mathlink_format_problems(link, m);
   }

   for (const auto& m: solvers) {
      MLPutRule(link, m.get_solver_name());
      mathlink_format_problems(link, m);
   }

   MLPutRule(link, "Observables");
   mathlink_format_problems(link, observables);

   MLEndPacket(link);
}

/******************************************************************/

void Model_data::put_warnings(MLINK link) const
{
   check_spectrum_pointer();
   const auto problems = spectrum->get_problems();
   const auto models = problems.get_model_problems();
   const auto n_models = models.size();
   const auto solvers = problems.get_bvp_solver_problems();
   const auto n_solvers = solvers.size();

   MLPutFunction(link, "List", n_models + n_solvers);

   for (const auto& m: models) {
      MLPutRule(link, m.get_model_name());
      mathlink_format_warnings(link, m);
   }

   for (const auto& m: solvers) {
      MLPutRule(link, m.get_solver_name());
      mathlink_format_problems(link, m);
   }

   MLEndPacket(link);
}

/******************************************************************/

THDMII_slha_io Model_data::get_slha_io() const
{
   check_spectrum_pointer();
   THDMII_slha_io slha_io;

   slha_io.set_settings(settings);
   slha_io.set_sminputs(qedqcd);
   slha_io.set_physical_input(physical_input);
   slha_io.set_modsel(modsel);
   slha_io.set_input(input);
   
   slha_io.set_print_imaginary_parts_of_majorana_mixings(
      settings.get(Spectrum_generator_settings::force_positive_masses));

   spectrum->fill_slha_io(slha_io, settings, flexibledecay_settings);

   return slha_io;
}

/******************************************************************/

void Model_data::put_slha(MLINK link) const
{
   const auto slha_io = get_slha_io();
   std::ostringstream ostr;

   slha_io.write_to_stream(ostr);

   MLPutString(link, ostr.str().c_str());
   MLEndPacket(link);
}

/******************************************************************/

template <typename Solver_type>
void put_spectrum(const standard_model::StandardModel<Solver_type>& model, MLINK link)
{
   MLPutRule(link, standard_model_info::model_name);
   MLPutFunction(link, "List", 46);

   MLPutRuleTo(link, MODELPARAMETER(MVG), "VG", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MHp), "Hp", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MFv), "Fv", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MAh), "Ah", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(Mhh), "hh", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MFd), "Fd", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MFu), "Fu", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MFe), "Fe", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MVWp), "VWp", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MVP), "VP", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MVZ), "VZ", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(Vd), "Vd");
   MLPutRuleTo(link, MODELPARAMETER(Ud), "Ud");
   MLPutRuleTo(link, MODELPARAMETER(Vu), "Vu");
   MLPutRuleTo(link, MODELPARAMETER(Uu), "Uu");
   MLPutRuleTo(link, MODELPARAMETER(Ve), "Ve");
   MLPutRuleTo(link, MODELPARAMETER(Ue), "Ue");
   MLPutRuleTo(link, MODELPARAMETER(ZZ), "ZZ");
   MLPutRuleTo(link, PHYSICALPARAMETER(MVG), "VG", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MHp), "Hp", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MFv), "Fv", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MAh), "Ah", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(Mhh), "hh", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MFd), "Fd", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MFu), "Fu", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MFe), "Fe", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MVWp), "VWp", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MVP), "VP", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MVZ), "VZ", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(Vd), "Vd", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(Ud), "Ud", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(Vu), "Vu", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(Uu), "Uu", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(Ve), "Ve", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(Ue), "Ue", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(ZZ), "ZZ", {"Pole"});
   MLPutRuleTo(link, MODELPARAMETER(g1), "g1");
   MLPutRuleTo(link, MODELPARAMETER(g2), "g2");
   MLPutRuleTo(link, MODELPARAMETER(g3), "g3");
   MLPutRuleTo(link, MODELPARAMETER(Lambdax), "\u03bb");
   MLPutRuleTo(link, MODELPARAMETER(Yu), "Yu");
   MLPutRuleTo(link, MODELPARAMETER(Yd), "Yd");
   MLPutRuleTo(link, MODELPARAMETER(Ye), "Ye");
   MLPutRuleTo(link, MODELPARAMETER(mu2), "mu2");
   MLPutRuleTo(link, MODELPARAMETER(v), "v");
   MLPutRuleTo(link, MODELPARAMETER(scale), "SCALE");
}

/******************************************************************/

void put_spectrum(const THDMII_slha& model, MLINK link)
{
   MLPutRule(link, THDMII_info::model_name);
   MLPutFunction(link, "List", 61);

   MLPutRuleTo(link, MODELPARAMETER(MVG), "VG", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MFv), "Fv", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(Mhh), "hh", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MAh), "Ah", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MHm), "Hm", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MFd), "Fd", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MFu), "Fu", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MFe), "Fe", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MVWm), "VWm", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MVP), "VP", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(MVZ), "VZ", {"M"});
   MLPutRuleTo(link, MODELPARAMETER(ZH), "ZH");
   MLPutRuleTo(link, MODELPARAMETER(ZA), "ZA");
   MLPutRuleTo(link, MODELPARAMETER(ZP), "ZP");
   MLPutRuleTo(link, MODELPARAMETER(Vd), "Vd");
   MLPutRuleTo(link, MODELPARAMETER(Ud), "Ud");
   MLPutRuleTo(link, MODELPARAMETER(Vu), "Vu");
   MLPutRuleTo(link, MODELPARAMETER(Uu), "Uu");
   MLPutRuleTo(link, MODELPARAMETER(Ve), "Ve");
   MLPutRuleTo(link, MODELPARAMETER(Ue), "Ue");
   MLPutRuleTo(link, MODELPARAMETER(ZZ), "ZZ");
   MLPutRuleTo(link, PHYSICALPARAMETER(MVG), "VG", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MFv), "Fv", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(Mhh), "hh", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MAh), "Ah", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MHm), "Hm", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MFd), "Fd", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MFu), "Fu", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MFe), "Fe", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MVWm), "VWm", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MVP), "VP", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(MVZ), "VZ", {"Pole", "M"});
   MLPutRuleTo(link, PHYSICALPARAMETER(ZH), "ZH", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(ZA), "ZA", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(ZP), "ZP", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(Vd), "Vd", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(Ud), "Ud", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(Vu), "Vu", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(Uu), "Uu", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(Ve), "Ve", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(Ue), "Ue", {"Pole"});
   MLPutRuleTo(link, PHYSICALPARAMETER(ZZ), "ZZ", {"Pole"});
   MLPutRuleTo(link, MODELPARAMETER(g1), "g1");
   MLPutRuleTo(link, MODELPARAMETER(g2), "g2");
   MLPutRuleTo(link, MODELPARAMETER(g3), "g3");
   MLPutRuleTo(link, MODELPARAMETER(Lambda6), "Lambda6");
   MLPutRuleTo(link, MODELPARAMETER(Lambda5), "Lambda5");
   MLPutRuleTo(link, MODELPARAMETER(Lambda7), "Lambda7");
   MLPutRuleTo(link, MODELPARAMETER(Lambda1), "Lambda1");
   MLPutRuleTo(link, MODELPARAMETER(Lambda4), "Lambda4");
   MLPutRuleTo(link, MODELPARAMETER(Lambda3), "Lambda3");
   MLPutRuleTo(link, MODELPARAMETER(Lambda2), "Lambda2");
   MLPutRuleTo(link, MODELPARAMETER(Yu), "Yu");
   MLPutRuleTo(link, MODELPARAMETER(Yd), "Yd");
   MLPutRuleTo(link, MODELPARAMETER(Ye), "Ye");
   MLPutRuleTo(link, MODELPARAMETER(M122), "M122");
   MLPutRuleTo(link, MODELPARAMETER(M112), "M112");
   MLPutRuleTo(link, MODELPARAMETER(M222), "M222");
   MLPutRuleTo(link, MODELPARAMETER(v1), "v1");
   MLPutRuleTo(link, MODELPARAMETER(v2), "v2");
   MLPutRuleTo(link, MODELPARAMETER(scale), "SCALE");

}

/******************************************************************/

template<typename... Ts>
void put_spectra(const std::tuple<Ts...>& models, MLINK link)
{
   MLPutFunction(link, "List", std::tuple_size<std::tuple<Ts...>>::value);

   const auto ps = [link] (const auto& model) { put_spectrum(model, link); };
   for_each_in_tuple(models, ps);

   MLEndPacket(link);
}

/******************************************************************/

template <typename Solver_type>
void THDMII_spectrum_impl<Solver_type>::put_model_spectra(MLINK link) const
{
   put_spectra(models, link);
}

/******************************************************************/

void Model_data::put_model_spectra(MLINK link) const
{
   check_spectrum_pointer();
   spectrum->put_model_spectra(link);
}

/******************************************************************/

template <typename Solver_type>
void THDMII_spectrum_impl<Solver_type>::calculate_spectrum(
   const Spectrum_generator_settings& settings,
   const SLHA_io::Modsel& modsel,
   const softsusy::QedQcd& qedqcd,
   const THDMII_input_parameters& input
)
{
   THDMII_spectrum_generator<Solver_type> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.set_parameter_output_scale(modsel.parameter_output_scale);
   spectrum_generator.run(qedqcd, input);

   models = spectrum_generator.get_models_slha();
   problems = spectrum_generator.get_problems();
   scales.HighScale = spectrum_generator.get_high_scale();
   scales.SUSYScale = spectrum_generator.get_susy_scale();
   scales.LowScale  = spectrum_generator.get_low_scale();
   scales.pole_mass_scale = spectrum_generator.get_pole_mass_scale();
}

/******************************************************************/

template <typename Solver_type>
void THDMII_spectrum_impl<Solver_type>::calculate_model_observables(
   const softsusy::QedQcd& qedqcd,
   
   const Physical_input& physical_input,
   const Spectrum_generator_settings& settings)
{
   observables = calculate_observables(
      std::get<0>(models),
      qedqcd,
      
      physical_input,
      settings,
      scales.pole_mass_scale);
}

/******************************************************************/

template <typename Solver_type>
void THDMII_spectrum_impl<Solver_type>::calculate_unitarity() {
   unitarityData = THDMII_unitarity::max_scattering_eigenvalue_infinite_s(std::get<0>(models));
}

/******************************************************************/

template <typename Solver_type>
void THDMII_spectrum_impl<Solver_type>::fill_slha_io(THDMII_slha_io& slha_io,
       const Spectrum_generator_settings& settings, const FlexibleDecay_settings& flexibledecay_settings) const
{
   const auto& problems = std::get<0>(models).get_problems();
   const auto force_output = std::get<0>(models).do_force_output();

   slha_io.set_spinfo(problems);
   if (!problems.have_problem() || force_output) {
      slha_io.set_spectrum(models);
      slha_io.set_extra(std::get<0>(models), scales, observables, settings);
   }

   

}

/******************************************************************/

void Model_data::put_observables(MLINK link) const
{
   check_spectrum_pointer();
   THDMII_observables observables = spectrum->get_observables();

   MLPutFunction(link, "List", 1);
   MLPutRule(link, THDMII_info::model_name);
   MLPutFunction(link, "List", 4);

   MLPutRuleTo(link, OBSERVABLE(amm_Fe_1), 2, {"FlexibleSUSYObservable`AMM", "Fe"});
   MLPutRuleTo(link, OBSERVABLE(amm_uncertainty_Fe_1), 2, {"FlexibleSUSYObservable`AMMUncertainty", "Fe"});
   MLPutRuleTo(link, OBSERVABLE(a_muon_gm2calc), "FlexibleSUSYObservable`aMuonGM2Calc");
   MLPutRuleTo(link, OBSERVABLE(a_muon_gm2calc_uncertainty), "FlexibleSUSYObservable`aMuonGM2CalcUncertainty");


   MLEndPacket(link);
}


/******************************************************************/

void Model_data::put_unitarity(MLINK link) const
{
   check_spectrum_pointer();
   const UnitarityInfiniteS unitarityData = spectrum->get_unitarity();
   MLPutFunction(link, "List", 1);
   MLPutRule(link, "THDMII");
   MLPutFunction(link, "List", 3);
   MLPutRuleTo(link, unitarityData.allowed, "FlexibleSUSYUnitarity`Allowed");
   MLPutRuleTo(link, unitarityData.renScale, "FlexibleSUSYUnitarity`RenormalizationScale");
   MLPutRuleTo(link, unitarityData.maxAbsReEigenval, "FlexibleSUSYUnitarity`MaxAbsReEigen");
   MLEndPacket(link);
}

/******************************************************************/

void Model_data::check_spectrum(MLINK link) const
{
   check_spectrum_pointer();
   const auto& problems = spectrum->get_problems();

   for (const auto& s: problems.get_problem_strings())
      put_message(link, "FSTHDMIICalculateSpectrum", "error", s);

   for (const auto& s: problems.get_warning_strings())
      put_message(link, "FSTHDMIICalculateSpectrum", "warning", s);

   if (problems.have_problem() &&
       !settings.get(Spectrum_generator_settings::force_output))
      throw EInvalidSpectrum();
}

/******************************************************************/

void Model_data::calculate_spectrum()
{
   const int solver_type = static_cast<int>(settings.get(
                                               Spectrum_generator_settings::solver));
   switch (solver_type) {
   case 0:
#ifdef ENABLE_TWO_SCALE_SOLVER
   case 1:
      spectrum.reset(new THDMII_spectrum_impl<Two_scale>());
      spectrum->calculate_spectrum(settings, modsel, qedqcd, input);
      if (!spectrum->get_problems().have_problem() || solver_type != 0) break;
#endif

   default:
      if (solver_type != 0) {
         throw SetupError("invalid solver type");
      }
   }
}

/******************************************************************/

void Model_data::calculate_model_observables()
{
   check_spectrum_pointer();
   spectrum->calculate_model_observables(qedqcd,
      
      physical_input,
      settings);
}

/******************************************************************/

void Model_data::calculate_unitarity()
{
   check_spectrum_pointer();
   spectrum->calculate_unitarity();
}

/******************************************************************/

template <typename Element_t>
Model_data make_data(const Dynamic_array_view<Element_t>& pars)
{
   using Index_t = typename Dynamic_array_view<Element_t>::Index_t;

   const Index_t n_settings = Spectrum_generator_settings::NUMBER_OF_OPTIONS,
      n_sm_parameters = softsusy::NUMBER_OF_LOW_ENERGY_INPUT_PARAMETERS
                        + Physical_input::NUMBER_OF_INPUT_PARAMETERS,
      n_input_pars = 10;
   const Index_t n_fd_settings = 0;
   const Index_t n_total = n_settings + n_sm_parameters + n_input_pars + n_fd_settings;

   if (pars.size() != n_total)
      throw EWrongNumberOfParameters(pars.size(), n_total);

   Index_t c = 0; // counter

   Spectrum_generator_settings settings;
   settings.set(Spectrum_generator_settings::precision, pars[c++]);
   settings.set(Spectrum_generator_settings::max_iterations, pars[c++]);
   settings.set(Spectrum_generator_settings::solver, pars[c++]);
   settings.set(Spectrum_generator_settings::calculate_sm_masses, pars[c++]);
   settings.set(Spectrum_generator_settings::pole_mass_loop_order, pars[c++]);
   settings.set(Spectrum_generator_settings::ewsb_loop_order, pars[c++]);
   settings.set(Spectrum_generator_settings::beta_loop_order, pars[c++]);
   settings.set(Spectrum_generator_settings::threshold_corrections_loop_order, pars[c++]);
   settings.set(Spectrum_generator_settings::higgs_2loop_correction_at_as, pars[c++]);
   settings.set(Spectrum_generator_settings::higgs_2loop_correction_ab_as, pars[c++]);
   settings.set(Spectrum_generator_settings::higgs_2loop_correction_at_at, pars[c++]);
   settings.set(Spectrum_generator_settings::higgs_2loop_correction_atau_atau, pars[c++]);
   settings.set(Spectrum_generator_settings::force_output, pars[c++]);
   settings.set(Spectrum_generator_settings::top_pole_qcd_corrections, pars[c++]);
   settings.set(Spectrum_generator_settings::beta_zero_threshold, pars[c++]);
   settings.set(Spectrum_generator_settings::force_positive_masses, pars[c++]);
   settings.set(Spectrum_generator_settings::pole_mass_scale, pars[c++]);
   settings.set(Spectrum_generator_settings::eft_pole_mass_scale, pars[c++]);
   settings.set(Spectrum_generator_settings::eft_matching_scale, pars[c++]);
   settings.set(Spectrum_generator_settings::eft_matching_loop_order_up, pars[c++]);
   settings.set(Spectrum_generator_settings::eft_matching_loop_order_down, pars[c++]);
   settings.set(Spectrum_generator_settings::eft_higgs_index, pars[c++]);
   settings.set(Spectrum_generator_settings::calculate_bsm_masses, pars[c++]);
   settings.set(Spectrum_generator_settings::threshold_corrections, pars[c++]);
   settings.set(Spectrum_generator_settings::higgs_3loop_ren_scheme_atb_as2, pars[c++]);
   settings.set(Spectrum_generator_settings::higgs_3loop_correction_at_as2, pars[c++]);
   settings.set(Spectrum_generator_settings::higgs_3loop_correction_ab_as2, pars[c++]);
   settings.set(Spectrum_generator_settings::higgs_3loop_correction_at2_as, pars[c++]);
   settings.set(Spectrum_generator_settings::higgs_3loop_correction_at3, pars[c++]);
   settings.set(Spectrum_generator_settings::higgs_4loop_correction_at_as3, pars[c++]);
   settings.set(Spectrum_generator_settings::loop_library, pars[c++]);
   settings.set(Spectrum_generator_settings::calculate_amm, pars[c++]);

   SLHA_io::Modsel modsel;
   modsel.parameter_output_scale = pars[c++];

   softsusy::QedQcd qedqcd;
   qedqcd.setAlpha(softsusy::ALPHA, pars[c]);
   qedqcd.setAlphaEmInput(pars[c++]);
   qedqcd.setFermiConstant(pars[c++]);
   qedqcd.setAlpha(softsusy::ALPHAS, pars[c]);
   qedqcd.setAlphaSInput(pars[c++]);
   qedqcd.setPoleMZ(pars[c]);
   qedqcd.set_scale(pars[c++]);
   qedqcd.setMass(softsusy::mBottom, pars[c]);
   qedqcd.setMbMb(pars[c++]);
   qedqcd.setPoleMt(pars[c++]);
   qedqcd.setMass(softsusy::mTau, pars[c]);
   qedqcd.setPoleMtau(pars[c++]);
   qedqcd.setNeutrinoPoleMass(3, pars[c++]);
   qedqcd.setPoleMW(pars[c++]);
   qedqcd.setMass(softsusy::mElectron, pars[c]);
   qedqcd.setPoleMel(pars[c++]);
   qedqcd.setNeutrinoPoleMass(1, pars[c++]);
   qedqcd.setMass(softsusy::mMuon, pars[c]);
   qedqcd.setPoleMmuon(pars[c++]);
   qedqcd.setNeutrinoPoleMass(2, pars[c++]);
   qedqcd.setMass(softsusy::mDown, pars[c]);
   qedqcd.setMd2GeV(pars[c++]);
   qedqcd.setMass(softsusy::mUp, pars[c]);
   qedqcd.setMu2GeV(pars[c++]);
   qedqcd.setMass(softsusy::mStrange, pars[c]);
   qedqcd.setMs2GeV(pars[c++]);
   qedqcd.setMass(softsusy::mCharm, pars[c]);
   qedqcd.setMcMc(pars[c++]);

   {
      flexiblesusy::CKM_parameters ckm;
      ckm.theta_12 = pars[c++];
      ckm.theta_13 = pars[c++];
      ckm.theta_23 = pars[c++];
      ckm.delta    = pars[c++];
      qedqcd.setCKM(ckm);
   }

   {
      flexiblesusy::PMNS_parameters pmns;
      pmns.theta_12 = pars[c++];
      pmns.theta_13 = pars[c++];
      pmns.theta_23 = pars[c++];
      pmns.delta    = pars[c++];
      pmns.alpha_1  = pars[c++];
      pmns.alpha_2  = pars[c++];
      qedqcd.setPMNS(pmns);
   }

   
   

   Physical_input physical_input;
   physical_input.set(Physical_input::alpha_em_0, pars[c++]);
   physical_input.set(Physical_input::mh_pole, pars[c++]);

   THDMII_input_parameters input;
   INPUTPARAMETER(Lambda1IN) = pars[c++];
   INPUTPARAMETER(Lambda2IN) = pars[c++];
   INPUTPARAMETER(Lambda3IN) = pars[c++];
   INPUTPARAMETER(Lambda4IN) = pars[c++];
   INPUTPARAMETER(Lambda5IN) = pars[c++];
   INPUTPARAMETER(Lambda6IN) = pars[c++];
   INPUTPARAMETER(Lambda7IN) = pars[c++];
   INPUTPARAMETER(M122IN) = pars[c++];
   INPUTPARAMETER(TanBeta) = pars[c++];
   INPUTPARAMETER(Qin) = pars[c++];


   Model_data data;
   data.set_settings(settings);
   data.set_modsel(modsel);
   data.set_sm_input_parameters(qedqcd);
   data.set_physical_input(physical_input);
   data.set_input_parameters(input);
   
   
   return data;
}

} // namespace THDMII_librarylink
} // namespace flexiblesusy

extern "C" {

/******************************************************************/

DLLEXPORT mint WolframLibrary_getVersion()
{
   return WolframLibraryVersion;
}

/******************************************************************/

DLLEXPORT int WolframLibrary_initialize(WolframLibraryData /* libData */)
{
   return LIBRARY_NO_ERROR;
}

/******************************************************************/

DLLEXPORT int FSTHDMIIGetSettings(WolframLibraryData /* libData */, MLINK link)
{
   using namespace flexiblesusy::THDMII_librarylink;

   if (!check_number_of_args(link, 1, "FSTHDMIIGetSettings"))
      return LIBRARY_TYPE_ERROR;

   const auto hid = get_handle_from(link);

   try {
      find_data(hid).put_settings(link);
   } catch (const flexiblesusy::Error& e) {
      std::cerr << e.what_detailed() << std::endl;
      put_error_output(link);
   }

   return LIBRARY_NO_ERROR;
}

/******************************************************************/

DLLEXPORT int FSTHDMIIGetSMInputParameters(WolframLibraryData /* libData */, MLINK link)
{
   using namespace flexiblesusy::THDMII_librarylink;

   if (!check_number_of_args(link, 1, "FSTHDMIIGetSMInputParameters"))
      return LIBRARY_TYPE_ERROR;

   const auto hid = get_handle_from(link);

   try {
      find_data(hid).put_sm_input_parameters(link);
   } catch (const flexiblesusy::Error& e) {
      std::cerr << e.what_detailed() << std::endl;
      put_error_output(link);
   }

   return LIBRARY_NO_ERROR;
}

/******************************************************************/

DLLEXPORT int FSTHDMIIGetInputParameters(WolframLibraryData /* libData */, MLINK link)
{
   using namespace flexiblesusy::THDMII_librarylink;

   if (!check_number_of_args(link, 1, "FSTHDMIIGetInputParameters"))
      return LIBRARY_TYPE_ERROR;

   const auto hid = get_handle_from(link);

   try {
      find_data(hid).put_input_parameters(link);
   } catch (const flexiblesusy::Error& e) {
      std::cerr << e.what_detailed() << std::endl;
      put_error_output(link);
   }

   return LIBRARY_NO_ERROR;
}

/******************************************************************/

DLLEXPORT int FSTHDMIIOpenHandle(
   WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
   using namespace flexiblesusy;
   using namespace flexiblesusy::THDMII_librarylink;

   if (Argc != 1)
      return LIBRARY_TYPE_ERROR;

   MTensor pars = MArgument_getMTensor(Args[0]);

   if (libData->MTensor_getType(pars) != MType_Real ||
       libData->MTensor_getRank(pars) != 1)
      return LIBRARY_TYPE_ERROR;

   try {
      auto data = make_data(make_dynamic_array_view(
                          libData->MTensor_getRealData(pars),
                          libData->MTensor_getDimensions(pars)[0]));

      const auto hid = get_new_handle();

      handles.insert(std::make_pair(hid, std::move(data)));

      MArgument_setInteger(Res, hid);
   } catch (const flexiblesusy::Error& e) {
      std::cerr << e.what_detailed() << std::endl;
      return LIBRARY_FUNCTION_ERROR;
   }

   return LIBRARY_NO_ERROR;
}

/******************************************************************/

DLLEXPORT int FSTHDMIICloseHandle(
   WolframLibraryData /* libData */, mint Argc, MArgument* Args, MArgument /* Res */)
{
   using namespace flexiblesusy::THDMII_librarylink;

   if (Argc != 1)
      return LIBRARY_TYPE_ERROR;

   const auto hid = get_handle_from(Args[0]);
   const auto handle = handles.find(hid);

   if (handle != handles.end())
      handles.erase(handle);

   return LIBRARY_NO_ERROR;
}

/******************************************************************/

DLLEXPORT int FSTHDMIISet(
   WolframLibraryData libData, mint Argc, MArgument* Args, MArgument /* Res */)
{
   using namespace flexiblesusy;
   using namespace flexiblesusy::THDMII_librarylink;

   if (Argc != 2)
      return LIBRARY_TYPE_ERROR;

   const auto hid = get_handle_from(Args[0]);
   MTensor pars = MArgument_getMTensor(Args[1]);

   if (libData->MTensor_getType(pars) != MType_Real ||
       libData->MTensor_getRank(pars) != 1)
      return LIBRARY_TYPE_ERROR;

   try {
      find_data(hid) =
         make_data(make_dynamic_array_view(
                      libData->MTensor_getRealData(pars),
                      libData->MTensor_getDimensions(pars)[0]));
   } catch (const flexiblesusy::Error& e) {
      std::cerr << e.what_detailed() << std::endl;
      return LIBRARY_FUNCTION_ERROR;
   }

   return LIBRARY_NO_ERROR;
}

/******************************************************************/

DLLEXPORT int FSTHDMIIGetProblems(
   WolframLibraryData /* libData */, MLINK link)
{
   using namespace flexiblesusy::THDMII_librarylink;

   if (!check_number_of_args(link, 1, "FSTHDMIIGetProblems"))
      return LIBRARY_TYPE_ERROR;

   const auto hid = get_handle_from(link);

   try {
      find_data(hid).put_problems(link);
   } catch (const flexiblesusy::Error& e) {
      std::cerr << e.what_detailed() << std::endl;
      put_error_output(link);
   }

   return LIBRARY_NO_ERROR;
}

/******************************************************************/

DLLEXPORT int FSTHDMIIToSLHA(WolframLibraryData /* libData */, MLINK link)
{
   using namespace flexiblesusy::THDMII_librarylink;

   if (!check_number_of_args(link, 1, "FSTHDMIIToSLHA"))
      return LIBRARY_TYPE_ERROR;

   const auto hid = get_handle_from(link);

   try {
      find_data(hid).put_slha(link);
   } catch (const flexiblesusy::Error& e) {
      std::cerr << e.what_detailed() << std::endl;
      put_error_output(link);
   }

   return LIBRARY_NO_ERROR;
}

/******************************************************************/

DLLEXPORT int FSTHDMIIGetWarnings(
   WolframLibraryData /* libData */, MLINK link)
{
   using namespace flexiblesusy::THDMII_librarylink;

   if (!check_number_of_args(link, 1, "FSTHDMIIGetWarnings"))
      return LIBRARY_TYPE_ERROR;

   const auto hid = get_handle_from(link);

   try {
      find_data(hid).put_warnings(link);
   } catch (const flexiblesusy::Error& e) {
      std::cerr << e.what_detailed() << std::endl;
      put_error_output(link);
   }

   return LIBRARY_NO_ERROR;
}

/******************************************************************/

DLLEXPORT int FSTHDMIICalculateSpectrum(
   WolframLibraryData /* libData */, MLINK link)
{
   using namespace flexiblesusy::THDMII_librarylink;

   if (!check_number_of_args(link, 1, "FSTHDMIICalculateSpectrum"))
      return LIBRARY_TYPE_ERROR;

   const auto hid = get_handle_from(link);

   try {
      auto& data = find_data(hid);

      {
         Redirect_output crd(link);
         data.calculate_spectrum();
      }

      data.check_spectrum(link);
      data.put_model_spectra(link);
   } catch (const flexiblesusy::Error&) {
      put_error_output(link);
   }

   return LIBRARY_NO_ERROR;
}

/******************************************************************/

DLLEXPORT int FSTHDMIICalculateObservables(
   WolframLibraryData /* libData */, MLINK link)
{
   using namespace flexiblesusy::THDMII_librarylink;

   if (!check_number_of_args(link, 1, "FSTHDMIICalculateObservables"))
      return LIBRARY_TYPE_ERROR;

   const auto hid = get_handle_from(link);

   try {
      auto& data = find_data(hid);

      if (data.get_model_scale() == 0.) {
         put_message(link,
            "FSTHDMIICalculateObservables", "warning",
            "Renormalization scale is 0.  Did you run "
            "FSTHDMIICalculateSpectrum[]?");
      }

      {
         Redirect_output crd(link);
         data.calculate_model_observables();
         auto setting = data.get_settings();
         setting.set(flexiblesusy::Spectrum_generator_settings::calculate_observables, 1.0);
         data.set_settings(setting);
      }

      data.put_observables(link);
   } catch (const flexiblesusy::Error& e) {
      put_message(link, "FSTHDMIICalculateObservables", "error", e.what_detailed());
      put_error_output(link);
   }

   return LIBRARY_NO_ERROR;
}


DLLEXPORT int FSTHDMIICalculateUnitarity(
   WolframLibraryData /* libData */, MLINK link)
{
   using namespace flexiblesusy::THDMII_librarylink;

   if (!check_number_of_args(link, 1, "FSTHDMIICalculateUnitarity"))
      return LIBRARY_TYPE_ERROR;

   const auto hid = get_handle_from(link);

   try {
      auto& data = find_data(hid);
      data.calculate_unitarity();
      data.put_unitarity(link);
   } catch (const flexiblesusy::Error& e) {
      put_message(link, "FSTHDMIICalculateUnitarity", "error", e.what());
      put_error_output(link);
   }

   return LIBRARY_NO_ERROR;
}

} // extern "C"
