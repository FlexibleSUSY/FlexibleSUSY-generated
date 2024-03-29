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


#include "HGTHDMIIMSSMBC_slha_io.hpp"
#include "HGTHDMIIMSSMBC_input_parameters.hpp"
#include "HGTHDMIIMSSMBC_mass_eigenstates.hpp"
#include "HGTHDMIIMSSMBC_model_slha.hpp"
#include "HGTHDMIIMSSMBC_observables.hpp"
#include "HGTHDMIIMSSMBC_physical.hpp"
#include "ew_input.hpp"
#include "logger.hpp"
#include "observable_problems.hpp"
#include "observable_problems_format_slha.hpp"
#include "numerics2.hpp"
#include "spectrum_generator_problems.hpp"
#include "standard_model.hpp"
#include "wrappers.hpp"
#include "config.h"
#include "spectrum_generator_settings.hpp"

#include <array>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iterator>
#include <string>

#define Pole(p) physical.p
#define PHYSICAL(p) model.get_physical().p
#define PHYSICAL_SLHA(p) model.get_physical_slha().p
#define LOCALPHYSICAL(p) physical.p
#define MODEL model
#define MODELPARAMETER(p) model.get_##p()
#define INPUTPARAMETER(p) input.p
#define EXTRAPARAMETER(p) model.get_##p()
#define OBSERVABLES observables
#define DEFINE_PHYSICAL_PARAMETER(p) decltype(LOCALPHYSICAL(p)) p;
#define LowEnergyConstant(p) Electroweak_constants::p
#define SCALES(p) scales.p

namespace flexiblesusy {

HGTHDMIIMSSMBC_slha_io::HGTHDMIIMSSMBC_slha_io()
   : slha_io()
   , print_imaginary_parts_of_majorana_mixings(false)
{
}

void HGTHDMIIMSSMBC_slha_io::clear()
{
   slha_io.clear();
}

void HGTHDMIIMSSMBC_slha_io::set_print_imaginary_parts_of_majorana_mixings(bool flag)
{
   print_imaginary_parts_of_majorana_mixings = flag;
}

/**
 * Reads DR-bar parameters, pole masses and mixing matrices from a
 * SLHA output file.
 */
void HGTHDMIIMSSMBC_slha_io::fill(HGTHDMIIMSSMBC_slha& model) const
{
   fill(static_cast<HGTHDMIIMSSMBC_mass_eigenstates&>(model));
   fill_physical(model.get_physical_slha());
}

/**
 * Stores the EXTPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void HGTHDMIIMSSMBC_slha_io::set_extpar(const HGTHDMIIMSSMBC_input_parameters& input)
{
   std::ostringstream extpar;

   extpar << "Block EXTPAR\n";
   extpar << FORMAT_ELEMENT(0, input.MSUSY, "MSUSY");
   extpar << FORMAT_ELEMENT(1, input.MEWSB, "MEWSB");
   extpar << FORMAT_ELEMENT(2, input.MuInput, "MuInput");
   extpar << FORMAT_ELEMENT(3, input.M1Input, "M1Input");
   extpar << FORMAT_ELEMENT(4, input.M2Input, "M2Input");
   extpar << FORMAT_ELEMENT(5, input.M3Input, "M3Input");
   extpar << FORMAT_ELEMENT(6, input.MAInput, "MAInput");
   extpar << FORMAT_ELEMENT(100, input.LambdaLoopOrder, "LambdaLoopOrder");
   slha_io.set_block(extpar);

}

/**
 * Stores the IMMINPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void HGTHDMIIMSSMBC_slha_io::set_imminpar(const HGTHDMIIMSSMBC_input_parameters& input)
{

}

/**
 * Stores the IMEXTPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void HGTHDMIIMSSMBC_slha_io::set_imextpar(const HGTHDMIIMSSMBC_input_parameters& input)
{

}

/**
 * Stores the MODSEL input parameters in the SLHA object.
 *
 * @param modsel struct of MODSEL parameters
 */
void HGTHDMIIMSSMBC_slha_io::set_modsel(const SLHA_io::Modsel& modsel)
{
   slha_io.set_modsel(modsel);
}

/**
 * Stores the MINPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void HGTHDMIIMSSMBC_slha_io::set_minpar(const HGTHDMIIMSSMBC_input_parameters& input)
{
   std::ostringstream minpar;

   minpar << "Block MINPAR\n";
   minpar << FORMAT_ELEMENT(3, input.TanBeta, "TanBeta");
   slha_io.set_block(minpar);

}

/**
 * Stores all input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void HGTHDMIIMSSMBC_slha_io::set_input(const HGTHDMIIMSSMBC_input_parameters& input)
{
   set_minpar(input);
   set_extpar(input);
   set_imminpar(input);
   set_imextpar(input);

   slha_io.set_block("AdIN", INPUTPARAMETER(AdInput), "AdInput");
   slha_io.set_block("AeIN", INPUTPARAMETER(AeInput), "AeInput");
   slha_io.set_block("AuIN", INPUTPARAMETER(AuInput), "AuInput");
   slha_io.set_block("MSDIN", INPUTPARAMETER(msdInput), "msdInput");
   slha_io.set_block("MSEIN", INPUTPARAMETER(mseInput), "mseInput");
   slha_io.set_block("MSLIN", INPUTPARAMETER(mslInput), "mslInput");
   slha_io.set_block("MSQIN", INPUTPARAMETER(msqInput), "msqInput");
   slha_io.set_block("MSUIN", INPUTPARAMETER(msuInput), "msuInput");

}

/**
 * Stores the additional physical input (FlexibleSUSYInput block) in
 * the SLHA object.
 *
 * @param input class of input
 */
void HGTHDMIIMSSMBC_slha_io::set_physical_input(const Physical_input& input)
{
   slha_io.set_physical_input(input);
}

/**
 * Stores the settings (FlexibleSUSY block) in the SLHA object.
 *
 * @param settings class of settings
 */
void HGTHDMIIMSSMBC_slha_io::set_settings(const Spectrum_generator_settings& settings)
{
   slha_io.set_settings(settings);
}

/**
 * Stores the settings (LToLConversion block) in the SLHA object.
 *
 * @param settings class of settings
 */
void HGTHDMIIMSSMBC_slha_io::set_LToLConversion_settings(const LToLConversion_settings& settings)
{
   slha_io.set_LToLConversion_settings(settings);
}

/**
 * Stores the settings (FlexibleSUSY block) in the SLHA object.
 *
 * @param settings class of settings
 */
void HGTHDMIIMSSMBC_slha_io::set_FlexibleDecay_settings(const FlexibleDecay_settings& settings)
{
   slha_io.set_FlexibleDecay_settings(settings);
}

/**
 * Stores the settings (FlexibleSUSYUnitarity block) in the SLHA object.
 */
void HGTHDMIIMSSMBC_slha_io::set_unitarity_infinite_s(
   const flexiblesusy::Spectrum_generator_settings& spectrum_generator_settings, UnitarityInfiniteS const& unitarity)
{
   slha_io.set_unitarity_infinite_s(spectrum_generator_settings, unitarity);
}

/**
 * Stores the SMINPUTS input parameters in the SLHA object.
 *
 * @param qedqcd class of Standard Model parameters
 */
void HGTHDMIIMSSMBC_slha_io::set_sminputs(const softsusy::QedQcd& qedqcd)
{
   slha_io.set_sminputs(qedqcd);
}

/**
 * Stores the spectrum generator information in the SPINFO block in
 * the SLHA object.
 *
 * @param problems struct with parameter point problems
 */
void HGTHDMIIMSSMBC_slha_io::set_spinfo(const Spectrum_generator_problems& problems)
{
   set_spinfo(problems.get_problem_strings(), problems.get_warning_strings());
}

/**
 * Stores the spectrum generator information in the SPINFO block in
 * the SLHA object.
 *
 * @param problems struct with parameter point problems
 */
void HGTHDMIIMSSMBC_slha_io::set_spinfo(const Problems& problems)
{
   set_spinfo(problems.get_problem_strings(), problems.get_warning_strings());
}

/**
 * Stores the given problems and warnings in the SPINFO block in the
 * SLHA object.
 *
 * @param problems vector of problem strings
 * @param warnings vector of warning strings
 */
void HGTHDMIIMSSMBC_slha_io::set_spinfo(
   const std::vector<std::string>& problems,
   const std::vector<std::string>& warnings)
{
   std::ostringstream spinfo;
   spinfo << "Block SPINFO\n"
          << FORMAT_SPINFO(1, PKGNAME)
          << FORMAT_SPINFO(2, FLEXIBLESUSY_VERSION);

   for (const auto& s: warnings)
      spinfo << FORMAT_SPINFO(3, s);

   for (const auto& s: problems)
      spinfo << FORMAT_SPINFO(4, s);

   spinfo << FORMAT_SPINFO(5, HGTHDMIIMSSMBC_info::model_name)
          << FORMAT_SPINFO(9, SARAH_VERSION);

   slha_io.set_block(spinfo, SLHA_io::front);
}

/**
 * Stores the particle masses in the SLHA object.
 *
 * @param physical struct of physical parameters
 *
 * @param write_sm_masses flag to indicate if Standard Model
 *    particle masses should be written as well
 */
void HGTHDMIIMSSMBC_slha_io::set_mass(const HGTHDMIIMSSMBC_physical& physical,
                                   bool write_sm_masses)
{
   std::ostringstream mass;

   mass << "Block MASS\n"
      << FORMAT_MASS(24, LOCALPHYSICAL(MVWm), "VWm")
      << FORMAT_MASS(1000021, LOCALPHYSICAL(MGlu), "Glu")
      << FORMAT_MASS(1000024, LOCALPHYSICAL(MCha(0)), "Cha(1)")
      << FORMAT_MASS(1000037, LOCALPHYSICAL(MCha(1)), "Cha(2)")
      << FORMAT_MASS(37, LOCALPHYSICAL(MHm(1)), "Hm(2)")
      << FORMAT_MASS(25, LOCALPHYSICAL(Mhh(0)), "hh(1)")
      << FORMAT_MASS(35, LOCALPHYSICAL(Mhh(1)), "hh(2)")
      << FORMAT_MASS(36, LOCALPHYSICAL(MAh(1)), "Ah(2)")
      << FORMAT_MASS(1000022, LOCALPHYSICAL(MChi(0)), "Chi(1)")
      << FORMAT_MASS(1000023, LOCALPHYSICAL(MChi(1)), "Chi(2)")
      << FORMAT_MASS(1000025, LOCALPHYSICAL(MChi(2)), "Chi(3)")
      << FORMAT_MASS(1000035, LOCALPHYSICAL(MChi(3)), "Chi(4)")
   ;

   if (write_sm_masses) {
      mass
         << FORMAT_MASS(21, LOCALPHYSICAL(MVG), "VG")
         << FORMAT_MASS(12, LOCALPHYSICAL(MFv(0)), "Fv(1)")
         << FORMAT_MASS(14, LOCALPHYSICAL(MFv(1)), "Fv(2)")
         << FORMAT_MASS(16, LOCALPHYSICAL(MFv(2)), "Fv(3)")
         << FORMAT_MASS(1, LOCALPHYSICAL(MFd(0)), "Fd(1)")
         << FORMAT_MASS(3, LOCALPHYSICAL(MFd(1)), "Fd(2)")
         << FORMAT_MASS(5, LOCALPHYSICAL(MFd(2)), "Fd(3)")
         << FORMAT_MASS(2, LOCALPHYSICAL(MFu(0)), "Fu(1)")
         << FORMAT_MASS(4, LOCALPHYSICAL(MFu(1)), "Fu(2)")
         << FORMAT_MASS(6, LOCALPHYSICAL(MFu(2)), "Fu(3)")
         << FORMAT_MASS(11, LOCALPHYSICAL(MFe(0)), "Fe(1)")
         << FORMAT_MASS(13, LOCALPHYSICAL(MFe(1)), "Fe(2)")
         << FORMAT_MASS(15, LOCALPHYSICAL(MFe(2)), "Fe(3)")
         << FORMAT_MASS(22, LOCALPHYSICAL(MVP), "VP")
         << FORMAT_MASS(23, LOCALPHYSICAL(MVZ), "VZ")
      ;
   }

   slha_io.set_block(mass);

}

/**
 * Stores the mixing matrices in the SLHA object.
 *
 * @param physical struct of physical parameters
 *
 * @param write_sm_mixing_matrics flag to indicate if Standard Model
 *    particle mixing matrices should be written as well
 */
void HGTHDMIIMSSMBC_slha_io::set_mixing_matrices(const HGTHDMIIMSSMBC_physical& physical,
                                              bool write_sm_mixing_matrics)
{
   slha_io.set_block("UMIX", LOCALPHYSICAL(UM), "UM");
   slha_io.set_block("VMIX", LOCALPHYSICAL(UP), "UP");
   slha_io.set_block("PSEUDOSCALARMIX", LOCALPHYSICAL(ZA), "ZA");
   slha_io.set_block("SCALARMIX", LOCALPHYSICAL(ZH), "ZH");
   slha_io.set_block("NMIX", LOCALPHYSICAL(ZN), "ZN");
   slha_io.set_block("CHARGEMIX", LOCALPHYSICAL(ZP), "ZP");

   if (write_sm_mixing_matrics) {
      slha_io.set_block("UULMIX", LOCALPHYSICAL(Vu), "Vu");
      slha_io.set_block("UDLMIX", LOCALPHYSICAL(Vd), "Vd");
      slha_io.set_block("UURMIX", LOCALPHYSICAL(Uu), "Uu");
      slha_io.set_block("UDRMIX", LOCALPHYSICAL(Ud), "Ud");
      slha_io.set_block("UELMIX", LOCALPHYSICAL(Ve), "Ve");
      slha_io.set_block("UERMIX", LOCALPHYSICAL(Ue), "Ue");
   }

   if (print_imaginary_parts_of_majorana_mixings) {
      slha_io.set_block_imag("IMNMIX", LOCALPHYSICAL(ZN), "ZN");
   }

}

void HGTHDMIIMSSMBC_slha_io::set_ckm(
   const Eigen::Matrix<std::complex<double>,3,3>& ckm_matrix,
   double scale)
{
   slha_io.set_block("VCKM"  , ckm_matrix.real(), "Re(CKM)", scale);
   slha_io.set_block("IMVCKM", ckm_matrix.imag(), "Im(CKM)", scale);
}

void HGTHDMIIMSSMBC_slha_io::set_pmns(
   const Eigen::Matrix<std::complex<double>,3,3>& pmns_matrix,
   double scale)
{
   slha_io.set_block("VPMNS"  , pmns_matrix.real(), "Re(PMNS)", scale);
   slha_io.set_block("IMVPMNS", pmns_matrix.imag(), "Im(PMNS)", scale);
}

/**
 * Stores the model (DR-bar) parameters in the SLHA object.
 *
 * @param model model class
 */
void HGTHDMIIMSSMBC_slha_io::set_model_parameters(const HGTHDMIIMSSMBC_slha& model)
{
   {
      std::ostringstream block;
      block << "Block gauge Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(g1) * 0.7745966692414834), "g1 * 0.7745966692414834")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(g2)), "g2")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(g3)), "g3")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("Yu", ToMatrix(MODELPARAMETER(Yu_slha)), "Yu", model.get_scale());
   slha_io.set_block("Yd", ToMatrix(MODELPARAMETER(Yd_slha)), "Yd", model.get_scale());
   slha_io.set_block("Ye", ToMatrix(MODELPARAMETER(Ye_slha)), "Ye", model.get_scale());
   {
      std::ostringstream block;
      block << "Block HMIX Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(31, (MODELPARAMETER(Lambda1)), "Lambda1")
            << FORMAT_ELEMENT(32, (MODELPARAMETER(Lambda2)), "Lambda2")
            << FORMAT_ELEMENT(33, (MODELPARAMETER(Lambda3)), "Lambda3")
            << FORMAT_ELEMENT(34, (MODELPARAMETER(Lambda4)), "Lambda4")
            << FORMAT_ELEMENT(35, (MODELPARAMETER(Lambda5)), "Lambda5")
            << FORMAT_ELEMENT(36, (MODELPARAMETER(Lambda6)), "Lambda6")
            << FORMAT_ELEMENT(37, (MODELPARAMETER(Lambda7)), "Lambda7")
            << FORMAT_ELEMENT(20, (MODELPARAMETER(M112)), "M112")
            << FORMAT_ELEMENT(21, (MODELPARAMETER(M222)), "M222")
            << FORMAT_ELEMENT(22, (MODELPARAMETER(M122)), "M122")
            << FORMAT_ELEMENT(102, (MODELPARAMETER(v1)), "v1")
            << FORMAT_ELEMENT(103, (MODELPARAMETER(v2)), "v2")
            << FORMAT_ELEMENT(1, (MODELPARAMETER(Mu)), "Mu")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block MSOFT Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(MassB)), "MassB")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(MassWB)), "MassWB")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(MassG)), "MassG")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block SplitMSSM Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(g2up)), "g2up")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(g2u)), "g2u")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(g1dp)), "g1dp")
            << FORMAT_ELEMENT(4, (MODELPARAMETER(g1d)), "g1d")
      ;
      slha_io.set_block(block);
   }

   {
      std::ostringstream block;
      block << "Block Phases Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (Re(MODELPARAMETER(PhaseGlu))), "Re(PhaseGlu)")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block IMPhases Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (Im(MODELPARAMETER(PhaseGlu))), "Im(PhaseGlu)")
      ;
      slha_io.set_block(block);
   }

}

void HGTHDMIIMSSMBC_slha_io::set_model_parameters(const standard_model::Standard_model& model)
{
   {
      std::ostringstream block;
      block << "Block SMGAUGE Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (model.get_g1() * standard_model_info::normalization_g1), "gY")
            << FORMAT_ELEMENT(2, (model.get_g2()), "g2")
            << FORMAT_ELEMENT(3, (model.get_g3()), "g3")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("SMYu", ToMatrix(model.get_Yu()), "Yu", model.get_scale());
   slha_io.set_block("SMYd", ToMatrix(model.get_Yd()), "Yd", model.get_scale());
   slha_io.set_block("SMYe", ToMatrix(model.get_Ye()), "Ye", model.get_scale());
   {
      std::ostringstream block;
      block << "Block SMSM Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (model.get_mu2()), "mu2")
            << FORMAT_ELEMENT(2, (model.get_Lambdax()), "Lambdax")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block SMHMIX Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(3, (model.get_v()), "v")
      ;
      slha_io.set_block(block);
   }
}

void HGTHDMIIMSSMBC_slha_io::set_mass(const standard_model::Standard_model_physical& physical)
{
   std::ostringstream mass;

   mass << "Block SMMASS\n"
      << FORMAT_MASS(24, physical.MVWp, "VWp")
      << FORMAT_MASS(21, physical.MVG, "VG")
      << FORMAT_MASS(12, physical.MFv(0), "Fv(1)")
      << FORMAT_MASS(14, physical.MFv(1), "Fv(2)")
      << FORMAT_MASS(16, physical.MFv(2), "Fv(3)")
      << FORMAT_MASS(25, physical.Mhh, "hh")
      << FORMAT_MASS(1, physical.MFd(0), "Fd(1)")
      << FORMAT_MASS(3, physical.MFd(1), "Fd(2)")
      << FORMAT_MASS(5, physical.MFd(2), "Fd(3)")
      << FORMAT_MASS(2, physical.MFu(0), "Fu(1)")
      << FORMAT_MASS(4, physical.MFu(1), "Fu(2)")
      << FORMAT_MASS(6, physical.MFu(2), "Fu(3)")
      << FORMAT_MASS(11, physical.MFe(0), "Fe(1)")
      << FORMAT_MASS(13, physical.MFe(1), "Fe(2)")
      << FORMAT_MASS(15, physical.MFe(2), "Fe(3)")
      << FORMAT_MASS(22, physical.MVP, "VP")
      << FORMAT_MASS(23, physical.MVZ, "VZ")
      ;

   slha_io.set_block(mass);
}

void HGTHDMIIMSSMBC_slha_io::set_mixing_matrices(const standard_model::Standard_model_physical& physical)
{
   slha_io.set_block("SMUULMIX", physical.Vu, "Vu");
   slha_io.set_block("SMUDLMIX", physical.Vd, "Vd");
   slha_io.set_block("SMUURMIX", physical.Uu, "Uu");
   slha_io.set_block("SMUDRMIX", physical.Ud, "Ud");
   slha_io.set_block("SMUELMIX", physical.Ve, "Ve");
   slha_io.set_block("SMUERMIX", physical.Ue, "Ue");
}

void HGTHDMIIMSSMBC_slha_io::set_spectrum(const standard_model::Standard_model& model)
{
   const auto& physical = model.get_physical();

   set_model_parameters(model);
   set_mass(physical);
   set_mixing_matrices(physical);
}

/**
 * Stores the decays calculation information in the DCINFO block
 * in the SLHA object.
 *
 * @param problems struct with decays calculation problems
 */
void HGTHDMIIMSSMBC_slha_io::set_dcinfo(const FlexibleDecay_problems& problems)
{
   set_dcinfo(problems.get_problem_strings(), problems.get_warning_strings());
}

/**
 * Stores given problems and warnings in the DCINFO block in the SLHA
 * object.
 *
 * @param problems vector of problem strings
 * @param warnings vector of warning strings
 */
void HGTHDMIIMSSMBC_slha_io::set_dcinfo(
   const std::vector<std::string>& problems,
   const std::vector<std::string>& warnings)
{
   std::ostringstream dcinfo;
   dcinfo << "Block DCINFO\n"
          << FORMAT_SPINFO(1, PKGNAME)
          << FORMAT_SPINFO(2, FLEXIBLESUSY_VERSION);

   for (const auto& s: warnings)
      dcinfo << FORMAT_SPINFO(3, s);

   for (const auto& s: problems)
      dcinfo << FORMAT_SPINFO(4, s);

   dcinfo << FORMAT_SPINFO(5, HGTHDMIIMSSMBC_info::model_name)
          << FORMAT_SPINFO(9, SARAH_VERSION);

   slha_io.set_block(dcinfo);
}

/**
 * Stores the branching ratios for a given particle in the SLHA
 * object.
 *
 * @param decays struct containing individual particle decays
 */
void HGTHDMIIMSSMBC_slha_io::set_decay_block(const Decays_list& decays_list, FlexibleDecay_settings const& flexibledecay_settings)
{
   const auto pdg = decays_list.get_particle_id();
   const auto width = decays_list.get_total_width();
   const std::string name = HGTHDMIIMSSMBC_info::get_particle_name_from_pdg(pdg);
   if (std::isnan(width) || std::isinf(width)) {
      throw std::runtime_error("Total width of " + name + " is " + std::to_string(width));
   }

   std::ostringstream decay;

   decay << "DECAY "
         << FORMAT_TOTAL_WIDTH(pdg, width, name + " decays");

   if (!is_zero(width, 1e-100)) {
      static constexpr double NEGATIVE_BR_TOLERANCE = 1e-11;
      const double MIN_BR_TO_PRINT = flexibledecay_settings.get(FlexibleDecay_settings::min_br_to_print);
      std::vector<Decay> sorted_decays_list = sort_decays_list(decays_list);
      for (const auto& channel : sorted_decays_list) {
         auto const partial_width = channel.get_width();
         auto branching_ratio = partial_width / width;
         if (partial_width < 0 && !is_zero(branching_ratio, NEGATIVE_BR_TOLERANCE)) {
            std::stringstream ss;
            ss << std::scientific << partial_width;
            throw std::runtime_error("Error in " + channel.get_proc_string() + ": partial width is negative (" + ss.str() + " GeV).");
         }
         else if (partial_width < 0 && is_zero(branching_ratio, NEGATIVE_BR_TOLERANCE)) {
            branching_ratio = 0;
         }
         if (branching_ratio < MIN_BR_TO_PRINT) continue;
         const auto final_state = channel.get_final_state_particle_ids();
         std::string comment = "BR(" + name + " ->";
         for (auto id : final_state) {
            comment += " " + HGTHDMIIMSSMBC_info::get_particle_name_from_pdg(id);
         }
         comment += ")";

         decay << format_decay(branching_ratio, final_state, comment);
      }
   }

   slha_io.set_block(decay);
}

void HGTHDMIIMSSMBC_slha_io::set_effectivecouplings_block(const std::vector<std::tuple<int, int, int, double, std::string>>& effCouplings)
{
   slha_io.set_effectivecouplings_block(effCouplings);
}




/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in SLHA convention
 */
void HGTHDMIIMSSMBC_slha_io::set_spectrum(const HGTHDMIIMSSMBC_slha& model)
{
   const HGTHDMIIMSSMBC_physical physical(model.get_physical_slha());
   const bool write_sm_masses = model.do_calculate_sm_pole_masses();

   set_model_parameters(model);
   set_mass(physical, write_sm_masses);
   set_mixing_matrices(physical, write_sm_masses);

   if (slha_io.get_modsel().quark_flavour_violated)
      set_ckm(model.get_ckm_matrix(), model.get_scale());

   if (slha_io.get_modsel().lepton_flavour_violated)
      set_pmns(model.get_pmns_matrix(), model.get_scale());
}

/**
 * Writes extra SLHA blocks
 *
 * @param model model class
 * @param scales struct of boundary condition scales
 * @param observables struct of observables
 */
void HGTHDMIIMSSMBC_slha_io::set_extra(
   const HGTHDMIIMSSMBC_slha& model,
   const HGTHDMIIMSSMBC_scales& scales,
   const HGTHDMIIMSSMBC_observables& observables,
   const flexiblesusy::Spectrum_generator_settings& spectrum_generator_settings)
{
   const HGTHDMIIMSSMBC_physical physical(model.get_physical_slha());

   if (observables.problems.have_problem()) {
      std::ostringstream block;
      block << "Block OBSINFO\n";
      slha_format_problems_and_warnings(observables.problems,
                                        std::ostream_iterator<std::string>(block));
      slha_io.set_block(block);
   }


}

/**
 * Write SLHA object to given output.  If output == "-", then the SLHA
 * object is written to std::cout.  Otherwise, output is interpreted
 * as a file name
 *
 * @param output "-" for cout, or file name
 */
void HGTHDMIIMSSMBC_slha_io::write_to(const std::string& output) const
{
   if (output == "-")
      write_to_stream(std::cout);
   else
      write_to_file(output);
}

void HGTHDMIIMSSMBC_slha_io::write_to_file(const std::string& file_name) const
{
   slha_io.write_to_file(file_name);
}

void HGTHDMIIMSSMBC_slha_io::write_to_stream() const
{
   write_to_stream(std::cout);
}

void HGTHDMIIMSSMBC_slha_io::write_to_stream(std::ostream& ostr) const
{
   slha_io.write_to_stream(ostr);
}

/**
 * Read (DR-bar) model parameter output scale from MODSEL entry 12
 */
double HGTHDMIIMSSMBC_slha_io::get_parameter_output_scale() const
{
   return slha_io.get_modsel().parameter_output_scale;
}

/**
 * Read SLHA object from file
 *
 * @param file_name file name
 */
void HGTHDMIIMSSMBC_slha_io::read_from_file(const std::string& file_name)
{
   slha_io.read_from_file(file_name);
}

/**
 * Read SLHA object from source
 *
 * calls SLHA_io::read_from_source()
 *
 * @param source source name
 */
void HGTHDMIIMSSMBC_slha_io::read_from_source(const std::string& source)
{
   slha_io.read_from_source(source);
}

/**
 * Read SLHA object from stream
 *
 * @param istr stream name
 */
void HGTHDMIIMSSMBC_slha_io::read_from_stream(std::istream& istr)
{
   slha_io.read_from_stream(istr);
}

/**
 * Fill struct of model input parameters from SLHA object (MINPAR,
 * EXTPAR and IMEXTPAR blocks)
 *
 * @param input struct of model input parameters
 */
void HGTHDMIIMSSMBC_slha_io::fill(HGTHDMIIMSSMBC_input_parameters& input) const
{
   SLHA_io::Tuple_processor minpar_processor = [&input] (int key, double value) {
      return fill_minpar_tuple(input, key, value);
   };

   SLHA_io::Tuple_processor extpar_processor = [&input] (int key, double value) {
      return fill_extpar_tuple(input, key, value);
   };

   SLHA_io::Tuple_processor imminpar_processor = [&input] (int key, double value) {
      return fill_imminpar_tuple(input, key, value);
   };

   SLHA_io::Tuple_processor imextpar_processor = [&input] (int key, double value) {
      return fill_imextpar_tuple(input, key, value);
   };

   slha_io.read_block("MINPAR", minpar_processor);
   slha_io.read_block("EXTPAR", extpar_processor);
   slha_io.read_block("IMMINPAR", imminpar_processor);
   slha_io.read_block("IMEXTPAR", imextpar_processor);

   slha_io.read_block("AdIN", input.AdInput);
   slha_io.read_block("AeIN", input.AeInput);
   slha_io.read_block("AuIN", input.AuInput);
   slha_io.read_block("MSDIN", input.msdInput);
   slha_io.read_block("MSEIN", input.mseInput);
   slha_io.read_block("MSLIN", input.mslInput);
   slha_io.read_block("MSQIN", input.msqInput);
   slha_io.read_block("MSUIN", input.msuInput);

}

/**
 * Reads DR-bar parameters from a SLHA output file.
 *
 * @param model model class to be filled
 */
void HGTHDMIIMSSMBC_slha_io::fill_drbar_parameters(HGTHDMIIMSSMBC_mass_eigenstates& model) const
{
   model.set_g1(slha_io.read_entry("gauge", 1) * 1.2909944487358056);
   model.set_g2(slha_io.read_entry("gauge", 2));
   model.set_g3(slha_io.read_entry("gauge", 3));
   {
      Eigen::Matrix<double,3,3> Yu;
      slha_io.read_block("Yu", Yu);
      model.set_Yu(Yu);
   }
   {
      Eigen::Matrix<double,3,3> Yd;
      slha_io.read_block("Yd", Yd);
      model.set_Yd(Yd);
   }
   {
      Eigen::Matrix<double,3,3> Ye;
      slha_io.read_block("Ye", Ye);
      model.set_Ye(Ye);
   }
   model.set_Lambda1(slha_io.read_entry("HMIX", 31));
   model.set_Lambda2(slha_io.read_entry("HMIX", 32));
   model.set_Lambda3(slha_io.read_entry("HMIX", 33));
   model.set_Lambda4(slha_io.read_entry("HMIX", 34));
   model.set_Lambda5(slha_io.read_entry("HMIX", 35));
   model.set_Lambda6(slha_io.read_entry("HMIX", 36));
   model.set_Lambda7(slha_io.read_entry("HMIX", 37));
   model.set_M112(slha_io.read_entry("HMIX", 20));
   model.set_M222(slha_io.read_entry("HMIX", 21));
   model.set_M122(slha_io.read_entry("HMIX", 22));
   model.set_v1(slha_io.read_entry("HMIX", 102));
   model.set_v2(slha_io.read_entry("HMIX", 103));
   model.set_Mu(slha_io.read_entry("HMIX", 1));
   model.set_MassB(slha_io.read_entry("MSOFT", 1));
   model.set_MassWB(slha_io.read_entry("MSOFT", 2));
   model.set_MassG(slha_io.read_entry("MSOFT", 3));
   model.set_g2up(slha_io.read_entry("SplitMSSM", 1));
   model.set_g2u(slha_io.read_entry("SplitMSSM", 2));
   model.set_g1dp(slha_io.read_entry("SplitMSSM", 3));
   model.set_g1d(slha_io.read_entry("SplitMSSM", 4));


   model.set_scale(read_scale());
}

/**
 * Reads DR-bar parameters, pole masses and mixing matrices (in
 * Haber-Kane convention) from a SLHA output file.
 *
 * @param model model class to be filled
 */
void HGTHDMIIMSSMBC_slha_io::fill(HGTHDMIIMSSMBC_mass_eigenstates& model) const
{
   fill_drbar_parameters(model);

   HGTHDMIIMSSMBC_physical physical_hk;
   fill_physical(physical_hk);
   physical_hk.convert_to_hk();
   model.get_physical() = physical_hk;
}

/**
 * Fill struct of extra physical input parameters from SLHA object
 * (FlexibleSUSYInput block)
 *
 * @param input struct of physical non-SLHA input parameters
 */
void HGTHDMIIMSSMBC_slha_io::fill(Physical_input& input) const
{
   slha_io.fill(input);
}

/**
 * Fill struct of spectrum generator settings from SLHA object
 * (FlexibleSUSY block)
 *
 * @param settings struct of spectrum generator settings to be filled
 */
void HGTHDMIIMSSMBC_slha_io::fill(Spectrum_generator_settings& settings) const
{
   slha_io.fill(settings);
}

/**
 * Fill struct of spectrum generator settings from SLHA object
 * (LToLConversion block)
 *
 * @param settings struct of spectrum generator settings to be filled
 */
void HGTHDMIIMSSMBC_slha_io::fill(LToLConversion_settings& settings) const
{
   slha_io.fill(settings);
}

/**
 * Fill struct of spectrum generator settings from SLHA object
 * (FlexibleSUSY block)
 *
 * @param settings struct of spectrum generator settings to be filled
 */
void HGTHDMIIMSSMBC_slha_io::fill(FlexibleDecay_settings& settings) const
{
   slha_io.fill(settings);
}

void HGTHDMIIMSSMBC_slha_io::fill_minpar_tuple(HGTHDMIIMSSMBC_input_parameters& input,
                                                int key, double value)
{
   switch (key) {
   case 3: input.TanBeta = value; break;
   default: WARNING("Unrecognized entry in block MINPAR: " << key); break;
   }

}

void HGTHDMIIMSSMBC_slha_io::fill_extpar_tuple(HGTHDMIIMSSMBC_input_parameters& input,
                                                int key, double value)
{
   switch (key) {
   case 0: input.MSUSY = value; break;
   case 1: input.MEWSB = value; break;
   case 2: input.MuInput = value; break;
   case 3: input.M1Input = value; break;
   case 4: input.M2Input = value; break;
   case 5: input.M3Input = value; break;
   case 6: input.MAInput = value; break;
   case 100: input.LambdaLoopOrder = value; break;
   default: WARNING("Unrecognized entry in block EXTPAR: " << key); break;
   }

}

void HGTHDMIIMSSMBC_slha_io::fill_imminpar_tuple(HGTHDMIIMSSMBC_input_parameters& input,
                                                int key, double value)
{
   switch (key) {
   default: WARNING("Unrecognized entry in block IMMINPAR: " << key); break;
   }

}

void HGTHDMIIMSSMBC_slha_io::fill_imextpar_tuple(HGTHDMIIMSSMBC_input_parameters& input,
                                                int key, double value)
{
   switch (key) {
   default: WARNING("Unrecognized entry in block IMEXTPAR: " << key); break;
   }

}

/**
 * Reads pole masses and mixing matrices from a SLHA output file to be filled.
 */
void HGTHDMIIMSSMBC_slha_io::fill_physical(HGTHDMIIMSSMBC_physical& physical) const
{
   {
      DEFINE_PHYSICAL_PARAMETER(ZH);
      slha_io.read_block("SCALARMIX", ZH);
      LOCALPHYSICAL(ZH) = ZH;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZA);
      slha_io.read_block("PSEUDOSCALARMIX", ZA);
      LOCALPHYSICAL(ZA) = ZA;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZP);
      slha_io.read_block("CHARGEMIX", ZP);
      LOCALPHYSICAL(ZP) = ZP;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(Vu);
      slha_io.read_block("UULMIX", Vu);
      LOCALPHYSICAL(Vu) = Vu;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(Vd);
      slha_io.read_block("UDLMIX", Vd);
      LOCALPHYSICAL(Vd) = Vd;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(Uu);
      slha_io.read_block("UURMIX", Uu);
      LOCALPHYSICAL(Uu) = Uu;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(Ud);
      slha_io.read_block("UDRMIX", Ud);
      LOCALPHYSICAL(Ud) = Ud;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(Ve);
      slha_io.read_block("UELMIX", Ve);
      LOCALPHYSICAL(Ve) = Ve;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(Ue);
      slha_io.read_block("UERMIX", Ue);
      LOCALPHYSICAL(Ue) = Ue;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZN);
      slha_io.read_block("NMIX", ZN);
      LOCALPHYSICAL(ZN) = ZN;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(UP);
      slha_io.read_block("VMIX", UP);
      LOCALPHYSICAL(UP) = UP;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(UM);
      slha_io.read_block("UMIX", UM);
      LOCALPHYSICAL(UM) = UM;
   }

   LOCALPHYSICAL(MVG) = slha_io.read_entry("MASS", 21);
   LOCALPHYSICAL(MFv)(0) = slha_io.read_entry("MASS", 12);
   LOCALPHYSICAL(MFv)(1) = slha_io.read_entry("MASS", 14);
   LOCALPHYSICAL(MFv)(2) = slha_io.read_entry("MASS", 16);
   LOCALPHYSICAL(MGlu) = slha_io.read_entry("MASS", 1000021);
   LOCALPHYSICAL(MVP) = slha_io.read_entry("MASS", 22);
   LOCALPHYSICAL(MVZ) = slha_io.read_entry("MASS", 23);
   LOCALPHYSICAL(Mhh)(0) = slha_io.read_entry("MASS", 25);
   LOCALPHYSICAL(Mhh)(1) = slha_io.read_entry("MASS", 35);
   LOCALPHYSICAL(MAh)(1) = slha_io.read_entry("MASS", 36);
   LOCALPHYSICAL(MHm)(1) = slha_io.read_entry("MASS", 37);
   LOCALPHYSICAL(MFd)(0) = slha_io.read_entry("MASS", 1);
   LOCALPHYSICAL(MFd)(1) = slha_io.read_entry("MASS", 3);
   LOCALPHYSICAL(MFd)(2) = slha_io.read_entry("MASS", 5);
   LOCALPHYSICAL(MFu)(0) = slha_io.read_entry("MASS", 2);
   LOCALPHYSICAL(MFu)(1) = slha_io.read_entry("MASS", 4);
   LOCALPHYSICAL(MFu)(2) = slha_io.read_entry("MASS", 6);
   LOCALPHYSICAL(MFe)(0) = slha_io.read_entry("MASS", 11);
   LOCALPHYSICAL(MFe)(1) = slha_io.read_entry("MASS", 13);
   LOCALPHYSICAL(MFe)(2) = slha_io.read_entry("MASS", 15);
   LOCALPHYSICAL(MChi)(0) = slha_io.read_entry("MASS", 1000022);
   LOCALPHYSICAL(MChi)(1) = slha_io.read_entry("MASS", 1000023);
   LOCALPHYSICAL(MChi)(2) = slha_io.read_entry("MASS", 1000025);
   LOCALPHYSICAL(MChi)(3) = slha_io.read_entry("MASS", 1000035);
   LOCALPHYSICAL(MCha)(0) = slha_io.read_entry("MASS", 1000024);
   LOCALPHYSICAL(MCha)(1) = slha_io.read_entry("MASS", 1000037);
   LOCALPHYSICAL(MVWm) = slha_io.read_entry("MASS", 24);

}

/**
 * Reads the renormalization scales from all DR-bar parameter blocks.
 * If blocks with different scales are found the last scale is
 * returned and a warning is printed.
 *
 * @return common renormalization scale
 */
double HGTHDMIIMSSMBC_slha_io::read_scale() const
{
   static const std::array<std::string, 7> drbar_blocks =
      { "gauge", "Yu", "Yd", "Ye", "HMIX", "MSOFT", "SplitMSSM" };

   double scale = 0.;

   for (const auto& block: drbar_blocks) {
      const double block_scale = slha_io.read_scale(block);
      if (!is_zero(block_scale)) {
         if (!is_zero(scale) && !is_equal(scale, block_scale))
            WARNING("DR-bar parameters defined at different scales");
         scale = block_scale;
      }
   }

   return scale;
}

} // namespace flexiblesusy
