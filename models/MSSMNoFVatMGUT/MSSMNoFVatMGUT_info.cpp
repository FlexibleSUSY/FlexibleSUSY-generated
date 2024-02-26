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


#include "MSSMNoFVatMGUT_info.hpp"

#include "error.hpp"

#include <iostream>
#include <vector>

namespace flexiblesusy {

namespace MSSMNoFVatMGUT_info {
   const double normalization_g1 = 0.7745966692414834;
   const double normalization_g2 = 1;
   const double normalization_g3 = 1;

   const std::array<int, NUMBER_OF_PARTICLES> particle_multiplicities = {1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2
      , 4, 2, 1, 1, 1};

   const std::array<std::string, NUMBER_OF_PARTICLES> particle_names = {"VG",
      "Glu", "Fd", "Fs", "Fb", "Fu", "Fc", "Ft", "Fve", "Fvm", "Fvt", "Fe", "Fm",
      "Ftau", "SveL", "SvmL", "SvtL", "Sd", "Su", "Se", "Sm", "Stau", "Ss", "Sc",
      "Sb", "St", "hh", "Ah", "Hpm", "Chi", "Cha", "VWm", "VP", "VZ"};

   const std::array<std::string, NUMBER_OF_PARTICLES> particle_latex_names = {
      "g", "\\tilde{g}", "d", "s", "b", "u", "c", "t", "\\nu_e", "\\nu_{\\mu}",
      "\\nu_{\\tau}", "e", "m", "\\tau", "\\tilde{\\nu}_e", "\\tilde{\\nu}_{\\mu}"
      , "\\tilde{\\nu}_{\\tau}", "\\tilde{d}", "\\tilde{u}", "\\tilde{e}",
      "\\tilde{\\mu}", "\\tilde{\\tau}", "\\tilde{s}", "\\tilde{c}", "\\tilde{b}",
      "\\tilde{t}", "h", "A^0", "H^-", "\\tilde{\\chi}^0", "\\tilde{\\chi}^-",
      "W^-", "\\gamma", "Z"};

   const std::array<std::string, NUMBER_OF_PARAMETERS> parameter_names = {
      "Yd(0,0)", "Yd(0,1)", "Yd(0,2)", "Yd(1,0)", "Yd(1,1)", "Yd(1,2)", "Yd(2,0)",
      "Yd(2,1)", "Yd(2,2)", "Ye(0,0)", "Ye(0,1)", "Ye(0,2)", "Ye(1,0)", "Ye(1,1)",
      "Ye(1,2)", "Ye(2,0)", "Ye(2,1)", "Ye(2,2)", "Yu(0,0)", "Yu(0,1)", "Yu(0,2)",
      "Yu(1,0)", "Yu(1,1)", "Yu(1,2)", "Yu(2,0)", "Yu(2,1)", "Yu(2,2)", "Mu", "g1"
      , "g2", "g3", "vd", "vu", "TYd(0,0)", "TYd(0,1)", "TYd(0,2)", "TYd(1,0)",
      "TYd(1,1)", "TYd(1,2)", "TYd(2,0)", "TYd(2,1)", "TYd(2,2)", "TYe(0,0)",
      "TYe(0,1)", "TYe(0,2)", "TYe(1,0)", "TYe(1,1)", "TYe(1,2)", "TYe(2,0)",
      "TYe(2,1)", "TYe(2,2)", "TYu(0,0)", "TYu(0,1)", "TYu(0,2)", "TYu(1,0)",
      "TYu(1,1)", "TYu(1,2)", "TYu(2,0)", "TYu(2,1)", "TYu(2,2)", "BMu",
      "mq2(0,0)", "mq2(0,1)", "mq2(0,2)", "mq2(1,0)", "mq2(1,1)", "mq2(1,2)",
      "mq2(2,0)", "mq2(2,1)", "mq2(2,2)", "ml2(0,0)", "ml2(0,1)", "ml2(0,2)",
      "ml2(1,0)", "ml2(1,1)", "ml2(1,2)", "ml2(2,0)", "ml2(2,1)", "ml2(2,2)",
      "mHd2", "mHu2", "md2(0,0)", "md2(0,1)", "md2(0,2)", "md2(1,0)", "md2(1,1)",
      "md2(1,2)", "md2(2,0)", "md2(2,1)", "md2(2,2)", "mu2(0,0)", "mu2(0,1)",
      "mu2(0,2)", "mu2(1,0)", "mu2(1,1)", "mu2(1,2)", "mu2(2,0)", "mu2(2,1)",
      "mu2(2,2)", "me2(0,0)", "me2(0,1)", "me2(0,2)", "me2(1,0)", "me2(1,1)",
      "me2(1,2)", "me2(2,0)", "me2(2,1)", "me2(2,2)", "MassB", "MassWB", "MassG"};

   const std::array<std::string, NUMBER_OF_MIXINGS> particle_mixing_names = {
      "ZD(0,0)", "ZD(0,1)", "ZD(1,0)", "ZD(1,1)", "ZU(0,0)", "ZU(0,1)", "ZU(1,0)",
      "ZU(1,1)", "ZE(0,0)", "ZE(0,1)", "ZE(1,0)", "ZE(1,1)", "ZM(0,0)", "ZM(0,1)",
      "ZM(1,0)", "ZM(1,1)", "ZTau(0,0)", "ZTau(0,1)", "ZTau(1,0)", "ZTau(1,1)",
      "ZS(0,0)", "ZS(0,1)", "ZS(1,0)", "ZS(1,1)", "ZC(0,0)", "ZC(0,1)", "ZC(1,0)",
      "ZC(1,1)", "ZB(0,0)", "ZB(0,1)", "ZB(1,0)", "ZB(1,1)", "ZT(0,0)", "ZT(0,1)",
      "ZT(1,0)", "ZT(1,1)", "ZH(0,0)", "ZH(0,1)", "ZH(1,0)", "ZH(1,1)", "ZA(0,0)",
      "ZA(0,1)", "ZA(1,0)", "ZA(1,1)", "ZP(0,0)", "ZP(0,1)", "ZP(1,0)", "ZP(1,1)",
      "Re(ZN(0,0))", "Im(ZN(0,0))", "Re(ZN(0,1))", "Im(ZN(0,1))", "Re(ZN(0,2))",
      "Im(ZN(0,2))", "Re(ZN(0,3))", "Im(ZN(0,3))", "Re(ZN(1,0))", "Im(ZN(1,0))",
      "Re(ZN(1,1))", "Im(ZN(1,1))", "Re(ZN(1,2))", "Im(ZN(1,2))", "Re(ZN(1,3))",
      "Im(ZN(1,3))", "Re(ZN(2,0))", "Im(ZN(2,0))", "Re(ZN(2,1))", "Im(ZN(2,1))",
      "Re(ZN(2,2))", "Im(ZN(2,2))", "Re(ZN(2,3))", "Im(ZN(2,3))", "Re(ZN(3,0))",
      "Im(ZN(3,0))", "Re(ZN(3,1))", "Im(ZN(3,1))", "Re(ZN(3,2))", "Im(ZN(3,2))",
      "Re(ZN(3,3))", "Im(ZN(3,3))", "Re(UM(0,0))", "Im(UM(0,0))", "Re(UM(0,1))",
      "Im(UM(0,1))", "Re(UM(1,0))", "Im(UM(1,0))", "Re(UM(1,1))", "Im(UM(1,1))",
      "Re(UP(0,0))", "Im(UP(0,0))", "Re(UP(0,1))", "Im(UP(0,1))", "Re(UP(1,0))",
      "Im(UP(1,0))", "Re(UP(1,1))", "Im(UP(1,1))", "ZZ(0,0)", "ZZ(0,1)", "ZZ(1,0)"
      , "ZZ(1,1)"};

   const std::array<std::string, NUMBER_OF_INPUT_PARAMETERS> input_parameter_names
       = {"TanBeta", "Sign(Mu)", "M1", "M2", "M3", "AtIN", "AbIN", "AtauIN",
      "AcIN", "AsIN", "AmuonIN", "AuIN", "AdIN", "AeIN", "mHd2IN", "mHu2IN",
      "ml11IN", "ml22IN", "ml33IN", "me11IN", "me22IN", "me33IN", "mq11IN",
      "mq22IN", "mq33IN", "mu11IN", "mu22IN", "mu33IN", "md11IN", "md22IN",
      "md33IN"};

   const std::array<std::string, NUMBER_OF_EXTRA_PARAMETERS> extra_parameter_names
       = {};

   const std::string model_name = "MSSMNoFVatMGUT";

int get_pdg_code_for_particle(Particles p)
{
   if (particle_multiplicities[p] > 1) {
      throw OutOfBoundsError(particle_names[p] + " must have a generation index");
   }

   int pdg = 0;
   switch (p) {

   case VG: pdg = 21; break;
   case Glu: pdg = 1000021; break;
   case Fd: pdg = 1; break;
   case Fs: pdg = 3; break;
   case Fb: pdg = 5; break;
   case Fu: pdg = 2; break;
   case Fc: pdg = 4; break;
   case Ft: pdg = 6; break;
   case Fve: pdg = 12; break;
   case Fvm: pdg = 14; break;
   case Fvt: pdg = 16; break;
   case Fe: pdg = 11; break;
   case Fm: pdg = 13; break;
   case Ftau: pdg = 15; break;
   case SveL: pdg = 1000012; break;
   case SvmL: pdg = 1000014; break;
   case SvtL: pdg = 1000016; break;
   case VWm: pdg = -24; break;
   case VP: pdg = 22; break;
   case VZ: pdg = 23; break;

   default: throw OutOfBoundsError("invalid particle " + std::to_string(p));
   }

   return pdg;
}

int get_pdg_code_for_particle(Particles p, int index)
{
   if (particle_multiplicities[p] == 1) {
      throw OutOfBoundsError(particle_names[p] + " does not carry an index");
   }

   std::vector<int> pdg_codes;
   switch (p) {

   case Sd: pdg_codes = {1000001, 2000001}; break;
   case Su: pdg_codes = {1000002, 2000002}; break;
   case Se: pdg_codes = {1000011, 2000011}; break;
   case Sm: pdg_codes = {1000013, 2000013}; break;
   case Stau: pdg_codes = {1000015, 2000015}; break;
   case Ss: pdg_codes = {1000003, 2000003}; break;
   case Sc: pdg_codes = {1000004, 2000004}; break;
   case Sb: pdg_codes = {1000005, 2000005}; break;
   case St: pdg_codes = {1000006, 2000006}; break;
   case hh: pdg_codes = {25, 35}; break;
   case Ah: pdg_codes = {0, 36}; break;
   case Hpm: pdg_codes = {0, -37}; break;
   case Chi: pdg_codes = {1000022, 1000023, 1000025, 1000035}; break;
   case Cha: pdg_codes = {-1000024, -1000037}; break;

   default: throw OutOfBoundsError("invalid particle " + std::to_string(p));
   }

   if (index < 0 || std::abs(index) >= pdg_codes.size()) {
      throw OutOfBoundsError("index " + std::to_string(index) + " out of bounds");
   }

   return pdg_codes[index];
}

std::pair<std::string, std::optional<unsigned int>> get_multiplet_and_index_from_pdg(int pdg)
{
   std::pair<std::string, std::optional<unsigned int>> name;

   switch (pdg) {

   case 21: name = {"VG", {}}; break;
   case 1000021: name = {"Glu", {}}; break;
   case 1: name = {"Fd", {}}; break;
   case 3: name = {"Fs", {}}; break;
   case 5: name = {"Fb", {}}; break;
   case 2: name = {"Fu", {}}; break;
   case 4: name = {"Fc", {}}; break;
   case 6: name = {"Ft", {}}; break;
   case 12: name = {"Fve", {}}; break;
   case 14: name = {"Fvm", {}}; break;
   case 16: name = {"Fvt", {}}; break;
   case 11: name = {"Fe", {}}; break;
   case 13: name = {"Fm", {}}; break;
   case 15: name = {"Ftau", {}}; break;
   case 1000012: name = {"SveL", {}}; break;
   case 1000014: name = {"SvmL", {}}; break;
   case 1000016: name = {"SvtL", {}}; break;
   case 1000001: name = {"Sd", 1}; break;
   case 2000001: name = {"Sd", 2}; break;
   case 1000002: name = {"Su", 1}; break;
   case 2000002: name = {"Su", 2}; break;
   case 1000011: name = {"Se", 1}; break;
   case 2000011: name = {"Se", 2}; break;
   case 1000013: name = {"Sm", 1}; break;
   case 2000013: name = {"Sm", 2}; break;
   case 1000015: name = {"Stau", 1}; break;
   case 2000015: name = {"Stau", 2}; break;
   case 1000003: name = {"Ss", 1}; break;
   case 2000003: name = {"Ss", 2}; break;
   case 1000004: name = {"Sc", 1}; break;
   case 2000004: name = {"Sc", 2}; break;
   case 1000005: name = {"Sb", 1}; break;
   case 2000005: name = {"Sb", 2}; break;
   case 1000006: name = {"St", 1}; break;
   case 2000006: name = {"St", 2}; break;
   case 25: name = {"hh", 1}; break;
   case 35: name = {"hh", 2}; break;
   case 36: name = {"Ah", 2}; break;
   case -37: name = {"Hpm", 2}; break;
   case 1000022: name = {"Chi", 1}; break;
   case 1000023: name = {"Chi", 2}; break;
   case 1000025: name = {"Chi", 3}; break;
   case 1000035: name = {"Chi", 4}; break;
   case -1000024: name = {"Cha", 1}; break;
   case -1000037: name = {"Cha", 2}; break;
   case -24: name = {"VWm", {}}; break;
   case 22: name = {"VP", {}}; break;
   case 23: name = {"VZ", {}}; break;
   case -1: name = {"barFd", {}}; break;
   case -3: name = {"barFs", {}}; break;
   case -5: name = {"barFb", {}}; break;
   case -2: name = {"barFu", {}}; break;
   case -4: name = {"barFc", {}}; break;
   case -6: name = {"barFt", {}}; break;
   case -12: name = {"barFve", {}}; break;
   case -14: name = {"barFvm", {}}; break;
   case -16: name = {"barFvt", {}}; break;
   case -11: name = {"barFe", {}}; break;
   case -13: name = {"barFm", {}}; break;
   case -15: name = {"barFtau", {}}; break;
   case -1000012: name = {"conjSveL", {}}; break;
   case -1000014: name = {"conjSvmL", {}}; break;
   case -1000016: name = {"conjSvtL", {}}; break;
   case -1000001: name = {"conjSd", 1}; break;
   case -2000001: name = {"conjSd", 2}; break;
   case -1000002: name = {"conjSu", 1}; break;
   case -2000002: name = {"conjSu", 2}; break;
   case -1000011: name = {"conjSe", 1}; break;
   case -2000011: name = {"conjSe", 2}; break;
   case -1000013: name = {"conjSm", 1}; break;
   case -2000013: name = {"conjSm", 2}; break;
   case -1000015: name = {"conjStau", 1}; break;
   case -2000015: name = {"conjStau", 2}; break;
   case -1000003: name = {"conjSs", 1}; break;
   case -2000003: name = {"conjSs", 2}; break;
   case -1000004: name = {"conjSc", 1}; break;
   case -2000004: name = {"conjSc", 2}; break;
   case -1000005: name = {"conjSb", 1}; break;
   case -2000005: name = {"conjSb", 2}; break;
   case -1000006: name = {"conjSt", 1}; break;
   case -2000006: name = {"conjSt", 2}; break;
   case 37: name = {"conjHpm", 2}; break;
   case 1000024: name = {"barCha", 1}; break;
   case 1000037: name = {"barCha", 2}; break;
   case 24: name = {"conjVWm", {}}; break;

   default: name = {"", {}};
   }

   return name;
}

std::string get_particle_name_from_pdg(int pdg)
{
   std::pair<std::string, std::optional<unsigned int>> const pair = get_multiplet_and_index_from_pdg(pdg);
   return pair.first + (pair.second.has_value() ? "(" + std::to_string(pair.second.value()) + ")" : "");
}

void print(std::ostream& ostr)
{
   ostr
      << "Model information\n"
      << "=================\n"
      << "Model name:                  " << model_name << '\n'
      << "Is a low-energy model:       "
      << (is_low_energy_model ? "yes" : "no") << '\n'
      << "Is a supersymmetric model:   "
      << (is_supersymmetric_model ? "yes" : "no") << '\n'
      << "Is a FlexibleEFTHiggs model: "
      << (is_FlexibleEFTHiggs ? "yes" : "no") << '\n'
      << "Number of multiplets:        " << NUMBER_OF_PARTICLES << '\n'
      << "Number of parameters:        " << NUMBER_OF_PARAMETERS << '\n'
      ;

   ostr << "\n"
      "Multiplets:                  ";
   for (int i = 0; i < NUMBER_OF_PARTICLES; i++) {
      ostr << particle_names[i]
           << '[' << particle_multiplicities[i] << ']';
      if (i + 1 < NUMBER_OF_PARTICLES)
         ostr << ", ";
   }

   ostr << "\n\n"
      "Parameters:                  ";
   for (int i = 0; i < NUMBER_OF_PARAMETERS; i++) {
      ostr << parameter_names[i];
      if (i + 1 < NUMBER_OF_PARAMETERS)
         ostr << ", ";
   }

   ostr << "\n\n"
      "Input parameters:            ";
   for (int i = 0; i < NUMBER_OF_INPUT_PARAMETERS; i++) {
      ostr << input_parameter_names[i];
      if (i + 1 < NUMBER_OF_INPUT_PARAMETERS)
         ostr << ", ";
   }

   ostr << '\n';
}

} // namespace MSSMNoFVatMGUT_info

} // namespace flexiblesusy

