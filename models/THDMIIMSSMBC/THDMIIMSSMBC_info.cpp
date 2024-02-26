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


#include "THDMIIMSSMBC_info.hpp"

#include "error.hpp"

#include <iostream>
#include <vector>

namespace flexiblesusy {

namespace THDMIIMSSMBC_info {
   const double normalization_g1 = 0.7745966692414834;
   const double normalization_g2 = 1;
   const double normalization_g3 = 1;

   const std::array<int, NUMBER_OF_PARTICLES> particle_multiplicities = {1, 3, 2,
      2, 2, 3, 3, 3, 1, 1, 1};

   const std::array<std::string, NUMBER_OF_PARTICLES> particle_names = {"VG", "Fv"
      , "hh", "Ah", "Hm", "Fd", "Fu", "Fe", "VWm", "VP", "VZ"};

   const std::array<std::string, NUMBER_OF_PARTICLES> particle_latex_names = {
      "g", "\\nu", "h", "A^0", "H^-", "d", "u", "e", "W^-", "\\gamma", "Z"};

   const std::array<std::string, NUMBER_OF_PARAMETERS> parameter_names = {"g1",
      "g2", "g3", "Lambda6", "Lambda5", "Lambda7", "Lambda1", "Lambda4", "Lambda3"
      , "Lambda2", "Yu(0,0)", "Yu(0,1)", "Yu(0,2)", "Yu(1,0)", "Yu(1,1)",
      "Yu(1,2)", "Yu(2,0)", "Yu(2,1)", "Yu(2,2)", "Yd(0,0)", "Yd(0,1)", "Yd(0,2)",
      "Yd(1,0)", "Yd(1,1)", "Yd(1,2)", "Yd(2,0)", "Yd(2,1)", "Yd(2,2)", "Ye(0,0)",
      "Ye(0,1)", "Ye(0,2)", "Ye(1,0)", "Ye(1,1)", "Ye(1,2)", "Ye(2,0)", "Ye(2,1)",
      "Ye(2,2)", "M122", "M112", "M222", "v1", "v2"};

   const std::array<std::string, NUMBER_OF_MIXINGS> particle_mixing_names = {
      "ZH(0,0)", "ZH(0,1)", "ZH(1,0)", "ZH(1,1)", "ZA(0,0)", "ZA(0,1)", "ZA(1,0)",
      "ZA(1,1)", "ZP(0,0)", "ZP(0,1)", "ZP(1,0)", "ZP(1,1)", "Re(Vd(0,0))",
      "Im(Vd(0,0))", "Re(Vd(0,1))", "Im(Vd(0,1))", "Re(Vd(0,2))", "Im(Vd(0,2))",
      "Re(Vd(1,0))", "Im(Vd(1,0))", "Re(Vd(1,1))", "Im(Vd(1,1))", "Re(Vd(1,2))",
      "Im(Vd(1,2))", "Re(Vd(2,0))", "Im(Vd(2,0))", "Re(Vd(2,1))", "Im(Vd(2,1))",
      "Re(Vd(2,2))", "Im(Vd(2,2))", "Re(Ud(0,0))", "Im(Ud(0,0))", "Re(Ud(0,1))",
      "Im(Ud(0,1))", "Re(Ud(0,2))", "Im(Ud(0,2))", "Re(Ud(1,0))", "Im(Ud(1,0))",
      "Re(Ud(1,1))", "Im(Ud(1,1))", "Re(Ud(1,2))", "Im(Ud(1,2))", "Re(Ud(2,0))",
      "Im(Ud(2,0))", "Re(Ud(2,1))", "Im(Ud(2,1))", "Re(Ud(2,2))", "Im(Ud(2,2))",
      "Re(Vu(0,0))", "Im(Vu(0,0))", "Re(Vu(0,1))", "Im(Vu(0,1))", "Re(Vu(0,2))",
      "Im(Vu(0,2))", "Re(Vu(1,0))", "Im(Vu(1,0))", "Re(Vu(1,1))", "Im(Vu(1,1))",
      "Re(Vu(1,2))", "Im(Vu(1,2))", "Re(Vu(2,0))", "Im(Vu(2,0))", "Re(Vu(2,1))",
      "Im(Vu(2,1))", "Re(Vu(2,2))", "Im(Vu(2,2))", "Re(Uu(0,0))", "Im(Uu(0,0))",
      "Re(Uu(0,1))", "Im(Uu(0,1))", "Re(Uu(0,2))", "Im(Uu(0,2))", "Re(Uu(1,0))",
      "Im(Uu(1,0))", "Re(Uu(1,1))", "Im(Uu(1,1))", "Re(Uu(1,2))", "Im(Uu(1,2))",
      "Re(Uu(2,0))", "Im(Uu(2,0))", "Re(Uu(2,1))", "Im(Uu(2,1))", "Re(Uu(2,2))",
      "Im(Uu(2,2))", "Re(Ve(0,0))", "Im(Ve(0,0))", "Re(Ve(0,1))", "Im(Ve(0,1))",
      "Re(Ve(0,2))", "Im(Ve(0,2))", "Re(Ve(1,0))", "Im(Ve(1,0))", "Re(Ve(1,1))",
      "Im(Ve(1,1))", "Re(Ve(1,2))", "Im(Ve(1,2))", "Re(Ve(2,0))", "Im(Ve(2,0))",
      "Re(Ve(2,1))", "Im(Ve(2,1))", "Re(Ve(2,2))", "Im(Ve(2,2))", "Re(Ue(0,0))",
      "Im(Ue(0,0))", "Re(Ue(0,1))", "Im(Ue(0,1))", "Re(Ue(0,2))", "Im(Ue(0,2))",
      "Re(Ue(1,0))", "Im(Ue(1,0))", "Re(Ue(1,1))", "Im(Ue(1,1))", "Re(Ue(1,2))",
      "Im(Ue(1,2))", "Re(Ue(2,0))", "Im(Ue(2,0))", "Re(Ue(2,1))", "Im(Ue(2,1))",
      "Re(Ue(2,2))", "Im(Ue(2,2))", "ZZ(0,0)", "ZZ(0,1)", "ZZ(1,0)", "ZZ(1,1)"};

   const std::array<std::string, NUMBER_OF_INPUT_PARAMETERS> input_parameter_names
       = {"TanBeta", "MSUSY", "MEWSB", "MuInput", "M1Input", "M2Input", "MAInput",
      "LambdaLoopOrder", "AeInput(0,0)", "AeInput(0,1)", "AeInput(0,2)",
      "AeInput(1,0)", "AeInput(1,1)", "AeInput(1,2)", "AeInput(2,0)",
      "AeInput(2,1)", "AeInput(2,2)", "AdInput(0,0)", "AdInput(0,1)",
      "AdInput(0,2)", "AdInput(1,0)", "AdInput(1,1)", "AdInput(1,2)",
      "AdInput(2,0)", "AdInput(2,1)", "AdInput(2,2)", "AuInput(0,0)",
      "AuInput(0,1)", "AuInput(0,2)", "AuInput(1,0)", "AuInput(1,1)",
      "AuInput(1,2)", "AuInput(2,0)", "AuInput(2,1)", "AuInput(2,2)",
      "mslInput(0)", "mslInput(1)", "mslInput(2)", "mseInput(0)", "mseInput(1)",
      "mseInput(2)", "msqInput(0)", "msqInput(1)", "msqInput(2)", "msdInput(0)",
      "msdInput(1)", "msdInput(2)", "msuInput(0)", "msuInput(1)", "msuInput(2)"};

   const std::array<std::string, NUMBER_OF_EXTRA_PARAMETERS> extra_parameter_names
       = {};

   const std::string model_name = "THDMIIMSSMBC";

int get_pdg_code_for_particle(Particles p)
{
   if (particle_multiplicities[p] > 1) {
      throw OutOfBoundsError(particle_names[p] + " must have a generation index");
   }

   int pdg = 0;
   switch (p) {

   case VG: pdg = 21; break;
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

   case Fv: pdg_codes = {12, 14, 16}; break;
   case hh: pdg_codes = {25, 35}; break;
   case Ah: pdg_codes = {0, 36}; break;
   case Hm: pdg_codes = {0, -37}; break;
   case Fd: pdg_codes = {1, 3, 5}; break;
   case Fu: pdg_codes = {2, 4, 6}; break;
   case Fe: pdg_codes = {11, 13, 15}; break;

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
   case 12: name = {"Fv", 1}; break;
   case 14: name = {"Fv", 2}; break;
   case 16: name = {"Fv", 3}; break;
   case 25: name = {"hh", 1}; break;
   case 35: name = {"hh", 2}; break;
   case 36: name = {"Ah", 2}; break;
   case -37: name = {"Hm", 2}; break;
   case 1: name = {"Fd", 1}; break;
   case 3: name = {"Fd", 2}; break;
   case 5: name = {"Fd", 3}; break;
   case 2: name = {"Fu", 1}; break;
   case 4: name = {"Fu", 2}; break;
   case 6: name = {"Fu", 3}; break;
   case 11: name = {"Fe", 1}; break;
   case 13: name = {"Fe", 2}; break;
   case 15: name = {"Fe", 3}; break;
   case -24: name = {"VWm", {}}; break;
   case 22: name = {"VP", {}}; break;
   case 23: name = {"VZ", {}}; break;
   case -12: name = {"barFv", 1}; break;
   case -14: name = {"barFv", 2}; break;
   case -16: name = {"barFv", 3}; break;
   case 37: name = {"conjHm", 2}; break;
   case -1: name = {"barFd", 1}; break;
   case -3: name = {"barFd", 2}; break;
   case -5: name = {"barFd", 3}; break;
   case -2: name = {"barFu", 1}; break;
   case -4: name = {"barFu", 2}; break;
   case -6: name = {"barFu", 3}; break;
   case -11: name = {"barFe", 1}; break;
   case -13: name = {"barFe", 2}; break;
   case -15: name = {"barFe", 3}; break;
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

} // namespace THDMIIMSSMBC_info

} // namespace flexiblesusy

