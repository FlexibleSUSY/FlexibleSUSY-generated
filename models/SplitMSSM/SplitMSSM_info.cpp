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


#include "SplitMSSM_info.hpp"

#include "error.hpp"

#include <iostream>
#include <vector>

namespace flexiblesusy {

namespace SplitMSSM_info {
   const double normalization_g1 = 0.7745966692414834;
   const double normalization_g2 = 1;
   const double normalization_g3 = 1;

   const std::array<int, NUMBER_OF_PARTICLES> particle_multiplicities = {1, 1, 3,
      1, 1, 1, 3, 3, 3, 4, 2, 1, 1, 1};

   const std::array<std::string, NUMBER_OF_PARTICLES> particle_names = {"VG", "Hp"
      , "Fv", "Glu", "Ah", "hh", "Fd", "Fu", "Fe", "Chi", "Cha", "VWp", "VP", "VZ"
      };

   const std::array<std::string, NUMBER_OF_PARTICLES> particle_latex_names = {
      "g", "H^+", "\\nu", "\\tilde{g}", "A^0", "h", "d", "u", "e",
      "\\tilde{\\chi}^0", "\\tilde{\\chi}^-", "W^+", "\\gamma", "Z"};

   const std::array<std::string, NUMBER_OF_PARAMETERS> parameter_names = {"g1",
      "g2", "g3", "Lambdax", "Yu(0,0)", "Yu(0,1)", "Yu(0,2)", "Yu(1,0)", "Yu(1,1)"
      , "Yu(1,2)", "Yu(2,0)", "Yu(2,1)", "Yu(2,2)", "Yd(0,0)", "Yd(0,1)",
      "Yd(0,2)", "Yd(1,0)", "Yd(1,1)", "Yd(1,2)", "Yd(2,0)", "Yd(2,1)", "Yd(2,2)",
      "Ye(0,0)", "Ye(0,1)", "Ye(0,2)", "Ye(1,0)", "Ye(1,1)", "Ye(1,2)", "Ye(2,0)",
      "Ye(2,1)", "Ye(2,2)", "gYd", "g2d", "gYu", "g2u", "MassB", "MassG", "MassWB"
      , "Mu", "mu2", "v"};

   const std::array<std::string, NUMBER_OF_MIXINGS> particle_mixing_names = {
      "Re(Vd(0,0))", "Im(Vd(0,0))", "Re(Vd(0,1))", "Im(Vd(0,1))", "Re(Vd(0,2))",
      "Im(Vd(0,2))", "Re(Vd(1,0))", "Im(Vd(1,0))", "Re(Vd(1,1))", "Im(Vd(1,1))",
      "Re(Vd(1,2))", "Im(Vd(1,2))", "Re(Vd(2,0))", "Im(Vd(2,0))", "Re(Vd(2,1))",
      "Im(Vd(2,1))", "Re(Vd(2,2))", "Im(Vd(2,2))", "Re(Ud(0,0))", "Im(Ud(0,0))",
      "Re(Ud(0,1))", "Im(Ud(0,1))", "Re(Ud(0,2))", "Im(Ud(0,2))", "Re(Ud(1,0))",
      "Im(Ud(1,0))", "Re(Ud(1,1))", "Im(Ud(1,1))", "Re(Ud(1,2))", "Im(Ud(1,2))",
      "Re(Ud(2,0))", "Im(Ud(2,0))", "Re(Ud(2,1))", "Im(Ud(2,1))", "Re(Ud(2,2))",
      "Im(Ud(2,2))", "Re(Vu(0,0))", "Im(Vu(0,0))", "Re(Vu(0,1))", "Im(Vu(0,1))",
      "Re(Vu(0,2))", "Im(Vu(0,2))", "Re(Vu(1,0))", "Im(Vu(1,0))", "Re(Vu(1,1))",
      "Im(Vu(1,1))", "Re(Vu(1,2))", "Im(Vu(1,2))", "Re(Vu(2,0))", "Im(Vu(2,0))",
      "Re(Vu(2,1))", "Im(Vu(2,1))", "Re(Vu(2,2))", "Im(Vu(2,2))", "Re(Uu(0,0))",
      "Im(Uu(0,0))", "Re(Uu(0,1))", "Im(Uu(0,1))", "Re(Uu(0,2))", "Im(Uu(0,2))",
      "Re(Uu(1,0))", "Im(Uu(1,0))", "Re(Uu(1,1))", "Im(Uu(1,1))", "Re(Uu(1,2))",
      "Im(Uu(1,2))", "Re(Uu(2,0))", "Im(Uu(2,0))", "Re(Uu(2,1))", "Im(Uu(2,1))",
      "Re(Uu(2,2))", "Im(Uu(2,2))", "Re(Ve(0,0))", "Im(Ve(0,0))", "Re(Ve(0,1))",
      "Im(Ve(0,1))", "Re(Ve(0,2))", "Im(Ve(0,2))", "Re(Ve(1,0))", "Im(Ve(1,0))",
      "Re(Ve(1,1))", "Im(Ve(1,1))", "Re(Ve(1,2))", "Im(Ve(1,2))", "Re(Ve(2,0))",
      "Im(Ve(2,0))", "Re(Ve(2,1))", "Im(Ve(2,1))", "Re(Ve(2,2))", "Im(Ve(2,2))",
      "Re(Ue(0,0))", "Im(Ue(0,0))", "Re(Ue(0,1))", "Im(Ue(0,1))", "Re(Ue(0,2))",
      "Im(Ue(0,2))", "Re(Ue(1,0))", "Im(Ue(1,0))", "Re(Ue(1,1))", "Im(Ue(1,1))",
      "Re(Ue(1,2))", "Im(Ue(1,2))", "Re(Ue(2,0))", "Im(Ue(2,0))", "Re(Ue(2,1))",
      "Im(Ue(2,1))", "Re(Ue(2,2))", "Im(Ue(2,2))", "Re(ZN(0,0))", "Im(ZN(0,0))",
      "Re(ZN(0,1))", "Im(ZN(0,1))", "Re(ZN(0,2))", "Im(ZN(0,2))", "Re(ZN(0,3))",
      "Im(ZN(0,3))", "Re(ZN(1,0))", "Im(ZN(1,0))", "Re(ZN(1,1))", "Im(ZN(1,1))",
      "Re(ZN(1,2))", "Im(ZN(1,2))", "Re(ZN(1,3))", "Im(ZN(1,3))", "Re(ZN(2,0))",
      "Im(ZN(2,0))", "Re(ZN(2,1))", "Im(ZN(2,1))", "Re(ZN(2,2))", "Im(ZN(2,2))",
      "Re(ZN(2,3))", "Im(ZN(2,3))", "Re(ZN(3,0))", "Im(ZN(3,0))", "Re(ZN(3,1))",
      "Im(ZN(3,1))", "Re(ZN(3,2))", "Im(ZN(3,2))", "Re(ZN(3,3))", "Im(ZN(3,3))",
      "Re(UM(0,0))", "Im(UM(0,0))", "Re(UM(0,1))", "Im(UM(0,1))", "Re(UM(1,0))",
      "Im(UM(1,0))", "Re(UM(1,1))", "Im(UM(1,1))", "Re(UP(0,0))", "Im(UP(0,0))",
      "Re(UP(0,1))", "Im(UP(0,1))", "Re(UP(1,0))", "Im(UP(1,0))", "Re(UP(1,1))",
      "Im(UP(1,1))", "ZZ(0,0)", "ZZ(0,1)", "ZZ(1,0)", "ZZ(1,1)"};

   const std::array<std::string, NUMBER_OF_INPUT_PARAMETERS> input_parameter_names
       = {"MSUSY", "M1Input", "M2Input", "M3Input", "MuInput", "mAInput", "MEWSB",
      "AtInput", "TanBeta", "LambdaLoopOrder", "msq2(0,0)", "msq2(0,1)",
      "msq2(0,2)", "msq2(1,0)", "msq2(1,1)", "msq2(1,2)", "msq2(2,0)", "msq2(2,1)"
      , "msq2(2,2)", "msu2(0,0)", "msu2(0,1)", "msu2(0,2)", "msu2(1,0)",
      "msu2(1,1)", "msu2(1,2)", "msu2(2,0)", "msu2(2,1)", "msu2(2,2)", "msd2(0,0)"
      , "msd2(0,1)", "msd2(0,2)", "msd2(1,0)", "msd2(1,1)", "msd2(1,2)",
      "msd2(2,0)", "msd2(2,1)", "msd2(2,2)", "msl2(0,0)", "msl2(0,1)", "msl2(0,2)"
      , "msl2(1,0)", "msl2(1,1)", "msl2(1,2)", "msl2(2,0)", "msl2(2,1)",
      "msl2(2,2)", "mse2(0,0)", "mse2(0,1)", "mse2(0,2)", "mse2(1,0)", "mse2(1,1)"
      , "mse2(1,2)", "mse2(2,0)", "mse2(2,1)", "mse2(2,2)"};

   const std::array<std::string, NUMBER_OF_EXTRA_PARAMETERS> extra_parameter_names
       = {};

   const std::string model_name = "SplitMSSM";

int get_pdg_code_for_particle(Particles p)
{
   if (particle_multiplicities[p] > 1) {
      throw OutOfBoundsError(particle_names[p] + " must have a generation index");
   }

   int pdg = 0;
   switch (p) {

   case VG: pdg = 21; break;
   case Hp: pdg = 0; break;
   case Glu: pdg = 1000021; break;
   case Ah: pdg = 0; break;
   case hh: pdg = 25; break;
   case VWp: pdg = 24; break;
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
   case Fd: pdg_codes = {1, 3, 5}; break;
   case Fu: pdg_codes = {2, 4, 6}; break;
   case Fe: pdg_codes = {11, 13, 15}; break;
   case Chi: pdg_codes = {1000022, 1000023, 1000025, 1000035}; break;
   case Cha: pdg_codes = {-1000024, -1000037}; break;

   default: throw OutOfBoundsError("invalid particle " + std::to_string(p));
   }

   const int n_pdgs = pdg_codes.size();
   if (index < 0 || index >= n_pdgs) {
      throw OutOfBoundsError("index " + std::to_string(index) + " out of bounds");
   }

   return pdg_codes[index];
}

std::string get_particle_name_from_pdg(int pdg)
{
   std::string name;

   switch (pdg) {

   case 21: name = "VG"; break;
   case 12: name = "Fv(1)"; break;
   case 14: name = "Fv(2)"; break;
   case 16: name = "Fv(3)"; break;
   case 1000021: name = "Glu"; break;
   case 25: name = "hh"; break;
   case 1: name = "Fd(1)"; break;
   case 3: name = "Fd(2)"; break;
   case 5: name = "Fd(3)"; break;
   case 2: name = "Fu(1)"; break;
   case 4: name = "Fu(2)"; break;
   case 6: name = "Fu(3)"; break;
   case 11: name = "Fe(1)"; break;
   case 13: name = "Fe(2)"; break;
   case 15: name = "Fe(3)"; break;
   case 1000022: name = "Chi(1)"; break;
   case 1000023: name = "Chi(2)"; break;
   case 1000025: name = "Chi(3)"; break;
   case 1000035: name = "Chi(4)"; break;
   case -1000024: name = "Cha(1)"; break;
   case -1000037: name = "Cha(2)"; break;
   case 24: name = "VWp"; break;
   case 22: name = "VP"; break;
   case 23: name = "VZ"; break;
   case -12: name = "barFv(1)"; break;
   case -14: name = "barFv(2)"; break;
   case -16: name = "barFv(3)"; break;
   case -1: name = "barFd(1)"; break;
   case -3: name = "barFd(2)"; break;
   case -5: name = "barFd(3)"; break;
   case -2: name = "barFu(1)"; break;
   case -4: name = "barFu(2)"; break;
   case -6: name = "barFu(3)"; break;
   case -11: name = "barFe(1)"; break;
   case -13: name = "barFe(2)"; break;
   case -15: name = "barFe(3)"; break;
   case 1000024: name = "barCha(1)"; break;
   case 1000037: name = "barCha(2)"; break;
   case -24: name = "conjVWp"; break;

   default: name = "";
   }

   return name;
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

} // namespace SplitMSSM_info

} // namespace flexiblesusy

