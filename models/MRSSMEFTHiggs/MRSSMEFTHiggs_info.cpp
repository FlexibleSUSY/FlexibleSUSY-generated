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

// File generated at Tue 10 Oct 2017 21:39:58

#include "MRSSMEFTHiggs_info.hpp"

#include <iostream>

namespace flexiblesusy {

namespace MRSSMEFTHiggs_info {
   const double normalization_g1 = 0.7745966692414834;
   const double normalization_g2 = 1;
   const double normalization_g3 = 1;

   const std::array<int, NUMBER_OF_PARTICLES> particle_multiplicities = {1, 1,
      3, 1, 1, 1, 1, 6, 3, 6, 6, 4, 4, 2, 4, 4, 2, 2, 3, 3, 3, 1, 1, 1};

   const std::array<std::string, NUMBER_OF_PARTICLES> particle_names = {"VG",
      "Glu", "Fv", "SRdp", "SRum", "sigmaO", "phiO", "Sd", "Sv", "Su", "Se", "hh",
      "Ah", "Rh", "Hpm", "Chi", "Cha1", "Cha2", "Fe", "Fd", "Fu", "VWm", "VP",
      "VZ"};

   const std::array<std::string, NUMBER_OF_PARTICLES> particle_latex_names = {
      "g", "\\tilde{g}", "\\nu", "R_d^+", "R_u^-", "\\sigma_o", "\\phi_o",
      "\\tilde{d}", "\\tilde{\\nu}", "\\tilde{u}", "\\tilde{e}", "h", "A^0", "R^h"
      , "H^-", "\\tilde{\\chi}^0", "\\tilde{\\chi}^+", "\\tilde{\\rho}^-", "e",
      "d", "u", "W^-", "\\gamma", "Z"};

   const std::array<std::string, NUMBER_OF_PARAMETERS> parameter_names = {
      "Yd(0,0)", "Yd(0,1)", "Yd(0,2)", "Yd(1,0)", "Yd(1,1)", "Yd(1,2)", "Yd(2,0)",
      "Yd(2,1)", "Yd(2,2)", "Ye(0,0)", "Ye(0,1)", "Ye(0,2)", "Ye(1,0)", "Ye(1,1)"
      , "Ye(1,2)", "Ye(2,0)", "Ye(2,1)", "Ye(2,2)", "LamTD", "LamTU", "LamSD",
      "LamSU", "Yu(0,0)", "Yu(0,1)", "Yu(0,2)", "Yu(1,0)", "Yu(1,1)", "Yu(1,2)",
      "Yu(2,0)", "Yu(2,1)", "Yu(2,2)", "Mu", "MuD", "MuU", "g1", "g2", "g3", "vd",
      "vu", "vT", "vS", "BMu", "BMuD", "BMuU", "mq2(0,0)", "mq2(0,1)", "mq2(0,2)"
      , "mq2(1,0)", "mq2(1,1)", "mq2(1,2)", "mq2(2,0)", "mq2(2,1)", "mq2(2,2)",
      "ml2(0,0)", "ml2(0,1)", "ml2(0,2)", "ml2(1,0)", "ml2(1,1)", "ml2(1,2)",
      "ml2(2,0)", "ml2(2,1)", "ml2(2,2)", "mHd2", "mHu2", "md2(0,0)", "md2(0,1)",
      "md2(0,2)", "md2(1,0)", "md2(1,1)", "md2(1,2)", "md2(2,0)", "md2(2,1)",
      "md2(2,2)", "mu2(0,0)", "mu2(0,1)", "mu2(0,2)", "mu2(1,0)", "mu2(1,1)",
      "mu2(1,2)", "mu2(2,0)", "mu2(2,1)", "mu2(2,2)", "me2(0,0)", "me2(0,1)",
      "me2(0,2)", "me2(1,0)", "me2(1,1)", "me2(1,2)", "me2(2,0)", "me2(2,1)",
      "me2(2,2)", "mS2", "mT2", "moc2", "mRd2", "mRu2", "MDBS", "MDWBT", "MDGoc"};

   const std::array<std::string, NUMBER_OF_MIXINGS> particle_mixing_names = {
      "ZD(0,0)", "ZD(0,1)", "ZD(0,2)", "ZD(0,3)", "ZD(0,4)", "ZD(0,5)", "ZD(1,0)"
      , "ZD(1,1)", "ZD(1,2)", "ZD(1,3)", "ZD(1,4)", "ZD(1,5)", "ZD(2,0)",
      "ZD(2,1)", "ZD(2,2)", "ZD(2,3)", "ZD(2,4)", "ZD(2,5)", "ZD(3,0)", "ZD(3,1)",
      "ZD(3,2)", "ZD(3,3)", "ZD(3,4)", "ZD(3,5)", "ZD(4,0)", "ZD(4,1)", "ZD(4,2)"
      , "ZD(4,3)", "ZD(4,4)", "ZD(4,5)", "ZD(5,0)", "ZD(5,1)", "ZD(5,2)",
      "ZD(5,3)", "ZD(5,4)", "ZD(5,5)", "ZV(0,0)", "ZV(0,1)", "ZV(0,2)", "ZV(1,0)",
      "ZV(1,1)", "ZV(1,2)", "ZV(2,0)", "ZV(2,1)", "ZV(2,2)", "ZU(0,0)", "ZU(0,1)"
      , "ZU(0,2)", "ZU(0,3)", "ZU(0,4)", "ZU(0,5)", "ZU(1,0)", "ZU(1,1)",
      "ZU(1,2)", "ZU(1,3)", "ZU(1,4)", "ZU(1,5)", "ZU(2,0)", "ZU(2,1)", "ZU(2,2)",
      "ZU(2,3)", "ZU(2,4)", "ZU(2,5)", "ZU(3,0)", "ZU(3,1)", "ZU(3,2)", "ZU(3,3)"
      , "ZU(3,4)", "ZU(3,5)", "ZU(4,0)", "ZU(4,1)", "ZU(4,2)", "ZU(4,3)",
      "ZU(4,4)", "ZU(4,5)", "ZU(5,0)", "ZU(5,1)", "ZU(5,2)", "ZU(5,3)", "ZU(5,4)",
      "ZU(5,5)", "ZE(0,0)", "ZE(0,1)", "ZE(0,2)", "ZE(0,3)", "ZE(0,4)", "ZE(0,5)"
      , "ZE(1,0)", "ZE(1,1)", "ZE(1,2)", "ZE(1,3)", "ZE(1,4)", "ZE(1,5)",
      "ZE(2,0)", "ZE(2,1)", "ZE(2,2)", "ZE(2,3)", "ZE(2,4)", "ZE(2,5)", "ZE(3,0)",
      "ZE(3,1)", "ZE(3,2)", "ZE(3,3)", "ZE(3,4)", "ZE(3,5)", "ZE(4,0)", "ZE(4,1)"
      , "ZE(4,2)", "ZE(4,3)", "ZE(4,4)", "ZE(4,5)", "ZE(5,0)", "ZE(5,1)",
      "ZE(5,2)", "ZE(5,3)", "ZE(5,4)", "ZE(5,5)", "ZH(0,0)", "ZH(0,1)", "ZH(0,2)",
      "ZH(0,3)", "ZH(1,0)", "ZH(1,1)", "ZH(1,2)", "ZH(1,3)", "ZH(2,0)", "ZH(2,1)"
      , "ZH(2,2)", "ZH(2,3)", "ZH(3,0)", "ZH(3,1)", "ZH(3,2)", "ZH(3,3)",
      "ZA(0,0)", "ZA(0,1)", "ZA(0,2)", "ZA(0,3)", "ZA(1,0)", "ZA(1,1)", "ZA(1,2)",
      "ZA(1,3)", "ZA(2,0)", "ZA(2,1)", "ZA(2,2)", "ZA(2,3)", "ZA(3,0)", "ZA(3,1)"
      , "ZA(3,2)", "ZA(3,3)", "ZHR(0,0)", "ZHR(0,1)", "ZHR(1,0)", "ZHR(1,1)",
      "ZP(0,0)", "ZP(0,1)", "ZP(0,2)", "ZP(0,3)", "ZP(1,0)", "ZP(1,1)", "ZP(1,2)",
      "ZP(1,3)", "ZP(2,0)", "ZP(2,1)", "ZP(2,2)", "ZP(2,3)", "ZP(3,0)", "ZP(3,1)"
      , "ZP(3,2)", "ZP(3,3)", "Re(ZN1(0,0))", "Im(ZN1(0,0))", "Re(ZN1(0,1))",
      "Im(ZN1(0,1))", "Re(ZN1(0,2))", "Im(ZN1(0,2))", "Re(ZN1(0,3))",
      "Im(ZN1(0,3))", "Re(ZN1(1,0))", "Im(ZN1(1,0))", "Re(ZN1(1,1))",
      "Im(ZN1(1,1))", "Re(ZN1(1,2))", "Im(ZN1(1,2))", "Re(ZN1(1,3))",
      "Im(ZN1(1,3))", "Re(ZN1(2,0))", "Im(ZN1(2,0))", "Re(ZN1(2,1))",
      "Im(ZN1(2,1))", "Re(ZN1(2,2))", "Im(ZN1(2,2))", "Re(ZN1(2,3))",
      "Im(ZN1(2,3))", "Re(ZN1(3,0))", "Im(ZN1(3,0))", "Re(ZN1(3,1))",
      "Im(ZN1(3,1))", "Re(ZN1(3,2))", "Im(ZN1(3,2))", "Re(ZN1(3,3))",
      "Im(ZN1(3,3))", "Re(ZN2(0,0))", "Im(ZN2(0,0))", "Re(ZN2(0,1))",
      "Im(ZN2(0,1))", "Re(ZN2(0,2))", "Im(ZN2(0,2))", "Re(ZN2(0,3))",
      "Im(ZN2(0,3))", "Re(ZN2(1,0))", "Im(ZN2(1,0))", "Re(ZN2(1,1))",
      "Im(ZN2(1,1))", "Re(ZN2(1,2))", "Im(ZN2(1,2))", "Re(ZN2(1,3))",
      "Im(ZN2(1,3))", "Re(ZN2(2,0))", "Im(ZN2(2,0))", "Re(ZN2(2,1))",
      "Im(ZN2(2,1))", "Re(ZN2(2,2))", "Im(ZN2(2,2))", "Re(ZN2(2,3))",
      "Im(ZN2(2,3))", "Re(ZN2(3,0))", "Im(ZN2(3,0))", "Re(ZN2(3,1))",
      "Im(ZN2(3,1))", "Re(ZN2(3,2))", "Im(ZN2(3,2))", "Re(ZN2(3,3))",
      "Im(ZN2(3,3))", "Re(UM1(0,0))", "Im(UM1(0,0))", "Re(UM1(0,1))",
      "Im(UM1(0,1))", "Re(UM1(1,0))", "Im(UM1(1,0))", "Re(UM1(1,1))",
      "Im(UM1(1,1))", "Re(UP1(0,0))", "Im(UP1(0,0))", "Re(UP1(0,1))",
      "Im(UP1(0,1))", "Re(UP1(1,0))", "Im(UP1(1,0))", "Re(UP1(1,1))",
      "Im(UP1(1,1))", "Re(UM2(0,0))", "Im(UM2(0,0))", "Re(UM2(0,1))",
      "Im(UM2(0,1))", "Re(UM2(1,0))", "Im(UM2(1,0))", "Re(UM2(1,1))",
      "Im(UM2(1,1))", "Re(UP2(0,0))", "Im(UP2(0,0))", "Re(UP2(0,1))",
      "Im(UP2(0,1))", "Re(UP2(1,0))", "Im(UP2(1,0))", "Re(UP2(1,1))",
      "Im(UP2(1,1))", "Re(ZEL(0,0))", "Im(ZEL(0,0))", "Re(ZEL(0,1))",
      "Im(ZEL(0,1))", "Re(ZEL(0,2))", "Im(ZEL(0,2))", "Re(ZEL(1,0))",
      "Im(ZEL(1,0))", "Re(ZEL(1,1))", "Im(ZEL(1,1))", "Re(ZEL(1,2))",
      "Im(ZEL(1,2))", "Re(ZEL(2,0))", "Im(ZEL(2,0))", "Re(ZEL(2,1))",
      "Im(ZEL(2,1))", "Re(ZEL(2,2))", "Im(ZEL(2,2))", "Re(ZER(0,0))",
      "Im(ZER(0,0))", "Re(ZER(0,1))", "Im(ZER(0,1))", "Re(ZER(0,2))",
      "Im(ZER(0,2))", "Re(ZER(1,0))", "Im(ZER(1,0))", "Re(ZER(1,1))",
      "Im(ZER(1,1))", "Re(ZER(1,2))", "Im(ZER(1,2))", "Re(ZER(2,0))",
      "Im(ZER(2,0))", "Re(ZER(2,1))", "Im(ZER(2,1))", "Re(ZER(2,2))",
      "Im(ZER(2,2))", "Re(ZDL(0,0))", "Im(ZDL(0,0))", "Re(ZDL(0,1))",
      "Im(ZDL(0,1))", "Re(ZDL(0,2))", "Im(ZDL(0,2))", "Re(ZDL(1,0))",
      "Im(ZDL(1,0))", "Re(ZDL(1,1))", "Im(ZDL(1,1))", "Re(ZDL(1,2))",
      "Im(ZDL(1,2))", "Re(ZDL(2,0))", "Im(ZDL(2,0))", "Re(ZDL(2,1))",
      "Im(ZDL(2,1))", "Re(ZDL(2,2))", "Im(ZDL(2,2))", "Re(ZDR(0,0))",
      "Im(ZDR(0,0))", "Re(ZDR(0,1))", "Im(ZDR(0,1))", "Re(ZDR(0,2))",
      "Im(ZDR(0,2))", "Re(ZDR(1,0))", "Im(ZDR(1,0))", "Re(ZDR(1,1))",
      "Im(ZDR(1,1))", "Re(ZDR(1,2))", "Im(ZDR(1,2))", "Re(ZDR(2,0))",
      "Im(ZDR(2,0))", "Re(ZDR(2,1))", "Im(ZDR(2,1))", "Re(ZDR(2,2))",
      "Im(ZDR(2,2))", "Re(ZUL(0,0))", "Im(ZUL(0,0))", "Re(ZUL(0,1))",
      "Im(ZUL(0,1))", "Re(ZUL(0,2))", "Im(ZUL(0,2))", "Re(ZUL(1,0))",
      "Im(ZUL(1,0))", "Re(ZUL(1,1))", "Im(ZUL(1,1))", "Re(ZUL(1,2))",
      "Im(ZUL(1,2))", "Re(ZUL(2,0))", "Im(ZUL(2,0))", "Re(ZUL(2,1))",
      "Im(ZUL(2,1))", "Re(ZUL(2,2))", "Im(ZUL(2,2))", "Re(ZUR(0,0))",
      "Im(ZUR(0,0))", "Re(ZUR(0,1))", "Im(ZUR(0,1))", "Re(ZUR(0,2))",
      "Im(ZUR(0,2))", "Re(ZUR(1,0))", "Im(ZUR(1,0))", "Re(ZUR(1,1))",
      "Im(ZUR(1,1))", "Re(ZUR(1,2))", "Im(ZUR(1,2))", "Re(ZUR(2,0))",
      "Im(ZUR(2,0))", "Re(ZUR(2,1))", "Im(ZUR(2,1))", "Re(ZUR(2,2))",
      "Im(ZUR(2,2))", "ZZ(0,0)", "ZZ(0,1)", "ZZ(1,0)", "ZZ(1,1)"};

   const std::array<std::string, NUMBER_OF_INPUT_PARAMETERS>
      input_parameter_names = {"TanBeta", "MS", "LamTDInput", "LamTUInput",
      "LamSDInput", "LamSUInput", "MuDInput", "MuUInput", "BMuInput",
      "mq2Input(0,0)", "mq2Input(0,1)", "mq2Input(0,2)", "mq2Input(1,0)",
      "mq2Input(1,1)", "mq2Input(1,2)", "mq2Input(2,0)", "mq2Input(2,1)",
      "mq2Input(2,2)", "ml2Input(0,0)", "ml2Input(0,1)", "ml2Input(0,2)",
      "ml2Input(1,0)", "ml2Input(1,1)", "ml2Input(1,2)", "ml2Input(2,0)",
      "ml2Input(2,1)", "ml2Input(2,2)", "md2Input(0,0)", "md2Input(0,1)",
      "md2Input(0,2)", "md2Input(1,0)", "md2Input(1,1)", "md2Input(1,2)",
      "md2Input(2,0)", "md2Input(2,1)", "md2Input(2,2)", "mu2Input(0,0)",
      "mu2Input(0,1)", "mu2Input(0,2)", "mu2Input(1,0)", "mu2Input(1,1)",
      "mu2Input(1,2)", "mu2Input(2,0)", "mu2Input(2,1)", "mu2Input(2,2)",
      "me2Input(0,0)", "me2Input(0,1)", "me2Input(0,2)", "me2Input(1,0)",
      "me2Input(1,1)", "me2Input(1,2)", "me2Input(2,0)", "me2Input(2,1)",
      "me2Input(2,2)", "mS2Input", "mT2Input", "moc2Input", "mRd2Input",
      "mRu2Input", "MDBSInput", "MDWBTInput", "MDGocInput"};

   const std::array<std::string, NUMBER_OF_EXTRA_PARAMETERS>
      extra_parameter_names = {};

   const std::string model_name = "MRSSMEFTHiggs";

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

} // namespace MRSSMEFTHiggs_info

} // namespace flexiblesusy
