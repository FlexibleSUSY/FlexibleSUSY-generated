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

// File generated at Sun 18 Oct 2015 11:51:19

#include "MRSSM_info.hpp"

#include <iostream>

namespace flexiblesusy {

namespace MRSSM_info {
   const double normalization_g1 = 0.7745966692414834;
   const double normalization_g2 = 1;
   const double normalization_g3 = 1;

   const unsigned particle_multiplicities[NUMBER_OF_PARTICLES] = {1, 3, 1, 1, 1
      , 1, 1, 6, 3, 6, 6, 4, 4, 2, 4, 4, 2, 2, 3, 3, 3, 1, 1, 1};

   const char* particle_names[NUMBER_OF_PARTICLES] = {"Glu", "Fv", "SRdp",
      "SRum", "sigmaO", "phiO", "VZ", "Sd", "Sv", "Su", "Se", "hh", "Ah", "Rh",
      "Hpm", "Chi", "Cha1", "Cha2", "Fe", "Fd", "Fu", "VG", "VP", "VWm"};

   const char* particle_latex_names[NUMBER_OF_PARTICLES] = {   "\\tilde{g}",
      "\\nu", "R_d^+", "R_u^-", "\\sigma_o", "\\phi_o", "Z", "\\tilde{d}",
      "\\tilde{\\nu}", "\\tilde{u}", "\\tilde{e}", "h", "A^0", "R^h", "H^-",
      "\\tilde{\\chi}^0", "\\tilde{\\chi}^+", "\\tilde{\\rho}^-", "e", "d", "u",
      "g", "\\gamma", "W^-"};

   const char* parameter_names[NUMBER_OF_PARAMETERS] = {"Yd(0,0)", "Yd(0,1)",
      "Yd(0,2)", "Yd(1,0)", "Yd(1,1)", "Yd(1,2)", "Yd(2,0)", "Yd(2,1)", "Yd(2,2)",
      "Ye(0,0)", "Ye(0,1)", "Ye(0,2)", "Ye(1,0)", "Ye(1,1)", "Ye(1,2)", "Ye(2,0)"
      , "Ye(2,1)", "Ye(2,2)", "LamTD", "LamTU", "LamSD", "LamSU", "Yu(0,0)",
      "Yu(0,1)", "Yu(0,2)", "Yu(1,0)", "Yu(1,1)", "Yu(1,2)", "Yu(2,0)", "Yu(2,1)",
      "Yu(2,2)", "Mu", "MuD", "MuU", "g1", "g2", "g3", "vd", "vu", "vT", "vS",
      "BMu", "BMuD", "BMuU", "mq2(0,0)", "mq2(0,1)", "mq2(0,2)", "mq2(1,0)",
      "mq2(1,1)", "mq2(1,2)", "mq2(2,0)", "mq2(2,1)", "mq2(2,2)", "ml2(0,0)",
      "ml2(0,1)", "ml2(0,2)", "ml2(1,0)", "ml2(1,1)", "ml2(1,2)", "ml2(2,0)",
      "ml2(2,1)", "ml2(2,2)", "mHd2", "mHu2", "md2(0,0)", "md2(0,1)", "md2(0,2)",
      "md2(1,0)", "md2(1,1)", "md2(1,2)", "md2(2,0)", "md2(2,1)", "md2(2,2)",
      "mu2(0,0)", "mu2(0,1)", "mu2(0,2)", "mu2(1,0)", "mu2(1,1)", "mu2(1,2)",
      "mu2(2,0)", "mu2(2,1)", "mu2(2,2)", "me2(0,0)", "me2(0,1)", "me2(0,2)",
      "me2(1,0)", "me2(1,1)", "me2(1,2)", "me2(2,0)", "me2(2,1)", "me2(2,2)",
      "mS2", "mT2", "moc2", "mRd2", "mRu2", "MDBS", "MDWBT", "MDGoc"};

   const char* model_name = "MRSSM";
   const bool is_low_energy_model = true;
   const bool is_supersymmetric_model = true;

void print(std::ostream& ostr)
{
   ostr
      << "Model information\n"
      << "=================\n"
      << "Model name:                " << model_name << '\n'
      << "Is a low-energy model:     "
      << (is_low_energy_model ? "yes" : "no") << '\n'
      << "Is a supersymmetric model: "
      << (is_supersymmetric_model ? "yes" : "no") << '\n'
      << "Number of multiplets:      " << NUMBER_OF_PARTICLES << '\n'
      << "Number of parameters:      " << NUMBER_OF_PARAMETERS << '\n'
      ;

   ostr << "\n"
      "Multiplets:                ";
   for (unsigned i = 0; i < NUMBER_OF_PARTICLES; i++) {
      ostr << particle_names[i]
           << '[' << particle_multiplicities[i] << ']';
      if (i + 1 < NUMBER_OF_PARTICLES)
         ostr << ", ";
   }

   ostr << "\n\n"
      "Parameters:                ";
   for (unsigned i = 0; i < NUMBER_OF_PARAMETERS; i++) {
      ostr << parameter_names[i];
      if (i + 1 < NUMBER_OF_PARAMETERS)
         ostr << ", ";
   }
   ostr << '\n';
}

} // namespace MRSSM_info

} // namespace flexiblesusy

