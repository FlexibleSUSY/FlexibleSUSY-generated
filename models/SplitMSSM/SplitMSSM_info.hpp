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

// File generated at Fri 8 Jan 2016 15:07:05

#ifndef SplitMSSM_INFO_H
#define SplitMSSM_INFO_H

#include <iosfwd>

namespace flexiblesusy {

namespace SplitMSSM_info {
   enum Particles : unsigned {VG, Hp, Fv, Glu, Ah, hh, VP, VZ, Fd, Fu, Fe, Chi,
      Cha, VWp, NUMBER_OF_PARTICLES};

   enum Parameters : unsigned {g1, g2, g3, Lambdax, Yu00, Yu01, Yu02, Yu10,
      Yu11, Yu12, Yu20, Yu21, Yu22, Yd00, Yd01, Yd02, Yd10, Yd11, Yd12, Yd20, Yd21
      , Yd22, Ye00, Ye01, Ye02, Ye10, Ye11, Ye12, Ye20, Ye21, Ye22, gYd, g2d, gYu,
      g2u, MassB, MassG, MassWB, Mu, mu2, v, NUMBER_OF_PARAMETERS};

   enum Mixings : unsigned {ReVd00, ImVd00, ReVd01, ImVd01, ReVd02, ImVd02,
      ReVd10, ImVd10, ReVd11, ImVd11, ReVd12, ImVd12, ReVd20, ImVd20, ReVd21,
      ImVd21, ReVd22, ImVd22, ReUd00, ImUd00, ReUd01, ImUd01, ReUd02, ImUd02,
      ReUd10, ImUd10, ReUd11, ImUd11, ReUd12, ImUd12, ReUd20, ImUd20, ReUd21,
      ImUd21, ReUd22, ImUd22, ReVu00, ImVu00, ReVu01, ImVu01, ReVu02, ImVu02,
      ReVu10, ImVu10, ReVu11, ImVu11, ReVu12, ImVu12, ReVu20, ImVu20, ReVu21,
      ImVu21, ReVu22, ImVu22, ReUu00, ImUu00, ReUu01, ImUu01, ReUu02, ImUu02,
      ReUu10, ImUu10, ReUu11, ImUu11, ReUu12, ImUu12, ReUu20, ImUu20, ReUu21,
      ImUu21, ReUu22, ImUu22, ReVe00, ImVe00, ReVe01, ImVe01, ReVe02, ImVe02,
      ReVe10, ImVe10, ReVe11, ImVe11, ReVe12, ImVe12, ReVe20, ImVe20, ReVe21,
      ImVe21, ReVe22, ImVe22, ReUe00, ImUe00, ReUe01, ImUe01, ReUe02, ImUe02,
      ReUe10, ImUe10, ReUe11, ImUe11, ReUe12, ImUe12, ReUe20, ImUe20, ReUe21,
      ImUe21, ReUe22, ImUe22, ReZN00, ImZN00, ReZN01, ImZN01, ReZN02, ImZN02,
      ReZN03, ImZN03, ReZN10, ImZN10, ReZN11, ImZN11, ReZN12, ImZN12, ReZN13,
      ImZN13, ReZN20, ImZN20, ReZN21, ImZN21, ReZN22, ImZN22, ReZN23, ImZN23,
      ReZN30, ImZN30, ReZN31, ImZN31, ReZN32, ImZN32, ReZN33, ImZN33, ReUM00,
      ImUM00, ReUM01, ImUM01, ReUM10, ImUM10, ReUM11, ImUM11, ReUP00, ImUP00,
      ReUP01, ImUP01, ReUP10, ImUP10, ReUP11, ImUP11, NUMBER_OF_MIXINGS};

   enum Input_parameters : unsigned {MSUSY, M1Input, M2Input, M3Input, MuInput,
      mAInput, MEWSB, AtInput, TanBeta, msq200, msq201, msq202, msq210, msq211,
      msq212, msq220, msq221, msq222, msu200, msu201, msu202, msu210, msu211,
      msu212, msu220, msu221, msu222, msd200, msd201, msd202, msd210, msd211,
      msd212, msd220, msd221, msd222, msl200, msl201, msl202, msl210, msl211,
      msl212, msl220, msl221, msl222, mse200, mse201, mse202, mse210, mse211,
      mse212, mse220, mse221, mse222, NUMBER_OF_INPUT_PARAMETERS};

   extern const double normalization_g1;
   extern const double normalization_g2;
   extern const double normalization_g3;

   extern const unsigned particle_multiplicities[NUMBER_OF_PARTICLES];
   extern const char* particle_names[NUMBER_OF_PARTICLES];
   extern const char* particle_latex_names[NUMBER_OF_PARTICLES];
   extern const char* parameter_names[NUMBER_OF_PARAMETERS];
   extern const char* particle_mixing_names[NUMBER_OF_MIXINGS];
   extern const char* input_parameter_names[NUMBER_OF_INPUT_PARAMETERS];
   extern const char* model_name;
   extern const bool is_low_energy_model;
   extern const bool is_supersymmetric_model;

   void print(std::ostream&);
}

} // namespace flexiblesusy

#endif
