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

// File generated at Sun 26 Aug 2018 14:07:31

#ifndef HGTHDMIIMSSMBC_INFO_H
#define HGTHDMIIMSSMBC_INFO_H

#include "problems.hpp"

#include <array>
#include <iosfwd>
#include <string>

namespace flexiblesusy {

namespace HGTHDMIIMSSMBC_info {
   enum Particles : int { VG, Fv, Glu, hh, Ah, Hm, Fd, Fu, Fe, Chi, Cha, VWm, VP,
      VZ, NUMBER_OF_PARTICLES };

   enum Masses : int { MVG, MFv_1, MFv_2, MFv_3, MGlu, Mhh_1, Mhh_2, MAh_1, MAh_2,
      MHm_1, MHm_2, MFd_1, MFd_2, MFd_3, MFu_1, MFu_2, MFu_3, MFe_1, MFe_2, MFe_3,
      MChi_1, MChi_2, MChi_3, MChi_4, MCha_1, MCha_2, MVWm, MVP, MVZ,
      NUMBER_OF_MASSES };

   enum Parameters : int { g1, g2, g3, Lambda6, Lambda5, Lambda7, Lambda1, Lambda4
      , Lambda3, Lambda2, Yu0_0, Yu0_1, Yu0_2, Yu1_0, Yu1_1, Yu1_2, Yu2_0, Yu2_1,
      Yu2_2, Yd0_0, Yd0_1, Yd0_2, Yd1_0, Yd1_1, Yd1_2, Yd2_0, Yd2_1, Yd2_2, Ye0_0,
      Ye0_1, Ye0_2, Ye1_0, Ye1_1, Ye1_2, Ye2_0, Ye2_1, Ye2_2, g1dp, g1d, g2up, g2u
      , MassB, MassG, MassWB, Mu, M122, M112, M222, v1, v2, NUMBER_OF_PARAMETERS }
      ;

   enum Mixings : int { ZH0_0, ZH0_1, ZH1_0, ZH1_1, ZA0_0, ZA0_1, ZA1_0, ZA1_1,
      ZP0_0, ZP0_1, ZP1_0, ZP1_1, ReVd0_0, ImVd0_0, ReVd0_1, ImVd0_1, ReVd0_2,
      ImVd0_2, ReVd1_0, ImVd1_0, ReVd1_1, ImVd1_1, ReVd1_2, ImVd1_2, ReVd2_0,
      ImVd2_0, ReVd2_1, ImVd2_1, ReVd2_2, ImVd2_2, ReUd0_0, ImUd0_0, ReUd0_1,
      ImUd0_1, ReUd0_2, ImUd0_2, ReUd1_0, ImUd1_0, ReUd1_1, ImUd1_1, ReUd1_2,
      ImUd1_2, ReUd2_0, ImUd2_0, ReUd2_1, ImUd2_1, ReUd2_2, ImUd2_2, ReVu0_0,
      ImVu0_0, ReVu0_1, ImVu0_1, ReVu0_2, ImVu0_2, ReVu1_0, ImVu1_0, ReVu1_1,
      ImVu1_1, ReVu1_2, ImVu1_2, ReVu2_0, ImVu2_0, ReVu2_1, ImVu2_1, ReVu2_2,
      ImVu2_2, ReUu0_0, ImUu0_0, ReUu0_1, ImUu0_1, ReUu0_2, ImUu0_2, ReUu1_0,
      ImUu1_0, ReUu1_1, ImUu1_1, ReUu1_2, ImUu1_2, ReUu2_0, ImUu2_0, ReUu2_1,
      ImUu2_1, ReUu2_2, ImUu2_2, ReVe0_0, ImVe0_0, ReVe0_1, ImVe0_1, ReVe0_2,
      ImVe0_2, ReVe1_0, ImVe1_0, ReVe1_1, ImVe1_1, ReVe1_2, ImVe1_2, ReVe2_0,
      ImVe2_0, ReVe2_1, ImVe2_1, ReVe2_2, ImVe2_2, ReUe0_0, ImUe0_0, ReUe0_1,
      ImUe0_1, ReUe0_2, ImUe0_2, ReUe1_0, ImUe1_0, ReUe1_1, ImUe1_1, ReUe1_2,
      ImUe1_2, ReUe2_0, ImUe2_0, ReUe2_1, ImUe2_1, ReUe2_2, ImUe2_2, ReZN0_0,
      ImZN0_0, ReZN0_1, ImZN0_1, ReZN0_2, ImZN0_2, ReZN0_3, ImZN0_3, ReZN1_0,
      ImZN1_0, ReZN1_1, ImZN1_1, ReZN1_2, ImZN1_2, ReZN1_3, ImZN1_3, ReZN2_0,
      ImZN2_0, ReZN2_1, ImZN2_1, ReZN2_2, ImZN2_2, ReZN2_3, ImZN2_3, ReZN3_0,
      ImZN3_0, ReZN3_1, ImZN3_1, ReZN3_2, ImZN3_2, ReZN3_3, ImZN3_3, ReUM0_0,
      ImUM0_0, ReUM0_1, ImUM0_1, ReUM1_0, ImUM1_0, ReUM1_1, ImUM1_1, ReUP0_0,
      ImUP0_0, ReUP0_1, ImUP0_1, ReUP1_0, ImUP1_0, ReUP1_1, ImUP1_1, ZZ0_0, ZZ0_1,
      ZZ1_0, ZZ1_1, NUMBER_OF_MIXINGS };

   enum Input_parameters : int { TanBeta, MSUSY, MEWSB, MuInput, M1Input, M2Input,
      M3Input, MAInput, AtInput, AbInput, AtauInput, LambdaLoopOrder,
      NUMBER_OF_INPUT_PARAMETERS };

   enum Extra_parameters : int { NUMBER_OF_EXTRA_PARAMETERS };

   extern const double normalization_g1;
   extern const double normalization_g2;
   extern const double normalization_g3;

   extern const std::array<int, NUMBER_OF_PARTICLES> particle_multiplicities;
   extern const std::array<std::string, NUMBER_OF_PARTICLES> particle_names;
   extern const std::array<std::string, NUMBER_OF_PARTICLES> particle_latex_names;
   extern const std::array<std::string, NUMBER_OF_PARAMETERS> parameter_names;
   extern const std::array<std::string, NUMBER_OF_MIXINGS> particle_mixing_names;
   extern const std::array<std::string, NUMBER_OF_INPUT_PARAMETERS> input_parameter_names;
   extern const std::array<std::string, NUMBER_OF_EXTRA_PARAMETERS> extra_parameter_names;
   extern const std::string model_name;
   constexpr bool is_low_energy_model = false;
   constexpr bool is_supersymmetric_model = false;
   constexpr bool is_FlexibleEFTHiggs = false;

   void print(std::ostream&);

   class HGTHDMIIMSSMBC_particle_names : public Names {
   public:
      virtual ~HGTHDMIIMSSMBC_particle_names() = default;
      virtual const std::string& get(int index) const override {
         return particle_names[index];
      }
      virtual int size() const override {
         return NUMBER_OF_PARTICLES;
      }
   };

   class HGTHDMIIMSSMBC_parameter_names : public Names {
   public:
      virtual ~HGTHDMIIMSSMBC_parameter_names() = default;
      virtual const std::string& get(int index) const override {
         return parameter_names[index];
      }
      virtual int size() const override {
         return NUMBER_OF_PARAMETERS;
      }
   };

   const HGTHDMIIMSSMBC_particle_names  particle_names_getter{};
   const HGTHDMIIMSSMBC_parameter_names parameter_names_getter{};

} // namespace HGTHDMIIMSSMBC_info

} // namespace flexiblesusy

#endif
