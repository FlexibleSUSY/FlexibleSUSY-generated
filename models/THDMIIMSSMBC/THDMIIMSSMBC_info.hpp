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


#ifndef THDMIIMSSMBC_INFO_H
#define THDMIIMSSMBC_INFO_H

#include "names.hpp"

#include <array>
#include <iosfwd>
#include <optional>
#include <string>
#include <utility>

namespace flexiblesusy {

namespace THDMIIMSSMBC_info {
   enum Particles : int { VG = 0, Fv, hh, Ah, Hm, Fd, Fu, Fe, VWm, VP, VZ,
      NUMBER_OF_PARTICLES };

   enum Masses : int { MVG = 0, MFv_1, MFv_2, MFv_3, Mhh_1, Mhh_2, MAh_1, MAh_2,
      MHm_1, MHm_2, MFd_1, MFd_2, MFd_3, MFu_1, MFu_2, MFu_3, MFe_1, MFe_2, MFe_3,
      MVWm, MVP, MVZ, NUMBER_OF_MASSES };

   enum Parameters : int { g1, g2, g3, Lambda6, Lambda5, Lambda7, Lambda1, Lambda4
      , Lambda3, Lambda2, Yu0_0, Yu0_1, Yu0_2, Yu1_0, Yu1_1, Yu1_2, Yu2_0, Yu2_1,
      Yu2_2, Yd0_0, Yd0_1, Yd0_2, Yd1_0, Yd1_1, Yd1_2, Yd2_0, Yd2_1, Yd2_2, Ye0_0,
      Ye0_1, Ye0_2, Ye1_0, Ye1_1, Ye1_2, Ye2_0, Ye2_1, Ye2_2, M122, M112, M222, v1
      , v2, NUMBER_OF_PARAMETERS };

   enum Mixings : int { ZH0_0 = 0, ZH0_1, ZH1_0, ZH1_1, ZA0_0, ZA0_1, ZA1_0, ZA1_1
      , ZP0_0, ZP0_1, ZP1_0, ZP1_1, ReVd0_0, ImVd0_0, ReVd0_1, ImVd0_1, ReVd0_2,
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
      ImUe1_2, ReUe2_0, ImUe2_0, ReUe2_1, ImUe2_1, ReUe2_2, ImUe2_2, ZZ0_0, ZZ0_1,
      ZZ1_0, ZZ1_1, NUMBER_OF_MIXINGS };

   enum Input_parameters : int { TanBeta, MSUSY, MEWSB, MuInput, M1Input, M2Input,
      MAInput, LambdaLoopOrder, AeInput0_0, AeInput0_1, AeInput0_2, AeInput1_0,
      AeInput1_1, AeInput1_2, AeInput2_0, AeInput2_1, AeInput2_2, AdInput0_0,
      AdInput0_1, AdInput0_2, AdInput1_0, AdInput1_1, AdInput1_2, AdInput2_0,
      AdInput2_1, AdInput2_2, AuInput0_0, AuInput0_1, AuInput0_2, AuInput1_0,
      AuInput1_1, AuInput1_2, AuInput2_0, AuInput2_1, AuInput2_2, mslInput0,
      mslInput1, mslInput2, mseInput0, mseInput1, mseInput2, msqInput0, msqInput1,
      msqInput2, msdInput0, msdInput1, msdInput2, msuInput0, msuInput1, msuInput2,
      NUMBER_OF_INPUT_PARAMETERS };

   enum Extra_parameters : int { NUMBER_OF_EXTRA_PARAMETERS };

   extern const double normalization_g1;
   extern const double normalization_g2;
   extern const double normalization_g3;

   static constexpr int number_of_neutral_goldstones = 1;
   static constexpr int number_of_charged_goldstones = 1;

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
   constexpr bool is_CP_violating_Higgs_sector {false};

   int get_pdg_code_for_particle(Particles);
   int get_pdg_code_for_particle(Particles, int);
   std::string get_particle_name_from_pdg(int);
   std::pair<std::string, std::optional<unsigned int>> get_multiplet_and_index_from_pdg(int);
   void print(std::ostream&);

   class THDMIIMSSMBC_particle_names : public Names {
   public:
      virtual ~THDMIIMSSMBC_particle_names() = default;
      virtual const std::string& get(int index) const override {
         return particle_names[index];
      }
      virtual int size() const override {
         return NUMBER_OF_PARTICLES;
      }
   };

   class THDMIIMSSMBC_parameter_names : public Names {
   public:
      virtual ~THDMIIMSSMBC_parameter_names() = default;
      virtual const std::string& get(int index) const override {
         return parameter_names[index];
      }
      virtual int size() const override {
         return NUMBER_OF_PARAMETERS;
      }
   };

   const THDMIIMSSMBC_particle_names  particle_names_getter{};
   const THDMIIMSSMBC_parameter_names parameter_names_getter{};

} // namespace THDMIIMSSMBC_info

} // namespace flexiblesusy

#endif
