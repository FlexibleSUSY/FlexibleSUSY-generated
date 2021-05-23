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


#ifndef MRSSM2_INFO_H
#define MRSSM2_INFO_H

#include "names.hpp"

#include <array>
#include <iosfwd>
#include <string>

namespace flexiblesusy {

namespace MRSSM2_info {
   enum Particles : int { VG = 0, Glu, Fv, SRdp, SRum, sigmaO, phiO, Sd, Sv, Su,
      Se, hh, Ah, Rh, Hpm, Chi, Cha1, Cha2, Fe, Fd, Fu, VWm, VP, VZ,
      NUMBER_OF_PARTICLES };

   enum Masses : int { MVG = 0, MGlu, MFv_1, MFv_2, MFv_3, MSRdp, MSRum, MsigmaO,
      MphiO, MSd_1, MSd_2, MSd_3, MSd_4, MSd_5, MSd_6, MSv_1, MSv_2, MSv_3, MSu_1,
      MSu_2, MSu_3, MSu_4, MSu_5, MSu_6, MSe_1, MSe_2, MSe_3, MSe_4, MSe_5, MSe_6,
      Mhh_1, Mhh_2, Mhh_3, Mhh_4, MAh_1, MAh_2, MAh_3, MAh_4, MRh_1, MRh_2, MHpm_1
      , MHpm_2, MHpm_3, MHpm_4, MChi_1, MChi_2, MChi_3, MChi_4, MCha1_1, MCha1_2,
      MCha2_1, MCha2_2, MFe_1, MFe_2, MFe_3, MFd_1, MFd_2, MFd_3, MFu_1, MFu_2,
      MFu_3, MVWm, MVP, MVZ, NUMBER_OF_MASSES };

   enum Parameters : int { Yd0_0, Yd0_1, Yd0_2, Yd1_0, Yd1_1, Yd1_2, Yd2_0, Yd2_1,
      Yd2_2, Ye0_0, Ye0_1, Ye0_2, Ye1_0, Ye1_1, Ye1_2, Ye2_0, Ye2_1, Ye2_2, LamTD,
      LamTU, LamSD, LamSU, Yu0_0, Yu0_1, Yu0_2, Yu1_0, Yu1_1, Yu1_2, Yu2_0, Yu2_1,
      Yu2_2, Mu, MuD, MuU, g1, g2, g3, vd, vu, vT, vS, BMu, BMuD, BMuU, mq20_0,
      mq20_1, mq20_2, mq21_0, mq21_1, mq21_2, mq22_0, mq22_1, mq22_2, ml20_0,
      ml20_1, ml20_2, ml21_0, ml21_1, ml21_2, ml22_0, ml22_1, ml22_2, mHd2, mHu2,
      md20_0, md20_1, md20_2, md21_0, md21_1, md21_2, md22_0, md22_1, md22_2,
      mu20_0, mu20_1, mu20_2, mu21_0, mu21_1, mu21_2, mu22_0, mu22_1, mu22_2,
      me20_0, me20_1, me20_2, me21_0, me21_1, me21_2, me22_0, me22_1, me22_2, mS2,
      mT2, moc2, mRd2, mRu2, MDBS, MDWBT, MDGoc, NUMBER_OF_PARAMETERS };

   enum Mixings : int { ZD0_0 = 0, ZD0_1, ZD0_2, ZD0_3, ZD0_4, ZD0_5, ZD1_0, ZD1_1
      , ZD1_2, ZD1_3, ZD1_4, ZD1_5, ZD2_0, ZD2_1, ZD2_2, ZD2_3, ZD2_4, ZD2_5,
      ZD3_0, ZD3_1, ZD3_2, ZD3_3, ZD3_4, ZD3_5, ZD4_0, ZD4_1, ZD4_2, ZD4_3, ZD4_4,
      ZD4_5, ZD5_0, ZD5_1, ZD5_2, ZD5_3, ZD5_4, ZD5_5, ZV0_0, ZV0_1, ZV0_2, ZV1_0,
      ZV1_1, ZV1_2, ZV2_0, ZV2_1, ZV2_2, ZU0_0, ZU0_1, ZU0_2, ZU0_3, ZU0_4, ZU0_5,
      ZU1_0, ZU1_1, ZU1_2, ZU1_3, ZU1_4, ZU1_5, ZU2_0, ZU2_1, ZU2_2, ZU2_3, ZU2_4,
      ZU2_5, ZU3_0, ZU3_1, ZU3_2, ZU3_3, ZU3_4, ZU3_5, ZU4_0, ZU4_1, ZU4_2, ZU4_3,
      ZU4_4, ZU4_5, ZU5_0, ZU5_1, ZU5_2, ZU5_3, ZU5_4, ZU5_5, ZE0_0, ZE0_1, ZE0_2,
      ZE0_3, ZE0_4, ZE0_5, ZE1_0, ZE1_1, ZE1_2, ZE1_3, ZE1_4, ZE1_5, ZE2_0, ZE2_1,
      ZE2_2, ZE2_3, ZE2_4, ZE2_5, ZE3_0, ZE3_1, ZE3_2, ZE3_3, ZE3_4, ZE3_5, ZE4_0,
      ZE4_1, ZE4_2, ZE4_3, ZE4_4, ZE4_5, ZE5_0, ZE5_1, ZE5_2, ZE5_3, ZE5_4, ZE5_5,
      ZH0_0, ZH0_1, ZH0_2, ZH0_3, ZH1_0, ZH1_1, ZH1_2, ZH1_3, ZH2_0, ZH2_1, ZH2_2,
      ZH2_3, ZH3_0, ZH3_1, ZH3_2, ZH3_3, ZA0_0, ZA0_1, ZA0_2, ZA0_3, ZA1_0, ZA1_1,
      ZA1_2, ZA1_3, ZA2_0, ZA2_1, ZA2_2, ZA2_3, ZA3_0, ZA3_1, ZA3_2, ZA3_3, ZHR0_0
      , ZHR0_1, ZHR1_0, ZHR1_1, ZP0_0, ZP0_1, ZP0_2, ZP0_3, ZP1_0, ZP1_1, ZP1_2,
      ZP1_3, ZP2_0, ZP2_1, ZP2_2, ZP2_3, ZP3_0, ZP3_1, ZP3_2, ZP3_3, ReZN10_0,
      ImZN10_0, ReZN10_1, ImZN10_1, ReZN10_2, ImZN10_2, ReZN10_3, ImZN10_3,
      ReZN11_0, ImZN11_0, ReZN11_1, ImZN11_1, ReZN11_2, ImZN11_2, ReZN11_3,
      ImZN11_3, ReZN12_0, ImZN12_0, ReZN12_1, ImZN12_1, ReZN12_2, ImZN12_2,
      ReZN12_3, ImZN12_3, ReZN13_0, ImZN13_0, ReZN13_1, ImZN13_1, ReZN13_2,
      ImZN13_2, ReZN13_3, ImZN13_3, ReZN20_0, ImZN20_0, ReZN20_1, ImZN20_1,
      ReZN20_2, ImZN20_2, ReZN20_3, ImZN20_3, ReZN21_0, ImZN21_0, ReZN21_1,
      ImZN21_1, ReZN21_2, ImZN21_2, ReZN21_3, ImZN21_3, ReZN22_0, ImZN22_0,
      ReZN22_1, ImZN22_1, ReZN22_2, ImZN22_2, ReZN22_3, ImZN22_3, ReZN23_0,
      ImZN23_0, ReZN23_1, ImZN23_1, ReZN23_2, ImZN23_2, ReZN23_3, ImZN23_3,
      ReUM10_0, ImUM10_0, ReUM10_1, ImUM10_1, ReUM11_0, ImUM11_0, ReUM11_1,
      ImUM11_1, ReUP10_0, ImUP10_0, ReUP10_1, ImUP10_1, ReUP11_0, ImUP11_0,
      ReUP11_1, ImUP11_1, ReUM20_0, ImUM20_0, ReUM20_1, ImUM20_1, ReUM21_0,
      ImUM21_0, ReUM21_1, ImUM21_1, ReUP20_0, ImUP20_0, ReUP20_1, ImUP20_1,
      ReUP21_0, ImUP21_0, ReUP21_1, ImUP21_1, ReZEL0_0, ImZEL0_0, ReZEL0_1,
      ImZEL0_1, ReZEL0_2, ImZEL0_2, ReZEL1_0, ImZEL1_0, ReZEL1_1, ImZEL1_1,
      ReZEL1_2, ImZEL1_2, ReZEL2_0, ImZEL2_0, ReZEL2_1, ImZEL2_1, ReZEL2_2,
      ImZEL2_2, ReZER0_0, ImZER0_0, ReZER0_1, ImZER0_1, ReZER0_2, ImZER0_2,
      ReZER1_0, ImZER1_0, ReZER1_1, ImZER1_1, ReZER1_2, ImZER1_2, ReZER2_0,
      ImZER2_0, ReZER2_1, ImZER2_1, ReZER2_2, ImZER2_2, ReZDL0_0, ImZDL0_0,
      ReZDL0_1, ImZDL0_1, ReZDL0_2, ImZDL0_2, ReZDL1_0, ImZDL1_0, ReZDL1_1,
      ImZDL1_1, ReZDL1_2, ImZDL1_2, ReZDL2_0, ImZDL2_0, ReZDL2_1, ImZDL2_1,
      ReZDL2_2, ImZDL2_2, ReZDR0_0, ImZDR0_0, ReZDR0_1, ImZDR0_1, ReZDR0_2,
      ImZDR0_2, ReZDR1_0, ImZDR1_0, ReZDR1_1, ImZDR1_1, ReZDR1_2, ImZDR1_2,
      ReZDR2_0, ImZDR2_0, ReZDR2_1, ImZDR2_1, ReZDR2_2, ImZDR2_2, ReZUL0_0,
      ImZUL0_0, ReZUL0_1, ImZUL0_1, ReZUL0_2, ImZUL0_2, ReZUL1_0, ImZUL1_0,
      ReZUL1_1, ImZUL1_1, ReZUL1_2, ImZUL1_2, ReZUL2_0, ImZUL2_0, ReZUL2_1,
      ImZUL2_1, ReZUL2_2, ImZUL2_2, ReZUR0_0, ImZUR0_0, ReZUR0_1, ImZUR0_1,
      ReZUR0_2, ImZUR0_2, ReZUR1_0, ImZUR1_0, ReZUR1_1, ImZUR1_1, ReZUR1_2,
      ImZUR1_2, ReZUR2_0, ImZUR2_0, ReZUR2_1, ImZUR2_1, ReZUR2_2, ImZUR2_2, ZZ0_0,
      ZZ0_1, ZZ1_0, ZZ1_1, NUMBER_OF_MIXINGS };

   enum Input_parameters : int { TanBeta, Ms, LamTDInput, LamTUInput, LamSDInput,
      LamSUInput, MuDInput, MuUInput, BMuInput, mq2Input0_0, mq2Input0_1,
      mq2Input0_2, mq2Input1_0, mq2Input1_1, mq2Input1_2, mq2Input2_0, mq2Input2_1
      , mq2Input2_2, ml2Input0_0, ml2Input0_1, ml2Input0_2, ml2Input1_0,
      ml2Input1_1, ml2Input1_2, ml2Input2_0, ml2Input2_1, ml2Input2_2, md2Input0_0
      , md2Input0_1, md2Input0_2, md2Input1_0, md2Input1_1, md2Input1_2,
      md2Input2_0, md2Input2_1, md2Input2_2, mu2Input0_0, mu2Input0_1, mu2Input0_2
      , mu2Input1_0, mu2Input1_1, mu2Input1_2, mu2Input2_0, mu2Input2_1,
      mu2Input2_2, me2Input0_0, me2Input0_1, me2Input0_2, me2Input1_0, me2Input1_1
      , me2Input1_2, me2Input2_0, me2Input2_1, me2Input2_2, mS2Input, mT2Input,
      moc2Input, mRd2Input, mRu2Input, MDBSInput, MDWBTInput, MDGocInput,
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
   constexpr bool is_low_energy_model = true;
   constexpr bool is_supersymmetric_model = true;
   constexpr bool is_FlexibleEFTHiggs = false;
   constexpr bool is_CP_violating_Higgs_sector {false};

   int get_pdg_code_for_particle(Particles);
   int get_pdg_code_for_particle(Particles, int);
   std::string get_particle_name_from_pdg(int);
   void print(std::ostream&);

   class MRSSM2_particle_names : public Names {
   public:
      virtual ~MRSSM2_particle_names() = default;
      virtual const std::string& get(int index) const override {
         return particle_names[index];
      }
      virtual int size() const override {
         return NUMBER_OF_PARTICLES;
      }
   };

   class MRSSM2_parameter_names : public Names {
   public:
      virtual ~MRSSM2_parameter_names() = default;
      virtual const std::string& get(int index) const override {
         return parameter_names[index];
      }
      virtual int size() const override {
         return NUMBER_OF_PARAMETERS;
      }
   };

   const MRSSM2_particle_names  particle_names_getter{};
   const MRSSM2_parameter_names parameter_names_getter{};

} // namespace MRSSM2_info

} // namespace flexiblesusy

#endif
