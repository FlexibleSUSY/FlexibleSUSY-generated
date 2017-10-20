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

// File generated at Fri 20 Oct 2017 09:00:44

#ifndef E6SSM_INFO_H
#define E6SSM_INFO_H

#include "problems.hpp"

#include <array>
#include <iosfwd>
#include <string>

namespace flexiblesusy {

namespace E6SSM_info {
   enum Particles : int { VG, Glu, Fv, ChaP, Sd, Sv, Su, Se, SDX, hh, Ah, Hpm,
      Chi, Cha, Fe, Fd, Fu, FDX, SHI0, SHIp, ChaI, ChiI, SSI0, FSI, SHp0, SHpp,
      ChiP, VWm, VP, VZ, VZp, NUMBER_OF_PARTICLES };

   enum Masses : int { MVG, MGlu, MFv_1, MFv_2, MFv_3, MChaP, MSd_1, MSd_2,
      MSd_3, MSd_4, MSd_5, MSd_6, MSv_1, MSv_2, MSv_3, MSu_1, MSu_2, MSu_3, MSu_4,
      MSu_5, MSu_6, MSe_1, MSe_2, MSe_3, MSe_4, MSe_5, MSe_6, MSDX_1, MSDX_2,
      MSDX_3, MSDX_4, MSDX_5, MSDX_6, Mhh_1, Mhh_2, Mhh_3, MAh_1, MAh_2, MAh_3,
      MHpm_1, MHpm_2, MChi_1, MChi_2, MChi_3, MChi_4, MChi_5, MChi_6, MCha_1,
      MCha_2, MFe_1, MFe_2, MFe_3, MFd_1, MFd_2, MFd_3, MFu_1, MFu_2, MFu_3,
      MFDX_1, MFDX_2, MFDX_3, MSHI0_1, MSHI0_2, MSHI0_3, MSHI0_4, MSHIp_1, MSHIp_2
      , MSHIp_3, MSHIp_4, MChaI_1, MChaI_2, MChiI_1, MChiI_2, MChiI_3, MChiI_4,
      MSSI0_1, MSSI0_2, MFSI_1, MFSI_2, MSHp0_1, MSHp0_2, MSHpp_1, MSHpp_2,
      MChiP_1, MChiP_2, MVWm, MVP, MVZ, MVZp, NUMBER_OF_MASSES };

   enum Parameters : int { Yd0_0, Yd0_1, Yd0_2, Yd1_0, Yd1_1, Yd1_2, Yd2_0,
      Yd2_1, Yd2_2, Ye0_0, Ye0_1, Ye0_2, Ye1_0, Ye1_1, Ye1_2, Ye2_0, Ye2_1, Ye2_2,
      Kappa0_0, Kappa0_1, Kappa0_2, Kappa1_0, Kappa1_1, Kappa1_2, Kappa2_0,
      Kappa2_1, Kappa2_2, Lambda120_0, Lambda120_1, Lambda121_0, Lambda121_1,
      Lambdax, Yu0_0, Yu0_1, Yu0_2, Yu1_0, Yu1_1, Yu1_2, Yu2_0, Yu2_1, Yu2_2, MuPr
      , g1, g2, g3, gN, vd, vu, vs, TYd0_0, TYd0_1, TYd0_2, TYd1_0, TYd1_1, TYd1_2
      , TYd2_0, TYd2_1, TYd2_2, TYe0_0, TYe0_1, TYe0_2, TYe1_0, TYe1_1, TYe1_2,
      TYe2_0, TYe2_1, TYe2_2, TKappa0_0, TKappa0_1, TKappa0_2, TKappa1_0,
      TKappa1_1, TKappa1_2, TKappa2_0, TKappa2_1, TKappa2_2, TLambda120_0,
      TLambda120_1, TLambda121_0, TLambda121_1, TLambdax, TYu0_0, TYu0_1, TYu0_2,
      TYu1_0, TYu1_1, TYu1_2, TYu2_0, TYu2_1, TYu2_2, BMuPr, mq20_0, mq20_1,
      mq20_2, mq21_0, mq21_1, mq21_2, mq22_0, mq22_1, mq22_2, ml20_0, ml20_1,
      ml20_2, ml21_0, ml21_1, ml21_2, ml22_0, ml22_1, ml22_2, mHd2, mHu2, md20_0,
      md20_1, md20_2, md21_0, md21_1, md21_2, md22_0, md22_1, md22_2, mu20_0,
      mu20_1, mu20_2, mu21_0, mu21_1, mu21_2, mu22_0, mu22_1, mu22_2, me20_0,
      me20_1, me20_2, me21_0, me21_1, me21_2, me22_0, me22_1, me22_2, ms2,
      mH1I20_0, mH1I20_1, mH1I21_0, mH1I21_1, mH2I20_0, mH2I20_1, mH2I21_0,
      mH2I21_1, msI20_0, msI20_1, msI21_0, msI21_1, mDx20_0, mDx20_1, mDx20_2,
      mDx21_0, mDx21_1, mDx21_2, mDx22_0, mDx22_1, mDx22_2, mDxbar20_0, mDxbar20_1
      , mDxbar20_2, mDxbar21_0, mDxbar21_1, mDxbar21_2, mDxbar22_0, mDxbar22_1,
      mDxbar22_2, mHp2, mHpbar2, MassB, MassWB, MassG, MassBp,
      NUMBER_OF_PARAMETERS };

   enum Mixings : int { ZD0_0, ZD0_1, ZD0_2, ZD0_3, ZD0_4, ZD0_5, ZD1_0, ZD1_1,
      ZD1_2, ZD1_3, ZD1_4, ZD1_5, ZD2_0, ZD2_1, ZD2_2, ZD2_3, ZD2_4, ZD2_5, ZD3_0
      , ZD3_1, ZD3_2, ZD3_3, ZD3_4, ZD3_5, ZD4_0, ZD4_1, ZD4_2, ZD4_3, ZD4_4,
      ZD4_5, ZD5_0, ZD5_1, ZD5_2, ZD5_3, ZD5_4, ZD5_5, ZV0_0, ZV0_1, ZV0_2, ZV1_0,
      ZV1_1, ZV1_2, ZV2_0, ZV2_1, ZV2_2, ZU0_0, ZU0_1, ZU0_2, ZU0_3, ZU0_4, ZU0_5
      , ZU1_0, ZU1_1, ZU1_2, ZU1_3, ZU1_4, ZU1_5, ZU2_0, ZU2_1, ZU2_2, ZU2_3,
      ZU2_4, ZU2_5, ZU3_0, ZU3_1, ZU3_2, ZU3_3, ZU3_4, ZU3_5, ZU4_0, ZU4_1, ZU4_2,
      ZU4_3, ZU4_4, ZU4_5, ZU5_0, ZU5_1, ZU5_2, ZU5_3, ZU5_4, ZU5_5, ZE0_0, ZE0_1
      , ZE0_2, ZE0_3, ZE0_4, ZE0_5, ZE1_0, ZE1_1, ZE1_2, ZE1_3, ZE1_4, ZE1_5,
      ZE2_0, ZE2_1, ZE2_2, ZE2_3, ZE2_4, ZE2_5, ZE3_0, ZE3_1, ZE3_2, ZE3_3, ZE3_4,
      ZE3_5, ZE4_0, ZE4_1, ZE4_2, ZE4_3, ZE4_4, ZE4_5, ZE5_0, ZE5_1, ZE5_2, ZE5_3
      , ZE5_4, ZE5_5, ZDX0_0, ZDX0_1, ZDX0_2, ZDX0_3, ZDX0_4, ZDX0_5, ZDX1_0,
      ZDX1_1, ZDX1_2, ZDX1_3, ZDX1_4, ZDX1_5, ZDX2_0, ZDX2_1, ZDX2_2, ZDX2_3,
      ZDX2_4, ZDX2_5, ZDX3_0, ZDX3_1, ZDX3_2, ZDX3_3, ZDX3_4, ZDX3_5, ZDX4_0,
      ZDX4_1, ZDX4_2, ZDX4_3, ZDX4_4, ZDX4_5, ZDX5_0, ZDX5_1, ZDX5_2, ZDX5_3,
      ZDX5_4, ZDX5_5, ZH0_0, ZH0_1, ZH0_2, ZH1_0, ZH1_1, ZH1_2, ZH2_0, ZH2_1,
      ZH2_2, ZA0_0, ZA0_1, ZA0_2, ZA1_0, ZA1_1, ZA1_2, ZA2_0, ZA2_1, ZA2_2, ZP0_0,
      ZP0_1, ZP1_0, ZP1_1, ReZN0_0, ImZN0_0, ReZN0_1, ImZN0_1, ReZN0_2, ImZN0_2,
      ReZN0_3, ImZN0_3, ReZN0_4, ImZN0_4, ReZN0_5, ImZN0_5, ReZN1_0, ImZN1_0,
      ReZN1_1, ImZN1_1, ReZN1_2, ImZN1_2, ReZN1_3, ImZN1_3, ReZN1_4, ImZN1_4,
      ReZN1_5, ImZN1_5, ReZN2_0, ImZN2_0, ReZN2_1, ImZN2_1, ReZN2_2, ImZN2_2,
      ReZN2_3, ImZN2_3, ReZN2_4, ImZN2_4, ReZN2_5, ImZN2_5, ReZN3_0, ImZN3_0,
      ReZN3_1, ImZN3_1, ReZN3_2, ImZN3_2, ReZN3_3, ImZN3_3, ReZN3_4, ImZN3_4,
      ReZN3_5, ImZN3_5, ReZN4_0, ImZN4_0, ReZN4_1, ImZN4_1, ReZN4_2, ImZN4_2,
      ReZN4_3, ImZN4_3, ReZN4_4, ImZN4_4, ReZN4_5, ImZN4_5, ReZN5_0, ImZN5_0,
      ReZN5_1, ImZN5_1, ReZN5_2, ImZN5_2, ReZN5_3, ImZN5_3, ReZN5_4, ImZN5_4,
      ReZN5_5, ImZN5_5, ReUM0_0, ImUM0_0, ReUM0_1, ImUM0_1, ReUM1_0, ImUM1_0,
      ReUM1_1, ImUM1_1, ReUP0_0, ImUP0_0, ReUP0_1, ImUP0_1, ReUP1_0, ImUP1_0,
      ReUP1_1, ImUP1_1, ReZEL0_0, ImZEL0_0, ReZEL0_1, ImZEL0_1, ReZEL0_2, ImZEL0_2
      , ReZEL1_0, ImZEL1_0, ReZEL1_1, ImZEL1_1, ReZEL1_2, ImZEL1_2, ReZEL2_0,
      ImZEL2_0, ReZEL2_1, ImZEL2_1, ReZEL2_2, ImZEL2_2, ReZER0_0, ImZER0_0,
      ReZER0_1, ImZER0_1, ReZER0_2, ImZER0_2, ReZER1_0, ImZER1_0, ReZER1_1,
      ImZER1_1, ReZER1_2, ImZER1_2, ReZER2_0, ImZER2_0, ReZER2_1, ImZER2_1,
      ReZER2_2, ImZER2_2, ReZDL0_0, ImZDL0_0, ReZDL0_1, ImZDL0_1, ReZDL0_2,
      ImZDL0_2, ReZDL1_0, ImZDL1_0, ReZDL1_1, ImZDL1_1, ReZDL1_2, ImZDL1_2,
      ReZDL2_0, ImZDL2_0, ReZDL2_1, ImZDL2_1, ReZDL2_2, ImZDL2_2, ReZDR0_0,
      ImZDR0_0, ReZDR0_1, ImZDR0_1, ReZDR0_2, ImZDR0_2, ReZDR1_0, ImZDR1_0,
      ReZDR1_1, ImZDR1_1, ReZDR1_2, ImZDR1_2, ReZDR2_0, ImZDR2_0, ReZDR2_1,
      ImZDR2_1, ReZDR2_2, ImZDR2_2, ReZUL0_0, ImZUL0_0, ReZUL0_1, ImZUL0_1,
      ReZUL0_2, ImZUL0_2, ReZUL1_0, ImZUL1_0, ReZUL1_1, ImZUL1_1, ReZUL1_2,
      ImZUL1_2, ReZUL2_0, ImZUL2_0, ReZUL2_1, ImZUL2_1, ReZUL2_2, ImZUL2_2,
      ReZUR0_0, ImZUR0_0, ReZUR0_1, ImZUR0_1, ReZUR0_2, ImZUR0_2, ReZUR1_0,
      ImZUR1_0, ReZUR1_1, ImZUR1_1, ReZUR1_2, ImZUR1_2, ReZUR2_0, ImZUR2_0,
      ReZUR2_1, ImZUR2_1, ReZUR2_2, ImZUR2_2, ReZDXL0_0, ImZDXL0_0, ReZDXL0_1,
      ImZDXL0_1, ReZDXL0_2, ImZDXL0_2, ReZDXL1_0, ImZDXL1_0, ReZDXL1_1, ImZDXL1_1,
      ReZDXL1_2, ImZDXL1_2, ReZDXL2_0, ImZDXL2_0, ReZDXL2_1, ImZDXL2_1, ReZDXL2_2
      , ImZDXL2_2, ReZDXR0_0, ImZDXR0_0, ReZDXR0_1, ImZDXR0_1, ReZDXR0_2,
      ImZDXR0_2, ReZDXR1_0, ImZDXR1_0, ReZDXR1_1, ImZDXR1_1, ReZDXR1_2, ImZDXR1_2,
      ReZDXR2_0, ImZDXR2_0, ReZDXR2_1, ImZDXR2_1, ReZDXR2_2, ImZDXR2_2, UHI00_0,
      UHI00_1, UHI00_2, UHI00_3, UHI01_0, UHI01_1, UHI01_2, UHI01_3, UHI02_0,
      UHI02_1, UHI02_2, UHI02_3, UHI03_0, UHI03_1, UHI03_2, UHI03_3, UHIp0_0,
      UHIp0_1, UHIp0_2, UHIp0_3, UHIp1_0, UHIp1_1, UHIp1_2, UHIp1_3, UHIp2_0,
      UHIp2_1, UHIp2_2, UHIp2_3, UHIp3_0, UHIp3_1, UHIp3_2, UHIp3_3, ReZMI0_0,
      ImZMI0_0, ReZMI0_1, ImZMI0_1, ReZMI1_0, ImZMI1_0, ReZMI1_1, ImZMI1_1,
      ReZPI0_0, ImZPI0_0, ReZPI0_1, ImZPI0_1, ReZPI1_0, ImZPI1_0, ReZPI1_1,
      ImZPI1_1, ReZNI0_0, ImZNI0_0, ReZNI0_1, ImZNI0_1, ReZNI0_2, ImZNI0_2,
      ReZNI0_3, ImZNI0_3, ReZNI1_0, ImZNI1_0, ReZNI1_1, ImZNI1_1, ReZNI1_2,
      ImZNI1_2, ReZNI1_3, ImZNI1_3, ReZNI2_0, ImZNI2_0, ReZNI2_1, ImZNI2_1,
      ReZNI2_2, ImZNI2_2, ReZNI2_3, ImZNI2_3, ReZNI3_0, ImZNI3_0, ReZNI3_1,
      ImZNI3_1, ReZNI3_2, ImZNI3_2, ReZNI3_3, ImZNI3_3, ZSSI0_0, ZSSI0_1, ZSSI1_0,
      ZSSI1_1, ReZFSI0_0, ImZFSI0_0, ReZFSI0_1, ImZFSI0_1, ReZFSI1_0, ImZFSI1_0,
      ReZFSI1_1, ImZFSI1_1, UHp00_0, UHp00_1, UHp01_0, UHp01_1, UHpp0_0, UHpp0_1,
      UHpp1_0, UHpp1_1, ReZNp0_0, ImZNp0_0, ReZNp0_1, ImZNp0_1, ReZNp1_0, ImZNp1_0
      , ReZNp1_1, ImZNp1_1, ZZ0_0, ZZ0_1, ZZ0_2, ZZ1_0, ZZ1_1, ZZ1_2, ZZ2_0, ZZ2_1
      , ZZ2_2, NUMBER_OF_MIXINGS };

   enum Input_parameters : int { m0, m12, TanBeta, Azero, LambdaInput,
      KappaInput, muPrimeInput, BmuPrimeInput, vSInput, Lambda12Input,
      NUMBER_OF_INPUT_PARAMETERS };

   enum Extra_parameters : int { NUMBER_OF_EXTRA_PARAMETERS };

   extern const double normalization_g1;
   extern const double normalization_g2;
   extern const double normalization_g3;
   extern const double normalization_gN;

   extern const std::array<int, NUMBER_OF_PARTICLES> particle_multiplicities;
   extern const std::array<std::string, NUMBER_OF_PARTICLES> particle_names;
   extern const std::array<std::string, NUMBER_OF_PARTICLES> particle_latex_names;
   extern const std::array<std::string, NUMBER_OF_PARAMETERS> parameter_names;
   extern const std::array<std::string, NUMBER_OF_MIXINGS> particle_mixing_names;
   extern const std::array<std::string, NUMBER_OF_INPUT_PARAMETERS> input_parameter_names;
   extern const std::array<std::string, NUMBER_OF_EXTRA_PARAMETERS> extra_parameter_names;
   extern const std::string model_name;
   constexpr bool is_low_energy_model = false;
   constexpr bool is_supersymmetric_model = true;
   constexpr bool is_FlexibleEFTHiggs = false;

   void print(std::ostream&);

   class E6SSM_particle_names : public Names {
   public:
      virtual ~E6SSM_particle_names() = default;
      virtual const std::string& get(int index) const override {
         return particle_names[index];
      }
      virtual int size() const override {
         return NUMBER_OF_PARTICLES;
      }
   };

   class E6SSM_parameter_names : public Names {
   public:
      virtual ~E6SSM_parameter_names() = default;
      virtual const std::string& get(int index) const override {
         return parameter_names[index];
      }
      virtual int size() const override {
         return NUMBER_OF_PARAMETERS;
      }
   };

   const E6SSM_particle_names  particle_names_getter{};
   const E6SSM_parameter_names parameter_names_getter{};

} // namespace E6SSM_info

} // namespace flexiblesusy

#endif
