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

// File generated at Sun 10 Jan 2016 15:36:37

#ifndef E6SSM_INFO_H
#define E6SSM_INFO_H

#include <iosfwd>

namespace flexiblesusy {

namespace E6SSM_info {
   enum Particles : unsigned {VG, Glu, Fv, ChaP, VP, VZ, VZp, Sd, Sv, Su, Se,
      SDX, hh, Ah, Hpm, Chi, Cha, Fe, Fd, Fu, FDX, SHI0, SHIp, ChaI, ChiI, SSI0,
      FSI, SHp0, SHpp, ChiP, VWm, NUMBER_OF_PARTICLES};

   enum Parameters : unsigned {Yd00, Yd01, Yd02, Yd10, Yd11, Yd12, Yd20, Yd21,
      Yd22, Ye00, Ye01, Ye02, Ye10, Ye11, Ye12, Ye20, Ye21, Ye22, Kappa00, Kappa01
      , Kappa02, Kappa10, Kappa11, Kappa12, Kappa20, Kappa21, Kappa22, Lambda1200,
      Lambda1201, Lambda1210, Lambda1211, Lambdax, Yu00, Yu01, Yu02, Yu10, Yu11,
      Yu12, Yu20, Yu21, Yu22, MuPr, g1, g2, g3, gN, vd, vu, vs, TYd00, TYd01,
      TYd02, TYd10, TYd11, TYd12, TYd20, TYd21, TYd22, TYe00, TYe01, TYe02, TYe10,
      TYe11, TYe12, TYe20, TYe21, TYe22, TKappa00, TKappa01, TKappa02, TKappa10,
      TKappa11, TKappa12, TKappa20, TKappa21, TKappa22, TLambda1200, TLambda1201,
      TLambda1210, TLambda1211, TLambdax, TYu00, TYu01, TYu02, TYu10, TYu11, TYu12
      , TYu20, TYu21, TYu22, BMuPr, mq200, mq201, mq202, mq210, mq211, mq212,
      mq220, mq221, mq222, ml200, ml201, ml202, ml210, ml211, ml212, ml220, ml221,
      ml222, mHd2, mHu2, md200, md201, md202, md210, md211, md212, md220, md221,
      md222, mu200, mu201, mu202, mu210, mu211, mu212, mu220, mu221, mu222, me200,
      me201, me202, me210, me211, me212, me220, me221, me222, ms2, mH1I200,
      mH1I201, mH1I210, mH1I211, mH2I200, mH2I201, mH2I210, mH2I211, msI200,
      msI201, msI210, msI211, mDx200, mDx201, mDx202, mDx210, mDx211, mDx212,
      mDx220, mDx221, mDx222, mDxbar200, mDxbar201, mDxbar202, mDxbar210,
      mDxbar211, mDxbar212, mDxbar220, mDxbar221, mDxbar222, mHp2, mHpbar2, MassB,
      MassWB, MassG, MassBp, NUMBER_OF_PARAMETERS};

   enum Mixings : unsigned {ZD00, ZD01, ZD02, ZD03, ZD04, ZD05, ZD10, ZD11,
      ZD12, ZD13, ZD14, ZD15, ZD20, ZD21, ZD22, ZD23, ZD24, ZD25, ZD30, ZD31, ZD32
      , ZD33, ZD34, ZD35, ZD40, ZD41, ZD42, ZD43, ZD44, ZD45, ZD50, ZD51, ZD52,
      ZD53, ZD54, ZD55, ZV00, ZV01, ZV02, ZV10, ZV11, ZV12, ZV20, ZV21, ZV22, ZU00
      , ZU01, ZU02, ZU03, ZU04, ZU05, ZU10, ZU11, ZU12, ZU13, ZU14, ZU15, ZU20,
      ZU21, ZU22, ZU23, ZU24, ZU25, ZU30, ZU31, ZU32, ZU33, ZU34, ZU35, ZU40, ZU41
      , ZU42, ZU43, ZU44, ZU45, ZU50, ZU51, ZU52, ZU53, ZU54, ZU55, ZE00, ZE01,
      ZE02, ZE03, ZE04, ZE05, ZE10, ZE11, ZE12, ZE13, ZE14, ZE15, ZE20, ZE21, ZE22
      , ZE23, ZE24, ZE25, ZE30, ZE31, ZE32, ZE33, ZE34, ZE35, ZE40, ZE41, ZE42,
      ZE43, ZE44, ZE45, ZE50, ZE51, ZE52, ZE53, ZE54, ZE55, ZDX00, ZDX01, ZDX02,
      ZDX03, ZDX04, ZDX05, ZDX10, ZDX11, ZDX12, ZDX13, ZDX14, ZDX15, ZDX20, ZDX21,
      ZDX22, ZDX23, ZDX24, ZDX25, ZDX30, ZDX31, ZDX32, ZDX33, ZDX34, ZDX35, ZDX40
      , ZDX41, ZDX42, ZDX43, ZDX44, ZDX45, ZDX50, ZDX51, ZDX52, ZDX53, ZDX54,
      ZDX55, ZH00, ZH01, ZH02, ZH10, ZH11, ZH12, ZH20, ZH21, ZH22, ZA00, ZA01,
      ZA02, ZA10, ZA11, ZA12, ZA20, ZA21, ZA22, ZP00, ZP01, ZP10, ZP11, ReZN00,
      ImZN00, ReZN01, ImZN01, ReZN02, ImZN02, ReZN03, ImZN03, ReZN04, ImZN04,
      ReZN05, ImZN05, ReZN10, ImZN10, ReZN11, ImZN11, ReZN12, ImZN12, ReZN13,
      ImZN13, ReZN14, ImZN14, ReZN15, ImZN15, ReZN20, ImZN20, ReZN21, ImZN21,
      ReZN22, ImZN22, ReZN23, ImZN23, ReZN24, ImZN24, ReZN25, ImZN25, ReZN30,
      ImZN30, ReZN31, ImZN31, ReZN32, ImZN32, ReZN33, ImZN33, ReZN34, ImZN34,
      ReZN35, ImZN35, ReZN40, ImZN40, ReZN41, ImZN41, ReZN42, ImZN42, ReZN43,
      ImZN43, ReZN44, ImZN44, ReZN45, ImZN45, ReZN50, ImZN50, ReZN51, ImZN51,
      ReZN52, ImZN52, ReZN53, ImZN53, ReZN54, ImZN54, ReZN55, ImZN55, ReUM00,
      ImUM00, ReUM01, ImUM01, ReUM10, ImUM10, ReUM11, ImUM11, ReUP00, ImUP00,
      ReUP01, ImUP01, ReUP10, ImUP10, ReUP11, ImUP11, ReZEL00, ImZEL00, ReZEL01,
      ImZEL01, ReZEL02, ImZEL02, ReZEL10, ImZEL10, ReZEL11, ImZEL11, ReZEL12,
      ImZEL12, ReZEL20, ImZEL20, ReZEL21, ImZEL21, ReZEL22, ImZEL22, ReZER00,
      ImZER00, ReZER01, ImZER01, ReZER02, ImZER02, ReZER10, ImZER10, ReZER11,
      ImZER11, ReZER12, ImZER12, ReZER20, ImZER20, ReZER21, ImZER21, ReZER22,
      ImZER22, ReZDL00, ImZDL00, ReZDL01, ImZDL01, ReZDL02, ImZDL02, ReZDL10,
      ImZDL10, ReZDL11, ImZDL11, ReZDL12, ImZDL12, ReZDL20, ImZDL20, ReZDL21,
      ImZDL21, ReZDL22, ImZDL22, ReZDR00, ImZDR00, ReZDR01, ImZDR01, ReZDR02,
      ImZDR02, ReZDR10, ImZDR10, ReZDR11, ImZDR11, ReZDR12, ImZDR12, ReZDR20,
      ImZDR20, ReZDR21, ImZDR21, ReZDR22, ImZDR22, ReZUL00, ImZUL00, ReZUL01,
      ImZUL01, ReZUL02, ImZUL02, ReZUL10, ImZUL10, ReZUL11, ImZUL11, ReZUL12,
      ImZUL12, ReZUL20, ImZUL20, ReZUL21, ImZUL21, ReZUL22, ImZUL22, ReZUR00,
      ImZUR00, ReZUR01, ImZUR01, ReZUR02, ImZUR02, ReZUR10, ImZUR10, ReZUR11,
      ImZUR11, ReZUR12, ImZUR12, ReZUR20, ImZUR20, ReZUR21, ImZUR21, ReZUR22,
      ImZUR22, ReZDXL00, ImZDXL00, ReZDXL01, ImZDXL01, ReZDXL02, ImZDXL02,
      ReZDXL10, ImZDXL10, ReZDXL11, ImZDXL11, ReZDXL12, ImZDXL12, ReZDXL20,
      ImZDXL20, ReZDXL21, ImZDXL21, ReZDXL22, ImZDXL22, ReZDXR00, ImZDXR00,
      ReZDXR01, ImZDXR01, ReZDXR02, ImZDXR02, ReZDXR10, ImZDXR10, ReZDXR11,
      ImZDXR11, ReZDXR12, ImZDXR12, ReZDXR20, ImZDXR20, ReZDXR21, ImZDXR21,
      ReZDXR22, ImZDXR22, UHI000, UHI001, UHI002, UHI003, UHI010, UHI011, UHI012,
      UHI013, UHI020, UHI021, UHI022, UHI023, UHI030, UHI031, UHI032, UHI033,
      UHIp00, UHIp01, UHIp02, UHIp03, UHIp10, UHIp11, UHIp12, UHIp13, UHIp20,
      UHIp21, UHIp22, UHIp23, UHIp30, UHIp31, UHIp32, UHIp33, ReZMI00, ImZMI00,
      ReZMI01, ImZMI01, ReZMI10, ImZMI10, ReZMI11, ImZMI11, ReZPI00, ImZPI00,
      ReZPI01, ImZPI01, ReZPI10, ImZPI10, ReZPI11, ImZPI11, ReZNI00, ImZNI00,
      ReZNI01, ImZNI01, ReZNI02, ImZNI02, ReZNI03, ImZNI03, ReZNI10, ImZNI10,
      ReZNI11, ImZNI11, ReZNI12, ImZNI12, ReZNI13, ImZNI13, ReZNI20, ImZNI20,
      ReZNI21, ImZNI21, ReZNI22, ImZNI22, ReZNI23, ImZNI23, ReZNI30, ImZNI30,
      ReZNI31, ImZNI31, ReZNI32, ImZNI32, ReZNI33, ImZNI33, ZSSI00, ZSSI01, ZSSI10
      , ZSSI11, ReZFSI00, ImZFSI00, ReZFSI01, ImZFSI01, ReZFSI10, ImZFSI10,
      ReZFSI11, ImZFSI11, UHp000, UHp001, UHp010, UHp011, UHpp00, UHpp01, UHpp10,
      UHpp11, ReZNp00, ImZNp00, ReZNp01, ImZNp01, ReZNp10, ImZNp10, ReZNp11,
      ImZNp11, NUMBER_OF_MIXINGS};

   enum Input_parameters : unsigned {m0, m12, TanBeta, Azero, LambdaInput,
      KappaInput, muPrimeInput, BmuPrimeInput, vSInput, Lambda12Input,
      NUMBER_OF_INPUT_PARAMETERS};

   extern const double normalization_g1;
   extern const double normalization_g2;
   extern const double normalization_g3;
   extern const double normalization_gN;

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
