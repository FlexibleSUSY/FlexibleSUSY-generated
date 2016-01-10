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

// File generated at Sun 10 Jan 2016 15:32:45

#ifndef MRSSM_INFO_H
#define MRSSM_INFO_H

#include <iosfwd>

namespace flexiblesusy {

namespace MRSSM_info {
   enum Particles : unsigned {VG, Glu, Fv, SRdp, SRum, sigmaO, phiO, VP, VZ, Sd
      , Sv, Su, Se, hh, Ah, Rh, Hpm, Chi, Cha1, Cha2, Fe, Fd, Fu, VWm,
      NUMBER_OF_PARTICLES};

   enum Parameters : unsigned {Yd00, Yd01, Yd02, Yd10, Yd11, Yd12, Yd20, Yd21,
      Yd22, Ye00, Ye01, Ye02, Ye10, Ye11, Ye12, Ye20, Ye21, Ye22, LamTD, LamTU,
      LamSD, LamSU, Yu00, Yu01, Yu02, Yu10, Yu11, Yu12, Yu20, Yu21, Yu22, Mu, MuD,
      MuU, g1, g2, g3, vd, vu, vT, vS, BMu, BMuD, BMuU, mq200, mq201, mq202,
      mq210, mq211, mq212, mq220, mq221, mq222, ml200, ml201, ml202, ml210, ml211,
      ml212, ml220, ml221, ml222, mHd2, mHu2, md200, md201, md202, md210, md211,
      md212, md220, md221, md222, mu200, mu201, mu202, mu210, mu211, mu212, mu220,
      mu221, mu222, me200, me201, me202, me210, me211, me212, me220, me221, me222
      , mS2, mT2, moc2, mRd2, mRu2, MDBS, MDWBT, MDGoc, NUMBER_OF_PARAMETERS};

   enum Mixings : unsigned {ZD00, ZD01, ZD02, ZD03, ZD04, ZD05, ZD10, ZD11,
      ZD12, ZD13, ZD14, ZD15, ZD20, ZD21, ZD22, ZD23, ZD24, ZD25, ZD30, ZD31, ZD32
      , ZD33, ZD34, ZD35, ZD40, ZD41, ZD42, ZD43, ZD44, ZD45, ZD50, ZD51, ZD52,
      ZD53, ZD54, ZD55, ZV00, ZV01, ZV02, ZV10, ZV11, ZV12, ZV20, ZV21, ZV22, ZU00
      , ZU01, ZU02, ZU03, ZU04, ZU05, ZU10, ZU11, ZU12, ZU13, ZU14, ZU15, ZU20,
      ZU21, ZU22, ZU23, ZU24, ZU25, ZU30, ZU31, ZU32, ZU33, ZU34, ZU35, ZU40, ZU41
      , ZU42, ZU43, ZU44, ZU45, ZU50, ZU51, ZU52, ZU53, ZU54, ZU55, ZE00, ZE01,
      ZE02, ZE03, ZE04, ZE05, ZE10, ZE11, ZE12, ZE13, ZE14, ZE15, ZE20, ZE21, ZE22
      , ZE23, ZE24, ZE25, ZE30, ZE31, ZE32, ZE33, ZE34, ZE35, ZE40, ZE41, ZE42,
      ZE43, ZE44, ZE45, ZE50, ZE51, ZE52, ZE53, ZE54, ZE55, ZH00, ZH01, ZH02, ZH03
      , ZH10, ZH11, ZH12, ZH13, ZH20, ZH21, ZH22, ZH23, ZH30, ZH31, ZH32, ZH33,
      ZA00, ZA01, ZA02, ZA03, ZA10, ZA11, ZA12, ZA13, ZA20, ZA21, ZA22, ZA23, ZA30
      , ZA31, ZA32, ZA33, ZHR00, ZHR01, ZHR10, ZHR11, ZP00, ZP01, ZP02, ZP03, ZP10
      , ZP11, ZP12, ZP13, ZP20, ZP21, ZP22, ZP23, ZP30, ZP31, ZP32, ZP33, ReZN100,
      ImZN100, ReZN101, ImZN101, ReZN102, ImZN102, ReZN103, ImZN103, ReZN110,
      ImZN110, ReZN111, ImZN111, ReZN112, ImZN112, ReZN113, ImZN113, ReZN120,
      ImZN120, ReZN121, ImZN121, ReZN122, ImZN122, ReZN123, ImZN123, ReZN130,
      ImZN130, ReZN131, ImZN131, ReZN132, ImZN132, ReZN133, ImZN133, ReZN200,
      ImZN200, ReZN201, ImZN201, ReZN202, ImZN202, ReZN203, ImZN203, ReZN210,
      ImZN210, ReZN211, ImZN211, ReZN212, ImZN212, ReZN213, ImZN213, ReZN220,
      ImZN220, ReZN221, ImZN221, ReZN222, ImZN222, ReZN223, ImZN223, ReZN230,
      ImZN230, ReZN231, ImZN231, ReZN232, ImZN232, ReZN233, ImZN233, ReUM100,
      ImUM100, ReUM101, ImUM101, ReUM110, ImUM110, ReUM111, ImUM111, ReUP100,
      ImUP100, ReUP101, ImUP101, ReUP110, ImUP110, ReUP111, ImUP111, ReUM200,
      ImUM200, ReUM201, ImUM201, ReUM210, ImUM210, ReUM211, ImUM211, ReUP200,
      ImUP200, ReUP201, ImUP201, ReUP210, ImUP210, ReUP211, ImUP211, ReZEL00,
      ImZEL00, ReZEL01, ImZEL01, ReZEL02, ImZEL02, ReZEL10, ImZEL10, ReZEL11,
      ImZEL11, ReZEL12, ImZEL12, ReZEL20, ImZEL20, ReZEL21, ImZEL21, ReZEL22,
      ImZEL22, ReZER00, ImZER00, ReZER01, ImZER01, ReZER02, ImZER02, ReZER10,
      ImZER10, ReZER11, ImZER11, ReZER12, ImZER12, ReZER20, ImZER20, ReZER21,
      ImZER21, ReZER22, ImZER22, ReZDL00, ImZDL00, ReZDL01, ImZDL01, ReZDL02,
      ImZDL02, ReZDL10, ImZDL10, ReZDL11, ImZDL11, ReZDL12, ImZDL12, ReZDL20,
      ImZDL20, ReZDL21, ImZDL21, ReZDL22, ImZDL22, ReZDR00, ImZDR00, ReZDR01,
      ImZDR01, ReZDR02, ImZDR02, ReZDR10, ImZDR10, ReZDR11, ImZDR11, ReZDR12,
      ImZDR12, ReZDR20, ImZDR20, ReZDR21, ImZDR21, ReZDR22, ImZDR22, ReZUL00,
      ImZUL00, ReZUL01, ImZUL01, ReZUL02, ImZUL02, ReZUL10, ImZUL10, ReZUL11,
      ImZUL11, ReZUL12, ImZUL12, ReZUL20, ImZUL20, ReZUL21, ImZUL21, ReZUL22,
      ImZUL22, ReZUR00, ImZUR00, ReZUR01, ImZUR01, ReZUR02, ImZUR02, ReZUR10,
      ImZUR10, ReZUR11, ImZUR11, ReZUR12, ImZUR12, ReZUR20, ImZUR20, ReZUR21,
      ImZUR21, ReZUR22, ImZUR22, NUMBER_OF_MIXINGS};

   enum Input_parameters : unsigned {TanBeta, LamTDInput, LamTUInput,
      LamSDInput, LamSUInput, MuInput, MuDInput, MuUInput, vTInput, vSInput,
      BMuInput, BMuDInput, BMuUInput, mq2Input00, mq2Input01, mq2Input02,
      mq2Input10, mq2Input11, mq2Input12, mq2Input20, mq2Input21, mq2Input22,
      ml2Input00, ml2Input01, ml2Input02, ml2Input10, ml2Input11, ml2Input12,
      ml2Input20, ml2Input21, ml2Input22, md2Input00, md2Input01, md2Input02,
      md2Input10, md2Input11, md2Input12, md2Input20, md2Input21, md2Input22,
      mu2Input00, mu2Input01, mu2Input02, mu2Input10, mu2Input11, mu2Input12,
      mu2Input20, mu2Input21, mu2Input22, me2Input00, me2Input01, me2Input02,
      me2Input10, me2Input11, me2Input12, me2Input20, me2Input21, me2Input22,
      moc2Input, mRd2Input, mRu2Input, MDBSInput, MDWBTInput, MDGocInput,
      NUMBER_OF_INPUT_PARAMETERS};

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
