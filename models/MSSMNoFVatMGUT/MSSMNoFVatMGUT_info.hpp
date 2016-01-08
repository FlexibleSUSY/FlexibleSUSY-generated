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

// File generated at Fri 8 Jan 2016 15:26:01

#ifndef MSSMNoFVatMGUT_INFO_H
#define MSSMNoFVatMGUT_INFO_H

#include <iosfwd>

namespace flexiblesusy {

namespace MSSMNoFVatMGUT_info {
   enum Particles : unsigned {VG, Glu, VP, VZ, Fd, Fs, Fb, Fu, Fc, Ft, Fve, Fvm
      , Fvt, Fe, Fm, Ftau, SveL, SvmL, SvtL, Sd, Su, Se, Sm, Stau, Ss, Sc, Sb, St,
      hh, Ah, Hpm, Chi, Cha, VWm, NUMBER_OF_PARTICLES};

   enum Parameters : unsigned {Yd00, Yd01, Yd02, Yd10, Yd11, Yd12, Yd20, Yd21,
      Yd22, Ye00, Ye01, Ye02, Ye10, Ye11, Ye12, Ye20, Ye21, Ye22, Yu00, Yu01, Yu02
      , Yu10, Yu11, Yu12, Yu20, Yu21, Yu22, Mu, g1, g2, g3, vd, vu, TYd00, TYd01,
      TYd02, TYd10, TYd11, TYd12, TYd20, TYd21, TYd22, TYe00, TYe01, TYe02, TYe10,
      TYe11, TYe12, TYe20, TYe21, TYe22, TYu00, TYu01, TYu02, TYu10, TYu11, TYu12
      , TYu20, TYu21, TYu22, BMu, mq200, mq201, mq202, mq210, mq211, mq212, mq220,
      mq221, mq222, ml200, ml201, ml202, ml210, ml211, ml212, ml220, ml221, ml222
      , mHd2, mHu2, md200, md201, md202, md210, md211, md212, md220, md221, md222,
      mu200, mu201, mu202, mu210, mu211, mu212, mu220, mu221, mu222, me200, me201
      , me202, me210, me211, me212, me220, me221, me222, MassB, MassWB, MassG,
      NUMBER_OF_PARAMETERS};

   enum Mixings : unsigned {ZD00, ZD01, ZD10, ZD11, ZU00, ZU01, ZU10, ZU11,
      ZE00, ZE01, ZE10, ZE11, ZM00, ZM01, ZM10, ZM11, ZTau00, ZTau01, ZTau10,
      ZTau11, ZS00, ZS01, ZS10, ZS11, ZC00, ZC01, ZC10, ZC11, ZB00, ZB01, ZB10,
      ZB11, ZT00, ZT01, ZT10, ZT11, ZH00, ZH01, ZH10, ZH11, ZA00, ZA01, ZA10, ZA11
      , ZP00, ZP01, ZP10, ZP11, ReZN00, ImZN00, ReZN01, ImZN01, ReZN02, ImZN02,
      ReZN03, ImZN03, ReZN10, ImZN10, ReZN11, ImZN11, ReZN12, ImZN12, ReZN13,
      ImZN13, ReZN20, ImZN20, ReZN21, ImZN21, ReZN22, ImZN22, ReZN23, ImZN23,
      ReZN30, ImZN30, ReZN31, ImZN31, ReZN32, ImZN32, ReZN33, ImZN33, ReUM00,
      ImUM00, ReUM01, ImUM01, ReUM10, ImUM10, ReUM11, ImUM11, ReUP00, ImUP00,
      ReUP01, ImUP01, ReUP10, ImUP10, ReUP11, ImUP11, NUMBER_OF_MIXINGS};

   enum Input_parameters : unsigned {TanBeta, SignMu, M1, M2, M3, AtIN, AbIN,
      AtauIN, AcIN, AsIN, AmuonIN, AuIN, AdIN, AeIN, mHd2IN, mHu2IN, ml11IN,
      ml22IN, ml33IN, me11IN, me22IN, me33IN, mq11IN, mq22IN, mq33IN, mu11IN,
      mu22IN, mu33IN, md11IN, md22IN, md33IN, NUMBER_OF_INPUT_PARAMETERS};

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
