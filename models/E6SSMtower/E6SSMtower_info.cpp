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

// File generated at Mon 19 Sep 2016 09:41:25

#include "E6SSMtower_info.hpp"

#include <iostream>

namespace flexiblesusy {

namespace E6SSMtower_info {
   const double normalization_g1 = 0.7745966692414834;
   const double normalization_g2 = 1;
   const double normalization_g3 = 1;
   const double normalization_gN = 1;

   const unsigned particle_multiplicities[NUMBER_OF_PARTICLES] = {1, 1, 3, 1, 6
      , 3, 6, 6, 6, 3, 3, 2, 6, 2, 3, 3, 3, 3, 4, 4, 2, 4, 2, 2, 2, 2, 2, 1, 1, 1,
      1};

   const char* particle_names[NUMBER_OF_PARTICLES] = {"VG", "Glu", "Fv", "ChaP"
      , "Sd", "Sv", "Su", "Se", "SDX", "hh", "Ah", "Hpm", "Chi", "Cha", "Fe", "Fd"
      , "Fu", "FDX", "SHI0", "SHIp", "ChaI", "ChiI", "SSI0", "FSI", "SHp0", "SHpp"
      , "ChiP", "VWm", "VP", "VZ", "VZp"};

   const char* particle_latex_names[NUMBER_OF_PARTICLES] = {   "g",
      "\\tilde{g}", "\\nu", "\\tilde{\\chi}^{'-}", "\\tilde{d}", "\\tilde{\\nu}",
      "\\tilde{u}", "\\tilde{e}", "\\tilde{x}", "h", "A^0", "H^-",
      "\\tilde{\\chi}^0", "\\tilde{\\chi}^-", "e", "d", "u", "x", "h^{0,Inert}",
      "h^{-,Inert}", "\\tilde{\\chi}^{-,Inert}", "\\tilde{\\chi}^{0,Inert}",
      "s^{Inert}", "\\tilde{S}^{Inert}", "H^{'0}", "H^{'-}", "\\tilde{\\chi}^{'0}"
      , "W^-", "\\gamma", "Z", "{Z'}"};

   const char* parameter_names[NUMBER_OF_PARAMETERS] = {"Yd(0,0)", "Yd(0,1)",
      "Yd(0,2)", "Yd(1,0)", "Yd(1,1)", "Yd(1,2)", "Yd(2,0)", "Yd(2,1)", "Yd(2,2)",
      "Ye(0,0)", "Ye(0,1)", "Ye(0,2)", "Ye(1,0)", "Ye(1,1)", "Ye(1,2)", "Ye(2,0)"
      , "Ye(2,1)", "Ye(2,2)", "Kappa(0,0)", "Kappa(0,1)", "Kappa(0,2)",
      "Kappa(1,0)", "Kappa(1,1)", "Kappa(1,2)", "Kappa(2,0)", "Kappa(2,1)",
      "Kappa(2,2)", "Lambda12(0,0)", "Lambda12(0,1)", "Lambda12(1,0)",
      "Lambda12(1,1)", "Lambdax", "Yu(0,0)", "Yu(0,1)", "Yu(0,2)", "Yu(1,0)",
      "Yu(1,1)", "Yu(1,2)", "Yu(2,0)", "Yu(2,1)", "Yu(2,2)", "MuPr", "g1", "g2",
      "g3", "gN", "vd", "vu", "vs", "TYd(0,0)", "TYd(0,1)", "TYd(0,2)", "TYd(1,0)"
      , "TYd(1,1)", "TYd(1,2)", "TYd(2,0)", "TYd(2,1)", "TYd(2,2)", "TYe(0,0)",
      "TYe(0,1)", "TYe(0,2)", "TYe(1,0)", "TYe(1,1)", "TYe(1,2)", "TYe(2,0)",
      "TYe(2,1)", "TYe(2,2)", "TKappa(0,0)", "TKappa(0,1)", "TKappa(0,2)",
      "TKappa(1,0)", "TKappa(1,1)", "TKappa(1,2)", "TKappa(2,0)", "TKappa(2,1)",
      "TKappa(2,2)", "TLambda12(0,0)", "TLambda12(0,1)", "TLambda12(1,0)",
      "TLambda12(1,1)", "TLambdax", "TYu(0,0)", "TYu(0,1)", "TYu(0,2)", "TYu(1,0)"
      , "TYu(1,1)", "TYu(1,2)", "TYu(2,0)", "TYu(2,1)", "TYu(2,2)", "BMuPr",
      "mq2(0,0)", "mq2(0,1)", "mq2(0,2)", "mq2(1,0)", "mq2(1,1)", "mq2(1,2)",
      "mq2(2,0)", "mq2(2,1)", "mq2(2,2)", "ml2(0,0)", "ml2(0,1)", "ml2(0,2)",
      "ml2(1,0)", "ml2(1,1)", "ml2(1,2)", "ml2(2,0)", "ml2(2,1)", "ml2(2,2)",
      "mHd2", "mHu2", "md2(0,0)", "md2(0,1)", "md2(0,2)", "md2(1,0)", "md2(1,1)",
      "md2(1,2)", "md2(2,0)", "md2(2,1)", "md2(2,2)", "mu2(0,0)", "mu2(0,1)",
      "mu2(0,2)", "mu2(1,0)", "mu2(1,1)", "mu2(1,2)", "mu2(2,0)", "mu2(2,1)",
      "mu2(2,2)", "me2(0,0)", "me2(0,1)", "me2(0,2)", "me2(1,0)", "me2(1,1)",
      "me2(1,2)", "me2(2,0)", "me2(2,1)", "me2(2,2)", "ms2", "mH1I2(0,0)",
      "mH1I2(0,1)", "mH1I2(1,0)", "mH1I2(1,1)", "mH2I2(0,0)", "mH2I2(0,1)",
      "mH2I2(1,0)", "mH2I2(1,1)", "msI2(0,0)", "msI2(0,1)", "msI2(1,0)",
      "msI2(1,1)", "mDx2(0,0)", "mDx2(0,1)", "mDx2(0,2)", "mDx2(1,0)", "mDx2(1,1)"
      , "mDx2(1,2)", "mDx2(2,0)", "mDx2(2,1)", "mDx2(2,2)", "mDxbar2(0,0)",
      "mDxbar2(0,1)", "mDxbar2(0,2)", "mDxbar2(1,0)", "mDxbar2(1,1)",
      "mDxbar2(1,2)", "mDxbar2(2,0)", "mDxbar2(2,1)", "mDxbar2(2,2)", "mHp2",
      "mHpbar2", "MassB", "MassWB", "MassG", "MassBp"};

   const char* particle_mixing_names[NUMBER_OF_MIXINGS] = {   "ZD(0,0)",
      "ZD(0,1)", "ZD(0,2)", "ZD(0,3)", "ZD(0,4)", "ZD(0,5)", "ZD(1,0)", "ZD(1,1)",
      "ZD(1,2)", "ZD(1,3)", "ZD(1,4)", "ZD(1,5)", "ZD(2,0)", "ZD(2,1)", "ZD(2,2)"
      , "ZD(2,3)", "ZD(2,4)", "ZD(2,5)", "ZD(3,0)", "ZD(3,1)", "ZD(3,2)",
      "ZD(3,3)", "ZD(3,4)", "ZD(3,5)", "ZD(4,0)", "ZD(4,1)", "ZD(4,2)", "ZD(4,3)",
      "ZD(4,4)", "ZD(4,5)", "ZD(5,0)", "ZD(5,1)", "ZD(5,2)", "ZD(5,3)", "ZD(5,4)"
      , "ZD(5,5)", "ZV(0,0)", "ZV(0,1)", "ZV(0,2)", "ZV(1,0)", "ZV(1,1)",
      "ZV(1,2)", "ZV(2,0)", "ZV(2,1)", "ZV(2,2)", "ZU(0,0)", "ZU(0,1)", "ZU(0,2)",
      "ZU(0,3)", "ZU(0,4)", "ZU(0,5)", "ZU(1,0)", "ZU(1,1)", "ZU(1,2)", "ZU(1,3)"
      , "ZU(1,4)", "ZU(1,5)", "ZU(2,0)", "ZU(2,1)", "ZU(2,2)", "ZU(2,3)",
      "ZU(2,4)", "ZU(2,5)", "ZU(3,0)", "ZU(3,1)", "ZU(3,2)", "ZU(3,3)", "ZU(3,4)",
      "ZU(3,5)", "ZU(4,0)", "ZU(4,1)", "ZU(4,2)", "ZU(4,3)", "ZU(4,4)", "ZU(4,5)"
      , "ZU(5,0)", "ZU(5,1)", "ZU(5,2)", "ZU(5,3)", "ZU(5,4)", "ZU(5,5)",
      "ZE(0,0)", "ZE(0,1)", "ZE(0,2)", "ZE(0,3)", "ZE(0,4)", "ZE(0,5)", "ZE(1,0)",
      "ZE(1,1)", "ZE(1,2)", "ZE(1,3)", "ZE(1,4)", "ZE(1,5)", "ZE(2,0)", "ZE(2,1)"
      , "ZE(2,2)", "ZE(2,3)", "ZE(2,4)", "ZE(2,5)", "ZE(3,0)", "ZE(3,1)",
      "ZE(3,2)", "ZE(3,3)", "ZE(3,4)", "ZE(3,5)", "ZE(4,0)", "ZE(4,1)", "ZE(4,2)",
      "ZE(4,3)", "ZE(4,4)", "ZE(4,5)", "ZE(5,0)", "ZE(5,1)", "ZE(5,2)", "ZE(5,3)"
      , "ZE(5,4)", "ZE(5,5)", "ZDX(0,0)", "ZDX(0,1)", "ZDX(0,2)", "ZDX(0,3)",
      "ZDX(0,4)", "ZDX(0,5)", "ZDX(1,0)", "ZDX(1,1)", "ZDX(1,2)", "ZDX(1,3)",
      "ZDX(1,4)", "ZDX(1,5)", "ZDX(2,0)", "ZDX(2,1)", "ZDX(2,2)", "ZDX(2,3)",
      "ZDX(2,4)", "ZDX(2,5)", "ZDX(3,0)", "ZDX(3,1)", "ZDX(3,2)", "ZDX(3,3)",
      "ZDX(3,4)", "ZDX(3,5)", "ZDX(4,0)", "ZDX(4,1)", "ZDX(4,2)", "ZDX(4,3)",
      "ZDX(4,4)", "ZDX(4,5)", "ZDX(5,0)", "ZDX(5,1)", "ZDX(5,2)", "ZDX(5,3)",
      "ZDX(5,4)", "ZDX(5,5)", "ZH(0,0)", "ZH(0,1)", "ZH(0,2)", "ZH(1,0)",
      "ZH(1,1)", "ZH(1,2)", "ZH(2,0)", "ZH(2,1)", "ZH(2,2)", "ZA(0,0)", "ZA(0,1)",
      "ZA(0,2)", "ZA(1,0)", "ZA(1,1)", "ZA(1,2)", "ZA(2,0)", "ZA(2,1)", "ZA(2,2)"
      , "ZP(0,0)", "ZP(0,1)", "ZP(1,0)", "ZP(1,1)", "Re(ZN(0,0))", "Im(ZN(0,0))",
      "Re(ZN(0,1))", "Im(ZN(0,1))", "Re(ZN(0,2))", "Im(ZN(0,2))", "Re(ZN(0,3))",
      "Im(ZN(0,3))", "Re(ZN(0,4))", "Im(ZN(0,4))", "Re(ZN(0,5))", "Im(ZN(0,5))",
      "Re(ZN(1,0))", "Im(ZN(1,0))", "Re(ZN(1,1))", "Im(ZN(1,1))", "Re(ZN(1,2))",
      "Im(ZN(1,2))", "Re(ZN(1,3))", "Im(ZN(1,3))", "Re(ZN(1,4))", "Im(ZN(1,4))",
      "Re(ZN(1,5))", "Im(ZN(1,5))", "Re(ZN(2,0))", "Im(ZN(2,0))", "Re(ZN(2,1))",
      "Im(ZN(2,1))", "Re(ZN(2,2))", "Im(ZN(2,2))", "Re(ZN(2,3))", "Im(ZN(2,3))",
      "Re(ZN(2,4))", "Im(ZN(2,4))", "Re(ZN(2,5))", "Im(ZN(2,5))", "Re(ZN(3,0))",
      "Im(ZN(3,0))", "Re(ZN(3,1))", "Im(ZN(3,1))", "Re(ZN(3,2))", "Im(ZN(3,2))",
      "Re(ZN(3,3))", "Im(ZN(3,3))", "Re(ZN(3,4))", "Im(ZN(3,4))", "Re(ZN(3,5))",
      "Im(ZN(3,5))", "Re(ZN(4,0))", "Im(ZN(4,0))", "Re(ZN(4,1))", "Im(ZN(4,1))",
      "Re(ZN(4,2))", "Im(ZN(4,2))", "Re(ZN(4,3))", "Im(ZN(4,3))", "Re(ZN(4,4))",
      "Im(ZN(4,4))", "Re(ZN(4,5))", "Im(ZN(4,5))", "Re(ZN(5,0))", "Im(ZN(5,0))",
      "Re(ZN(5,1))", "Im(ZN(5,1))", "Re(ZN(5,2))", "Im(ZN(5,2))", "Re(ZN(5,3))",
      "Im(ZN(5,3))", "Re(ZN(5,4))", "Im(ZN(5,4))", "Re(ZN(5,5))", "Im(ZN(5,5))",
      "Re(UM(0,0))", "Im(UM(0,0))", "Re(UM(0,1))", "Im(UM(0,1))", "Re(UM(1,0))",
      "Im(UM(1,0))", "Re(UM(1,1))", "Im(UM(1,1))", "Re(UP(0,0))", "Im(UP(0,0))",
      "Re(UP(0,1))", "Im(UP(0,1))", "Re(UP(1,0))", "Im(UP(1,0))", "Re(UP(1,1))",
      "Im(UP(1,1))", "Re(ZEL(0,0))", "Im(ZEL(0,0))", "Re(ZEL(0,1))",
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
      "Im(ZUR(2,2))", "Re(ZDXL(0,0))", "Im(ZDXL(0,0))", "Re(ZDXL(0,1))",
      "Im(ZDXL(0,1))", "Re(ZDXL(0,2))", "Im(ZDXL(0,2))", "Re(ZDXL(1,0))",
      "Im(ZDXL(1,0))", "Re(ZDXL(1,1))", "Im(ZDXL(1,1))", "Re(ZDXL(1,2))",
      "Im(ZDXL(1,2))", "Re(ZDXL(2,0))", "Im(ZDXL(2,0))", "Re(ZDXL(2,1))",
      "Im(ZDXL(2,1))", "Re(ZDXL(2,2))", "Im(ZDXL(2,2))", "Re(ZDXR(0,0))",
      "Im(ZDXR(0,0))", "Re(ZDXR(0,1))", "Im(ZDXR(0,1))", "Re(ZDXR(0,2))",
      "Im(ZDXR(0,2))", "Re(ZDXR(1,0))", "Im(ZDXR(1,0))", "Re(ZDXR(1,1))",
      "Im(ZDXR(1,1))", "Re(ZDXR(1,2))", "Im(ZDXR(1,2))", "Re(ZDXR(2,0))",
      "Im(ZDXR(2,0))", "Re(ZDXR(2,1))", "Im(ZDXR(2,1))", "Re(ZDXR(2,2))",
      "Im(ZDXR(2,2))", "UHI0(0,0)", "UHI0(0,1)", "UHI0(0,2)", "UHI0(0,3)",
      "UHI0(1,0)", "UHI0(1,1)", "UHI0(1,2)", "UHI0(1,3)", "UHI0(2,0)", "UHI0(2,1)"
      , "UHI0(2,2)", "UHI0(2,3)", "UHI0(3,0)", "UHI0(3,1)", "UHI0(3,2)",
      "UHI0(3,3)", "UHIp(0,0)", "UHIp(0,1)", "UHIp(0,2)", "UHIp(0,3)", "UHIp(1,0)"
      , "UHIp(1,1)", "UHIp(1,2)", "UHIp(1,3)", "UHIp(2,0)", "UHIp(2,1)",
      "UHIp(2,2)", "UHIp(2,3)", "UHIp(3,0)", "UHIp(3,1)", "UHIp(3,2)", "UHIp(3,3)"
      , "Re(ZMI(0,0))", "Im(ZMI(0,0))", "Re(ZMI(0,1))", "Im(ZMI(0,1))",
      "Re(ZMI(1,0))", "Im(ZMI(1,0))", "Re(ZMI(1,1))", "Im(ZMI(1,1))",
      "Re(ZPI(0,0))", "Im(ZPI(0,0))", "Re(ZPI(0,1))", "Im(ZPI(0,1))",
      "Re(ZPI(1,0))", "Im(ZPI(1,0))", "Re(ZPI(1,1))", "Im(ZPI(1,1))",
      "Re(ZNI(0,0))", "Im(ZNI(0,0))", "Re(ZNI(0,1))", "Im(ZNI(0,1))",
      "Re(ZNI(0,2))", "Im(ZNI(0,2))", "Re(ZNI(0,3))", "Im(ZNI(0,3))",
      "Re(ZNI(1,0))", "Im(ZNI(1,0))", "Re(ZNI(1,1))", "Im(ZNI(1,1))",
      "Re(ZNI(1,2))", "Im(ZNI(1,2))", "Re(ZNI(1,3))", "Im(ZNI(1,3))",
      "Re(ZNI(2,0))", "Im(ZNI(2,0))", "Re(ZNI(2,1))", "Im(ZNI(2,1))",
      "Re(ZNI(2,2))", "Im(ZNI(2,2))", "Re(ZNI(2,3))", "Im(ZNI(2,3))",
      "Re(ZNI(3,0))", "Im(ZNI(3,0))", "Re(ZNI(3,1))", "Im(ZNI(3,1))",
      "Re(ZNI(3,2))", "Im(ZNI(3,2))", "Re(ZNI(3,3))", "Im(ZNI(3,3))", "ZSSI(0,0)",
      "ZSSI(0,1)", "ZSSI(1,0)", "ZSSI(1,1)", "Re(ZFSI(0,0))", "Im(ZFSI(0,0))",
      "Re(ZFSI(0,1))", "Im(ZFSI(0,1))", "Re(ZFSI(1,0))", "Im(ZFSI(1,0))",
      "Re(ZFSI(1,1))", "Im(ZFSI(1,1))", "UHp0(0,0)", "UHp0(0,1)", "UHp0(1,0)",
      "UHp0(1,1)", "UHpp(0,0)", "UHpp(0,1)", "UHpp(1,0)", "UHpp(1,1)",
      "Re(ZNp(0,0))", "Im(ZNp(0,0))", "Re(ZNp(0,1))", "Im(ZNp(0,1))",
      "Re(ZNp(1,0))", "Im(ZNp(1,0))", "Re(ZNp(1,1))", "Im(ZNp(1,1))", "ZZ(0,0)",
      "ZZ(0,1)", "ZZ(0,2)", "ZZ(1,0)", "ZZ(1,1)", "ZZ(1,2)", "ZZ(2,0)", "ZZ(2,1)",
      "ZZ(2,2)"};

   const char* input_parameter_names[NUMBER_OF_INPUT_PARAMETERS] = {"MSUSY",
      "M1Input", "M2Input", "M3Input", "MuInput", "mAInput", "TanBeta",
      "LambdaInput", "gNInput", "M4Input", "mHp2Input", "mHpbar2Input",
      "MuPrInput", "BMuPrInput", "mq2Input(0,0)", "mq2Input(0,1)", "mq2Input(0,2)"
      , "mq2Input(1,0)", "mq2Input(1,1)", "mq2Input(1,2)", "mq2Input(2,0)",
      "mq2Input(2,1)", "mq2Input(2,2)", "mu2Input(0,0)", "mu2Input(0,1)",
      "mu2Input(0,2)", "mu2Input(1,0)", "mu2Input(1,1)", "mu2Input(1,2)",
      "mu2Input(2,0)", "mu2Input(2,1)", "mu2Input(2,2)", "md2Input(0,0)",
      "md2Input(0,1)", "md2Input(0,2)", "md2Input(1,0)", "md2Input(1,1)",
      "md2Input(1,2)", "md2Input(2,0)", "md2Input(2,1)", "md2Input(2,2)",
      "ml2Input(0,0)", "ml2Input(0,1)", "ml2Input(0,2)", "ml2Input(1,0)",
      "ml2Input(1,1)", "ml2Input(1,2)", "ml2Input(2,0)", "ml2Input(2,1)",
      "ml2Input(2,2)", "me2Input(0,0)", "me2Input(0,1)", "me2Input(0,2)",
      "me2Input(1,0)", "me2Input(1,1)", "me2Input(1,2)", "me2Input(2,0)",
      "me2Input(2,1)", "me2Input(2,2)", "AuInput(0,0)", "AuInput(0,1)",
      "AuInput(0,2)", "AuInput(1,0)", "AuInput(1,1)", "AuInput(1,2)",
      "AuInput(2,0)", "AuInput(2,1)", "AuInput(2,2)", "AdInput(0,0)",
      "AdInput(0,1)", "AdInput(0,2)", "AdInput(1,0)", "AdInput(1,1)",
      "AdInput(1,2)", "AdInput(2,0)", "AdInput(2,1)", "AdInput(2,2)",
      "AeInput(0,0)", "AeInput(0,1)", "AeInput(0,2)", "AeInput(1,0)",
      "AeInput(1,1)", "AeInput(1,2)", "AeInput(2,0)", "AeInput(2,1)",
      "AeInput(2,2)", "Lambda12Input(0,0)", "Lambda12Input(0,1)",
      "Lambda12Input(1,0)", "Lambda12Input(1,1)", "ALambda12Input(0,0)",
      "ALambda12Input(0,1)", "ALambda12Input(1,0)", "ALambda12Input(1,1)",
      "KappaInput(0,0)", "KappaInput(0,1)", "KappaInput(0,2)", "KappaInput(1,0)",
      "KappaInput(1,1)", "KappaInput(1,2)", "KappaInput(2,0)", "KappaInput(2,1)",
      "KappaInput(2,2)", "AKappaInput(0,0)", "AKappaInput(0,1)",
      "AKappaInput(0,2)", "AKappaInput(1,0)", "AKappaInput(1,1)",
      "AKappaInput(1,2)", "AKappaInput(2,0)", "AKappaInput(2,1)",
      "AKappaInput(2,2)", "mDx2Input(0,0)", "mDx2Input(0,1)", "mDx2Input(0,2)",
      "mDx2Input(1,0)", "mDx2Input(1,1)", "mDx2Input(1,2)", "mDx2Input(2,0)",
      "mDx2Input(2,1)", "mDx2Input(2,2)", "mDxbar2Input(0,0)", "mDxbar2Input(0,1)"
      , "mDxbar2Input(0,2)", "mDxbar2Input(1,0)", "mDxbar2Input(1,1)",
      "mDxbar2Input(1,2)", "mDxbar2Input(2,0)", "mDxbar2Input(2,1)",
      "mDxbar2Input(2,2)", "mH1I2Input(0,0)", "mH1I2Input(0,1)", "mH1I2Input(1,0)"
      , "mH1I2Input(1,1)", "mH2I2Input(0,0)", "mH2I2Input(0,1)", "mH2I2Input(1,0)"
      , "mH2I2Input(1,1)", "msI2Input(0,0)", "msI2Input(0,1)", "msI2Input(1,0)",
      "msI2Input(1,1)"};

   const char* model_name = "E6SSMtower";
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

} // namespace E6SSMtower_info

} // namespace flexiblesusy
