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


#include "CE6SSM_info.hpp"

#include "error.hpp"

#include <iostream>
#include <vector>

namespace flexiblesusy {

namespace CE6SSM_info {
   const double normalization_g1 = 0.7745966692414834;
   const double normalization_g2 = 1;
   const double normalization_g3 = 1;
   const double normalization_gN = 1;

   const std::array<int, NUMBER_OF_PARTICLES> particle_multiplicities = {1, 1, 3,
      1, 6, 3, 6, 6, 6, 3, 3, 2, 6, 2, 3, 3, 3, 3, 4, 4, 2, 4, 2, 2, 2, 2, 2, 1, 1
      , 1, 1};

   const std::array<std::string, NUMBER_OF_PARTICLES> particle_names = {"VG",
      "Glu", "Fv", "ChaP", "Sd", "Sv", "Su", "Se", "SDX", "hh", "Ah", "Hpm", "Chi"
      , "Cha", "Fe", "Fd", "Fu", "FDX", "SHI0", "SHIp", "ChaI", "ChiI", "SSI0",
      "FSI", "SHp0", "SHpp", "ChiP", "VWm", "VP", "VZ", "VZp"};

   const std::array<std::string, NUMBER_OF_PARTICLES> particle_latex_names = {
      "g", "\\tilde{g}", "\\nu", "\\tilde{\\chi}^{'-}", "\\tilde{d}",
      "\\tilde{\\nu}", "\\tilde{u}", "\\tilde{e}", "\\tilde{x}", "h", "A^0", "H^-"
      , "\\tilde{\\chi}^0", "\\tilde{\\chi}^-", "e", "d", "u", "x", "h^{0,Inert}",
      "h^{-,Inert}", "\\tilde{\\chi}^{-,Inert}", "\\tilde{\\chi}^{0,Inert}",
      "s^{Inert}", "\\tilde{S}^{Inert}", "H^{'0}", "H^{'-}", "\\tilde{\\chi}^{'0}"
      , "W^-", "\\gamma", "Z", "{Z'}"};

   const std::array<std::string, NUMBER_OF_PARAMETERS> parameter_names = {
      "Yd(0,0)", "Yd(0,1)", "Yd(0,2)", "Yd(1,0)", "Yd(1,1)", "Yd(1,2)", "Yd(2,0)",
      "Yd(2,1)", "Yd(2,2)", "Ye(0,0)", "Ye(0,1)", "Ye(0,2)", "Ye(1,0)", "Ye(1,1)",
      "Ye(1,2)", "Ye(2,0)", "Ye(2,1)", "Ye(2,2)", "Kappa(0,0)", "Kappa(0,1)",
      "Kappa(0,2)", "Kappa(1,0)", "Kappa(1,1)", "Kappa(1,2)", "Kappa(2,0)",
      "Kappa(2,1)", "Kappa(2,2)", "Lambda12(0,0)", "Lambda12(0,1)",
      "Lambda12(1,0)", "Lambda12(1,1)", "Lambdax", "Yu(0,0)", "Yu(0,1)", "Yu(0,2)"
      , "Yu(1,0)", "Yu(1,1)", "Yu(1,2)", "Yu(2,0)", "Yu(2,1)", "Yu(2,2)", "MuPr",
      "g1", "g2", "g3", "gN", "vd", "vu", "vs", "TYd(0,0)", "TYd(0,1)", "TYd(0,2)"
      , "TYd(1,0)", "TYd(1,1)", "TYd(1,2)", "TYd(2,0)", "TYd(2,1)", "TYd(2,2)",
      "TYe(0,0)", "TYe(0,1)", "TYe(0,2)", "TYe(1,0)", "TYe(1,1)", "TYe(1,2)",
      "TYe(2,0)", "TYe(2,1)", "TYe(2,2)", "TKappa(0,0)", "TKappa(0,1)",
      "TKappa(0,2)", "TKappa(1,0)", "TKappa(1,1)", "TKappa(1,2)", "TKappa(2,0)",
      "TKappa(2,1)", "TKappa(2,2)", "TLambda12(0,0)", "TLambda12(0,1)",
      "TLambda12(1,0)", "TLambda12(1,1)", "TLambdax", "TYu(0,0)", "TYu(0,1)",
      "TYu(0,2)", "TYu(1,0)", "TYu(1,1)", "TYu(1,2)", "TYu(2,0)", "TYu(2,1)",
      "TYu(2,2)", "BMuPr", "mq2(0,0)", "mq2(0,1)", "mq2(0,2)", "mq2(1,0)",
      "mq2(1,1)", "mq2(1,2)", "mq2(2,0)", "mq2(2,1)", "mq2(2,2)", "ml2(0,0)",
      "ml2(0,1)", "ml2(0,2)", "ml2(1,0)", "ml2(1,1)", "ml2(1,2)", "ml2(2,0)",
      "ml2(2,1)", "ml2(2,2)", "mHd2", "mHu2", "md2(0,0)", "md2(0,1)", "md2(0,2)",
      "md2(1,0)", "md2(1,1)", "md2(1,2)", "md2(2,0)", "md2(2,1)", "md2(2,2)",
      "mu2(0,0)", "mu2(0,1)", "mu2(0,2)", "mu2(1,0)", "mu2(1,1)", "mu2(1,2)",
      "mu2(2,0)", "mu2(2,1)", "mu2(2,2)", "me2(0,0)", "me2(0,1)", "me2(0,2)",
      "me2(1,0)", "me2(1,1)", "me2(1,2)", "me2(2,0)", "me2(2,1)", "me2(2,2)",
      "ms2", "mH1I2(0,0)", "mH1I2(0,1)", "mH1I2(1,0)", "mH1I2(1,1)", "mH2I2(0,0)",
      "mH2I2(0,1)", "mH2I2(1,0)", "mH2I2(1,1)", "msI2(0,0)", "msI2(0,1)",
      "msI2(1,0)", "msI2(1,1)", "mDx2(0,0)", "mDx2(0,1)", "mDx2(0,2)", "mDx2(1,0)"
      , "mDx2(1,1)", "mDx2(1,2)", "mDx2(2,0)", "mDx2(2,1)", "mDx2(2,2)",
      "mDxbar2(0,0)", "mDxbar2(0,1)", "mDxbar2(0,2)", "mDxbar2(1,0)",
      "mDxbar2(1,1)", "mDxbar2(1,2)", "mDxbar2(2,0)", "mDxbar2(2,1)",
      "mDxbar2(2,2)", "mHp2", "mHpbar2", "MassB", "MassWB", "MassG", "MassBp"};

   const std::array<std::string, NUMBER_OF_MIXINGS> particle_mixing_names = {
      "ZD(0,0)", "ZD(0,1)", "ZD(0,2)", "ZD(0,3)", "ZD(0,4)", "ZD(0,5)", "ZD(1,0)",
      "ZD(1,1)", "ZD(1,2)", "ZD(1,3)", "ZD(1,4)", "ZD(1,5)", "ZD(2,0)", "ZD(2,1)",
      "ZD(2,2)", "ZD(2,3)", "ZD(2,4)", "ZD(2,5)", "ZD(3,0)", "ZD(3,1)", "ZD(3,2)",
      "ZD(3,3)", "ZD(3,4)", "ZD(3,5)", "ZD(4,0)", "ZD(4,1)", "ZD(4,2)", "ZD(4,3)",
      "ZD(4,4)", "ZD(4,5)", "ZD(5,0)", "ZD(5,1)", "ZD(5,2)", "ZD(5,3)", "ZD(5,4)",
      "ZD(5,5)", "ZV(0,0)", "ZV(0,1)", "ZV(0,2)", "ZV(1,0)", "ZV(1,1)", "ZV(1,2)",
      "ZV(2,0)", "ZV(2,1)", "ZV(2,2)", "ZU(0,0)", "ZU(0,1)", "ZU(0,2)", "ZU(0,3)",
      "ZU(0,4)", "ZU(0,5)", "ZU(1,0)", "ZU(1,1)", "ZU(1,2)", "ZU(1,3)", "ZU(1,4)",
      "ZU(1,5)", "ZU(2,0)", "ZU(2,1)", "ZU(2,2)", "ZU(2,3)", "ZU(2,4)", "ZU(2,5)",
      "ZU(3,0)", "ZU(3,1)", "ZU(3,2)", "ZU(3,3)", "ZU(3,4)", "ZU(3,5)", "ZU(4,0)",
      "ZU(4,1)", "ZU(4,2)", "ZU(4,3)", "ZU(4,4)", "ZU(4,5)", "ZU(5,0)", "ZU(5,1)",
      "ZU(5,2)", "ZU(5,3)", "ZU(5,4)", "ZU(5,5)", "ZE(0,0)", "ZE(0,1)", "ZE(0,2)",
      "ZE(0,3)", "ZE(0,4)", "ZE(0,5)", "ZE(1,0)", "ZE(1,1)", "ZE(1,2)", "ZE(1,3)",
      "ZE(1,4)", "ZE(1,5)", "ZE(2,0)", "ZE(2,1)", "ZE(2,2)", "ZE(2,3)", "ZE(2,4)",
      "ZE(2,5)", "ZE(3,0)", "ZE(3,1)", "ZE(3,2)", "ZE(3,3)", "ZE(3,4)", "ZE(3,5)",
      "ZE(4,0)", "ZE(4,1)", "ZE(4,2)", "ZE(4,3)", "ZE(4,4)", "ZE(4,5)", "ZE(5,0)",
      "ZE(5,1)", "ZE(5,2)", "ZE(5,3)", "ZE(5,4)", "ZE(5,5)", "ZDX(0,0)",
      "ZDX(0,1)", "ZDX(0,2)", "ZDX(0,3)", "ZDX(0,4)", "ZDX(0,5)", "ZDX(1,0)",
      "ZDX(1,1)", "ZDX(1,2)", "ZDX(1,3)", "ZDX(1,4)", "ZDX(1,5)", "ZDX(2,0)",
      "ZDX(2,1)", "ZDX(2,2)", "ZDX(2,3)", "ZDX(2,4)", "ZDX(2,5)", "ZDX(3,0)",
      "ZDX(3,1)", "ZDX(3,2)", "ZDX(3,3)", "ZDX(3,4)", "ZDX(3,5)", "ZDX(4,0)",
      "ZDX(4,1)", "ZDX(4,2)", "ZDX(4,3)", "ZDX(4,4)", "ZDX(4,5)", "ZDX(5,0)",
      "ZDX(5,1)", "ZDX(5,2)", "ZDX(5,3)", "ZDX(5,4)", "ZDX(5,5)", "ZH(0,0)",
      "ZH(0,1)", "ZH(0,2)", "ZH(1,0)", "ZH(1,1)", "ZH(1,2)", "ZH(2,0)", "ZH(2,1)",
      "ZH(2,2)", "ZA(0,0)", "ZA(0,1)", "ZA(0,2)", "ZA(1,0)", "ZA(1,1)", "ZA(1,2)",
      "ZA(2,0)", "ZA(2,1)", "ZA(2,2)", "ZP(0,0)", "ZP(0,1)", "ZP(1,0)", "ZP(1,1)",
      "Re(ZN(0,0))", "Im(ZN(0,0))", "Re(ZN(0,1))", "Im(ZN(0,1))", "Re(ZN(0,2))",
      "Im(ZN(0,2))", "Re(ZN(0,3))", "Im(ZN(0,3))", "Re(ZN(0,4))", "Im(ZN(0,4))",
      "Re(ZN(0,5))", "Im(ZN(0,5))", "Re(ZN(1,0))", "Im(ZN(1,0))", "Re(ZN(1,1))",
      "Im(ZN(1,1))", "Re(ZN(1,2))", "Im(ZN(1,2))", "Re(ZN(1,3))", "Im(ZN(1,3))",
      "Re(ZN(1,4))", "Im(ZN(1,4))", "Re(ZN(1,5))", "Im(ZN(1,5))", "Re(ZN(2,0))",
      "Im(ZN(2,0))", "Re(ZN(2,1))", "Im(ZN(2,1))", "Re(ZN(2,2))", "Im(ZN(2,2))",
      "Re(ZN(2,3))", "Im(ZN(2,3))", "Re(ZN(2,4))", "Im(ZN(2,4))", "Re(ZN(2,5))",
      "Im(ZN(2,5))", "Re(ZN(3,0))", "Im(ZN(3,0))", "Re(ZN(3,1))", "Im(ZN(3,1))",
      "Re(ZN(3,2))", "Im(ZN(3,2))", "Re(ZN(3,3))", "Im(ZN(3,3))", "Re(ZN(3,4))",
      "Im(ZN(3,4))", "Re(ZN(3,5))", "Im(ZN(3,5))", "Re(ZN(4,0))", "Im(ZN(4,0))",
      "Re(ZN(4,1))", "Im(ZN(4,1))", "Re(ZN(4,2))", "Im(ZN(4,2))", "Re(ZN(4,3))",
      "Im(ZN(4,3))", "Re(ZN(4,4))", "Im(ZN(4,4))", "Re(ZN(4,5))", "Im(ZN(4,5))",
      "Re(ZN(5,0))", "Im(ZN(5,0))", "Re(ZN(5,1))", "Im(ZN(5,1))", "Re(ZN(5,2))",
      "Im(ZN(5,2))", "Re(ZN(5,3))", "Im(ZN(5,3))", "Re(ZN(5,4))", "Im(ZN(5,4))",
      "Re(ZN(5,5))", "Im(ZN(5,5))", "Re(UM(0,0))", "Im(UM(0,0))", "Re(UM(0,1))",
      "Im(UM(0,1))", "Re(UM(1,0))", "Im(UM(1,0))", "Re(UM(1,1))", "Im(UM(1,1))",
      "Re(UP(0,0))", "Im(UP(0,0))", "Re(UP(0,1))", "Im(UP(0,1))", "Re(UP(1,0))",
      "Im(UP(1,0))", "Re(UP(1,1))", "Im(UP(1,1))", "Re(ZEL(0,0))", "Im(ZEL(0,0))",
      "Re(ZEL(0,1))", "Im(ZEL(0,1))", "Re(ZEL(0,2))", "Im(ZEL(0,2))",
      "Re(ZEL(1,0))", "Im(ZEL(1,0))", "Re(ZEL(1,1))", "Im(ZEL(1,1))",
      "Re(ZEL(1,2))", "Im(ZEL(1,2))", "Re(ZEL(2,0))", "Im(ZEL(2,0))",
      "Re(ZEL(2,1))", "Im(ZEL(2,1))", "Re(ZEL(2,2))", "Im(ZEL(2,2))",
      "Re(ZER(0,0))", "Im(ZER(0,0))", "Re(ZER(0,1))", "Im(ZER(0,1))",
      "Re(ZER(0,2))", "Im(ZER(0,2))", "Re(ZER(1,0))", "Im(ZER(1,0))",
      "Re(ZER(1,1))", "Im(ZER(1,1))", "Re(ZER(1,2))", "Im(ZER(1,2))",
      "Re(ZER(2,0))", "Im(ZER(2,0))", "Re(ZER(2,1))", "Im(ZER(2,1))",
      "Re(ZER(2,2))", "Im(ZER(2,2))", "Re(ZDL(0,0))", "Im(ZDL(0,0))",
      "Re(ZDL(0,1))", "Im(ZDL(0,1))", "Re(ZDL(0,2))", "Im(ZDL(0,2))",
      "Re(ZDL(1,0))", "Im(ZDL(1,0))", "Re(ZDL(1,1))", "Im(ZDL(1,1))",
      "Re(ZDL(1,2))", "Im(ZDL(1,2))", "Re(ZDL(2,0))", "Im(ZDL(2,0))",
      "Re(ZDL(2,1))", "Im(ZDL(2,1))", "Re(ZDL(2,2))", "Im(ZDL(2,2))",
      "Re(ZDR(0,0))", "Im(ZDR(0,0))", "Re(ZDR(0,1))", "Im(ZDR(0,1))",
      "Re(ZDR(0,2))", "Im(ZDR(0,2))", "Re(ZDR(1,0))", "Im(ZDR(1,0))",
      "Re(ZDR(1,1))", "Im(ZDR(1,1))", "Re(ZDR(1,2))", "Im(ZDR(1,2))",
      "Re(ZDR(2,0))", "Im(ZDR(2,0))", "Re(ZDR(2,1))", "Im(ZDR(2,1))",
      "Re(ZDR(2,2))", "Im(ZDR(2,2))", "Re(ZUL(0,0))", "Im(ZUL(0,0))",
      "Re(ZUL(0,1))", "Im(ZUL(0,1))", "Re(ZUL(0,2))", "Im(ZUL(0,2))",
      "Re(ZUL(1,0))", "Im(ZUL(1,0))", "Re(ZUL(1,1))", "Im(ZUL(1,1))",
      "Re(ZUL(1,2))", "Im(ZUL(1,2))", "Re(ZUL(2,0))", "Im(ZUL(2,0))",
      "Re(ZUL(2,1))", "Im(ZUL(2,1))", "Re(ZUL(2,2))", "Im(ZUL(2,2))",
      "Re(ZUR(0,0))", "Im(ZUR(0,0))", "Re(ZUR(0,1))", "Im(ZUR(0,1))",
      "Re(ZUR(0,2))", "Im(ZUR(0,2))", "Re(ZUR(1,0))", "Im(ZUR(1,0))",
      "Re(ZUR(1,1))", "Im(ZUR(1,1))", "Re(ZUR(1,2))", "Im(ZUR(1,2))",
      "Re(ZUR(2,0))", "Im(ZUR(2,0))", "Re(ZUR(2,1))", "Im(ZUR(2,1))",
      "Re(ZUR(2,2))", "Im(ZUR(2,2))", "Re(ZDXL(0,0))", "Im(ZDXL(0,0))",
      "Re(ZDXL(0,1))", "Im(ZDXL(0,1))", "Re(ZDXL(0,2))", "Im(ZDXL(0,2))",
      "Re(ZDXL(1,0))", "Im(ZDXL(1,0))", "Re(ZDXL(1,1))", "Im(ZDXL(1,1))",
      "Re(ZDXL(1,2))", "Im(ZDXL(1,2))", "Re(ZDXL(2,0))", "Im(ZDXL(2,0))",
      "Re(ZDXL(2,1))", "Im(ZDXL(2,1))", "Re(ZDXL(2,2))", "Im(ZDXL(2,2))",
      "Re(ZDXR(0,0))", "Im(ZDXR(0,0))", "Re(ZDXR(0,1))", "Im(ZDXR(0,1))",
      "Re(ZDXR(0,2))", "Im(ZDXR(0,2))", "Re(ZDXR(1,0))", "Im(ZDXR(1,0))",
      "Re(ZDXR(1,1))", "Im(ZDXR(1,1))", "Re(ZDXR(1,2))", "Im(ZDXR(1,2))",
      "Re(ZDXR(2,0))", "Im(ZDXR(2,0))", "Re(ZDXR(2,1))", "Im(ZDXR(2,1))",
      "Re(ZDXR(2,2))", "Im(ZDXR(2,2))", "UHI0(0,0)", "UHI0(0,1)", "UHI0(0,2)",
      "UHI0(0,3)", "UHI0(1,0)", "UHI0(1,1)", "UHI0(1,2)", "UHI0(1,3)", "UHI0(2,0)"
      , "UHI0(2,1)", "UHI0(2,2)", "UHI0(2,3)", "UHI0(3,0)", "UHI0(3,1)",
      "UHI0(3,2)", "UHI0(3,3)", "UHIp(0,0)", "UHIp(0,1)", "UHIp(0,2)", "UHIp(0,3)"
      , "UHIp(1,0)", "UHIp(1,1)", "UHIp(1,2)", "UHIp(1,3)", "UHIp(2,0)",
      "UHIp(2,1)", "UHIp(2,2)", "UHIp(2,3)", "UHIp(3,0)", "UHIp(3,1)", "UHIp(3,2)"
      , "UHIp(3,3)", "Re(ZMI(0,0))", "Im(ZMI(0,0))", "Re(ZMI(0,1))",
      "Im(ZMI(0,1))", "Re(ZMI(1,0))", "Im(ZMI(1,0))", "Re(ZMI(1,1))",
      "Im(ZMI(1,1))", "Re(ZPI(0,0))", "Im(ZPI(0,0))", "Re(ZPI(0,1))",
      "Im(ZPI(0,1))", "Re(ZPI(1,0))", "Im(ZPI(1,0))", "Re(ZPI(1,1))",
      "Im(ZPI(1,1))", "Re(ZNI(0,0))", "Im(ZNI(0,0))", "Re(ZNI(0,1))",
      "Im(ZNI(0,1))", "Re(ZNI(0,2))", "Im(ZNI(0,2))", "Re(ZNI(0,3))",
      "Im(ZNI(0,3))", "Re(ZNI(1,0))", "Im(ZNI(1,0))", "Re(ZNI(1,1))",
      "Im(ZNI(1,1))", "Re(ZNI(1,2))", "Im(ZNI(1,2))", "Re(ZNI(1,3))",
      "Im(ZNI(1,3))", "Re(ZNI(2,0))", "Im(ZNI(2,0))", "Re(ZNI(2,1))",
      "Im(ZNI(2,1))", "Re(ZNI(2,2))", "Im(ZNI(2,2))", "Re(ZNI(2,3))",
      "Im(ZNI(2,3))", "Re(ZNI(3,0))", "Im(ZNI(3,0))", "Re(ZNI(3,1))",
      "Im(ZNI(3,1))", "Re(ZNI(3,2))", "Im(ZNI(3,2))", "Re(ZNI(3,3))",
      "Im(ZNI(3,3))", "ZSSI(0,0)", "ZSSI(0,1)", "ZSSI(1,0)", "ZSSI(1,1)",
      "Re(ZFSI(0,0))", "Im(ZFSI(0,0))", "Re(ZFSI(0,1))", "Im(ZFSI(0,1))",
      "Re(ZFSI(1,0))", "Im(ZFSI(1,0))", "Re(ZFSI(1,1))", "Im(ZFSI(1,1))",
      "UHp0(0,0)", "UHp0(0,1)", "UHp0(1,0)", "UHp0(1,1)", "UHpp(0,0)", "UHpp(0,1)"
      , "UHpp(1,0)", "UHpp(1,1)", "Re(ZNp(0,0))", "Im(ZNp(0,0))", "Re(ZNp(0,1))",
      "Im(ZNp(0,1))", "Re(ZNp(1,0))", "Im(ZNp(1,0))", "Re(ZNp(1,1))",
      "Im(ZNp(1,1))", "ZZ(0,0)", "ZZ(0,1)", "ZZ(0,2)", "ZZ(1,0)", "ZZ(1,1)",
      "ZZ(1,2)", "ZZ(2,0)", "ZZ(2,1)", "ZZ(2,2)"};

   const std::array<std::string, NUMBER_OF_INPUT_PARAMETERS> input_parameter_names
       = {"TanBeta", "m0SqGuess", "m12Guess", "AzeroGuess", "LambdaInput",
      "KappaInput", "MuPrimeInput", "BMuPrimeInput", "vsInput", "Lambda12Input"};

   const std::array<std::string, NUMBER_OF_EXTRA_PARAMETERS> extra_parameter_names
       = {"m0Sq", "m12", "Azero", "MuPrBV"};

   const std::string model_name = "CE6SSM";

int get_pdg_code_for_particle(Particles p)
{
   if (particle_multiplicities[p] > 1) {
      throw OutOfBoundsError(particle_names[p] + " must have a generation index");
   }

   int pdg = 0;
   switch (p) {

   case VG: pdg = 21; break;
   case Glu: pdg = 1000021; break;
   case ChaP: pdg = 1000091; break;
   case VWm: pdg = -24; break;
   case VP: pdg = 22; break;
   case VZ: pdg = 23; break;
   case VZp: pdg = 31; break;

   default: throw OutOfBoundsError("invalid particle " + std::to_string(p));
   }

   return pdg;
}

int get_pdg_code_for_particle(Particles p, int index)
{
   if (particle_multiplicities[p] == 1) {
      throw OutOfBoundsError(particle_names[p] + " does not carry an index");
   }

   std::vector<int> pdg_codes;
   switch (p) {

   case Fv: pdg_codes = {12, 14, 16}; break;
   case Sd: pdg_codes = {1000001, 1000003, 1000005, 2000001, 2000003, 2000005}; break;
   case Sv: pdg_codes = {1000012, 1000014, 1000016}; break;
   case Su: pdg_codes = {1000002, 1000004, 1000006, 2000002, 2000004, 2000006}; break;
   case Se: pdg_codes = {1000011, 1000013, 1000015, 2000011, 2000013, 2000015}; break;
   case SDX: pdg_codes = {1000051, 2000051, 1000052, 2000052, 1000053, 2000053}; break;
   case hh: pdg_codes = {25, 35, 45}; break;
   case Ah: pdg_codes = {0, 0, 36}; break;
   case Hpm: pdg_codes = {0, -37}; break;
   case Chi: pdg_codes = {1000022, 1000023, 1000025, 1000035, 1000045, 1000055}; break;
   case Cha: pdg_codes = {-1000024, -1000037}; break;
   case Fe: pdg_codes = {11, 13, 15}; break;
   case Fd: pdg_codes = {1, 3, 5}; break;
   case Fu: pdg_codes = {2, 4, 6}; break;
   case FDX: pdg_codes = {51, 52, 53}; break;
   case SHI0: pdg_codes = {82, 86, 84, 88}; break;
   case SHIp: pdg_codes = {81, 85, 83, 87}; break;
   case ChaI: pdg_codes = {1000085, 1000086}; break;
   case ChiI: pdg_codes = {1000081, 1000082, 1000083, 1000084}; break;
   case SSI0: pdg_codes = {89, 90}; break;
   case FSI: pdg_codes = {1000089, 1000090}; break;
   case SHp0: pdg_codes = {92, 94}; break;
   case SHpp: pdg_codes = {91, 93}; break;
   case ChiP: pdg_codes = {1000092, 1000094}; break;

   default: throw OutOfBoundsError("invalid particle " + std::to_string(p));
   }

   if (index < 0 || std::abs(index) >= pdg_codes.size()) {
      throw OutOfBoundsError("index " + std::to_string(index) + " out of bounds");
   }

   return pdg_codes[index];
}

std::pair<std::string, std::optional<unsigned int>> get_multiplet_and_index_from_pdg(int pdg)
{
   std::pair<std::string, std::optional<unsigned int>> name;

   switch (pdg) {

   case 21: name = {"VG", {}}; break;
   case 1000021: name = {"Glu", {}}; break;
   case 12: name = {"Fv", 1}; break;
   case 14: name = {"Fv", 2}; break;
   case 16: name = {"Fv", 3}; break;
   case 1000091: name = {"ChaP", {}}; break;
   case 1000001: name = {"Sd", 1}; break;
   case 1000003: name = {"Sd", 2}; break;
   case 1000005: name = {"Sd", 3}; break;
   case 2000001: name = {"Sd", 4}; break;
   case 2000003: name = {"Sd", 5}; break;
   case 2000005: name = {"Sd", 6}; break;
   case 1000012: name = {"Sv", 1}; break;
   case 1000014: name = {"Sv", 2}; break;
   case 1000016: name = {"Sv", 3}; break;
   case 1000002: name = {"Su", 1}; break;
   case 1000004: name = {"Su", 2}; break;
   case 1000006: name = {"Su", 3}; break;
   case 2000002: name = {"Su", 4}; break;
   case 2000004: name = {"Su", 5}; break;
   case 2000006: name = {"Su", 6}; break;
   case 1000011: name = {"Se", 1}; break;
   case 1000013: name = {"Se", 2}; break;
   case 1000015: name = {"Se", 3}; break;
   case 2000011: name = {"Se", 4}; break;
   case 2000013: name = {"Se", 5}; break;
   case 2000015: name = {"Se", 6}; break;
   case 1000051: name = {"SDX", 1}; break;
   case 2000051: name = {"SDX", 2}; break;
   case 1000052: name = {"SDX", 3}; break;
   case 2000052: name = {"SDX", 4}; break;
   case 1000053: name = {"SDX", 5}; break;
   case 2000053: name = {"SDX", 6}; break;
   case 25: name = {"hh", 1}; break;
   case 35: name = {"hh", 2}; break;
   case 45: name = {"hh", 3}; break;
   case 36: name = {"Ah", 3}; break;
   case -37: name = {"Hpm", 2}; break;
   case 1000022: name = {"Chi", 1}; break;
   case 1000023: name = {"Chi", 2}; break;
   case 1000025: name = {"Chi", 3}; break;
   case 1000035: name = {"Chi", 4}; break;
   case 1000045: name = {"Chi", 5}; break;
   case 1000055: name = {"Chi", 6}; break;
   case -1000024: name = {"Cha", 1}; break;
   case -1000037: name = {"Cha", 2}; break;
   case 11: name = {"Fe", 1}; break;
   case 13: name = {"Fe", 2}; break;
   case 15: name = {"Fe", 3}; break;
   case 1: name = {"Fd", 1}; break;
   case 3: name = {"Fd", 2}; break;
   case 5: name = {"Fd", 3}; break;
   case 2: name = {"Fu", 1}; break;
   case 4: name = {"Fu", 2}; break;
   case 6: name = {"Fu", 3}; break;
   case 51: name = {"FDX", 1}; break;
   case 52: name = {"FDX", 2}; break;
   case 53: name = {"FDX", 3}; break;
   case 82: name = {"SHI0", 1}; break;
   case 86: name = {"SHI0", 2}; break;
   case 84: name = {"SHI0", 3}; break;
   case 88: name = {"SHI0", 4}; break;
   case 81: name = {"SHIp", 1}; break;
   case 85: name = {"SHIp", 2}; break;
   case 83: name = {"SHIp", 3}; break;
   case 87: name = {"SHIp", 4}; break;
   case 1000085: name = {"ChaI", 1}; break;
   case 1000086: name = {"ChaI", 2}; break;
   case 1000081: name = {"ChiI", 1}; break;
   case 1000082: name = {"ChiI", 2}; break;
   case 1000083: name = {"ChiI", 3}; break;
   case 1000084: name = {"ChiI", 4}; break;
   case 89: name = {"SSI0", 1}; break;
   case 90: name = {"SSI0", 2}; break;
   case 1000089: name = {"FSI", 1}; break;
   case 1000090: name = {"FSI", 2}; break;
   case 92: name = {"SHp0", 1}; break;
   case 94: name = {"SHp0", 2}; break;
   case 91: name = {"SHpp", 1}; break;
   case 93: name = {"SHpp", 2}; break;
   case 1000092: name = {"ChiP", 1}; break;
   case 1000094: name = {"ChiP", 2}; break;
   case -24: name = {"VWm", {}}; break;
   case 22: name = {"VP", {}}; break;
   case 23: name = {"VZ", {}}; break;
   case 31: name = {"VZp", {}}; break;
   case -12: name = {"barFv", 1}; break;
   case -14: name = {"barFv", 2}; break;
   case -16: name = {"barFv", 3}; break;
   case -1000091: name = {"barChaP", {}}; break;
   case -1000001: name = {"conjSd", 1}; break;
   case -1000003: name = {"conjSd", 2}; break;
   case -1000005: name = {"conjSd", 3}; break;
   case -2000001: name = {"conjSd", 4}; break;
   case -2000003: name = {"conjSd", 5}; break;
   case -2000005: name = {"conjSd", 6}; break;
   case -1000012: name = {"conjSv", 1}; break;
   case -1000014: name = {"conjSv", 2}; break;
   case -1000016: name = {"conjSv", 3}; break;
   case -1000002: name = {"conjSu", 1}; break;
   case -1000004: name = {"conjSu", 2}; break;
   case -1000006: name = {"conjSu", 3}; break;
   case -2000002: name = {"conjSu", 4}; break;
   case -2000004: name = {"conjSu", 5}; break;
   case -2000006: name = {"conjSu", 6}; break;
   case -1000011: name = {"conjSe", 1}; break;
   case -1000013: name = {"conjSe", 2}; break;
   case -1000015: name = {"conjSe", 3}; break;
   case -2000011: name = {"conjSe", 4}; break;
   case -2000013: name = {"conjSe", 5}; break;
   case -2000015: name = {"conjSe", 6}; break;
   case -1000051: name = {"conjSDX", 1}; break;
   case -2000051: name = {"conjSDX", 2}; break;
   case -1000052: name = {"conjSDX", 3}; break;
   case -2000052: name = {"conjSDX", 4}; break;
   case -1000053: name = {"conjSDX", 5}; break;
   case -2000053: name = {"conjSDX", 6}; break;
   case 37: name = {"conjHpm", 2}; break;
   case 1000024: name = {"barCha", 1}; break;
   case 1000037: name = {"barCha", 2}; break;
   case -11: name = {"barFe", 1}; break;
   case -13: name = {"barFe", 2}; break;
   case -15: name = {"barFe", 3}; break;
   case -1: name = {"barFd", 1}; break;
   case -3: name = {"barFd", 2}; break;
   case -5: name = {"barFd", 3}; break;
   case -2: name = {"barFu", 1}; break;
   case -4: name = {"barFu", 2}; break;
   case -6: name = {"barFu", 3}; break;
   case -51: name = {"barFDX", 1}; break;
   case -52: name = {"barFDX", 2}; break;
   case -53: name = {"barFDX", 3}; break;
   case -82: name = {"conjSHI0", 1}; break;
   case -86: name = {"conjSHI0", 2}; break;
   case -84: name = {"conjSHI0", 3}; break;
   case -88: name = {"conjSHI0", 4}; break;
   case -81: name = {"conjSHIp", 1}; break;
   case -85: name = {"conjSHIp", 2}; break;
   case -83: name = {"conjSHIp", 3}; break;
   case -87: name = {"conjSHIp", 4}; break;
   case -1000085: name = {"barChaI", 1}; break;
   case -1000086: name = {"barChaI", 2}; break;
   case -89: name = {"conjSSI0", 1}; break;
   case -90: name = {"conjSSI0", 2}; break;
   case -92: name = {"conjSHp0", 1}; break;
   case -94: name = {"conjSHp0", 2}; break;
   case -91: name = {"conjSHpp", 1}; break;
   case -93: name = {"conjSHpp", 2}; break;
   case 24: name = {"conjVWm", {}}; break;

   default: name = {"", {}};
   }

   return name;
}

std::string get_particle_name_from_pdg(int pdg)
{
   std::pair<std::string, std::optional<unsigned int>> const pair = get_multiplet_and_index_from_pdg(pdg);
   return pair.first + (pair.second.has_value() ? "(" + std::to_string(pair.second.value()) + ")" : "");
}

void print(std::ostream& ostr)
{
   ostr
      << "Model information\n"
      << "=================\n"
      << "Model name:                  " << model_name << '\n'
      << "Is a low-energy model:       "
      << (is_low_energy_model ? "yes" : "no") << '\n'
      << "Is a supersymmetric model:   "
      << (is_supersymmetric_model ? "yes" : "no") << '\n'
      << "Is a FlexibleEFTHiggs model: "
      << (is_FlexibleEFTHiggs ? "yes" : "no") << '\n'
      << "Number of multiplets:        " << NUMBER_OF_PARTICLES << '\n'
      << "Number of parameters:        " << NUMBER_OF_PARAMETERS << '\n'
      ;

   ostr << "\n"
      "Multiplets:                  ";
   for (int i = 0; i < NUMBER_OF_PARTICLES; i++) {
      ostr << particle_names[i]
           << '[' << particle_multiplicities[i] << ']';
      if (i + 1 < NUMBER_OF_PARTICLES)
         ostr << ", ";
   }

   ostr << "\n\n"
      "Parameters:                  ";
   for (int i = 0; i < NUMBER_OF_PARAMETERS; i++) {
      ostr << parameter_names[i];
      if (i + 1 < NUMBER_OF_PARAMETERS)
         ostr << ", ";
   }

   ostr << "\n\n"
      "Input parameters:            ";
   for (int i = 0; i < NUMBER_OF_INPUT_PARAMETERS; i++) {
      ostr << input_parameter_names[i];
      if (i + 1 < NUMBER_OF_INPUT_PARAMETERS)
         ostr << ", ";
   }

   ostr << '\n';
}

} // namespace CE6SSM_info

} // namespace flexiblesusy

