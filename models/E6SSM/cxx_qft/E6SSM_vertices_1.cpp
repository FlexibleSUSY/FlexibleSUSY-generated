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


/**
 * @file cxx_qft/E6SSM_vertices.cpp
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#include "E6SSM_context_base.hpp"
#include "E6SSM_input_parameters.hpp"
#include "E6SSM_vertices.hpp"

#include "concatenate.hpp"
#include "wrappers.hpp"

#define INPUTPARAMETER(p) context.model.get_input_parameters().p
#define MODELPARAMETER(p) context.model.get_##p()
#define DERIVEDPARAMETER(p) context.model.p()
#define PHASE(p) context.model.get_##p()

namespace flexiblesusy::E6SSM_cxx_diagrams::detail {

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::Ah, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::ChaI>::type, E6SSM_cxx_diagrams::fields::ChaI>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto ZMI = MODELPARAMETER(ZMI);
   const auto ZPI = MODELPARAMETER(ZPI);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*SUM(j1,0,1,Conj(ZMI(gt2,j1))*Conj(ZPI(gt1,j1))*Lambda12(j1,j1))*ZA(gt3,2);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*SUM(j1,0,1,Conj(Lambda12(j1,j1))*ZMI(gt1,j1)*ZPI(gt2,j1))*ZA(gt3,2);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::Ah, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Cha>::type, E6SSM_cxx_diagrams::fields::Cha>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto UM = MODELPARAMETER(UM);
   const auto UP = MODELPARAMETER(UP);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,0.7071067811865475)*(g2*Conj(UM(gt2,0))*Conj(UP(gt1,1))*ZA(gt3,1) + Conj(UM(gt2,1))*(g2*Conj(UP(gt1,0))*ZA(gt3,0) - Conj(UP(gt1,1))*Lambdax*ZA(gt3,2)));

   const std::complex<double> right = std::complex<double>(0.,-0.7071067811865475)*(g2*UM(gt1,0)*UP(gt2,1)*ZA(gt3,1) + UM(gt1,1)*(g2*UP(gt2,0)*ZA(gt3,0) - Conj(Lambdax)*UP(gt2,1)*ZA(gt3,2)));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::Ah, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fd>::type, E6SSM_cxx_diagrams::fields::Fd>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto ZDR = MODELPARAMETER(ZDR);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*SUM(j1,0,2,Conj(ZDL(gt2,j1))*Conj(ZDR(gt1,j1))*Yd(j1,j1))*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*SUM(j1,0,2,Conj(Yd(j1,j1))*ZDL(gt1,j1)*ZDR(gt2,j1))*ZA(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::Ah, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::FDX>::type, E6SSM_cxx_diagrams::fields::FDX>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Kappa = MODELPARAMETER(Kappa);
   const auto ZDXL = MODELPARAMETER(ZDXL);
   const auto ZDXR = MODELPARAMETER(ZDXR);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*SUM(j1,0,2,Conj(ZDXL(gt2,j1))*Conj(ZDXR(gt1,j1))*Kappa(j1,j1))*ZA(gt3,2);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*SUM(j1,0,2,Conj(Kappa(j1,j1))*ZDXL(gt1,j1)*ZDXR(gt2,j1))*ZA(gt3,2);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::Ah, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fe>::type, E6SSM_cxx_diagrams::fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*SUM(j1,0,2,Conj(ZEL(gt2,j1))*Conj(ZER(gt1,j1))*Ye(j1,j1))*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*SUM(j1,0,2,Conj(Ye(j1,j1))*ZEL(gt1,j1)*ZER(gt2,j1))*ZA(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::Ah, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fu>::type, E6SSM_cxx_diagrams::fields::Fu>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto ZUR = MODELPARAMETER(ZUR);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*SUM(j1,0,2,Conj(ZUL(gt2,j1))*Conj(ZUR(gt1,j1))*Yu(j1,j1))*ZA(gt3,1);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*SUM(j1,0,2,Conj(Yu(j1,j1))*ZUL(gt1,j1)*ZUR(gt2,j1))*ZA(gt3,1);

   return {left, right};
}

cxx_diagrams::ScalarVertex VertexImpl<E6SSM_cxx_diagrams::fields::Ah, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Hpm>::type, E6SSM_cxx_diagrams::fields::Hpm>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TLambdax = MODELPARAMETER(TLambdax);
   const auto vu = MODELPARAMETER(vu);
   const auto g2 = MODELPARAMETER(g2);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto vd = MODELPARAMETER(vd);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,-0.25)*(vu*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZA(gt1,0)*(ZP(gt2,1)*ZP(gt3,0) - ZP(gt2,0)*ZP(gt3,1)) + vd*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZA(gt1,1)*(ZP(gt2,1)*ZP(gt3,0) - ZP(gt2,0)*ZP(gt3,1)) + 2.8284271247461903*ZA(gt1,2)*(-(Conj(TLambdax)*ZP(gt2,1)*ZP(gt3,0)) + TLambdax*ZP(gt2,0)*ZP(gt3,1)));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<E6SSM_cxx_diagrams::fields::Ah, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Sd>::type, E6SSM_cxx_diagrams::fields::Sd>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TYd = MODELPARAMETER(TYd);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto vs = MODELPARAMETER(vs);
   const auto vu = MODELPARAMETER(vu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,0.5)*(1.4142135623730951*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j1))*ZD(gt3,j1))*ZA(gt1,0) - 1.4142135623730951*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,3 + j1)*TYd(j1,j1))*ZA(gt1,0) + (Lambdax*SUM(j1,0,2,Conj(Yd(j1,j1))*Conj(ZD(gt2,3 + j1))*ZD(gt3,j1)) - Conj(Lambdax)*SUM(j1,0,2,Conj(ZD(gt2,j1))*Yd(j1,j1)*ZD(gt3,3 + j1)))*(vs*ZA(gt1,1) + vu*ZA(gt1,2)));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<E6SSM_cxx_diagrams::fields::Ah, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SDX>::type, E6SSM_cxx_diagrams::fields::SDX>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TKappa = MODELPARAMETER(TKappa);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto vu = MODELPARAMETER(vu);
   const auto vd = MODELPARAMETER(vd);
   const auto Kappa = MODELPARAMETER(Kappa);
   const auto ZDX = MODELPARAMETER(ZDX);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,0.5)*(Lambdax*SUM(j1,0,2,Conj(ZDX(gt2,3 + j1))*Conj(Kappa(j1,j1))*ZDX(gt3,j1))*(vu*ZA(gt1,0) + vd*ZA(gt1,1)) - Conj(Lambdax)*SUM(j1,0,2,Conj(ZDX(gt2,j1))*ZDX(gt3,3 + j1)*Kappa(j1,j1))*(vu*ZA(gt1,0) + vd*ZA(gt1,1)) + 1.4142135623730951*(SUM(j1,0,2,Conj(ZDX(gt2,3 + j1))*Conj(TKappa(j1,j1))*ZDX(gt3,j1)) - SUM(j1,0,2,Conj(ZDX(gt2,j1))*ZDX(gt3,3 + j1)*TKappa(j1,j1)))*ZA(gt1,2));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<E6SSM_cxx_diagrams::fields::Ah, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Se>::type, E6SSM_cxx_diagrams::fields::Se>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TYe = MODELPARAMETER(TYe);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto vs = MODELPARAMETER(vs);
   const auto vu = MODELPARAMETER(vu);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,0.5)*(1.4142135623730951*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j1))*ZE(gt3,j1))*ZA(gt1,0) - 1.4142135623730951*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,3 + j1)*TYe(j1,j1))*ZA(gt1,0) + (Lambdax*SUM(j1,0,2,Conj(Ye(j1,j1))*Conj(ZE(gt2,3 + j1))*ZE(gt3,j1)) - Conj(Lambdax)*SUM(j1,0,2,Conj(ZE(gt2,j1))*Ye(j1,j1)*ZE(gt3,3 + j1)))*(vs*ZA(gt1,1) + vu*ZA(gt1,2)));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<E6SSM_cxx_diagrams::fields::Ah, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHIp>::type, E6SSM_cxx_diagrams::fields::SHIp>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TLambda12 = MODELPARAMETER(TLambda12);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto vu = MODELPARAMETER(vu);
   const auto vd = MODELPARAMETER(vd);
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,0.5)*(Lambdax*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*Conj(Lambda12(j1,j1))*UHIp(gt3,j1))*(vu*ZA(gt1,0) + vd*ZA(gt1,1)) - Conj(Lambdax)*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt3,2 + j1)*Lambda12(j1,j1))*(vu*ZA(gt1,0) + vd*ZA(gt1,1)) + 1.4142135623730951*(SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*Conj(TLambda12(j1,j1))*UHIp(gt3,j1)) - SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt3,2 + j1)*TLambda12(j1,j1)))*ZA(gt1,2));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<E6SSM_cxx_diagrams::fields::Ah, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Su>::type, E6SSM_cxx_diagrams::fields::Su>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TYu = MODELPARAMETER(TYu);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto vs = MODELPARAMETER(vs);
   const auto vd = MODELPARAMETER(vd);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,0.5)*(1.4142135623730951*(SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j1))*ZU(gt3,j1)) - SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,3 + j1)*TYu(j1,j1)))*ZA(gt1,1) + Lambdax*SUM(j1,0,2,Conj(Yu(j1,j1))*Conj(ZU(gt2,3 + j1))*ZU(gt3,j1))*(vs*ZA(gt1,0) + vd*ZA(gt1,2)) - Conj(Lambdax)*SUM(j1,0,2,Conj(ZU(gt2,j1))*Yu(j1,j1)*ZU(gt3,3 + j1))*(vs*ZA(gt1,0) + vd*ZA(gt1,2)));

   return {result};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::Chi, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Se>::type, E6SSM_cxx_diagrams::fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZER = MODELPARAMETER(ZER);

   const std::complex<double> left = 0.5477225575051661*g1*Conj(ZN(gt1,0))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZE(gt3,j1)) + 0.7071067811865475*g2*Conj(ZN(gt1,1))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZE(gt3,j1)) - 0.4472135954999579*gN*Conj(ZN(gt1,5))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZE(gt3,j1)) - Conj(ZN(gt1,2))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*Ye(j1,j1)*ZE(gt3,3 + j1));

   const std::complex<double> right = -(SUM(j1,0,2,Conj(Ye(j1,j1))*ZE(gt3,j1)*ZER(gt2,j1))*ZN(gt1,2)) - 0.22360679774997896*SUM(j1,0,2,ZE(gt3,3 + j1)*ZER(gt2,j1))*(4.898979485566356*g1*ZN(gt1,0) + gN*ZN(gt1,5));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::Fe, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fe>::type, E6SSM_cxx_diagrams::fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.7745966692414834*g1*Cos(ThetaW)*KroneckerDelta(gt1,gt2);

   const std::complex<double> right = -0.5*KroneckerDelta(gt1,gt2)*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::Fe, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fe>::type, E6SSM_cxx_diagrams::fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.05*KroneckerDelta(gt1,gt2)*(15.491933384829668*g1*Cos(ThetaWp)*Sin(ThetaW) - 3.1622776601683795*gN*Sin(ThetaWp));

   const std::complex<double> right = -0.1*KroneckerDelta(gt1,gt2)*(5*g2*Cos(ThetaW)*Cos(ThetaWp) - 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) - 3.1622776601683795*gN*Sin(ThetaWp));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::hh, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::ChaI>::type, E6SSM_cxx_diagrams::fields::ChaI>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto ZMI = MODELPARAMETER(ZMI);
   const auto ZPI = MODELPARAMETER(ZPI);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*SUM(j1,0,1,Conj(ZMI(gt2,j1))*Conj(ZPI(gt1,j1))*Lambda12(j1,j1))*ZH(gt3,2);

   const std::complex<double> right = -0.7071067811865475*SUM(j1,0,1,Conj(Lambda12(j1,j1))*ZMI(gt1,j1)*ZPI(gt2,j1))*ZH(gt3,2);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::hh, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Cha>::type, E6SSM_cxx_diagrams::fields::Cha>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto UM = MODELPARAMETER(UM);
   const auto UP = MODELPARAMETER(UP);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*(g2*Conj(UM(gt2,0))*Conj(UP(gt1,1))*ZH(gt3,1) + Conj(UM(gt2,1))*(g2*Conj(UP(gt1,0))*ZH(gt3,0) + Conj(UP(gt1,1))*Lambdax*ZH(gt3,2)));

   const std::complex<double> right = -0.7071067811865475*(g2*UM(gt1,0)*UP(gt2,1)*ZH(gt3,1) + UM(gt1,1)*(g2*UP(gt2,0)*ZH(gt3,0) + Conj(Lambdax)*UP(gt2,1)*ZH(gt3,2)));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::hh, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fd>::type, E6SSM_cxx_diagrams::fields::Fd>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto ZDR = MODELPARAMETER(ZDR);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*SUM(j1,0,2,Conj(ZDL(gt2,j1))*Conj(ZDR(gt1,j1))*Yd(j1,j1))*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*SUM(j1,0,2,Conj(Yd(j1,j1))*ZDL(gt1,j1)*ZDR(gt2,j1))*ZH(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::hh, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::FDX>::type, E6SSM_cxx_diagrams::fields::FDX>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Kappa = MODELPARAMETER(Kappa);
   const auto ZDXL = MODELPARAMETER(ZDXL);
   const auto ZDXR = MODELPARAMETER(ZDXR);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*SUM(j1,0,2,Conj(ZDXL(gt2,j1))*Conj(ZDXR(gt1,j1))*Kappa(j1,j1))*ZH(gt3,2);

   const std::complex<double> right = -0.7071067811865475*SUM(j1,0,2,Conj(Kappa(j1,j1))*ZDXL(gt1,j1)*ZDXR(gt2,j1))*ZH(gt3,2);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::hh, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fe>::type, E6SSM_cxx_diagrams::fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*SUM(j1,0,2,Conj(ZEL(gt2,j1))*Conj(ZER(gt1,j1))*Ye(j1,j1))*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*SUM(j1,0,2,Conj(Ye(j1,j1))*ZEL(gt1,j1)*ZER(gt2,j1))*ZH(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::hh, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fu>::type, E6SSM_cxx_diagrams::fields::Fu>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto ZUR = MODELPARAMETER(ZUR);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*SUM(j1,0,2,Conj(ZUL(gt2,j1))*Conj(ZUR(gt1,j1))*Yu(j1,j1))*ZH(gt3,1);

   const std::complex<double> right = -0.7071067811865475*SUM(j1,0,2,Conj(Yu(j1,j1))*ZUL(gt1,j1)*ZUR(gt2,j1))*ZH(gt3,1);

   return {left, right};
}

cxx_diagrams::ScalarVertex VertexImpl<E6SSM_cxx_diagrams::fields::hh, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Hpm>::type, E6SSM_cxx_diagrams::fields::Hpm>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TLambdax = MODELPARAMETER(TLambdax);
   const auto vs = MODELPARAMETER(vs);
   const auto gN = MODELPARAMETER(gN);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto vd = MODELPARAMETER(vd);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.025*(2*ZH(gt1,1)*(ZP(gt2,0)*(vu*(3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZP(gt3,0) - 5*vd*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gt3,1)) - ZP(gt2,1)*(5*vd*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gt3,0) + vu*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZP(gt3,1))) - ZH(gt1,0)*(ZP(gt2,0)*(vd*(6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*ZP(gt3,0) + 10*vu*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gt3,1)) + 2*ZP(gt2,1)*(5*vu*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gt3,0) + vd*(-3*Sqr(g1) + 5*Sqr(g2) + 3*Sqr(gN))*ZP(gt3,1))) + 5*ZH(gt1,2)*(2*ZP(gt2,1)*(-2.8284271247461903*Conj(TLambdax)*ZP(gt3,0) + vs*(-4*AbsSqr(Lambdax) + Sqr(gN))*ZP(gt3,1)) + ZP(gt2,0)*(vs*(-8*AbsSqr(Lambdax) + 3*Sqr(gN))*ZP(gt3,0) - 5.656854249492381*TLambdax*ZP(gt3,1))));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<E6SSM_cxx_diagrams::fields::hh, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Sd>::type, E6SSM_cxx_diagrams::fields::Sd>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TYd = MODELPARAMETER(TYd);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto vs = MODELPARAMETER(vs);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.025*(-40*vd*SUM(j1,0,2,AbsSqr(Yd(j1,j1))*Conj(ZD(gt2,j1))*ZD(gt3,j1))*ZH(gt1,0) - 28.284271247461902*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j1))*ZD(gt3,j1))*ZH(gt1,0) + 4*vd*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*ZH(gt1,0) + 6*vd*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*ZH(gt1,0) - 40*vd*SUM(j1,0,2,AbsSqr(Yd(j1,j1))*Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*ZH(gt1,0) - 28.284271247461902*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,3 + j1)*TYd(j1,j1))*ZH(gt1,0) + 20*vs*Lambdax*SUM(j1,0,2,Conj(Yd(j1,j1))*Conj(ZD(gt2,3 + j1))*ZD(gt3,j1))*ZH(gt1,1) - 4*vu*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*ZH(gt1,1) + 4*vu*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*ZH(gt1,1) + 20*vs*Conj(Lambdax)*SUM(j1,0,2,Conj(ZD(gt2,j1))*Yd(j1,j1)*ZD(gt3,3 + j1))*ZH(gt1,1) + 20*vu*Lambdax*SUM(j1,0,2,Conj(Yd(j1,j1))*Conj(ZD(gt2,3 + j1))*ZD(gt3,j1))*ZH(gt1,2) - 10*vs*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*ZH(gt1,2) + 20*vu*Conj(Lambdax)*SUM(j1,0,2,Conj(ZD(gt2,j1))*Yd(j1,j1)*ZD(gt3,3 + j1))*ZH(gt1,2) + SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*(vd*(2*Sqr(g1) + 10*Sqr(g2) + 3*Sqr(gN))*ZH(gt1,0) - 2*vu*(Sqr(g1) + 5*Sqr(g2) - Sqr(gN))*ZH(gt1,1) - 5*vs*Sqr(gN)*ZH(gt1,2)));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<E6SSM_cxx_diagrams::fields::hh, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SDX>::type, E6SSM_cxx_diagrams::fields::SDX>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TKappa = MODELPARAMETER(TKappa);
   const auto g1 = MODELPARAMETER(g1);
   const auto vd = MODELPARAMETER(vd);
   const auto gN = MODELPARAMETER(gN);
   const auto vu = MODELPARAMETER(vu);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto vs = MODELPARAMETER(vs);
   const auto Kappa = MODELPARAMETER(Kappa);
   const auto ZDX = MODELPARAMETER(ZDX);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.025*(4*vd*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gt2,3 + j1))*ZDX(gt3,3 + j1))*ZH(gt1,0) - 9*vd*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gt2,3 + j1))*ZDX(gt3,3 + j1))*ZH(gt1,0) + 20*vu*Conj(Lambdax)*SUM(j1,0,2,Conj(ZDX(gt2,j1))*ZDX(gt3,3 + j1)*Kappa(j1,j1))*ZH(gt1,0) - 4*vu*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gt2,3 + j1))*ZDX(gt3,3 + j1))*ZH(gt1,1) - 6*vu*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gt2,3 + j1))*ZDX(gt3,3 + j1))*ZH(gt1,1) + 20*vd*Conj(Lambdax)*SUM(j1,0,2,Conj(ZDX(gt2,j1))*ZDX(gt3,3 + j1)*Kappa(j1,j1))*ZH(gt1,1) + 20*Lambdax*SUM(j1,0,2,Conj(ZDX(gt2,3 + j1))*Conj(Kappa(j1,j1))*ZDX(gt3,j1))*(vu*ZH(gt1,0) + vd*ZH(gt1,1)) - 40*vs*SUM(j1,0,2,AbsSqr(Kappa(j1,j1))*Conj(ZDX(gt2,j1))*ZDX(gt3,j1))*ZH(gt1,2) - 28.284271247461902*SUM(j1,0,2,Conj(ZDX(gt2,3 + j1))*Conj(TKappa(j1,j1))*ZDX(gt3,j1))*ZH(gt1,2) + 15*vs*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gt2,3 + j1))*ZDX(gt3,3 + j1))*ZH(gt1,2) - 40*vs*SUM(j1,0,2,AbsSqr(Kappa(j1,j1))*Conj(ZDX(gt2,3 + j1))*ZDX(gt3,3 + j1))*ZH(gt1,2) - 28.284271247461902*SUM(j1,0,2,Conj(ZDX(gt2,j1))*ZDX(gt3,3 + j1)*TKappa(j1,j1))*ZH(gt1,2) - 2*SUM(j1,0,2,Conj(ZDX(gt2,j1))*ZDX(gt3,j1))*(vd*(2*Sqr(g1) + 3*Sqr(gN))*ZH(gt1,0) + 2*vu*(-Sqr(g1) + Sqr(gN))*ZH(gt1,1) - 5*vs*Sqr(gN)*ZH(gt1,2)));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<E6SSM_cxx_diagrams::fields::hh, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Se>::type, E6SSM_cxx_diagrams::fields::Se>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TYe = MODELPARAMETER(TYe);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto vs = MODELPARAMETER(vs);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto vu = MODELPARAMETER(vu);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.025*(-40*vd*SUM(j1,0,2,AbsSqr(Ye(j1,j1))*Conj(ZE(gt2,j1))*ZE(gt3,j1))*ZH(gt1,0) - 28.284271247461902*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j1))*ZE(gt3,j1))*ZH(gt1,0) + 12*vd*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*ZH(gt1,0) + 3*vd*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*ZH(gt1,0) - 40*vd*SUM(j1,0,2,AbsSqr(Ye(j1,j1))*Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*ZH(gt1,0) - 28.284271247461902*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,3 + j1)*TYe(j1,j1))*ZH(gt1,0) + 20*vs*Lambdax*SUM(j1,0,2,Conj(Ye(j1,j1))*Conj(ZE(gt2,3 + j1))*ZE(gt3,j1))*ZH(gt1,1) - 12*vu*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*ZH(gt1,1) + 2*vu*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*ZH(gt1,1) + 20*vs*Conj(Lambdax)*SUM(j1,0,2,Conj(ZE(gt2,j1))*Ye(j1,j1)*ZE(gt3,3 + j1))*ZH(gt1,1) + 20*vu*Lambdax*SUM(j1,0,2,Conj(Ye(j1,j1))*Conj(ZE(gt2,3 + j1))*ZE(gt3,j1))*ZH(gt1,2) - 5*vs*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*ZH(gt1,2) + 20*vu*Conj(Lambdax)*SUM(j1,0,2,Conj(ZE(gt2,j1))*Ye(j1,j1)*ZE(gt3,3 + j1))*ZH(gt1,2) - 2*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*(vd*(3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,0) + vu*(-3*Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN))*ZH(gt1,1) + 5*vs*Sqr(gN)*ZH(gt1,2)));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<E6SSM_cxx_diagrams::fields::hh, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHIp>::type, E6SSM_cxx_diagrams::fields::SHIp>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TLambda12 = MODELPARAMETER(TLambda12);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto vs = MODELPARAMETER(vs);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.025*(20*Lambdax*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*Conj(Lambda12(j1,j1))*UHIp(gt3,j1))*(vu*ZH(gt1,0) + vd*ZH(gt1,1)) + 2*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*UHIp(gt3,2 + j1))*(vd*(3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,0) + vu*(-3*Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN))*ZH(gt1,1) + 5*vs*Sqr(gN)*ZH(gt1,2)) + SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt3,j1))*(vd*(-6*Sqr(g1) + 10*Sqr(g2) - 9*Sqr(gN))*ZH(gt1,0) + 2*vu*(3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,1) + 15*vs*Sqr(gN)*ZH(gt1,2)) + 20*(Conj(Lambdax)*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt3,2 + j1)*Lambda12(j1,j1))*(vu*ZH(gt1,0) + vd*ZH(gt1,1)) - (2*vs*SUM(j1,0,1,AbsSqr(Lambda12(j1,j1))*Conj(UHIp(gt2,j1))*UHIp(gt3,j1)) + 1.4142135623730951*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*Conj(TLambda12(j1,j1))*UHIp(gt3,j1)) + 2*vs*SUM(j1,0,1,AbsSqr(Lambda12(j1,j1))*Conj(UHIp(gt2,2 + j1))*UHIp(gt3,2 + j1)) + 1.4142135623730951*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt3,2 + j1)*TLambda12(j1,j1)))*ZH(gt1,2)));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<E6SSM_cxx_diagrams::fields::hh, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHpp>::type, E6SSM_cxx_diagrams::fields::SHpp>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto vs = MODELPARAMETER(vs);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = -0.05*(Conj(UHpp(gt2,0))*UHpp(gt3,0) - Conj(UHpp(gt2,1))*UHpp(gt3,1))*(vd*(3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,0) + vu*(-3*Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN))*ZH(gt1,1) + 5*vs*Sqr(gN)*ZH(gt1,2));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<E6SSM_cxx_diagrams::fields::hh, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Su>::type, E6SSM_cxx_diagrams::fields::Su>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TYu = MODELPARAMETER(TYu);
   const auto g1 = MODELPARAMETER(g1);
   const auto vd = MODELPARAMETER(vd);
   const auto gN = MODELPARAMETER(gN);
   const auto vs = MODELPARAMETER(vs);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto vu = MODELPARAMETER(vu);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.025*(-8*vd*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*ZH(gt1,0) + 3*vd*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*ZH(gt1,0) + 20*vs*Conj(Lambdax)*SUM(j1,0,2,Conj(ZU(gt2,j1))*Yu(j1,j1)*ZU(gt3,3 + j1))*ZH(gt1,0) - 40*vu*SUM(j1,0,2,AbsSqr(Yu(j1,j1))*Conj(ZU(gt2,j1))*ZU(gt3,j1))*ZH(gt1,1) - 28.284271247461902*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j1))*ZU(gt3,j1))*ZH(gt1,1) + 8*vu*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*ZH(gt1,1) + 2*vu*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*ZH(gt1,1) - 40*vu*SUM(j1,0,2,AbsSqr(Yu(j1,j1))*Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*ZH(gt1,1) - 28.284271247461902*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,3 + j1)*TYu(j1,j1))*ZH(gt1,1) - 5*vs*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*ZH(gt1,2) + 20*vd*Conj(Lambdax)*SUM(j1,0,2,Conj(ZU(gt2,j1))*Yu(j1,j1)*ZU(gt3,3 + j1))*ZH(gt1,2) + 20*Lambdax*SUM(j1,0,2,Conj(Yu(j1,j1))*Conj(ZU(gt2,3 + j1))*ZU(gt3,j1))*(vs*ZH(gt1,0) + vd*ZH(gt1,2)) + SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*(vd*(2*Sqr(g1) - 10*Sqr(g2) + 3*Sqr(gN))*ZH(gt1,0) + 2*vu*(-Sqr(g1) + 5*Sqr(g2) + Sqr(gN))*ZH(gt1,1) - 5*vs*Sqr(gN)*ZH(gt1,2)));

   return {result};
}

cxx_diagrams::InverseMetricVertex VertexImpl<E6SSM_cxx_diagrams::fields::hh, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::VWm>::type, E6SSM_cxx_diagrams::fields::VWm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.5*Sqr(g2)*(vd*ZH(gt1,0) + vu*ZH(gt1,1));

   return {result};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::VP, E6SSM_cxx_diagrams::fields::ChaI, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::ChaI>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.5*KroneckerDelta(gt1,gt2)*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   const std::complex<double> right = -0.5*KroneckerDelta(gt1,gt2)*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::VP, E6SSM_cxx_diagrams::fields::Cha, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Cha>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UP = MODELPARAMETER(UP);
   const auto UM = MODELPARAMETER(UM);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -(g2*Conj(UP(gt1,0))*Sin(ThetaW)*UP(gt2,0)) - 0.5*Conj(UP(gt1,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UP(gt2,1);

   const std::complex<double> right = -(g2*Conj(UM(gt2,0))*Sin(ThetaW)*UM(gt1,0)) - 0.5*Conj(UM(gt2,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UM(gt1,1);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::VP, E6SSM_cxx_diagrams::fields::Fd, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fd>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.2581988897471611*g1*Cos(ThetaW)*KroneckerDelta(gt1,gt2);

   const std::complex<double> right = 0.16666666666666666*KroneckerDelta(gt1,gt2)*(0.7745966692414834*g1*Cos(ThetaW) - 3*g2*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::VP, E6SSM_cxx_diagrams::fields::FDX, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::FDX>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.2581988897471611*g1*Cos(ThetaW)*KroneckerDelta(gt1,gt2);

   const std::complex<double> right = -0.2581988897471611*g1*Cos(ThetaW)*KroneckerDelta(gt1,gt2);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::VP, E6SSM_cxx_diagrams::fields::Fe, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fe>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.7745966692414834*g1*Cos(ThetaW)*KroneckerDelta(gt1,gt2);

   const std::complex<double> right = -0.5*KroneckerDelta(gt1,gt2)*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::VP, E6SSM_cxx_diagrams::fields::Fu, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fu>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.5163977794943222*g1*Cos(ThetaW)*KroneckerDelta(gt1,gt2);

   const std::complex<double> right = 0.16666666666666666*KroneckerDelta(gt1,gt2)*(0.7745966692414834*g1*Cos(ThetaW) + 3*g2*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<E6SSM_cxx_diagrams::fields::VP, E6SSM_cxx_diagrams::fields::Hpm, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Hpm>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*(ZP(gt1,0)*ZP(gt2,0) + ZP(gt1,1)*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<E6SSM_cxx_diagrams::fields::VP, E6SSM_cxx_diagrams::fields::Sd, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Sd>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.03333333333333333*(3.872983346207417*g1*Cos(ThetaW) - 15*g2*Sin(ThetaW))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1)) + 0.2581988897471611*g1*Cos(ThetaW)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1));

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<E6SSM_cxx_diagrams::fields::VP, E6SSM_cxx_diagrams::fields::SDX, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SDX>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.2581988897471611*g1*Cos(ThetaW)*KroneckerDelta(gt1,gt2);

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<E6SSM_cxx_diagrams::fields::VP, E6SSM_cxx_diagrams::fields::Se, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Se>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) + 0.7745966692414834*g1*Cos(ThetaW)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1));

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<E6SSM_cxx_diagrams::fields::VP, E6SSM_cxx_diagrams::fields::SHIp, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHIp>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*KroneckerDelta(gt1,gt2)*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<E6SSM_cxx_diagrams::fields::VP, E6SSM_cxx_diagrams::fields::SHpp, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHpp>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*(Conj(UHpp(gt1,0))*UHpp(gt2,0) + Conj(UHpp(gt1,1))*UHpp(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<E6SSM_cxx_diagrams::fields::VP, E6SSM_cxx_diagrams::fields::Su, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Su>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.16666666666666666*(0.7745966692414834*g1*Cos(ThetaW) + 3*g2*Sin(ThetaW))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1)) - 0.5163977794943222*g1*Cos(ThetaW)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1));

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::TripleVectorVertex VertexImpl<E6SSM_cxx_diagrams::fields::VP, E6SSM_cxx_diagrams::fields::VWm, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::VWm>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Sin(ThetaW);

   return {result, cxx_diagrams::TripleVectorVertex::odd_permutation{}};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::VP, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fe>::type, E6SSM_cxx_diagrams::fields::Fe>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.5*KroneckerDelta(gt1,gt2)*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   const std::complex<double> right = 0.7745966692414834*g1*Cos(ThetaW)*KroneckerDelta(gt1,gt2);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::VZ, E6SSM_cxx_diagrams::fields::ChaI, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::ChaI>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = -0.1*KroneckerDelta(gt1,gt2)*(5*g2*Cos(ThetaW)*Cos(ThetaWp) - 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) - 3.1622776601683795*gN*Sin(ThetaWp));

   const std::complex<double> right = -0.05*KroneckerDelta(gt1,gt2)*(10*g2*Cos(ThetaW)*Cos(ThetaWp) - 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::VZ, E6SSM_cxx_diagrams::fields::Cha, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Cha>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto UP = MODELPARAMETER(UP);
   const auto UM = MODELPARAMETER(UM);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = -(g2*Conj(UP(gt1,0))*Cos(ThetaW)*Cos(ThetaWp)*UP(gt2,0)) + 0.1*Conj(UP(gt1,1))*(-5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*UP(gt2,1);

   const std::complex<double> right = -(g2*Conj(UM(gt2,0))*Cos(ThetaW)*Cos(ThetaWp)*UM(gt1,0)) - 0.05*Conj(UM(gt2,1))*(10*g2*Cos(ThetaW)*Cos(ThetaWp) - 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp))*UM(gt1,1);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::VZ, E6SSM_cxx_diagrams::fields::Fd, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fd>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = KroneckerDelta(gt1,gt2)*(0.2581988897471611*g1*Cos(ThetaWp)*Sin(ThetaW) - 0.31622776601683794*gN*Sin(ThetaWp));

   const std::complex<double> right = -0.016666666666666666*KroneckerDelta(gt1,gt2)*(30*g2*Cos(ThetaW)*Cos(ThetaWp) + 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::VZ, E6SSM_cxx_diagrams::fields::FDX, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::FDX>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.016666666666666666*KroneckerDelta(gt1,gt2)*(15.491933384829668*g1*Cos(ThetaWp)*Sin(ThetaW) + 28.460498941515414*gN*Sin(ThetaWp));

   const std::complex<double> right = KroneckerDelta(gt1,gt2)*(0.2581988897471611*g1*Cos(ThetaWp)*Sin(ThetaW) - 0.31622776601683794*gN*Sin(ThetaWp));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::VZ, E6SSM_cxx_diagrams::fields::Fe, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fe>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.05*KroneckerDelta(gt1,gt2)*(15.491933384829668*g1*Cos(ThetaWp)*Sin(ThetaW) - 3.1622776601683795*gN*Sin(ThetaWp));

   const std::complex<double> right = -0.1*KroneckerDelta(gt1,gt2)*(5*g2*Cos(ThetaW)*Cos(ThetaWp) - 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) - 3.1622776601683795*gN*Sin(ThetaWp));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::VZ, E6SSM_cxx_diagrams::fields::Fu, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fu>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.016666666666666666*KroneckerDelta(gt1,gt2)*(30.983866769659336*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp));

   const std::complex<double> right = 0.016666666666666666*KroneckerDelta(gt1,gt2)*(30*g2*Cos(ThetaW)*Cos(ThetaWp) - 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::VZp, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fe>::type, E6SSM_cxx_diagrams::fields::Fe>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.1*KroneckerDelta(gt1,gt2)*(3.1622776601683795*gN*Cos(ThetaWp) + (5*g2*Cos(ThetaW) - 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp));

   const std::complex<double> right = 0.05*KroneckerDelta(gt1,gt2)*(3.1622776601683795*gN*Cos(ThetaWp) + 15.491933384829668*g1*Sin(ThetaW)*Sin(ThetaWp));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<E6SSM_cxx_diagrams::fields::VZ, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fe>::type, E6SSM_cxx_diagrams::fields::Fe>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = 0.1*KroneckerDelta(gt1,gt2)*(5*g2*Cos(ThetaW)*Cos(ThetaWp) - 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) - 3.1622776601683795*gN*Sin(ThetaWp));

   const std::complex<double> right = -0.05*KroneckerDelta(gt1,gt2)*(15.491933384829668*g1*Cos(ThetaWp)*Sin(ThetaW) - 3.1622776601683795*gN*Sin(ThetaWp));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Cha>::type, E6SSM_cxx_diagrams::fields::Cha, E6SSM_cxx_diagrams::fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UM = MODELPARAMETER(UM);
   const auto UP = MODELPARAMETER(UP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = g2*Conj(UM(gt2,0))*Sin(ThetaW)*UM(gt1,0) + 0.5*Conj(UM(gt2,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UM(gt1,1);

   const std::complex<double> right = g2*Conj(UP(gt1,0))*Sin(ThetaW)*UP(gt2,0) + 0.5*Conj(UP(gt1,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UP(gt2,1);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fe>::type, E6SSM_cxx_diagrams::fields::Cha, E6SSM_cxx_diagrams::fields::Sv>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto UM = MODELPARAMETER(UM);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto UP = MODELPARAMETER(UP);

   const std::complex<double> left = Conj(UM(gt2,1))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Conj(ZV(gt3,j1))*Ye(j1,j1));

   const std::complex<double> right = -(g2*SUM(j1,0,2,Conj(ZV(gt3,j1))*ZEL(gt1,j1))*UP(gt2,0));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fe>::type, E6SSM_cxx_diagrams::fields::Fe, E6SSM_cxx_diagrams::fields::Ah>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*SUM(j1,0,2,Conj(ZEL(gt2,j1))*Conj(ZER(gt1,j1))*Ye(j1,j1))*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*SUM(j1,0,2,Conj(Ye(j1,j1))*ZEL(gt1,j1)*ZER(gt2,j1))*ZA(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fe>::type, E6SSM_cxx_diagrams::fields::Fe, E6SSM_cxx_diagrams::fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*SUM(j1,0,2,Conj(ZEL(gt2,j1))*Conj(ZER(gt1,j1))*Ye(j1,j1))*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*SUM(j1,0,2,Conj(Ye(j1,j1))*ZEL(gt1,j1)*ZER(gt2,j1))*ZH(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fe>::type, E6SSM_cxx_diagrams::fields::Fe, E6SSM_cxx_diagrams::fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.5*KroneckerDelta(gt1,gt2)*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   const std::complex<double> right = 0.7745966692414834*g1*Cos(ThetaW)*KroneckerDelta(gt1,gt2);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fe>::type, E6SSM_cxx_diagrams::fields::Fe, E6SSM_cxx_diagrams::fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.1*KroneckerDelta(gt1,gt2)*(3.1622776601683795*gN*Cos(ThetaWp) + (5*g2*Cos(ThetaW) - 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp));

   const std::complex<double> right = 0.05*KroneckerDelta(gt1,gt2)*(3.1622776601683795*gN*Cos(ThetaWp) + 15.491933384829668*g1*Sin(ThetaW)*Sin(ThetaWp));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fe>::type, E6SSM_cxx_diagrams::fields::Fe, E6SSM_cxx_diagrams::fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = 0.1*KroneckerDelta(gt1,gt2)*(5*g2*Cos(ThetaW)*Cos(ThetaWp) - 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) - 3.1622776601683795*gN*Sin(ThetaWp));

   const std::complex<double> right = -0.05*KroneckerDelta(gt1,gt2)*(15.491933384829668*g1*Cos(ThetaWp)*Sin(ThetaW) - 3.1622776601683795*gN*Sin(ThetaWp));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fe>::type, E6SSM_cxx_diagrams::fields::Hpm, E6SSM_cxx_diagrams::fields::Fv>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> left = IF(gt2 < 3,Conj(ZER(gt1,gt2))*Ye(gt2,gt2)*ZP(gt3,0),0);

   const std::complex<double> right = 0;

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fe>::type, E6SSM_cxx_diagrams::fields::Se, E6SSM_cxx_diagrams::fields::Chi>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZEL = MODELPARAMETER(ZEL);

   const std::complex<double> left = -1.0954451150103321*g1*Conj(ZN(gt2,0))*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*Conj(ZER(gt1,j1))) - 0.22360679774997896*gN*Conj(ZN(gt2,5))*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*Conj(ZER(gt1,j1))) - Conj(ZN(gt2,2))*SUM(j1,0,2,Conj(ZE(gt3,j1))*Conj(ZER(gt1,j1))*Ye(j1,j1));

   const std::complex<double> right = -(SUM(j1,0,2,Conj(Ye(j1,j1))*Conj(ZE(gt3,3 + j1))*ZEL(gt1,j1))*ZN(gt2,2)) + SUM(j1,0,2,Conj(ZE(gt3,j1))*ZEL(gt1,j1))*(0.5477225575051661*g1*ZN(gt2,0) + 0.7071067811865475*g2*ZN(gt2,1) - 0.4472135954999579*gN*ZN(gt2,5));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fe>::type, E6SSM_cxx_diagrams::fields::VWm, E6SSM_cxx_diagrams::fields::Fv>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZEL = MODELPARAMETER(ZEL);

   const std::complex<double> left = IF(gt2 < 3,-0.7071067811865475*g2*ZEL(gt1,gt2),0);

   const std::complex<double> right = 0;

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fv>::type, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Hpm>::type, E6SSM_cxx_diagrams::fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> left = 0;

   const std::complex<double> right = IF(gt1 < 3,Conj(Ye(gt1,gt1))*ZER(gt2,gt1)*ZP(gt3,0),0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fv>::type, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::VWm>::type, E6SSM_cxx_diagrams::fields::Fe>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZEL = MODELPARAMETER(ZEL);

   const std::complex<double> left = IF(gt1 < 3,-0.7071067811865475*g2*Conj(ZEL(gt2,gt1)),0);

   const std::complex<double> right = 0;

   return {left, right};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Hpm>::type, E6SSM_cxx_diagrams::fields::Hpm, E6SSM_cxx_diagrams::fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 0;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*(ZP(gt1,0)*ZP(gt2,0) + ZP(gt1,1)*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::InverseMetricVertex VertexImpl<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Hpm>::type, E6SSM_cxx_diagrams::fields::VWm, E6SSM_cxx_diagrams::fields::VP>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.3872983346207417*g1*g2*Cos(ThetaW)*(vd*ZP(gt1,0) - vu*ZP(gt1,1));

   return {result};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Se>::type, E6SSM_cxx_diagrams::fields::Se, E6SSM_cxx_diagrams::fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 0;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) + 0.7745966692414834*g1*Cos(ThetaW)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1));

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::ChiralVertex VertexImpl<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Sv>::type, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Cha>::type, E6SSM_cxx_diagrams::fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto UP = MODELPARAMETER(UP);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZER = MODELPARAMETER(ZER);
   const auto UM = MODELPARAMETER(UM);

   const std::complex<double> left = -(g2*Conj(UP(gt1,0))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZV(gt3,j1)));

   const std::complex<double> right = SUM(j1,0,2,Conj(Ye(j1,j1))*ZER(gt2,j1)*ZV(gt3,j1))*UM(gt1,1);

   return {left, right};
}

cxx_diagrams::InverseMetricVertex VertexImpl<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::VWm>::type, E6SSM_cxx_diagrams::fields::Hpm, E6SSM_cxx_diagrams::fields::VP>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.3872983346207417*g1*g2*Cos(ThetaW)*(vd*ZP(gt1,0) - vu*ZP(gt1,1));

   return {result};
}

cxx_diagrams::TripleVectorVertex VertexImpl<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::VWm>::type, E6SSM_cxx_diagrams::fields::VWm, E6SSM_cxx_diagrams::fields::VP>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Sin(ThetaW));

   return {result, cxx_diagrams::TripleVectorVertex::odd_permutation{}};
}

} // namespace flexiblesusy::E6SSM_cxx_diagrams::detail
