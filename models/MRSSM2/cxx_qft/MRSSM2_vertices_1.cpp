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
 * @file cxx_qft/MRSSM2_vertices.cpp
 *
 * This file was generated with FlexibleSUSY 2.6.2 and SARAH 4.14.5 .
 */

#include "MRSSM2_context_base.hpp"
#include "MRSSM2_input_parameters.hpp"
#include "MRSSM2_vertices.hpp"

#include "concatenate.hpp"
#include "wrappers.hpp"

#define INPUTPARAMETER(p) context.model.get_input_parameters().p
#define MODELPARAMETER(p) context.model.get_##p()
#define DERIVEDPARAMETER(p) context.model.p()
#define PHASE(p) context.model.get_##p()

namespace flexiblesusy {
namespace MRSSM2_cxx_diagrams {
namespace detail {

ChiralVertex VertexImpl<fields::Fe, typename fields::bar<fields::Fe>::type, fields::VP>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fe, fields::Ah>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,j2)))*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZEL(gt1,j2))*ZA(gt3,0);

   return {left, right};
}

ChiralVertex VertexImpl<fields::VP, fields::Cha1, typename fields::bar<fields::Cha1>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = g2*Conj(UM1(gt1,0))*Sin(ThetaW)*UM1(gt2,0) + 0.5*Conj(UM1(gt1,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UM1(gt2,1);

   const std::complex<double> right = g2*Conj(UP1(gt2,0))*Sin(ThetaW)*UP1(gt1,0) + 0.5*Conj(UP1(gt2,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UP1(gt1,1);

   return {left, right};
}

ChiralVertex VertexImpl<fields::Ah, typename fields::bar<fields::Cha1>::type, fields::Cha1>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0,0.5)*(Conj(UM1(gt1,0))*(-1.4142135623730951*LamTD*Conj(UP1(gt2,1))*ZA(gt3,0) + 2*g2*Conj(UP1(gt2,0))*ZA(gt3,3)) + Conj(UM1(gt1,1))*(1.4142135623730951*g2*Conj(UP1(gt2,0))*ZA(gt3,0) + Conj(UP1(gt2,1))*(-1.4142135623730951*LamSD*ZA(gt3,2) + LamTD*ZA(gt3,3))));

   const std::complex<double> right = std::complex<double>(0,-0.5)*(UM1(gt2,0)*(-1.4142135623730951*Conj(LamTD)*UP1(gt1,1)*ZA(gt3,0) + 2*g2*UP1(gt1,0)*ZA(gt3,3)) + UM1(gt2,1)*(1.4142135623730951*g2*UP1(gt1,0)*ZA(gt3,0) + UP1(gt1,1)*(-1.4142135623730951*Conj(LamSD)*ZA(gt3,2) + Conj(LamTD)*ZA(gt3,3))));

   return {left, right};
}

ChiralVertex VertexImpl<fields::VP, fields::Cha2, typename fields::bar<fields::Cha2>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UP2 = MODELPARAMETER(UP2);
   const auto UM2 = MODELPARAMETER(UM2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -(g2*Conj(UP2(gt1,0))*Sin(ThetaW)*UP2(gt2,0)) - 0.5*Conj(UP2(gt1,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UP2(gt2,1);

   const std::complex<double> right = -(g2*Conj(UM2(gt2,0))*Sin(ThetaW)*UM2(gt1,0)) - 0.5*Conj(UM2(gt2,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UM2(gt1,1);

   return {left, right};
}

ChiralVertex VertexImpl<fields::Ah, typename fields::bar<fields::Cha2>::type, fields::Cha2>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto UM2 = MODELPARAMETER(UM2);
   const auto UP2 = MODELPARAMETER(UP2);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0,0.5)*(g2*Conj(UM2(gt2,0))*(1.4142135623730951*Conj(UP2(gt1,1))*ZA(gt3,1) - 2*Conj(UP2(gt1,0))*ZA(gt3,3)) + Conj(UM2(gt2,1))*(1.4142135623730951*LamTU*Conj(UP2(gt1,0))*ZA(gt3,1) + Conj(UP2(gt1,1))*(1.4142135623730951*LamSU*ZA(gt3,2) + LamTU*ZA(gt3,3))));

   const std::complex<double> right = std::complex<double>(0,-0.5)*(1.4142135623730951*Conj(LamSU)*UM2(gt1,1)*UP2(gt2,1)*ZA(gt3,2) + g2*UM2(gt1,0)*(1.4142135623730951*UP2(gt2,1)*ZA(gt3,1) - 2*UP2(gt2,0)*ZA(gt3,3)) + Conj(LamTU)*UM2(gt1,1)*(1.4142135623730951*UP2(gt2,0)*ZA(gt3,1) + UP2(gt2,1)*ZA(gt3,3)));

   return {left, right};
}

ChiralVertex VertexImpl<fields::VP, fields::Fd, typename fields::bar<fields::Fd>::type>::evaluate(
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

ChiralVertex VertexImpl<fields::Ah, typename fields::bar<fields::Fd>::type, fields::Fd>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto ZDR = MODELPARAMETER(ZDR);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,Conj(ZDL(gt2,j2))*SUM(j1,0,2,Conj(ZDR(gt1,j1))*Yd(j1,j2)))*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*ZDR(gt2,j1))*ZDL(gt1,j2))*ZA(gt3,0);

   return {left, right};
}

ChiralVertex VertexImpl<fields::VP, fields::Fe, typename fields::bar<fields::Fe>::type>::evaluate(
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

ChiralVertex VertexImpl<fields::Ah, typename fields::bar<fields::Fe>::type, fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,j2)))*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZEL(gt1,j2))*ZA(gt3,0);

   return {left, right};
}

ChiralVertex VertexImpl<fields::VP, fields::Fu, typename fields::bar<fields::Fu>::type>::evaluate(
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

ChiralVertex VertexImpl<fields::Ah, typename fields::bar<fields::Fu>::type, fields::Fu>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto ZUR = MODELPARAMETER(ZUR);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,Conj(ZUL(gt2,j2))*SUM(j1,0,2,Conj(ZUR(gt1,j1))*Yu(j1,j2)))*ZA(gt3,1);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*ZUR(gt2,j1))*ZUL(gt1,j2))*ZA(gt3,1);

   return {left, right};
}

MomentumDifferenceVertex VertexImpl<fields::VP, fields::Hpm, typename fields::conj<fields::Hpm>::type>::evaluate(
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

   const std::complex<double> result = 0.5*((0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*ZP(gt1,0)*ZP(gt2,0) + (0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*ZP(gt1,1)*ZP(gt2,1) + 2*g2*Sin(ThetaW)*(ZP(gt1,2)*ZP(gt2,2) + ZP(gt1,3)*ZP(gt2,3)));

   return {result, minuend_index, subtrahend_index};
}

ScalarVertex VertexImpl<fields::Ah, typename fields::conj<fields::Hpm>::type, fields::Hpm>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto vu = MODELPARAMETER(vu);
   const auto vd = MODELPARAMETER(vd);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto vT = MODELPARAMETER(vT);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto vS = MODELPARAMETER(vS);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto MuD = MODELPARAMETER(MuD);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto MuU = MODELPARAMETER(MuU);
   const auto g1 = MODELPARAMETER(g1);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,0.05)*(ZA(gt1,2)*(10*LamTD*vd*Conj(LamSD)*ZP(gt2,2)*ZP(gt3,0) - 10*LamSD*vd*Conj(LamTD)*ZP(gt2,3)*ZP(gt3,0) - 7.745966692414834*g1*MDBS*ZP(gt2,1)*ZP(gt3,1) + 14.142135623730951*MuU*Conj(LamSU)*ZP(gt2,1)*ZP(gt3,1) + 7.0710678118654755*LamTU*vT*Conj(LamSU)*ZP(gt2,1)*ZP(gt3,1) - 7.0710678118654755*LamSU*vT*Conj(LamTU)*ZP(gt2,1)*ZP(gt3,1) + 7.745966692414834*g1*Conj(MDBS)*ZP(gt2,1)*ZP(gt3,1) - 14.142135623730951*LamSU*Conj(MuU)*ZP(gt2,1)*ZP(gt3,1) + 10*LamTU*vu*Conj(LamSU)*ZP(gt2,2)*ZP(gt3,1) - 10*LamSU*vu*Conj(LamTU)*ZP(gt2,3)*ZP(gt3,1) - 10*LamSU*vu*Conj(LamTU)*ZP(gt2,1)*ZP(gt3,2) + 10*LamTU*vu*Conj(LamSU)*ZP(gt2,1)*ZP(gt3,3) + ZP(gt2,0)*((7.745966692414834*g1*MDBS + 7.0710678118654755*(2*MuD - LamTD*vT)*Conj(LamSD) + 7.0710678118654755*LamSD*vT*Conj(LamTD) - 7.745966692414834*g1*Conj(MDBS) - 14.142135623730951*LamSD*Conj(MuD))*ZP(gt3,0) + 10*vd*(-(LamSD*Conj(LamTD)*ZP(gt3,2)) + LamTD*Conj(LamSD)*ZP(gt3,3)))) - 5*(ZA(gt1,0)*(vu*Sqr(g2)*ZP(gt2,1)*ZP(gt3,0) + (2*LamTD*vS*Conj(LamSD) + 1.4142135623730951*(2*g2*MDWBT - vT*AbsSqr(LamTD) + 2*LamTD*Conj(MuD) + vT*Sqr(g2)))*ZP(gt2,2)*ZP(gt3,0) + 1.4142135623730951*vT*AbsSqr(LamTD)*ZP(gt2,3)*ZP(gt3,0) + 2.8284271247461903*MuD*Conj(LamTD)*ZP(gt2,3)*ZP(gt3,0) + 2*LamSD*vS*Conj(LamTD)*ZP(gt2,3)*ZP(gt3,0) + 2.8284271247461903*g2*Conj(MDWBT)*ZP(gt2,3)*ZP(gt3,0) - 1.4142135623730951*vT*Sqr(g2)*ZP(gt2,3)*ZP(gt3,0) - vu*Sqr(g2)*ZP(gt2,0)*ZP(gt3,1) + 1.4142135623730951*vT*AbsSqr(LamTD)*ZP(gt2,0)*ZP(gt3,2) - 2.8284271247461903*MuD*Conj(LamTD)*ZP(gt2,0)*ZP(gt3,2) - 2*LamSD*vS*Conj(LamTD)*ZP(gt2,0)*ZP(gt3,2) - 2.8284271247461903*g2*Conj(MDWBT)*ZP(gt2,0)*ZP(gt3,2) - 1.4142135623730951*vT*Sqr(g2)*ZP(gt2,0)*ZP(gt3,2) - 2.8284271247461903*g2*MDWBT*ZP(gt2,0)*ZP(gt3,3) - 1.4142135623730951*vT*AbsSqr(LamTD)*ZP(gt2,0)*ZP(gt3,3) - 2*LamTD*vS*Conj(LamSD)*ZP(gt2,0)*ZP(gt3,3) - 2.8284271247461903*LamTD*Conj(MuD)*ZP(gt2,0)*ZP(gt3,3) + 1.4142135623730951*vT*Sqr(g2)*ZP(gt2,0)*ZP(gt3,3)) + ZA(gt1,1)*(-((vd*Sqr(g2)*ZP(gt2,0) + (2*LamTU*vS*Conj(LamSU) + 1.4142135623730951*(2*g2*MDWBT - vT*AbsSqr(LamTU) + 2*LamTU*Conj(MuU) + vT*Sqr(g2)))*ZP(gt2,2) + ((2.8284271247461903*MuU + 2*LamSU*vS + 1.4142135623730951*LamTU*vT)*Conj(LamTU) + 1.4142135623730951*g2*(-(g2*vT) + 2*Conj(MDWBT)))*ZP(gt2,3))*ZP(gt3,1)) + ZP(gt2,1)*(vd*Sqr(g2)*ZP(gt3,0) + ((2.8284271247461903*MuU + 2*LamSU*vS - 1.4142135623730951*LamTU*vT)*Conj(LamTU) + 1.4142135623730951*g2*(g2*vT + 2*Conj(MDWBT)))*ZP(gt3,2) + (2*LamTU*vS*Conj(LamSU) + 1.4142135623730951*(2*g2*MDWBT + vT*AbsSqr(LamTU) + 2*LamTU*Conj(MuU) - vT*Sqr(g2)))*ZP(gt3,3))) + ZA(gt1,3)*(1.4142135623730951*vd*AbsSqr(LamTD)*ZP(gt2,3)*ZP(gt3,0) - 1.4142135623730951*vd*Sqr(g2)*ZP(gt2,3)*ZP(gt3,0) + 2*g2*MDWBT*ZP(gt2,1)*ZP(gt3,1) + 1.4142135623730951*LamTU*vS*Conj(LamSU)*ZP(gt2,1)*ZP(gt3,1) - 2*MuU*Conj(LamTU)*ZP(gt2,1)*ZP(gt3,1) - 1.4142135623730951*LamSU*vS*Conj(LamTU)*ZP(gt2,1)*ZP(gt3,1) - 2*g2*Conj(MDWBT)*ZP(gt2,1)*ZP(gt3,1) + 2*LamTU*Conj(MuU)*ZP(gt2,1)*ZP(gt3,1) + 1.4142135623730951*vu*AbsSqr(LamTU)*ZP(gt2,3)*ZP(gt3,1) - 1.4142135623730951*vu*Sqr(g2)*ZP(gt2,3)*ZP(gt3,1) - 1.4142135623730951*vu*AbsSqr(LamTU)*ZP(gt2,1)*ZP(gt3,2) + 1.4142135623730951*vu*Sqr(g2)*ZP(gt2,1)*ZP(gt3,2) - 4*vT*Sqr(g2)*ZP(gt2,3)*ZP(gt3,2) - 1.4142135623730951*vu*AbsSqr(LamTU)*ZP(gt2,1)*ZP(gt3,3) + 1.4142135623730951*vu*Sqr(g2)*ZP(gt2,1)*ZP(gt3,3) + ZP(gt2,2)*(1.4142135623730951*vd*(AbsSqr(LamTD) - Sqr(g2))*ZP(gt3,0) + 1.4142135623730951*vu*(AbsSqr(LamTU) - Sqr(g2))*ZP(gt3,1) + 4*vT*Sqr(g2)*ZP(gt3,3)) + ZP(gt2,0)*((-2*g2*MDWBT - 1.4142135623730951*LamTD*vS*Conj(LamSD) + 2*MuD*Conj(LamTD) + 1.4142135623730951*LamSD*vS*Conj(LamTD) + 2*g2*Conj(MDWBT) - 2*LamTD*Conj(MuD))*ZP(gt3,0) + 1.4142135623730951*vd*(-AbsSqr(LamTD) + Sqr(g2))*(ZP(gt3,2) + ZP(gt3,3))))));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::VP, fields::Sd, typename fields::conj<fields::Sd>::type>::evaluate(
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

ScalarVertex VertexImpl<fields::Ah, typename fields::conj<fields::Sd>::type, fields::Sd>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Mu = MODELPARAMETER(Mu);
   const auto g1 = MODELPARAMETER(g1);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,-0.16666666666666666)*(4.242640687119286*Conj(Mu)*SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1)))*ZA(gt1,1) - 4.242640687119286*Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1)))*ZD(gt3,j2))*ZA(gt1,1) + 0.7745966692414834*g1*MDBS*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*ZA(gt1,2) - 0.7745966692414834*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*ZA(gt1,2) + 1.5491933384829668*g1*MDBS*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*ZA(gt1,2) - 1.5491933384829668*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*ZA(gt1,2) - 3*g2*MDWBT*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*ZA(gt1,3) + 3*g2*Conj(MDWBT)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*ZA(gt1,3));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::VP, fields::Se, typename fields::conj<fields::Se>::type>::evaluate(
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

ScalarVertex VertexImpl<fields::Ah, typename fields::conj<fields::Se>::type, fields::Se>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Mu = MODELPARAMETER(Mu);
   const auto g1 = MODELPARAMETER(g1);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,-0.1)*(7.0710678118654755*Conj(Mu)*SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)))*ZA(gt1,1) - 7.0710678118654755*Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1)))*ZE(gt3,j2))*ZA(gt1,1) - 3.872983346207417*g1*MDBS*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*ZA(gt1,2) + 3.872983346207417*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*ZA(gt1,2) + 7.745966692414834*g1*MDBS*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*ZA(gt1,2) - 7.745966692414834*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*ZA(gt1,2) - 5*g2*MDWBT*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*ZA(gt1,3) + 5*g2*Conj(MDWBT)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*ZA(gt1,3));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::VP, fields::SRdp, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.1*(-3.872983346207417*g1*Cos(ThetaW) - 5*g2*Sin(ThetaW));

   return {result, minuend_index, subtrahend_index};
}

ScalarVertex VertexImpl<fields::Ah, typename fields::conj<fields::SRdp>::type, fields::SRdp>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MuD = MODELPARAMETER(MuD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto vT = MODELPARAMETER(vT);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto vS = MODELPARAMETER(vS);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,-0.25)*((1.5491933384829668*g1*MDBS + 1.4142135623730951*(-2*MuD + LamTD*vT)*Conj(LamSD) - 1.4142135623730951*LamSD*vT*Conj(LamTD) - 1.5491933384829668*g1*Conj(MDBS) + 2.8284271247461903*LamSD*Conj(MuD))*ZA(gt1,2) + (2*g2*MDWBT - 1.4142135623730951*LamTD*vS*Conj(LamSD) + (2*MuD + 1.4142135623730951*LamSD*vS)*Conj(LamTD) - 2*g2*Conj(MDWBT) - 2*LamTD*Conj(MuD))*ZA(gt1,3));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::VP, fields::SRum, typename fields::conj<fields::SRum>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   return {result, minuend_index, subtrahend_index};
}

ScalarVertex VertexImpl<fields::Ah, typename fields::conj<fields::SRum>::type, fields::SRum>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MuU = MODELPARAMETER(MuU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto vT = MODELPARAMETER(vT);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto vS = MODELPARAMETER(vS);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,-0.25)*((-1.5491933384829668*g1*MDBS - 1.4142135623730951*(2*MuU + LamTU*vT)*Conj(LamSU) + 1.4142135623730951*LamSU*vT*Conj(LamTU) + 1.5491933384829668*g1*Conj(MDBS) + 2.8284271247461903*LamSU*Conj(MuU))*ZA(gt1,2) + (-2*g2*MDWBT + 1.4142135623730951*LamTU*vS*Conj(LamSU) - (2*MuU + 1.4142135623730951*LamSU*vS)*Conj(LamTU) + 2*g2*Conj(MDWBT) + 2*LamTU*Conj(MuU))*ZA(gt1,3));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::VP, fields::Su, typename fields::conj<fields::Su>::type>::evaluate(
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

ScalarVertex VertexImpl<fields::Ah, typename fields::conj<fields::Su>::type, fields::Su>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Mu = MODELPARAMETER(Mu);
   const auto g1 = MODELPARAMETER(g1);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,-0.16666666666666666)*(4.242640687119286*Conj(Mu)*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt3,3 + j1)))*ZA(gt1,0) - 4.242640687119286*Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1)))*ZU(gt3,j2))*ZA(gt1,0) + 0.7745966692414834*g1*MDBS*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*ZA(gt1,2) - 0.7745966692414834*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*ZA(gt1,2) - 3.0983866769659336*g1*MDBS*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*ZA(gt1,2) + 3.0983866769659336*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*ZA(gt1,2) + 3*g2*MDWBT*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*ZA(gt1,3) - 3*g2*Conj(MDWBT)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*ZA(gt1,3));

   return {result};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fe, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,j2)))*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZEL(gt1,j2))*ZH(gt3,0);

   return {left, right};
}

ChiralVertex VertexImpl<fields::hh, typename fields::bar<fields::Cha1>::type, fields::Cha1>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = 0.5*(-(Conj(UM1(gt1,0))*(1.4142135623730951*LamTD*Conj(UP1(gt2,1))*ZH(gt3,0) + 2*g2*Conj(UP1(gt2,0))*ZH(gt3,3))) - Conj(UM1(gt1,1))*(1.4142135623730951*g2*Conj(UP1(gt2,0))*ZH(gt3,0) + Conj(UP1(gt2,1))*(1.4142135623730951*LamSD*ZH(gt3,2) - LamTD*ZH(gt3,3))));

   const std::complex<double> right = 0.5*(-(UM1(gt2,0)*(1.4142135623730951*Conj(LamTD)*UP1(gt1,1)*ZH(gt3,0) + 2*g2*UP1(gt1,0)*ZH(gt3,3))) - UM1(gt2,1)*(1.4142135623730951*g2*UP1(gt1,0)*ZH(gt3,0) + UP1(gt1,1)*(1.4142135623730951*Conj(LamSD)*ZH(gt3,2) - Conj(LamTD)*ZH(gt3,3))));

   return {left, right};
}

ChiralVertex VertexImpl<fields::hh, typename fields::bar<fields::Cha2>::type, fields::Cha2>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto UM2 = MODELPARAMETER(UM2);
   const auto UP2 = MODELPARAMETER(UP2);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = 0.5*(g2*Conj(UM2(gt2,0))*(-1.4142135623730951*Conj(UP2(gt1,1))*ZH(gt3,1) + 2*Conj(UP2(gt1,0))*ZH(gt3,3)) + Conj(UM2(gt2,1))*(1.4142135623730951*LamTU*Conj(UP2(gt1,0))*ZH(gt3,1) + Conj(UP2(gt1,1))*(1.4142135623730951*LamSU*ZH(gt3,2) + LamTU*ZH(gt3,3))));

   const std::complex<double> right = 0.5*(1.4142135623730951*Conj(LamSU)*UM2(gt1,1)*UP2(gt2,1)*ZH(gt3,2) + UM2(gt1,0)*(-1.4142135623730951*g2*UP2(gt2,1)*ZH(gt3,1) + 2*g2*UP2(gt2,0)*ZH(gt3,3)) + Conj(LamTU)*UM2(gt1,1)*(1.4142135623730951*UP2(gt2,0)*ZH(gt3,1) + UP2(gt2,1)*ZH(gt3,3)));

   return {left, right};
}

ChiralVertex VertexImpl<fields::hh, typename fields::bar<fields::Fd>::type, fields::Fd>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto ZDR = MODELPARAMETER(ZDR);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(ZDL(gt2,j2))*SUM(j1,0,2,Conj(ZDR(gt1,j1))*Yd(j1,j2)))*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*ZDR(gt2,j1))*ZDL(gt1,j2))*ZH(gt3,0);

   return {left, right};
}

ChiralVertex VertexImpl<fields::hh, typename fields::bar<fields::Fe>::type, fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,j2)))*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZEL(gt1,j2))*ZH(gt3,0);

   return {left, right};
}

ChiralVertex VertexImpl<fields::hh, typename fields::bar<fields::Fu>::type, fields::Fu>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto ZUR = MODELPARAMETER(ZUR);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(ZUL(gt2,j2))*SUM(j1,0,2,Conj(ZUR(gt1,j1))*Yu(j1,j2)))*ZH(gt3,1);

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*ZUR(gt2,j1))*ZUL(gt1,j2))*ZH(gt3,1);

   return {left, right};
}

ScalarVertex VertexImpl<fields::hh, typename fields::conj<fields::Hpm>::type, fields::Hpm>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MuD = MODELPARAMETER(MuD);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto vS = MODELPARAMETER(vS);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto vT = MODELPARAMETER(vT);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto vd = MODELPARAMETER(vd);
   const auto MuU = MODELPARAMETER(MuU);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto vu = MODELPARAMETER(vu);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.05*(7.745966692414834*g1*MDBS*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) - 20*vS*AbsSqr(LamSD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) - 14.142135623730951*MuD*Conj(LamSD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) + 7.0710678118654755*LamTD*vT*Conj(LamSD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) + 7.0710678118654755*LamSD*vT*Conj(LamTD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) + 7.745966692414834*g1*Conj(MDBS)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) - 14.142135623730951*LamSD*Conj(MuD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) + 10*g2*MDWBT*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) - 10*vT*AbsSqr(LamTD)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) + 7.0710678118654755*LamTD*vS*Conj(LamSD)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) + 10*MuD*Conj(LamTD)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) + 7.0710678118654755*LamSD*vS*Conj(LamTD)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) + 10*g2*Conj(MDWBT)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) + 10*LamTD*Conj(MuD)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) - 10*LamTD*vd*Conj(LamSD)*ZH(gt1,2)*ZP(gt2,2)*ZP(gt3,0) + 7.0710678118654755*vd*AbsSqr(LamTD)*ZH(gt1,3)*ZP(gt2,2)*ZP(gt3,0) - 7.0710678118654755*vd*Sqr(g2)*ZH(gt1,3)*ZP(gt2,2)*ZP(gt3,0) - 10*LamSD*vd*Conj(LamTD)*ZH(gt1,2)*ZP(gt2,3)*ZP(gt3,0) - 7.0710678118654755*vd*AbsSqr(LamTD)*ZH(gt1,3)*ZP(gt2,3)*ZP(gt3,0) + 7.0710678118654755*vd*Sqr(g2)*ZH(gt1,3)*ZP(gt2,3)*ZP(gt3,0) - 7.745966692414834*g1*MDBS*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) - 20*vS*AbsSqr(LamSU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) - 14.142135623730951*MuU*Conj(LamSU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) - 7.0710678118654755*LamTU*vT*Conj(LamSU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) - 7.0710678118654755*LamSU*vT*Conj(LamTU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) - 7.745966692414834*g1*Conj(MDBS)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) - 14.142135623730951*LamSU*Conj(MuU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) - 10*g2*MDWBT*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,1) - 10*vT*AbsSqr(LamTU)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,1) - 7.0710678118654755*LamTU*vS*Conj(LamSU)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,1) - 10*MuU*Conj(LamTU)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,1) - 7.0710678118654755*LamSU*vS*Conj(LamTU)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,1) - 10*g2*Conj(MDWBT)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,1) - 10*LamTU*Conj(MuU)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,1) - 10*LamTU*vu*Conj(LamSU)*ZH(gt1,2)*ZP(gt2,2)*ZP(gt3,1) + 7.0710678118654755*vu*AbsSqr(LamTU)*ZH(gt1,3)*ZP(gt2,2)*ZP(gt3,1) - 7.0710678118654755*vu*Sqr(g2)*ZH(gt1,3)*ZP(gt2,2)*ZP(gt3,1) - 10*LamSU*vu*Conj(LamTU)*ZH(gt1,2)*ZP(gt2,3)*ZP(gt3,1) - 7.0710678118654755*vu*AbsSqr(LamTU)*ZH(gt1,3)*ZP(gt2,3)*ZP(gt3,1) + 7.0710678118654755*vu*Sqr(g2)*ZH(gt1,3)*ZP(gt2,3)*ZP(gt3,1) - 10*LamSD*vd*Conj(LamTD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,2) + 7.0710678118654755*vd*AbsSqr(LamTD)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,2) - 7.0710678118654755*vd*Sqr(g2)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,2) - 10*LamSU*vu*Conj(LamTU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,2) + 7.0710678118654755*vu*AbsSqr(LamTU)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,2) - 7.0710678118654755*vu*Sqr(g2)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,2) - 20*vT*Sqr(g2)*ZH(gt1,3)*ZP(gt2,2)*ZP(gt3,2) + 20*vT*Sqr(g2)*ZH(gt1,3)*ZP(gt2,3)*ZP(gt3,2) - 10*LamTD*vd*Conj(LamSD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,3) - 7.0710678118654755*vd*AbsSqr(LamTD)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,3) + 7.0710678118654755*vd*Sqr(g2)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,3) - 10*LamTU*vu*Conj(LamSU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,3) - 7.0710678118654755*vu*AbsSqr(LamTU)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,3) + 7.0710678118654755*vu*Sqr(g2)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,3) + 20*vT*Sqr(g2)*ZH(gt1,3)*ZP(gt2,2)*ZP(gt3,3) - 20*vT*Sqr(g2)*ZH(gt1,3)*ZP(gt2,3)*ZP(gt3,3) - ZH(gt1,1)*(ZP(gt2,0)*((-3*vu*Sqr(g1) + 5*vu*Sqr(g2))*ZP(gt3,0) + 5*vd*Sqr(g2)*ZP(gt3,1)) + 5*ZP(gt2,2)*((2*LamTU*vS*Conj(LamSU) + 1.4142135623730951*(2*g2*MDWBT - vT*AbsSqr(LamTU) + 2*LamTU*Conj(MuU) + vT*Sqr(g2)))*ZP(gt3,1) + 2*vu*Sqr(g2)*ZP(gt3,2)) + 5*ZP(gt2,3)*(((2.8284271247461903*MuU + 2*LamSU*vS + 1.4142135623730951*LamTU*vT)*Conj(LamTU) + 1.4142135623730951*g2*(-(g2*vT) + 2*Conj(MDWBT)))*ZP(gt3,1) - 2*vu*(-2*AbsSqr(LamTU) + Sqr(g2))*ZP(gt3,3)) + ZP(gt2,1)*(5*vd*Sqr(g2)*ZP(gt3,0) + vu*(3*Sqr(g1) + 5*Sqr(g2))*ZP(gt3,1) + 5*((2.8284271247461903*MuU + 2*LamSU*vS - 1.4142135623730951*LamTU*vT)*Conj(LamTU) + 1.4142135623730951*g2*(g2*vT + 2*Conj(MDWBT)))*ZP(gt3,2) + 5*(2*LamTU*vS*Conj(LamSU) + 1.4142135623730951*(2*g2*MDWBT + vT*AbsSqr(LamTU) + 2*LamTU*Conj(MuU) - vT*Sqr(g2)))*ZP(gt3,3))) - ZH(gt1,0)*(ZP(gt2,1)*(5*vu*Sqr(g2)*ZP(gt3,0) + vd*(-3*Sqr(g1) + 5*Sqr(g2))*ZP(gt3,1)) + 5*(ZP(gt2,2)*((2*LamTD*vS*Conj(LamSD) + 1.4142135623730951*(2*g2*MDWBT - vT*AbsSqr(LamTD) + 2*LamTD*Conj(MuD) + vT*Sqr(g2)))*ZP(gt3,0) - 2*vd*(-2*AbsSqr(LamTD) + Sqr(g2))*ZP(gt3,2)) + ZP(gt2,3)*(((2.8284271247461903*MuD + 2*LamSD*vS + 1.4142135623730951*LamTD*vT)*Conj(LamTD) + 1.4142135623730951*g2*(-(g2*vT) + 2*Conj(MDWBT)))*ZP(gt3,0) + 2*vd*Sqr(g2)*ZP(gt3,3))) + ZP(gt2,0)*(vd*(3*Sqr(g1) + 5*Sqr(g2))*ZP(gt3,0) + 5*(vu*Sqr(g2)*ZP(gt3,1) + ((2.8284271247461903*MuD + 2*LamSD*vS - 1.4142135623730951*LamTD*vT)*Conj(LamTD) + 1.4142135623730951*g2*(g2*vT + 2*Conj(MDWBT)))*ZP(gt3,2) + (2*LamTD*vS*Conj(LamSD) + 1.4142135623730951*(2*g2*MDWBT + vT*AbsSqr(LamTD) + 2*LamTD*Conj(MuD) - vT*Sqr(g2)))*ZP(gt3,3)))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, typename fields::conj<fields::Sd>::type, fields::Sd>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto vd = MODELPARAMETER(vd);
   const auto Mu = MODELPARAMETER(Mu);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.016666666666666666*(30*(-2*vd*SUM(j3,0,2,Conj(ZD(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,j1))*ZD(gt3,3 + j2)))*ZH(gt1,0) - 2*vd*SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gt3,j3))*ZH(gt1,0) + 1.4142135623730951*(Conj(Mu)*SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1))) + Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1)))*ZD(gt3,j2)))*ZH(gt1,1)) + 2*g1*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*(3*g1*vd*ZH(gt1,0) - 3*g1*vu*ZH(gt1,1) - 7.745966692414834*(MDBS + Conj(MDBS))*ZH(gt1,2)) + SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*(3*vd*(Sqr(g1) + 5*Sqr(g2))*ZH(gt1,0) - 3*vu*(Sqr(g1) + 5*Sqr(g2))*ZH(gt1,1) - 7.745966692414834*g1*MDBS*ZH(gt1,2) - 7.745966692414834*g1*Conj(MDBS)*ZH(gt1,2) + 30*g2*MDWBT*ZH(gt1,3) + 30*g2*Conj(MDWBT)*ZH(gt1,3)));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, typename fields::conj<fields::Se>::type, fields::Se>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto vd = MODELPARAMETER(vd);
   const auto Mu = MODELPARAMETER(Mu);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.05*(2*(5*(-2*vd*SUM(j3,0,2,Conj(ZE(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Ye(j2,j1))*ZE(gt3,3 + j2)))*ZH(gt1,0) - 2*vd*SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZE(gt3,j3))*ZH(gt1,0) + 1.4142135623730951*(Conj(Mu)*SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1))) + Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1)))*ZE(gt3,j2)))*ZH(gt1,1)) + g1*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*(3*g1*vd*ZH(gt1,0) - 3*g1*vu*ZH(gt1,1) - 7.745966692414834*(MDBS + Conj(MDBS))*ZH(gt1,2))) + SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*((-3*vd*Sqr(g1) + 5*vd*Sqr(g2))*ZH(gt1,0) + vu*(3*Sqr(g1) - 5*Sqr(g2))*ZH(gt1,1) + 7.745966692414834*g1*(MDBS + Conj(MDBS))*ZH(gt1,2) + 10*g2*(MDWBT + Conj(MDWBT))*ZH(gt1,3)));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, typename fields::conj<fields::SRdp>::type, fields::SRdp>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto vu = MODELPARAMETER(vu);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MuD = MODELPARAMETER(MuD);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto vS = MODELPARAMETER(vS);
   const auto vT = MODELPARAMETER(vT);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.25*(vd*(-4*AbsSqr(LamTD) + 0.6*Sqr(g1) - Sqr(g2))*ZH(gt1,0) + vu*(-0.6*Sqr(g1) + Sqr(g2))*ZH(gt1,1) - 1.5491933384829668*g1*MDBS*ZH(gt1,2) - 4*vS*AbsSqr(LamSD)*ZH(gt1,2) - 2.8284271247461903*MuD*Conj(LamSD)*ZH(gt1,2) + 1.4142135623730951*LamTD*vT*Conj(LamSD)*ZH(gt1,2) + 1.4142135623730951*LamSD*vT*Conj(LamTD)*ZH(gt1,2) - 1.5491933384829668*g1*Conj(MDBS)*ZH(gt1,2) - 2.8284271247461903*LamSD*Conj(MuD)*ZH(gt1,2) - 2*g2*MDWBT*ZH(gt1,3) - 2*vT*AbsSqr(LamTD)*ZH(gt1,3) + 1.4142135623730951*LamTD*vS*Conj(LamSD)*ZH(gt1,3) + 2*MuD*Conj(LamTD)*ZH(gt1,3) + 1.4142135623730951*LamSD*vS*Conj(LamTD)*ZH(gt1,3) - 2*g2*Conj(MDWBT)*ZH(gt1,3) + 2*LamTD*Conj(MuD)*ZH(gt1,3));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, typename fields::conj<fields::SRum>::type, fields::SRum>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MuU = MODELPARAMETER(MuU);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto vS = MODELPARAMETER(vS);
   const auto vT = MODELPARAMETER(vT);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.25*(-0.2*vd*(3*Sqr(g1) - 5*Sqr(g2))*ZH(gt1,0) - vu*(4*AbsSqr(LamTU) - 0.6*Sqr(g1) + Sqr(g2))*ZH(gt1,1) + 1.5491933384829668*g1*MDBS*ZH(gt1,2) - 4*vS*AbsSqr(LamSU)*ZH(gt1,2) - 2.8284271247461903*MuU*Conj(LamSU)*ZH(gt1,2) - 1.4142135623730951*LamTU*vT*Conj(LamSU)*ZH(gt1,2) - 1.4142135623730951*LamSU*vT*Conj(LamTU)*ZH(gt1,2) + 1.5491933384829668*g1*Conj(MDBS)*ZH(gt1,2) - 2.8284271247461903*LamSU*Conj(MuU)*ZH(gt1,2) + 2*g2*MDWBT*ZH(gt1,3) - 2*vT*AbsSqr(LamTU)*ZH(gt1,3) - 1.4142135623730951*LamTU*vS*Conj(LamSU)*ZH(gt1,3) - 2*MuU*Conj(LamTU)*ZH(gt1,3) - 1.4142135623730951*LamSU*vS*Conj(LamTU)*ZH(gt1,3) + 2*g2*Conj(MDWBT)*ZH(gt1,3) - 2*LamTU*Conj(MuU)*ZH(gt1,3));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, typename fields::conj<fields::Su>::type, fields::Su>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Mu = MODELPARAMETER(Mu);
   const auto vu = MODELPARAMETER(vu);
   const auto g1 = MODELPARAMETER(g1);
   const auto vd = MODELPARAMETER(vd);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.016666666666666666*(30*(1.4142135623730951*Conj(Mu)*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt3,3 + j1)))*ZH(gt1,0) + 1.4142135623730951*Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1)))*ZU(gt3,j2))*ZH(gt1,0) - 2*vu*(SUM(j3,0,2,Conj(ZU(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))*ZU(gt3,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gt3,j3)))*ZH(gt1,1)) + 4*g1*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*(-3*g1*vd*ZH(gt1,0) + 3*g1*vu*ZH(gt1,1) + 7.745966692414834*(MDBS + Conj(MDBS))*ZH(gt1,2)) + SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*(3*vd*(Sqr(g1) - 5*Sqr(g2))*ZH(gt1,0) - 3*vu*(Sqr(g1) - 5*Sqr(g2))*ZH(gt1,1) - 7.745966692414834*g1*(MDBS + Conj(MDBS))*ZH(gt1,2) - 30*g2*(MDWBT + Conj(MDWBT))*ZH(gt1,3)));

   return {result};
}

TripleVectorVertex VertexImpl<fields::VP, fields::VWm, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Sin(ThetaW);

   return {result, TripleVectorVertex::odd_permutation{}};
}

InverseMetricVertex VertexImpl<fields::hh, typename fields::conj<fields::VWm>::type, fields::VWm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto vT = MODELPARAMETER(vT);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.5*Sqr(g2)*(vd*ZH(gt1,0) + vu*ZH(gt1,1) + 4*vT*ZH(gt1,3));

   return {result};
}

ChiralVertex VertexImpl<fields::Fe, typename fields::bar<fields::Fe>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.7745966692414834*g1*KroneckerDelta(gt1,gt2)*Sin(ThetaW);

   const std::complex<double> right = 0.1*KroneckerDelta(gt1,gt2)*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW));

   return {left, right};
}

ChiralVertex VertexImpl<fields::VZ, fields::Cha1, typename fields::bar<fields::Cha1>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = g2*Conj(UM1(gt1,0))*Cos(ThetaW)*UM1(gt2,0) + 0.5*Conj(UM1(gt1,1))*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*UM1(gt2,1);

   const std::complex<double> right = g2*Conj(UP1(gt2,0))*Cos(ThetaW)*UP1(gt1,0) + 0.5*Conj(UP1(gt2,1))*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*UP1(gt1,1);

   return {left, right};
}

ChiralVertex VertexImpl<fields::VZ, fields::Cha2, typename fields::bar<fields::Cha2>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UP2 = MODELPARAMETER(UP2);
   const auto UM2 = MODELPARAMETER(UM2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -(g2*Conj(UP2(gt1,0))*Cos(ThetaW)*UP2(gt2,0)) + 0.1*Conj(UP2(gt1,1))*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*UP2(gt2,1);

   const std::complex<double> right = -(g2*Conj(UM2(gt2,0))*Cos(ThetaW)*UM2(gt1,0)) + 0.1*Conj(UM2(gt2,1))*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*UM2(gt1,1);

   return {left, right};
}

ChiralVertex VertexImpl<fields::VZ, fields::Fd, typename fields::bar<fields::Fd>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.2581988897471611*g1*KroneckerDelta(gt1,gt2)*Sin(ThetaW);

   const std::complex<double> right = -0.16666666666666666*KroneckerDelta(gt1,gt2)*(3*g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW));

   return {left, right};
}

ChiralVertex VertexImpl<fields::VZ, fields::Fe, typename fields::bar<fields::Fe>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.7745966692414834*g1*KroneckerDelta(gt1,gt2)*Sin(ThetaW);

   const std::complex<double> right = 0.1*KroneckerDelta(gt1,gt2)*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW));

   return {left, right};
}

ChiralVertex VertexImpl<fields::VZ, fields::Fu, typename fields::bar<fields::Fu>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.5163977794943222*g1*KroneckerDelta(gt1,gt2)*Sin(ThetaW);

   const std::complex<double> right = 0.16666666666666666*KroneckerDelta(gt1,gt2)*(3*g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Fe, fields::Ah, typename fields::bar<fields::Fe>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt3 = indices[1];
   const int gt1 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,j2)))*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZEL(gt1,j2))*ZA(gt3,0);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Ah, fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,j2)))*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZEL(gt1,j2))*ZA(gt3,0);

   return {left, right};
}

ChiralVertex VertexImpl<fields::Fe, fields::Chi, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZN1 = MODELPARAMETER(ZN1);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZN2 = MODELPARAMETER(ZN2);

   const std::complex<double> left = 0.7071067811865475*(0.7745966692414834*g1*Conj(ZN1(gt1,0)) + g2*Conj(ZN1(gt1,1)))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZE(gt3,j1));

   const std::complex<double> right = -(SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZE(gt3,j2))*ZN2(gt1,2));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, typename fields::bar<fields::Chi>::type, fields::Se>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZN2 = MODELPARAMETER(ZN2);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZN1 = MODELPARAMETER(ZN1);

   const std::complex<double> left = -(Conj(ZN2(gt1,2))*SUM(j2,0,2,Conj(ZE(gt3,j2))*SUM(j1,0,2,Conj(ZER(gt2,j1))*Ye(j1,j2))));

   const std::complex<double> right = 0.7071067811865475*SUM(j1,0,2,Conj(ZE(gt3,j1))*ZEL(gt2,j1))*(0.7745966692414834*g1*ZN1(gt1,0) + g2*ZN1(gt1,1));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Fe, fields::hh, typename fields::bar<fields::Fe>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt3 = indices[1];
   const int gt1 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,j2)))*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZEL(gt1,j2))*ZH(gt3,0);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::hh, fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,j2)))*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZEL(gt1,j2))*ZH(gt3,0);

   return {left, right};
}

ChiralVertex VertexImpl<fields::Fe, typename fields::bar<fields::Chi>::type, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZN2 = MODELPARAMETER(ZN2);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZN1 = MODELPARAMETER(ZN1);

   const std::complex<double> left = -(Conj(ZN2(gt1,2))*SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1))));

   const std::complex<double> right = -1.0954451150103321*g1*SUM(j1,0,2,ZE(gt3,3 + j1)*ZER(gt2,j1))*ZN1(gt1,0);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Chi, fields::Se>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZN1 = MODELPARAMETER(ZN1);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZN2 = MODELPARAMETER(ZN2);

   const std::complex<double> left = -1.0954451150103321*g1*Conj(ZN1(gt2,0))*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*Conj(ZER(gt1,j1)));

   const std::complex<double> right = -(SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt3,3 + j1)))*ZEL(gt1,j2))*ZN2(gt2,2));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Fe, typename fields::bar<fields::Fv>::type, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const int gt3 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> left = 0;

   const std::complex<double> right = SUM(j1,0,2,Conj(Ye(j1,gt1))*ZER(gt2,j1))*ZP(gt3,0);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fv, fields::Hpm>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> left = SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,gt2))*ZP(gt3,0);

   const std::complex<double> right = 0;

   return {left, right};
}

ChiralVertex VertexImpl<fields::Fe, typename fields::conj<fields::Sv>::type, fields::Cha1>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt3 = indices[1];
   const int gt1 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZER = MODELPARAMETER(ZER);
   const auto UM1 = MODELPARAMETER(UM1);

   const std::complex<double> left = -(g2*Conj(UP1(gt1,0))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZV(gt3,j1)));

   const std::complex<double> right = SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZV(gt3,j2))*UM1(gt1,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Sv, typename fields::bar<fields::Cha1>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt3 = indices[1];
   const int gt1 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto UP1 = MODELPARAMETER(UP1);

   const std::complex<double> left = Conj(UM1(gt1,1))*SUM(j2,0,2,Conj(ZV(gt3,j2))*SUM(j1,0,2,Conj(ZER(gt2,j1))*Ye(j1,j2)));

   const std::complex<double> right = -(g2*SUM(j1,0,2,Conj(ZV(gt3,j1))*ZEL(gt2,j1))*UP1(gt1,0));

   return {left, right};
}

ChiralVertex VertexImpl<fields::VP, typename fields::bar<fields::Cha1>::type, fields::Cha1>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -(g2*Conj(UP1(gt2,0))*Sin(ThetaW)*UP1(gt1,0)) - 0.5*Conj(UP1(gt2,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UP1(gt1,1);

   const std::complex<double> right = -(g2*Conj(UM1(gt1,0))*Sin(ThetaW)*UM1(gt2,0)) - 0.5*Conj(UM1(gt1,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UM1(gt2,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fe, fields::VP>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Hpm, fields::Fv>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> left = SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,gt2))*ZP(gt3,0);

   const std::complex<double> right = 0;

   return {left, right};
}

MomentumDifferenceVertex VertexImpl<typename fields::conj<fields::Hpm>::type, fields::Hpm, fields::VP>::evaluate(
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

   const std::complex<double> result = 0.5*((0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*ZP(gt1,0)*ZP(gt2,0) + (0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*ZP(gt1,1)*ZP(gt2,1) + 2*g2*Sin(ThetaW)*(ZP(gt1,2)*ZP(gt2,2) + ZP(gt1,3)*ZP(gt2,3)));

   return {result, minuend_index, subtrahend_index};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fv>::type, typename fields::conj<fields::Hpm>::type, fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> left = 0;

   const std::complex<double> right = SUM(j1,0,2,Conj(Ye(j1,gt1))*ZER(gt2,j1))*ZP(gt3,0);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Se, fields::Chi>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZN1 = MODELPARAMETER(ZN1);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZN2 = MODELPARAMETER(ZN2);

   const std::complex<double> left = -1.0954451150103321*g1*Conj(ZN1(gt2,0))*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*Conj(ZER(gt1,j1)));

   const std::complex<double> right = -(SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt3,3 + j1)))*ZEL(gt1,j2))*ZN2(gt2,2));

   return {left, right};
}

MomentumDifferenceVertex VertexImpl<typename fields::conj<fields::Se>::type, fields::Se, fields::VP>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Chi>::type, typename fields::conj<fields::Se>::type, fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZN2 = MODELPARAMETER(ZN2);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZN1 = MODELPARAMETER(ZN1);

   const std::complex<double> left = -(Conj(ZN2(gt1,2))*SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1))));

   const std::complex<double> right = -1.0954451150103321*g1*SUM(j1,0,2,ZE(gt3,3 + j1)*ZER(gt2,j1))*ZN1(gt1,0);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Se, typename fields::bar<fields::Chi>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt3 = indices[1];
   const int gt1 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZN2 = MODELPARAMETER(ZN2);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZN1 = MODELPARAMETER(ZN1);

   const std::complex<double> left = -(Conj(ZN2(gt1,2))*SUM(j2,0,2,Conj(ZE(gt3,j2))*SUM(j1,0,2,Conj(ZER(gt2,j1))*Ye(j1,j2))));

   const std::complex<double> right = 0.7071067811865475*SUM(j1,0,2,Conj(ZE(gt3,j1))*ZEL(gt2,j1))*(0.7745966692414834*g1*ZN1(gt1,0) + g2*ZN1(gt1,1));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Chi, typename fields::conj<fields::Se>::type, fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZN1 = MODELPARAMETER(ZN1);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZN2 = MODELPARAMETER(ZN2);

   const std::complex<double> left = 0.7071067811865475*(0.7745966692414834*g1*Conj(ZN1(gt1,0)) + g2*Conj(ZN1(gt1,1)))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZE(gt3,j1));

   const std::complex<double> right = -(SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZE(gt3,j2))*ZN2(gt1,2));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, typename fields::bar<fields::Cha1>::type, fields::Sv>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto UP1 = MODELPARAMETER(UP1);

   const std::complex<double> left = Conj(UM1(gt1,1))*SUM(j2,0,2,Conj(ZV(gt3,j2))*SUM(j1,0,2,Conj(ZER(gt2,j1))*Ye(j1,j2)));

   const std::complex<double> right = -(g2*SUM(j1,0,2,Conj(ZV(gt3,j1))*ZEL(gt2,j1))*UP1(gt1,0));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Cha1, typename fields::bar<fields::Cha1>::type, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = g2*Conj(UM1(gt1,0))*Sin(ThetaW)*UM1(gt2,0) + 0.5*Conj(UM1(gt1,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UM1(gt2,1);

   const std::complex<double> right = g2*Conj(UP1(gt2,0))*Sin(ThetaW)*UP1(gt1,0) + 0.5*Conj(UP1(gt2,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UP1(gt1,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::conj<fields::Sv>::type, fields::Cha1, fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZER = MODELPARAMETER(ZER);
   const auto UM1 = MODELPARAMETER(UM1);

   const std::complex<double> left = -(g2*Conj(UP1(gt1,0))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZV(gt3,j1)));

   const std::complex<double> right = SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZV(gt3,j2))*UM1(gt1,1);

   return {left, right};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto vd = MODELPARAMETER(vd);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto vu = MODELPARAMETER(vu);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MuD = MODELPARAMETER(MuD);
   const auto vS = MODELPARAMETER(vS);
   const auto vT = MODELPARAMETER(vT);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto MuU = MODELPARAMETER(MuU);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.05*(-5*(vd*Conj(LamSD)*(1.4142135623730951*LamTD*ZA(gt1,3)*ZA(gt2,2) + ZA(gt1,2)*(4*LamSD*ZA(gt2,2) + 1.4142135623730951*LamTD*ZA(gt2,3)))*ZH(gt3,0) + vd*Conj(LamTD)*(1.4142135623730951*LamSD*ZA(gt1,2)*ZA(gt2,3) + ZA(gt1,3)*(1.4142135623730951*LamSD*ZA(gt2,2) + 2*LamTD*ZA(gt2,3)))*ZH(gt3,0) + vu*(-(Conj(LamTU)*(1.4142135623730951*LamSU*ZA(gt1,2)*ZA(gt2,3) + ZA(gt1,3)*(1.4142135623730951*LamSU*ZA(gt2,2) - 2*LamTU*ZA(gt2,3)))) - Conj(LamSU)*(1.4142135623730951*LamTU*ZA(gt1,3)*ZA(gt2,2) + ZA(gt1,2)*(-4*LamSU*ZA(gt2,2) + 1.4142135623730951*LamTU*ZA(gt2,3))))*ZH(gt3,1)) - ZA(gt1,0)*ZA(gt2,0)*(vd*(3*Sqr(g1) + 5*Sqr(g2))*ZH(gt3,0) - vu*(3*Sqr(g1) + 5*Sqr(g2))*ZH(gt3,1) - 7.745966692414834*g1*MDBS*ZH(gt3,2) + 20*vS*AbsSqr(LamSD)*ZH(gt3,2) + 14.142135623730951*MuD*Conj(LamSD)*ZH(gt3,2) + 7.0710678118654755*LamTD*vT*Conj(LamSD)*ZH(gt3,2) + 7.0710678118654755*LamSD*vT*Conj(LamTD)*ZH(gt3,2) - 7.745966692414834*g1*Conj(MDBS)*ZH(gt3,2) + 14.142135623730951*LamSD*Conj(MuD)*ZH(gt3,2) + 10*g2*MDWBT*ZH(gt3,3) + 10*vT*AbsSqr(LamTD)*ZH(gt3,3) + 7.0710678118654755*LamTD*vS*Conj(LamSD)*ZH(gt3,3) + 10*MuD*Conj(LamTD)*ZH(gt3,3) + 7.0710678118654755*LamSD*vS*Conj(LamTD)*ZH(gt3,3) + 10*g2*Conj(MDWBT)*ZH(gt3,3) + 10*LamTD*Conj(MuD)*ZH(gt3,3)) + ZA(gt1,1)*ZA(gt2,1)*(vd*(3*Sqr(g1) + 5*Sqr(g2))*ZH(gt3,0) - vu*(3*Sqr(g1) + 5*Sqr(g2))*ZH(gt3,1) - 7.745966692414834*g1*MDBS*ZH(gt3,2) - 20*vS*AbsSqr(LamSU)*ZH(gt3,2) - 14.142135623730951*MuU*Conj(LamSU)*ZH(gt3,2) + 7.0710678118654755*LamTU*vT*Conj(LamSU)*ZH(gt3,2) + 7.0710678118654755*LamSU*vT*Conj(LamTU)*ZH(gt3,2) - 7.745966692414834*g1*Conj(MDBS)*ZH(gt3,2) - 14.142135623730951*LamSU*Conj(MuU)*ZH(gt3,2) + 10*g2*MDWBT*ZH(gt3,3) - 10*vT*AbsSqr(LamTU)*ZH(gt3,3) + 7.0710678118654755*LamTU*vS*Conj(LamSU)*ZH(gt3,3) + 10*MuU*Conj(LamTU)*ZH(gt3,3) + 7.0710678118654755*LamSU*vS*Conj(LamTU)*ZH(gt3,3) + 10*g2*Conj(MDWBT)*ZH(gt3,3) + 10*LamTU*Conj(MuU)*ZH(gt3,3)));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Mu = MODELPARAMETER(Mu);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,0.5)*Mu*(1.4142135623730951*Conj(LamSD)*ZA(gt1,2)*ZP(gt2,1) + Conj(LamTD)*(-(ZA(gt1,3)*ZP(gt2,1)) + 1.4142135623730951*ZA(gt1,1)*ZP(gt2,2)));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Hpm, typename fields::conj<fields::SRum>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Mu = MODELPARAMETER(Mu);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,-0.5)*Mu*(1.4142135623730951*Conj(LamSU)*ZA(gt1,2)*ZP(gt2,0) + Conj(LamTU)*(ZA(gt1,3)*ZP(gt2,0) + 1.4142135623730951*ZA(gt1,0)*ZP(gt2,3)));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::hh, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto vd = MODELPARAMETER(vd);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto vS = MODELPARAMETER(vS);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto MuU = MODELPARAMETER(MuU);
   const auto vu = MODELPARAMETER(vu);
   const auto MuD = MODELPARAMETER(MuD);
   const auto g1 = MODELPARAMETER(g1);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto vT = MODELPARAMETER(vT);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = std::complex<double>(0,-0.25)*(ZA(gt1,3)*(1.4142135623730951*LamTD*vd*Conj(LamSD)*ZH(gt2,2)*ZH(gt3,0) - 1.4142135623730951*LamSD*vd*Conj(LamTD)*ZH(gt2,2)*ZH(gt3,0) - 2*g2*MDWBT*ZH(gt2,1)*ZH(gt3,1) - 1.4142135623730951*LamTU*vS*Conj(LamSU)*ZH(gt2,1)*ZH(gt3,1) + 2*MuU*Conj(LamTU)*ZH(gt2,1)*ZH(gt3,1) + 1.4142135623730951*LamSU*vS*Conj(LamTU)*ZH(gt2,1)*ZH(gt3,1) + 2*g2*Conj(MDWBT)*ZH(gt2,1)*ZH(gt3,1) - 2*LamTU*Conj(MuU)*ZH(gt2,1)*ZH(gt3,1) - 1.4142135623730951*LamTU*vu*Conj(LamSU)*ZH(gt2,2)*ZH(gt3,1) + 1.4142135623730951*LamSU*vu*Conj(LamTU)*ZH(gt2,2)*ZH(gt3,1) - 1.4142135623730951*LamTU*vu*Conj(LamSU)*ZH(gt2,1)*ZH(gt3,2) + 1.4142135623730951*LamSU*vu*Conj(LamTU)*ZH(gt2,1)*ZH(gt3,2) + ZH(gt2,0)*((2*g2*MDWBT + 1.4142135623730951*LamTD*vS*Conj(LamSD) - 2*MuD*Conj(LamTD) - 1.4142135623730951*LamSD*vS*Conj(LamTD) - 2*g2*Conj(MDWBT) + 2*LamTD*Conj(MuD))*ZH(gt3,0) + 1.4142135623730951*vd*(LamTD*Conj(LamSD) - LamSD*Conj(LamTD))*ZH(gt3,2))) + ZA(gt1,2)*(-1.4142135623730951*LamTD*vd*Conj(LamSD)*ZH(gt2,3)*ZH(gt3,0) + 1.4142135623730951*LamSD*vd*Conj(LamTD)*ZH(gt2,3)*ZH(gt3,0) + 1.5491933384829668*g1*MDBS*ZH(gt2,1)*ZH(gt3,1) - 2.8284271247461903*MuU*Conj(LamSU)*ZH(gt2,1)*ZH(gt3,1) + 1.4142135623730951*LamTU*vT*Conj(LamSU)*ZH(gt2,1)*ZH(gt3,1) - 1.4142135623730951*LamSU*vT*Conj(LamTU)*ZH(gt2,1)*ZH(gt3,1) - 1.5491933384829668*g1*Conj(MDBS)*ZH(gt2,1)*ZH(gt3,1) + 2.8284271247461903*LamSU*Conj(MuU)*ZH(gt2,1)*ZH(gt3,1) + 1.4142135623730951*LamTU*vu*Conj(LamSU)*ZH(gt2,3)*ZH(gt3,1) - 1.4142135623730951*LamSU*vu*Conj(LamTU)*ZH(gt2,3)*ZH(gt3,1) + 1.4142135623730951*LamTU*vu*Conj(LamSU)*ZH(gt2,1)*ZH(gt3,3) - 1.4142135623730951*LamSU*vu*Conj(LamTU)*ZH(gt2,1)*ZH(gt3,3) + ZH(gt2,0)*((-1.5491933384829668*g1*MDBS - 1.4142135623730951*(2*MuD + LamTD*vT)*Conj(LamSD) + 1.4142135623730951*LamSD*vT*Conj(LamTD) + 1.5491933384829668*g1*Conj(MDBS) + 2.8284271247461903*LamSD*Conj(MuD))*ZH(gt3,0) + 1.4142135623730951*vd*(-(LamTD*Conj(LamSD)) + LamSD*Conj(LamTD))*ZH(gt3,3))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Hpm, typename fields::conj<fields::SRum>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Mu = MODELPARAMETER(Mu);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.5*Mu*(1.4142135623730951*Conj(LamSU)*ZH(gt1,2)*ZP(gt2,0) + Conj(LamTU)*(ZH(gt1,3)*ZP(gt2,0) - 1.4142135623730951*ZH(gt1,0)*ZP(gt2,3)));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::hh, typename fields::conj<fields::Rh>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Mu = MODELPARAMETER(Mu);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = std::complex<double>(0,0.25)*Mu*(2*Conj(LamSD)*(ZA(gt1,2)*ZH(gt2,1) - ZA(gt1,1)*ZH(gt2,2))*ZHR(gt3,0) + 1.4142135623730951*Conj(LamTD)*(ZA(gt1,3)*ZH(gt2,1) - ZA(gt1,1)*ZH(gt2,3))*ZHR(gt3,0) + (Conj(LamSU)*(-2*ZA(gt1,2)*ZH(gt2,0) + 2*ZA(gt1,0)*ZH(gt2,2)) + 1.4142135623730951*Conj(LamTU)*(ZA(gt1,3)*ZH(gt2,0) - ZA(gt1,0)*ZH(gt2,3)))*ZHR(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::Rh, typename fields::conj<fields::SRum>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto vd = MODELPARAMETER(vd);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto g2 = MODELPARAMETER(g2);
   const auto vu = MODELPARAMETER(vu);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto vT = MODELPARAMETER(vT);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto vS = MODELPARAMETER(vS);
   const auto MuU = MODELPARAMETER(MuU);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-2.8284271247461903*LamSD*vd*Conj(LamSU)*ZHR(gt2,0)*ZP(gt1,1) - ZHR(gt2,1)*(1.4142135623730951*vd*Sqr(g2)*ZP(gt1,0) + 1.4142135623730951*vu*(-2*AbsSqr(LamSU) + Sqr(g2))*ZP(gt1,1) + 4*g2*MDWBT*ZP(gt1,2) - 2.8284271247461903*LamTU*vS*Conj(LamSU)*ZP(gt1,2) - 4*LamTU*Conj(MuU)*ZP(gt1,2) + 2*vT*Sqr(g2)*ZP(gt1,2) + 4*g2*Conj(MDWBT)*ZP(gt1,3) - 2*vT*Sqr(g2)*ZP(gt1,3)) - Conj(LamTU)*(1.4142135623730951*LamTD*ZHR(gt2,0)*(2*vu*ZP(gt1,0) + vd*ZP(gt1,1)) + ZHR(gt2,1)*(1.4142135623730951*LamTU*vu*ZP(gt1,1) - 2*(LamTU*vT*ZP(gt1,2) + (2*MuU + 1.4142135623730951*LamSU*vS - LamTU*vT)*ZP(gt1,3)))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Mu = MODELPARAMETER(Mu);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.5*Mu*(-1.4142135623730951*Conj(LamSD)*ZH(gt1,2)*ZP(gt2,1) + Conj(LamTD)*(ZH(gt1,3)*ZP(gt2,1) + 1.4142135623730951*ZH(gt1,1)*ZP(gt2,2)));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MuD = MODELPARAMETER(MuD);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto vS = MODELPARAMETER(vS);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto vT = MODELPARAMETER(vT);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto vd = MODELPARAMETER(vd);
   const auto MuU = MODELPARAMETER(MuU);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto vu = MODELPARAMETER(vu);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.25*(1.5491933384829668*g1*MDBS*ZH(gt1,2)*ZH(gt2,0)*ZH(gt3,0) - 4*vS*AbsSqr(LamSD)*ZH(gt1,2)*ZH(gt2,0)*ZH(gt3,0) - 2.8284271247461903*MuD*Conj(LamSD)*ZH(gt1,2)*ZH(gt2,0)*ZH(gt3,0) - 1.4142135623730951*LamTD*vT*Conj(LamSD)*ZH(gt1,2)*ZH(gt2,0)*ZH(gt3,0) - 1.4142135623730951*LamSD*vT*Conj(LamTD)*ZH(gt1,2)*ZH(gt2,0)*ZH(gt3,0) + 1.5491933384829668*g1*Conj(MDBS)*ZH(gt1,2)*ZH(gt2,0)*ZH(gt3,0) - 2.8284271247461903*LamSD*Conj(MuD)*ZH(gt1,2)*ZH(gt2,0)*ZH(gt3,0) - 2*g2*MDWBT*ZH(gt1,3)*ZH(gt2,0)*ZH(gt3,0) - 2*vT*AbsSqr(LamTD)*ZH(gt1,3)*ZH(gt2,0)*ZH(gt3,0) - 1.4142135623730951*LamTD*vS*Conj(LamSD)*ZH(gt1,3)*ZH(gt2,0)*ZH(gt3,0) - 2*MuD*Conj(LamTD)*ZH(gt1,3)*ZH(gt2,0)*ZH(gt3,0) - 1.4142135623730951*LamSD*vS*Conj(LamTD)*ZH(gt1,3)*ZH(gt2,0)*ZH(gt3,0) - 2*g2*Conj(MDWBT)*ZH(gt1,3)*ZH(gt2,0)*ZH(gt3,0) - 2*LamTD*Conj(MuD)*ZH(gt1,3)*ZH(gt2,0)*ZH(gt3,0) - 4*vd*AbsSqr(LamSD)*ZH(gt1,2)*ZH(gt2,2)*ZH(gt3,0) - 1.4142135623730951*LamTD*vd*Conj(LamSD)*ZH(gt1,3)*ZH(gt2,2)*ZH(gt3,0) - 1.4142135623730951*LamSD*vd*Conj(LamTD)*ZH(gt1,3)*ZH(gt2,2)*ZH(gt3,0) - 1.4142135623730951*LamTD*vd*Conj(LamSD)*ZH(gt1,2)*ZH(gt2,3)*ZH(gt3,0) - 1.4142135623730951*LamSD*vd*Conj(LamTD)*ZH(gt1,2)*ZH(gt2,3)*ZH(gt3,0) - 2*vd*AbsSqr(LamTD)*ZH(gt1,3)*ZH(gt2,3)*ZH(gt3,0) - 1.5491933384829668*g1*MDBS*ZH(gt1,2)*ZH(gt2,1)*ZH(gt3,1) - 4*vS*AbsSqr(LamSU)*ZH(gt1,2)*ZH(gt2,1)*ZH(gt3,1) - 2.8284271247461903*MuU*Conj(LamSU)*ZH(gt1,2)*ZH(gt2,1)*ZH(gt3,1) + 1.4142135623730951*LamTU*vT*Conj(LamSU)*ZH(gt1,2)*ZH(gt2,1)*ZH(gt3,1) + 1.4142135623730951*LamSU*vT*Conj(LamTU)*ZH(gt1,2)*ZH(gt2,1)*ZH(gt3,1) - 1.5491933384829668*g1*Conj(MDBS)*ZH(gt1,2)*ZH(gt2,1)*ZH(gt3,1) - 2.8284271247461903*LamSU*Conj(MuU)*ZH(gt1,2)*ZH(gt2,1)*ZH(gt3,1) + 2*g2*MDWBT*ZH(gt1,3)*ZH(gt2,1)*ZH(gt3,1) - 2*vT*AbsSqr(LamTU)*ZH(gt1,3)*ZH(gt2,1)*ZH(gt3,1) + 1.4142135623730951*LamTU*vS*Conj(LamSU)*ZH(gt1,3)*ZH(gt2,1)*ZH(gt3,1) + 2*MuU*Conj(LamTU)*ZH(gt1,3)*ZH(gt2,1)*ZH(gt3,1) + 1.4142135623730951*LamSU*vS*Conj(LamTU)*ZH(gt1,3)*ZH(gt2,1)*ZH(gt3,1) + 2*g2*Conj(MDWBT)*ZH(gt1,3)*ZH(gt2,1)*ZH(gt3,1) + 2*LamTU*Conj(MuU)*ZH(gt1,3)*ZH(gt2,1)*ZH(gt3,1) - 4*vu*AbsSqr(LamSU)*ZH(gt1,2)*ZH(gt2,2)*ZH(gt3,1) + 1.4142135623730951*LamTU*vu*Conj(LamSU)*ZH(gt1,3)*ZH(gt2,2)*ZH(gt3,1) + 1.4142135623730951*LamSU*vu*Conj(LamTU)*ZH(gt1,3)*ZH(gt2,2)*ZH(gt3,1) + 1.4142135623730951*LamTU*vu*Conj(LamSU)*ZH(gt1,2)*ZH(gt2,3)*ZH(gt3,1) + 1.4142135623730951*LamSU*vu*Conj(LamTU)*ZH(gt1,2)*ZH(gt2,3)*ZH(gt3,1) - 2*vu*AbsSqr(LamTU)*ZH(gt1,3)*ZH(gt2,3)*ZH(gt3,1) - 4*vd*AbsSqr(LamSD)*ZH(gt1,2)*ZH(gt2,0)*ZH(gt3,2) - 1.4142135623730951*LamTD*vd*Conj(LamSD)*ZH(gt1,3)*ZH(gt2,0)*ZH(gt3,2) - 1.4142135623730951*LamSD*vd*Conj(LamTD)*ZH(gt1,3)*ZH(gt2,0)*ZH(gt3,2) - 4*vu*AbsSqr(LamSU)*ZH(gt1,2)*ZH(gt2,1)*ZH(gt3,2) + 1.4142135623730951*LamTU*vu*Conj(LamSU)*ZH(gt1,3)*ZH(gt2,1)*ZH(gt3,2) + 1.4142135623730951*LamSU*vu*Conj(LamTU)*ZH(gt1,3)*ZH(gt2,1)*ZH(gt3,2) - 1.4142135623730951*LamTD*vd*Conj(LamSD)*ZH(gt1,2)*ZH(gt2,0)*ZH(gt3,3) - 1.4142135623730951*LamSD*vd*Conj(LamTD)*ZH(gt1,2)*ZH(gt2,0)*ZH(gt3,3) - 2*vd*AbsSqr(LamTD)*ZH(gt1,3)*ZH(gt2,0)*ZH(gt3,3) + 1.4142135623730951*LamTU*vu*Conj(LamSU)*ZH(gt1,2)*ZH(gt2,1)*ZH(gt3,3) + 1.4142135623730951*LamSU*vu*Conj(LamTU)*ZH(gt1,2)*ZH(gt2,1)*ZH(gt3,3) - 2*vu*AbsSqr(LamTU)*ZH(gt1,3)*ZH(gt2,1)*ZH(gt3,3) - ZH(gt1,0)*(-1.5491933384829668*g1*MDBS*ZH(gt2,2)*ZH(gt3,0) + 4*vS*AbsSqr(LamSD)*ZH(gt2,2)*ZH(gt3,0) + 2.8284271247461903*MuD*Conj(LamSD)*ZH(gt2,2)*ZH(gt3,0) + 1.4142135623730951*LamTD*vT*Conj(LamSD)*ZH(gt2,2)*ZH(gt3,0) + 1.4142135623730951*LamSD*vT*Conj(LamTD)*ZH(gt2,2)*ZH(gt3,0) - 1.5491933384829668*g1*Conj(MDBS)*ZH(gt2,2)*ZH(gt3,0) + 2.8284271247461903*LamSD*Conj(MuD)*ZH(gt2,2)*ZH(gt3,0) + 2*g2*MDWBT*ZH(gt2,3)*ZH(gt3,0) + 2*vT*AbsSqr(LamTD)*ZH(gt2,3)*ZH(gt3,0) + 1.4142135623730951*LamTD*vS*Conj(LamSD)*ZH(gt2,3)*ZH(gt3,0) + 2*MuD*Conj(LamTD)*ZH(gt2,3)*ZH(gt3,0) + 1.4142135623730951*LamSD*vS*Conj(LamTD)*ZH(gt2,3)*ZH(gt3,0) + 2*g2*Conj(MDWBT)*ZH(gt2,3)*ZH(gt3,0) + 2*LamTD*Conj(MuD)*ZH(gt2,3)*ZH(gt3,0) - (0.6*Sqr(g1) + Sqr(g2))*ZH(gt2,1)*(vu*ZH(gt3,0) + vd*ZH(gt3,1)) + 4*vd*AbsSqr(LamSD)*ZH(gt2,2)*ZH(gt3,2) + 1.4142135623730951*LamTD*vd*Conj(LamSD)*ZH(gt2,3)*ZH(gt3,2) + 1.4142135623730951*LamSD*vd*Conj(LamTD)*ZH(gt2,3)*ZH(gt3,2) + 1.4142135623730951*LamTD*vd*Conj(LamSD)*ZH(gt2,2)*ZH(gt3,3) + 1.4142135623730951*LamSD*vd*Conj(LamTD)*ZH(gt2,2)*ZH(gt3,3) + 2*vd*AbsSqr(LamTD)*ZH(gt2,3)*ZH(gt3,3) + ZH(gt2,0)*(3*vd*(0.6*Sqr(g1) + Sqr(g2))*ZH(gt3,0) - vu*(0.6*Sqr(g1) + Sqr(g2))*ZH(gt3,1) - 1.5491933384829668*g1*MDBS*ZH(gt3,2) + 4*vS*AbsSqr(LamSD)*ZH(gt3,2) + 2.8284271247461903*MuD*Conj(LamSD)*ZH(gt3,2) + 1.4142135623730951*LamTD*vT*Conj(LamSD)*ZH(gt3,2) + 1.4142135623730951*LamSD*vT*Conj(LamTD)*ZH(gt3,2) - 1.5491933384829668*g1*Conj(MDBS)*ZH(gt3,2) + 2.8284271247461903*LamSD*Conj(MuD)*ZH(gt3,2) + 2*g2*MDWBT*ZH(gt3,3) + 2*vT*AbsSqr(LamTD)*ZH(gt3,3) + 1.4142135623730951*LamTD*vS*Conj(LamSD)*ZH(gt3,3) + 2*MuD*Conj(LamTD)*ZH(gt3,3) + 1.4142135623730951*LamSD*vS*Conj(LamTD)*ZH(gt3,3) + 2*g2*Conj(MDWBT)*ZH(gt3,3) + 2*LamTD*Conj(MuD)*ZH(gt3,3))) + ZH(gt1,1)*(-1.5491933384829668*g1*MDBS*ZH(gt2,2)*ZH(gt3,1) - 4*vS*AbsSqr(LamSU)*ZH(gt2,2)*ZH(gt3,1) - 2.8284271247461903*MuU*Conj(LamSU)*ZH(gt2,2)*ZH(gt3,1) + 1.4142135623730951*LamTU*vT*Conj(LamSU)*ZH(gt2,2)*ZH(gt3,1) + 1.4142135623730951*LamSU*vT*Conj(LamTU)*ZH(gt2,2)*ZH(gt3,1) - 1.5491933384829668*g1*Conj(MDBS)*ZH(gt2,2)*ZH(gt3,1) - 2.8284271247461903*LamSU*Conj(MuU)*ZH(gt2,2)*ZH(gt3,1) + 2*g2*MDWBT*ZH(gt2,3)*ZH(gt3,1) - 2*vT*AbsSqr(LamTU)*ZH(gt2,3)*ZH(gt3,1) + 1.4142135623730951*LamTU*vS*Conj(LamSU)*ZH(gt2,3)*ZH(gt3,1) + 2*MuU*Conj(LamTU)*ZH(gt2,3)*ZH(gt3,1) + 1.4142135623730951*LamSU*vS*Conj(LamTU)*ZH(gt2,3)*ZH(gt3,1) + 2*g2*Conj(MDWBT)*ZH(gt2,3)*ZH(gt3,1) + 2*LamTU*Conj(MuU)*ZH(gt2,3)*ZH(gt3,1) + (0.6*Sqr(g1) + Sqr(g2))*ZH(gt2,0)*(vu*ZH(gt3,0) + vd*ZH(gt3,1)) - 4*vu*AbsSqr(LamSU)*ZH(gt2,2)*ZH(gt3,2) + 1.4142135623730951*LamTU*vu*Conj(LamSU)*ZH(gt2,3)*ZH(gt3,2) + 1.4142135623730951*LamSU*vu*Conj(LamTU)*ZH(gt2,3)*ZH(gt3,2) + 1.4142135623730951*LamTU*vu*Conj(LamSU)*ZH(gt2,2)*ZH(gt3,3) + 1.4142135623730951*LamSU*vu*Conj(LamTU)*ZH(gt2,2)*ZH(gt3,3) - 2*vu*AbsSqr(LamTU)*ZH(gt2,3)*ZH(gt3,3) + ZH(gt2,1)*(vd*(0.6*Sqr(g1) + Sqr(g2))*ZH(gt3,0) - 3*vu*(0.6*Sqr(g1) + Sqr(g2))*ZH(gt3,1) - 1.5491933384829668*g1*MDBS*ZH(gt3,2) - 4*vS*AbsSqr(LamSU)*ZH(gt3,2) - 2.8284271247461903*MuU*Conj(LamSU)*ZH(gt3,2) + 1.4142135623730951*LamTU*vT*Conj(LamSU)*ZH(gt3,2) + 1.4142135623730951*LamSU*vT*Conj(LamTU)*ZH(gt3,2) - 1.5491933384829668*g1*Conj(MDBS)*ZH(gt3,2) - 2.8284271247461903*LamSU*Conj(MuU)*ZH(gt3,2) + 2*g2*MDWBT*ZH(gt3,3) - 2*vT*AbsSqr(LamTU)*ZH(gt3,3) + 1.4142135623730951*LamTU*vS*Conj(LamSU)*ZH(gt3,3) + 2*MuU*Conj(LamTU)*ZH(gt3,3) + 1.4142135623730951*LamSU*vS*Conj(LamTU)*ZH(gt3,3) + 2*g2*Conj(MDWBT)*ZH(gt3,3) + 2*LamTU*Conj(MuU)*ZH(gt3,3))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, typename fields::conj<fields::Rh>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Mu = MODELPARAMETER(Mu);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.25*Mu*(2*Conj(LamSD)*(ZH(gt1,2)*ZH(gt2,1) + ZH(gt1,1)*ZH(gt2,2))*ZHR(gt3,0) + 1.4142135623730951*Conj(LamTD)*(ZH(gt1,3)*ZH(gt2,1) + ZH(gt1,1)*ZH(gt2,3))*ZHR(gt3,0) + (-2*Conj(LamSU)*(ZH(gt1,2)*ZH(gt2,0) + ZH(gt1,0)*ZH(gt2,2)) + 1.4142135623730951*Conj(LamTU)*(ZH(gt1,3)*ZH(gt2,0) + ZH(gt1,0)*ZH(gt2,3)))*ZHR(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Hpm, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MuD = MODELPARAMETER(MuD);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto vS = MODELPARAMETER(vS);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto vT = MODELPARAMETER(vT);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto vd = MODELPARAMETER(vd);
   const auto MuU = MODELPARAMETER(MuU);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto vu = MODELPARAMETER(vu);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.05*(7.745966692414834*g1*MDBS*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) - 20*vS*AbsSqr(LamSD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) - 14.142135623730951*MuD*Conj(LamSD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) + 7.0710678118654755*LamTD*vT*Conj(LamSD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) + 7.0710678118654755*LamSD*vT*Conj(LamTD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) + 7.745966692414834*g1*Conj(MDBS)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) - 14.142135623730951*LamSD*Conj(MuD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) + 10*g2*MDWBT*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) - 10*vT*AbsSqr(LamTD)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) + 7.0710678118654755*LamTD*vS*Conj(LamSD)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) + 10*MuD*Conj(LamTD)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) + 7.0710678118654755*LamSD*vS*Conj(LamTD)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) + 10*g2*Conj(MDWBT)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) + 10*LamTD*Conj(MuD)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) - 10*LamTD*vd*Conj(LamSD)*ZH(gt1,2)*ZP(gt2,2)*ZP(gt3,0) + 7.0710678118654755*vd*AbsSqr(LamTD)*ZH(gt1,3)*ZP(gt2,2)*ZP(gt3,0) - 7.0710678118654755*vd*Sqr(g2)*ZH(gt1,3)*ZP(gt2,2)*ZP(gt3,0) - 10*LamSD*vd*Conj(LamTD)*ZH(gt1,2)*ZP(gt2,3)*ZP(gt3,0) - 7.0710678118654755*vd*AbsSqr(LamTD)*ZH(gt1,3)*ZP(gt2,3)*ZP(gt3,0) + 7.0710678118654755*vd*Sqr(g2)*ZH(gt1,3)*ZP(gt2,3)*ZP(gt3,0) - 7.745966692414834*g1*MDBS*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) - 20*vS*AbsSqr(LamSU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) - 14.142135623730951*MuU*Conj(LamSU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) - 7.0710678118654755*LamTU*vT*Conj(LamSU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) - 7.0710678118654755*LamSU*vT*Conj(LamTU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) - 7.745966692414834*g1*Conj(MDBS)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) - 14.142135623730951*LamSU*Conj(MuU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) - 10*g2*MDWBT*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,1) - 10*vT*AbsSqr(LamTU)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,1) - 7.0710678118654755*LamTU*vS*Conj(LamSU)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,1) - 10*MuU*Conj(LamTU)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,1) - 7.0710678118654755*LamSU*vS*Conj(LamTU)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,1) - 10*g2*Conj(MDWBT)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,1) - 10*LamTU*Conj(MuU)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,1) - 10*LamTU*vu*Conj(LamSU)*ZH(gt1,2)*ZP(gt2,2)*ZP(gt3,1) + 7.0710678118654755*vu*AbsSqr(LamTU)*ZH(gt1,3)*ZP(gt2,2)*ZP(gt3,1) - 7.0710678118654755*vu*Sqr(g2)*ZH(gt1,3)*ZP(gt2,2)*ZP(gt3,1) - 10*LamSU*vu*Conj(LamTU)*ZH(gt1,2)*ZP(gt2,3)*ZP(gt3,1) - 7.0710678118654755*vu*AbsSqr(LamTU)*ZH(gt1,3)*ZP(gt2,3)*ZP(gt3,1) + 7.0710678118654755*vu*Sqr(g2)*ZH(gt1,3)*ZP(gt2,3)*ZP(gt3,1) - 10*LamSD*vd*Conj(LamTD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,2) + 7.0710678118654755*vd*AbsSqr(LamTD)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,2) - 7.0710678118654755*vd*Sqr(g2)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,2) - 10*LamSU*vu*Conj(LamTU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,2) + 7.0710678118654755*vu*AbsSqr(LamTU)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,2) - 7.0710678118654755*vu*Sqr(g2)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,2) - 20*vT*Sqr(g2)*ZH(gt1,3)*ZP(gt2,2)*ZP(gt3,2) + 20*vT*Sqr(g2)*ZH(gt1,3)*ZP(gt2,3)*ZP(gt3,2) - 10*LamTD*vd*Conj(LamSD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,3) - 7.0710678118654755*vd*AbsSqr(LamTD)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,3) + 7.0710678118654755*vd*Sqr(g2)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,3) - 10*LamTU*vu*Conj(LamSU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,3) - 7.0710678118654755*vu*AbsSqr(LamTU)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,3) + 7.0710678118654755*vu*Sqr(g2)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,3) + 20*vT*Sqr(g2)*ZH(gt1,3)*ZP(gt2,2)*ZP(gt3,3) - 20*vT*Sqr(g2)*ZH(gt1,3)*ZP(gt2,3)*ZP(gt3,3) - ZH(gt1,1)*(ZP(gt2,0)*((-3*vu*Sqr(g1) + 5*vu*Sqr(g2))*ZP(gt3,0) + 5*vd*Sqr(g2)*ZP(gt3,1)) + 5*ZP(gt2,2)*((2*LamTU*vS*Conj(LamSU) + 1.4142135623730951*(2*g2*MDWBT - vT*AbsSqr(LamTU) + 2*LamTU*Conj(MuU) + vT*Sqr(g2)))*ZP(gt3,1) + 2*vu*Sqr(g2)*ZP(gt3,2)) + 5*ZP(gt2,3)*(((2.8284271247461903*MuU + 2*LamSU*vS + 1.4142135623730951*LamTU*vT)*Conj(LamTU) + 1.4142135623730951*g2*(-(g2*vT) + 2*Conj(MDWBT)))*ZP(gt3,1) - 2*vu*(-2*AbsSqr(LamTU) + Sqr(g2))*ZP(gt3,3)) + ZP(gt2,1)*(5*vd*Sqr(g2)*ZP(gt3,0) + vu*(3*Sqr(g1) + 5*Sqr(g2))*ZP(gt3,1) + 5*((2.8284271247461903*MuU + 2*LamSU*vS - 1.4142135623730951*LamTU*vT)*Conj(LamTU) + 1.4142135623730951*g2*(g2*vT + 2*Conj(MDWBT)))*ZP(gt3,2) + 5*(2*LamTU*vS*Conj(LamSU) + 1.4142135623730951*(2*g2*MDWBT + vT*AbsSqr(LamTU) + 2*LamTU*Conj(MuU) - vT*Sqr(g2)))*ZP(gt3,3))) - ZH(gt1,0)*(ZP(gt2,1)*(5*vu*Sqr(g2)*ZP(gt3,0) + vd*(-3*Sqr(g1) + 5*Sqr(g2))*ZP(gt3,1)) + 5*(ZP(gt2,2)*((2*LamTD*vS*Conj(LamSD) + 1.4142135623730951*(2*g2*MDWBT - vT*AbsSqr(LamTD) + 2*LamTD*Conj(MuD) + vT*Sqr(g2)))*ZP(gt3,0) - 2*vd*(-2*AbsSqr(LamTD) + Sqr(g2))*ZP(gt3,2)) + ZP(gt2,3)*(((2.8284271247461903*MuD + 2*LamSD*vS + 1.4142135623730951*LamTD*vT)*Conj(LamTD) + 1.4142135623730951*g2*(-(g2*vT) + 2*Conj(MDWBT)))*ZP(gt3,0) + 2*vd*Sqr(g2)*ZP(gt3,3))) + ZP(gt2,0)*(vd*(3*Sqr(g1) + 5*Sqr(g2))*ZP(gt3,0) + 5*(vu*Sqr(g2)*ZP(gt3,1) + ((2.8284271247461903*MuD + 2*LamSD*vS - 1.4142135623730951*LamTD*vT)*Conj(LamTD) + 1.4142135623730951*g2*(g2*vT + 2*Conj(MDWBT)))*ZP(gt3,2) + (2*LamTD*vS*Conj(LamSD) + 1.4142135623730951*(2*g2*MDWBT + vT*AbsSqr(LamTD) + 2*LamTD*Conj(MuD) - vT*Sqr(g2)))*ZP(gt3,3)))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::SRum, typename fields::conj<fields::SRum>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MuU = MODELPARAMETER(MuU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto vT = MODELPARAMETER(vT);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto vS = MODELPARAMETER(vS);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,-0.25)*((-1.5491933384829668*g1*MDBS - 1.4142135623730951*(2*MuU + LamTU*vT)*Conj(LamSU) + 1.4142135623730951*LamSU*vT*Conj(LamTU) + 1.5491933384829668*g1*Conj(MDBS) + 2.8284271247461903*LamSU*Conj(MuU))*ZA(gt1,2) + (-2*g2*MDWBT + 1.4142135623730951*LamTU*vS*Conj(LamSU) - (2*MuU + 1.4142135623730951*LamSU*vS)*Conj(LamTU) + 2*g2*Conj(MDWBT) + 2*LamTU*Conj(MuU))*ZA(gt1,3));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::SRum, typename fields::conj<fields::SRum>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MuU = MODELPARAMETER(MuU);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto vS = MODELPARAMETER(vS);
   const auto vT = MODELPARAMETER(vT);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.25*(-0.2*vd*(3*Sqr(g1) - 5*Sqr(g2))*ZH(gt1,0) - vu*(4*AbsSqr(LamTU) - 0.6*Sqr(g1) + Sqr(g2))*ZH(gt1,1) + 1.5491933384829668*g1*MDBS*ZH(gt1,2) - 4*vS*AbsSqr(LamSU)*ZH(gt1,2) - 2.8284271247461903*MuU*Conj(LamSU)*ZH(gt1,2) - 1.4142135623730951*LamTU*vT*Conj(LamSU)*ZH(gt1,2) - 1.4142135623730951*LamSU*vT*Conj(LamTU)*ZH(gt1,2) + 1.5491933384829668*g1*Conj(MDBS)*ZH(gt1,2) - 2.8284271247461903*LamSU*Conj(MuU)*ZH(gt1,2) + 2*g2*MDWBT*ZH(gt1,3) - 2*vT*AbsSqr(LamTU)*ZH(gt1,3) - 1.4142135623730951*LamTU*vS*Conj(LamSU)*ZH(gt1,3) - 2*MuU*Conj(LamTU)*ZH(gt1,3) - 1.4142135623730951*LamSU*vS*Conj(LamTU)*ZH(gt1,3) + 2*g2*Conj(MDWBT)*ZH(gt1,3) - 2*LamTU*Conj(MuU)*ZH(gt1,3));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Sd, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto vd = MODELPARAMETER(vd);
   const auto Mu = MODELPARAMETER(Mu);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.016666666666666666*(30*(-2*vd*SUM(j3,0,2,Conj(ZD(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,j1))*ZD(gt3,3 + j2)))*ZH(gt1,0) - 2*vd*SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gt3,j3))*ZH(gt1,0) + 1.4142135623730951*(Conj(Mu)*SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1))) + Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1)))*ZD(gt3,j2)))*ZH(gt1,1)) + 2*g1*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*(3*g1*vd*ZH(gt1,0) - 3*g1*vu*ZH(gt1,1) - 7.745966692414834*(MDBS + Conj(MDBS))*ZH(gt1,2)) + SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*(3*vd*(Sqr(g1) + 5*Sqr(g2))*ZH(gt1,0) - 3*vu*(Sqr(g1) + 5*Sqr(g2))*ZH(gt1,1) - 7.745966692414834*g1*MDBS*ZH(gt1,2) - 7.745966692414834*g1*Conj(MDBS)*ZH(gt1,2) + 30*g2*MDWBT*ZH(gt1,3) + 30*g2*Conj(MDWBT)*ZH(gt1,3)));

   return {result};
}

ScalarVertex VertexImpl<fields::Su, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto vS = MODELPARAMETER(vS);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto vT = MODELPARAMETER(vT);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto MuD = MODELPARAMETER(MuD);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = 0.5*(1.4142135623730951*vS*Conj(LamSD) - vT*Conj(LamTD) + 2*Conj(MuD))*SUM(j2,0,2,Conj(ZU(gt1,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt2,3 + j1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Sd, typename fields::conj<fields::SRum>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const auto vS = MODELPARAMETER(vS);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto vT = MODELPARAMETER(vT);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto MuU = MODELPARAMETER(MuU);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);

   const std::complex<double> result = -0.5*(1.4142135623730951*vS*Conj(LamSU) + vT*Conj(LamTU) + 2*Conj(MuU))*SUM(j2,0,2,Conj(ZD(gt1,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt3,3 + j1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Rh, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto vu = MODELPARAMETER(vu);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto vd = MODELPARAMETER(vd);
   const auto g2 = MODELPARAMETER(g2);
   const auto vT = MODELPARAMETER(vT);
   const auto MuD = MODELPARAMETER(MuD);
   const auto vS = MODELPARAMETER(vS);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-1.4142135623730951*ZHR(gt1,1)*(2*LamSU*vu*Conj(LamSD)*ZP(gt2,0) + LamTU*Conj(LamTD)*(vu*ZP(gt2,0) + 2*vd*ZP(gt2,1))) - ZHR(gt1,0)*(1.4142135623730951*vd*(-2*AbsSqr(LamSD) + AbsSqr(LamTD) + Sqr(g2))*ZP(gt2,0) + 1.4142135623730951*vu*Sqr(g2)*ZP(gt2,1) - 2*vT*AbsSqr(LamTD)*ZP(gt2,2) - 4*MuD*Conj(LamTD)*ZP(gt2,2) - 2.8284271247461903*LamSD*vS*Conj(LamTD)*ZP(gt2,2) + 4*g2*Conj(MDWBT)*ZP(gt2,2) + 2*vT*Sqr(g2)*ZP(gt2,2) + 4*g2*MDWBT*ZP(gt2,3) + 2*vT*AbsSqr(LamTD)*ZP(gt2,3) - 2.8284271247461903*LamTD*vS*Conj(LamSD)*ZP(gt2,3) - 4*LamTD*Conj(MuD)*ZP(gt2,3) - 2*vT*Sqr(g2)*ZP(gt2,3)));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::SRdp, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MuD = MODELPARAMETER(MuD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto vT = MODELPARAMETER(vT);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto vS = MODELPARAMETER(vS);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,-0.25)*((1.5491933384829668*g1*MDBS + 1.4142135623730951*(-2*MuD + LamTD*vT)*Conj(LamSD) - 1.4142135623730951*LamSD*vT*Conj(LamTD) - 1.5491933384829668*g1*Conj(MDBS) + 2.8284271247461903*LamSD*Conj(MuD))*ZA(gt1,2) + (2*g2*MDWBT - 1.4142135623730951*LamTD*vS*Conj(LamSD) + (2*MuD + 1.4142135623730951*LamSD*vS)*Conj(LamTD) - 2*g2*Conj(MDWBT) - 2*LamTD*Conj(MuD))*ZA(gt1,3));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::SRdp, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto vu = MODELPARAMETER(vu);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MuD = MODELPARAMETER(MuD);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto vS = MODELPARAMETER(vS);
   const auto vT = MODELPARAMETER(vT);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.25*(vd*(-4*AbsSqr(LamTD) + 0.6*Sqr(g1) - Sqr(g2))*ZH(gt1,0) + vu*(-0.6*Sqr(g1) + Sqr(g2))*ZH(gt1,1) - 1.5491933384829668*g1*MDBS*ZH(gt1,2) - 4*vS*AbsSqr(LamSD)*ZH(gt1,2) - 2.8284271247461903*MuD*Conj(LamSD)*ZH(gt1,2) + 1.4142135623730951*LamTD*vT*Conj(LamSD)*ZH(gt1,2) + 1.4142135623730951*LamSD*vT*Conj(LamTD)*ZH(gt1,2) - 1.5491933384829668*g1*Conj(MDBS)*ZH(gt1,2) - 2.8284271247461903*LamSD*Conj(MuD)*ZH(gt1,2) - 2*g2*MDWBT*ZH(gt1,3) - 2*vT*AbsSqr(LamTD)*ZH(gt1,3) + 1.4142135623730951*LamTD*vS*Conj(LamSD)*ZH(gt1,3) + 2*MuD*Conj(LamTD)*ZH(gt1,3) + 1.4142135623730951*LamSD*vS*Conj(LamTD)*ZH(gt1,3) - 2*g2*Conj(MDWBT)*ZH(gt1,3) + 2*LamTD*Conj(MuD)*ZH(gt1,3));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Su, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Mu = MODELPARAMETER(Mu);
   const auto vu = MODELPARAMETER(vu);
   const auto g1 = MODELPARAMETER(g1);
   const auto vd = MODELPARAMETER(vd);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.016666666666666666*(30*(1.4142135623730951*Conj(Mu)*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt3,3 + j1)))*ZH(gt1,0) + 1.4142135623730951*Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1)))*ZU(gt3,j2))*ZH(gt1,0) - 2*vu*(SUM(j3,0,2,Conj(ZU(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))*ZU(gt3,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gt3,j3)))*ZH(gt1,1)) + 4*g1*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*(-3*g1*vd*ZH(gt1,0) + 3*g1*vu*ZH(gt1,1) + 7.745966692414834*(MDBS + Conj(MDBS))*ZH(gt1,2)) + SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*(3*vd*(Sqr(g1) - 5*Sqr(g2))*ZH(gt1,0) - 3*vu*(Sqr(g1) - 5*Sqr(g2))*ZH(gt1,1) - 7.745966692414834*g1*(MDBS + Conj(MDBS))*ZH(gt1,2) - 30*g2*(MDWBT + Conj(MDWBT))*ZH(gt1,3)));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::Su, typename fields::conj<fields::SRum>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt4 = indices[2];
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = Conj(LamTU)*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt4,3 + j1)))*ZP(gt1,3);

   return {result};
}

ScalarVertex VertexImpl<fields::Sd, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -(Conj(LamTD)*SUM(j2,0,2,Conj(ZD(gt1,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1)))*ZP(gt2,2));

   return {result};
}

ScalarVertex VertexImpl<fields::Se, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Se>::type, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -(Conj(LamTD)*SUM(j2,0,2,Conj(ZE(gt1,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)))*ZP(gt2,2));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Sd, typename fields::conj<fields::SRum>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt4 = indices[2];
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = -0.5*SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt4,3 + j1)))*(1.4142135623730951*Conj(LamSU)*ZH(gt1,2) + Conj(LamTU)*ZH(gt1,3));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Su, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.5*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1)))*(1.4142135623730951*Conj(LamSD)*ZH(gt1,2) - Conj(LamTD)*ZH(gt1,3));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Hpm, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto vu = MODELPARAMETER(vu);
   const auto vd = MODELPARAMETER(vd);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto vT = MODELPARAMETER(vT);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto vS = MODELPARAMETER(vS);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto MuD = MODELPARAMETER(MuD);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto MuU = MODELPARAMETER(MuU);
   const auto g1 = MODELPARAMETER(g1);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,0.05)*(ZA(gt1,2)*(10*LamTD*vd*Conj(LamSD)*ZP(gt2,2)*ZP(gt3,0) - 10*LamSD*vd*Conj(LamTD)*ZP(gt2,3)*ZP(gt3,0) - 7.745966692414834*g1*MDBS*ZP(gt2,1)*ZP(gt3,1) + 14.142135623730951*MuU*Conj(LamSU)*ZP(gt2,1)*ZP(gt3,1) + 7.0710678118654755*LamTU*vT*Conj(LamSU)*ZP(gt2,1)*ZP(gt3,1) - 7.0710678118654755*LamSU*vT*Conj(LamTU)*ZP(gt2,1)*ZP(gt3,1) + 7.745966692414834*g1*Conj(MDBS)*ZP(gt2,1)*ZP(gt3,1) - 14.142135623730951*LamSU*Conj(MuU)*ZP(gt2,1)*ZP(gt3,1) + 10*LamTU*vu*Conj(LamSU)*ZP(gt2,2)*ZP(gt3,1) - 10*LamSU*vu*Conj(LamTU)*ZP(gt2,3)*ZP(gt3,1) - 10*LamSU*vu*Conj(LamTU)*ZP(gt2,1)*ZP(gt3,2) + 10*LamTU*vu*Conj(LamSU)*ZP(gt2,1)*ZP(gt3,3) + ZP(gt2,0)*((7.745966692414834*g1*MDBS + 7.0710678118654755*(2*MuD - LamTD*vT)*Conj(LamSD) + 7.0710678118654755*LamSD*vT*Conj(LamTD) - 7.745966692414834*g1*Conj(MDBS) - 14.142135623730951*LamSD*Conj(MuD))*ZP(gt3,0) + 10*vd*(-(LamSD*Conj(LamTD)*ZP(gt3,2)) + LamTD*Conj(LamSD)*ZP(gt3,3)))) - 5*(ZA(gt1,0)*(vu*Sqr(g2)*ZP(gt2,1)*ZP(gt3,0) + (2*LamTD*vS*Conj(LamSD) + 1.4142135623730951*(2*g2*MDWBT - vT*AbsSqr(LamTD) + 2*LamTD*Conj(MuD) + vT*Sqr(g2)))*ZP(gt2,2)*ZP(gt3,0) + 1.4142135623730951*vT*AbsSqr(LamTD)*ZP(gt2,3)*ZP(gt3,0) + 2.8284271247461903*MuD*Conj(LamTD)*ZP(gt2,3)*ZP(gt3,0) + 2*LamSD*vS*Conj(LamTD)*ZP(gt2,3)*ZP(gt3,0) + 2.8284271247461903*g2*Conj(MDWBT)*ZP(gt2,3)*ZP(gt3,0) - 1.4142135623730951*vT*Sqr(g2)*ZP(gt2,3)*ZP(gt3,0) - vu*Sqr(g2)*ZP(gt2,0)*ZP(gt3,1) + 1.4142135623730951*vT*AbsSqr(LamTD)*ZP(gt2,0)*ZP(gt3,2) - 2.8284271247461903*MuD*Conj(LamTD)*ZP(gt2,0)*ZP(gt3,2) - 2*LamSD*vS*Conj(LamTD)*ZP(gt2,0)*ZP(gt3,2) - 2.8284271247461903*g2*Conj(MDWBT)*ZP(gt2,0)*ZP(gt3,2) - 1.4142135623730951*vT*Sqr(g2)*ZP(gt2,0)*ZP(gt3,2) - 2.8284271247461903*g2*MDWBT*ZP(gt2,0)*ZP(gt3,3) - 1.4142135623730951*vT*AbsSqr(LamTD)*ZP(gt2,0)*ZP(gt3,3) - 2*LamTD*vS*Conj(LamSD)*ZP(gt2,0)*ZP(gt3,3) - 2.8284271247461903*LamTD*Conj(MuD)*ZP(gt2,0)*ZP(gt3,3) + 1.4142135623730951*vT*Sqr(g2)*ZP(gt2,0)*ZP(gt3,3)) + ZA(gt1,1)*(-((vd*Sqr(g2)*ZP(gt2,0) + (2*LamTU*vS*Conj(LamSU) + 1.4142135623730951*(2*g2*MDWBT - vT*AbsSqr(LamTU) + 2*LamTU*Conj(MuU) + vT*Sqr(g2)))*ZP(gt2,2) + ((2.8284271247461903*MuU + 2*LamSU*vS + 1.4142135623730951*LamTU*vT)*Conj(LamTU) + 1.4142135623730951*g2*(-(g2*vT) + 2*Conj(MDWBT)))*ZP(gt2,3))*ZP(gt3,1)) + ZP(gt2,1)*(vd*Sqr(g2)*ZP(gt3,0) + ((2.8284271247461903*MuU + 2*LamSU*vS - 1.4142135623730951*LamTU*vT)*Conj(LamTU) + 1.4142135623730951*g2*(g2*vT + 2*Conj(MDWBT)))*ZP(gt3,2) + (2*LamTU*vS*Conj(LamSU) + 1.4142135623730951*(2*g2*MDWBT + vT*AbsSqr(LamTU) + 2*LamTU*Conj(MuU) - vT*Sqr(g2)))*ZP(gt3,3))) + ZA(gt1,3)*(1.4142135623730951*vd*AbsSqr(LamTD)*ZP(gt2,3)*ZP(gt3,0) - 1.4142135623730951*vd*Sqr(g2)*ZP(gt2,3)*ZP(gt3,0) + 2*g2*MDWBT*ZP(gt2,1)*ZP(gt3,1) + 1.4142135623730951*LamTU*vS*Conj(LamSU)*ZP(gt2,1)*ZP(gt3,1) - 2*MuU*Conj(LamTU)*ZP(gt2,1)*ZP(gt3,1) - 1.4142135623730951*LamSU*vS*Conj(LamTU)*ZP(gt2,1)*ZP(gt3,1) - 2*g2*Conj(MDWBT)*ZP(gt2,1)*ZP(gt3,1) + 2*LamTU*Conj(MuU)*ZP(gt2,1)*ZP(gt3,1) + 1.4142135623730951*vu*AbsSqr(LamTU)*ZP(gt2,3)*ZP(gt3,1) - 1.4142135623730951*vu*Sqr(g2)*ZP(gt2,3)*ZP(gt3,1) - 1.4142135623730951*vu*AbsSqr(LamTU)*ZP(gt2,1)*ZP(gt3,2) + 1.4142135623730951*vu*Sqr(g2)*ZP(gt2,1)*ZP(gt3,2) - 4*vT*Sqr(g2)*ZP(gt2,3)*ZP(gt3,2) - 1.4142135623730951*vu*AbsSqr(LamTU)*ZP(gt2,1)*ZP(gt3,3) + 1.4142135623730951*vu*Sqr(g2)*ZP(gt2,1)*ZP(gt3,3) + ZP(gt2,2)*(1.4142135623730951*vd*(AbsSqr(LamTD) - Sqr(g2))*ZP(gt3,0) + 1.4142135623730951*vu*(AbsSqr(LamTU) - Sqr(g2))*ZP(gt3,1) + 4*vT*Sqr(g2)*ZP(gt3,3)) + ZP(gt2,0)*((-2*g2*MDWBT - 1.4142135623730951*LamTD*vS*Conj(LamSD) + 2*MuD*Conj(LamTD) + 1.4142135623730951*LamSD*vS*Conj(LamTD) + 2*g2*Conj(MDWBT) - 2*LamTD*Conj(MuD))*ZP(gt3,0) + 1.4142135623730951*vd*(-AbsSqr(LamTD) + Sqr(g2))*(ZP(gt3,2) + ZP(gt3,3))))));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Rh>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Mu = MODELPARAMETER(Mu);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = Mu*(Conj(LamTD)*ZHR(gt3,0)*ZP(gt1,3)*ZP(gt2,1) - Conj(LamTU)*ZHR(gt3,1)*ZP(gt1,0)*ZP(gt2,2));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::Su, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Mu = MODELPARAMETER(Mu);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto vT = MODELPARAMETER(vT);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(2*(2*Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1)))*ZD(gt3,j2))*ZP(gt1,0) + 1.4142135623730951*vd*SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gt3,j3))*ZP(gt1,0) + 2*Conj(Mu)*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1)))*ZP(gt1,1) + 1.4142135623730951*vu*SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZD(gt3,j3))*ZP(gt1,1) + 1.4142135623730951*SUM(j3,0,2,Conj(ZU(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yd(j2,j1))*ZD(gt3,3 + j2)))*(vu*ZP(gt1,0) + vd*ZP(gt1,1))) - g2*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZD(gt3,j1))*(1.4142135623730951*g2*vd*ZP(gt1,0) + 1.4142135623730951*g2*vu*ZP(gt1,1) + 4*MDWBT*ZP(gt1,2) + 2*g2*vT*ZP(gt1,2) - 2*g2*vT*ZP(gt1,3) + 4*Conj(MDWBT)*ZP(gt1,3)));

   return {result};
}

ScalarVertex VertexImpl<fields::Sd, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Mu = MODELPARAMETER(Mu);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto g2 = MODELPARAMETER(g2);
   const auto vT = MODELPARAMETER(vT);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(2*(2*Conj(Mu)*SUM(j2,0,2,Conj(ZD(gt1,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt3,3 + j1)))*ZP(gt2,0) + 1.4142135623730951*vd*SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt1,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZU(gt3,j3))*ZP(gt2,0) + 2*Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt1,3 + j1)))*ZU(gt3,j2))*ZP(gt2,1) + 1.4142135623730951*vu*SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt1,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gt3,j3))*ZP(gt2,1) + 1.4142135623730951*SUM(j3,0,2,Conj(ZD(gt1,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yu(j2,j1))*ZU(gt3,3 + j2)))*(vu*ZP(gt2,0) + vd*ZP(gt2,1))) - g2*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZU(gt3,j1))*(1.4142135623730951*g2*vd*ZP(gt2,0) + 1.4142135623730951*g2*vu*ZP(gt2,1) + 2*g2*vT*ZP(gt2,2) + 4*Conj(MDWBT)*ZP(gt2,2) + 4*MDWBT*ZP(gt2,3) - 2*g2*vT*ZP(gt2,3)));

   return {result};
}

ScalarVertex VertexImpl<fields::Sv, typename fields::conj<fields::Se>::type, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto vS = MODELPARAMETER(vS);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto vT = MODELPARAMETER(vT);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto MuD = MODELPARAMETER(MuD);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);

   const std::complex<double> result = 0.5*(1.4142135623730951*vS*Conj(LamSD) - vT*Conj(LamTD) + 2*Conj(MuD))*SUM(j2,0,2,Conj(ZV(gt1,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt2,3 + j1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Se, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto vd = MODELPARAMETER(vd);
   const auto Mu = MODELPARAMETER(Mu);
   const auto g2 = MODELPARAMETER(g2);
   const auto vu = MODELPARAMETER(vu);
   const auto vT = MODELPARAMETER(vT);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.7071067811865475*vd*SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gt1,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZV(gt3,j3))*ZP(gt2,0) + Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt1,3 + j1)))*ZV(gt3,j2))*ZP(gt2,1) - 0.25*g2*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZV(gt3,j1))*(1.4142135623730951*g2*vd*ZP(gt2,0) + 1.4142135623730951*g2*vu*ZP(gt2,1) + 2*g2*vT*ZP(gt2,2) + 4*Conj(MDWBT)*ZP(gt2,2) + 4*MDWBT*ZP(gt2,3) - 2*g2*vT*ZP(gt2,3));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::hh, fields::Rh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Mu = MODELPARAMETER(Mu);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = std::complex<double>(0,0.25)*Conj(Mu)*(2*LamSD*ZA(gt1,1)*ZH(gt2,2)*ZHR(gt3,0) + 1.4142135623730951*LamTD*ZA(gt1,1)*ZH(gt2,3)*ZHR(gt3,0) - 2*LamSU*ZA(gt1,0)*ZH(gt2,2)*ZHR(gt3,1) + 1.4142135623730951*LamTU*ZA(gt1,0)*ZH(gt2,3)*ZHR(gt3,1) + ZA(gt1,2)*(-2*LamSD*ZH(gt2,1)*ZHR(gt3,0) + 2*LamSU*ZH(gt2,0)*ZHR(gt3,1)) - 1.4142135623730951*ZA(gt1,3)*(LamTD*ZH(gt2,1)*ZHR(gt3,0) + LamTU*ZH(gt2,0)*ZHR(gt3,1)));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::Ah, fields::hh, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = std::complex<double>(0,-0.5)*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(ZA(gt1,0)*ZH(gt2,0) - ZA(gt1,1)*ZH(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*((g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*ZP(gt1,0)*ZP(gt2,0) + (g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*ZP(gt1,1)*ZP(gt2,1) + 2*g2*Cos(ThetaW)*(ZP(gt1,2)*ZP(gt2,2) + ZP(gt1,3)*ZP(gt2,3)));

   return {result, minuend_index, subtrahend_index};
}

ScalarVertex VertexImpl<fields::Hpm, fields::Rh, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Mu = MODELPARAMETER(Mu);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = Conj(Mu)*(-(LamTU*ZHR(gt2,1)*ZP(gt1,2)*ZP(gt3,0)) + LamTD*ZHR(gt2,0)*ZP(gt1,1)*ZP(gt3,3));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::Rh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Mu = MODELPARAMETER(Mu);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.25*Conj(Mu)*(2*LamSD*ZH(gt1,1)*ZH(gt2,2)*ZHR(gt3,0) + 1.4142135623730951*LamTD*ZH(gt1,1)*ZH(gt2,3)*ZHR(gt3,0) - 2*LamSU*ZH(gt1,0)*ZH(gt2,2)*ZHR(gt3,1) + 1.4142135623730951*LamTU*ZH(gt1,0)*ZH(gt2,3)*ZHR(gt3,1) + 2*ZH(gt1,2)*(LamSD*ZH(gt2,1)*ZHR(gt3,0) - LamSU*ZH(gt2,0)*ZHR(gt3,1)) + 1.4142135623730951*ZH(gt1,3)*(LamTD*ZH(gt2,1)*ZHR(gt3,0) + LamTU*ZH(gt2,0)*ZHR(gt3,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Hpm, fields::SRdp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Mu = MODELPARAMETER(Mu);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.5*Conj(Mu)*(-1.4142135623730951*LamSD*ZH(gt1,2)*ZP(gt2,1) + LamTD*ZH(gt1,3)*ZP(gt2,1) + 1.4142135623730951*LamTD*ZH(gt1,1)*ZP(gt2,2));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::SRum, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const auto Mu = MODELPARAMETER(Mu);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,0.5)*Conj(Mu)*(1.4142135623730951*LamSU*ZA(gt1,2)*ZP(gt3,0) + LamTU*ZA(gt1,3)*ZP(gt3,0) + 1.4142135623730951*LamTU*ZA(gt1,0)*ZP(gt3,3));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::SRum, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const auto Mu = MODELPARAMETER(Mu);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.5*Conj(Mu)*(1.4142135623730951*LamSU*ZH(gt1,2)*ZP(gt3,0) + LamTU*ZH(gt1,3)*ZP(gt3,0) - 1.4142135623730951*LamTU*ZH(gt1,0)*ZP(gt3,3));

   return {result};
}

ScalarVertex VertexImpl<fields::SRum, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Rh>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt3 = indices[1];
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto vd = MODELPARAMETER(vd);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto vu = MODELPARAMETER(vu);
   const auto g2 = MODELPARAMETER(g2);
   const auto MuU = MODELPARAMETER(MuU);
   const auto vS = MODELPARAMETER(vS);
   const auto vT = MODELPARAMETER(vT);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-2.8284271247461903*LamSU*vd*Conj(LamSD)*ZHR(gt3,0)*ZP(gt2,1) - 1.4142135623730951*LamTU*Conj(LamTD)*ZHR(gt3,0)*(2*vu*ZP(gt2,0) + vd*ZP(gt2,1)) - ZHR(gt3,1)*(1.4142135623730951*vd*Sqr(g2)*ZP(gt2,0) + 1.4142135623730951*vu*(-2*AbsSqr(LamSU) + AbsSqr(LamTU) + Sqr(g2))*ZP(gt2,1) + 2*(-((2*MuU + 1.4142135623730951*LamSU*vS + LamTU*vT)*Conj(LamTU)) + g2*(g2*vT + 2*Conj(MDWBT)))*ZP(gt2,2) - 2*(-2*g2*MDWBT - vT*AbsSqr(LamTU) + 1.4142135623730951*LamTU*vS*Conj(LamSU) + 2*LamTU*Conj(MuU) + vT*Sqr(g2))*ZP(gt2,3)));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::hh, fields::Hpm, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.5*g2*(ZH(gt1,0)*ZP(gt2,0) - ZH(gt1,1)*ZP(gt2,1) + 1.4142135623730951*ZH(gt1,3)*(ZP(gt2,2) + ZP(gt2,3)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::Ah, typename fields::conj<fields::Hpm>::type, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,-0.5)*g2*(ZA(gt1,0)*ZP(gt2,0) + ZA(gt1,1)*ZP(gt2,1) + 1.4142135623730951*ZA(gt1,3)*(ZP(gt2,2) - ZP(gt2,3)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::hh, typename fields::conj<fields::Hpm>::type, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -0.5*g2*(ZH(gt1,0)*ZP(gt2,0) - ZH(gt1,1)*ZP(gt2,1) + 1.4142135623730951*ZH(gt1,3)*(ZP(gt2,2) + ZP(gt2,3)));

   return {result, minuend_index, subtrahend_index};
}

ScalarVertex VertexImpl<fields::hh, fields::Se, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto vd = MODELPARAMETER(vd);
   const auto Mu = MODELPARAMETER(Mu);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.05*(2*(5*(-2*vd*SUM(j3,0,2,Conj(ZE(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Ye(j2,j1))*ZE(gt3,3 + j2)))*ZH(gt1,0) - 2*vd*SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZE(gt3,j3))*ZH(gt1,0) + 1.4142135623730951*(Conj(Mu)*SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1))) + Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1)))*ZE(gt3,j2)))*ZH(gt1,1)) + g1*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*(3*g1*vd*ZH(gt1,0) - 3*g1*vu*ZH(gt1,1) - 7.745966692414834*(MDBS + Conj(MDBS))*ZH(gt1,2))) + SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*((-3*vd*Sqr(g1) + 5*vd*Sqr(g2))*ZH(gt1,0) + vu*(3*Sqr(g1) - 5*Sqr(g2))*ZH(gt1,1) + 7.745966692414834*g1*(MDBS + Conj(MDBS))*ZH(gt1,2) + 10*g2*(MDWBT + Conj(MDWBT))*ZH(gt1,3)));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::Rh, typename fields::conj<fields::SRdp>::type, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.7071067811865475*g2*ZHR(gt1,0);

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::SRdp, typename fields::conj<fields::SRdp>::type, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.1*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW));

   return {result, minuend_index, subtrahend_index};
}

ScalarVertex VertexImpl<fields::hh, fields::Rh, typename fields::conj<fields::Rh>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MuD = MODELPARAMETER(MuD);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto vS = MODELPARAMETER(vS);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto vT = MODELPARAMETER(vT);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto MuU = MODELPARAMETER(MuU);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto vu = MODELPARAMETER(vu);
   const auto vd = MODELPARAMETER(vd);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = 0.25*(-1.5491933384829668*g1*MDBS*ZH(gt1,2)*ZHR(gt2,0)*ZHR(gt3,0) - 4*vS*AbsSqr(LamSD)*ZH(gt1,2)*ZHR(gt2,0)*ZHR(gt3,0) - 2.8284271247461903*MuD*Conj(LamSD)*ZH(gt1,2)*ZHR(gt2,0)*ZHR(gt3,0) - 1.4142135623730951*LamTD*vT*Conj(LamSD)*ZH(gt1,2)*ZHR(gt2,0)*ZHR(gt3,0) - 1.4142135623730951*LamSD*vT*Conj(LamTD)*ZH(gt1,2)*ZHR(gt2,0)*ZHR(gt3,0) - 1.5491933384829668*g1*Conj(MDBS)*ZH(gt1,2)*ZHR(gt2,0)*ZHR(gt3,0) - 2.8284271247461903*LamSD*Conj(MuD)*ZH(gt1,2)*ZHR(gt2,0)*ZHR(gt3,0) + 2*g2*MDWBT*ZH(gt1,3)*ZHR(gt2,0)*ZHR(gt3,0) - 2*vT*AbsSqr(LamTD)*ZH(gt1,3)*ZHR(gt2,0)*ZHR(gt3,0) - 1.4142135623730951*LamTD*vS*Conj(LamSD)*ZH(gt1,3)*ZHR(gt2,0)*ZHR(gt3,0) - 2*MuD*Conj(LamTD)*ZH(gt1,3)*ZHR(gt2,0)*ZHR(gt3,0) - 1.4142135623730951*LamSD*vS*Conj(LamTD)*ZH(gt1,3)*ZHR(gt2,0)*ZHR(gt3,0) + 2*g2*Conj(MDWBT)*ZH(gt1,3)*ZHR(gt2,0)*ZHR(gt3,0) - 2*LamTD*Conj(MuD)*ZH(gt1,3)*ZHR(gt2,0)*ZHR(gt3,0) + 1.5491933384829668*g1*MDBS*ZH(gt1,2)*ZHR(gt2,1)*ZHR(gt3,1) - 4*vS*AbsSqr(LamSU)*ZH(gt1,2)*ZHR(gt2,1)*ZHR(gt3,1) - 2.8284271247461903*MuU*Conj(LamSU)*ZH(gt1,2)*ZHR(gt2,1)*ZHR(gt3,1) + 1.4142135623730951*LamTU*vT*Conj(LamSU)*ZH(gt1,2)*ZHR(gt2,1)*ZHR(gt3,1) + 1.4142135623730951*LamSU*vT*Conj(LamTU)*ZH(gt1,2)*ZHR(gt2,1)*ZHR(gt3,1) + 1.5491933384829668*g1*Conj(MDBS)*ZH(gt1,2)*ZHR(gt2,1)*ZHR(gt3,1) - 2.8284271247461903*LamSU*Conj(MuU)*ZH(gt1,2)*ZHR(gt2,1)*ZHR(gt3,1) - 2*g2*MDWBT*ZH(gt1,3)*ZHR(gt2,1)*ZHR(gt3,1) - 2*vT*AbsSqr(LamTU)*ZH(gt1,3)*ZHR(gt2,1)*ZHR(gt3,1) + 1.4142135623730951*LamTU*vS*Conj(LamSU)*ZH(gt1,3)*ZHR(gt2,1)*ZHR(gt3,1) + 2*MuU*Conj(LamTU)*ZH(gt1,3)*ZHR(gt2,1)*ZHR(gt3,1) + 1.4142135623730951*LamSU*vS*Conj(LamTU)*ZH(gt1,3)*ZHR(gt2,1)*ZHR(gt3,1) - 2*g2*Conj(MDWBT)*ZH(gt1,3)*ZHR(gt2,1)*ZHR(gt3,1) + 2*LamTU*Conj(MuU)*ZH(gt1,3)*ZHR(gt2,1)*ZHR(gt3,1) + ZH(gt1,0)*(ZHR(gt2,0)*(vd*(-4*AbsSqr(LamSD) - 2*AbsSqr(LamTD) + 0.6*Sqr(g1) + Sqr(g2))*ZHR(gt3,0) + vu*(2*LamSD*Conj(LamSU) - LamTD*Conj(LamTU))*ZHR(gt3,1)) + ZHR(gt2,1)*(2*LamSU*vu*Conj(LamSD)*ZHR(gt3,0) - LamTU*vu*Conj(LamTD)*ZHR(gt3,0) - vd*(0.6*Sqr(g1) + Sqr(g2))*ZHR(gt3,1))) - ZH(gt1,1)*(ZHR(gt2,0)*(vu*(0.6*Sqr(g1) + Sqr(g2))*ZHR(gt3,0) + vd*(-2*LamSD*Conj(LamSU) + LamTD*Conj(LamTU))*ZHR(gt3,1)) + ZHR(gt2,1)*(-2*LamSU*vd*Conj(LamSD)*ZHR(gt3,0) + LamTU*vd*Conj(LamTD)*ZHR(gt3,0) - vu*(-4*AbsSqr(LamSU) - 2*AbsSqr(LamTU) + 0.6*Sqr(g1) + Sqr(g2))*ZHR(gt3,1))));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::SRdp, typename fields::conj<fields::SRdp>::type, fields::VP>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.1*(-3.872983346207417*g1*Cos(ThetaW) - 5*g2*Sin(ThetaW));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*((0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*ZP(gt1,0)*ZP(gt2,0) + (0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*ZP(gt1,1)*ZP(gt2,1) + 2*g2*Sin(ThetaW)*(ZP(gt1,2)*ZP(gt2,2) + ZP(gt1,3)*ZP(gt2,3)));

   return {result, minuend_index, subtrahend_index};
}

ScalarVertex VertexImpl<fields::hh, fields::Sv, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = -0.05*KroneckerDelta(gt2,gt3)*(vd*(3*Sqr(g1) + 5*Sqr(g2))*ZH(gt1,0) - vu*(3*Sqr(g1) + 5*Sqr(g2))*ZH(gt1,1) - 7.745966692414834*g1*(MDBS + Conj(MDBS))*ZH(gt1,2) + 10*g2*(MDWBT + Conj(MDWBT))*ZH(gt1,3));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Hpm, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.05*(ZA(gt1,0)*(-(ZA(gt2,0)*((3*Sqr(g1) + 5*Sqr(g2))*ZP(gt3,0)*ZP(gt4,0) + (-3*Sqr(g1) + 5*Sqr(g2))*ZP(gt3,1)*ZP(gt4,1) - 10*(-2*AbsSqr(LamTD) + Sqr(g2))*ZP(gt3,2)*ZP(gt4,2) + 10*Sqr(g2)*ZP(gt3,3)*ZP(gt4,3))) + 5*(1.4142135623730951*AbsSqr(LamTD)*ZA(gt2,3)*ZP(gt3,2)*ZP(gt4,0) - 1.4142135623730951*Sqr(g2)*ZA(gt2,3)*ZP(gt3,2)*ZP(gt4,0) + 2*LamSD*Conj(LamTD)*ZA(gt2,2)*ZP(gt3,3)*ZP(gt4,0) + 1.4142135623730951*AbsSqr(LamTD)*ZA(gt2,3)*ZP(gt3,3)*ZP(gt4,0) - 1.4142135623730951*Sqr(g2)*ZA(gt2,3)*ZP(gt3,3)*ZP(gt4,0) + Sqr(g2)*ZA(gt2,1)*(ZP(gt3,1)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,1)) - 2*LamSD*Conj(LamTD)*ZA(gt2,2)*ZP(gt3,0)*ZP(gt4,2) + 1.4142135623730951*AbsSqr(LamTD)*ZA(gt2,3)*ZP(gt3,0)*ZP(gt4,2) - 1.4142135623730951*Sqr(g2)*ZA(gt2,3)*ZP(gt3,0)*ZP(gt4,2) + 1.4142135623730951*AbsSqr(LamTD)*ZA(gt2,3)*ZP(gt3,0)*ZP(gt4,3) - 1.4142135623730951*Sqr(g2)*ZA(gt2,3)*ZP(gt3,0)*ZP(gt4,3) + 2*LamTD*Conj(LamSD)*ZA(gt2,2)*(-(ZP(gt3,2)*ZP(gt4,0)) + ZP(gt3,0)*ZP(gt4,3)))) + ZA(gt1,1)*(ZA(gt2,1)*((3*Sqr(g1) - 5*Sqr(g2))*ZP(gt3,0)*ZP(gt4,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZP(gt3,1)*ZP(gt4,1) - 10*Sqr(g2)*ZP(gt3,2)*ZP(gt4,2) + 10*(-2*AbsSqr(LamTU) + Sqr(g2))*ZP(gt3,3)*ZP(gt4,3)) + 5*(-1.4142135623730951*AbsSqr(LamTU)*ZA(gt2,3)*ZP(gt3,2)*ZP(gt4,1) + 1.4142135623730951*Sqr(g2)*ZA(gt2,3)*ZP(gt3,2)*ZP(gt4,1) - 2*LamSU*Conj(LamTU)*ZA(gt2,2)*ZP(gt3,3)*ZP(gt4,1) - 1.4142135623730951*AbsSqr(LamTU)*ZA(gt2,3)*ZP(gt3,3)*ZP(gt4,1) + 1.4142135623730951*Sqr(g2)*ZA(gt2,3)*ZP(gt3,3)*ZP(gt4,1) + Sqr(g2)*ZA(gt2,0)*(ZP(gt3,1)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,1)) + 2*LamSU*Conj(LamTU)*ZA(gt2,2)*ZP(gt3,1)*ZP(gt4,2) - 1.4142135623730951*AbsSqr(LamTU)*ZA(gt2,3)*ZP(gt3,1)*ZP(gt4,2) + 1.4142135623730951*Sqr(g2)*ZA(gt2,3)*ZP(gt3,1)*ZP(gt4,2) - 1.4142135623730951*AbsSqr(LamTU)*ZA(gt2,3)*ZP(gt3,1)*ZP(gt4,3) + 1.4142135623730951*Sqr(g2)*ZA(gt2,3)*ZP(gt3,1)*ZP(gt4,3) + 2*LamTU*Conj(LamSU)*ZA(gt2,2)*(ZP(gt3,2)*ZP(gt4,1) - ZP(gt3,1)*ZP(gt4,3)))) - 5*(1.4142135623730951*Sqr(g2)*ZA(gt1,3)*ZA(gt2,0)*ZP(gt3,2)*ZP(gt4,0) + 1.4142135623730951*Sqr(g2)*ZA(gt1,3)*ZA(gt2,0)*ZP(gt3,3)*ZP(gt4,0) + 4*AbsSqr(LamSU)*ZA(gt1,2)*ZA(gt2,2)*ZP(gt3,1)*ZP(gt4,1) + 1.4142135623730951*LamTU*Conj(LamSU)*ZA(gt1,3)*ZA(gt2,2)*ZP(gt3,1)*ZP(gt4,1) + 1.4142135623730951*LamSU*Conj(LamTU)*ZA(gt1,3)*ZA(gt2,2)*ZP(gt3,1)*ZP(gt4,1) + 1.4142135623730951*LamTU*Conj(LamSU)*ZA(gt1,2)*ZA(gt2,3)*ZP(gt3,1)*ZP(gt4,1) + 1.4142135623730951*LamSU*Conj(LamTU)*ZA(gt1,2)*ZA(gt2,3)*ZP(gt3,1)*ZP(gt4,1) + 2*AbsSqr(LamTU)*ZA(gt1,3)*ZA(gt2,3)*ZP(gt3,1)*ZP(gt4,1) - 2*LamTU*Conj(LamSU)*ZA(gt1,2)*ZA(gt2,1)*ZP(gt3,2)*ZP(gt4,1) + 1.4142135623730951*AbsSqr(LamTU)*ZA(gt1,3)*ZA(gt2,1)*ZP(gt3,2)*ZP(gt4,1) - 1.4142135623730951*Sqr(g2)*ZA(gt1,3)*ZA(gt2,1)*ZP(gt3,2)*ZP(gt4,1) + 2*LamSU*Conj(LamTU)*ZA(gt1,2)*ZA(gt2,1)*ZP(gt3,3)*ZP(gt4,1) + 1.4142135623730951*AbsSqr(LamTU)*ZA(gt1,3)*ZA(gt2,1)*ZP(gt3,3)*ZP(gt4,1) - 1.4142135623730951*Sqr(g2)*ZA(gt1,3)*ZA(gt2,1)*ZP(gt3,3)*ZP(gt4,1) + 1.4142135623730951*Sqr(g2)*ZA(gt1,3)*ZA(gt2,0)*ZP(gt3,0)*ZP(gt4,2) - 2*LamSU*Conj(LamTU)*ZA(gt1,2)*ZA(gt2,1)*ZP(gt3,1)*ZP(gt4,2) + 1.4142135623730951*AbsSqr(LamTU)*ZA(gt1,3)*ZA(gt2,1)*ZP(gt3,1)*ZP(gt4,2) - 1.4142135623730951*Sqr(g2)*ZA(gt1,3)*ZA(gt2,1)*ZP(gt3,1)*ZP(gt4,2) + 4*Sqr(g2)*ZA(gt1,3)*ZA(gt2,3)*ZP(gt3,2)*ZP(gt4,2) + 4*Sqr(g2)*ZA(gt1,3)*ZA(gt2,3)*ZP(gt3,3)*ZP(gt4,2) + 1.4142135623730951*Sqr(g2)*ZA(gt1,3)*ZA(gt2,0)*ZP(gt3,0)*ZP(gt4,3) + 2*LamTU*Conj(LamSU)*ZA(gt1,2)*ZA(gt2,1)*ZP(gt3,1)*ZP(gt4,3) + 1.4142135623730951*AbsSqr(LamTU)*ZA(gt1,3)*ZA(gt2,1)*ZP(gt3,1)*ZP(gt4,3) - 1.4142135623730951*Sqr(g2)*ZA(gt1,3)*ZA(gt2,1)*ZP(gt3,1)*ZP(gt4,3) + 4*Sqr(g2)*ZA(gt1,3)*ZA(gt2,3)*ZP(gt3,2)*ZP(gt4,3) + 4*Sqr(g2)*ZA(gt1,3)*ZA(gt2,3)*ZP(gt3,3)*ZP(gt4,3) - Conj(LamTD)*(LamSD*ZA(gt1,2)*(1.4142135623730951*ZA(gt2,3)*ZP(gt3,0)*ZP(gt4,0) + 2*ZA(gt2,0)*(ZP(gt3,3)*ZP(gt4,0) - ZP(gt3,0)*ZP(gt4,2))) + ZA(gt1,3)*(1.4142135623730951*LamSD*ZA(gt2,2)*ZP(gt3,0)*ZP(gt4,0) - 2*LamTD*ZA(gt2,3)*ZP(gt3,0)*ZP(gt4,0) + 1.4142135623730951*LamTD*ZA(gt2,0)*(ZP(gt3,2)*ZP(gt4,0) + ZP(gt3,3)*ZP(gt4,0) + ZP(gt3,0)*(ZP(gt4,2) + ZP(gt4,3))))) + Conj(LamSD)*(-1.4142135623730951*LamTD*ZA(gt1,3)*ZA(gt2,2)*ZP(gt3,0)*ZP(gt4,0) + ZA(gt1,2)*(4*LamSD*ZA(gt2,2)*ZP(gt3,0)*ZP(gt4,0) - LamTD*(1.4142135623730951*ZA(gt2,3)*ZP(gt3,0)*ZP(gt4,0) + 2*ZA(gt2,0)*(-(ZP(gt3,2)*ZP(gt4,0)) + ZP(gt3,0)*ZP(gt4,3)))))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::Hpm, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.05*(-(ZH(gt1,0)*(ZH(gt2,0)*((3*Sqr(g1) + 5*Sqr(g2))*ZP(gt3,0)*ZP(gt4,0) + (-3*Sqr(g1) + 5*Sqr(g2))*ZP(gt3,1)*ZP(gt4,1) - 10*(-2*AbsSqr(LamTD) + Sqr(g2))*ZP(gt3,2)*ZP(gt4,2) + 10*Sqr(g2)*ZP(gt3,3)*ZP(gt4,3)) + 5*(-1.4142135623730951*AbsSqr(LamTD)*ZH(gt2,3)*ZP(gt3,2)*ZP(gt4,0) + 1.4142135623730951*Sqr(g2)*ZH(gt2,3)*ZP(gt3,2)*ZP(gt4,0) + 2*LamSD*Conj(LamTD)*ZH(gt2,2)*ZP(gt3,3)*ZP(gt4,0) + 1.4142135623730951*AbsSqr(LamTD)*ZH(gt2,3)*ZP(gt3,3)*ZP(gt4,0) - 1.4142135623730951*Sqr(g2)*ZH(gt2,3)*ZP(gt3,3)*ZP(gt4,0) + Sqr(g2)*ZH(gt2,1)*(ZP(gt3,1)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,1)) + 2*LamSD*Conj(LamTD)*ZH(gt2,2)*ZP(gt3,0)*ZP(gt4,2) - 1.4142135623730951*AbsSqr(LamTD)*ZH(gt2,3)*ZP(gt3,0)*ZP(gt4,2) + 1.4142135623730951*Sqr(g2)*ZH(gt2,3)*ZP(gt3,0)*ZP(gt4,2) + 1.4142135623730951*AbsSqr(LamTD)*ZH(gt2,3)*ZP(gt3,0)*ZP(gt4,3) - 1.4142135623730951*Sqr(g2)*ZH(gt2,3)*ZP(gt3,0)*ZP(gt4,3) + 2*LamTD*Conj(LamSD)*ZH(gt2,2)*(ZP(gt3,2)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,3))))) + ZH(gt1,1)*(ZH(gt2,1)*((3*Sqr(g1) - 5*Sqr(g2))*ZP(gt3,0)*ZP(gt4,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZP(gt3,1)*ZP(gt4,1) - 10*Sqr(g2)*ZP(gt3,2)*ZP(gt4,2) + 10*(-2*AbsSqr(LamTU) + Sqr(g2))*ZP(gt3,3)*ZP(gt4,3)) - 5*(-1.4142135623730951*AbsSqr(LamTU)*ZH(gt2,3)*ZP(gt3,2)*ZP(gt4,1) + 1.4142135623730951*Sqr(g2)*ZH(gt2,3)*ZP(gt3,2)*ZP(gt4,1) + 2*LamSU*Conj(LamTU)*ZH(gt2,2)*ZP(gt3,3)*ZP(gt4,1) + 1.4142135623730951*AbsSqr(LamTU)*ZH(gt2,3)*ZP(gt3,3)*ZP(gt4,1) - 1.4142135623730951*Sqr(g2)*ZH(gt2,3)*ZP(gt3,3)*ZP(gt4,1) + Sqr(g2)*ZH(gt2,0)*(ZP(gt3,1)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,1)) + 2*LamSU*Conj(LamTU)*ZH(gt2,2)*ZP(gt3,1)*ZP(gt4,2) - 1.4142135623730951*AbsSqr(LamTU)*ZH(gt2,3)*ZP(gt3,1)*ZP(gt4,2) + 1.4142135623730951*Sqr(g2)*ZH(gt2,3)*ZP(gt3,1)*ZP(gt4,2) + 1.4142135623730951*AbsSqr(LamTU)*ZH(gt2,3)*ZP(gt3,1)*ZP(gt4,3) - 1.4142135623730951*Sqr(g2)*ZH(gt2,3)*ZP(gt3,1)*ZP(gt4,3) + 2*LamTU*Conj(LamSU)*ZH(gt2,2)*(ZP(gt3,2)*ZP(gt4,1) + ZP(gt3,1)*ZP(gt4,3)))) - 5*(1.4142135623730951*Sqr(g2)*ZH(gt1,3)*ZH(gt2,0)*ZP(gt3,2)*ZP(gt4,0) - 1.4142135623730951*Sqr(g2)*ZH(gt1,3)*ZH(gt2,0)*ZP(gt3,3)*ZP(gt4,0) + 4*AbsSqr(LamSU)*ZH(gt1,2)*ZH(gt2,2)*ZP(gt3,1)*ZP(gt4,1) + 1.4142135623730951*LamTU*Conj(LamSU)*ZH(gt1,3)*ZH(gt2,2)*ZP(gt3,1)*ZP(gt4,1) + 1.4142135623730951*LamSU*Conj(LamTU)*ZH(gt1,3)*ZH(gt2,2)*ZP(gt3,1)*ZP(gt4,1) + 1.4142135623730951*LamTU*Conj(LamSU)*ZH(gt1,2)*ZH(gt2,3)*ZP(gt3,1)*ZP(gt4,1) + 1.4142135623730951*LamSU*Conj(LamTU)*ZH(gt1,2)*ZH(gt2,3)*ZP(gt3,1)*ZP(gt4,1) + 2*AbsSqr(LamTU)*ZH(gt1,3)*ZH(gt2,3)*ZP(gt3,1)*ZP(gt4,1) + 2*LamTU*Conj(LamSU)*ZH(gt1,2)*ZH(gt2,1)*ZP(gt3,2)*ZP(gt4,1) - 1.4142135623730951*AbsSqr(LamTU)*ZH(gt1,3)*ZH(gt2,1)*ZP(gt3,2)*ZP(gt4,1) + 1.4142135623730951*Sqr(g2)*ZH(gt1,3)*ZH(gt2,1)*ZP(gt3,2)*ZP(gt4,1) + 2*LamSU*Conj(LamTU)*ZH(gt1,2)*ZH(gt2,1)*ZP(gt3,3)*ZP(gt4,1) + 1.4142135623730951*AbsSqr(LamTU)*ZH(gt1,3)*ZH(gt2,1)*ZP(gt3,3)*ZP(gt4,1) - 1.4142135623730951*Sqr(g2)*ZH(gt1,3)*ZH(gt2,1)*ZP(gt3,3)*ZP(gt4,1) + 1.4142135623730951*Sqr(g2)*ZH(gt1,3)*ZH(gt2,0)*ZP(gt3,0)*ZP(gt4,2) + 2*LamSU*Conj(LamTU)*ZH(gt1,2)*ZH(gt2,1)*ZP(gt3,1)*ZP(gt4,2) - 1.4142135623730951*AbsSqr(LamTU)*ZH(gt1,3)*ZH(gt2,1)*ZP(gt3,1)*ZP(gt4,2) + 1.4142135623730951*Sqr(g2)*ZH(gt1,3)*ZH(gt2,1)*ZP(gt3,1)*ZP(gt4,2) + 4*Sqr(g2)*ZH(gt1,3)*ZH(gt2,3)*ZP(gt3,2)*ZP(gt4,2) - 4*Sqr(g2)*ZH(gt1,3)*ZH(gt2,3)*ZP(gt3,3)*ZP(gt4,2) - 1.4142135623730951*Sqr(g2)*ZH(gt1,3)*ZH(gt2,0)*ZP(gt3,0)*ZP(gt4,3) + 2*LamTU*Conj(LamSU)*ZH(gt1,2)*ZH(gt2,1)*ZP(gt3,1)*ZP(gt4,3) + 1.4142135623730951*AbsSqr(LamTU)*ZH(gt1,3)*ZH(gt2,1)*ZP(gt3,1)*ZP(gt4,3) - 1.4142135623730951*Sqr(g2)*ZH(gt1,3)*ZH(gt2,1)*ZP(gt3,1)*ZP(gt4,3) - 4*Sqr(g2)*ZH(gt1,3)*ZH(gt2,3)*ZP(gt3,2)*ZP(gt4,3) + 4*Sqr(g2)*ZH(gt1,3)*ZH(gt2,3)*ZP(gt3,3)*ZP(gt4,3) + Conj(LamSD)*(-1.4142135623730951*LamTD*ZH(gt1,3)*ZH(gt2,2)*ZP(gt3,0)*ZP(gt4,0) + ZH(gt1,2)*(4*LamSD*ZH(gt2,2)*ZP(gt3,0)*ZP(gt4,0) - 1.4142135623730951*LamTD*ZH(gt2,3)*ZP(gt3,0)*ZP(gt4,0) + 2*LamTD*ZH(gt2,0)*(ZP(gt3,2)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,3)))) + Conj(LamTD)*(LamSD*ZH(gt1,2)*(-1.4142135623730951*ZH(gt2,3)*ZP(gt3,0)*ZP(gt4,0) + 2*ZH(gt2,0)*(ZP(gt3,3)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,2))) + ZH(gt1,3)*(-1.4142135623730951*LamSD*ZH(gt2,2)*ZP(gt3,0)*ZP(gt4,0) + 2*LamTD*ZH(gt2,3)*ZP(gt3,0)*ZP(gt4,0) + 1.4142135623730951*LamTD*ZH(gt2,0)*(-(ZP(gt3,2)*ZP(gt4,0)) + ZP(gt3,3)*ZP(gt4,0) + ZP(gt3,0)*(-ZP(gt4,2) + ZP(gt4,3)))))));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::Hpm, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(ZP(gt1,0)*(-2*(0.6*Sqr(g1) + Sqr(g2))*ZP(gt2,0)*ZP(gt3,0)*ZP(gt4,0) + (0.6*Sqr(g1) + Sqr(g2))*ZP(gt2,1)*(ZP(gt3,1)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,1)) - 2*Sqr(g2)*ZP(gt2,2)*(ZP(gt3,2)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,2)) + 2*(-2*AbsSqr(LamTD) + Sqr(g2))*ZP(gt2,3)*(ZP(gt3,3)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,3))) + ZP(gt1,1)*(-2*(0.6*Sqr(g1) + Sqr(g2))*ZP(gt2,1)*ZP(gt3,1)*ZP(gt4,1) + (0.6*Sqr(g1) + Sqr(g2))*ZP(gt2,0)*(ZP(gt3,1)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,1)) + 2*(-2*AbsSqr(LamTU) + Sqr(g2))*ZP(gt2,2)*(ZP(gt3,2)*ZP(gt4,1) + ZP(gt3,1)*ZP(gt4,2)) - 2*Sqr(g2)*ZP(gt2,3)*(ZP(gt3,3)*ZP(gt4,1) + ZP(gt3,1)*ZP(gt4,3))) - 2*(ZP(gt1,3)*(-((-2*AbsSqr(LamTD) + Sqr(g2))*ZP(gt2,0)*(ZP(gt3,3)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,3))) + Sqr(g2)*(4*ZP(gt2,3)*ZP(gt3,3)*ZP(gt4,3) + ZP(gt2,1)*(ZP(gt3,3)*ZP(gt4,1) + ZP(gt3,1)*ZP(gt4,3)) - 2*ZP(gt2,2)*(ZP(gt3,3)*ZP(gt4,2) + ZP(gt3,2)*ZP(gt4,3)))) + ZP(gt1,2)*(Sqr(g2)*ZP(gt2,0)*(ZP(gt3,2)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,2)) - (-2*AbsSqr(LamTU) + Sqr(g2))*ZP(gt2,1)*(ZP(gt3,2)*ZP(gt4,1) + ZP(gt3,1)*ZP(gt4,2)) + 2*Sqr(g2)*(2*ZP(gt2,2)*ZP(gt3,2)*ZP(gt4,2) - ZP(gt2,3)*(ZP(gt3,3)*ZP(gt4,2) + ZP(gt3,2)*ZP(gt4,3))))));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::Rh, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Rh>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.05*(-(ZHR(gt2,1)*ZHR(gt4,1)*((3*Sqr(g1) - 5*Sqr(g2))*ZP(gt1,0)*ZP(gt3,0) + (20*AbsSqr(LamTU) - 3*Sqr(g1) + 5*Sqr(g2))*ZP(gt1,1)*ZP(gt3,1) - 10*(-2*AbsSqr(LamTU) + Sqr(g2))*ZP(gt1,2)*ZP(gt3,2) + 10*Sqr(g2)*ZP(gt1,3)*ZP(gt3,3))) + ZHR(gt2,0)*ZHR(gt4,0)*((-20*AbsSqr(LamTD) + 3*Sqr(g1) - 5*Sqr(g2))*ZP(gt1,0)*ZP(gt3,0) + (-3*Sqr(g1) + 5*Sqr(g2))*ZP(gt1,1)*ZP(gt3,1) - 10*Sqr(g2)*ZP(gt1,2)*ZP(gt3,2) + 10*(-2*AbsSqr(LamTD) + Sqr(g2))*ZP(gt1,3)*ZP(gt3,3)));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::Sd, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.05*(2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt4,3 + j1))*(ZP(gt1,0)*ZP(gt3,0) - ZP(gt1,1)*ZP(gt3,1)) - 20*(SUM(j3,0,2,Conj(ZD(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,j1))*ZD(gt4,3 + j2)))*ZP(gt1,0)*ZP(gt3,0) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZD(gt4,j3))*ZP(gt1,1)*ZP(gt3,1)) + SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt4,j1))*((Sqr(g1) - 5*Sqr(g2))*ZP(gt1,0)*ZP(gt3,0) - (Sqr(g1) - 5*Sqr(g2))*ZP(gt1,1)*ZP(gt3,1) + 10*Sqr(g2)*(-(ZP(gt1,2)*ZP(gt3,2)) + ZP(gt1,3)*ZP(gt3,3))));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::Se, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.05*(-20*SUM(j3,0,2,Conj(ZE(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Ye(j2,j1))*ZE(gt4,3 + j2)))*ZP(gt1,0)*ZP(gt3,0) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt4,3 + j1))*(ZP(gt1,0)*ZP(gt3,0) - ZP(gt1,1)*ZP(gt3,1)) + SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt4,j1))*(-((3*Sqr(g1) + 5*Sqr(g2))*ZP(gt1,0)*ZP(gt3,0)) + (3*Sqr(g1) + 5*Sqr(g2))*ZP(gt1,1)*ZP(gt3,1) + 10*Sqr(g2)*(-(ZP(gt1,2)*ZP(gt3,2)) + ZP(gt1,3)*ZP(gt3,3))));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::SRdp, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*((-4*AbsSqr(LamSD) - 2*AbsSqr(LamTD) + 0.6*Sqr(g1) + Sqr(g2))*ZP(gt1,0)*ZP(gt3,0) - (0.6*Sqr(g1) + Sqr(g2))*ZP(gt1,1)*ZP(gt3,1) - 4*AbsSqr(LamTD)*ZP(gt1,2)*ZP(gt3,2) + 2*Sqr(g2)*ZP(gt1,2)*ZP(gt3,2) - 2*Sqr(g2)*ZP(gt1,3)*ZP(gt3,3));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::SRum, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::SRum>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-((0.6*Sqr(g1) + Sqr(g2))*ZP(gt1,0)*ZP(gt3,0)) + (-4*AbsSqr(LamSU) - 2*AbsSqr(LamTU) + 0.6*Sqr(g1) + Sqr(g2))*ZP(gt1,1)*ZP(gt3,1) - 2*Sqr(g2)*ZP(gt1,2)*ZP(gt3,2) - 4*AbsSqr(LamTU)*ZP(gt1,3)*ZP(gt3,3) + 2*Sqr(g2)*ZP(gt1,3)*ZP(gt3,3));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::Su, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.05*(-4*(Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt4,3 + j1))*(ZP(gt1,0)*ZP(gt3,0) - ZP(gt1,1)*ZP(gt3,1)) + 5*(SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZU(gt4,j3))*ZP(gt1,0)*ZP(gt3,0) + SUM(j3,0,2,Conj(ZU(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))*ZU(gt4,3 + j2)))*ZP(gt1,1)*ZP(gt3,1))) + SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt4,j1))*((Sqr(g1) + 5*Sqr(g2))*ZP(gt1,0)*ZP(gt3,0) - (Sqr(g1) + 5*Sqr(g2))*ZP(gt1,1)*ZP(gt3,1) + 10*Sqr(g2)*(ZP(gt1,2)*ZP(gt3,2) - ZP(gt1,3)*ZP(gt3,3))));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::Sv, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -(SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZV(gt4,j3))*ZP(gt1,0)*ZP(gt3,0)) + 0.05*KroneckerDelta(gt2,gt4)*((-3*Sqr(g1) + 5*Sqr(g2))*ZP(gt1,0)*ZP(gt3,0) + (3*Sqr(g1) - 5*Sqr(g2))*ZP(gt1,1)*ZP(gt3,1) + 10*Sqr(g2)*(ZP(gt1,2)*ZP(gt3,2) - ZP(gt1,3)*ZP(gt3,3)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, fields::VP, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(Sqr(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*ZP(gt1,0)*ZP(gt2,0) + Sqr(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*ZP(gt1,1)*ZP(gt2,1) + 4*Sqr(g2)*Sqr(Sin(ThetaW))*(ZP(gt1,2)*ZP(gt2,2) + ZP(gt1,3)*ZP(gt2,3)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::VWm>::type, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.5*Sqr(g2)*(ZP(gt1,0)*ZP(gt2,0) + ZP(gt1,1)*ZP(gt2,1) + 2*ZP(gt1,2)*ZP(gt2,2) + 2*ZP(gt1,3)*ZP(gt2,3));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(Sqr(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*ZP(gt1,0)*ZP(gt2,0) + Sqr(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*ZP(gt1,1)*ZP(gt2,1) + 4*Sqr(g2)*Sqr(Cos(ThetaW))*(ZP(gt1,2)*ZP(gt2,2) + ZP(gt1,3)*ZP(gt2,3)));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::SRdp, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = 0.05*((-20*AbsSqr(LamTD) + 3*Sqr(g1) - 5*Sqr(g2))*ZA(gt1,0)*ZA(gt2,0) + (-3*Sqr(g1) + 5*Sqr(g2))*ZA(gt1,1)*ZA(gt2,1) + 5*Conj(LamTD)*(1.4142135623730951*LamSD*ZA(gt1,2)*ZA(gt2,3) + ZA(gt1,3)*(1.4142135623730951*LamSD*ZA(gt2,2) - 2*LamTD*ZA(gt2,3))) + 5*Conj(LamSD)*(1.4142135623730951*LamTD*ZA(gt1,3)*ZA(gt2,2) + ZA(gt1,2)*(-4*LamSD*ZA(gt2,2) + 1.4142135623730951*LamTD*ZA(gt2,3))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::SRdp, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.05*((-20*AbsSqr(LamTD) + 3*Sqr(g1) - 5*Sqr(g2))*ZH(gt1,0)*ZH(gt2,0) + (-3*Sqr(g1) + 5*Sqr(g2))*ZH(gt1,1)*ZH(gt2,1) + 5*Conj(LamTD)*(1.4142135623730951*LamSD*ZH(gt1,2)*ZH(gt2,3) + ZH(gt1,3)*(1.4142135623730951*LamSD*ZH(gt2,2) - 2*LamTD*ZH(gt2,3))) + 5*Conj(LamSD)*(1.4142135623730951*LamTD*ZH(gt1,3)*ZH(gt2,2) + ZH(gt1,2)*(-4*LamSD*ZH(gt2,2) + 1.4142135623730951*LamTD*ZH(gt2,3))));

   return {result};
}

ScalarVertex VertexImpl<fields::Rh, fields::SRdp, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = 0.05*(-((3*Sqr(g1) + 5*Sqr(g2))*ZHR(gt1,0)*ZHR(gt3,0)) + (3*Sqr(g1) - 5*Sqr(g2))*ZHR(gt1,1)*ZHR(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Sd, fields::SRdp, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = 0.05*(-((Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))) - 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Se, fields::SRdp, typename fields::conj<fields::Se>::type, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZE = MODELPARAMETER(ZE);

   const std::complex<double> result = 0.05*((3*Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1)) - 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1)));

   return {result};
}

ScalarVertex VertexImpl<fields::SRdp, fields::SRdp, typename fields::conj<fields::SRdp>::type, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> result = 0.1*(-3*Sqr(g1) - 5*Sqr(g2));

   return {result};
}

ScalarVertex VertexImpl<fields::SRdp, fields::SRum, typename fields::conj<fields::SRdp>::type, typename fields::conj<fields::SRum>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> result = 0.25*(0.6*Sqr(g1) + Sqr(g2));

   return {result};
}

ScalarVertex VertexImpl<fields::SRdp, fields::Su, typename fields::conj<fields::SRdp>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt4 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZU = MODELPARAMETER(ZU);

   const std::complex<double> result = 0.05*(-((Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt4,j1))) + 4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt4,3 + j1)));

   return {result};
}

ScalarVertex VertexImpl<fields::SRdp, fields::Sv, typename fields::conj<fields::SRdp>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt4 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> result = 0.05*KroneckerDelta(gt2,gt4)*(3*Sqr(g1) - 5*Sqr(g2));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SRdp, typename fields::conj<fields::SRdp>::type, fields::VP, fields::VP>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*Sqr(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SRdp, typename fields::conj<fields::SRdp>::type, typename fields::conj<fields::VWm>::type, fields::VWm>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> result = 0.5*Sqr(g2);

   return {result};
}

InverseMetricVertex VertexImpl<fields::SRdp, typename fields::conj<fields::SRdp>::type, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*Sqr(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Rh, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,-0.25)*(1.4142135623730951*LamTU*Conj(LamTD)*ZA(gt1,1)*ZHR(gt2,1)*ZP(gt3,0) - 1.4142135623730951*Sqr(g2)*ZA(gt1,1)*ZHR(gt2,0)*ZP(gt3,1) + 1.4142135623730951*ZA(gt1,0)*((-2*AbsSqr(LamSD) + AbsSqr(LamTD) + Sqr(g2))*ZHR(gt2,0)*ZP(gt3,0) - 2*LamTU*Conj(LamTD)*ZHR(gt2,1)*ZP(gt3,1)) - 2.8284271247461903*LamSD*Conj(LamTD)*ZA(gt1,2)*ZHR(gt2,0)*ZP(gt3,2) - 2*AbsSqr(LamTD)*ZA(gt1,3)*ZHR(gt2,0)*ZP(gt3,2) + 2*Sqr(g2)*ZA(gt1,3)*ZHR(gt2,0)*ZP(gt3,2) - 2*AbsSqr(LamTD)*ZA(gt1,3)*ZHR(gt2,0)*ZP(gt3,3) + 2*Sqr(g2)*ZA(gt1,3)*ZHR(gt2,0)*ZP(gt3,3) + 2.8284271247461903*Conj(LamSD)*(LamSU*ZA(gt1,1)*ZHR(gt2,1)*ZP(gt3,0) + LamTD*ZA(gt1,2)*ZHR(gt2,0)*ZP(gt3,3)));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Rh, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-1.4142135623730951*LamTU*Conj(LamTD)*ZH(gt1,1)*ZHR(gt2,1)*ZP(gt3,0) - 1.4142135623730951*Sqr(g2)*ZH(gt1,1)*ZHR(gt2,0)*ZP(gt3,1) - 1.4142135623730951*ZH(gt1,0)*((-2*AbsSqr(LamSD) + AbsSqr(LamTD) + Sqr(g2))*ZHR(gt2,0)*ZP(gt3,0) + 2*LamTU*Conj(LamTD)*ZHR(gt2,1)*ZP(gt3,1)) + 2.8284271247461903*LamSD*Conj(LamTD)*ZH(gt1,2)*ZHR(gt2,0)*ZP(gt3,2) + 2*AbsSqr(LamTD)*ZH(gt1,3)*ZHR(gt2,0)*ZP(gt3,2) - 2*Sqr(g2)*ZH(gt1,3)*ZHR(gt2,0)*ZP(gt3,2) - 2*AbsSqr(LamTD)*ZH(gt1,3)*ZHR(gt2,0)*ZP(gt3,3) + 2*Sqr(g2)*ZH(gt1,3)*ZHR(gt2,0)*ZP(gt3,3) + 2.8284271247461903*Conj(LamSD)*(-(LamSU*ZH(gt1,1)*ZHR(gt2,1)*ZP(gt3,0)) + LamTD*ZH(gt1,2)*ZHR(gt2,0)*ZP(gt3,3)));

   return {result};
}

ScalarVertex VertexImpl<fields::SRum, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt3 = indices[1];
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.5*(2*LamSU*Conj(LamSD) - LamTU*Conj(LamTD))*(ZP(gt2,1)*ZP(gt3,0) + ZP(gt2,0)*ZP(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Ah, fields::hh>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = std::complex<double>(0.,-0.35355339059327373)*(LamTD*Conj(LamSD)*(ZA(gt1,3)*ZA(gt2,0)*ZA(gt3,0)*ZH(gt4,2) - ZA(gt1,2)*ZA(gt2,0)*ZA(gt3,0)*ZH(gt4,3) + ZA(gt1,0)*(ZA(gt2,3)*ZA(gt3,0)*ZH(gt4,2) - ZA(gt2,2)*ZA(gt3,0)*ZH(gt4,3) + ZA(gt2,0)*(ZA(gt3,3)*ZH(gt4,2) - ZA(gt3,2)*ZH(gt4,3)))) - (LamTU*Conj(LamSU) - LamSU*Conj(LamTU))*(ZA(gt1,3)*ZA(gt2,1)*ZA(gt3,1)*ZH(gt4,2) - ZA(gt1,2)*ZA(gt2,1)*ZA(gt3,1)*ZH(gt4,3) + ZA(gt1,1)*(ZA(gt2,3)*ZA(gt3,1)*ZH(gt4,2) - ZA(gt2,2)*ZA(gt3,1)*ZH(gt4,3) + ZA(gt2,1)*(ZA(gt3,3)*ZH(gt4,2) - ZA(gt3,2)*ZH(gt4,3)))) + LamSD*Conj(LamTD)*(-(ZA(gt1,3)*ZA(gt2,0)*ZA(gt3,0)*ZH(gt4,2)) + ZA(gt1,2)*ZA(gt2,0)*ZA(gt3,0)*ZH(gt4,3) + ZA(gt1,0)*(-(ZA(gt2,3)*ZA(gt3,0)*ZH(gt4,2)) + ZA(gt2,2)*ZA(gt3,0)*ZH(gt4,3) + ZA(gt2,0)*(-(ZA(gt3,3)*ZH(gt4,2)) + ZA(gt3,2)*ZH(gt4,3)))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::hh, fields::hh, fields::hh>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = std::complex<double>(0.,-0.35355339059327373)*(LamTD*Conj(LamSD)*(ZA(gt1,3)*(ZH(gt2,2)*ZH(gt3,0)*ZH(gt4,0) + ZH(gt2,0)*(ZH(gt3,2)*ZH(gt4,0) + ZH(gt3,0)*ZH(gt4,2))) - ZA(gt1,2)*(ZH(gt2,3)*ZH(gt3,0)*ZH(gt4,0) + ZH(gt2,0)*(ZH(gt3,3)*ZH(gt4,0) + ZH(gt3,0)*ZH(gt4,3)))) + LamSD*Conj(LamTD)*(-(ZA(gt1,3)*(ZH(gt2,2)*ZH(gt3,0)*ZH(gt4,0) + ZH(gt2,0)*(ZH(gt3,2)*ZH(gt4,0) + ZH(gt3,0)*ZH(gt4,2)))) + ZA(gt1,2)*(ZH(gt2,3)*ZH(gt3,0)*ZH(gt4,0) + ZH(gt2,0)*(ZH(gt3,3)*ZH(gt4,0) + ZH(gt3,0)*ZH(gt4,3)))) - (LamTU*Conj(LamSU) - LamSU*Conj(LamTU))*(ZA(gt1,3)*(ZH(gt2,2)*ZH(gt3,1)*ZH(gt4,1) + ZH(gt2,1)*(ZH(gt3,2)*ZH(gt4,1) + ZH(gt3,1)*ZH(gt4,2))) - ZA(gt1,2)*(ZH(gt2,3)*ZH(gt3,1)*ZH(gt4,1) + ZH(gt2,1)*(ZH(gt3,3)*ZH(gt4,1) + ZH(gt3,1)*ZH(gt4,3)))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::hh, fields::Hpm, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,0.25)*(-(Sqr(g2)*ZA(gt1,1)*ZH(gt2,0)*ZP(gt3,1)*ZP(gt4,0)) - Sqr(g2)*ZA(gt1,0)*ZH(gt2,1)*ZP(gt3,1)*ZP(gt4,0) + 1.4142135623730951*Sqr(g2)*ZA(gt1,3)*ZH(gt2,0)*ZP(gt3,2)*ZP(gt4,0) - 1.4142135623730951*Sqr(g2)*ZA(gt1,0)*ZH(gt2,3)*ZP(gt3,2)*ZP(gt4,0) + 1.4142135623730951*Sqr(g2)*ZA(gt1,3)*ZH(gt2,0)*ZP(gt3,3)*ZP(gt4,0) + 1.4142135623730951*Sqr(g2)*ZA(gt1,0)*ZH(gt2,3)*ZP(gt3,3)*ZP(gt4,0) + Sqr(g2)*ZA(gt1,1)*ZH(gt2,0)*ZP(gt3,0)*ZP(gt4,1) + Sqr(g2)*ZA(gt1,0)*ZH(gt2,1)*ZP(gt3,0)*ZP(gt4,1) - 1.4142135623730951*LamTU*Conj(LamSU)*ZA(gt1,3)*ZH(gt2,2)*ZP(gt3,1)*ZP(gt4,1) + 1.4142135623730951*LamSU*Conj(LamTU)*ZA(gt1,3)*ZH(gt2,2)*ZP(gt3,1)*ZP(gt4,1) + 1.4142135623730951*LamTU*Conj(LamSU)*ZA(gt1,2)*ZH(gt2,3)*ZP(gt3,1)*ZP(gt4,1) - 1.4142135623730951*LamSU*Conj(LamTU)*ZA(gt1,2)*ZH(gt2,3)*ZP(gt3,1)*ZP(gt4,1) + 2*LamTU*Conj(LamSU)*ZA(gt1,2)*ZH(gt2,1)*ZP(gt3,2)*ZP(gt4,1) - 1.4142135623730951*AbsSqr(LamTU)*ZA(gt1,3)*ZH(gt2,1)*ZP(gt3,2)*ZP(gt4,1) + 1.4142135623730951*Sqr(g2)*ZA(gt1,3)*ZH(gt2,1)*ZP(gt3,2)*ZP(gt4,1) + 2*LamTU*Conj(LamSU)*ZA(gt1,1)*ZH(gt2,2)*ZP(gt3,2)*ZP(gt4,1) - 1.4142135623730951*AbsSqr(LamTU)*ZA(gt1,1)*ZH(gt2,3)*ZP(gt3,2)*ZP(gt4,1) + 1.4142135623730951*Sqr(g2)*ZA(gt1,1)*ZH(gt2,3)*ZP(gt3,2)*ZP(gt4,1) - 2*LamSU*Conj(LamTU)*ZA(gt1,2)*ZH(gt2,1)*ZP(gt3,3)*ZP(gt4,1) - 1.4142135623730951*AbsSqr(LamTU)*ZA(gt1,3)*ZH(gt2,1)*ZP(gt3,3)*ZP(gt4,1) + 1.4142135623730951*Sqr(g2)*ZA(gt1,3)*ZH(gt2,1)*ZP(gt3,3)*ZP(gt4,1) + 2*LamSU*Conj(LamTU)*ZA(gt1,1)*ZH(gt2,2)*ZP(gt3,3)*ZP(gt4,1) + 1.4142135623730951*AbsSqr(LamTU)*ZA(gt1,1)*ZH(gt2,3)*ZP(gt3,3)*ZP(gt4,1) - 1.4142135623730951*Sqr(g2)*ZA(gt1,1)*ZH(gt2,3)*ZP(gt3,3)*ZP(gt4,1) - 1.4142135623730951*Sqr(g2)*ZA(gt1,3)*ZH(gt2,0)*ZP(gt3,0)*ZP(gt4,2) + 1.4142135623730951*Sqr(g2)*ZA(gt1,0)*ZH(gt2,3)*ZP(gt3,0)*ZP(gt4,2) - 2*LamSU*Conj(LamTU)*ZA(gt1,2)*ZH(gt2,1)*ZP(gt3,1)*ZP(gt4,2) + 1.4142135623730951*AbsSqr(LamTU)*ZA(gt1,3)*ZH(gt2,1)*ZP(gt3,1)*ZP(gt4,2) - 1.4142135623730951*Sqr(g2)*ZA(gt1,3)*ZH(gt2,1)*ZP(gt3,1)*ZP(gt4,2) - 2*LamSU*Conj(LamTU)*ZA(gt1,1)*ZH(gt2,2)*ZP(gt3,1)*ZP(gt4,2) + 1.4142135623730951*AbsSqr(LamTU)*ZA(gt1,1)*ZH(gt2,3)*ZP(gt3,1)*ZP(gt4,2) - 1.4142135623730951*Sqr(g2)*ZA(gt1,1)*ZH(gt2,3)*ZP(gt3,1)*ZP(gt4,2) + 4*Sqr(g2)*ZA(gt1,3)*ZH(gt2,3)*ZP(gt3,3)*ZP(gt4,2) - 1.4142135623730951*Sqr(g2)*ZA(gt1,3)*ZH(gt2,0)*ZP(gt3,0)*ZP(gt4,3) - 1.4142135623730951*Sqr(g2)*ZA(gt1,0)*ZH(gt2,3)*ZP(gt3,0)*ZP(gt4,3) + 2*LamTU*Conj(LamSU)*ZA(gt1,2)*ZH(gt2,1)*ZP(gt3,1)*ZP(gt4,3) + 1.4142135623730951*AbsSqr(LamTU)*ZA(gt1,3)*ZH(gt2,1)*ZP(gt3,1)*ZP(gt4,3) - 1.4142135623730951*Sqr(g2)*ZA(gt1,3)*ZH(gt2,1)*ZP(gt3,1)*ZP(gt4,3) - 2*LamTU*Conj(LamSU)*ZA(gt1,1)*ZH(gt2,2)*ZP(gt3,1)*ZP(gt4,3) - 1.4142135623730951*AbsSqr(LamTU)*ZA(gt1,1)*ZH(gt2,3)*ZP(gt3,1)*ZP(gt4,3) + 1.4142135623730951*Sqr(g2)*ZA(gt1,1)*ZH(gt2,3)*ZP(gt3,1)*ZP(gt4,3) - 4*Sqr(g2)*ZA(gt1,3)*ZH(gt2,3)*ZP(gt3,2)*ZP(gt4,3) + LamTD*Conj(LamSD)*(1.4142135623730951*ZA(gt1,3)*ZH(gt2,2)*ZP(gt3,0)*ZP(gt4,0) + 2*ZA(gt1,0)*ZH(gt2,2)*(-(ZP(gt3,2)*ZP(gt4,0)) + ZP(gt3,0)*ZP(gt4,3)) + ZA(gt1,2)*(-1.4142135623730951*ZH(gt2,3)*ZP(gt3,0)*ZP(gt4,0) + 2*ZH(gt2,0)*(ZP(gt3,2)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,3)))) - Conj(LamTD)*(LamSD*ZA(gt1,2)*(-1.4142135623730951*ZH(gt2,3)*ZP(gt3,0)*ZP(gt4,0) + 2*ZH(gt2,0)*(ZP(gt3,3)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,2))) + ZA(gt1,0)*(2*LamSD*ZH(gt2,2)*(ZP(gt3,3)*ZP(gt4,0) - ZP(gt3,0)*ZP(gt4,2)) + 1.4142135623730951*LamTD*ZH(gt2,3)*(-(ZP(gt3,2)*ZP(gt4,0)) + ZP(gt3,3)*ZP(gt4,0) + ZP(gt3,0)*(ZP(gt4,2) - ZP(gt4,3)))) + 1.4142135623730951*ZA(gt1,3)*(LamSD*ZH(gt2,2)*ZP(gt3,0)*ZP(gt4,0) + LamTD*ZH(gt2,0)*(ZP(gt3,2)*ZP(gt4,0) + ZP(gt3,3)*ZP(gt4,0) - ZP(gt3,0)*(ZP(gt4,2) + ZP(gt4,3))))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::hh, fields::Rh, typename fields::conj<fields::Rh>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = std::complex<double>(0,0.25)*(Conj(LamSD)*(-1.4142135623730951*LamTD*ZA(gt1,3)*ZH(gt2,2)*ZHR(gt3,0) + 1.4142135623730951*LamTD*ZA(gt1,2)*ZH(gt2,3)*ZHR(gt3,0) + 2*LamSU*(ZA(gt1,1)*ZH(gt2,0) - ZA(gt1,0)*ZH(gt2,1))*ZHR(gt3,1))*ZHR(gt4,0) + Conj(LamTD)*(1.4142135623730951*LamSD*ZA(gt1,3)*ZH(gt2,2)*ZHR(gt3,0) - 1.4142135623730951*LamSD*ZA(gt1,2)*ZH(gt2,3)*ZHR(gt3,0) + LamTU*(-(ZA(gt1,1)*ZH(gt2,0)) + ZA(gt1,0)*ZH(gt2,1))*ZHR(gt3,1))*ZHR(gt4,0) + (Conj(LamSU)*(-2*LamSD*ZA(gt1,1)*ZH(gt2,0)*ZHR(gt3,0) + 2*LamSD*ZA(gt1,0)*ZH(gt2,1)*ZHR(gt3,0) + 1.4142135623730951*LamTU*(ZA(gt1,3)*ZH(gt2,2) - ZA(gt1,2)*ZH(gt2,3))*ZHR(gt3,1)) + Conj(LamTU)*(LamTD*ZA(gt1,1)*ZH(gt2,0)*ZHR(gt3,0) - LamTD*ZA(gt1,0)*ZH(gt2,1)*ZHR(gt3,0) + 1.4142135623730951*LamSU*(-(ZA(gt1,3)*ZH(gt2,2)) + ZA(gt1,2)*ZH(gt2,3))*ZHR(gt3,1)))*ZHR(gt4,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::hh, fields::SRdp, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = std::complex<double>(0.,0.35355339059327373)*(LamTD*Conj(LamSD) - LamSD*Conj(LamTD))*(ZA(gt1,3)*ZH(gt2,2) - ZA(gt1,2)*ZH(gt2,3));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::hh, fields::SRum, typename fields::conj<fields::SRum>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = std::complex<double>(0.,-0.35355339059327373)*(LamTU*Conj(LamSU) - LamSU*Conj(LamTU))*(ZA(gt1,3)*ZH(gt2,2) - ZA(gt1,2)*ZH(gt2,3));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::hh, fields::hh>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.05*(-5*(Conj(LamSD)*(1.4142135623730951*LamTD*ZA(gt1,3)*ZA(gt2,2) + ZA(gt1,2)*(4*LamSD*ZA(gt2,2) + 1.4142135623730951*LamTD*ZA(gt2,3)))*ZH(gt3,0)*ZH(gt4,0) + Conj(LamTD)*(1.4142135623730951*LamSD*ZA(gt1,2)*ZA(gt2,3) + ZA(gt1,3)*(1.4142135623730951*LamSD*ZA(gt2,2) + 2*LamTD*ZA(gt2,3)))*ZH(gt3,0)*ZH(gt4,0) + (-(Conj(LamTU)*(1.4142135623730951*LamSU*ZA(gt1,2)*ZA(gt2,3) + ZA(gt1,3)*(1.4142135623730951*LamSU*ZA(gt2,2) - 2*LamTU*ZA(gt2,3)))) - Conj(LamSU)*(1.4142135623730951*LamTU*ZA(gt1,3)*ZA(gt2,2) + ZA(gt1,2)*(-4*LamSU*ZA(gt2,2) + 1.4142135623730951*LamTU*ZA(gt2,3))))*ZH(gt3,1)*ZH(gt4,1)) - ZA(gt1,0)*ZA(gt2,0)*((3*Sqr(g1) + 5*Sqr(g2))*ZH(gt3,0)*ZH(gt4,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZH(gt3,1)*ZH(gt4,1) + 5*Conj(LamSD)*(1.4142135623730951*LamTD*ZH(gt3,3)*ZH(gt4,2) + ZH(gt3,2)*(4*LamSD*ZH(gt4,2) + 1.4142135623730951*LamTD*ZH(gt4,3))) + 5*Conj(LamTD)*(1.4142135623730951*LamSD*ZH(gt3,2)*ZH(gt4,3) + ZH(gt3,3)*(1.4142135623730951*LamSD*ZH(gt4,2) + 2*LamTD*ZH(gt4,3)))) + ZA(gt1,1)*ZA(gt2,1)*((3*Sqr(g1) + 5*Sqr(g2))*ZH(gt3,0)*ZH(gt4,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZH(gt3,1)*ZH(gt4,1) + 5*Conj(LamTU)*(1.4142135623730951*LamSU*ZH(gt3,2)*ZH(gt4,3) + ZH(gt3,3)*(1.4142135623730951*LamSU*ZH(gt4,2) - 2*LamTU*ZH(gt4,3))) + 5*Conj(LamSU)*(1.4142135623730951*LamTU*ZH(gt3,3)*ZH(gt4,2) + ZH(gt3,2)*(-4*LamSU*ZH(gt4,2) + 1.4142135623730951*LamTU*ZH(gt4,3)))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::hh, fields::hh>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.05*(-(ZH(gt1,0)*(-((3*Sqr(g1) + 5*Sqr(g2))*ZH(gt2,1)*(ZH(gt3,1)*ZH(gt4,0) + ZH(gt3,0)*ZH(gt4,1))) + 5*Conj(LamSD)*(1.4142135623730951*LamTD*ZH(gt2,3)*(ZH(gt3,2)*ZH(gt4,0) + ZH(gt3,0)*ZH(gt4,2)) + ZH(gt2,2)*(4*LamSD*ZH(gt3,2)*ZH(gt4,0) + 1.4142135623730951*LamTD*ZH(gt3,3)*ZH(gt4,0) + ZH(gt3,0)*(4*LamSD*ZH(gt4,2) + 1.4142135623730951*LamTD*ZH(gt4,3)))) + 5*Conj(LamTD)*(1.4142135623730951*LamSD*ZH(gt2,2)*(ZH(gt3,3)*ZH(gt4,0) + ZH(gt3,0)*ZH(gt4,3)) + ZH(gt2,3)*(1.4142135623730951*LamSD*ZH(gt3,2)*ZH(gt4,0) + 2*LamTD*ZH(gt3,3)*ZH(gt4,0) + ZH(gt3,0)*(1.4142135623730951*LamSD*ZH(gt4,2) + 2*LamTD*ZH(gt4,3)))) + ZH(gt2,0)*(3*(3*Sqr(g1) + 5*Sqr(g2))*ZH(gt3,0)*ZH(gt4,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZH(gt3,1)*ZH(gt4,1) + 5*Conj(LamSD)*(1.4142135623730951*LamTD*ZH(gt3,3)*ZH(gt4,2) + ZH(gt3,2)*(4*LamSD*ZH(gt4,2) + 1.4142135623730951*LamTD*ZH(gt4,3))) + 5*Conj(LamTD)*(1.4142135623730951*LamSD*ZH(gt3,2)*ZH(gt4,3) + ZH(gt3,3)*(1.4142135623730951*LamSD*ZH(gt4,2) + 2*LamTD*ZH(gt4,3)))))) + ZH(gt1,1)*((3*Sqr(g1) + 5*Sqr(g2))*ZH(gt2,0)*(ZH(gt3,1)*ZH(gt4,0) + ZH(gt3,0)*ZH(gt4,1)) + 5*Conj(LamTU)*(1.4142135623730951*LamSU*ZH(gt2,2)*(ZH(gt3,3)*ZH(gt4,1) + ZH(gt3,1)*ZH(gt4,3)) + ZH(gt2,3)*(1.4142135623730951*LamSU*ZH(gt3,2)*ZH(gt4,1) - 2*LamTU*ZH(gt3,3)*ZH(gt4,1) + ZH(gt3,1)*(1.4142135623730951*LamSU*ZH(gt4,2) - 2*LamTU*ZH(gt4,3)))) + 5*Conj(LamSU)*(1.4142135623730951*LamTU*ZH(gt2,3)*(ZH(gt3,2)*ZH(gt4,1) + ZH(gt3,1)*ZH(gt4,2)) + ZH(gt2,2)*(-4*LamSU*ZH(gt3,2)*ZH(gt4,1) + 1.4142135623730951*LamTU*ZH(gt3,3)*ZH(gt4,1) + ZH(gt3,1)*(-4*LamSU*ZH(gt4,2) + 1.4142135623730951*LamTU*ZH(gt4,3)))) + ZH(gt2,1)*((3*Sqr(g1) + 5*Sqr(g2))*ZH(gt3,0)*ZH(gt4,0) - 3*(3*Sqr(g1) + 5*Sqr(g2))*ZH(gt3,1)*ZH(gt4,1) + 5*Conj(LamTU)*(1.4142135623730951*LamSU*ZH(gt3,2)*ZH(gt4,3) + ZH(gt3,3)*(1.4142135623730951*LamSU*ZH(gt4,2) - 2*LamTU*ZH(gt4,3))) + 5*Conj(LamSU)*(1.4142135623730951*LamTU*ZH(gt3,3)*ZH(gt4,2) + ZH(gt3,2)*(-4*LamSU*ZH(gt4,2) + 1.4142135623730951*LamTU*ZH(gt4,3))))) + 5*(-4*AbsSqr(LamSU)*ZH(gt1,2)*ZH(gt2,2)*ZH(gt3,1)*ZH(gt4,1) + 1.4142135623730951*LamTU*Conj(LamSU)*ZH(gt1,3)*ZH(gt2,2)*ZH(gt3,1)*ZH(gt4,1) + 1.4142135623730951*LamSU*Conj(LamTU)*ZH(gt1,3)*ZH(gt2,2)*ZH(gt3,1)*ZH(gt4,1) + 1.4142135623730951*LamTU*Conj(LamSU)*ZH(gt1,2)*ZH(gt2,3)*ZH(gt3,1)*ZH(gt4,1) + 1.4142135623730951*LamSU*Conj(LamTU)*ZH(gt1,2)*ZH(gt2,3)*ZH(gt3,1)*ZH(gt4,1) - 2*AbsSqr(LamTU)*ZH(gt1,3)*ZH(gt2,3)*ZH(gt3,1)*ZH(gt4,1) - 4*AbsSqr(LamSU)*ZH(gt1,2)*ZH(gt2,1)*ZH(gt3,2)*ZH(gt4,1) + 1.4142135623730951*LamTU*Conj(LamSU)*ZH(gt1,3)*ZH(gt2,1)*ZH(gt3,2)*ZH(gt4,1) + 1.4142135623730951*LamSU*Conj(LamTU)*ZH(gt1,3)*ZH(gt2,1)*ZH(gt3,2)*ZH(gt4,1) + 1.4142135623730951*LamTU*Conj(LamSU)*ZH(gt1,2)*ZH(gt2,1)*ZH(gt3,3)*ZH(gt4,1) + 1.4142135623730951*LamSU*Conj(LamTU)*ZH(gt1,2)*ZH(gt2,1)*ZH(gt3,3)*ZH(gt4,1) - 2*AbsSqr(LamTU)*ZH(gt1,3)*ZH(gt2,1)*ZH(gt3,3)*ZH(gt4,1) - 4*AbsSqr(LamSU)*ZH(gt1,2)*ZH(gt2,1)*ZH(gt3,1)*ZH(gt4,2) + 1.4142135623730951*LamTU*Conj(LamSU)*ZH(gt1,3)*ZH(gt2,1)*ZH(gt3,1)*ZH(gt4,2) + 1.4142135623730951*LamSU*Conj(LamTU)*ZH(gt1,3)*ZH(gt2,1)*ZH(gt3,1)*ZH(gt4,2) + 1.4142135623730951*LamTU*Conj(LamSU)*ZH(gt1,2)*ZH(gt2,1)*ZH(gt3,1)*ZH(gt4,3) + 1.4142135623730951*LamSU*Conj(LamTU)*ZH(gt1,2)*ZH(gt2,1)*ZH(gt3,1)*ZH(gt4,3) - 2*AbsSqr(LamTU)*ZH(gt1,3)*ZH(gt2,1)*ZH(gt3,1)*ZH(gt4,3) - Conj(LamSD)*(1.4142135623730951*LamTD*ZH(gt1,3)*(ZH(gt2,2)*ZH(gt3,0)*ZH(gt4,0) + ZH(gt2,0)*(ZH(gt3,2)*ZH(gt4,0) + ZH(gt3,0)*ZH(gt4,2))) + ZH(gt1,2)*(4*LamSD*ZH(gt2,2)*ZH(gt3,0)*ZH(gt4,0) + 1.4142135623730951*LamTD*ZH(gt2,3)*ZH(gt3,0)*ZH(gt4,0) + ZH(gt2,0)*(4*LamSD*ZH(gt3,2)*ZH(gt4,0) + 1.4142135623730951*LamTD*ZH(gt3,3)*ZH(gt4,0) + ZH(gt3,0)*(4*LamSD*ZH(gt4,2) + 1.4142135623730951*LamTD*ZH(gt4,3))))) - Conj(LamTD)*(1.4142135623730951*LamSD*ZH(gt1,2)*(ZH(gt2,3)*ZH(gt3,0)*ZH(gt4,0) + ZH(gt2,0)*(ZH(gt3,3)*ZH(gt4,0) + ZH(gt3,0)*ZH(gt4,3))) + ZH(gt1,3)*(1.4142135623730951*LamSD*ZH(gt2,2)*ZH(gt3,0)*ZH(gt4,0) + 2*LamTD*ZH(gt2,3)*ZH(gt3,0)*ZH(gt4,0) + ZH(gt2,0)*(1.4142135623730951*LamSD*ZH(gt3,2)*ZH(gt4,0) + 2*LamTD*ZH(gt3,3)*ZH(gt4,0) + ZH(gt3,0)*(1.4142135623730951*LamSD*ZH(gt4,2) + 2*LamTD*ZH(gt4,3)))))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::Rh, typename fields::conj<fields::Rh>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = 0.05*(-5*(Conj(LamSD)*(1.4142135623730951*LamTD*ZH(gt1,3)*ZH(gt2,2) + ZH(gt1,2)*(4*LamSD*ZH(gt2,2) + 1.4142135623730951*LamTD*ZH(gt2,3)))*ZHR(gt3,0)*ZHR(gt4,0) + Conj(LamTD)*(1.4142135623730951*LamSD*ZH(gt1,2)*ZH(gt2,3) + ZH(gt1,3)*(1.4142135623730951*LamSD*ZH(gt2,2) + 2*LamTD*ZH(gt2,3)))*ZHR(gt3,0)*ZHR(gt4,0) + (-(Conj(LamTU)*(1.4142135623730951*LamSU*ZH(gt1,2)*ZH(gt2,3) + ZH(gt1,3)*(1.4142135623730951*LamSU*ZH(gt2,2) - 2*LamTU*ZH(gt2,3)))) - Conj(LamSU)*(1.4142135623730951*LamTU*ZH(gt1,3)*ZH(gt2,2) + ZH(gt1,2)*(-4*LamSU*ZH(gt2,2) + 1.4142135623730951*LamTU*ZH(gt2,3))))*ZHR(gt3,1)*ZHR(gt4,1)) - ZH(gt1,1)*(5*ZH(gt2,0)*(-2*LamSU*Conj(LamSD)*ZHR(gt3,1)*ZHR(gt4,0) + LamTU*Conj(LamTD)*ZHR(gt3,1)*ZHR(gt4,0) + (-2*LamSD*Conj(LamSU) + LamTD*Conj(LamTU))*ZHR(gt3,0)*ZHR(gt4,1)) + ZH(gt2,1)*((3*Sqr(g1) + 5*Sqr(g2))*ZHR(gt3,0)*ZHR(gt4,0) + (20*AbsSqr(LamSU) + 10*AbsSqr(LamTU) - 3*Sqr(g1) - 5*Sqr(g2))*ZHR(gt3,1)*ZHR(gt4,1))) + ZH(gt1,0)*(5*ZH(gt2,1)*(2*LamSU*Conj(LamSD)*ZHR(gt3,1)*ZHR(gt4,0) - LamTU*Conj(LamTD)*ZHR(gt3,1)*ZHR(gt4,0) + (2*LamSD*Conj(LamSU) - LamTD*Conj(LamTU))*ZHR(gt3,0)*ZHR(gt4,1)) + ZH(gt2,0)*((-20*AbsSqr(LamSD) - 10*AbsSqr(LamTD) + 3*Sqr(g1) + 5*Sqr(g2))*ZHR(gt3,0)*ZHR(gt4,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZHR(gt3,1)*ZHR(gt4,1))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::Sd, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.05*(-20*(SUM(j3,0,2,Conj(ZD(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,j1))*ZD(gt4,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt3,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gt4,j3)))*ZH(gt1,0)*ZH(gt2,0) + (Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZD(gt3,j1))*ZD(gt4,j1))*(ZH(gt1,0)*ZH(gt2,0) - ZH(gt1,1)*ZH(gt2,1)) + 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt3,3 + j1))*ZD(gt4,3 + j1))*(ZH(gt1,0)*ZH(gt2,0) - ZH(gt1,1)*ZH(gt2,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::Se, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.05*(-20*(SUM(j3,0,2,Conj(ZE(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Ye(j2,j1))*ZE(gt4,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gt3,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZE(gt4,j3)))*ZH(gt1,0)*ZH(gt2,0) - (3*Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZE(gt3,j1))*ZE(gt4,j1))*(ZH(gt1,0)*ZH(gt2,0) - ZH(gt1,1)*ZH(gt2,1)) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*ZE(gt4,3 + j1))*(ZH(gt1,0)*ZH(gt2,0) - ZH(gt1,1)*ZH(gt2,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::SRum, typename fields::conj<fields::SRum>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.05*((-3*Sqr(g1) + 5*Sqr(g2))*ZH(gt1,0)*ZH(gt2,0) + (-20*AbsSqr(LamTU) + 3*Sqr(g1) - 5*Sqr(g2))*ZH(gt1,1)*ZH(gt2,1) - 5*(Conj(LamSU)*(1.4142135623730951*LamTU*ZH(gt1,3)*ZH(gt2,2) + ZH(gt1,2)*(4*LamSU*ZH(gt2,2) + 1.4142135623730951*LamTU*ZH(gt2,3))) + Conj(LamTU)*(1.4142135623730951*LamSU*ZH(gt1,2)*ZH(gt2,3) + ZH(gt1,3)*(1.4142135623730951*LamSU*ZH(gt2,2) + 2*LamTU*ZH(gt2,3)))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::Su, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.05*((Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZU(gt3,j1))*ZU(gt4,j1))*(ZH(gt1,0)*ZH(gt2,0) - ZH(gt1,1)*ZH(gt2,1)) - 4*(5*(SUM(j3,0,2,Conj(ZU(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))*ZU(gt4,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gt4,j3)))*ZH(gt1,1)*ZH(gt2,1) + Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt3,3 + j1))*ZU(gt4,3 + j1))*(ZH(gt1,0)*ZH(gt2,0) - ZH(gt1,1)*ZH(gt2,1))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::Sv, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = -0.25*KroneckerDelta(gt3,gt4)*(0.6*Sqr(g1) + Sqr(g2))*(ZH(gt1,0)*ZH(gt2,0) - ZH(gt1,1)*ZH(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::hh, typename fields::conj<fields::VWm>::type, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.5*Sqr(g2)*(ZH(gt1,0)*ZH(gt2,0) + ZH(gt1,1)*ZH(gt2,1) + 4*ZH(gt1,3)*ZH(gt2,3));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::hh, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*Sqr(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(ZH(gt1,0)*ZH(gt2,0) + ZH(gt1,1)*ZH(gt2,1));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Sd, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = 0.5*SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt4,3 + j1)))*(1.4142135623730951*Conj(LamSD)*ZH(gt1,2) + Conj(LamTD)*ZH(gt1,3))*ZHR(gt3,0);

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Se, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = 0.5*SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt4,3 + j1)))*(1.4142135623730951*Conj(LamSD)*ZH(gt1,2) + Conj(LamTD)*ZH(gt1,3))*ZHR(gt3,0);

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Su, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.5*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt4,3 + j1)))*(1.4142135623730951*Conj(LamSU)*ZH(gt1,2) - Conj(LamTU)*ZH(gt1,3))*ZHR(gt3,1);

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Sd, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.7071067811865475*(SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZU(gt4,j3))*ZH(gt1,0)*ZP(gt3,0) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gt4,j3))*ZH(gt1,1)*ZP(gt3,1) + SUM(j3,0,2,Conj(ZD(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yu(j2,j1))*ZU(gt4,3 + j2)))*(ZH(gt1,1)*ZP(gt3,0) + ZH(gt1,0)*ZP(gt3,1))) - 0.25*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZU(gt4,j1))*(1.4142135623730951*ZH(gt1,0)*ZP(gt3,0) + 1.4142135623730951*ZH(gt1,1)*ZP(gt3,1) + 2*ZH(gt1,3)*(ZP(gt3,2) - ZP(gt3,3)));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Se, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.7071067811865475*SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZV(gt4,j3))*ZH(gt1,0)*ZP(gt3,0) - 0.25*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZV(gt4,j1))*(1.4142135623730951*ZH(gt1,0)*ZP(gt3,0) + 1.4142135623730951*ZH(gt1,1)*ZP(gt3,1) + 2*ZH(gt1,3)*(ZP(gt3,2) - ZP(gt3,3)));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Sv, typename fields::conj<fields::Se>::type, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.5*SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)))*(1.4142135623730951*Conj(LamSD)*ZH(gt1,2) - Conj(LamTD)*ZH(gt1,3));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::Ah, fields::Hpm, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,-0.5)*g2*(ZA(gt1,0)*ZP(gt2,0) + ZA(gt1,1)*ZP(gt2,1) + 1.4142135623730951*ZA(gt1,3)*(ZP(gt2,2) - ZP(gt2,3)));

   return {result, minuend_index, subtrahend_index};
}

ScalarVertex VertexImpl<fields::Ah, fields::Hpm, fields::SRdp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Mu = MODELPARAMETER(Mu);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,-0.5)*Conj(Mu)*(1.4142135623730951*LamSD*ZA(gt1,2)*ZP(gt2,1) - LamTD*ZA(gt1,3)*ZP(gt2,1) + 1.4142135623730951*LamTD*ZA(gt1,1)*ZP(gt2,2));

   return {result};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Chi>::type, fields::Cha2, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto g1 = MODELPARAMETER(g1);
   const auto UM2 = MODELPARAMETER(UM2);
   const auto ZN2 = MODELPARAMETER(ZN2);
   const auto ZP = MODELPARAMETER(ZP);
   const auto UP2 = MODELPARAMETER(UP2);
   const auto ZN1 = MODELPARAMETER(ZN1);

   const std::complex<double> left = -(g2*Conj(UM2(gt2,0))*(Conj(ZN2(gt1,2))*ZP(gt3,0) + 1.4142135623730951*Conj(ZN2(gt1,1))*ZP(gt3,2))) + 0.5*Conj(UM2(gt2,1))*(2*LamSU*Conj(ZN2(gt1,0))*ZP(gt3,1) + 1.4142135623730951*LamTU*Conj(ZN2(gt1,1))*ZP(gt3,1) + 2*LamTU*Conj(ZN2(gt1,3))*ZP(gt3,3));

   const std::complex<double> right = Conj(LamTD)*UP2(gt2,0)*ZN1(gt1,2)*ZP(gt3,0) - 0.5*UP2(gt2,1)*(1.0954451150103321*g1*ZN1(gt1,0)*ZP(gt3,1) + 1.4142135623730951*g2*ZN1(gt1,1)*ZP(gt3,1) + 2*Conj(LamTU)*ZN1(gt1,3)*ZP(gt3,2)) - 1.4142135623730951*g2*UP2(gt2,0)*ZN1(gt1,1)*ZP(gt3,3);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha2>::type, fields::Chi, fields::Hpm>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto UP2 = MODELPARAMETER(UP2);
   const auto ZN1 = MODELPARAMETER(ZN1);
   const auto ZP = MODELPARAMETER(ZP);
   const auto UM2 = MODELPARAMETER(UM2);
   const auto ZN2 = MODELPARAMETER(ZN2);

   const std::complex<double> left = -0.5*Conj(UP2(gt1,1))*(1.0954451150103321*g1*Conj(ZN1(gt2,0))*ZP(gt3,1) + 1.4142135623730951*g2*Conj(ZN1(gt2,1))*ZP(gt3,1) + 2*LamTU*Conj(ZN1(gt2,3))*ZP(gt3,2)) + Conj(UP2(gt1,0))*(LamTD*Conj(ZN1(gt2,2))*ZP(gt3,0) - 1.4142135623730951*g2*Conj(ZN1(gt2,1))*ZP(gt3,3));

   const std::complex<double> right = -(g2*UM2(gt1,0)*(ZN2(gt2,2)*ZP(gt3,0) + 1.4142135623730951*ZN2(gt2,1)*ZP(gt3,2))) + 0.5*UM2(gt1,1)*(2*Conj(LamSU)*ZN2(gt2,0)*ZP(gt3,1) + Conj(LamTU)*(1.4142135623730951*ZN2(gt2,1)*ZP(gt3,1) + 2*ZN2(gt2,3)*ZP(gt3,3)));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha1>::type, fields::Chi, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto ZN1 = MODELPARAMETER(ZN1);
   const auto ZP = MODELPARAMETER(ZP);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto ZN2 = MODELPARAMETER(ZN2);

   const std::complex<double> left = Conj(UM1(gt1,0))*(-(LamTU*Conj(ZN1(gt2,3))*ZP(gt3,1)) + 1.4142135623730951*g2*Conj(ZN1(gt2,1))*ZP(gt3,2)) + Conj(UM1(gt1,1))*(0.5477225575051661*g1*Conj(ZN1(gt2,0))*ZP(gt3,0) + 0.7071067811865475*g2*Conj(ZN1(gt2,1))*ZP(gt3,0) + LamTD*Conj(ZN1(gt2,2))*ZP(gt3,3));

   const std::complex<double> right = -(Conj(LamSD)*UP1(gt1,1)*ZN2(gt2,0)*ZP(gt3,0)) + 0.5*Conj(LamTD)*UP1(gt1,1)*(1.4142135623730951*ZN2(gt2,1)*ZP(gt3,0) - 2*ZN2(gt2,2)*ZP(gt3,2)) + g2*UP1(gt1,0)*(-(ZN2(gt2,3)*ZP(gt3,1)) + 1.4142135623730951*ZN2(gt2,1)*ZP(gt3,3));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Chi>::type, fields::Cha1, fields::Hpm>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto g1 = MODELPARAMETER(g1);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto ZN2 = MODELPARAMETER(ZN2);
   const auto ZP = MODELPARAMETER(ZP);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto ZN1 = MODELPARAMETER(ZN1);

   const std::complex<double> left = Conj(UP1(gt2,1))*(-(LamSD*Conj(ZN2(gt1,0))*ZP(gt3,0)) + 0.7071067811865475*LamTD*Conj(ZN2(gt1,1))*ZP(gt3,0) - LamTD*Conj(ZN2(gt1,2))*ZP(gt3,2)) + g2*Conj(UP1(gt2,0))*(-(Conj(ZN2(gt1,3))*ZP(gt3,1)) + 1.4142135623730951*Conj(ZN2(gt1,1))*ZP(gt3,3));

   const std::complex<double> right = UM1(gt2,0)*(-(Conj(LamTU)*ZN1(gt1,3)*ZP(gt3,1)) + 1.4142135623730951*g2*ZN1(gt1,1)*ZP(gt3,2)) + UM1(gt2,1)*(0.5477225575051661*g1*ZN1(gt1,0)*ZP(gt3,0) + 0.7071067811865475*g2*ZN1(gt1,1)*ZP(gt3,0) + Conj(LamTD)*ZN1(gt1,2)*ZP(gt3,3));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fd, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto ZUR = MODELPARAMETER(ZUR);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ZDR = MODELPARAMETER(ZDR);
   const auto ZUL = MODELPARAMETER(ZUL);

   const std::complex<double> left = SUM(j2,0,2,Conj(ZDL(gt2,j2))*SUM(j1,0,2,Conj(ZUR(gt1,j1))*Yu(j1,j2)))*ZP(gt3,1);

   const std::complex<double> right = SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*ZDR(gt2,j1))*ZUL(gt1,j2))*ZP(gt3,0);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fu, fields::Hpm>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Yd = MODELPARAMETER(Yd);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto ZDR = MODELPARAMETER(ZDR);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ZUR = MODELPARAMETER(ZUR);
   const auto ZDL = MODELPARAMETER(ZDL);

   const std::complex<double> left = SUM(j2,0,2,Conj(ZUL(gt2,j2))*SUM(j1,0,2,Conj(ZDR(gt1,j1))*Yd(j1,j2)))*ZP(gt3,0);

   const std::complex<double> right = SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*ZUR(gt2,j1))*ZDL(gt1,j2))*ZP(gt3,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fv>::type, fields::Fe, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> left = 0;

   const std::complex<double> right = SUM(j1,0,2,Conj(Ye(j1,gt1))*ZER(gt2,j1))*ZP(gt3,0);

   return {left, right};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gWm, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.25*g2*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(vd*ZP(gt3,0) - vu*ZP(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gWm>::type, fields::gZ, fields::Hpm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto vT = MODELPARAMETER(vT);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.25*g2*(vd*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*ZP(gt3,0) + (-(g2*vu*Cos(ThetaW)) + 0.7745966692414834*g1*vu*Sin(ThetaW))*ZP(gt3,1) + 2.8284271247461903*g2*vT*Cos(ThetaW)*(ZP(gt3,2) + ZP(gt3,3)));

   return {result};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gZ, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto vT = MODELPARAMETER(vT);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.25*g2*(vd*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*ZP(gt3,0) + (-(g2*vu*Cos(ThetaW)) + 0.7745966692414834*g1*vu*Sin(ThetaW))*ZP(gt3,1) + 2.8284271247461903*g2*vT*Cos(ThetaW)*(ZP(gt3,2) + ZP(gt3,3)));

   return {result};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gWmC, fields::Hpm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.25*g2*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(vd*ZP(gt3,0) - vu*ZP(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::SRdp, typename fields::conj<fields::Rh>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto vu = MODELPARAMETER(vu);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto vd = MODELPARAMETER(vd);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto vT = MODELPARAMETER(vT);
   const auto vS = MODELPARAMETER(vS);
   const auto MuD = MODELPARAMETER(MuD);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-1.4142135623730951*ZHR(gt3,1)*(2*LamSD*vu*Conj(LamSU)*ZP(gt1,0) + LamTD*Conj(LamTU)*(vu*ZP(gt1,0) + 2*vd*ZP(gt1,1))) - ZHR(gt3,0)*(1.4142135623730951*vd*(-2*AbsSqr(LamSD) + AbsSqr(LamTD) + Sqr(g2))*ZP(gt1,0) + 1.4142135623730951*vu*Sqr(g2)*ZP(gt1,1) + 4*g2*MDWBT*ZP(gt1,2) - 2*vT*AbsSqr(LamTD)*ZP(gt1,2) - 2.8284271247461903*LamTD*vS*Conj(LamSD)*ZP(gt1,2) - 4*LamTD*Conj(MuD)*ZP(gt1,2) + 2*vT*Sqr(g2)*ZP(gt1,2) + 2*vT*AbsSqr(LamTD)*ZP(gt1,3) - 4*MuD*Conj(LamTD)*ZP(gt1,3) - 2.8284271247461903*LamSD*vS*Conj(LamTD)*ZP(gt1,3) + 4*g2*Conj(MDWBT)*ZP(gt1,3) - 2*vT*Sqr(g2)*ZP(gt1,3)));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::Sv, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto vd = MODELPARAMETER(vd);
   const auto Mu = MODELPARAMETER(Mu);
   const auto g2 = MODELPARAMETER(g2);
   const auto vu = MODELPARAMETER(vu);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto vT = MODELPARAMETER(vT);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.7071067811865475*vd*SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZE(gt3,j3))*ZP(gt1,0) + Conj(Mu)*SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)))*ZP(gt1,1) - 0.25*g2*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZE(gt3,j1))*(1.4142135623730951*g2*vd*ZP(gt1,0) + 1.4142135623730951*g2*vu*ZP(gt1,1) + 4*MDWBT*ZP(gt1,2) + 2*g2*vT*ZP(gt1,2) - 2*g2*vT*ZP(gt1,3) + 4*Conj(MDWBT)*ZP(gt1,3));

   return {result};
}

InverseMetricVertex VertexImpl<typename fields::conj<fields::Hpm>::type, fields::VP, fields::VWm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto vT = MODELPARAMETER(vT);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.5*g2*(0.7745966692414834*g1*vd*Cos(ThetaW)*ZP(gt1,0) - 0.7745966692414834*g1*vu*Cos(ThetaW)*ZP(gt1,1) + 1.4142135623730951*g2*vT*Sin(ThetaW)*(ZP(gt1,2) + ZP(gt1,3)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VP>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto vT = MODELPARAMETER(vT);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.5*g2*(0.7745966692414834*g1*vd*Cos(ThetaW)*ZP(gt1,0) - 0.7745966692414834*g1*vu*Cos(ThetaW)*ZP(gt1,1) + 1.4142135623730951*g2*vT*Sin(ThetaW)*(ZP(gt1,2) + ZP(gt1,3)));

   return {result};
}

InverseMetricVertex VertexImpl<typename fields::conj<fields::Hpm>::type, fields::VWm, fields::VZ>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto vT = MODELPARAMETER(vT);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.5*g2*(-0.7745966692414834*g1*vd*Sin(ThetaW)*ZP(gt1,0) + 0.7745966692414834*g1*vu*Sin(ThetaW)*ZP(gt1,1) + 1.4142135623730951*g2*vT*Cos(ThetaW)*(ZP(gt1,2) + ZP(gt1,3)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VZ>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto vT = MODELPARAMETER(vT);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.5*g2*(-0.7745966692414834*g1*vd*Sin(ThetaW)*ZP(gt1,0) + 0.7745966692414834*g1*vu*Sin(ThetaW)*ZP(gt1,1) + 1.4142135623730951*g2*vT*Cos(ThetaW)*(ZP(gt1,2) + ZP(gt1,3)));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::Sd, typename fields::conj<fields::Su>::type, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);

   const std::complex<double> result = -0.7071067811865475*g2*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZU(gt2,j1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::Se, typename fields::conj<fields::Sv>::type, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);

   const std::complex<double> result = -0.7071067811865475*g2*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZV(gt2,j1));

   return {result, minuend_index, subtrahend_index};
}

ChiralVertex VertexImpl<fields::Cha1, fields::Chi, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto ZN1 = MODELPARAMETER(ZN1);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto ZN2 = MODELPARAMETER(ZN2);

   const std::complex<double> left = -0.7071067811865475*Conj(UP1(gt1,1))*(0.7745966692414834*g1*Conj(ZN1(gt2,0)) + g2*Conj(ZN1(gt2,1))) - g2*Conj(UP1(gt1,0))*Conj(ZN1(gt2,2));

   const std::complex<double> right = -(Conj(LamSD)*UM1(gt1,1)*ZN2(gt2,0)) + Conj(LamTD)*(0.7071067811865475*UM1(gt1,1)*ZN2(gt2,1) - UM1(gt1,0)*ZN2(gt2,2));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha1>::type, typename fields::bar<fields::Chi>::type, fields::SRdp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto ZN2 = MODELPARAMETER(ZN2);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto ZN1 = MODELPARAMETER(ZN1);

   const std::complex<double> left = Conj(UM1(gt1,1))*(-(LamSD*Conj(ZN2(gt2,0))) + 0.7071067811865475*LamTD*Conj(ZN2(gt2,1))) - LamTD*Conj(UM1(gt1,0))*Conj(ZN2(gt2,2));

   const std::complex<double> right = -0.7071067811865475*UP1(gt1,1)*(0.7745966692414834*g1*ZN1(gt2,0) + g2*ZN1(gt2,1)) - g2*UP1(gt1,0)*ZN1(gt2,2);

   return {left, right};
}

MomentumDifferenceVertex VertexImpl<fields::SRdp, typename fields::conj<fields::Rh>::type, fields::VWm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.7071067811865475*g2*ZHR(gt2,0);

   return {result, minuend_index, subtrahend_index};
}

ScalarVertex VertexImpl<fields::Sd, fields::SRdp, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const auto MuD = MODELPARAMETER(MuD);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto vS = MODELPARAMETER(vS);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto vT = MODELPARAMETER(vT);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);

   const std::complex<double> result = 0.5*(2*MuD + 1.4142135623730951*LamSD*vS - LamTD*vT)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt1,3 + j1)))*ZU(gt3,j2));

   return {result};
}

ScalarVertex VertexImpl<fields::Se, fields::SRdp, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const auto MuD = MODELPARAMETER(MuD);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto vS = MODELPARAMETER(vS);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto vT = MODELPARAMETER(vT);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);

   const std::complex<double> result = 0.5*(2*MuD + 1.4142135623730951*LamSD*vS - LamTD*vT)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt1,3 + j1)))*ZV(gt3,j2));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Ah>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MuD = MODELPARAMETER(MuD);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto vT = MODELPARAMETER(vT);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto vS = MODELPARAMETER(vS);
   const auto MuU = MODELPARAMETER(MuU);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,-0.25)*(-1.5491933384829668*g1*MDBS*ZA(gt1,0)*ZA(gt2,2)*ZA(gt3,0) - 2.8284271247461903*MuD*Conj(LamSD)*ZA(gt1,0)*ZA(gt2,2)*ZA(gt3,0) - 1.4142135623730951*LamTD*vT*Conj(LamSD)*ZA(gt1,0)*ZA(gt2,2)*ZA(gt3,0) + 1.4142135623730951*LamSD*vT*Conj(LamTD)*ZA(gt1,0)*ZA(gt2,2)*ZA(gt3,0) + 1.5491933384829668*g1*Conj(MDBS)*ZA(gt1,0)*ZA(gt2,2)*ZA(gt3,0) + 2.8284271247461903*LamSD*Conj(MuD)*ZA(gt1,0)*ZA(gt2,2)*ZA(gt3,0) + 2*g2*MDWBT*ZA(gt1,0)*ZA(gt2,3)*ZA(gt3,0) + 1.4142135623730951*LamTD*vS*Conj(LamSD)*ZA(gt1,0)*ZA(gt2,3)*ZA(gt3,0) - 2*MuD*Conj(LamTD)*ZA(gt1,0)*ZA(gt2,3)*ZA(gt3,0) - 1.4142135623730951*LamSD*vS*Conj(LamTD)*ZA(gt1,0)*ZA(gt2,3)*ZA(gt3,0) - 2*g2*Conj(MDWBT)*ZA(gt1,0)*ZA(gt2,3)*ZA(gt3,0) + 2*LamTD*Conj(MuD)*ZA(gt1,0)*ZA(gt2,3)*ZA(gt3,0) + 1.5491933384829668*g1*MDBS*ZA(gt1,1)*ZA(gt2,2)*ZA(gt3,1) - 2.8284271247461903*MuU*Conj(LamSU)*ZA(gt1,1)*ZA(gt2,2)*ZA(gt3,1) + 1.4142135623730951*LamTU*vT*Conj(LamSU)*ZA(gt1,1)*ZA(gt2,2)*ZA(gt3,1) - 1.4142135623730951*LamSU*vT*Conj(LamTU)*ZA(gt1,1)*ZA(gt2,2)*ZA(gt3,1) - 1.5491933384829668*g1*Conj(MDBS)*ZA(gt1,1)*ZA(gt2,2)*ZA(gt3,1) + 2.8284271247461903*LamSU*Conj(MuU)*ZA(gt1,1)*ZA(gt2,2)*ZA(gt3,1) - 2*g2*MDWBT*ZA(gt1,1)*ZA(gt2,3)*ZA(gt3,1) - 1.4142135623730951*LamTU*vS*Conj(LamSU)*ZA(gt1,1)*ZA(gt2,3)*ZA(gt3,1) + 2*MuU*Conj(LamTU)*ZA(gt1,1)*ZA(gt2,3)*ZA(gt3,1) + 1.4142135623730951*LamSU*vS*Conj(LamTU)*ZA(gt1,1)*ZA(gt2,3)*ZA(gt3,1) + 2*g2*Conj(MDWBT)*ZA(gt1,1)*ZA(gt2,3)*ZA(gt3,1) - 2*LamTU*Conj(MuU)*ZA(gt1,1)*ZA(gt2,3)*ZA(gt3,1) + ZA(gt1,2)*((-1.5491933384829668*g1*MDBS - 1.4142135623730951*(2*MuD + LamTD*vT)*Conj(LamSD) + 1.4142135623730951*LamSD*vT*Conj(LamTD) + 1.5491933384829668*g1*Conj(MDBS) + 2.8284271247461903*LamSD*Conj(MuD))*ZA(gt2,0)*ZA(gt3,0) + (1.5491933384829668*g1*MDBS + 1.4142135623730951*(-2*MuU + LamTU*vT)*Conj(LamSU) - 1.4142135623730951*LamSU*vT*Conj(LamTU) - 1.5491933384829668*g1*Conj(MDBS) + 2.8284271247461903*LamSU*Conj(MuU))*ZA(gt2,1)*ZA(gt3,1)) + ZA(gt1,3)*((2*g2*MDWBT + 1.4142135623730951*LamTD*vS*Conj(LamSD) - 2*MuD*Conj(LamTD) - 1.4142135623730951*LamSD*vS*Conj(LamTD) - 2*g2*Conj(MDWBT) + 2*LamTD*Conj(MuD))*ZA(gt2,0)*ZA(gt3,0) + (-2*g2*MDWBT - 1.4142135623730951*LamTU*vS*Conj(LamSU) + 2*MuU*Conj(LamTU) + 1.4142135623730951*LamSU*vS*Conj(LamTU) + 2*g2*Conj(MDWBT) - 2*LamTU*Conj(MuU))*ZA(gt2,1)*ZA(gt3,1)) - 1.5491933384829668*g1*MDBS*ZA(gt1,0)*ZA(gt2,0)*ZA(gt3,2) - 2.8284271247461903*MuD*Conj(LamSD)*ZA(gt1,0)*ZA(gt2,0)*ZA(gt3,2) - 1.4142135623730951*LamTD*vT*Conj(LamSD)*ZA(gt1,0)*ZA(gt2,0)*ZA(gt3,2) + 1.4142135623730951*LamSD*vT*Conj(LamTD)*ZA(gt1,0)*ZA(gt2,0)*ZA(gt3,2) + 1.5491933384829668*g1*Conj(MDBS)*ZA(gt1,0)*ZA(gt2,0)*ZA(gt3,2) + 2.8284271247461903*LamSD*Conj(MuD)*ZA(gt1,0)*ZA(gt2,0)*ZA(gt3,2) + 1.5491933384829668*g1*MDBS*ZA(gt1,1)*ZA(gt2,1)*ZA(gt3,2) - 2.8284271247461903*MuU*Conj(LamSU)*ZA(gt1,1)*ZA(gt2,1)*ZA(gt3,2) + 1.4142135623730951*LamTU*vT*Conj(LamSU)*ZA(gt1,1)*ZA(gt2,1)*ZA(gt3,2) - 1.4142135623730951*LamSU*vT*Conj(LamTU)*ZA(gt1,1)*ZA(gt2,1)*ZA(gt3,2) - 1.5491933384829668*g1*Conj(MDBS)*ZA(gt1,1)*ZA(gt2,1)*ZA(gt3,2) + 2.8284271247461903*LamSU*Conj(MuU)*ZA(gt1,1)*ZA(gt2,1)*ZA(gt3,2) + 2*g2*MDWBT*ZA(gt1,0)*ZA(gt2,0)*ZA(gt3,3) + 1.4142135623730951*LamTD*vS*Conj(LamSD)*ZA(gt1,0)*ZA(gt2,0)*ZA(gt3,3) - 2*MuD*Conj(LamTD)*ZA(gt1,0)*ZA(gt2,0)*ZA(gt3,3) - 1.4142135623730951*LamSD*vS*Conj(LamTD)*ZA(gt1,0)*ZA(gt2,0)*ZA(gt3,3) - 2*g2*Conj(MDWBT)*ZA(gt1,0)*ZA(gt2,0)*ZA(gt3,3) + 2*LamTD*Conj(MuD)*ZA(gt1,0)*ZA(gt2,0)*ZA(gt3,3) - 2*g2*MDWBT*ZA(gt1,1)*ZA(gt2,1)*ZA(gt3,3) - 1.4142135623730951*LamTU*vS*Conj(LamSU)*ZA(gt1,1)*ZA(gt2,1)*ZA(gt3,3) + 2*MuU*Conj(LamTU)*ZA(gt1,1)*ZA(gt2,1)*ZA(gt3,3) + 1.4142135623730951*LamSU*vS*Conj(LamTU)*ZA(gt1,1)*ZA(gt2,1)*ZA(gt3,3) + 2*g2*Conj(MDWBT)*ZA(gt1,1)*ZA(gt2,1)*ZA(gt3,3) - 2*LamTU*Conj(MuU)*ZA(gt1,1)*ZA(gt2,1)*ZA(gt3,3));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, typename fields::conj<fields::Rh>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Mu = MODELPARAMETER(Mu);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.25*Mu*(2*Conj(LamSD)*(ZA(gt1,2)*ZA(gt2,1) + ZA(gt1,1)*ZA(gt2,2))*ZHR(gt3,0) + 1.4142135623730951*Conj(LamTD)*(ZA(gt1,3)*ZA(gt2,1) + ZA(gt1,1)*ZA(gt2,3))*ZHR(gt3,0) + (-2*Conj(LamSU)*(ZA(gt1,2)*ZA(gt2,0) + ZA(gt1,0)*ZA(gt2,2)) + 1.4142135623730951*Conj(LamTU)*(ZA(gt1,3)*ZA(gt2,0) + ZA(gt1,0)*ZA(gt2,3)))*ZHR(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Rh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Mu = MODELPARAMETER(Mu);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.25*Conj(Mu)*(2*LamSD*ZA(gt1,1)*ZA(gt2,2)*ZHR(gt3,0) + 1.4142135623730951*LamTD*ZA(gt1,1)*ZA(gt2,3)*ZHR(gt3,0) - 2*LamSU*ZA(gt1,0)*ZA(gt2,2)*ZHR(gt3,1) + 1.4142135623730951*LamTU*ZA(gt1,0)*ZA(gt2,3)*ZHR(gt3,1) + 2*ZA(gt1,2)*(LamSD*ZA(gt2,1)*ZHR(gt3,0) - LamSU*ZA(gt2,0)*ZHR(gt3,1)) + 1.4142135623730951*ZA(gt1,3)*(LamTD*ZA(gt2,1)*ZHR(gt3,0) + LamTU*ZA(gt2,0)*ZHR(gt3,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Rh, typename fields::conj<fields::Rh>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto vu = MODELPARAMETER(vu);
   const auto vd = MODELPARAMETER(vd);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto g1 = MODELPARAMETER(g1);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MuD = MODELPARAMETER(MuD);
   const auto vT = MODELPARAMETER(vT);
   const auto MuU = MODELPARAMETER(MuU);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto vS = MODELPARAMETER(vS);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = std::complex<double>(0,-0.25)*((vu*ZA(gt1,0) - vd*ZA(gt1,1))*(2*LamSU*Conj(LamSD)*ZHR(gt2,1)*ZHR(gt3,0) - LamTU*Conj(LamTD)*ZHR(gt2,1)*ZHR(gt3,0) + (-2*LamSD*Conj(LamSU) + LamTD*Conj(LamTU))*ZHR(gt2,0)*ZHR(gt3,1)) + ZA(gt1,2)*((1.5491933384829668*g1*MDBS - 1.4142135623730951*(2*MuD + LamTD*vT)*Conj(LamSD) + 1.4142135623730951*LamSD*vT*Conj(LamTD) - 1.5491933384829668*g1*Conj(MDBS) + 2.8284271247461903*LamSD*Conj(MuD))*ZHR(gt2,0)*ZHR(gt3,0) + (-1.5491933384829668*g1*MDBS + 1.4142135623730951*(-2*MuU + LamTU*vT)*Conj(LamSU) - 1.4142135623730951*LamSU*vT*Conj(LamTU) + 1.5491933384829668*g1*Conj(MDBS) + 2.8284271247461903*LamSU*Conj(MuU))*ZHR(gt2,1)*ZHR(gt3,1)) + ZA(gt1,3)*((-2*g2*MDWBT + 1.4142135623730951*LamTD*vS*Conj(LamSD) - (2*MuD + 1.4142135623730951*LamSD*vS)*Conj(LamTD) + 2*g2*Conj(MDWBT) + 2*LamTD*Conj(MuD))*ZHR(gt2,0)*ZHR(gt3,0) + (2*g2*MDWBT - 1.4142135623730951*LamTU*vS*Conj(LamSU) + (2*MuU + 1.4142135623730951*LamSU*vS)*Conj(LamTU) - 2*g2*Conj(MDWBT) - 2*LamTU*Conj(MuU))*ZHR(gt2,1)*ZHR(gt3,1)));

   return {result};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha1>::type, fields::Cha1, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = 0.5*(-(Conj(UM1(gt1,0))*(1.4142135623730951*LamTD*Conj(UP1(gt2,1))*ZH(gt3,0) + 2*g2*Conj(UP1(gt2,0))*ZH(gt3,3))) - Conj(UM1(gt1,1))*(1.4142135623730951*g2*Conj(UP1(gt2,0))*ZH(gt3,0) + Conj(UP1(gt2,1))*(1.4142135623730951*LamSD*ZH(gt3,2) - LamTD*ZH(gt3,3))));

   const std::complex<double> right = 0.5*(-(UM1(gt2,0)*(1.4142135623730951*Conj(LamTD)*UP1(gt1,1)*ZH(gt3,0) + 2*g2*UP1(gt1,0)*ZH(gt3,3))) - UM1(gt2,1)*(1.4142135623730951*g2*UP1(gt1,0)*ZH(gt3,0) + UP1(gt1,1)*(1.4142135623730951*Conj(LamSD)*ZH(gt3,2) - Conj(LamTD)*ZH(gt3,3))));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha1>::type, fields::Cha1, fields::Ah>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0,0.5)*(Conj(UM1(gt1,0))*(-1.4142135623730951*LamTD*Conj(UP1(gt2,1))*ZA(gt3,0) + 2*g2*Conj(UP1(gt2,0))*ZA(gt3,3)) + Conj(UM1(gt1,1))*(1.4142135623730951*g2*Conj(UP1(gt2,0))*ZA(gt3,0) + Conj(UP1(gt2,1))*(-1.4142135623730951*LamSD*ZA(gt3,2) + LamTD*ZA(gt3,3))));

   const std::complex<double> right = std::complex<double>(0,-0.5)*(UM1(gt2,0)*(-1.4142135623730951*Conj(LamTD)*UP1(gt1,1)*ZA(gt3,0) + 2*g2*UP1(gt1,0)*ZA(gt3,3)) + UM1(gt2,1)*(1.4142135623730951*g2*UP1(gt1,0)*ZA(gt3,0) + UP1(gt1,1)*(-1.4142135623730951*Conj(LamSD)*ZA(gt3,2) + Conj(LamTD)*ZA(gt3,3))));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha2>::type, fields::Cha2, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto UM2 = MODELPARAMETER(UM2);
   const auto UP2 = MODELPARAMETER(UP2);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = 0.5*(g2*Conj(UM2(gt2,0))*(-1.4142135623730951*Conj(UP2(gt1,1))*ZH(gt3,1) + 2*Conj(UP2(gt1,0))*ZH(gt3,3)) + Conj(UM2(gt2,1))*(1.4142135623730951*LamTU*Conj(UP2(gt1,0))*ZH(gt3,1) + Conj(UP2(gt1,1))*(1.4142135623730951*LamSU*ZH(gt3,2) + LamTU*ZH(gt3,3))));

   const std::complex<double> right = 0.5*(1.4142135623730951*Conj(LamSU)*UM2(gt1,1)*UP2(gt2,1)*ZH(gt3,2) + UM2(gt1,0)*(-1.4142135623730951*g2*UP2(gt2,1)*ZH(gt3,1) + 2*g2*UP2(gt2,0)*ZH(gt3,3)) + Conj(LamTU)*UM2(gt1,1)*(1.4142135623730951*UP2(gt2,0)*ZH(gt3,1) + UP2(gt2,1)*ZH(gt3,3)));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha2>::type, fields::Cha2, fields::Ah>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto UM2 = MODELPARAMETER(UM2);
   const auto UP2 = MODELPARAMETER(UP2);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0,0.5)*(g2*Conj(UM2(gt2,0))*(1.4142135623730951*Conj(UP2(gt1,1))*ZA(gt3,1) - 2*Conj(UP2(gt1,0))*ZA(gt3,3)) + Conj(UM2(gt2,1))*(1.4142135623730951*LamTU*Conj(UP2(gt1,0))*ZA(gt3,1) + Conj(UP2(gt1,1))*(1.4142135623730951*LamSU*ZA(gt3,2) + LamTU*ZA(gt3,3))));

   const std::complex<double> right = std::complex<double>(0,-0.5)*(1.4142135623730951*Conj(LamSU)*UM2(gt1,1)*UP2(gt2,1)*ZA(gt3,2) + g2*UM2(gt1,0)*(1.4142135623730951*UP2(gt2,1)*ZA(gt3,1) - 2*UP2(gt2,0)*ZA(gt3,3)) + Conj(LamTU)*UM2(gt1,1)*(1.4142135623730951*UP2(gt2,0)*ZA(gt3,1) + UP2(gt2,1)*ZA(gt3,3)));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Chi>::type, fields::Chi, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto g1 = MODELPARAMETER(g1);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto ZN1 = MODELPARAMETER(ZN1);
   const auto ZN2 = MODELPARAMETER(ZN2);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = 0.1*(3.872983346207417*g1*Conj(ZN1(gt2,0))*(Conj(ZN2(gt1,2))*ZH(gt3,0) - Conj(ZN2(gt1,3))*ZH(gt3,1)) + 5*Conj(ZN1(gt2,1))*(-(g2*Conj(ZN2(gt1,2))*ZH(gt3,0)) + g2*Conj(ZN2(gt1,3))*ZH(gt3,1)) + 5*Conj(ZN1(gt2,2))*(1.4142135623730951*LamSD*Conj(ZN2(gt1,0))*ZH(gt3,0) + LamTD*Conj(ZN2(gt1,1))*ZH(gt3,0) + Conj(ZN2(gt1,2))*(1.4142135623730951*LamSD*ZH(gt3,2) + LamTD*ZH(gt3,3))) + 5*Conj(ZN1(gt2,3))*(-1.4142135623730951*LamSU*Conj(ZN2(gt1,0))*ZH(gt3,1) + LamTU*Conj(ZN2(gt1,1))*ZH(gt3,1) + Conj(ZN2(gt1,3))*(-1.4142135623730951*LamSU*ZH(gt3,2) + LamTU*ZH(gt3,3))));

   const std::complex<double> right = 0.5*(Conj(LamTD)*ZH(gt3,0)*ZN1(gt1,2)*ZN2(gt2,1) + Conj(LamTU)*ZH(gt3,1)*ZN1(gt1,3)*ZN2(gt2,1) + 0.7745966692414834*g1*ZH(gt3,0)*ZN1(gt1,0)*ZN2(gt2,2) - g2*ZH(gt3,0)*ZN1(gt1,1)*ZN2(gt2,2) + Conj(LamTD)*ZH(gt3,3)*ZN1(gt1,2)*ZN2(gt2,2) + 1.4142135623730951*Conj(LamSD)*ZN1(gt1,2)*(ZH(gt3,0)*ZN2(gt2,0) + ZH(gt3,2)*ZN2(gt2,2)) - 0.7745966692414834*g1*ZH(gt3,1)*ZN1(gt1,0)*ZN2(gt2,3) + g2*ZH(gt3,1)*ZN1(gt1,1)*ZN2(gt2,3) + Conj(LamTU)*ZH(gt3,3)*ZN1(gt1,3)*ZN2(gt2,3) - 1.4142135623730951*Conj(LamSU)*ZN1(gt1,3)*(ZH(gt3,1)*ZN2(gt2,0) + ZH(gt3,2)*ZN2(gt2,3)));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Chi>::type, fields::Chi, fields::Ah>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto g1 = MODELPARAMETER(g1);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto ZN1 = MODELPARAMETER(ZN1);
   const auto ZN2 = MODELPARAMETER(ZN2);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0,0.1)*(5*g2*Conj(ZN1(gt2,1))*(Conj(ZN2(gt1,2))*ZA(gt3,0) - Conj(ZN2(gt1,3))*ZA(gt3,1)) + 3.872983346207417*g1*Conj(ZN1(gt2,0))*(-(Conj(ZN2(gt1,2))*ZA(gt3,0)) + Conj(ZN2(gt1,3))*ZA(gt3,1)) + 5*Conj(ZN1(gt2,2))*(1.4142135623730951*LamSD*Conj(ZN2(gt1,0))*ZA(gt3,0) + LamTD*Conj(ZN2(gt1,1))*ZA(gt3,0) + Conj(ZN2(gt1,2))*(1.4142135623730951*LamSD*ZA(gt3,2) + LamTD*ZA(gt3,3))) + 5*Conj(ZN1(gt2,3))*(-1.4142135623730951*LamSU*Conj(ZN2(gt1,0))*ZA(gt3,1) + LamTU*Conj(ZN2(gt1,1))*ZA(gt3,1) + Conj(ZN2(gt1,3))*(-1.4142135623730951*LamSU*ZA(gt3,2) + LamTU*ZA(gt3,3))));

   const std::complex<double> right = std::complex<double>(0,-0.5)*(Conj(LamTD)*ZA(gt3,0)*ZN1(gt1,2)*ZN2(gt2,1) + Conj(LamTU)*ZA(gt3,1)*ZN1(gt1,3)*ZN2(gt2,1) - 0.7745966692414834*g1*ZA(gt3,0)*ZN1(gt1,0)*ZN2(gt2,2) + g2*ZA(gt3,0)*ZN1(gt1,1)*ZN2(gt2,2) + Conj(LamTD)*ZA(gt3,3)*ZN1(gt1,2)*ZN2(gt2,2) + 1.4142135623730951*Conj(LamSD)*ZN1(gt1,2)*(ZA(gt3,0)*ZN2(gt2,0) + ZA(gt3,2)*ZN2(gt2,2)) + 0.7745966692414834*g1*ZA(gt3,1)*ZN1(gt1,0)*ZN2(gt2,3) - g2*ZA(gt3,1)*ZN1(gt1,1)*ZN2(gt2,3) + Conj(LamTU)*ZA(gt3,3)*ZN1(gt1,3)*ZN2(gt2,3) - 1.4142135623730951*Conj(LamSU)*ZN1(gt1,3)*(ZA(gt3,1)*ZN2(gt2,0) + ZA(gt3,2)*ZN2(gt2,3)));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fd, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto ZDR = MODELPARAMETER(ZDR);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(ZDL(gt2,j2))*SUM(j1,0,2,Conj(ZDR(gt1,j1))*Yd(j1,j2)))*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*ZDR(gt2,j1))*ZDL(gt1,j2))*ZH(gt3,0);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fd, fields::Ah>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto ZDR = MODELPARAMETER(ZDR);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,Conj(ZDL(gt2,j2))*SUM(j1,0,2,Conj(ZDR(gt1,j1))*Yd(j1,j2)))*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*ZDR(gt2,j1))*ZDL(gt1,j2))*ZA(gt3,0);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fu, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto ZUR = MODELPARAMETER(ZUR);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(ZUL(gt2,j2))*SUM(j1,0,2,Conj(ZUR(gt1,j1))*Yu(j1,j2)))*ZH(gt3,1);

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*ZUR(gt2,j1))*ZUL(gt1,j2))*ZH(gt3,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fu, fields::Ah>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto ZUR = MODELPARAMETER(ZUR);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,Conj(ZUL(gt2,j2))*SUM(j1,0,2,Conj(ZUR(gt1,j1))*Yu(j1,j2)))*ZA(gt3,1);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*ZUR(gt2,j1))*ZUL(gt1,j2))*ZA(gt3,1);

   return {left, right};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gWm>::type, fields::gWm, fields::hh>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto vT = MODELPARAMETER(vT);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = -0.25*Sqr(g2)*(vd*ZH(gt3,0) + vu*ZH(gt3,1) + 4*vT*ZH(gt3,3));

   return {result};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gWm>::type, fields::gWm, fields::Ah>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,-0.25)*Sqr(g2)*(vd*ZA(gt3,0) - vu*ZA(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gWmC, fields::hh>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto vT = MODELPARAMETER(vT);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = -0.25*Sqr(g2)*(vd*ZH(gt3,0) + vu*ZH(gt3,1) + 4*vT*ZH(gt3,3));

   return {result};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gWmC, fields::Ah>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,0.25)*Sqr(g2)*(vd*ZA(gt3,0) - vu*ZA(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gZ, fields::hh>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.25*Sqr(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(vd*ZH(gt3,0) + vu*ZH(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Sd, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Mu = MODELPARAMETER(Mu);
   const auto g1 = MODELPARAMETER(g1);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,-0.16666666666666666)*(4.242640687119286*Conj(Mu)*SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1)))*ZA(gt1,1) - 4.242640687119286*Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1)))*ZD(gt3,j2))*ZA(gt1,1) + 0.7745966692414834*g1*MDBS*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*ZA(gt1,2) - 0.7745966692414834*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*ZA(gt1,2) + 1.5491933384829668*g1*MDBS*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*ZA(gt1,2) - 1.5491933384829668*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*ZA(gt1,2) - 3*g2*MDWBT*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*ZA(gt1,3) + 3*g2*Conj(MDWBT)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*ZA(gt1,3));

   return {result};
}

ScalarVertex VertexImpl<fields::Sd, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto vS = MODELPARAMETER(vS);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto vT = MODELPARAMETER(vT);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto MuD = MODELPARAMETER(MuD);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = 0.5*(1.4142135623730951*vS*Conj(LamSD) + vT*Conj(LamTD) + 2*Conj(MuD))*SUM(j2,0,2,Conj(ZD(gt1,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1)))*ZHR(gt2,0);

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Se, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Mu = MODELPARAMETER(Mu);
   const auto g1 = MODELPARAMETER(g1);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,-0.1)*(7.0710678118654755*Conj(Mu)*SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)))*ZA(gt1,1) - 7.0710678118654755*Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1)))*ZE(gt3,j2))*ZA(gt1,1) - 3.872983346207417*g1*MDBS*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*ZA(gt1,2) + 3.872983346207417*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*ZA(gt1,2) + 7.745966692414834*g1*MDBS*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*ZA(gt1,2) - 7.745966692414834*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*ZA(gt1,2) - 5*g2*MDWBT*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*ZA(gt1,3) + 5*g2*Conj(MDWBT)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*ZA(gt1,3));

   return {result};
}

ScalarVertex VertexImpl<fields::Se, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto vS = MODELPARAMETER(vS);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto vT = MODELPARAMETER(vT);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto MuD = MODELPARAMETER(MuD);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = 0.5*(1.4142135623730951*vS*Conj(LamSD) + vT*Conj(LamTD) + 2*Conj(MuD))*SUM(j2,0,2,Conj(ZE(gt1,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)))*ZHR(gt2,0);

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Su, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Mu = MODELPARAMETER(Mu);
   const auto g1 = MODELPARAMETER(g1);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,-0.16666666666666666)*(4.242640687119286*Conj(Mu)*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt3,3 + j1)))*ZA(gt1,0) - 4.242640687119286*Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1)))*ZU(gt3,j2))*ZA(gt1,0) + 0.7745966692414834*g1*MDBS*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*ZA(gt1,2) - 0.7745966692414834*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*ZA(gt1,2) - 3.0983866769659336*g1*MDBS*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*ZA(gt1,2) + 3.0983866769659336*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*ZA(gt1,2) + 3*g2*MDWBT*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*ZA(gt1,3) - 3*g2*Conj(MDWBT)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*ZA(gt1,3));

   return {result};
}

ScalarVertex VertexImpl<fields::Su, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto vS = MODELPARAMETER(vS);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto vT = MODELPARAMETER(vT);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto MuU = MODELPARAMETER(MuU);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.5*(1.4142135623730951*vS*Conj(LamSU) - vT*Conj(LamTU) + 2*Conj(MuU))*SUM(j2,0,2,Conj(ZU(gt1,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt3,3 + j1)))*ZHR(gt2,1);

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Sv, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,-0.5)*KroneckerDelta(gt2,gt3)*(0.7745966692414834*g1*(-MDBS + Conj(MDBS))*ZA(gt1,2) + g2*(MDWBT - Conj(MDWBT))*ZA(gt1,3));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*Sqr(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(vd*ZH(gt1,0) + vu*ZH(gt1,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::Sd, fields::SRdp, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt4 = indices[2];
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -(LamTD*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1)))*ZD(gt4,j2))*ZP(gt1,2));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::Se, fields::SRdp, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt4 = indices[2];
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -(LamTD*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1)))*ZE(gt4,j2))*ZP(gt1,2));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Hpm, fields::SRdp, typename fields::conj<fields::Rh>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt4 = indices[2];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-2.8284271247461903*LamSD*Conj(LamSU)*ZH(gt1,1)*ZHR(gt4,1)*ZP(gt2,0) - 1.4142135623730951*LamTD*Conj(LamTU)*ZH(gt1,1)*ZHR(gt4,1)*ZP(gt2,0) - 1.4142135623730951*Sqr(g2)*ZH(gt1,1)*ZHR(gt4,0)*ZP(gt2,1) - 1.4142135623730951*ZH(gt1,0)*((-2*AbsSqr(LamSD) + AbsSqr(LamTD) + Sqr(g2))*ZHR(gt4,0)*ZP(gt2,0) + 2*LamTD*Conj(LamTU)*ZHR(gt4,1)*ZP(gt2,1)) + 2.8284271247461903*LamTD*Conj(LamSD)*ZH(gt1,2)*ZHR(gt4,0)*ZP(gt2,2) + 2*AbsSqr(LamTD)*ZH(gt1,3)*ZHR(gt4,0)*ZP(gt2,2) - 2*Sqr(g2)*ZH(gt1,3)*ZHR(gt4,0)*ZP(gt2,2) + 2.8284271247461903*LamSD*Conj(LamTD)*ZH(gt1,2)*ZHR(gt4,0)*ZP(gt2,3) - 2*AbsSqr(LamTD)*ZH(gt1,3)*ZHR(gt4,0)*ZP(gt2,3) + 2*Sqr(g2)*ZH(gt1,3)*ZHR(gt4,0)*ZP(gt2,3));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Sd, fields::SRdp, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt4 = indices[2];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.5*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1)))*ZU(gt4,j2))*(1.4142135623730951*LamSD*ZH(gt1,2) - LamTD*ZH(gt1,3));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Se, fields::SRdp, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt4 = indices[2];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.5*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1)))*ZV(gt4,j2))*(1.4142135623730951*LamSD*ZH(gt1,2) - LamTD*ZH(gt1,3));

   return {result};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha1>::type, fields::Cha1, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -(g2*Conj(UP1(gt2,0))*Sin(ThetaW)*UP1(gt1,0)) - 0.5*Conj(UP1(gt2,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UP1(gt1,1);

   const std::complex<double> right = -(g2*Conj(UM1(gt1,0))*Sin(ThetaW)*UM1(gt2,0)) - 0.5*Conj(UM1(gt1,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UM1(gt2,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha1>::type, fields::Cha1, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -(g2*Conj(UP1(gt2,0))*Cos(ThetaW)*UP1(gt1,0)) + 0.1*Conj(UP1(gt2,1))*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*UP1(gt1,1);

   const std::complex<double> right = -(g2*Conj(UM1(gt1,0))*Cos(ThetaW)*UM1(gt2,0)) + 0.1*Conj(UM1(gt1,1))*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*UM1(gt2,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha2>::type, fields::Cha2, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UM2 = MODELPARAMETER(UM2);
   const auto UP2 = MODELPARAMETER(UP2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = g2*Conj(UM2(gt2,0))*Sin(ThetaW)*UM2(gt1,0) + 0.5*Conj(UM2(gt2,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UM2(gt1,1);

   const std::complex<double> right = g2*Conj(UP2(gt1,0))*Sin(ThetaW)*UP2(gt2,0) + 0.5*Conj(UP2(gt1,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UP2(gt2,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha2>::type, fields::Cha2, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UM2 = MODELPARAMETER(UM2);
   const auto UP2 = MODELPARAMETER(UP2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = g2*Conj(UM2(gt2,0))*Cos(ThetaW)*UM2(gt1,0) + 0.5*Conj(UM2(gt2,1))*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*UM2(gt1,1);

   const std::complex<double> right = g2*Conj(UP2(gt1,0))*Cos(ThetaW)*UP2(gt2,0) + 0.5*Conj(UP2(gt1,1))*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*UP2(gt2,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Chi>::type, fields::Chi, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZN1 = MODELPARAMETER(ZN1);
   const auto ZN2 = MODELPARAMETER(ZN2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.5*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(Conj(ZN1(gt2,2))*ZN1(gt1,2) - Conj(ZN1(gt2,3))*ZN1(gt1,3));

   const std::complex<double> right = 0.5*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(Conj(ZN2(gt1,2))*ZN2(gt2,2) - Conj(ZN2(gt1,3))*ZN2(gt2,3));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fd, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.16666666666666666*KroneckerDelta(gt1,gt2)*(0.7745966692414834*g1*Cos(ThetaW) - 3*g2*Sin(ThetaW));

   const std::complex<double> right = 0.2581988897471611*g1*Cos(ThetaW)*KroneckerDelta(gt1,gt2);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fd, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.16666666666666666*KroneckerDelta(gt1,gt2)*(3*g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW));

   const std::complex<double> right = -0.2581988897471611*g1*KroneckerDelta(gt1,gt2)*Sin(ThetaW);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fe, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.5*KroneckerDelta(gt1,gt2)*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW));

   const std::complex<double> right = -0.7745966692414834*g1*KroneckerDelta(gt1,gt2)*Sin(ThetaW);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fu, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.16666666666666666*KroneckerDelta(gt1,gt2)*(0.7745966692414834*g1*Cos(ThetaW) + 3*g2*Sin(ThetaW));

   const std::complex<double> right = -0.5163977794943222*g1*Cos(ThetaW)*KroneckerDelta(gt1,gt2);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fu, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.03333333333333333*KroneckerDelta(gt1,gt2)*(-15*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW));

   const std::complex<double> right = 0.5163977794943222*g1*KroneckerDelta(gt1,gt2)*Sin(ThetaW);

   return {left, right};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWm>::type, fields::gWm, fields::VP>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Sin(ThetaW));

   return {result, 1};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWm>::type, fields::gWm, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Cos(ThetaW));

   return {result, 1};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gWmC, fields::VP>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Sin(ThetaW);

   return {result, 1};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gWmC, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Cos(ThetaW);

   return {result, 1};
}

MomentumDifferenceVertex VertexImpl<fields::Rh, typename fields::conj<fields::Rh>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(ZHR(gt1,0)*ZHR(gt2,0) - ZHR(gt1,1)*ZHR(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.03333333333333333*(3.872983346207417*g1*Cos(ThetaW) - 15*g2*Sin(ThetaW))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1)) + 0.2581988897471611*g1*Cos(ThetaW)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.16666666666666666*(3*g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1)) - 0.2581988897471611*g1*Sin(ThetaW)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::Se, typename fields::conj<fields::Se>::type, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) + 0.7745966692414834*g1*Cos(ThetaW)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::Se, typename fields::conj<fields::Se>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) - 0.7745966692414834*g1*Sin(ThetaW)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::SRum, typename fields::conj<fields::SRum>::type, fields::VP>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::SRum, typename fields::conj<fields::SRum>::type, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.16666666666666666*(0.7745966692414834*g1*Cos(ThetaW) + 3*g2*Sin(ThetaW))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1)) - 0.5163977794943222*g1*Cos(ThetaW)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.03333333333333333*((-15*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1)) + 15.491933384829668*g1*Sin(ThetaW)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::Sv, typename fields::conj<fields::Sv>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.5*KroneckerDelta(gt1,gt2)*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW));

   return {result, minuend_index, subtrahend_index};
}

TripleVectorVertex VertexImpl<typename fields::conj<fields::VWm>::type, fields::VP, fields::VWm>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Sin(ThetaW);

   return {result, TripleVectorVertex::odd_permutation{}};
}

TripleVectorVertex VertexImpl<typename fields::conj<fields::VWm>::type, fields::VWm, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Cos(ThetaW));

   return {result, TripleVectorVertex::odd_permutation{}};
}

ScalarVertex VertexImpl<fields::SRum, fields::Su, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt3 = indices[1];
   const auto MuU = MODELPARAMETER(MuU);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto vS = MODELPARAMETER(vS);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto vT = MODELPARAMETER(vT);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = -0.5*(2*MuU + 1.4142135623730951*LamSU*vS + LamTU*vT)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1)))*ZD(gt3,j2));

   return {result};
}

ChiralVertex VertexImpl<fields::Cha2, fields::Chi, typename fields::conj<fields::SRum>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto UM2 = MODELPARAMETER(UM2);
   const auto ZN1 = MODELPARAMETER(ZN1);
   const auto UP2 = MODELPARAMETER(UP2);
   const auto ZN2 = MODELPARAMETER(ZN2);

   const std::complex<double> left = 0.7071067811865475*Conj(UM2(gt1,1))*(0.7745966692414834*g1*Conj(ZN1(gt2,0)) + g2*Conj(ZN1(gt2,1))) - g2*Conj(UM2(gt1,0))*Conj(ZN1(gt2,3));

   const std::complex<double> right = Conj(LamSU)*UP2(gt1,1)*ZN2(gt2,0) + Conj(LamTU)*(0.7071067811865475*UP2(gt1,1)*ZN2(gt2,1) + UP2(gt1,0)*ZN2(gt2,3));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha2>::type, typename fields::bar<fields::Chi>::type, fields::SRum>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto UP2 = MODELPARAMETER(UP2);
   const auto ZN2 = MODELPARAMETER(ZN2);
   const auto UM2 = MODELPARAMETER(UM2);
   const auto ZN1 = MODELPARAMETER(ZN1);

   const std::complex<double> left = Conj(UP2(gt1,1))*(LamSU*Conj(ZN2(gt2,0)) + 0.7071067811865475*LamTU*Conj(ZN2(gt2,1))) + LamTU*Conj(UP2(gt1,0))*Conj(ZN2(gt2,3));

   const std::complex<double> right = UM2(gt1,1)*(0.5477225575051661*g1*ZN1(gt2,0) + 0.7071067811865475*g2*ZN1(gt2,1)) - g2*UM2(gt1,0)*ZN1(gt2,3);

   return {left, right};
}

MomentumDifferenceVertex VertexImpl<fields::SRum, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.7071067811865475*g2*ZHR(gt2,1);

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::Rh, typename fields::conj<fields::SRum>::type, fields::VWm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.7071067811865475*g2*ZHR(gt1,1);

   return {result, minuend_index, subtrahend_index};
}

ScalarVertex VertexImpl<fields::SRum, fields::Su, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt3 = indices[1];
   const int gt4 = indices[2];
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = LamTU*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1)))*ZU(gt4,j2))*ZP(gt3,3);

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::SRum, typename fields::conj<fields::SRum>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = 0.05*((-3*Sqr(g1) + 5*Sqr(g2))*ZA(gt1,0)*ZA(gt2,0) + (-20*AbsSqr(LamTU) + 3*Sqr(g1) - 5*Sqr(g2))*ZA(gt1,1)*ZA(gt2,1) - 5*(Conj(LamSU)*(1.4142135623730951*LamTU*ZA(gt1,3)*ZA(gt2,2) + ZA(gt1,2)*(4*LamSU*ZA(gt2,2) + 1.4142135623730951*LamTU*ZA(gt2,3))) + Conj(LamTU)*(1.4142135623730951*LamSU*ZA(gt1,2)*ZA(gt2,3) + ZA(gt1,3)*(1.4142135623730951*LamSU*ZA(gt2,2) + 2*LamTU*ZA(gt2,3)))));

   return {result};
}

ScalarVertex VertexImpl<fields::Rh, fields::SRum, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::SRum>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = 0.05*((3*Sqr(g1) - 5*Sqr(g2))*ZHR(gt1,0)*ZHR(gt3,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZHR(gt1,1)*ZHR(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Sd, fields::SRum, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::SRum>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = 0.05*((Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1)) + 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Se, fields::SRum, typename fields::conj<fields::Se>::type, typename fields::conj<fields::SRum>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZE = MODELPARAMETER(ZE);

   const std::complex<double> result = 0.05*(-((3*Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1)));

   return {result};
}

ScalarVertex VertexImpl<fields::SRum, fields::SRum, typename fields::conj<fields::SRum>::type, typename fields::conj<fields::SRum>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> result = 0.1*(-3*Sqr(g1) - 5*Sqr(g2));

   return {result};
}

ScalarVertex VertexImpl<fields::SRum, fields::Su, typename fields::conj<fields::SRum>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt4 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZU = MODELPARAMETER(ZU);

   const std::complex<double> result = 0.05*((Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt4,j1)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt4,3 + j1)));

   return {result};
}

ScalarVertex VertexImpl<fields::SRum, fields::Sv, typename fields::conj<fields::SRum>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt4 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> result = -0.05*KroneckerDelta(gt2,gt4)*(3*Sqr(g1) - 5*Sqr(g2));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SRum, typename fields::conj<fields::SRum>::type, fields::VP, fields::VP>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*Sqr(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SRum, typename fields::conj<fields::SRum>::type, typename fields::conj<fields::VWm>::type, fields::VWm>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> result = 0.5*Sqr(g2);

   return {result};
}

InverseMetricVertex VertexImpl<fields::SRum, typename fields::conj<fields::SRum>::type, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*Sqr(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::SRum, fields::Su, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt4 = indices[2];
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = -0.5*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt3,3 + j1)))*ZD(gt4,j2))*(1.4142135623730951*LamSU*ZH(gt1,2) + LamTU*ZH(gt1,3));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::SRum, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Rh>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt4 = indices[2];
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-1.4142135623730951*LamTU*Conj(LamTD)*ZHR(gt4,0)*(2*ZH(gt1,1)*ZP(gt3,0) + ZH(gt1,0)*ZP(gt3,1)) - 1.4142135623730951*ZH(gt1,0)*(Sqr(g2)*ZHR(gt4,1)*ZP(gt3,0) + 2*LamSU*Conj(LamSD)*ZHR(gt4,0)*ZP(gt3,1)) + ZHR(gt4,1)*(-1.4142135623730951*(-2*AbsSqr(LamSU) + AbsSqr(LamTU) + Sqr(g2))*ZH(gt1,1)*ZP(gt3,1) + 2*(Conj(LamTU)*(1.4142135623730951*LamSU*ZH(gt1,2)*ZP(gt3,2) + LamTU*ZH(gt1,3)*(ZP(gt3,2) - ZP(gt3,3))) + 1.4142135623730951*LamTU*Conj(LamSU)*ZH(gt1,2)*ZP(gt3,3) + Sqr(g2)*ZH(gt1,3)*(-ZP(gt3,2) + ZP(gt3,3)))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Hpm, fields::Rh, typename fields::conj<fields::SRum>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-1.4142135623730951*ZH(gt1,0)*(Sqr(g2)*ZHR(gt3,1)*ZP(gt2,0) + 2*LamSD*Conj(LamSU)*ZHR(gt3,0)*ZP(gt2,1)) + ZHR(gt3,1)*(-1.4142135623730951*(-2*AbsSqr(LamSU) + Sqr(g2))*ZH(gt1,1)*ZP(gt2,1) + 2.8284271247461903*LamTU*Conj(LamSU)*ZH(gt1,2)*ZP(gt2,2) + 2*Sqr(g2)*ZH(gt1,3)*(-ZP(gt2,2) + ZP(gt2,3))) - Conj(LamTU)*(1.4142135623730951*LamTD*ZH(gt1,0)*ZHR(gt3,0)*ZP(gt2,1) + 1.4142135623730951*ZH(gt1,1)*(2*LamTD*ZHR(gt3,0)*ZP(gt2,0) + LamTU*ZHR(gt3,1)*ZP(gt2,1)) - 2*ZHR(gt3,1)*(LamTU*ZH(gt1,3)*(ZP(gt2,2) - ZP(gt2,3)) + 1.4142135623730951*LamSU*ZH(gt1,2)*ZP(gt2,3))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Hpm, fields::Rh, typename fields::conj<fields::SRum>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,-0.25)*(1.4142135623730951*ZA(gt1,0)*(-(Sqr(g2)*ZHR(gt3,1)*ZP(gt2,0)) + 2*LamSD*Conj(LamSU)*ZHR(gt3,0)*ZP(gt2,1)) + ZHR(gt3,1)*(1.4142135623730951*(-2*AbsSqr(LamSU) + Sqr(g2))*ZA(gt1,1)*ZP(gt2,1) + 2.8284271247461903*LamTU*Conj(LamSU)*ZA(gt1,2)*ZP(gt2,2) - 2*Sqr(g2)*ZA(gt1,3)*(ZP(gt2,2) + ZP(gt2,3))) + Conj(LamTU)*(1.4142135623730951*LamTD*ZA(gt1,0)*ZHR(gt3,0)*ZP(gt2,1) + 1.4142135623730951*ZA(gt1,1)*(-2*LamTD*ZHR(gt3,0)*ZP(gt2,0) + LamTU*ZHR(gt3,1)*ZP(gt2,1)) + 2*ZHR(gt3,1)*(-1.4142135623730951*LamSU*ZA(gt1,2)*ZP(gt2,3) + LamTU*ZA(gt1,3)*(ZP(gt2,2) + ZP(gt2,3)))));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::Hpm, fields::SRdp, typename fields::conj<fields::SRum>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.5*(2*LamSD*Conj(LamSU) - LamTD*Conj(LamTU))*(ZP(gt1,1)*ZP(gt2,0) + ZP(gt1,0)*ZP(gt2,1));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Hpm, fields::Su, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.7071067811865475*(SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gt4,j3))*ZH(gt1,0)*ZP(gt2,0) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZD(gt4,j3))*ZH(gt1,1)*ZP(gt2,1) + SUM(j3,0,2,Conj(ZU(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yd(j2,j1))*ZD(gt4,3 + j2)))*(ZH(gt1,1)*ZP(gt2,0) + ZH(gt1,0)*ZP(gt2,1))) - 0.25*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt3,j1))*ZD(gt4,j1))*(1.4142135623730951*ZH(gt1,0)*ZP(gt2,0) + 1.4142135623730951*ZH(gt1,1)*ZP(gt2,1) + 2*ZH(gt1,3)*(ZP(gt2,2) - ZP(gt2,3)));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::Su, typename fields::conj<fields::Sd>::type, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = -0.7071067811865475*g2*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZD(gt2,j1));

   return {result, minuend_index, subtrahend_index};
}

ScalarVertex VertexImpl<fields::Sd, fields::sigmaO, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const auto g3 = MODELPARAMETER(g3);
   const auto MDGoc = MODELPARAMETER(MDGoc);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = std::complex<double>(0,-0.5)*g3*(MDGoc - Conj(MDGoc))*(SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1)) - SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1)));

   return {result};
}

ScalarVertex VertexImpl<fields::sigmaO, fields::Su, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt3 = indices[1];
   const auto g3 = MODELPARAMETER(g3);
   const auto MDGoc = MODELPARAMETER(MDGoc);
   const auto ZU = MODELPARAMETER(ZU);

   const std::complex<double> result = std::complex<double>(0,-0.5)*g3*(MDGoc - Conj(MDGoc))*(SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1)) - SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1)));

   return {result};
}

ScalarVertex VertexImpl<fields::phiO, fields::Sd, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt3 = indices[1];
   const auto g3 = MODELPARAMETER(g3);
   const auto MDGoc = MODELPARAMETER(MDGoc);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = -0.5*g3*(MDGoc + Conj(MDGoc))*(SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1)) - SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1)));

   return {result};
}

ScalarVertex VertexImpl<fields::phiO, fields::Su, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt3 = indices[1];
   const auto g3 = MODELPARAMETER(g3);
   const auto MDGoc = MODELPARAMETER(MDGoc);
   const auto ZU = MODELPARAMETER(ZU);

   const std::complex<double> result = -0.5*g3*(MDGoc + Conj(MDGoc))*(SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1)) - SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Rh, fields::Sd, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto MuD = MODELPARAMETER(MuD);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto vS = MODELPARAMETER(vS);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto vT = MODELPARAMETER(vT);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = 0.5*(2*MuD + 1.4142135623730951*LamSD*vS + LamTD*vT)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1)))*ZD(gt3,j2))*ZHR(gt1,0);

   return {result};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha1>::type, fields::Fu, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Yd = MODELPARAMETER(Yd);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> left = Conj(UM1(gt1,1))*SUM(j2,0,2,Conj(ZUL(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1)));

   const std::complex<double> right = 0;

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, fields::Cha1, fields::Sd>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto UM1 = MODELPARAMETER(UM1);

   const std::complex<double> left = 0;

   const std::complex<double> right = SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt3,3 + j1)))*ZUL(gt1,j2))*UM1(gt2,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Chi>::type, fields::Fd, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZN2 = MODELPARAMETER(ZN2);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZDR = MODELPARAMETER(ZDR);
   const auto ZN1 = MODELPARAMETER(ZN1);

   const std::complex<double> left = -(Conj(ZN2(gt1,2))*SUM(j2,0,2,Conj(ZDL(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1))));

   const std::complex<double> right = -0.3651483716701107*g1*SUM(j1,0,2,ZD(gt3,3 + j1)*ZDR(gt2,j1))*ZN1(gt1,0);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Chi, fields::Sd>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZN1 = MODELPARAMETER(ZN1);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZDR = MODELPARAMETER(ZDR);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto ZN2 = MODELPARAMETER(ZN2);

   const std::complex<double> left = -0.3651483716701107*g1*Conj(ZN1(gt2,0))*SUM(j1,0,2,Conj(ZD(gt3,3 + j1))*Conj(ZDR(gt1,j1)));

   const std::complex<double> right = -(SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt3,3 + j1)))*ZDL(gt1,j2))*ZN2(gt2,2));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Cha2, fields::Fu, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Yu = MODELPARAMETER(Yu);
   const auto UM2 = MODELPARAMETER(UM2);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZUR = MODELPARAMETER(ZUR);
   const auto UP2 = MODELPARAMETER(UP2);

   const std::complex<double> left = -(g2*Conj(UM2(gt1,0))*SUM(j1,0,2,Conj(ZUL(gt2,j1))*ZD(gt3,j1)));

   const std::complex<double> right = SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*ZUR(gt2,j1))*ZD(gt3,j2))*UP2(gt1,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha2>::type, typename fields::bar<fields::Fu>::type, fields::Sd>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Yu = MODELPARAMETER(Yu);
   const auto UP2 = MODELPARAMETER(UP2);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZUR = MODELPARAMETER(ZUR);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto UM2 = MODELPARAMETER(UM2);

   const std::complex<double> left = Conj(UP2(gt1,1))*SUM(j2,0,2,Conj(ZD(gt3,j2))*SUM(j1,0,2,Conj(ZUR(gt2,j1))*Yu(j1,j2)));

   const std::complex<double> right = -(g2*SUM(j1,0,2,Conj(ZD(gt3,j1))*ZUL(gt2,j1))*UM2(gt1,0));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Chi, fields::Fd, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZN1 = MODELPARAMETER(ZN1);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZDR = MODELPARAMETER(ZDR);
   const auto ZN2 = MODELPARAMETER(ZN2);

   const std::complex<double> left = -0.2357022603955158*(0.7745966692414834*g1*Conj(ZN1(gt1,0)) - 3*g2*Conj(ZN1(gt1,1)))*SUM(j1,0,2,Conj(ZDL(gt2,j1))*ZD(gt3,j1));

   const std::complex<double> right = -(SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*ZDR(gt2,j1))*ZD(gt3,j2))*ZN2(gt1,2));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Chi>::type, typename fields::bar<fields::Fd>::type, fields::Sd>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZN2 = MODELPARAMETER(ZN2);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZDR = MODELPARAMETER(ZDR);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto ZN1 = MODELPARAMETER(ZN1);

   const std::complex<double> left = -(Conj(ZN2(gt1,2))*SUM(j2,0,2,Conj(ZD(gt3,j2))*SUM(j1,0,2,Conj(ZDR(gt2,j1))*Yd(j1,j2))));

   const std::complex<double> right = -0.2357022603955158*SUM(j1,0,2,Conj(ZD(gt3,j1))*ZDL(gt2,j1))*(0.7745966692414834*g1*ZN1(gt1,0) - 3*g2*ZN1(gt1,1));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Fd, fields::Glu, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const auto g3 = MODELPARAMETER(g3);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> left = -0.7071067811865475*g3*SUM(j1,0,2,Conj(ZDL(gt1,j1))*ZD(gt3,j1));

   const std::complex<double> right = 0;

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, typename fields::bar<fields::Glu>::type, fields::Sd>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const auto g3 = MODELPARAMETER(g3);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZDL = MODELPARAMETER(ZDL);

   const std::complex<double> left = 0;

   const std::complex<double> right = -0.7071067811865475*g3*SUM(j1,0,2,Conj(ZD(gt3,j1))*ZDL(gt1,j1));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Glu>::type, fields::Fd, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt3 = indices[1];
   const auto g3 = MODELPARAMETER(g3);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZDR = MODELPARAMETER(ZDR);

   const std::complex<double> left = 0;

   const std::complex<double> right = 0.7071067811865475*g3*SUM(j1,0,2,ZD(gt3,3 + j1)*ZDR(gt2,j1));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Glu, fields::Sd>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const auto g3 = MODELPARAMETER(g3);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZDR = MODELPARAMETER(ZDR);

   const std::complex<double> left = 0.7071067811865475*g3*SUM(j1,0,2,Conj(ZD(gt3,3 + j1))*Conj(ZDR(gt1,j1)));

   const std::complex<double> right = 0;

   return {left, right};
}

MomentumDifferenceVertex VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VG>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g3 = MODELPARAMETER(g3);

   const std::complex<double> result = -0.5*g3*KroneckerDelta(gt1,gt2);

   return {result, minuend_index, subtrahend_index};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Sd, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = 0.05*(-20*(SUM(j3,0,2,Conj(ZD(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,j1))*ZD(gt4,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt3,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gt4,j3)))*ZA(gt1,0)*ZA(gt2,0) + (Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZD(gt3,j1))*ZD(gt4,j1))*(ZA(gt1,0)*ZA(gt2,0) - ZA(gt1,1)*ZA(gt2,1)) + 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt3,3 + j1))*ZD(gt4,3 + j1))*(ZA(gt1,0)*ZA(gt2,0) - ZA(gt1,1)*ZA(gt2,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Rh, fields::Sd, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.05*((Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt4,j1)) + 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt4,3 + j1)))*(ZHR(gt1,0)*ZHR(gt3,0) - ZHR(gt1,1)*ZHR(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Sd, fields::Sd, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto g3 = MODELPARAMETER(g3);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = 0.008333333333333333*(-(Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt2,j2))*ZD(gt3,j2))) - 15*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt2,j2))*ZD(gt3,j2)) + 10*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt2,j2))*ZD(gt3,j2)) - 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt4,3 + j1))*SUM(j2,0,2,Conj(ZD(gt2,j2))*ZD(gt3,j2)) - 10*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt4,3 + j1))*SUM(j2,0,2,Conj(ZD(gt2,j2))*ZD(gt3,j2)) - 30*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt4,j1))*(SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) - SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2))) + 30*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt4,3 + j1))*(SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) - SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2))) - 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt2,3 + j2))*ZD(gt3,3 + j2)) - 10*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt2,3 + j2))*ZD(gt3,3 + j2)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt4,3 + j1))*SUM(j2,0,2,Conj(ZD(gt2,3 + j2))*ZD(gt3,3 + j2)) + 10*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt4,3 + j1))*SUM(j2,0,2,Conj(ZD(gt2,3 + j2))*ZD(gt3,3 + j2)) - Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt4,j2)) - 15*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt4,j2)) + 10*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt4,j2)) - 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt4,j2)) - 10*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt4,j2)) - 30*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt2,j2))*ZD(gt4,j2)) + 30*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZD(gt2,j2))*ZD(gt4,j2)) - 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt4,3 + j2)) - 10*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt4,3 + j2)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt4,3 + j2)) + 10*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt4,3 + j2)) + 30*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt2,3 + j2))*ZD(gt4,3 + j2)) - 30*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZD(gt2,3 + j2))*ZD(gt4,3 + j2)) - 120*SUM(j2,0,2,Conj(ZD(gt1,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt4,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yd(j3,j4))*Conj(ZD(gt2,3 + j3)))*ZD(gt3,j4)) - 120*SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yd(j3,j4))*Conj(ZD(gt1,3 + j3)))*ZD(gt4,j4)));

   return {result};
}

ScalarVertex VertexImpl<fields::Sd, fields::Se, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = 0.025*(-2*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt4,3 + j1))*(SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) + 2*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2))) + SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt4,j1))*((Sqr(g1) - 5*Sqr(g2))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) + 2*Sqr(g1)*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2))) + Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) - 5*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) + 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) - 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt4,3 + j2)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt4,3 + j2)) - 40*SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt4,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yd(j3,j4))*Conj(ZD(gt1,3 + j3)))*ZD(gt3,j4)) - 40*SUM(j2,0,2,Conj(ZD(gt1,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Ye(j3,j4))*Conj(ZE(gt2,3 + j3)))*ZE(gt4,j4)));

   return {result};
}

ScalarVertex VertexImpl<fields::Sd, fields::Su, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto g3 = MODELPARAMETER(g3);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = 0.008333333333333333*(-(SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt4,j1))*((Sqr(g1) - 15*Sqr(g2) - 10*Sqr(g3))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) + 2*(Sqr(g1) + 5*Sqr(g3))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2)))) + 2*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt4,3 + j1))*((2*Sqr(g1) - 5*Sqr(g3))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) + (4*Sqr(g1) + 5*Sqr(g3))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2))) - Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) + 15*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) + 10*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) - 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) - 10*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) + 4*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)) - 10*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)) + 8*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)) + 10*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)));

   return {result};
}

ScalarVertex VertexImpl<fields::Sd, fields::Sv, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = 0.05*KroneckerDelta(gt2,gt4)*((Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1)) + 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VG, fields::VG>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g3 = MODELPARAMETER(g3);

   const std::complex<double> result = 0.25*KroneckerDelta(gt1,gt2)*Sqr(g3);

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VP, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05555555555555555*(Sqr(0.7745966692414834*g1*Cos(ThetaW) - 3*g2*Sin(ThetaW))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1)) + 2.4*Sqr(g1)*Sqr(Cos(ThetaW))*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::VWm>::type, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = 0.5*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05555555555555555*(Sqr(3*g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1)) + 2.4*Sqr(g1)*Sqr(Sin(ThetaW))*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Sd, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = std::complex<double>(0,-0.5)*SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt4,3 + j1)))*(1.4142135623730951*Conj(LamSD)*ZA(gt1,2) + Conj(LamTD)*ZA(gt1,3))*ZHR(gt3,0);

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Rh, fields::Sd, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = std::complex<double>(0,0.5)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt3,3 + j1)))*ZD(gt4,j2))*(1.4142135623730951*LamSD*ZA(gt1,2) + LamTD*ZA(gt1,3))*ZHR(gt2,0);

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Rh, fields::Sd, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = 0.5*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt3,3 + j1)))*ZD(gt4,j2))*(1.4142135623730951*LamSD*ZH(gt1,2) + LamTD*ZH(gt1,3))*ZHR(gt2,0);

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Rh, fields::Se, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = 0.5*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt3,3 + j1)))*ZE(gt4,j2))*(1.4142135623730951*LamSD*ZH(gt1,2) + LamTD*ZH(gt1,3))*ZHR(gt2,0);

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Rh, fields::Su, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.5*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt3,3 + j1)))*ZU(gt4,j2))*(1.4142135623730951*LamSU*ZH(gt1,2) - LamTU*ZH(gt1,3))*ZHR(gt2,1);

   return {result};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fd, fields::VG>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g3 = MODELPARAMETER(g3);

   const std::complex<double> left = -0.5*g3*KroneckerDelta(gt1,gt2);

   const std::complex<double> right = -0.5*g3*KroneckerDelta(gt1,gt2);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fu, fields::VG>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g3 = MODELPARAMETER(g3);

   const std::complex<double> left = -0.5*g3*KroneckerDelta(gt1,gt2);

   const std::complex<double> right = -0.5*g3*KroneckerDelta(gt1,gt2);

   return {left, right};
}

ScalarVertex VertexImpl<fields::Rh, fields::Se, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto MuD = MODELPARAMETER(MuD);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto vS = MODELPARAMETER(vS);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto vT = MODELPARAMETER(vT);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = 0.5*(2*MuD + 1.4142135623730951*LamSD*vS + LamTD*vT)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1)))*ZE(gt3,j2))*ZHR(gt1,0);

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VG>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g3 = MODELPARAMETER(g3);

   const std::complex<double> result = -0.5*g3*KroneckerDelta(gt1,gt2);

   return {result, minuend_index, subtrahend_index};
}

ScalarVertex VertexImpl<fields::Rh, fields::Su, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto MuU = MODELPARAMETER(MuU);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto vS = MODELPARAMETER(vS);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto vT = MODELPARAMETER(vT);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.5*(2*MuU + 1.4142135623730951*LamSU*vS - LamTU*vT)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1)))*ZU(gt3,j2))*ZHR(gt1,1);

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::Sv, typename fields::conj<fields::Se>::type, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);

   const std::complex<double> result = -0.7071067811865475*g2*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZE(gt2,j1));

   return {result, minuend_index, subtrahend_index};
}

ChiralVertex VertexImpl<fields::Cha1, fields::Fe, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZER = MODELPARAMETER(ZER);
   const auto UM1 = MODELPARAMETER(UM1);

   const std::complex<double> left = -(g2*Conj(UP1(gt1,0))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZV(gt3,j1)));

   const std::complex<double> right = SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZV(gt3,j2))*UM1(gt1,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha1>::type, typename fields::bar<fields::Fe>::type, fields::Sv>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto UP1 = MODELPARAMETER(UP1);

   const std::complex<double> left = Conj(UM1(gt1,1))*SUM(j2,0,2,Conj(ZV(gt3,j2))*SUM(j1,0,2,Conj(ZER(gt2,j1))*Ye(j1,j2)));

   const std::complex<double> right = -(g2*SUM(j1,0,2,Conj(ZV(gt3,j1))*ZEL(gt2,j1))*UP1(gt1,0));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Chi, fields::Fv, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZN1 = MODELPARAMETER(ZN1);
   const auto ZV = MODELPARAMETER(ZV);

   const std::complex<double> left = IF(gt2 < 3,0.5477225575051661*g1*Conj(ZN1(gt1,0))*ZV(gt3,gt2),0) + IF(gt2 < 3,-0.7071067811865475*g2*Conj(ZN1(gt1,1))*ZV(gt3,gt2),0);

   const std::complex<double> right = 0;

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Chi>::type, typename fields::bar<fields::Fv>::type, fields::Sv>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZN1 = MODELPARAMETER(ZN1);

   const std::complex<double> left = 0;

   const std::complex<double> right = IF(gt2 < 3,0.5477225575051661*g1*Conj(ZV(gt3,gt2))*ZN1(gt1,0),0) + IF(gt2 < 3,-0.7071067811865475*g2*Conj(ZV(gt3,gt2))*ZN1(gt1,1),0);

   return {left, right};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Sv, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = -0.25*KroneckerDelta(gt3,gt4)*(0.6*Sqr(g1) + Sqr(g2))*(ZA(gt1,0)*ZA(gt2,0) - ZA(gt1,1)*ZA(gt2,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Rh, fields::Sv, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = 0.25*KroneckerDelta(gt2,gt4)*(0.6*Sqr(g1) + Sqr(g2))*(ZHR(gt1,0)*ZHR(gt3,0) - ZHR(gt1,1)*ZHR(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Se, fields::Sv, typename fields::conj<fields::Se>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);

   const std::complex<double> result = 0.05*(KroneckerDelta(gt2,gt4)*((-3*Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1)) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1))) - 5*(Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZV(gt4,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZE(gt3,j2)) + Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZV(gt4,j2)) + 4*SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Ye(j3,j4))*Conj(ZE(gt1,3 + j3)))*ZV(gt4,j4))));

   return {result};
}

ScalarVertex VertexImpl<fields::Su, fields::Sv, typename fields::conj<fields::Su>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZU = MODELPARAMETER(ZU);

   const std::complex<double> result = 0.05*KroneckerDelta(gt2,gt4)*((Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt3,j1)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt3,3 + j1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Sv, fields::Sv, typename fields::conj<fields::Sv>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> result = -0.25*(KroneckerDelta(gt1,gt4)*KroneckerDelta(gt2,gt3) + KroneckerDelta(gt1,gt3)*KroneckerDelta(gt2,gt4))*(0.6*Sqr(g1) + Sqr(g2));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sv, typename fields::conj<fields::Sv>::type, typename fields::conj<fields::VWm>::type, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> result = 0.5*KroneckerDelta(gt1,gt2)*Sqr(g2);

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sv, typename fields::conj<fields::Sv>::type, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*KroneckerDelta(gt1,gt2)*Sqr(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Hpm, fields::Sv, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.7071067811865475*SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gt3,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZE(gt4,j3))*ZH(gt1,0)*ZP(gt2,0) - 0.25*Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt3,j1))*ZE(gt4,j1))*(1.4142135623730951*ZH(gt1,0)*ZP(gt2,0) + 1.4142135623730951*ZH(gt1,1)*ZP(gt2,1) + 2*ZH(gt1,3)*(ZP(gt2,2) - ZP(gt2,3)));

   return {result};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha2>::type, fields::Fd, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Yu = MODELPARAMETER(Yu);
   const auto UP2 = MODELPARAMETER(UP2);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto ZU = MODELPARAMETER(ZU);

   const std::complex<double> left = Conj(UP2(gt1,1))*SUM(j2,0,2,Conj(ZDL(gt2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt3,3 + j1)));

   const std::complex<double> right = 0;

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Cha2, fields::Su>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto UP2 = MODELPARAMETER(UP2);

   const std::complex<double> left = 0;

   const std::complex<double> right = SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt3,3 + j1)))*ZDL(gt1,j2))*UP2(gt2,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Chi>::type, fields::Fu, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZN2 = MODELPARAMETER(ZN2);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZUR = MODELPARAMETER(ZUR);
   const auto ZN1 = MODELPARAMETER(ZN1);

   const std::complex<double> left = -(Conj(ZN2(gt1,3))*SUM(j2,0,2,Conj(ZUL(gt2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt3,3 + j1))));

   const std::complex<double> right = 0.7302967433402214*g1*SUM(j1,0,2,ZU(gt3,3 + j1)*ZUR(gt2,j1))*ZN1(gt1,0);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, fields::Chi, fields::Su>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZN1 = MODELPARAMETER(ZN1);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZUR = MODELPARAMETER(ZUR);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto ZN2 = MODELPARAMETER(ZN2);

   const std::complex<double> left = 0.7302967433402214*g1*Conj(ZN1(gt2,0))*SUM(j1,0,2,Conj(ZU(gt3,3 + j1))*Conj(ZUR(gt1,j1)));

   const std::complex<double> right = -(SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt3,3 + j1)))*ZUL(gt1,j2))*ZN2(gt2,3));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Cha1, fields::Fd, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Yd = MODELPARAMETER(Yd);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZDR = MODELPARAMETER(ZDR);
   const auto UM1 = MODELPARAMETER(UM1);

   const std::complex<double> left = -(g2*Conj(UP1(gt1,0))*SUM(j1,0,2,Conj(ZDL(gt2,j1))*ZU(gt3,j1)));

   const std::complex<double> right = SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*ZDR(gt2,j1))*ZU(gt3,j2))*UM1(gt1,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha1>::type, typename fields::bar<fields::Fd>::type, fields::Su>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Yd = MODELPARAMETER(Yd);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZDR = MODELPARAMETER(ZDR);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto UP1 = MODELPARAMETER(UP1);

   const std::complex<double> left = Conj(UM1(gt1,1))*SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Conj(ZDR(gt2,j1))*Yd(j1,j2)));

   const std::complex<double> right = -(g2*SUM(j1,0,2,Conj(ZU(gt3,j1))*ZDL(gt2,j1))*UP1(gt1,0));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Chi, fields::Fu, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZN1 = MODELPARAMETER(ZN1);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZUR = MODELPARAMETER(ZUR);
   const auto ZN2 = MODELPARAMETER(ZN2);

   const std::complex<double> left = -0.2357022603955158*(0.7745966692414834*g1*Conj(ZN1(gt1,0)) + 3*g2*Conj(ZN1(gt1,1)))*SUM(j1,0,2,Conj(ZUL(gt2,j1))*ZU(gt3,j1));

   const std::complex<double> right = -(SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*ZUR(gt2,j1))*ZU(gt3,j2))*ZN2(gt1,3));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Chi>::type, typename fields::bar<fields::Fu>::type, fields::Su>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZN2 = MODELPARAMETER(ZN2);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZUR = MODELPARAMETER(ZUR);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto ZN1 = MODELPARAMETER(ZN1);

   const std::complex<double> left = -(Conj(ZN2(gt1,3))*SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Conj(ZUR(gt2,j1))*Yu(j1,j2))));

   const std::complex<double> right = -0.2357022603955158*SUM(j1,0,2,Conj(ZU(gt3,j1))*ZUL(gt2,j1))*(0.7745966692414834*g1*ZN1(gt1,0) + 3*g2*ZN1(gt1,1));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Fu, fields::Glu, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const auto g3 = MODELPARAMETER(g3);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto ZU = MODELPARAMETER(ZU);

   const std::complex<double> left = -0.7071067811865475*g3*SUM(j1,0,2,Conj(ZUL(gt1,j1))*ZU(gt3,j1));

   const std::complex<double> right = 0;

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, typename fields::bar<fields::Glu>::type, fields::Su>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const auto g3 = MODELPARAMETER(g3);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZUL = MODELPARAMETER(ZUL);

   const std::complex<double> left = 0;

   const std::complex<double> right = -0.7071067811865475*g3*SUM(j1,0,2,Conj(ZU(gt3,j1))*ZUL(gt1,j1));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Glu>::type, fields::Fu, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt3 = indices[1];
   const auto g3 = MODELPARAMETER(g3);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZUR = MODELPARAMETER(ZUR);

   const std::complex<double> left = 0;

   const std::complex<double> right = 0.7071067811865475*g3*SUM(j1,0,2,ZU(gt3,3 + j1)*ZUR(gt2,j1));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, fields::Glu, fields::Su>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const auto g3 = MODELPARAMETER(g3);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZUR = MODELPARAMETER(ZUR);

   const std::complex<double> left = 0.7071067811865475*g3*SUM(j1,0,2,Conj(ZU(gt3,3 + j1))*Conj(ZUR(gt1,j1)));

   const std::complex<double> right = 0;

   return {left, right};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Su, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = 0.05*((Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZU(gt3,j1))*ZU(gt4,j1))*(ZA(gt1,0)*ZA(gt2,0) - ZA(gt1,1)*ZA(gt2,1)) - 4*(5*(SUM(j3,0,2,Conj(ZU(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))*ZU(gt4,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gt4,j3)))*ZA(gt1,1)*ZA(gt2,1) + Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt3,3 + j1))*ZU(gt4,3 + j1))*(ZA(gt1,0)*ZA(gt2,0) - ZA(gt1,1)*ZA(gt2,1))));

   return {result};
}

ScalarVertex VertexImpl<fields::Rh, fields::Su, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.05*((Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt4,j1)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt4,3 + j1)))*(ZHR(gt1,0)*ZHR(gt3,0) - ZHR(gt1,1)*ZHR(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Se, fields::Su, typename fields::conj<fields::Se>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZE = MODELPARAMETER(ZE);

   const std::complex<double> result = 0.025*(-4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt4,3 + j1))*(SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt3,j2)) - 2*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt3,3 + j2))) + SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt4,j1))*((Sqr(g1) + 5*Sqr(g2))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt3,j2)) - 2*Sqr(g1)*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt3,3 + j2))) + Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) + 5*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) - 2*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)) + 8*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)));

   return {result};
}

ScalarVertex VertexImpl<fields::Su, fields::Su, typename fields::conj<fields::Su>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto g3 = MODELPARAMETER(g3);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);

   const std::complex<double> result = 0.008333333333333333*(-(Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt4,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt3,j2))) - 15*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt4,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt3,j2)) + 10*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt4,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt3,j2)) + 4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt4,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt3,j2)) - 10*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt4,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt3,j2)) - 30*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt4,j1))*(SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt3,j2)) - SUM(j2,0,2,Conj(ZU(gt1,3 + j2))*ZU(gt3,3 + j2))) + 30*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt4,3 + j1))*(SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt3,j2)) - SUM(j2,0,2,Conj(ZU(gt1,3 + j2))*ZU(gt3,3 + j2))) + 4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt4,j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt3,3 + j2)) - 10*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt4,j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt3,3 + j2)) - 16*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt4,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt3,3 + j2)) + 10*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt4,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt3,3 + j2)) - Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt4,j2)) - 15*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt4,j2)) + 10*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt4,j2)) + 4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt4,j2)) - 10*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt4,j2)) - 30*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) + 30*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) + 4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt1,3 + j2))*ZU(gt4,3 + j2)) - 10*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt1,3 + j2))*ZU(gt4,3 + j2)) - 16*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt1,3 + j2))*ZU(gt4,3 + j2)) + 10*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt1,3 + j2))*ZU(gt4,3 + j2)) + 30*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)) - 30*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)) - 120*SUM(j2,0,2,Conj(ZU(gt1,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt4,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yu(j3,j4))*Conj(ZU(gt2,3 + j3)))*ZU(gt3,j4)) - 120*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt3,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yu(j3,j4))*Conj(ZU(gt1,3 + j3)))*ZU(gt4,j4)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VG, fields::VG>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g3 = MODELPARAMETER(g3);

   const std::complex<double> result = 0.25*KroneckerDelta(gt1,gt2)*Sqr(g3);

   return {result};
}

InverseMetricVertex VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VP, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05555555555555555*(Sqr(0.7745966692414834*g1*Cos(ThetaW) + 3*g2*Sin(ThetaW))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1)) + 9.6*Sqr(g1)*Sqr(Cos(ThetaW))*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, typename fields::conj<fields::VWm>::type, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZU = MODELPARAMETER(ZU);

   const std::complex<double> result = 0.5*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05555555555555555*(Sqr(-3*g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1)) + 9.6*Sqr(g1)*Sqr(Sin(ThetaW))*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Su, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = std::complex<double>(0,0.5)*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt4,3 + j1)))*(1.4142135623730951*Conj(LamSU)*ZA(gt1,2) - Conj(LamTU)*ZA(gt1,3))*ZHR(gt3,1);

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Rh, fields::Su, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = std::complex<double>(0,-0.5)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt3,3 + j1)))*ZU(gt4,j2))*(1.4142135623730951*LamSU*ZA(gt1,2) - LamTU*ZA(gt1,3))*ZHR(gt2,1);

   return {result};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha1>::type, fields::Fv, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto ZE = MODELPARAMETER(ZE);

   const std::complex<double> left = Conj(UM1(gt1,1))*SUM(j1,0,2,Ye(j1,gt2)*ZE(gt3,3 + j1));

   const std::complex<double> right = 0;

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fv>::type, fields::Cha1, fields::Se>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto UM1 = MODELPARAMETER(UM1);

   const std::complex<double> left = 0;

   const std::complex<double> right = SUM(j1,0,2,Conj(Ye(j1,gt1))*Conj(ZE(gt3,3 + j1)))*UM1(gt2,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Chi>::type, fields::Fe, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZN2 = MODELPARAMETER(ZN2);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZN1 = MODELPARAMETER(ZN1);

   const std::complex<double> left = -(Conj(ZN2(gt1,2))*SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1))));

   const std::complex<double> right = -1.0954451150103321*g1*SUM(j1,0,2,ZE(gt3,3 + j1)*ZER(gt2,j1))*ZN1(gt1,0);

   return {left, right};
}

ChiralVertex VertexImpl<fields::Cha2, fields::Fv, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto UM2 = MODELPARAMETER(UM2);
   const auto ZE = MODELPARAMETER(ZE);

   const std::complex<double> left = IF(gt2 < 3,-(g2*Conj(UM2(gt1,0))*ZE(gt3,gt2)),0);

   const std::complex<double> right = 0;

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha2>::type, typename fields::bar<fields::Fv>::type, fields::Se>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZE = MODELPARAMETER(ZE);
   const auto UM2 = MODELPARAMETER(UM2);

   const std::complex<double> left = 0;

   const std::complex<double> right = IF(gt2 < 3,-(g2*Conj(ZE(gt3,gt2))*UM2(gt1,0)),0);

   return {left, right};
}

ChiralVertex VertexImpl<fields::Chi, fields::Fe, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZN1 = MODELPARAMETER(ZN1);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZN2 = MODELPARAMETER(ZN2);

   const std::complex<double> left = 0.7071067811865475*(0.7745966692414834*g1*Conj(ZN1(gt1,0)) + g2*Conj(ZN1(gt1,1)))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZE(gt3,j1));

   const std::complex<double> right = -(SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZE(gt3,j2))*ZN2(gt1,2));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Chi>::type, typename fields::bar<fields::Fe>::type, fields::Se>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZN2 = MODELPARAMETER(ZN2);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZN1 = MODELPARAMETER(ZN1);

   const std::complex<double> left = -(Conj(ZN2(gt1,2))*SUM(j2,0,2,Conj(ZE(gt3,j2))*SUM(j1,0,2,Conj(ZER(gt2,j1))*Ye(j1,j2))));

   const std::complex<double> right = 0.7071067811865475*SUM(j1,0,2,Conj(ZE(gt3,j1))*ZEL(gt2,j1))*(0.7745966692414834*g1*ZN1(gt1,0) + g2*ZN1(gt1,1));

   return {left, right};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Se, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = 0.05*(-20*(SUM(j3,0,2,Conj(ZE(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Ye(j2,j1))*ZE(gt4,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gt3,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZE(gt4,j3)))*ZA(gt1,0)*ZA(gt2,0) - (3*Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZE(gt3,j1))*ZE(gt4,j1))*(ZA(gt1,0)*ZA(gt2,0) - ZA(gt1,1)*ZA(gt2,1)) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*ZE(gt4,3 + j1))*(ZA(gt1,0)*ZA(gt2,0) - ZA(gt1,1)*ZA(gt2,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Rh, fields::Se, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = 0.05*((3*Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt4,j1)) - 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt4,3 + j1)))*(ZHR(gt1,0)*ZHR(gt3,0) - ZHR(gt1,1)*ZHR(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Se, fields::Se, typename fields::conj<fields::Se>::type, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);

   const std::complex<double> result = 0.025*(-3*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt4,j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt3,j2)) - 5*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt4,j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt3,j2)) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt4,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt3,j2)) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt4,3 + j1))*(SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt3,j2)) - 2*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt3,3 + j2))) + SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt4,j1))*(-((3*Sqr(g1) + 5*Sqr(g2))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt3,j2))) + 6*Sqr(g1)*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt3,3 + j2))) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt4,j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt3,3 + j2)) - 12*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt4,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt3,3 + j2)) - 3*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt4,j2)) - 5*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt4,j2)) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt4,j2)) - 3*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) - 5*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt4,3 + j2)) - 12*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt4,3 + j2)) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt4,3 + j2)) - 12*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt4,3 + j2)) - 40*SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt4,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Ye(j3,j4))*Conj(ZE(gt1,3 + j3)))*ZE(gt3,j4)) - 40*SUM(j2,0,2,Conj(ZE(gt1,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt4,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Ye(j3,j4))*Conj(ZE(gt2,3 + j3)))*ZE(gt3,j4)) - 40*SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Ye(j3,j4))*Conj(ZE(gt1,3 + j3)))*ZE(gt4,j4)) - 40*SUM(j2,0,2,Conj(ZE(gt1,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Ye(j3,j4))*Conj(ZE(gt2,3 + j3)))*ZE(gt4,j4)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Se, typename fields::conj<fields::Se>::type, fields::VP, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(Sqr(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) + 2.4*Sqr(g1)*Sqr(Cos(ThetaW))*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Se, typename fields::conj<fields::Se>::type, typename fields::conj<fields::VWm>::type, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZE = MODELPARAMETER(ZE);

   const std::complex<double> result = 0.5*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Se, typename fields::conj<fields::Se>::type, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(Sqr(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) + 2.4*Sqr(g1)*Sqr(Sin(ThetaW))*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Se, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = std::complex<double>(0,-0.5)*SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt4,3 + j1)))*(1.4142135623730951*Conj(LamSD)*ZA(gt1,2) + Conj(LamTD)*ZA(gt1,3))*ZHR(gt3,0);

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Rh, fields::Se, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = std::complex<double>(0,0.5)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt3,3 + j1)))*ZE(gt4,j2))*(1.4142135623730951*LamSD*ZA(gt1,2) + LamTD*ZA(gt1,3))*ZHR(gt2,0);

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Ah, fields::Ah>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = 0.05*(-(ZA(gt1,0)*(-((3*Sqr(g1) + 5*Sqr(g2))*ZA(gt2,1)*(ZA(gt3,1)*ZA(gt4,0) + ZA(gt3,0)*ZA(gt4,1))) + 5*Conj(LamSD)*(1.4142135623730951*LamTD*ZA(gt2,3)*(ZA(gt3,2)*ZA(gt4,0) + ZA(gt3,0)*ZA(gt4,2)) + ZA(gt2,2)*(4*LamSD*ZA(gt3,2)*ZA(gt4,0) + 1.4142135623730951*LamTD*ZA(gt3,3)*ZA(gt4,0) + ZA(gt3,0)*(4*LamSD*ZA(gt4,2) + 1.4142135623730951*LamTD*ZA(gt4,3)))) + 5*Conj(LamTD)*(1.4142135623730951*LamSD*ZA(gt2,2)*(ZA(gt3,3)*ZA(gt4,0) + ZA(gt3,0)*ZA(gt4,3)) + ZA(gt2,3)*(1.4142135623730951*LamSD*ZA(gt3,2)*ZA(gt4,0) + 2*LamTD*ZA(gt3,3)*ZA(gt4,0) + ZA(gt3,0)*(1.4142135623730951*LamSD*ZA(gt4,2) + 2*LamTD*ZA(gt4,3)))) + ZA(gt2,0)*(3*(3*Sqr(g1) + 5*Sqr(g2))*ZA(gt3,0)*ZA(gt4,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZA(gt3,1)*ZA(gt4,1) + 5*Conj(LamSD)*(1.4142135623730951*LamTD*ZA(gt3,3)*ZA(gt4,2) + ZA(gt3,2)*(4*LamSD*ZA(gt4,2) + 1.4142135623730951*LamTD*ZA(gt4,3))) + 5*Conj(LamTD)*(1.4142135623730951*LamSD*ZA(gt3,2)*ZA(gt4,3) + ZA(gt3,3)*(1.4142135623730951*LamSD*ZA(gt4,2) + 2*LamTD*ZA(gt4,3)))))) + ZA(gt1,1)*((3*Sqr(g1) + 5*Sqr(g2))*ZA(gt2,0)*(ZA(gt3,1)*ZA(gt4,0) + ZA(gt3,0)*ZA(gt4,1)) + 5*Conj(LamTU)*(1.4142135623730951*LamSU*ZA(gt2,2)*(ZA(gt3,3)*ZA(gt4,1) + ZA(gt3,1)*ZA(gt4,3)) + ZA(gt2,3)*(1.4142135623730951*LamSU*ZA(gt3,2)*ZA(gt4,1) - 2*LamTU*ZA(gt3,3)*ZA(gt4,1) + ZA(gt3,1)*(1.4142135623730951*LamSU*ZA(gt4,2) - 2*LamTU*ZA(gt4,3)))) + 5*Conj(LamSU)*(1.4142135623730951*LamTU*ZA(gt2,3)*(ZA(gt3,2)*ZA(gt4,1) + ZA(gt3,1)*ZA(gt4,2)) + ZA(gt2,2)*(-4*LamSU*ZA(gt3,2)*ZA(gt4,1) + 1.4142135623730951*LamTU*ZA(gt3,3)*ZA(gt4,1) + ZA(gt3,1)*(-4*LamSU*ZA(gt4,2) + 1.4142135623730951*LamTU*ZA(gt4,3)))) + ZA(gt2,1)*((3*Sqr(g1) + 5*Sqr(g2))*ZA(gt3,0)*ZA(gt4,0) - 3*(3*Sqr(g1) + 5*Sqr(g2))*ZA(gt3,1)*ZA(gt4,1) + 5*Conj(LamTU)*(1.4142135623730951*LamSU*ZA(gt3,2)*ZA(gt4,3) + ZA(gt3,3)*(1.4142135623730951*LamSU*ZA(gt4,2) - 2*LamTU*ZA(gt4,3))) + 5*Conj(LamSU)*(1.4142135623730951*LamTU*ZA(gt3,3)*ZA(gt4,2) + ZA(gt3,2)*(-4*LamSU*ZA(gt4,2) + 1.4142135623730951*LamTU*ZA(gt4,3))))) + 5*(-4*AbsSqr(LamSU)*ZA(gt1,2)*ZA(gt2,2)*ZA(gt3,1)*ZA(gt4,1) + 1.4142135623730951*LamTU*Conj(LamSU)*ZA(gt1,3)*ZA(gt2,2)*ZA(gt3,1)*ZA(gt4,1) + 1.4142135623730951*LamSU*Conj(LamTU)*ZA(gt1,3)*ZA(gt2,2)*ZA(gt3,1)*ZA(gt4,1) + 1.4142135623730951*LamTU*Conj(LamSU)*ZA(gt1,2)*ZA(gt2,3)*ZA(gt3,1)*ZA(gt4,1) + 1.4142135623730951*LamSU*Conj(LamTU)*ZA(gt1,2)*ZA(gt2,3)*ZA(gt3,1)*ZA(gt4,1) - 2*AbsSqr(LamTU)*ZA(gt1,3)*ZA(gt2,3)*ZA(gt3,1)*ZA(gt4,1) - 4*AbsSqr(LamSU)*ZA(gt1,2)*ZA(gt2,1)*ZA(gt3,2)*ZA(gt4,1) + 1.4142135623730951*LamTU*Conj(LamSU)*ZA(gt1,3)*ZA(gt2,1)*ZA(gt3,2)*ZA(gt4,1) + 1.4142135623730951*LamSU*Conj(LamTU)*ZA(gt1,3)*ZA(gt2,1)*ZA(gt3,2)*ZA(gt4,1) + 1.4142135623730951*LamTU*Conj(LamSU)*ZA(gt1,2)*ZA(gt2,1)*ZA(gt3,3)*ZA(gt4,1) + 1.4142135623730951*LamSU*Conj(LamTU)*ZA(gt1,2)*ZA(gt2,1)*ZA(gt3,3)*ZA(gt4,1) - 2*AbsSqr(LamTU)*ZA(gt1,3)*ZA(gt2,1)*ZA(gt3,3)*ZA(gt4,1) - 4*AbsSqr(LamSU)*ZA(gt1,2)*ZA(gt2,1)*ZA(gt3,1)*ZA(gt4,2) + 1.4142135623730951*LamTU*Conj(LamSU)*ZA(gt1,3)*ZA(gt2,1)*ZA(gt3,1)*ZA(gt4,2) + 1.4142135623730951*LamSU*Conj(LamTU)*ZA(gt1,3)*ZA(gt2,1)*ZA(gt3,1)*ZA(gt4,2) + 1.4142135623730951*LamTU*Conj(LamSU)*ZA(gt1,2)*ZA(gt2,1)*ZA(gt3,1)*ZA(gt4,3) + 1.4142135623730951*LamSU*Conj(LamTU)*ZA(gt1,2)*ZA(gt2,1)*ZA(gt3,1)*ZA(gt4,3) - 2*AbsSqr(LamTU)*ZA(gt1,3)*ZA(gt2,1)*ZA(gt3,1)*ZA(gt4,3) - Conj(LamSD)*(1.4142135623730951*LamTD*ZA(gt1,3)*(ZA(gt2,2)*ZA(gt3,0)*ZA(gt4,0) + ZA(gt2,0)*(ZA(gt3,2)*ZA(gt4,0) + ZA(gt3,0)*ZA(gt4,2))) + ZA(gt1,2)*(4*LamSD*ZA(gt2,2)*ZA(gt3,0)*ZA(gt4,0) + 1.4142135623730951*LamTD*ZA(gt2,3)*ZA(gt3,0)*ZA(gt4,0) + ZA(gt2,0)*(4*LamSD*ZA(gt3,2)*ZA(gt4,0) + 1.4142135623730951*LamTD*ZA(gt3,3)*ZA(gt4,0) + ZA(gt3,0)*(4*LamSD*ZA(gt4,2) + 1.4142135623730951*LamTD*ZA(gt4,3))))) - Conj(LamTD)*(1.4142135623730951*LamSD*ZA(gt1,2)*(ZA(gt2,3)*ZA(gt3,0)*ZA(gt4,0) + ZA(gt2,0)*(ZA(gt3,3)*ZA(gt4,0) + ZA(gt3,0)*ZA(gt4,3))) + ZA(gt1,3)*(1.4142135623730951*LamSD*ZA(gt2,2)*ZA(gt3,0)*ZA(gt4,0) + 2*LamTD*ZA(gt2,3)*ZA(gt3,0)*ZA(gt4,0) + ZA(gt2,0)*(1.4142135623730951*LamSD*ZA(gt3,2)*ZA(gt4,0) + 2*LamTD*ZA(gt3,3)*ZA(gt4,0) + ZA(gt3,0)*(1.4142135623730951*LamSD*ZA(gt4,2) + 2*LamTD*ZA(gt4,3)))))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Rh, typename fields::conj<fields::Rh>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = 0.05*(-5*(Conj(LamSD)*(1.4142135623730951*LamTD*ZA(gt1,3)*ZA(gt2,2) + ZA(gt1,2)*(4*LamSD*ZA(gt2,2) + 1.4142135623730951*LamTD*ZA(gt2,3)))*ZHR(gt3,0)*ZHR(gt4,0) + Conj(LamTD)*(1.4142135623730951*LamSD*ZA(gt1,2)*ZA(gt2,3) + ZA(gt1,3)*(1.4142135623730951*LamSD*ZA(gt2,2) + 2*LamTD*ZA(gt2,3)))*ZHR(gt3,0)*ZHR(gt4,0) + (-(Conj(LamTU)*(1.4142135623730951*LamSU*ZA(gt1,2)*ZA(gt2,3) + ZA(gt1,3)*(1.4142135623730951*LamSU*ZA(gt2,2) - 2*LamTU*ZA(gt2,3)))) - Conj(LamSU)*(1.4142135623730951*LamTU*ZA(gt1,3)*ZA(gt2,2) + ZA(gt1,2)*(-4*LamSU*ZA(gt2,2) + 1.4142135623730951*LamTU*ZA(gt2,3))))*ZHR(gt3,1)*ZHR(gt4,1)) - ZA(gt1,1)*(5*ZA(gt2,0)*(-2*LamSU*Conj(LamSD)*ZHR(gt3,1)*ZHR(gt4,0) + LamTU*Conj(LamTD)*ZHR(gt3,1)*ZHR(gt4,0) + (-2*LamSD*Conj(LamSU) + LamTD*Conj(LamTU))*ZHR(gt3,0)*ZHR(gt4,1)) + ZA(gt2,1)*((3*Sqr(g1) + 5*Sqr(g2))*ZHR(gt3,0)*ZHR(gt4,0) + (20*AbsSqr(LamSU) + 10*AbsSqr(LamTU) - 3*Sqr(g1) - 5*Sqr(g2))*ZHR(gt3,1)*ZHR(gt4,1))) + ZA(gt1,0)*(5*ZA(gt2,1)*(2*LamSU*Conj(LamSD)*ZHR(gt3,1)*ZHR(gt4,0) - LamTU*Conj(LamTD)*ZHR(gt3,1)*ZHR(gt4,0) + (2*LamSD*Conj(LamSU) - LamTD*Conj(LamTU))*ZHR(gt3,0)*ZHR(gt4,1)) + ZA(gt2,0)*((-20*AbsSqr(LamSD) - 10*AbsSqr(LamTD) + 3*Sqr(g1) + 5*Sqr(g2))*ZHR(gt3,0)*ZHR(gt4,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZHR(gt3,1)*ZHR(gt4,1))));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Ah, fields::Ah, typename fields::conj<fields::VWm>::type, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = 0.5*Sqr(g2)*(ZA(gt1,0)*ZA(gt2,0) + ZA(gt1,1)*ZA(gt2,1) + 4*ZA(gt1,3)*ZA(gt2,3));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Ah, fields::Ah, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*Sqr(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(ZA(gt1,0)*ZA(gt2,0) + ZA(gt1,1)*ZA(gt2,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Rh, fields::Rh, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::Rh>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.25*(0.6*Sqr(g1) + Sqr(g2))*(-(ZHR(gt1,1)*(-2*ZHR(gt2,1)*ZHR(gt3,1)*ZHR(gt4,1) + ZHR(gt2,0)*(ZHR(gt3,1)*ZHR(gt4,0) + ZHR(gt3,0)*ZHR(gt4,1)))) + ZHR(gt1,0)*(2*ZHR(gt2,0)*ZHR(gt3,0)*ZHR(gt4,0) - ZHR(gt2,1)*(ZHR(gt3,1)*ZHR(gt4,0) + ZHR(gt3,0)*ZHR(gt4,1))));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Rh, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::VWm>::type, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = 0.5*Sqr(g2)*(ZHR(gt1,0)*ZHR(gt2,0) + ZHR(gt1,1)*ZHR(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Rh, typename fields::conj<fields::Rh>::type, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*Sqr(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(ZHR(gt1,0)*ZHR(gt2,0) + ZHR(gt1,1)*ZHR(gt2,1));

   return {result};
}

ChiralVertex VertexImpl<fields::Cha1, fields::Cha2, typename fields::conj<fields::Rh>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto UM2 = MODELPARAMETER(UM2);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto UP2 = MODELPARAMETER(UP2);

   const std::complex<double> left = -(g2*(Conj(UM2(gt2,0))*Conj(UP1(gt1,1))*ZHR(gt3,0) + Conj(UM2(gt2,1))*Conj(UP1(gt1,0))*ZHR(gt3,1)));

   const std::complex<double> right = Conj(LamTD)*UM1(gt1,1)*UP2(gt2,0)*ZHR(gt3,0) - Conj(LamTU)*UM1(gt1,0)*UP2(gt2,1)*ZHR(gt3,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha1>::type, typename fields::bar<fields::Cha2>::type, fields::Rh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto g2 = MODELPARAMETER(g2);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto UP2 = MODELPARAMETER(UP2);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto UM2 = MODELPARAMETER(UM2);
   const auto UP1 = MODELPARAMETER(UP1);

   const std::complex<double> left = LamTD*Conj(UM1(gt1,1))*Conj(UP2(gt2,0))*ZHR(gt3,0) - LamTU*Conj(UM1(gt1,0))*Conj(UP2(gt2,1))*ZHR(gt3,1);

   const std::complex<double> right = -(g2*(UM2(gt2,0)*UP1(gt1,1)*ZHR(gt3,0) + UM2(gt2,1)*UP1(gt1,0)*ZHR(gt3,1)));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Chi, fields::Chi, typename fields::conj<fields::Rh>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto ZN1 = MODELPARAMETER(ZN1);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ZN2 = MODELPARAMETER(ZN2);

   const std::complex<double> left = -0.7071067811865475*(Conj(ZN1(gt1,2))*(0.7745966692414834*g1*Conj(ZN1(gt2,0)) - g2*Conj(ZN1(gt2,1)))*ZHR(gt3,0) - g2*Conj(ZN1(gt1,1))*Conj(ZN1(gt2,2))*ZHR(gt3,0) - 0.7745966692414834*g1*Conj(ZN1(gt1,3))*Conj(ZN1(gt2,0))*ZHR(gt3,1) + g2*Conj(ZN1(gt1,3))*Conj(ZN1(gt2,1))*ZHR(gt3,1) + g2*Conj(ZN1(gt1,1))*Conj(ZN1(gt2,3))*ZHR(gt3,1) + 0.7745966692414834*g1*Conj(ZN1(gt1,0))*(Conj(ZN1(gt2,2))*ZHR(gt3,0) - Conj(ZN1(gt2,3))*ZHR(gt3,1)));

   const std::complex<double> right = Conj(LamSD)*ZHR(gt3,0)*(ZN2(gt1,2)*ZN2(gt2,0) + ZN2(gt1,0)*ZN2(gt2,2)) - Conj(LamSU)*ZHR(gt3,1)*(ZN2(gt1,3)*ZN2(gt2,0) + ZN2(gt1,0)*ZN2(gt2,3)) + 0.7071067811865475*(Conj(LamTD)*ZHR(gt3,0)*(ZN2(gt1,2)*ZN2(gt2,1) + ZN2(gt1,1)*ZN2(gt2,2)) + Conj(LamTU)*ZHR(gt3,1)*(ZN2(gt1,3)*ZN2(gt2,1) + ZN2(gt1,1)*ZN2(gt2,3)));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Chi>::type, typename fields::bar<fields::Chi>::type, fields::Rh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZN2 = MODELPARAMETER(ZN2);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ZN1 = MODELPARAMETER(ZN1);

   const std::complex<double> left = 0.5*(Conj(ZN2(gt1,2))*(2*LamSD*Conj(ZN2(gt2,0)) + 1.4142135623730951*LamTD*Conj(ZN2(gt2,1)))*ZHR(gt3,0) + 1.4142135623730951*LamTD*Conj(ZN2(gt1,1))*Conj(ZN2(gt2,2))*ZHR(gt3,0) - 2*LamSU*Conj(ZN2(gt1,3))*Conj(ZN2(gt2,0))*ZHR(gt3,1) + 1.4142135623730951*LamTU*Conj(ZN2(gt1,3))*Conj(ZN2(gt2,1))*ZHR(gt3,1) + 1.4142135623730951*LamTU*Conj(ZN2(gt1,1))*Conj(ZN2(gt2,3))*ZHR(gt3,1) + 2*Conj(ZN2(gt1,0))*(LamSD*Conj(ZN2(gt2,2))*ZHR(gt3,0) - LamSU*Conj(ZN2(gt2,3))*ZHR(gt3,1)));

   const std::complex<double> right = 0.1414213562373095*(-(ZHR(gt3,0)*(ZN1(gt1,2)*(3.872983346207417*g1*ZN1(gt2,0) - 5*g2*ZN1(gt2,1)) + (3.872983346207417*g1*ZN1(gt1,0) - 5*g2*ZN1(gt1,1))*ZN1(gt2,2))) + ZHR(gt3,1)*(ZN1(gt1,3)*(3.872983346207417*g1*ZN1(gt2,0) - 5*g2*ZN1(gt2,1)) + (3.872983346207417*g1*ZN1(gt1,0) - 5*g2*ZN1(gt1,1))*ZN1(gt2,3)));

   return {left, right};
}

ScalarVertex VertexImpl<fields::Ah, fields::SRum, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Rh>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt4 = indices[2];
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,-0.25)*(1.4142135623730951*LamTU*Conj(LamTD)*ZHR(gt4,0)*(2*ZA(gt1,1)*ZP(gt3,0) - ZA(gt1,0)*ZP(gt3,1)) + 1.4142135623730951*ZA(gt1,0)*(Sqr(g2)*ZHR(gt4,1)*ZP(gt3,0) - 2*LamSU*Conj(LamSD)*ZHR(gt4,0)*ZP(gt3,1)) - ZHR(gt4,1)*(1.4142135623730951*(-2*AbsSqr(LamSU) + AbsSqr(LamTU) + Sqr(g2))*ZA(gt1,1)*ZP(gt3,1) - 2*Sqr(g2)*ZA(gt1,3)*ZP(gt3,2) - 2.8284271247461903*LamTU*Conj(LamSU)*ZA(gt1,2)*ZP(gt3,3) - 2*Sqr(g2)*ZA(gt1,3)*ZP(gt3,3) + 2*Conj(LamTU)*(1.4142135623730951*LamSU*ZA(gt1,2)*ZP(gt3,2) + LamTU*ZA(gt1,3)*(ZP(gt3,2) + ZP(gt3,3)))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Hpm, fields::SRdp, typename fields::conj<fields::Rh>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt4 = indices[2];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,0.25)*(2.8284271247461903*LamSD*Conj(LamSU)*ZA(gt1,1)*ZHR(gt4,1)*ZP(gt2,0) + 1.4142135623730951*LamTD*Conj(LamTU)*ZA(gt1,1)*ZHR(gt4,1)*ZP(gt2,0) - 1.4142135623730951*Sqr(g2)*ZA(gt1,1)*ZHR(gt4,0)*ZP(gt2,1) + 1.4142135623730951*ZA(gt1,0)*((-2*AbsSqr(LamSD) + AbsSqr(LamTD) + Sqr(g2))*ZHR(gt4,0)*ZP(gt2,0) - 2*LamTD*Conj(LamTU)*ZHR(gt4,1)*ZP(gt2,1)) - 2.8284271247461903*LamTD*Conj(LamSD)*ZA(gt1,2)*ZHR(gt4,0)*ZP(gt2,2) - 2*AbsSqr(LamTD)*ZA(gt1,3)*ZHR(gt4,0)*ZP(gt2,2) + 2*Sqr(g2)*ZA(gt1,3)*ZHR(gt4,0)*ZP(gt2,2) + 2.8284271247461903*LamSD*Conj(LamTD)*ZA(gt1,2)*ZHR(gt4,0)*ZP(gt2,3) - 2*AbsSqr(LamTD)*ZA(gt1,3)*ZHR(gt4,0)*ZP(gt2,3) + 2*Sqr(g2)*ZA(gt1,3)*ZHR(gt4,0)*ZP(gt2,3));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.5*g2*(0.7745966692414834*g1*Cos(ThetaW)*ZH(gt1,0)*ZP(gt2,0) - 0.7745966692414834*g1*Cos(ThetaW)*ZH(gt1,1)*ZP(gt2,1) + 1.4142135623730951*g2*Sin(ThetaW)*ZH(gt1,3)*(ZP(gt2,2) + ZP(gt2,3)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.5*g2*(-0.7745966692414834*g1*Sin(ThetaW)*ZH(gt1,0)*ZP(gt2,0) + 0.7745966692414834*g1*Sin(ThetaW)*ZH(gt1,1)*ZP(gt2,1) + 1.4142135623730951*g2*Cos(ThetaW)*ZH(gt1,3)*(ZP(gt2,2) + ZP(gt2,3)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, typename fields::conj<fields::Hpm>::type, fields::VP, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.5*g2*(0.7745966692414834*g1*Cos(ThetaW)*ZH(gt1,0)*ZP(gt2,0) - 0.7745966692414834*g1*Cos(ThetaW)*ZH(gt1,1)*ZP(gt2,1) + 1.4142135623730951*g2*Sin(ThetaW)*ZH(gt1,3)*(ZP(gt2,2) + ZP(gt2,3)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, typename fields::conj<fields::Hpm>::type, fields::VWm, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.5*g2*(-0.7745966692414834*g1*Sin(ThetaW)*ZH(gt1,0)*ZP(gt2,0) + 0.7745966692414834*g1*Sin(ThetaW)*ZH(gt1,1)*ZP(gt2,1) + 1.4142135623730951*g2*Cos(ThetaW)*ZH(gt1,3)*(ZP(gt2,2) + ZP(gt2,3)));

   return {result};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha1>::type, fields::Chi, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZN1 = MODELPARAMETER(ZN1);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto ZN2 = MODELPARAMETER(ZN2);

   const std::complex<double> left = g2*Conj(ZN1(gt2,1))*UP1(gt1,0) - 0.7071067811865475*g2*Conj(ZN1(gt2,2))*UP1(gt1,1);

   const std::complex<double> right = g2*Conj(UM1(gt1,0))*ZN2(gt2,1) + 0.7071067811865475*g2*Conj(UM1(gt1,1))*ZN2(gt2,2);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Chi>::type, fields::Cha2, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto UM2 = MODELPARAMETER(UM2);
   const auto ZN1 = MODELPARAMETER(ZN1);
   const auto ZN2 = MODELPARAMETER(ZN2);
   const auto UP2 = MODELPARAMETER(UP2);

   const std::complex<double> left = -0.5*g2*(2*Conj(UM2(gt2,0))*ZN1(gt1,1) + 1.4142135623730951*Conj(UM2(gt2,1))*ZN1(gt1,3));

   const std::complex<double> right = -(g2*Conj(ZN2(gt1,1))*UP2(gt2,0)) + 0.7071067811865475*g2*Conj(ZN2(gt1,3))*UP2(gt2,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fd, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto ZUL = MODELPARAMETER(ZUL);

   const std::complex<double> left = -0.7071067811865475*g2*SUM(j1,0,2,Conj(ZDL(gt2,j1))*ZUL(gt1,j1));

   const std::complex<double> right = 0;

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fv>::type, fields::Fe, typename fields::conj<fields::VWm>::type>::evaluate(
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

ScalarVertex VertexImpl<typename fields::bar<fields::gWm>::type, fields::gP, fields::Hpm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto vT = MODELPARAMETER(vT);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.25*g2*(vd*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*ZP(gt3,0) - vu*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*ZP(gt3,1) + 2.8284271247461903*g2*vT*Sin(ThetaW)*(ZP(gt3,2) + ZP(gt3,3)));

   return {result};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gP>::type, fields::gWm, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Sin(ThetaW);

   return {result, 1};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gZ, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Cos(ThetaW));

   return {result, 1};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gWm, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Cos(ThetaW);

   return {result, 1};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha2>::type, fields::Chi, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZN1 = MODELPARAMETER(ZN1);
   const auto UM2 = MODELPARAMETER(UM2);
   const auto UP2 = MODELPARAMETER(UP2);
   const auto ZN2 = MODELPARAMETER(ZN2);

   const std::complex<double> left = -0.5*g2*(2*Conj(ZN1(gt2,1))*UM2(gt1,0) + 1.4142135623730951*Conj(ZN1(gt2,3))*UM2(gt1,1));

   const std::complex<double> right = -(g2*Conj(UP2(gt1,0))*ZN2(gt2,1)) + 0.7071067811865475*g2*Conj(UP2(gt1,1))*ZN2(gt2,3);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Chi>::type, fields::Cha1, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto ZN1 = MODELPARAMETER(ZN1);
   const auto ZN2 = MODELPARAMETER(ZN2);
   const auto UM1 = MODELPARAMETER(UM1);

   const std::complex<double> left = g2*Conj(UP1(gt2,0))*ZN1(gt1,1) - 0.7071067811865475*g2*Conj(UP1(gt2,1))*ZN1(gt1,2);

   const std::complex<double> right = g2*Conj(ZN2(gt1,1))*UM1(gt2,0) + 0.7071067811865475*g2*Conj(ZN2(gt1,2))*UM1(gt2,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fu, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto ZDL = MODELPARAMETER(ZDL);

   const std::complex<double> left = -0.7071067811865475*g2*SUM(j1,0,2,Conj(ZUL(gt2,j1))*ZDL(gt1,j1));

   const std::complex<double> right = 0;

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fv, fields::VWm>::evaluate(
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

ScalarVertex VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gP, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto vT = MODELPARAMETER(vT);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.25*g2*(vd*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*ZP(gt3,0) - vu*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*ZP(gt3,1) + 2.8284271247461903*g2*vT*Sin(ThetaW)*(ZP(gt3,2) + ZP(gt3,3)));

   return {result};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gP>::type, fields::gWmC, fields::VWm>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Sin(ThetaW));

   return {result, 1};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWm>::type, fields::gZ, fields::VWm>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Cos(ThetaW);

   return {result, 1};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gWmC, fields::VWm>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Cos(ThetaW));

   return {result, 1};
}

MomentumDifferenceVertex VertexImpl<fields::Ah, fields::hh, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt3385 = indices[0];
   const int gt3386 = indices[1];

   const std::complex<double> result = 0;

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, fields::VP, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*((7.745966692414834*g1*g2*Cos(2*ThetaW) + Sin(2*ThetaW)*(-3*Sqr(g1) + 5*Sqr(g2)))*ZP(gt1,0)*ZP(gt2,0) + (7.745966692414834*g1*g2*Cos(2*ThetaW) + Sin(2*ThetaW)*(-3*Sqr(g1) + 5*Sqr(g2)))*ZP(gt1,1)*ZP(gt2,1) + 20*Sin(2*ThetaW)*Sqr(g2)*(ZP(gt1,2)*ZP(gt2,2) + ZP(gt1,3)*ZP(gt2,3)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VP, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.016666666666666666*(-((7.745966692414834*g1*g2*Cos(2*ThetaW) + Sin(2*ThetaW)*(Sqr(g1) - 15*Sqr(g2)))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1))) - 4*Sin(2*ThetaW)*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Se, typename fields::conj<fields::Se>::type, fields::VP, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*((7.745966692414834*g1*g2*Cos(2*ThetaW) + Sin(2*ThetaW)*(-3*Sqr(g1) + 5*Sqr(g2)))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) - 12*Sin(2*ThetaW)*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SRdp, typename fields::conj<fields::SRdp>::type, fields::VP, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*(7.745966692414834*g1*g2*Cos(2*ThetaW) + Sin(2*ThetaW)*(-3*Sqr(g1) + 5*Sqr(g2)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SRum, typename fields::conj<fields::SRum>::type, fields::VP, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*(7.745966692414834*g1*g2*Cos(2*ThetaW) + Sin(2*ThetaW)*(-3*Sqr(g1) + 5*Sqr(g2)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VP, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.016666666666666666*((7.745966692414834*g1*g2*Cos(2*ThetaW) - Sin(2*ThetaW)*(Sqr(g1) - 15*Sqr(g2)))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1)) - 16*Sin(2*ThetaW)*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1)));

   return {result};
}

QuadrupleVectorVertex VertexImpl<typename fields::conj<fields::VWm>::type, fields::VP, fields::VWm, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> part1 = Cos(ThetaW)*Sin(ThetaW)*Sqr(g2);

   const std::complex<double> part2 = -(Sin(2*ThetaW)*Sqr(g2));

   const std::complex<double> part3 = Cos(ThetaW)*Sin(ThetaW)*Sqr(g2);

   return {part1, part2, part3};
}

InverseMetricVertex VertexImpl<fields::Ah, typename fields::conj<fields::Hpm>::type, fields::VP, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = std::complex<double>(0,-0.5)*g2*(0.7745966692414834*g1*Cos(ThetaW)*ZA(gt1,0)*ZP(gt2,0) + 0.7745966692414834*g1*Cos(ThetaW)*ZA(gt1,1)*ZP(gt2,1) + 1.4142135623730951*g2*Sin(ThetaW)*ZA(gt1,3)*(ZP(gt2,2) - ZP(gt2,3)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Ah, fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = std::complex<double>(0,0.5)*g2*(0.7745966692414834*g1*Cos(ThetaW)*ZA(gt1,0)*ZP(gt2,0) + 0.7745966692414834*g1*Cos(ThetaW)*ZA(gt1,1)*ZP(gt2,1) + 1.4142135623730951*g2*Sin(ThetaW)*ZA(gt1,3)*(ZP(gt2,2) - ZP(gt2,3)));

   return {result};
}

QuadrupleVectorVertex VertexImpl<typename fields::conj<fields::VWm>::type, fields::VWm, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> part1 = -2*Sqr(g2)*Sqr(Cos(ThetaW));

   const std::complex<double> part2 = Sqr(g2)*Sqr(Cos(ThetaW));

   const std::complex<double> part3 = Sqr(g2)*Sqr(Cos(ThetaW));

   return {part1, part2, part3};
}

InverseMetricVertex VertexImpl<fields::Ah, typename fields::conj<fields::Hpm>::type, fields::VWm, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = std::complex<double>(0,-0.5)*g2*(-0.7745966692414834*g1*Sin(ThetaW)*ZA(gt1,0)*ZP(gt2,0) - 0.7745966692414834*g1*Sin(ThetaW)*ZA(gt1,1)*ZP(gt2,1) + 1.4142135623730951*g2*Cos(ThetaW)*ZA(gt1,3)*(ZP(gt2,2) - ZP(gt2,3)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Ah, fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = std::complex<double>(0,-0.5)*g2*(0.7745966692414834*g1*Sin(ThetaW)*ZA(gt1,0)*ZP(gt2,0) + 0.7745966692414834*g1*Sin(ThetaW)*ZA(gt1,1)*ZP(gt2,1) + 1.4142135623730951*g2*Cos(ThetaW)*ZA(gt1,3)*(-ZP(gt2,2) + ZP(gt2,3)));

   return {result};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fv>::type, fields::Fv, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.5*KroneckerDelta(gt1,gt2)*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW));

   const std::complex<double> right = 0;

   return {left, right};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gP, fields::hh>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.025*(7.745966692414834*g1*g2*Cos(2*ThetaW) + Sin(2*ThetaW)*(3*Sqr(g1) - 5*Sqr(g2)))*(vd*ZH(gt3,0) + vu*ZH(gt3,1));

   return {result};
}

QuadrupleVectorVertex VertexImpl<typename fields::conj<fields::VWm>::type, fields::VP, fields::VP, fields::VWm>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> part1 = Sqr(g2)*Sqr(Sin(ThetaW));

   const std::complex<double> part2 = Sqr(g2)*Sqr(Sin(ThetaW));

   const std::complex<double> part3 = -2*Sqr(g2)*Sqr(Sin(ThetaW));

   return {part1, part2, part3};
}

QuadrupleVectorVertex VertexImpl<typename fields::conj<fields::VWm>::type, typename fields::conj<fields::VWm>::type, fields::VWm, fields::VWm>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> part1 = 2*Sqr(g2);

   const std::complex<double> part2 = -Sqr(g2);

   const std::complex<double> part3 = -Sqr(g2);

   return {part1, part2, part3};
}

InverseMetricVertex VertexImpl<fields::Hpm, fields::Hpm, typename fields::conj<fields::VWm>::type, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -2*Sqr(g2)*(ZP(gt1,3)*ZP(gt2,2) + ZP(gt1,2)*ZP(gt2,3));

   return {result};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gP, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Sin(ThetaW));

   return {result, 1};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWm>::type, fields::gP, fields::VWm>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Sin(ThetaW);

   return {result, 1};
}

InverseMetricVertex VertexImpl<typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Hpm>::type, fields::VWm, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -2*Sqr(g2)*(ZP(gt1,3)*ZP(gt2,2) + ZP(gt1,2)*ZP(gt2,3));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Sd, typename fields::conj<fields::SRum>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt4 = indices[2];
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,0.5)*SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt4,3 + j1)))*(1.4142135623730951*Conj(LamSU)*ZA(gt1,2) + Conj(LamTU)*ZA(gt1,3));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Su, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,-0.5)*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1)))*(1.4142135623730951*Conj(LamSD)*ZA(gt1,2) - Conj(LamTD)*ZA(gt1,3));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Sd, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,-0.25)*(2.8284271247461903*(-(SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZU(gt4,j3))*ZA(gt1,0)*ZP(gt3,0)) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gt4,j3))*ZA(gt1,1)*ZP(gt3,1) + SUM(j3,0,2,Conj(ZD(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yu(j2,j1))*ZU(gt4,3 + j2)))*(-(ZA(gt1,1)*ZP(gt3,0)) + ZA(gt1,0)*ZP(gt3,1))) + Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZU(gt4,j1))*(1.4142135623730951*ZA(gt1,0)*ZP(gt3,0) - 1.4142135623730951*ZA(gt1,1)*ZP(gt3,1) + 2*ZA(gt1,3)*(ZP(gt3,2) + ZP(gt3,3))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Se, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,-0.25)*(-2.8284271247461903*SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZV(gt4,j3))*ZA(gt1,0)*ZP(gt3,0) + Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZV(gt4,j1))*(1.4142135623730951*ZA(gt1,0)*ZP(gt3,0) - 1.4142135623730951*ZA(gt1,1)*ZP(gt3,1) + 2*ZA(gt1,3)*(ZP(gt3,2) + ZP(gt3,3))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Sv, typename fields::conj<fields::Se>::type, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,-0.5)*SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)))*(1.4142135623730951*Conj(LamSD)*ZA(gt1,2) - Conj(LamTD)*ZA(gt1,3));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Sd, fields::SRdp, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt4 = indices[2];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,0.5)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1)))*ZU(gt4,j2))*(1.4142135623730951*LamSD*ZA(gt1,2) - LamTD*ZA(gt1,3));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Se, fields::SRdp, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt4 = indices[2];
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,0.5)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1)))*ZV(gt4,j2))*(1.4142135623730951*LamSD*ZA(gt1,2) - LamTD*ZA(gt1,3));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::SRum, fields::Su, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt4 = indices[2];
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,-0.5)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt3,3 + j1)))*ZD(gt4,j2))*(1.4142135623730951*LamSU*ZA(gt1,2) + LamTU*ZA(gt1,3));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Hpm, fields::Su, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,0.25)*(2.8284271247461903*(-(SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gt4,j3))*ZA(gt1,0)*ZP(gt2,0)) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZD(gt4,j3))*ZA(gt1,1)*ZP(gt2,1) + SUM(j3,0,2,Conj(ZU(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yd(j2,j1))*ZD(gt4,3 + j2)))*(-(ZA(gt1,1)*ZP(gt2,0)) + ZA(gt1,0)*ZP(gt2,1))) + Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt3,j1))*ZD(gt4,j1))*(1.4142135623730951*ZA(gt1,0)*ZP(gt2,0) - 1.4142135623730951*ZA(gt1,1)*ZP(gt2,1) + 2*ZA(gt1,3)*(ZP(gt2,2) + ZP(gt2,3))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Hpm, fields::Sv, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,0.25)*(-2.8284271247461903*SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gt3,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZE(gt4,j3))*ZA(gt1,0)*ZP(gt2,0) + Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt3,j1))*ZE(gt4,j1))*(1.4142135623730951*ZA(gt1,0)*ZP(gt2,0) - 1.4142135623730951*ZA(gt1,1)*ZP(gt2,1) + 2*ZA(gt1,3)*(ZP(gt2,2) + ZP(gt2,3))));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::Su, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -(Conj(LamTD)*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt4,3 + j1)))*ZHR(gt3,0)*ZP(gt1,3));

   return {result};
}

ScalarVertex VertexImpl<fields::Rh, fields::Rh, typename fields::conj<fields::SRdp>::type, typename fields::conj<fields::SRum>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.5*Sqr(g2)*(ZHR(gt1,1)*ZHR(gt2,0) + ZHR(gt1,0)*ZHR(gt2,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Rh, fields::Sd, typename fields::conj<fields::SRum>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt4 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.5*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZU(gt4,j1))*ZHR(gt1,1);

   return {result};
}

ScalarVertex VertexImpl<fields::Rh, fields::Se, typename fields::conj<fields::SRum>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt4 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.5*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZV(gt4,j1))*ZHR(gt1,1);

   return {result};
}

InverseMetricVertex VertexImpl<fields::Rh, typename fields::conj<fields::SRum>::type, fields::VP, fields::VWm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.5477225575051661*g1*g2*Cos(ThetaW)*ZHR(gt1,1);

   return {result};
}

InverseMetricVertex VertexImpl<fields::Rh, typename fields::conj<fields::SRum>::type, fields::VWm, fields::VZ>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5477225575051661*g1*g2*Sin(ThetaW)*ZHR(gt1,1);

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::Rh, fields::Su, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = LamTU*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt3,3 + j1)))*ZD(gt4,j2))*ZHR(gt2,1)*ZP(gt1,2);

   return {result};
}

ScalarVertex VertexImpl<fields::SRum, fields::Su, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt3 = indices[1];
   const int gt4 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.5*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZD(gt4,j1))*ZHR(gt3,1);

   return {result};
}

ScalarVertex VertexImpl<fields::Rh, fields::Su, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.5*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZD(gt3,j1))*ZHR(gt1,0);

   return {result};
}

ScalarVertex VertexImpl<fields::Se, fields::Su, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = 0.25*(-(Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZV(gt4,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZD(gt3,j2))) - Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZV(gt4,j2)) - 4*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Ye(j3,j4))*Conj(ZE(gt1,3 + j3)))*ZV(gt4,j4)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Su, typename fields::conj<fields::Sd>::type, fields::VP, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.18257418583505536*g1*g2*Cos(ThetaW)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZD(gt2,j1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Su, typename fields::conj<fields::Sd>::type, fields::VWm, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.18257418583505536*g1*g2*Sin(ThetaW)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZD(gt2,j1));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::Sv, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -(Conj(LamTD)*SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt4,3 + j1)))*ZHR(gt3,0)*ZP(gt1,3));

   return {result};
}

ScalarVertex VertexImpl<fields::SRum, fields::Sv, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt3 = indices[1];
   const int gt4 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.5*Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZE(gt4,j1))*ZHR(gt3,1);

   return {result};
}

ScalarVertex VertexImpl<fields::Rh, fields::Sv, typename fields::conj<fields::Se>::type, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.5*Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZE(gt3,j1))*ZHR(gt1,0);

   return {result};
}

ScalarVertex VertexImpl<fields::Sd, fields::Sv, typename fields::conj<fields::Se>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);

   const std::complex<double> result = 0.25*(-(Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZU(gt4,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZE(gt3,j2))) - Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZU(gt4,j2)) - 4*SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yd(j3,j4))*Conj(ZD(gt1,3 + j3)))*ZU(gt4,j4)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sv, typename fields::conj<fields::Se>::type, fields::VP, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.5477225575051661*g1*g2*Cos(ThetaW)*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZE(gt2,j1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sv, typename fields::conj<fields::Se>::type, fields::VWm, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5477225575051661*g1*g2*Sin(ThetaW)*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZE(gt2,j1));

   return {result};
}

ScalarVertex VertexImpl<fields::Sd, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = Conj(LamTU)*SUM(j2,0,2,Conj(ZD(gt1,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt4,3 + j1)))*ZHR(gt3,1)*ZP(gt2,2);

   return {result};
}

ScalarVertex VertexImpl<fields::SRdp, fields::SRum, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::Rh>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt4 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.5*Sqr(g2)*(ZHR(gt3,1)*ZHR(gt4,0) + ZHR(gt3,0)*ZHR(gt4,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Sd, fields::SRdp, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt4 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.5*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZU(gt4,j1))*ZHR(gt3,0);

   return {result};
}

ScalarVertex VertexImpl<fields::Se, fields::SRdp, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt4 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.5*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZV(gt4,j1))*ZHR(gt3,0);

   return {result};
}

InverseMetricVertex VertexImpl<fields::SRdp, typename fields::conj<fields::Rh>::type, fields::VP, fields::VWm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5477225575051661*g1*g2*Cos(ThetaW)*ZHR(gt2,0);

   return {result};
}

InverseMetricVertex VertexImpl<fields::SRdp, typename fields::conj<fields::Rh>::type, fields::VWm, fields::VZ>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.5477225575051661*g1*g2*Sin(ThetaW)*ZHR(gt2,0);

   return {result};
}

ScalarVertex VertexImpl<fields::Rh, fields::Sd, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -(LamTD*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1)))*ZU(gt4,j2))*ZHR(gt1,0)*ZP(gt3,3));

   return {result};
}

ScalarVertex VertexImpl<fields::Rh, fields::Se, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -(LamTD*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1)))*ZV(gt4,j2))*ZHR(gt1,0)*ZP(gt3,3));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SRum, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::VWm>::type, fields::VP>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.5477225575051661*g1*g2*Cos(ThetaW)*ZHR(gt2,1);

   return {result};
}

InverseMetricVertex VertexImpl<fields::Rh, typename fields::conj<fields::SRdp>::type, typename fields::conj<fields::VWm>::type, fields::VP>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5477225575051661*g1*g2*Cos(ThetaW)*ZHR(gt1,0);

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sd, typename fields::conj<fields::Su>::type, typename fields::conj<fields::VWm>::type, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.18257418583505536*g1*g2*Cos(ThetaW)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZU(gt2,j1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Se, typename fields::conj<fields::Sv>::type, typename fields::conj<fields::VWm>::type, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.5477225575051661*g1*g2*Cos(ThetaW)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZV(gt2,j1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SRum, typename fields::conj<fields::Rh>::type, typename fields::conj<fields::VWm>::type, fields::VZ>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5477225575051661*g1*g2*Sin(ThetaW)*ZHR(gt2,1);

   return {result};
}

InverseMetricVertex VertexImpl<fields::Rh, typename fields::conj<fields::SRdp>::type, typename fields::conj<fields::VWm>::type, fields::VZ>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.5477225575051661*g1*g2*Sin(ThetaW)*ZHR(gt1,0);

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sd, typename fields::conj<fields::Su>::type, typename fields::conj<fields::VWm>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.18257418583505536*g1*g2*Sin(ThetaW)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZU(gt2,j1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Se, typename fields::conj<fields::Sv>::type, typename fields::conj<fields::VWm>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5477225575051661*g1*g2*Sin(ThetaW)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZV(gt2,j1));

   return {result};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Glu>::type, fields::Glu, fields::sigmaO>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g3 = MODELPARAMETER(g3);

   const std::complex<double> left = g3;

   const std::complex<double> right = -g3;

   return {left, right};
}

MomentumDifferenceVertex VertexImpl<fields::sigmaO, fields::sigmaO, fields::VG>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const auto g3 = MODELPARAMETER(g3);

   const std::complex<double> result = std::complex<double>(0,-1)*g3;

   return {result, minuend_index, subtrahend_index};
}

ScalarVertex VertexImpl<fields::phiO, fields::Sd, fields::sigmaO, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt4 = indices[1];
   const auto g3 = MODELPARAMETER(g3);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = -0.5*Sqr(g3)*(SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt4,j1)) - SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt4,3 + j1)));

   return {result};
}

ScalarVertex VertexImpl<fields::phiO, fields::sigmaO, fields::Su, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt4 = indices[1];
   const auto g3 = MODELPARAMETER(g3);
   const auto ZU = MODELPARAMETER(ZU);

   const std::complex<double> result = -0.5*Sqr(g3)*(SUM(j1,0,2,Conj(ZU(gt3,j1))*ZU(gt4,j1)) - SUM(j1,0,2,Conj(ZU(gt3,3 + j1))*ZU(gt4,3 + j1)));

   return {result};
}

ScalarVertex VertexImpl<fields::phiO, fields::phiO, fields::sigmaO, fields::sigmaO>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g3 = MODELPARAMETER(g3);

   const std::complex<double> result = -Sqr(g3);

   return {result};
}

InverseMetricVertex VertexImpl<fields::sigmaO, fields::sigmaO, fields::VG, fields::VG>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g3 = MODELPARAMETER(g3);

   const std::complex<double> result = Sqr(g3);

   return {result};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Glu>::type, fields::Glu, fields::phiO>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g3 = MODELPARAMETER(g3);

   const std::complex<double> left = std::complex<double>(0,1)*g3;

   const std::complex<double> right = std::complex<double>(0,1)*g3;

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Glu>::type, fields::Glu, fields::VG>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g3 = MODELPARAMETER(g3);

   const std::complex<double> left = std::complex<double>(0,1)*g3;

   const std::complex<double> right = std::complex<double>(0,1)*g3;

   return {left, right};
}

InverseMetricVertex VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VG, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g3 = MODELPARAMETER(g3);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.16666666666666666*g3*((0.7745966692414834*g1*Cos(ThetaW) - 3*g2*Sin(ThetaW))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1)) - 1.5491933384829668*g1*Cos(ThetaW)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VG, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g3 = MODELPARAMETER(g3);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.16666666666666666*g3*((0.7745966692414834*g1*Cos(ThetaW) + 3*g2*Sin(ThetaW))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1)) + 3.0983866769659336*g1*Cos(ThetaW)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VG, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g3 = MODELPARAMETER(g3);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.16666666666666666*g3*(3*g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1)) + 0.2581988897471611*g1*g3*Sin(ThetaW)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VG, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g3 = MODELPARAMETER(g3);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.16666666666666666*g3*((3*g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1)) - 3.0983866769659336*g1*Sin(ThetaW)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1)));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::phiO, fields::phiO, fields::VG>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const auto g3 = MODELPARAMETER(g3);

   const std::complex<double> result = std::complex<double>(0,-1)*g3;

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::phiO, fields::phiO, fields::VG, fields::VG>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g3 = MODELPARAMETER(g3);

   const std::complex<double> result = Sqr(g3);

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, fields::Ah, fields::hh>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = std::complex<double>(0,-0.5)*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(ZA(gt1,0)*ZH(gt2,0) - ZA(gt1,1)*ZH(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

ChiralVertex VertexImpl<fields::VZ, fields::Chi, typename fields::bar<fields::Chi>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZN2 = MODELPARAMETER(ZN2);
   const auto ZN1 = MODELPARAMETER(ZN1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.5*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(Conj(ZN2(gt1,2))*ZN2(gt2,2) - Conj(ZN2(gt1,3))*ZN2(gt2,3));

   const std::complex<double> right = -0.5*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(Conj(ZN1(gt2,2))*ZN1(gt1,2) - Conj(ZN1(gt2,3))*ZN1(gt1,3));

   return {left, right};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, fields::hh, fields::Ah>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = std::complex<double>(0,-0.5)*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(ZA(gt1,0)*ZH(gt2,0) - ZA(gt1,1)*ZH(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, fields::Hpm, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*((g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*ZP(gt1,0)*ZP(gt2,0) + (g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*ZP(gt1,1)*ZP(gt2,1) + 2*g2*Cos(ThetaW)*(ZP(gt1,2)*ZP(gt2,2) + ZP(gt1,3)*ZP(gt2,3)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, fields::Rh, typename fields::conj<fields::Rh>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(ZHR(gt1,0)*ZHR(gt2,0) - ZHR(gt1,1)*ZHR(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, fields::Sd, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.16666666666666666*(3*g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1)) - 0.2581988897471611*g1*Sin(ThetaW)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, fields::Se, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) - 0.7745966692414834*g1*Sin(ThetaW)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, fields::SRdp, typename fields::conj<fields::SRdp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.1*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, fields::SRum, typename fields::conj<fields::SRum>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, fields::Su, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.03333333333333333*((-15*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1)) + 15.491933384829668*g1*Sin(ThetaW)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, fields::Sv, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.5*KroneckerDelta(gt1,gt2)*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW));

   return {result, minuend_index, subtrahend_index};
}

ChiralVertex VertexImpl<fields::VZ, typename fields::bar<fields::Cha1>::type, fields::Cha1>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -(g2*Conj(UP1(gt2,0))*Cos(ThetaW)*UP1(gt1,0)) + 0.1*Conj(UP1(gt2,1))*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*UP1(gt1,1);

   const std::complex<double> right = -(g2*Conj(UM1(gt1,0))*Cos(ThetaW)*UM1(gt2,0)) + 0.1*Conj(UM1(gt1,1))*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*UM1(gt2,1);

   return {left, right};
}

ChiralVertex VertexImpl<fields::VZ, typename fields::bar<fields::Cha2>::type, fields::Cha2>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UM2 = MODELPARAMETER(UM2);
   const auto UP2 = MODELPARAMETER(UP2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = g2*Conj(UM2(gt2,0))*Cos(ThetaW)*UM2(gt1,0) + 0.5*Conj(UM2(gt2,1))*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*UM2(gt1,1);

   const std::complex<double> right = g2*Conj(UP2(gt1,0))*Cos(ThetaW)*UP2(gt2,0) + 0.5*Conj(UP2(gt1,1))*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*UP2(gt2,1);

   return {left, right};
}

ChiralVertex VertexImpl<fields::VZ, typename fields::bar<fields::Chi>::type, fields::Chi>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZN1 = MODELPARAMETER(ZN1);
   const auto ZN2 = MODELPARAMETER(ZN2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.5*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(Conj(ZN1(gt2,2))*ZN1(gt1,2) - Conj(ZN1(gt2,3))*ZN1(gt1,3));

   const std::complex<double> right = 0.5*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(Conj(ZN2(gt1,2))*ZN2(gt2,2) - Conj(ZN2(gt1,3))*ZN2(gt2,3));

   return {left, right};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, typename fields::conj<fields::Hpm>::type, fields::Hpm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*((g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*ZP(gt1,0)*ZP(gt2,0) + (g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*ZP(gt1,1)*ZP(gt2,1) + 2*g2*Cos(ThetaW)*(ZP(gt1,2)*ZP(gt2,2) + ZP(gt1,3)*ZP(gt2,3)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, typename fields::conj<fields::Rh>::type, fields::Rh>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZHR = MODELPARAMETER(ZHR);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(ZHR(gt1,0)*ZHR(gt2,0) - ZHR(gt1,1)*ZHR(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, typename fields::conj<fields::Sd>::type, fields::Sd>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.16666666666666666*(3*g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1)) - 0.2581988897471611*g1*Sin(ThetaW)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, typename fields::conj<fields::Se>::type, fields::Se>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) - 0.7745966692414834*g1*Sin(ThetaW)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, typename fields::conj<fields::SRdp>::type, fields::SRdp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.1*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, typename fields::conj<fields::SRum>::type, fields::SRum>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, typename fields::conj<fields::Su>::type, fields::Su>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.03333333333333333*((-15*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1)) + 15.491933384829668*g1*Sin(ThetaW)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, typename fields::conj<fields::Sv>::type, fields::Sv>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.5*KroneckerDelta(gt1,gt2)*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW));

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::VZ, fields::hh, fields::VP>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3394 = indices[0];

   const std::complex<double> result = 0;

   return {result};
}

InverseMetricVertex VertexImpl<fields::VZ, fields::hh, fields::VZ>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*Sqr(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(vd*ZH(gt1,0) + vu*ZH(gt1,1));

   return {result};
}

} // namespace detail
} // namespace MRSSM2_cxx_diagrams
} // namespace flexiblesusy
