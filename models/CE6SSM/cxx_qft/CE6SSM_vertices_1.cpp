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
 * @file cxx_qft/CE6SSM_vertices.cpp
 *
 * This file was generated with FlexibleSUSY 2.6.0 and SARAH 4.14.5 .
 */

#include "CE6SSM_context_base.hpp"
#include "CE6SSM_input_parameters.hpp"
#include "CE6SSM_vertices.hpp"

#include "concatenate.hpp"
#include "wrappers.hpp"

#define INPUTPARAMETER(p) context.model.get_input_parameters().p
#define MODELPARAMETER(p) context.model.get_##p()
#define DERIVEDPARAMETER(p) context.model.p()
#define PHASE(p) context.model.get_##p()

namespace flexiblesusy {
namespace CE6SSM_cxx_diagrams {
namespace detail {

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Cha, fields::Sv>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Cha>::type, fields::Cha, fields::VP>::evaluate(
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

ChiralVertex VertexImpl<typename fields::conj<fields::Sv>::type, typename fields::bar<fields::Cha>::type, fields::Fe>::evaluate(
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

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*SUM(j1,0,2,Conj(ZEL(gt2,j1))*Conj(ZER(gt1,j1))*Ye(j1,j1))*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*SUM(j1,0,2,Conj(Ye(j1,j1))*ZEL(gt1,j1)*ZER(gt2,j1))*ZA(gt3,0);

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

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*SUM(j1,0,2,Conj(ZEL(gt2,j1))*Conj(ZER(gt1,j1))*Ye(j1,j1))*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*SUM(j1,0,2,Conj(Ye(j1,j1))*ZEL(gt1,j1)*ZER(gt2,j1))*ZA(gt3,0);

   return {left, right};
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

   const std::complex<double> left = -0.7071067811865475*SUM(j1,0,2,Conj(ZEL(gt2,j1))*Conj(ZER(gt1,j1))*Ye(j1,j1))*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*SUM(j1,0,2,Conj(Ye(j1,j1))*ZEL(gt1,j1)*ZER(gt2,j1))*ZH(gt3,0);

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

   const std::complex<double> left = -0.7071067811865475*SUM(j1,0,2,Conj(ZEL(gt2,j1))*Conj(ZER(gt1,j1))*Ye(j1,j1))*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*SUM(j1,0,2,Conj(Ye(j1,j1))*ZEL(gt1,j1)*ZER(gt2,j1))*ZH(gt3,0);

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

   const std::complex<double> left = IF(gt2 < 3,Conj(ZER(gt1,gt2))*Ye(gt2,gt2)*ZP(gt3,0),0);

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

   const std::complex<double> result = 0.5*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*(ZP(gt1,0)*ZP(gt2,0) + ZP(gt1,1)*ZP(gt2,1));

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

   const std::complex<double> right = IF(gt1 < 3,Conj(Ye(gt1,gt1))*ZER(gt2,gt1)*ZP(gt3,0),0);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Se, fields::Chi>::evaluate(
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

ChiralVertex VertexImpl<fields::Chi, typename fields::conj<fields::Se>::type, fields::Fe>::evaluate(
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

ScalarVertex VertexImpl<fields::hh, fields::Sd, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
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

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto TLambdax = MODELPARAMETER(TLambdax);
   const auto gN = MODELPARAMETER(gN);
   const auto vd = MODELPARAMETER(vd);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto vu = MODELPARAMETER(vu);
   const auto vs = MODELPARAMETER(vs);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.025*(-5*ZA(gt1,2)*(2.8284271247461903*Conj(TLambdax)*(ZA(gt2,1)*ZH(gt3,0) + ZA(gt2,0)*ZH(gt3,1)) + 2.8284271247461903*TLambdax*(ZA(gt2,1)*ZH(gt3,0) + ZA(gt2,0)*ZH(gt3,1)) + ZA(gt2,2)*((8*vd*AbsSqr(Lambdax) - 3*vd*Sqr(gN))*ZH(gt3,0) - 2*vu*(-4*AbsSqr(Lambdax) + Sqr(gN))*ZH(gt3,1) + 5*vs*Sqr(gN)*ZH(gt3,2))) + 2*ZA(gt1,1)*(ZA(gt2,1)*(vd*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt3,0) - vu*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZH(gt3,1) + 5*vs*(-4*AbsSqr(Lambdax) + Sqr(gN))*ZH(gt3,2)) - 7.0710678118654755*(Conj(TLambdax) + TLambdax)*(ZA(gt2,2)*ZH(gt3,0) + ZA(gt2,0)*ZH(gt3,2))) - ZA(gt1,0)*(ZA(gt2,0)*(vd*(6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*ZH(gt3,0) + 2*vu*(20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2) + 3*Sqr(gN))*ZH(gt3,1) + 5*vs*(8*AbsSqr(Lambdax) - 3*Sqr(gN))*ZH(gt3,2)) + 14.142135623730951*(Conj(TLambdax) + TLambdax)*(ZA(gt2,2)*ZH(gt3,1) + ZA(gt2,1)*ZH(gt3,2))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Sd, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
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

ScalarVertex VertexImpl<fields::Ah, fields::hh, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto TLambdax = MODELPARAMETER(TLambdax);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = std::complex<double>(0.,-0.35355339059327373)*(Conj(TLambdax) - TLambdax)*(ZA(gt1,2)*(ZH(gt2,1)*ZH(gt3,0) + ZH(gt2,0)*ZH(gt3,1)) + ZA(gt1,1)*(ZH(gt2,2)*ZH(gt3,0) + ZH(gt2,0)*ZH(gt3,2)) + ZA(gt1,0)*(ZH(gt2,2)*ZH(gt3,1) + ZH(gt2,1)*ZH(gt3,2)));

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
   const auto gN = MODELPARAMETER(gN);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = std::complex<double>(0,-0.05)*((10*g2*Cos(ThetaW)*Cos(ThetaWp) + 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp))*ZA(gt1,0)*ZH(gt2,0) - 2*(5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*ZA(gt1,1)*ZH(gt2,1) + 15.811388300841898*gN*Sin(ThetaWp)*ZA(gt1,2)*ZH(gt2,2));

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
   const auto gN = MODELPARAMETER(gN);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.016666666666666666*((30*g2*Cos(ThetaW)*Cos(ThetaWp) + 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1)) + 2*(-7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::Ah, fields::hh, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = std::complex<double>(0,0.05)*((9.486832980505138*gN*Cos(ThetaWp) + 2*(5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*ZA(gt1,0)*ZH(gt2,0) + 2*(3.1622776601683795*gN*Cos(ThetaWp) - (5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*ZA(gt1,1)*ZH(gt2,1) - 15.811388300841898*gN*Cos(ThetaWp)*ZA(gt1,2)*ZH(gt2,2));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.016666666666666666*(-((9.486832980505138*gN*Cos(ThetaWp) + 2*(15*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1))) + 2*(9.486832980505138*gN*Cos(ThetaWp) + 7.745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1)));

   return {result, minuend_index, subtrahend_index};
}

ChiralVertex VertexImpl<fields::Chi, fields::Chi, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = 0.05*(-10*g2*Conj(ZN(gt1,1))*Conj(ZN(gt2,2))*ZH(gt3,0) + 9.486832980505138*gN*Conj(ZN(gt1,5))*Conj(ZN(gt2,2))*ZH(gt3,0) + 14.142135623730951*Conj(ZN(gt1,4))*Conj(ZN(gt2,3))*Lambdax*ZH(gt3,0) + 14.142135623730951*Conj(ZN(gt1,3))*Conj(ZN(gt2,4))*Lambdax*ZH(gt3,0) - 7.745966692414834*g1*Conj(ZN(gt1,3))*Conj(ZN(gt2,0))*ZH(gt3,1) + 10*g2*Conj(ZN(gt1,3))*Conj(ZN(gt2,1))*ZH(gt3,1) + 10*g2*Conj(ZN(gt1,1))*Conj(ZN(gt2,3))*ZH(gt3,1) + 6.324555320336759*gN*Conj(ZN(gt1,5))*Conj(ZN(gt2,3))*ZH(gt3,1) + 6.324555320336759*gN*Conj(ZN(gt1,3))*Conj(ZN(gt2,5))*ZH(gt3,1) + 14.142135623730951*Conj(ZN(gt1,4))*Conj(ZN(gt2,2))*Lambdax*ZH(gt3,1) + 7.745966692414834*g1*Conj(ZN(gt1,0))*(Conj(ZN(gt2,2))*ZH(gt3,0) - Conj(ZN(gt2,3))*ZH(gt3,1)) - 15.811388300841898*gN*Conj(ZN(gt1,5))*Conj(ZN(gt2,4))*ZH(gt3,2) - 15.811388300841898*gN*Conj(ZN(gt1,4))*Conj(ZN(gt2,5))*ZH(gt3,2) + 14.142135623730951*Conj(ZN(gt1,3))*Conj(ZN(gt2,2))*Lambdax*ZH(gt3,2) + Conj(ZN(gt1,2))*(7.745966692414834*g1*Conj(ZN(gt2,0))*ZH(gt3,0) - 10*g2*Conj(ZN(gt2,1))*ZH(gt3,0) + 9.486832980505138*gN*Conj(ZN(gt2,5))*ZH(gt3,0) + 14.142135623730951*Lambdax*(Conj(ZN(gt2,4))*ZH(gt3,1) + Conj(ZN(gt2,3))*ZH(gt3,2))));

   const std::complex<double> right = 0.05*(2*ZH(gt3,1)*((-3.872983346207417*g1*ZN(gt1,0) + 5*g2*ZN(gt1,1) + 3.1622776601683795*gN*ZN(gt1,5))*ZN(gt2,3) + 7.0710678118654755*Conj(Lambdax)*(ZN(gt1,4)*ZN(gt2,2) + ZN(gt1,2)*ZN(gt2,4)) + ZN(gt1,3)*(-3.872983346207417*g1*ZN(gt2,0) + 5*g2*ZN(gt2,1) + 3.1622776601683795*gN*ZN(gt2,5))) + ZH(gt3,0)*(7.745966692414834*g1*ZN(gt1,0)*ZN(gt2,2) - 10*g2*ZN(gt1,1)*ZN(gt2,2) + 9.486832980505138*gN*ZN(gt1,5)*ZN(gt2,2) + 14.142135623730951*Conj(Lambdax)*ZN(gt1,4)*ZN(gt2,3) + 14.142135623730951*Conj(Lambdax)*ZN(gt1,3)*ZN(gt2,4) + ZN(gt1,2)*(7.745966692414834*g1*ZN(gt2,0) - 10*g2*ZN(gt2,1) + 9.486832980505138*gN*ZN(gt2,5))) + 7.0710678118654755*ZH(gt3,2)*(2*Conj(Lambdax)*(ZN(gt1,3)*ZN(gt2,2) + ZN(gt1,2)*ZN(gt2,3)) - 2.23606797749979*gN*(ZN(gt1,5)*ZN(gt2,4) + ZN(gt1,4)*ZN(gt2,5))));

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
   const auto gN = MODELPARAMETER(gN);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZDR = MODELPARAMETER(ZDR);

   const std::complex<double> left = -0.18257418583505536*g1*Conj(ZN(gt1,0))*SUM(j1,0,2,Conj(ZDL(gt2,j1))*ZD(gt3,j1)) + 0.7071067811865475*g2*Conj(ZN(gt1,1))*SUM(j1,0,2,Conj(ZDL(gt2,j1))*ZD(gt3,j1)) - 0.22360679774997896*gN*Conj(ZN(gt1,5))*SUM(j1,0,2,Conj(ZDL(gt2,j1))*ZD(gt3,j1)) - Conj(ZN(gt1,2))*SUM(j1,0,2,Conj(ZDL(gt2,j1))*Yd(j1,j1)*ZD(gt3,3 + j1));

   const std::complex<double> right = -(SUM(j1,0,2,Conj(Yd(j1,j1))*ZD(gt3,j1)*ZDR(gt2,j1))*ZN(gt1,2)) - 0.14907119849998596*SUM(j1,0,2,ZD(gt3,3 + j1)*ZDR(gt2,j1))*(2.449489742783178*g1*ZN(gt1,0) + 3*gN*ZN(gt1,5));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Chi, fields::Sd>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZDR = MODELPARAMETER(ZDR);
   const auto ZDL = MODELPARAMETER(ZDL);

   const std::complex<double> left = -0.3651483716701107*g1*Conj(ZN(gt2,0))*SUM(j1,0,2,Conj(ZD(gt3,3 + j1))*Conj(ZDR(gt1,j1))) - 0.4472135954999579*gN*Conj(ZN(gt2,5))*SUM(j1,0,2,Conj(ZD(gt3,3 + j1))*Conj(ZDR(gt1,j1))) - Conj(ZN(gt2,2))*SUM(j1,0,2,Conj(ZD(gt3,j1))*Conj(ZDR(gt1,j1))*Yd(j1,j1));

   const std::complex<double> right = -(SUM(j1,0,2,Conj(Yd(j1,j1))*Conj(ZD(gt3,3 + j1))*ZDL(gt1,j1))*ZN(gt2,2)) - 0.03333333333333333*SUM(j1,0,2,Conj(ZD(gt3,j1))*ZDL(gt1,j1))*(5.477225575051661*g1*ZN(gt2,0) - 21.213203435596427*g2*ZN(gt2,1) + 6.708203932499369*gN*ZN(gt2,5));

   return {left, right};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto TLambdax = MODELPARAMETER(TLambdax);
   const auto gN = MODELPARAMETER(gN);
   const auto vd = MODELPARAMETER(vd);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto vs = MODELPARAMETER(vs);
   const auto vu = MODELPARAMETER(vu);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.025*(5*ZH(gt1,2)*(2.8284271247461903*Conj(TLambdax)*ZH(gt2,1)*ZH(gt3,0) + 2.8284271247461903*TLambdax*ZH(gt2,1)*ZH(gt3,0) - 8*vd*AbsSqr(Lambdax)*ZH(gt2,2)*ZH(gt3,0) + 3*vd*Sqr(gN)*ZH(gt2,2)*ZH(gt3,0) - 8*vs*AbsSqr(Lambdax)*ZH(gt2,1)*ZH(gt3,1) + 2*vs*Sqr(gN)*ZH(gt2,1)*ZH(gt3,1) - 8*vu*AbsSqr(Lambdax)*ZH(gt2,2)*ZH(gt3,1) + 2*vu*Sqr(gN)*ZH(gt2,2)*ZH(gt3,1) - 8*vu*AbsSqr(Lambdax)*ZH(gt2,1)*ZH(gt3,2) + 2*vu*Sqr(gN)*ZH(gt2,1)*ZH(gt3,2) - 15*vs*Sqr(gN)*ZH(gt2,2)*ZH(gt3,2) + ZH(gt2,0)*(vs*(-8*AbsSqr(Lambdax) + 3*Sqr(gN))*ZH(gt3,0) + 2.8284271247461903*Conj(TLambdax)*ZH(gt3,1) + 2.8284271247461903*TLambdax*ZH(gt3,1) - 8*vd*AbsSqr(Lambdax)*ZH(gt3,2) + 3*vd*Sqr(gN)*ZH(gt3,2))) + ZH(gt1,0)*(-(ZH(gt2,0)*(3*vd*(6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*ZH(gt3,0) + 2*vu*(20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2) + 3*Sqr(gN))*ZH(gt3,1) + 5*vs*(8*AbsSqr(Lambdax) - 3*Sqr(gN))*ZH(gt3,2))) + 5*ZH(gt2,2)*(vs*(-8*AbsSqr(Lambdax) + 3*Sqr(gN))*ZH(gt3,0) + 2.8284271247461903*Conj(TLambdax)*ZH(gt3,1) + 2.8284271247461903*TLambdax*ZH(gt3,1) - 8*vd*AbsSqr(Lambdax)*ZH(gt3,2) + 3*vd*Sqr(gN)*ZH(gt3,2)) + 2*ZH(gt2,1)*(vu*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt3,0) + vd*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt3,1) + 7.0710678118654755*(Conj(TLambdax) + TLambdax)*ZH(gt3,2))) + 2*ZH(gt1,1)*(ZH(gt2,1)*(vd*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt3,0) - 3*vu*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZH(gt3,1) + 5*vs*(-4*AbsSqr(Lambdax) + Sqr(gN))*ZH(gt3,2)) + ZH(gt2,0)*(vu*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt3,0) + vd*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt3,1) + 7.0710678118654755*(Conj(TLambdax) + TLambdax)*ZH(gt3,2)) + 5*ZH(gt2,2)*(1.4142135623730951*Conj(TLambdax)*ZH(gt3,0) + 1.4142135623730951*TLambdax*ZH(gt3,0) + (-4*AbsSqr(Lambdax) + Sqr(gN))*(vs*ZH(gt3,1) + vu*ZH(gt3,2)))));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto vd = MODELPARAMETER(vd);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto vu = MODELPARAMETER(vu);
   const auto vs = MODELPARAMETER(vs);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05*(vd*(-14.696938456699067*g1*gN*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) + 10*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Cos(ThetaWp)) + Cos(ThetaW)*(-18.973665961010276*g2*gN*Cos(ThetaWp)*Sin(ThetaWp) + 15.491933384829668*g1*g2*Sin(ThetaW)*Sqr(Cos(ThetaWp))) + 6*Sqr(g1)*Sqr(Cos(ThetaWp))*Sqr(Sin(ThetaW)) + 9*Sqr(gN)*Sqr(Sin(ThetaWp)))*ZH(gt1,0) + 2*vu*(3.1622776601683795*g2*gN*Cos(ThetaW)*Sin(2*ThetaWp) + g1*Sin(ThetaW)*(7.745966692414834*g2*Cos(ThetaW) + 3*g1*Sin(ThetaW))*Sqr(Cos(ThetaWp)) + 5*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Cos(ThetaWp)) + gN*(2.449489742783178*g1*Sin(ThetaW)*Sin(2*ThetaWp) + 2*gN*Sqr(Sin(ThetaWp))))*ZH(gt1,1) + 25*vs*Sqr(gN)*Sqr(Sin(ThetaWp))*ZH(gt1,2));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::VZ, fields::VZp>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto vd = MODELPARAMETER(vd);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto vs = MODELPARAMETER(vs);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05*(-(vd*(9.486832980505138*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) - 9*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN) + 5*Sin(2*ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 7.348469228349534*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + g1*(3.872983346207417*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) + 3*g1*Sin(2*ThetaWp)*Sqr(Sin(ThetaW)) - 7.348469228349534*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*ZH(gt1,0)) + vu*(6.324555320336759*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) + 4*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN) - 5*Sin(2*ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 4.898979485566356*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) - g1*(3.872983346207417*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) + 3*g1*Sin(2*ThetaWp)*Sqr(Sin(ThetaW)) + 4.898979485566356*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*ZH(gt1,1) + 25*vs*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN)*ZH(gt1,2));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::VZp, fields::VZp>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto vd = MODELPARAMETER(vd);
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto vs = MODELPARAMETER(vs);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*(vd*(9*Sqr(gN)*Sqr(Cos(ThetaWp)) + (g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(9.486832980505138*gN*Sin(2*ThetaWp) + 10*g2*Cos(ThetaW)*Sqr(Sin(ThetaWp)) + 7.745966692414834*g1*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*ZH(gt1,0) + 2*vu*(-3.1622776601683795*gN*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*Sin(2*ThetaWp) + 2*Sqr(gN)*Sqr(Cos(ThetaWp)) + 5*Sqr(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*Sqr(Sin(ThetaWp)))*ZH(gt1,1) + 25*vs*Sqr(gN)*Sqr(Cos(ThetaWp))*ZH(gt1,2));

   return {result};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha>::type, fields::Cha, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto UM = MODELPARAMETER(UM);
   const auto UP = MODELPARAMETER(UP);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*(g2*Conj(UM(gt2,0))*Conj(UP(gt1,1))*ZH(gt3,1) + Conj(UM(gt2,1))*(g2*Conj(UP(gt1,0))*ZH(gt3,0) + Conj(UP(gt1,1))*Lambdax*ZH(gt3,2)));

   const std::complex<double> right = -0.7071067811865475*(g2*UM(gt1,0)*UP(gt2,1)*ZH(gt3,1) + UM(gt1,1)*(g2*UP(gt2,0)*ZH(gt3,0) + Conj(Lambdax)*UP(gt2,1)*ZH(gt3,2)));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Cha, fields::Fu, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Yu = MODELPARAMETER(Yu);
   const auto UM = MODELPARAMETER(UM);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZUR = MODELPARAMETER(ZUR);
   const auto UP = MODELPARAMETER(UP);

   const std::complex<double> left = -(g2*Conj(UM(gt1,0))*SUM(j1,0,2,Conj(ZUL(gt2,j1))*ZD(gt3,j1))) + Conj(UM(gt1,1))*SUM(j1,0,2,Conj(ZUL(gt2,j1))*Yd(j1,j1)*ZD(gt3,3 + j1));

   const std::complex<double> right = SUM(j1,0,2,Conj(Yu(j1,j1))*ZD(gt3,j1)*ZUR(gt2,j1))*UP(gt1,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha>::type, typename fields::bar<fields::Fu>::type, fields::Sd>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto UP = MODELPARAMETER(UP);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZUR = MODELPARAMETER(ZUR);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto UM = MODELPARAMETER(UM);

   const std::complex<double> left = Conj(UP(gt1,1))*SUM(j1,0,2,Conj(ZD(gt3,j1))*Conj(ZUR(gt2,j1))*Yu(j1,j1));

   const std::complex<double> right = -(g2*SUM(j1,0,2,Conj(ZD(gt3,j1))*ZUL(gt2,j1))*UM(gt1,0)) + SUM(j1,0,2,Conj(Yd(j1,j1))*Conj(ZD(gt3,3 + j1))*ZUL(gt2,j1))*UM(gt1,1);

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

   const std::complex<double> left = -0.7071067811865475*SUM(j1,0,2,Conj(ZDL(gt2,j1))*Conj(ZDR(gt1,j1))*Yd(j1,j1))*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*SUM(j1,0,2,Conj(Yd(j1,j1))*ZDL(gt1,j1)*ZDR(gt2,j1))*ZH(gt3,0);

   return {left, right};
}

ChiralVertex VertexImpl<fields::Glu, fields::Fd, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt3 = indices[1];
   const auto g3 = MODELPARAMETER(g3);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZDR = MODELPARAMETER(ZDR);
   const auto PhaseGlu = PHASE(PhaseGlu);

   const std::complex<double> left = -0.7071067811865475*g3*PhaseGlu*SUM(j1,0,2,Conj(ZDL(gt2,j1))*ZD(gt3,j1));

   const std::complex<double> right = 0.7071067811865475*g3*Conj(PhaseGlu)*SUM(j1,0,2,ZD(gt3,3 + j1)*ZDR(gt2,j1));

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
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto PhaseGlu = PHASE(PhaseGlu);

   const std::complex<double> left = 0.7071067811865475*g3*PhaseGlu*SUM(j1,0,2,Conj(ZD(gt3,3 + j1))*Conj(ZDR(gt1,j1)));

   const std::complex<double> right = -0.7071067811865475*g3*Conj(PhaseGlu)*SUM(j1,0,2,Conj(ZD(gt3,j1))*ZDL(gt1,j1));

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

   const std::complex<double> left = -0.7071067811865475*SUM(j1,0,2,Conj(ZUL(gt2,j1))*Conj(ZUR(gt1,j1))*Yu(j1,j1))*ZH(gt3,1);

   const std::complex<double> right = -0.7071067811865475*SUM(j1,0,2,Conj(Yu(j1,j1))*ZUL(gt1,j1)*ZUR(gt2,j1))*ZH(gt3,1);

   return {left, right};
}

ScalarVertex VertexImpl<fields::hh, fields::Hpm, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
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

ScalarVertex VertexImpl<fields::Hpm, fields::Su, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto TYd = MODELPARAMETER(TYd);
   const auto TYu = MODELPARAMETER(TYu);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto vs = MODELPARAMETER(vs);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-1.4142135623730951*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZD(gt3,j1))*(vd*ZP(gt1,0) + vu*ZP(gt1,1)) + 2*(1.4142135623730951*vd*SUM(j1,0,2,AbsSqr(Yd(j1,j1))*Conj(ZU(gt2,j1))*ZD(gt3,j1))*ZP(gt1,0) + 1.4142135623730951*vs*Lambdax*SUM(j1,0,2,Conj(Yu(j1,j1))*Conj(ZU(gt2,3 + j1))*ZD(gt3,j1))*ZP(gt1,0) + 1.4142135623730951*vu*SUM(j1,0,2,Conj(Yu(j1,j1))*Conj(ZU(gt2,3 + j1))*Yd(j1,j1)*ZD(gt3,3 + j1))*ZP(gt1,0) + 2*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZD(gt3,3 + j1)*TYd(j1,j1))*ZP(gt1,0) + 1.4142135623730951*vu*SUM(j1,0,2,AbsSqr(Yu(j1,j1))*Conj(ZU(gt2,j1))*ZD(gt3,j1))*ZP(gt1,1) + 2*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j1))*ZD(gt3,j1))*ZP(gt1,1) + 1.4142135623730951*vs*Conj(Lambdax)*SUM(j1,0,2,Conj(ZU(gt2,j1))*Yd(j1,j1)*ZD(gt3,3 + j1))*ZP(gt1,1) + 1.4142135623730951*vd*SUM(j1,0,2,Conj(Yu(j1,j1))*Conj(ZU(gt2,3 + j1))*Yd(j1,j1)*ZD(gt3,3 + j1))*ZP(gt1,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Sd, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto TYd = MODELPARAMETER(TYd);
   const auto TYu = MODELPARAMETER(TYu);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto vs = MODELPARAMETER(vs);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-1.4142135623730951*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZU(gt3,j1))*(vd*ZP(gt2,0) + vu*ZP(gt2,1)) + 2*(1.4142135623730951*vd*SUM(j1,0,2,AbsSqr(Yd(j1,j1))*Conj(ZD(gt1,j1))*ZU(gt3,j1))*ZP(gt2,0) + 2*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*Conj(TYd(j1,j1))*ZU(gt3,j1))*ZP(gt2,0) + 1.4142135623730951*vs*Conj(Lambdax)*SUM(j1,0,2,Conj(ZD(gt1,j1))*Yu(j1,j1)*ZU(gt3,3 + j1))*ZP(gt2,0) + 1.4142135623730951*vu*SUM(j1,0,2,Conj(Yd(j1,j1))*Conj(ZD(gt1,3 + j1))*Yu(j1,j1)*ZU(gt3,3 + j1))*ZP(gt2,0) + 1.4142135623730951*vu*SUM(j1,0,2,AbsSqr(Yu(j1,j1))*Conj(ZD(gt1,j1))*ZU(gt3,j1))*ZP(gt2,1) + 1.4142135623730951*vs*Lambdax*SUM(j1,0,2,Conj(Yd(j1,j1))*Conj(ZD(gt1,3 + j1))*ZU(gt3,j1))*ZP(gt2,1) + 1.4142135623730951*vd*SUM(j1,0,2,Conj(Yd(j1,j1))*Conj(ZD(gt1,3 + j1))*Yu(j1,j1)*ZU(gt3,3 + j1))*ZP(gt2,1) + 2*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZU(gt3,3 + j1)*TYu(j1,j1))*ZP(gt2,1)));

   return {result};
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

   const std::complex<double> result = -0.5*g2*(ZH(gt1,0)*ZP(gt2,0) - ZH(gt1,1)*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
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

ScalarVertex VertexImpl<fields::hh, fields::Su, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
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

   const std::complex<double> result = 0.5*g2*(ZH(gt1,0)*ZP(gt2,0) - ZH(gt1,1)*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::hh, typename fields::conj<fields::VWm>::type, fields::VWm>::evaluate(
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

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Sd, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = 0.025*(SUM(j1,0,2,Conj(ZD(gt3,j1))*ZD(gt4,j1))*((2*Sqr(g1) + 10*Sqr(g2) + 3*Sqr(gN))*ZA(gt1,0)*ZA(gt2,0) - 2*(Sqr(g1) + 5*Sqr(g2) - Sqr(gN))*ZA(gt1,1)*ZA(gt2,1) - 5*Sqr(gN)*ZA(gt1,2)*ZA(gt2,2)) - 2*(20*SUM(j1,0,2,AbsSqr(Yd(j1,j1))*Conj(ZD(gt3,j1))*ZD(gt4,j1))*ZA(gt1,0)*ZA(gt2,0) + SUM(j1,0,2,Conj(ZD(gt3,3 + j1))*ZD(gt4,3 + j1))*(-((2*Sqr(g1) + 3*Sqr(gN))*ZA(gt1,0)*ZA(gt2,0)) + 2*(Sqr(g1) - Sqr(gN))*ZA(gt1,1)*ZA(gt2,1) + 5*Sqr(gN)*ZA(gt1,2)*ZA(gt2,2)) + 10*(2*SUM(j1,0,2,AbsSqr(Yd(j1,j1))*Conj(ZD(gt3,3 + j1))*ZD(gt4,3 + j1))*ZA(gt1,0)*ZA(gt2,0) + (Lambdax*SUM(j1,0,2,Conj(Yd(j1,j1))*Conj(ZD(gt3,3 + j1))*ZD(gt4,j1)) + Conj(Lambdax)*SUM(j1,0,2,Conj(ZD(gt3,j1))*Yd(j1,j1)*ZD(gt4,3 + j1)))*(ZA(gt1,2)*ZA(gt2,1) + ZA(gt1,1)*ZA(gt2,2)))));

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
   const auto gN = MODELPARAMETER(gN);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.008333333333333333*(SUM(j1,0,2,Conj(ZD(gt3,j1))*ZD(gt4,j1))*((6*Sqr(g1) + 30*Sqr(g2) + 9*Sqr(gN))*ZH(gt1,0)*ZH(gt2,0) - 6*(Sqr(g1) + 5*Sqr(g2) - Sqr(gN))*ZH(gt1,1)*ZH(gt2,1) - 15*Sqr(gN)*ZH(gt1,2)*ZH(gt2,2)) + 2*(-60*SUM(j1,0,2,AbsSqr(Yd(j1,j1))*Conj(ZD(gt3,j1))*ZD(gt4,j1))*ZH(gt1,0)*ZH(gt2,0) + SUM(j1,0,2,Conj(ZD(gt3,3 + j1))*ZD(gt4,3 + j1))*(3*(2*Sqr(g1) + 3*Sqr(gN))*ZH(gt1,0)*ZH(gt2,0) + 6*(-Sqr(g1) + Sqr(gN))*ZH(gt1,1)*ZH(gt2,1) - 15*Sqr(gN)*ZH(gt1,2)*ZH(gt2,2)) + 30*(-2*SUM(j1,0,2,AbsSqr(Yd(j1,j1))*Conj(ZD(gt3,3 + j1))*ZD(gt4,3 + j1))*ZH(gt1,0)*ZH(gt2,0) + (Lambdax*SUM(j1,0,2,Conj(Yd(j1,j1))*Conj(ZD(gt3,3 + j1))*ZD(gt4,j1)) + Conj(Lambdax)*SUM(j1,0,2,Conj(ZD(gt3,j1))*Yd(j1,j1)*ZD(gt4,3 + j1)))*(ZH(gt1,2)*ZH(gt2,1) + ZH(gt1,1)*ZH(gt2,2)))));

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
   const auto gN = MODELPARAMETER(gN);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.025*(SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt4,3 + j1))*(2*(2*Sqr(g1) + 3*Sqr(gN))*ZP(gt1,0)*ZP(gt3,0) + 4*(-Sqr(g1) + Sqr(gN))*ZP(gt1,1)*ZP(gt3,1)) + SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt4,j1))*((2*Sqr(g1) - 10*Sqr(g2) + 3*Sqr(gN))*ZP(gt1,0)*ZP(gt3,0) + 2*(-Sqr(g1) + 5*Sqr(g2) + Sqr(gN))*ZP(gt1,1)*ZP(gt3,1)) - 40*(SUM(j1,0,2,AbsSqr(Yd(j1,j1))*Conj(ZD(gt2,3 + j1))*ZD(gt4,3 + j1))*ZP(gt1,0)*ZP(gt3,0) + SUM(j1,0,2,AbsSqr(Yu(j1,j1))*Conj(ZD(gt2,j1))*ZD(gt4,j1))*ZP(gt1,1)*ZP(gt3,1)));

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
   const auto gN = MODELPARAMETER(gN);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = 0.004166666666666667*(-2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt2,j2))*ZD(gt3,j2)) - 30*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt2,j2))*ZD(gt3,j2)) + 20*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt2,j2))*ZD(gt3,j2)) - 3*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt2,j2))*ZD(gt3,j2)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt4,3 + j1))*SUM(j2,0,2,Conj(ZD(gt2,j2))*ZD(gt3,j2)) - 20*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt4,3 + j1))*SUM(j2,0,2,Conj(ZD(gt2,j2))*ZD(gt3,j2)) - 6*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt4,3 + j1))*SUM(j2,0,2,Conj(ZD(gt2,j2))*ZD(gt3,j2)) - 240*SUM(j1,0,2,Conj(ZD(gt1,j1))*Yd(j1,j1)*ZD(gt4,3 + j1))*SUM(j2,0,2,Conj(Yd(j2,j2))*Conj(ZD(gt2,3 + j2))*ZD(gt3,j2)) - 60*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt4,j1))*(SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) - SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2))) + 60*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt4,3 + j1))*(SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) - SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2))) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt2,3 + j2))*ZD(gt3,3 + j2)) - 20*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt2,3 + j2))*ZD(gt3,3 + j2)) - 6*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt2,3 + j2))*ZD(gt3,3 + j2)) - 8*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt4,3 + j1))*SUM(j2,0,2,Conj(ZD(gt2,3 + j2))*ZD(gt3,3 + j2)) + 20*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt4,3 + j1))*SUM(j2,0,2,Conj(ZD(gt2,3 + j2))*ZD(gt3,3 + j2)) - 12*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt4,3 + j1))*SUM(j2,0,2,Conj(ZD(gt2,3 + j2))*ZD(gt3,3 + j2)) - 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt4,j2)) - 30*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt4,j2)) + 20*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt4,j2)) - 3*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt4,j2)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt4,j2)) - 20*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt4,j2)) - 6*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt4,j2)) - 240*SUM(j1,0,2,Conj(ZD(gt2,j1))*Yd(j1,j1)*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(Yd(j2,j2))*Conj(ZD(gt1,3 + j2))*ZD(gt4,j2)) - 60*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt2,j2))*ZD(gt4,j2)) + 60*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZD(gt2,j2))*ZD(gt4,j2)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt4,3 + j2)) - 20*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt4,3 + j2)) - 6*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt4,3 + j2)) - 8*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt4,3 + j2)) + 20*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt4,3 + j2)) - 12*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt4,3 + j2)) + 60*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt2,3 + j2))*ZD(gt4,3 + j2)) - 60*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZD(gt2,3 + j2))*ZD(gt4,3 + j2)));

   return {result};
}

ScalarVertex VertexImpl<fields::Sd, fields::SDX, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::SDX>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g3 = MODELPARAMETER(g3);
   const auto ZDX = MODELPARAMETER(ZDX);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = -0.25*Sqr(g3)*(SUM(j1,0,2,Conj(ZDX(gt2,j1))*ZDX(gt4,j1))*(SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) - SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2))) + SUM(j1,0,2,Conj(ZDX(gt2,3 + j1))*ZDX(gt4,3 + j1))*(-SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) + SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2))) + (SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1)) - SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1)))*(SUM(j2,0,2,Conj(ZDX(gt2,j2))*ZDX(gt4,j2)) - SUM(j2,0,2,Conj(ZDX(gt2,3 + j2))*ZDX(gt4,3 + j2))));

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
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = 0.0125*(-80*SUM(j1,0,2,Conj(ZE(gt2,j1))*Ye(j1,j1)*ZE(gt4,3 + j1))*SUM(j2,0,2,Conj(Yd(j2,j2))*Conj(ZD(gt1,3 + j2))*ZD(gt3,j2)) - (4*Sqr(g1) + Sqr(gN))*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt4,3 + j1))*(SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) + 2*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2))) + 2*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt4,j1))*((Sqr(g1) - 5*Sqr(g2) - Sqr(gN))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) + 2*(Sqr(g1) - Sqr(gN))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2))) + 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) - 10*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) - 2*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) + 4*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) - 4*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) - 80*SUM(j1,0,2,Conj(ZD(gt1,j1))*Yd(j1,j1)*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(Ye(j2,j2))*Conj(ZE(gt2,3 + j2))*ZE(gt4,j2)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt4,3 + j2)) - Sqr(gN)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt4,3 + j2)) - 8*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt4,3 + j2)) - 2*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt4,3 + j2)));

   return {result};
}

ScalarVertex VertexImpl<fields::Sd, fields::SHI0, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::SHI0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto ZD = MODELPARAMETER(ZD);
   const auto UHI0 = MODELPARAMETER(UHI0);

   const std::complex<double> result = 0.0125*(SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*((2*Sqr(g1) + 10*Sqr(g2) + 3*Sqr(gN))*SUM(j2,0,1,Conj(UHI0(gt2,j2))*UHI0(gt4,j2)) - 2*(Sqr(g1) + 5*Sqr(g2) - Sqr(gN))*SUM(j2,0,1,Conj(UHI0(gt2,2 + j2))*UHI0(gt4,2 + j2))) + SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*((4*Sqr(g1) + 6*Sqr(gN))*SUM(j2,0,1,Conj(UHI0(gt2,j2))*UHI0(gt4,j2)) + 4*(-Sqr(g1) + Sqr(gN))*SUM(j2,0,1,Conj(UHI0(gt2,2 + j2))*UHI0(gt4,2 + j2))) + 2*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) + 10*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) + 3*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) - 2*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*UHI0(gt4,2 + j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) - 10*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*UHI0(gt4,2 + j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) + 2*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*UHI0(gt4,2 + j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) + 4*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2)) + 6*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2)) - 4*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*UHI0(gt4,2 + j1))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2)) + 4*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*UHI0(gt4,2 + j1))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2)));

   return {result};
}

ScalarVertex VertexImpl<fields::Sd, fields::SHIp, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::SHIp>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZD = MODELPARAMETER(ZD);
   const auto UHIp = MODELPARAMETER(UHIp);

   const std::complex<double> result = 0.0125*(SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*((4*Sqr(g1) + 6*Sqr(gN))*SUM(j2,0,1,Conj(UHIp(gt2,j2))*UHIp(gt4,j2)) + 4*(-Sqr(g1) + Sqr(gN))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*UHIp(gt4,2 + j2))) + SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*((2*Sqr(g1) - 10*Sqr(g2) + 3*Sqr(gN))*SUM(j2,0,1,Conj(UHIp(gt2,j2))*UHIp(gt4,j2)) + 2*(-Sqr(g1) + 5*Sqr(g2) + Sqr(gN))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*UHIp(gt4,2 + j2))) + 2*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) - 10*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) + 3*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) - 2*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*UHIp(gt4,2 + j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) + 10*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*UHIp(gt4,2 + j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) + 2*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*UHIp(gt4,2 + j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) + 4*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2)) + 6*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2)) - 4*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*UHIp(gt4,2 + j1))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2)) + 4*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*UHIp(gt4,2 + j1))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2)));

   return {result};
}

ScalarVertex VertexImpl<fields::Sd, fields::SHp0, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::SHp0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto ZD = MODELPARAMETER(ZD);
   const auto UHp0 = MODELPARAMETER(UHp0);

   const std::complex<double> result = 0.05*((Sqr(g1) + 5*Sqr(g2) - Sqr(gN))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1)) + 2*(Sqr(g1) - Sqr(gN))*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1)))*(Conj(UHp0(gt2,0))*UHp0(gt4,0) - Conj(UHp0(gt2,1))*UHp0(gt4,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Sd, fields::SHpp, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::SHpp>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto ZD = MODELPARAMETER(ZD);
   const auto UHpp = MODELPARAMETER(UHpp);

   const std::complex<double> result = 0.05*((Sqr(g1) - 5*Sqr(g2) - Sqr(gN))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1)) + 2*(Sqr(g1) - Sqr(gN))*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1)))*(Conj(UHpp(gt2,0))*UHpp(gt4,0) - Conj(UHpp(gt2,1))*UHpp(gt4,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Sd, fields::SSI0, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::SSI0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto gN = MODELPARAMETER(gN);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = -0.125*KroneckerDelta(gt2,gt4)*Sqr(gN)*(SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1)) + 2*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Sd, fields::Su, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto g3 = MODELPARAMETER(g3);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);

   const std::complex<double> result = 0.25*(-(Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZU(gt4,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZD(gt3,j2))) - 4*SUM(j1,0,2,Conj(ZD(gt1,j1))*Yu(j1,j1)*ZU(gt4,3 + j1))*SUM(j2,0,2,Conj(Yu(j2,j2))*Conj(ZU(gt2,3 + j2))*ZD(gt3,j2)) + Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt4,3 + j1))*(SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) - SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2))) + Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt4,j1))*(-SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) + SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2))) - Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZU(gt4,j2)) - 4*SUM(j1,0,2,Conj(ZU(gt2,j1))*Yd(j1,j1)*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(Yd(j2,j2))*Conj(ZD(gt1,3 + j2))*ZU(gt4,j2)) - Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) + Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) + Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)) - Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)));

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
   const auto gN = MODELPARAMETER(gN);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = 0.05*KroneckerDelta(gt2,gt4)*((Sqr(g1) + 5*Sqr(g2) - Sqr(gN))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1)) + 2*(Sqr(g1) - Sqr(gN))*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1)));

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
   const auto gN = MODELPARAMETER(gN);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.016666666666666666*((-4.898979485566356*g1*gN*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) + 30*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Cos(ThetaWp)) + Cos(ThetaW)*(-9.486832980505138*g2*gN*Sin(2*ThetaWp) + 15.491933384829668*g1*g2*Sin(ThetaW)*Sqr(Cos(ThetaWp))) + 2*Sqr(g1)*Sqr(Cos(ThetaWp))*Sqr(Sin(ThetaW)) + 3*Sqr(gN)*Sqr(Sin(ThetaWp)))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1)) + 4*(-4.898979485566356*g1*gN*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) + 2*Sqr(g1)*Sqr(Cos(ThetaWp))*Sqr(Sin(ThetaW)) + 3*Sqr(gN)*Sqr(Sin(ThetaWp)))*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VZp, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.016666666666666666*((3*Sqr(gN)*Sqr(Cos(ThetaWp)) + 30*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Sin(ThetaWp)) + g1*Sin(ThetaW)*(2.449489742783178*gN*Sin(2*ThetaWp) + 2*g1*Sin(ThetaW)*Sqr(Sin(ThetaWp))) + Cos(ThetaW)*(9.486832980505138*g2*gN*Sin(2*ThetaWp) + 15.491933384829668*g1*g2*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1)) + 4*(3*Sqr(gN)*Sqr(Cos(ThetaWp)) + g1*Sin(ThetaW)*(2.449489742783178*gN*Sin(2*ThetaWp) + 2*g1*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::hh, fields::Sd, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = std::complex<double>(0,0.5)*(Lambdax*SUM(j1,0,2,Conj(Yd(j1,j1))*Conj(ZD(gt3,3 + j1))*ZD(gt4,j1)) - Conj(Lambdax)*SUM(j1,0,2,Conj(ZD(gt3,j1))*Yd(j1,j1)*ZD(gt4,3 + j1)))*(ZA(gt1,2)*ZH(gt2,1) + ZA(gt1,1)*ZH(gt2,2));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::SDX, typename fields::conj<fields::SDX>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
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

ScalarVertex VertexImpl<fields::hh, fields::Se, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
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

ScalarVertex VertexImpl<fields::hh, fields::SHI0, typename fields::conj<fields::SHI0>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto TLambda12 = MODELPARAMETER(TLambda12);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto vu = MODELPARAMETER(vu);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto vs = MODELPARAMETER(vs);
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.025*(SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt3,j1))*(-(vd*(6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*ZH(gt1,0)) + 2*vu*(3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,1) + 15*vs*Sqr(gN)*ZH(gt1,2)) - 2*(10*Lambdax*SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*Conj(Lambda12(j1,j1))*UHI0(gt3,j1))*(vu*ZH(gt1,0) + vd*ZH(gt1,1)) + SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*UHI0(gt3,2 + j1))*(vd*(-3*Sqr(g1) - 5*Sqr(g2) + 3*Sqr(gN))*ZH(gt1,0) + vu*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZH(gt1,1) - 5*vs*Sqr(gN)*ZH(gt1,2)) + 10*(Conj(Lambdax)*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt3,2 + j1)*Lambda12(j1,j1))*(vu*ZH(gt1,0) + vd*ZH(gt1,1)) - (-2*vs*SUM(j1,0,1,AbsSqr(Lambda12(j1,j1))*Conj(UHI0(gt2,j1))*UHI0(gt3,j1)) + 1.4142135623730951*SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*Conj(TLambda12(j1,j1))*UHI0(gt3,j1)) - 2*vs*SUM(j1,0,1,AbsSqr(Lambda12(j1,j1))*Conj(UHI0(gt2,2 + j1))*UHI0(gt3,2 + j1)) + 1.4142135623730951*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt3,2 + j1)*TLambda12(j1,j1)))*ZH(gt1,2))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::SHIp, typename fields::conj<fields::SHIp>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
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

ScalarVertex VertexImpl<fields::hh, fields::SHp0, typename fields::conj<fields::SHp0>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto vs = MODELPARAMETER(vs);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = -0.05*(Conj(UHp0(gt2,0))*UHp0(gt3,0) - Conj(UHp0(gt2,1))*UHp0(gt3,1))*(vd*(3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,0) - vu*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZH(gt1,1) + 5*vs*Sqr(gN)*ZH(gt1,2));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::SHpp, typename fields::conj<fields::SHpp>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
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

ScalarVertex VertexImpl<fields::hh, fields::SSI0, typename fields::conj<fields::SSI0>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto gN = MODELPARAMETER(gN);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto vs = MODELPARAMETER(vs);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.125*KroneckerDelta(gt2,gt3)*Sqr(gN)*(3*vd*ZH(gt1,0) + 2*vu*ZH(gt1,1) - 5*vs*ZH(gt1,2));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Sv, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto vs = MODELPARAMETER(vs);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = -0.05*KroneckerDelta(gt2,gt3)*(vd*(3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,0) - vu*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZH(gt1,1) + 5*vs*Sqr(gN)*ZH(gt1,2));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VZ, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.016666666666666666*(-((9.486832980505138*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) - 3*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN) + 30*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 2.449489742783178*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + g1*(3.872983346207417*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) + g1*Sin(2*ThetaWp)*Sqr(Sin(ThetaW)) - 2.449489742783178*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1))) + 4*(Cos(ThetaWp)*Sin(ThetaWp)*(-Sqr(g1) + Cos(2*ThetaW)*Sqr(g1) + 3*Sqr(gN)) - 2.449489742783178*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + 2.449489742783178*g1*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp)))*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1)));

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
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,-0.25)*(-2*AbsSqr(Lambdax) + Sqr(g2))*(ZA(gt1,1)*ZH(gt2,0) + ZA(gt1,0)*ZH(gt2,1))*(ZP(gt3,1)*ZP(gt4,0) - ZP(gt3,0)*ZP(gt4,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::hh, fields::SDX, typename fields::conj<fields::SDX>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Kappa = MODELPARAMETER(Kappa);
   const auto ZDX = MODELPARAMETER(ZDX);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = std::complex<double>(0,0.5)*(Lambdax*SUM(j1,0,2,Conj(ZDX(gt3,3 + j1))*Conj(Kappa(j1,j1))*ZDX(gt4,j1)) - Conj(Lambdax)*SUM(j1,0,2,Conj(ZDX(gt3,j1))*ZDX(gt4,3 + j1)*Kappa(j1,j1)))*(ZA(gt1,1)*ZH(gt2,0) + ZA(gt1,0)*ZH(gt2,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::hh, fields::Se, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = std::complex<double>(0,0.5)*(Lambdax*SUM(j1,0,2,Conj(Ye(j1,j1))*Conj(ZE(gt3,3 + j1))*ZE(gt4,j1)) - Conj(Lambdax)*SUM(j1,0,2,Conj(ZE(gt3,j1))*Ye(j1,j1)*ZE(gt4,3 + j1)))*(ZA(gt1,2)*ZH(gt2,1) + ZA(gt1,1)*ZH(gt2,2));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::hh, fields::SHI0, typename fields::conj<fields::SHI0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = std::complex<double>(0,-0.5)*(Lambdax*SUM(j1,0,1,Conj(UHI0(gt3,2 + j1))*Conj(Lambda12(j1,j1))*UHI0(gt4,j1)) - Conj(Lambdax)*SUM(j1,0,1,Conj(UHI0(gt3,j1))*UHI0(gt4,2 + j1)*Lambda12(j1,j1)))*(ZA(gt1,1)*ZH(gt2,0) + ZA(gt1,0)*ZH(gt2,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::hh, fields::SHIp, typename fields::conj<fields::SHIp>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = std::complex<double>(0,0.5)*(Lambdax*SUM(j1,0,1,Conj(UHIp(gt3,2 + j1))*Conj(Lambda12(j1,j1))*UHIp(gt4,j1)) - Conj(Lambdax)*SUM(j1,0,1,Conj(UHIp(gt3,j1))*UHIp(gt4,2 + j1)*Lambda12(j1,j1)))*(ZA(gt1,1)*ZH(gt2,0) + ZA(gt1,0)*ZH(gt2,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::hh, fields::Su, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = std::complex<double>(0,0.5)*(Lambdax*SUM(j1,0,2,Conj(Yu(j1,j1))*Conj(ZU(gt3,3 + j1))*ZU(gt4,j1)) - Conj(Lambdax)*SUM(j1,0,2,Conj(ZU(gt3,j1))*Yu(j1,j1)*ZU(gt4,3 + j1)))*(ZA(gt1,2)*ZH(gt2,0) + ZA(gt1,0)*ZH(gt2,2));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::hh, fields::hh>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto gN = MODELPARAMETER(gN);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.025*(-(ZA(gt1,0)*ZA(gt2,0)*((6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*ZH(gt3,0)*ZH(gt4,0) + 2*(20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2) + 3*Sqr(gN))*ZH(gt3,1)*ZH(gt4,1) + 5*(8*AbsSqr(Lambdax) - 3*Sqr(gN))*ZH(gt3,2)*ZH(gt4,2))) + 5*ZA(gt1,2)*ZA(gt2,2)*((-8*AbsSqr(Lambdax) + 3*Sqr(gN))*ZH(gt3,0)*ZH(gt4,0) + 2*(-4*AbsSqr(Lambdax) + Sqr(gN))*ZH(gt3,1)*ZH(gt4,1) - 5*Sqr(gN)*ZH(gt3,2)*ZH(gt4,2)) + 2*ZA(gt1,1)*ZA(gt2,1)*((-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt3,0)*ZH(gt4,0) - (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZH(gt3,1)*ZH(gt4,1) + 5*(-4*AbsSqr(Lambdax) + Sqr(gN))*ZH(gt3,2)*ZH(gt4,2)));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::hh, fields::hh>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto gN = MODELPARAMETER(gN);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.025*(-(ZH(gt1,0)*(-2*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt2,1)*(ZH(gt3,1)*ZH(gt4,0) + ZH(gt3,0)*ZH(gt4,1)) - 5*(-8*AbsSqr(Lambdax) + 3*Sqr(gN))*ZH(gt2,2)*(ZH(gt3,2)*ZH(gt4,0) + ZH(gt3,0)*ZH(gt4,2)) + ZH(gt2,0)*(3*(6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*ZH(gt3,0)*ZH(gt4,0) + 2*(20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2) + 3*Sqr(gN))*ZH(gt3,1)*ZH(gt4,1) + 5*(8*AbsSqr(Lambdax) - 3*Sqr(gN))*ZH(gt3,2)*ZH(gt4,2)))) + 5*ZH(gt1,2)*((-8*AbsSqr(Lambdax) + 3*Sqr(gN))*ZH(gt2,0)*(ZH(gt3,2)*ZH(gt4,0) + ZH(gt3,0)*ZH(gt4,2)) + 2*(-4*AbsSqr(Lambdax) + Sqr(gN))*ZH(gt2,1)*(ZH(gt3,2)*ZH(gt4,1) + ZH(gt3,1)*ZH(gt4,2)) + ZH(gt2,2)*((-8*AbsSqr(Lambdax) + 3*Sqr(gN))*ZH(gt3,0)*ZH(gt4,0) + 2*(-4*AbsSqr(Lambdax) + Sqr(gN))*ZH(gt3,1)*ZH(gt4,1) - 15*Sqr(gN)*ZH(gt3,2)*ZH(gt4,2))) + 2*ZH(gt1,1)*((-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt2,0)*(ZH(gt3,1)*ZH(gt4,0) + ZH(gt3,0)*ZH(gt4,1)) + 5*(-4*AbsSqr(Lambdax) + Sqr(gN))*ZH(gt2,2)*(ZH(gt3,2)*ZH(gt4,1) + ZH(gt3,1)*ZH(gt4,2)) + ZH(gt2,1)*((-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt3,0)*ZH(gt4,0) - 3*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZH(gt3,1)*ZH(gt4,1) + 5*(-4*AbsSqr(Lambdax) + Sqr(gN))*ZH(gt3,2)*ZH(gt4,2))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::Hpm, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto gN = MODELPARAMETER(gN);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.025*(5*ZH(gt1,2)*ZH(gt2,2)*((-8*AbsSqr(Lambdax) + 3*Sqr(gN))*ZP(gt3,0)*ZP(gt4,0) + 2*(-4*AbsSqr(Lambdax) + Sqr(gN))*ZP(gt3,1)*ZP(gt4,1)) + 2*ZH(gt1,1)*(-5*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZH(gt2,0)*(ZP(gt3,1)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,1)) + ZH(gt2,1)*((3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZP(gt3,0)*ZP(gt4,0) - (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZP(gt3,1)*ZP(gt4,1))) - ZH(gt1,0)*(10*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZH(gt2,1)*(ZP(gt3,1)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,1)) + ZH(gt2,0)*((6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*ZP(gt3,0)*ZP(gt4,0) + 2*(-3*Sqr(g1) + 5*Sqr(g2) + 3*Sqr(gN))*ZP(gt3,1)*ZP(gt4,1))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::SDX, typename fields::conj<fields::SDX>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Kappa = MODELPARAMETER(Kappa);
   const auto ZDX = MODELPARAMETER(ZDX);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.025*(-2*SUM(j1,0,2,Conj(ZDX(gt3,j1))*ZDX(gt4,j1))*((2*Sqr(g1) + 3*Sqr(gN))*ZH(gt1,0)*ZH(gt2,0) + 2*(-Sqr(g1) + Sqr(gN))*ZH(gt1,1)*ZH(gt2,1) - 5*Sqr(gN)*ZH(gt1,2)*ZH(gt2,2)) + SUM(j1,0,2,Conj(ZDX(gt3,3 + j1))*ZDX(gt4,3 + j1))*((4*Sqr(g1) - 9*Sqr(gN))*ZH(gt1,0)*ZH(gt2,0) - 2*(2*Sqr(g1) + 3*Sqr(gN))*ZH(gt1,1)*ZH(gt2,1) + 15*Sqr(gN)*ZH(gt1,2)*ZH(gt2,2)) + 20*(Lambdax*SUM(j1,0,2,Conj(ZDX(gt3,3 + j1))*Conj(Kappa(j1,j1))*ZDX(gt4,j1))*(ZH(gt1,1)*ZH(gt2,0) + ZH(gt1,0)*ZH(gt2,1)) + Conj(Lambdax)*SUM(j1,0,2,Conj(ZDX(gt3,j1))*ZDX(gt4,3 + j1)*Kappa(j1,j1))*(ZH(gt1,1)*ZH(gt2,0) + ZH(gt1,0)*ZH(gt2,1)) - 2*(SUM(j1,0,2,AbsSqr(Kappa(j1,j1))*Conj(ZDX(gt3,j1))*ZDX(gt4,j1)) + SUM(j1,0,2,AbsSqr(Kappa(j1,j1))*Conj(ZDX(gt3,3 + j1))*ZDX(gt4,3 + j1)))*ZH(gt1,2)*ZH(gt2,2)));

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
   const auto gN = MODELPARAMETER(gN);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.025*(-40*SUM(j1,0,2,AbsSqr(Ye(j1,j1))*Conj(ZE(gt3,j1))*ZE(gt4,j1))*ZH(gt1,0)*ZH(gt2,0) + 12*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*ZE(gt4,3 + j1))*ZH(gt1,0)*ZH(gt2,0) + 3*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*ZE(gt4,3 + j1))*ZH(gt1,0)*ZH(gt2,0) - 40*SUM(j1,0,2,AbsSqr(Ye(j1,j1))*Conj(ZE(gt3,3 + j1))*ZE(gt4,3 + j1))*ZH(gt1,0)*ZH(gt2,0) - 12*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*ZE(gt4,3 + j1))*ZH(gt1,1)*ZH(gt2,1) + 2*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*ZE(gt4,3 + j1))*ZH(gt1,1)*ZH(gt2,1) + 20*Lambdax*SUM(j1,0,2,Conj(Ye(j1,j1))*Conj(ZE(gt3,3 + j1))*ZE(gt4,j1))*ZH(gt1,2)*ZH(gt2,1) + 20*Conj(Lambdax)*SUM(j1,0,2,Conj(ZE(gt3,j1))*Ye(j1,j1)*ZE(gt4,3 + j1))*ZH(gt1,2)*ZH(gt2,1) + 20*Lambdax*SUM(j1,0,2,Conj(Ye(j1,j1))*Conj(ZE(gt3,3 + j1))*ZE(gt4,j1))*ZH(gt1,1)*ZH(gt2,2) + 20*Conj(Lambdax)*SUM(j1,0,2,Conj(ZE(gt3,j1))*Ye(j1,j1)*ZE(gt4,3 + j1))*ZH(gt1,1)*ZH(gt2,2) - 5*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*ZE(gt4,3 + j1))*ZH(gt1,2)*ZH(gt2,2) - 2*SUM(j1,0,2,Conj(ZE(gt3,j1))*ZE(gt4,j1))*((3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,0)*ZH(gt2,0) + (-3*Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN))*ZH(gt1,1)*ZH(gt2,1) + 5*Sqr(gN)*ZH(gt1,2)*ZH(gt2,2)));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::SHI0, typename fields::conj<fields::SHI0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.025*(2*SUM(j1,0,1,Conj(UHI0(gt3,2 + j1))*UHI0(gt4,2 + j1))*((3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,0)*ZH(gt2,0) - (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZH(gt1,1)*ZH(gt2,1) + 5*Sqr(gN)*ZH(gt1,2)*ZH(gt2,2)) + SUM(j1,0,1,Conj(UHI0(gt3,j1))*UHI0(gt4,j1))*(-((6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*ZH(gt1,0)*ZH(gt2,0)) + 2*(3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,1)*ZH(gt2,1) + 15*Sqr(gN)*ZH(gt1,2)*ZH(gt2,2)) - 20*(Lambdax*SUM(j1,0,1,Conj(UHI0(gt3,2 + j1))*Conj(Lambda12(j1,j1))*UHI0(gt4,j1))*(ZH(gt1,1)*ZH(gt2,0) + ZH(gt1,0)*ZH(gt2,1)) + Conj(Lambdax)*SUM(j1,0,1,Conj(UHI0(gt3,j1))*UHI0(gt4,2 + j1)*Lambda12(j1,j1))*(ZH(gt1,1)*ZH(gt2,0) + ZH(gt1,0)*ZH(gt2,1)) + 2*(SUM(j1,0,1,AbsSqr(Lambda12(j1,j1))*Conj(UHI0(gt3,j1))*UHI0(gt4,j1)) + SUM(j1,0,1,AbsSqr(Lambda12(j1,j1))*Conj(UHI0(gt3,2 + j1))*UHI0(gt4,2 + j1)))*ZH(gt1,2)*ZH(gt2,2)));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::SHIp, typename fields::conj<fields::SHIp>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.025*(2*SUM(j1,0,1,Conj(UHIp(gt3,2 + j1))*UHIp(gt4,2 + j1))*((3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,0)*ZH(gt2,0) + (-3*Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN))*ZH(gt1,1)*ZH(gt2,1) + 5*Sqr(gN)*ZH(gt1,2)*ZH(gt2,2)) + SUM(j1,0,1,Conj(UHIp(gt3,j1))*UHIp(gt4,j1))*((-6*Sqr(g1) + 10*Sqr(g2) - 9*Sqr(gN))*ZH(gt1,0)*ZH(gt2,0) + 2*(3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,1)*ZH(gt2,1) + 15*Sqr(gN)*ZH(gt1,2)*ZH(gt2,2)) + 20*(Lambdax*SUM(j1,0,1,Conj(UHIp(gt3,2 + j1))*Conj(Lambda12(j1,j1))*UHIp(gt4,j1))*(ZH(gt1,1)*ZH(gt2,0) + ZH(gt1,0)*ZH(gt2,1)) + Conj(Lambdax)*SUM(j1,0,1,Conj(UHIp(gt3,j1))*UHIp(gt4,2 + j1)*Lambda12(j1,j1))*(ZH(gt1,1)*ZH(gt2,0) + ZH(gt1,0)*ZH(gt2,1)) - 2*(SUM(j1,0,1,AbsSqr(Lambda12(j1,j1))*Conj(UHIp(gt3,j1))*UHIp(gt4,j1)) + SUM(j1,0,1,AbsSqr(Lambda12(j1,j1))*Conj(UHIp(gt3,2 + j1))*UHIp(gt4,2 + j1)))*ZH(gt1,2)*ZH(gt2,2)));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::SHp0, typename fields::conj<fields::SHp0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = -0.05*(Conj(UHp0(gt3,0))*UHp0(gt4,0) - Conj(UHp0(gt3,1))*UHp0(gt4,1))*((3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,0)*ZH(gt2,0) - (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZH(gt1,1)*ZH(gt2,1) + 5*Sqr(gN)*ZH(gt1,2)*ZH(gt2,2));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::SHpp, typename fields::conj<fields::SHpp>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = -0.05*(Conj(UHpp(gt3,0))*UHpp(gt4,0) - Conj(UHpp(gt3,1))*UHpp(gt4,1))*((3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,0)*ZH(gt2,0) + (-3*Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN))*ZH(gt1,1)*ZH(gt2,1) + 5*Sqr(gN)*ZH(gt1,2)*ZH(gt2,2));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::SSI0, typename fields::conj<fields::SSI0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto gN = MODELPARAMETER(gN);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.125*KroneckerDelta(gt3,gt4)*Sqr(gN)*(3*ZH(gt1,0)*ZH(gt2,0) + 2*ZH(gt1,1)*ZH(gt2,1) - 5*ZH(gt1,2)*ZH(gt2,2));

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
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.025*(SUM(j1,0,2,Conj(ZU(gt3,3 + j1))*ZU(gt4,3 + j1))*((-8*Sqr(g1) + 3*Sqr(gN))*ZH(gt1,0)*ZH(gt2,0) + 2*(4*Sqr(g1) + Sqr(gN))*ZH(gt1,1)*ZH(gt2,1) - 5*Sqr(gN)*ZH(gt1,2)*ZH(gt2,2)) + SUM(j1,0,2,Conj(ZU(gt3,j1))*ZU(gt4,j1))*((2*Sqr(g1) - 10*Sqr(g2) + 3*Sqr(gN))*ZH(gt1,0)*ZH(gt2,0) + 2*(-Sqr(g1) + 5*Sqr(g2) + Sqr(gN))*ZH(gt1,1)*ZH(gt2,1) - 5*Sqr(gN)*ZH(gt1,2)*ZH(gt2,2)) + 20*(-2*(SUM(j1,0,2,AbsSqr(Yu(j1,j1))*Conj(ZU(gt3,j1))*ZU(gt4,j1)) + SUM(j1,0,2,AbsSqr(Yu(j1,j1))*Conj(ZU(gt3,3 + j1))*ZU(gt4,3 + j1)))*ZH(gt1,1)*ZH(gt2,1) + Lambdax*SUM(j1,0,2,Conj(Yu(j1,j1))*Conj(ZU(gt3,3 + j1))*ZU(gt4,j1))*(ZH(gt1,2)*ZH(gt2,0) + ZH(gt1,0)*ZH(gt2,2)) + Conj(Lambdax)*SUM(j1,0,2,Conj(ZU(gt3,j1))*Yu(j1,j1)*ZU(gt4,3 + j1))*(ZH(gt1,2)*ZH(gt2,0) + ZH(gt1,0)*ZH(gt2,2))));

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
   const auto gN = MODELPARAMETER(gN);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = -0.05*KroneckerDelta(gt3,gt4)*((3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZH(gt1,0)*ZH(gt2,0) - (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZH(gt1,1)*ZH(gt2,1) + 5*Sqr(gN)*ZH(gt1,2)*ZH(gt2,2));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::hh, typename fields::conj<fields::VWm>::type, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.5*Sqr(g2)*(ZH(gt1,0)*ZH(gt2,0) + ZH(gt1,1)*ZH(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::hh, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05*((-14.696938456699067*g1*gN*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) + 10*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Cos(ThetaWp)) + Cos(ThetaW)*(-9.486832980505138*g2*gN*Sin(2*ThetaWp) + 15.491933384829668*g1*g2*Sin(ThetaW)*Sqr(Cos(ThetaWp))) + 6*Sqr(g1)*Sqr(Cos(ThetaWp))*Sqr(Sin(ThetaW)) + 9*Sqr(gN)*Sqr(Sin(ThetaWp)))*ZH(gt1,0)*ZH(gt2,0) + 2*(3.1622776601683795*g2*gN*Cos(ThetaW)*Sin(2*ThetaWp) + g1*Sin(ThetaW)*(7.745966692414834*g2*Cos(ThetaW) + 3*g1*Sin(ThetaW))*Sqr(Cos(ThetaWp)) + 5*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Cos(ThetaWp)) + gN*(2.449489742783178*g1*Sin(ThetaW)*Sin(2*ThetaWp) + 2*gN*Sqr(Sin(ThetaWp))))*ZH(gt1,1)*ZH(gt2,1) + 25*Sqr(gN)*Sqr(Sin(ThetaWp))*ZH(gt1,2)*ZH(gt2,2));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::hh, fields::VZp, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*((9*Sqr(gN)*Sqr(Cos(ThetaWp)) + (g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(9.486832980505138*gN*Sin(2*ThetaWp) + 10*g2*Cos(ThetaW)*Sqr(Sin(ThetaWp)) + 7.745966692414834*g1*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*ZH(gt1,0)*ZH(gt2,0) + 2*(-3.1622776601683795*gN*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*Sin(2*ThetaWp) + 2*Sqr(gN)*Sqr(Cos(ThetaWp)) + 5*Sqr(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*Sqr(Sin(ThetaWp)))*ZH(gt1,1)*ZH(gt2,1) + 25*Sqr(gN)*Sqr(Cos(ThetaWp))*ZH(gt1,2)*ZH(gt2,2));

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
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -0.35355339059327373*(Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZU(gt4,j1))*(ZH(gt1,0)*ZP(gt3,0) + ZH(gt1,1)*ZP(gt3,1)) - 2*(SUM(j1,0,2,AbsSqr(Yd(j1,j1))*Conj(ZD(gt2,j1))*ZU(gt4,j1))*ZH(gt1,0)*ZP(gt3,0) + Conj(Lambdax)*SUM(j1,0,2,Conj(ZD(gt2,j1))*Yu(j1,j1)*ZU(gt4,3 + j1))*ZH(gt1,2)*ZP(gt3,0) + SUM(j1,0,2,AbsSqr(Yu(j1,j1))*Conj(ZD(gt2,j1))*ZU(gt4,j1))*ZH(gt1,1)*ZP(gt3,1) + Lambdax*SUM(j1,0,2,Conj(Yd(j1,j1))*Conj(ZD(gt2,3 + j1))*ZU(gt4,j1))*ZH(gt1,2)*ZP(gt3,1) + SUM(j1,0,2,Conj(Yd(j1,j1))*Conj(ZD(gt2,3 + j1))*Yu(j1,j1)*ZU(gt4,3 + j1))*(ZH(gt1,1)*ZP(gt3,0) + ZH(gt1,0)*ZP(gt3,1))));

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
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -0.35355339059327373*(Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt3,j1))*ZD(gt4,j1))*(ZH(gt1,0)*ZP(gt2,0) + ZH(gt1,1)*ZP(gt2,1)) - 2*(SUM(j1,0,2,AbsSqr(Yd(j1,j1))*Conj(ZU(gt3,j1))*ZD(gt4,j1))*ZH(gt1,0)*ZP(gt2,0) + Lambdax*SUM(j1,0,2,Conj(Yu(j1,j1))*Conj(ZU(gt3,3 + j1))*ZD(gt4,j1))*ZH(gt1,2)*ZP(gt2,0) + SUM(j1,0,2,AbsSqr(Yu(j1,j1))*Conj(ZU(gt3,j1))*ZD(gt4,j1))*ZH(gt1,1)*ZP(gt2,1) + Conj(Lambdax)*SUM(j1,0,2,Conj(ZU(gt3,j1))*Yd(j1,j1)*ZD(gt4,3 + j1))*ZH(gt1,2)*ZP(gt2,1) + SUM(j1,0,2,Conj(Yu(j1,j1))*Conj(ZU(gt3,3 + j1))*Yd(j1,j1)*ZD(gt4,3 + j1))*(ZH(gt1,1)*ZP(gt2,0) + ZH(gt1,0)*ZP(gt2,1))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Ah>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto TLambdax = MODELPARAMETER(TLambdax);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0.,0.35355339059327373)*(Conj(TLambdax) - TLambdax)*(ZA(gt1,2)*(ZA(gt2,1)*ZA(gt3,0) + ZA(gt2,0)*ZA(gt3,1)) + ZA(gt1,1)*(ZA(gt2,2)*ZA(gt3,0) + ZA(gt2,0)*ZA(gt3,2)) + ZA(gt1,0)*(ZA(gt2,2)*ZA(gt3,1) + ZA(gt2,1)*ZA(gt3,2)));

   return {result};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha>::type, fields::Cha, fields::Ah>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto UM = MODELPARAMETER(UM);
   const auto UP = MODELPARAMETER(UP);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,0.7071067811865475)*(g2*Conj(UM(gt2,0))*Conj(UP(gt1,1))*ZA(gt3,1) + Conj(UM(gt2,1))*(g2*Conj(UP(gt1,0))*ZA(gt3,0) - Conj(UP(gt1,1))*Lambdax*ZA(gt3,2)));

   const std::complex<double> right = std::complex<double>(0.,-0.7071067811865475)*(g2*UM(gt1,0)*UP(gt2,1)*ZA(gt3,1) + UM(gt1,1)*(g2*UP(gt2,0)*ZA(gt3,0) - Conj(Lambdax)*UP(gt2,1)*ZA(gt3,2)));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha>::type, fields::Cha, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto UM = MODELPARAMETER(UM);
   const auto UP = MODELPARAMETER(UP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = g2*Conj(UM(gt2,0))*Cos(ThetaW)*Cos(ThetaWp)*UM(gt1,0) + 0.05*Conj(UM(gt2,1))*(10*g2*Cos(ThetaW)*Cos(ThetaWp) - 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp))*UM(gt1,1);

   const std::complex<double> right = g2*Conj(UP(gt1,0))*Cos(ThetaW)*Cos(ThetaWp)*UP(gt2,0) + 0.1*Conj(UP(gt1,1))*(5*g2*Cos(ThetaW)*Cos(ThetaWp) - 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) - 3.1622776601683795*gN*Sin(ThetaWp))*UP(gt2,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha>::type, fields::Cha, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto UM = MODELPARAMETER(UM);
   const auto UP = MODELPARAMETER(UP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = -(g2*Conj(UM(gt2,0))*Cos(ThetaW)*Sin(ThetaWp)*UM(gt1,0)) + 0.05*Conj(UM(gt2,1))*(9.486832980505138*gN*Cos(ThetaWp) + 2*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*UM(gt1,1);

   const std::complex<double> right = -(g2*Conj(UP(gt1,0))*Cos(ThetaW)*Sin(ThetaWp)*UP(gt2,0)) - 0.1*Conj(UP(gt1,1))*(3.1622776601683795*gN*Cos(ThetaWp) + (5*g2*Cos(ThetaW) - 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*UP(gt2,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::ChaI>::type, fields::ChaI, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto ZMI = MODELPARAMETER(ZMI);
   const auto ZPI = MODELPARAMETER(ZPI);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*SUM(j1,0,1,Conj(ZMI(gt2,j1))*Conj(ZPI(gt1,j1))*Lambda12(j1,j1))*ZH(gt3,2);

   const std::complex<double> right = -0.7071067811865475*SUM(j1,0,1,Conj(Lambda12(j1,j1))*ZMI(gt1,j1)*ZPI(gt2,j1))*ZH(gt3,2);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::ChaI>::type, fields::ChaI, fields::Ah>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto ZMI = MODELPARAMETER(ZMI);
   const auto ZPI = MODELPARAMETER(ZPI);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*SUM(j1,0,1,Conj(ZMI(gt2,j1))*Conj(ZPI(gt1,j1))*Lambda12(j1,j1))*ZA(gt3,2);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*SUM(j1,0,1,Conj(Lambda12(j1,j1))*ZMI(gt1,j1)*ZPI(gt2,j1))*ZA(gt3,2);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::ChaI>::type, fields::ChaI, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.5*KroneckerDelta(gt1,gt2)*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   const std::complex<double> right = 0.5*KroneckerDelta(gt1,gt2)*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::ChaI>::type, fields::ChaI, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = 0.05*KroneckerDelta(gt1,gt2)*(10*g2*Cos(ThetaW)*Cos(ThetaWp) - 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp));

   const std::complex<double> right = 0.1*KroneckerDelta(gt1,gt2)*(5*g2*Cos(ThetaW)*Cos(ThetaWp) - 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) - 3.1622776601683795*gN*Sin(ThetaWp));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::ChaI>::type, fields::ChaI, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.05*KroneckerDelta(gt1,gt2)*(9.486832980505138*gN*Cos(ThetaWp) + 2*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp));

   const std::complex<double> right = -0.1*KroneckerDelta(gt1,gt2)*(3.1622776601683795*gN*Cos(ThetaWp) + (5*g2*Cos(ThetaW) - 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Chi, fields::Chi, fields::Ah>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0,0.05)*(10*g2*Conj(ZN(gt1,1))*Conj(ZN(gt2,2))*ZA(gt3,0) - 9.486832980505138*gN*Conj(ZN(gt1,5))*Conj(ZN(gt2,2))*ZA(gt3,0) + 14.142135623730951*Conj(ZN(gt1,4))*Conj(ZN(gt2,3))*Lambdax*ZA(gt3,0) + 14.142135623730951*Conj(ZN(gt1,3))*Conj(ZN(gt2,4))*Lambdax*ZA(gt3,0) + 7.745966692414834*g1*Conj(ZN(gt1,3))*Conj(ZN(gt2,0))*ZA(gt3,1) - 10*g2*Conj(ZN(gt1,3))*Conj(ZN(gt2,1))*ZA(gt3,1) - 10*g2*Conj(ZN(gt1,1))*Conj(ZN(gt2,3))*ZA(gt3,1) - 6.324555320336759*gN*Conj(ZN(gt1,5))*Conj(ZN(gt2,3))*ZA(gt3,1) - 6.324555320336759*gN*Conj(ZN(gt1,3))*Conj(ZN(gt2,5))*ZA(gt3,1) + 14.142135623730951*Conj(ZN(gt1,4))*Conj(ZN(gt2,2))*Lambdax*ZA(gt3,1) + 7.745966692414834*g1*Conj(ZN(gt1,0))*(-(Conj(ZN(gt2,2))*ZA(gt3,0)) + Conj(ZN(gt2,3))*ZA(gt3,1)) + 15.811388300841898*gN*Conj(ZN(gt1,5))*Conj(ZN(gt2,4))*ZA(gt3,2) + 15.811388300841898*gN*Conj(ZN(gt1,4))*Conj(ZN(gt2,5))*ZA(gt3,2) + 14.142135623730951*Conj(ZN(gt1,3))*Conj(ZN(gt2,2))*Lambdax*ZA(gt3,2) + Conj(ZN(gt1,2))*(-7.745966692414834*g1*Conj(ZN(gt2,0))*ZA(gt3,0) + 10*g2*Conj(ZN(gt2,1))*ZA(gt3,0) - 9.486832980505138*gN*Conj(ZN(gt2,5))*ZA(gt3,0) + 14.142135623730951*Lambdax*(Conj(ZN(gt2,4))*ZA(gt3,1) + Conj(ZN(gt2,3))*ZA(gt3,2))));

   const std::complex<double> right = std::complex<double>(0,-0.05)*(ZA(gt3,0)*(-7.745966692414834*g1*ZN(gt1,0)*ZN(gt2,2) + 10*g2*ZN(gt1,1)*ZN(gt2,2) - 9.486832980505138*gN*ZN(gt1,5)*ZN(gt2,2) + 14.142135623730951*Conj(Lambdax)*ZN(gt1,4)*ZN(gt2,3) + 14.142135623730951*Conj(Lambdax)*ZN(gt1,3)*ZN(gt2,4) + ZN(gt1,2)*(-7.745966692414834*g1*ZN(gt2,0) + 10*g2*ZN(gt2,1) - 9.486832980505138*gN*ZN(gt2,5))) + 2*ZA(gt3,1)*((3.872983346207417*g1*ZN(gt1,0) - 5*g2*ZN(gt1,1) - 3.1622776601683795*gN*ZN(gt1,5))*ZN(gt2,3) + 7.0710678118654755*Conj(Lambdax)*(ZN(gt1,4)*ZN(gt2,2) + ZN(gt1,2)*ZN(gt2,4)) + ZN(gt1,3)*(3.872983346207417*g1*ZN(gt2,0) - 5*g2*ZN(gt2,1) - 3.1622776601683795*gN*ZN(gt2,5))) + 7.0710678118654755*ZA(gt3,2)*(2*Conj(Lambdax)*(ZN(gt1,3)*ZN(gt2,2) + ZN(gt1,2)*ZN(gt2,3)) + 2.23606797749979*gN*(ZN(gt1,5)*ZN(gt2,4) + ZN(gt1,4)*ZN(gt2,5))));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Chi, fields::Chi, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = 0.05*(Conj(ZN(gt2,2))*(-10*g2*Cos(ThetaW)*Cos(ThetaWp) - 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp))*ZN(gt1,2) + 2*Conj(ZN(gt2,3))*(5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*ZN(gt1,3) - 15.811388300841898*gN*Conj(ZN(gt2,4))*Sin(ThetaWp)*ZN(gt1,4));

   const std::complex<double> right = 0.05*(Conj(ZN(gt1,2))*(10*g2*Cos(ThetaW)*Cos(ThetaWp) + 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp))*ZN(gt2,2) - 2*Conj(ZN(gt1,3))*(5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*ZN(gt2,3) + 15.811388300841898*gN*Conj(ZN(gt1,4))*Sin(ThetaWp)*ZN(gt2,4));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Chi, fields::Chi, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.05*(Conj(ZN(gt2,2))*(9.486832980505138*gN*Cos(ThetaWp) + 2*(5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*ZN(gt1,2) + 2*Conj(ZN(gt2,3))*(3.1622776601683795*gN*Cos(ThetaWp) - (5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*ZN(gt1,3) - 15.811388300841898*gN*Conj(ZN(gt2,4))*Cos(ThetaWp)*ZN(gt1,4));

   const std::complex<double> right = 0.05*(-(Conj(ZN(gt1,2))*(9.486832980505138*gN*Cos(ThetaWp) + 2*(5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*ZN(gt2,2)) + 2*Conj(ZN(gt1,3))*(-3.1622776601683795*gN*Cos(ThetaWp) + 5*g2*Cos(ThetaW)*Sin(ThetaWp) + 3.872983346207417*g1*Sin(ThetaW)*Sin(ThetaWp))*ZN(gt2,3) + 15.811388300841898*gN*Conj(ZN(gt1,4))*Cos(ThetaWp)*ZN(gt2,4));

   return {left, right};
}

ChiralVertex VertexImpl<fields::ChiI, fields::ChiI, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto ZNI = MODELPARAMETER(ZNI);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = 0.7071067811865475*(SUM(j1,0,1,Conj(ZNI(gt1,2 + j1))*Conj(ZNI(gt2,j1))*Lambda12(j1,j1)) + SUM(j1,0,1,Conj(ZNI(gt1,j1))*Conj(ZNI(gt2,2 + j1))*Lambda12(j1,j1)))*ZH(gt3,2);

   const std::complex<double> right = 0.7071067811865475*(SUM(j1,0,1,Conj(Lambda12(j1,j1))*ZNI(gt1,2 + j1)*ZNI(gt2,j1)) + SUM(j1,0,1,Conj(Lambda12(j1,j1))*ZNI(gt1,j1)*ZNI(gt2,2 + j1)))*ZH(gt3,2);

   return {left, right};
}

ChiralVertex VertexImpl<fields::ChiI, fields::ChiI, fields::Ah>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto ZNI = MODELPARAMETER(ZNI);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,0.7071067811865475)*(SUM(j1,0,1,Conj(ZNI(gt1,2 + j1))*Conj(ZNI(gt2,j1))*Lambda12(j1,j1)) + SUM(j1,0,1,Conj(ZNI(gt1,j1))*Conj(ZNI(gt2,2 + j1))*Lambda12(j1,j1)))*ZA(gt3,2);

   const std::complex<double> right = std::complex<double>(0.,-0.7071067811865475)*(SUM(j1,0,1,Conj(Lambda12(j1,j1))*ZNI(gt1,2 + j1)*ZNI(gt2,j1)) + SUM(j1,0,1,Conj(Lambda12(j1,j1))*ZNI(gt1,j1)*ZNI(gt2,2 + j1)))*ZA(gt3,2);

   return {left, right};
}

ChiralVertex VertexImpl<fields::ChiI, fields::ChiI, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZNI = MODELPARAMETER(ZNI);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = 0.05*((-10*g2*Cos(ThetaW)*Cos(ThetaWp) - 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp))*SUM(j1,0,1,Conj(ZNI(gt2,j1))*ZNI(gt1,j1)) + 2*(5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*SUM(j1,0,1,Conj(ZNI(gt2,2 + j1))*ZNI(gt1,2 + j1)));

   const std::complex<double> right = 0.05*((10*g2*Cos(ThetaW)*Cos(ThetaWp) + 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp))*SUM(j1,0,1,Conj(ZNI(gt1,j1))*ZNI(gt2,j1)) - 2*(5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*SUM(j1,0,1,Conj(ZNI(gt1,2 + j1))*ZNI(gt2,2 + j1)));

   return {left, right};
}

ChiralVertex VertexImpl<fields::ChiI, fields::ChiI, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZNI = MODELPARAMETER(ZNI);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.05*((9.486832980505138*gN*Cos(ThetaWp) + 2*(5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*SUM(j1,0,1,Conj(ZNI(gt2,j1))*ZNI(gt1,j1)) + 2*(3.1622776601683795*gN*Cos(ThetaWp) - (5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*SUM(j1,0,1,Conj(ZNI(gt2,2 + j1))*ZNI(gt1,2 + j1)));

   const std::complex<double> right = 0.05*(-((9.486832980505138*gN*Cos(ThetaWp) + 2*(5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*SUM(j1,0,1,Conj(ZNI(gt1,j1))*ZNI(gt2,j1))) + 2*(-3.1622776601683795*gN*Cos(ThetaWp) + (5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*SUM(j1,0,1,Conj(ZNI(gt1,2 + j1))*ZNI(gt2,2 + j1)));

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

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*SUM(j1,0,2,Conj(ZDL(gt2,j1))*Conj(ZDR(gt1,j1))*Yd(j1,j1))*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*SUM(j1,0,2,Conj(Yd(j1,j1))*ZDL(gt1,j1)*ZDR(gt2,j1))*ZA(gt3,0);

   return {left, right};
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
   const auto gN = MODELPARAMETER(gN);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = 0.016666666666666666*KroneckerDelta(gt1,gt2)*(30*g2*Cos(ThetaW)*Cos(ThetaWp) + 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp));

   const std::complex<double> right = KroneckerDelta(gt1,gt2)*(-0.2581988897471611*g1*Cos(ThetaWp)*Sin(ThetaW) + 0.31622776601683794*gN*Sin(ThetaWp));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fd, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.016666666666666666*KroneckerDelta(gt1,gt2)*(9.486832980505138*gN*Cos(ThetaWp) + 2*(15*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp));

   const std::complex<double> right = KroneckerDelta(gt1,gt2)*(0.31622776601683794*gN*Cos(ThetaWp) + 0.2581988897471611*g1*Sin(ThetaW)*Sin(ThetaWp));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::FDX>::type, fields::FDX, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Kappa = MODELPARAMETER(Kappa);
   const auto ZDXL = MODELPARAMETER(ZDXL);
   const auto ZDXR = MODELPARAMETER(ZDXR);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*SUM(j1,0,2,Conj(ZDXL(gt2,j1))*Conj(ZDXR(gt1,j1))*Kappa(j1,j1))*ZH(gt3,2);

   const std::complex<double> right = -0.7071067811865475*SUM(j1,0,2,Conj(Kappa(j1,j1))*ZDXL(gt1,j1)*ZDXR(gt2,j1))*ZH(gt3,2);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::FDX>::type, fields::FDX, fields::Ah>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Kappa = MODELPARAMETER(Kappa);
   const auto ZDXL = MODELPARAMETER(ZDXL);
   const auto ZDXR = MODELPARAMETER(ZDXR);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*SUM(j1,0,2,Conj(ZDXL(gt2,j1))*Conj(ZDXR(gt1,j1))*Kappa(j1,j1))*ZA(gt3,2);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*SUM(j1,0,2,Conj(Kappa(j1,j1))*ZDXL(gt1,j1)*ZDXR(gt2,j1))*ZA(gt3,2);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::FDX>::type, fields::FDX, fields::VG>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g3 = MODELPARAMETER(g3);

   const std::complex<double> left = -0.5*g3*KroneckerDelta(gt1,gt2);

   const std::complex<double> right = -0.5*g3*KroneckerDelta(gt1,gt2);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::FDX>::type, fields::FDX, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.2581988897471611*g1*Cos(ThetaW)*KroneckerDelta(gt1,gt2);

   const std::complex<double> right = 0.2581988897471611*g1*Cos(ThetaW)*KroneckerDelta(gt1,gt2);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::FDX>::type, fields::FDX, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = KroneckerDelta(gt1,gt2)*(-0.2581988897471611*g1*Cos(ThetaWp)*Sin(ThetaW) + 0.31622776601683794*gN*Sin(ThetaWp));

   const std::complex<double> right = -0.016666666666666666*KroneckerDelta(gt1,gt2)*(15.491933384829668*g1*Cos(ThetaWp)*Sin(ThetaW) + 28.460498941515414*gN*Sin(ThetaWp));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::FDX>::type, fields::FDX, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = KroneckerDelta(gt1,gt2)*(0.31622776601683794*gN*Cos(ThetaWp) + 0.2581988897471611*g1*Sin(ThetaW)*Sin(ThetaWp));

   const std::complex<double> right = -0.016666666666666666*KroneckerDelta(gt1,gt2)*(28.460498941515414*gN*Cos(ThetaWp) - 15.491933384829668*g1*Sin(ThetaW)*Sin(ThetaWp));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fe, fields::VZ>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fe, fields::VZp>::evaluate(
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

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*SUM(j1,0,2,Conj(ZUL(gt2,j1))*Conj(ZUR(gt1,j1))*Yu(j1,j1))*ZA(gt3,1);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*SUM(j1,0,2,Conj(Yu(j1,j1))*ZUL(gt1,j1)*ZUR(gt2,j1))*ZA(gt3,1);

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
   const auto gN = MODELPARAMETER(gN);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = -0.016666666666666666*KroneckerDelta(gt1,gt2)*(30*g2*Cos(ThetaW)*Cos(ThetaWp) - 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp));

   const std::complex<double> right = 0.016666666666666666*KroneckerDelta(gt1,gt2)*(30.983866769659336*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fu, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.016666666666666666*KroneckerDelta(gt1,gt2)*(9.486832980505138*gN*Cos(ThetaWp) + 2*(-15*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp));

   const std::complex<double> right = 0.016666666666666666*KroneckerDelta(gt1,gt2)*(9.486832980505138*gN*Cos(ThetaWp) - 30.983866769659336*g1*Sin(ThetaW)*Sin(ThetaWp));

   return {left, right};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gWm>::type, fields::gWm, fields::hh>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = -0.25*Sqr(g2)*(vd*ZH(gt3,0) + vu*ZH(gt3,1));

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
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -(g2*Cos(ThetaW)*Cos(ThetaWp));

   return {result, 1};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWm>::type, fields::gWm, fields::VZp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = g2*Cos(ThetaW)*Sin(ThetaWp);

   return {result, 1};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gWmC, fields::hh>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = -0.25*Sqr(g2)*(vd*ZH(gt3,0) + vu*ZH(gt3,1));

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
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = g2*Cos(ThetaW)*Cos(ThetaWp);

   return {result, 1};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gWmC, fields::VZp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -(g2*Cos(ThetaW)*Sin(ThetaWp));

   return {result, 1};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gZ, fields::hh>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto vd = MODELPARAMETER(vd);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto vu = MODELPARAMETER(vu);
   const auto vs = MODELPARAMETER(vs);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.025*(-(vd*(-14.696938456699067*g1*gN*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) + 10*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Cos(ThetaWp)) + Cos(ThetaW)*(-18.973665961010276*g2*gN*Cos(ThetaWp)*Sin(ThetaWp) + 15.491933384829668*g1*g2*Sin(ThetaW)*Sqr(Cos(ThetaWp))) + 6*Sqr(g1)*Sqr(Cos(ThetaWp))*Sqr(Sin(ThetaW)) + 9*Sqr(gN)*Sqr(Sin(ThetaWp)))*ZH(gt3,0)) - 2*vu*(3.1622776601683795*g2*gN*Cos(ThetaW)*Sin(2*ThetaWp) + g1*Sin(ThetaW)*(7.745966692414834*g2*Cos(ThetaW) + 3*g1*Sin(ThetaW))*Sqr(Cos(ThetaWp)) + 5*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Cos(ThetaWp)) + gN*(2.449489742783178*g1*Sin(ThetaW)*Sin(2*ThetaWp) + 2*gN*Sqr(Sin(ThetaWp))))*ZH(gt3,1) - 25*vs*Sqr(gN)*Sqr(Sin(ThetaWp))*ZH(gt3,2));

   return {result};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gZp>::type, fields::gZ, fields::hh>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto vd = MODELPARAMETER(vd);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto vs = MODELPARAMETER(vs);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.025*(vd*(9.486832980505138*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) - 9*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN) + 5*Sin(2*ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 7.348469228349534*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + g1*(3.872983346207417*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) + 3*g1*Sin(2*ThetaWp)*Sqr(Sin(ThetaW)) - 7.348469228349534*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*ZH(gt3,0) + vu*(-6.324555320336759*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) - 4*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN) + 5*Sin(2*ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) - 4.898979485566356*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + g1*(3.872983346207417*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) + 3*g1*Sin(2*ThetaWp)*Sqr(Sin(ThetaW)) + 4.898979485566356*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*ZH(gt3,1) - 25*vs*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN)*ZH(gt3,2));

   return {result};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gZp, fields::hh>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto vd = MODELPARAMETER(vd);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto vs = MODELPARAMETER(vs);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.025*(vd*(9.486832980505138*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) - 9*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN) + 5*Sin(2*ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 7.348469228349534*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + g1*(3.872983346207417*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) + 3*g1*Sin(2*ThetaWp)*Sqr(Sin(ThetaW)) - 7.348469228349534*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*ZH(gt3,0) + vu*(-6.324555320336759*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) - 4*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN) + 5*Sin(2*ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) - 4.898979485566356*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + g1*(3.872983346207417*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) + 3*g1*Sin(2*ThetaWp)*Sqr(Sin(ThetaW)) + 4.898979485566356*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*ZH(gt3,1) - 25*vs*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN)*ZH(gt3,2));

   return {result};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gZp>::type, fields::gZp, fields::hh>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto vd = MODELPARAMETER(vd);
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto vs = MODELPARAMETER(vs);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.025*(-(vd*(9*Sqr(gN)*Sqr(Cos(ThetaWp)) + (g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(9.486832980505138*gN*Sin(2*ThetaWp) + 10*g2*Cos(ThetaW)*Sqr(Sin(ThetaWp)) + 7.745966692414834*g1*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*ZH(gt3,0)) - 2*vu*(-3.1622776601683795*g2*gN*Cos(ThetaW)*Sin(2*ThetaWp) + 2*Sqr(gN)*Sqr(Cos(ThetaWp)) + 5*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Sin(ThetaWp)) + g1*(-2.449489742783178*gN*Sin(ThetaW)*Sin(2*ThetaWp) + 3.872983346207417*g2*Sin(2*ThetaW)*Sqr(Sin(ThetaWp)) + 3*g1*Sqr(Sin(ThetaW))*Sqr(Sin(ThetaWp))))*ZH(gt3,1) - 25*vs*Sqr(gN)*Sqr(Cos(ThetaWp))*ZH(gt3,2));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Hpm, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
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

   const std::complex<double> result = 0.5*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*(ZP(gt1,0)*ZP(gt2,0) + ZP(gt1,1)*ZP(gt2,1));

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
   const auto gN = MODELPARAMETER(gN);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05*((10*g2*Cos(ThetaW)*Cos(ThetaWp) - 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp))*ZP(gt1,0)*ZP(gt2,0) + 2*(5*g2*Cos(ThetaW)*Cos(ThetaWp) - 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) - 3.1622776601683795*gN*Sin(ThetaWp))*ZP(gt1,1)*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*((9.486832980505138*gN*Cos(ThetaWp) + 2*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*ZP(gt1,0)*ZP(gt2,0) - 2*(3.1622776601683795*gN*Cos(ThetaWp) + (5*g2*Cos(ThetaW) - 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*ZP(gt1,1)*ZP(gt2,1));

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

   const std::complex<double> result = std::complex<double>(0,-0.5)*g2*(ZA(gt1,0)*ZP(gt2,0) + ZA(gt1,1)*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<typename fields::conj<fields::Hpm>::type, fields::VP, fields::VWm>::evaluate(
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

InverseMetricVertex VertexImpl<typename fields::conj<fields::Hpm>::type, fields::VWm, fields::VZ>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*g2*(vd*(7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp))*ZP(gt1,0) - 2*vu*(3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*ZP(gt1,1));

   return {result};
}

InverseMetricVertex VertexImpl<typename fields::conj<fields::Hpm>::type, fields::VWm, fields::VZp>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.05*g2*(vd*(9.486832980505138*gN*Cos(ThetaWp) + 7.745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt1,0) + 2*vu*(3.1622776601683795*gN*Cos(ThetaWp) - 3.872983346207417*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt1,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::SDX, typename fields::conj<fields::SDX>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
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

MomentumDifferenceVertex VertexImpl<fields::SDX, typename fields::conj<fields::SDX>::type, fields::VG>::evaluate(
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

MomentumDifferenceVertex VertexImpl<fields::SDX, typename fields::conj<fields::SDX>::type, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.2581988897471611*g1*Cos(ThetaW)*KroneckerDelta(gt1,gt2);

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::SDX, typename fields::conj<fields::SDX>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZDX = MODELPARAMETER(ZDX);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.016666666666666666*(-2*(7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt2,j1)) - (15.491933384829668*g1*Cos(ThetaWp)*Sin(ThetaW) + 28.460498941515414*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt2,3 + j1)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::SDX, typename fields::conj<fields::SDX>::type, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZDX = MODELPARAMETER(ZDX);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.016666666666666666*((18.973665961010276*gN*Cos(ThetaWp) + 15.491933384829668*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt2,j1)) + (-28.460498941515414*gN*Cos(ThetaWp) + 15.491933384829668*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt2,3 + j1)));

   return {result, minuend_index, subtrahend_index};
}

ScalarVertex VertexImpl<fields::Ah, fields::Se, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
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
   const auto gN = MODELPARAMETER(gN);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05*(2*(5*g2*Cos(ThetaW)*Cos(ThetaWp) - 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) - 3.1622776601683795*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) + (-15.491933384829668*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::Se, typename fields::conj<fields::Se>::type, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*(-2*(3.1622776601683795*gN*Cos(ThetaWp) + (5*g2*Cos(ThetaW) - 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) + (3.1622776601683795*gN*Cos(ThetaWp) + 15.491933384829668*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1)));

   return {result, minuend_index, subtrahend_index};
}

ScalarVertex VertexImpl<fields::Ah, fields::SHI0, typename fields::conj<fields::SHI0>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto TLambda12 = MODELPARAMETER(TLambda12);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto vu = MODELPARAMETER(vu);
   const auto vd = MODELPARAMETER(vd);
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,-0.5)*(Lambdax*SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*Conj(Lambda12(j1,j1))*UHI0(gt3,j1))*(vu*ZA(gt1,0) + vd*ZA(gt1,1)) - Conj(Lambdax)*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt3,2 + j1)*Lambda12(j1,j1))*(vu*ZA(gt1,0) + vd*ZA(gt1,1)) + 1.4142135623730951*(SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*Conj(TLambda12(j1,j1))*UHI0(gt3,j1)) - SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt3,2 + j1)*TLambda12(j1,j1)))*ZA(gt1,2));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::SHI0, typename fields::conj<fields::SHI0>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05*((-10*g2*Cos(ThetaW)*Cos(ThetaWp) - 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp))*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt2,j1)) - 2*(5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt2,2 + j1)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::SHI0, typename fields::conj<fields::SHI0>::type, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*((9.486832980505138*gN*Cos(ThetaWp) + 2*(5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt2,j1)) + 2*(-3.1622776601683795*gN*Cos(ThetaWp) + (5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt2,2 + j1)));

   return {result, minuend_index, subtrahend_index};
}

ScalarVertex VertexImpl<fields::Ah, fields::SHIp, typename fields::conj<fields::SHIp>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
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

MomentumDifferenceVertex VertexImpl<fields::SHIp, typename fields::conj<fields::SHIp>::type, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*KroneckerDelta(gt1,gt2)*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::SHIp, typename fields::conj<fields::SHIp>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05*((10*g2*Cos(ThetaW)*Cos(ThetaWp) - 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp))*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt2,j1)) + 2*(5*g2*Cos(ThetaW)*Cos(ThetaWp) - 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) - 3.1622776601683795*gN*Sin(ThetaWp))*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt2,2 + j1)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::SHIp, typename fields::conj<fields::SHIp>::type, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*((9.486832980505138*gN*Cos(ThetaWp) + 2*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt2,j1)) - 2*(3.1622776601683795*gN*Cos(ThetaWp) + (5*g2*Cos(ThetaW) - 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt2,2 + j1)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::SHp0, typename fields::conj<fields::SHp0>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -0.1*(5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*(Conj(UHp0(gt1,0))*UHp0(gt2,0) + Conj(UHp0(gt1,1))*UHp0(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::SHp0, typename fields::conj<fields::SHp0>::type, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.1*(3.1622776601683795*gN*Cos(ThetaWp) - (5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*(Conj(UHp0(gt1,0))*UHp0(gt2,0) + Conj(UHp0(gt1,1))*UHp0(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::SHpp, typename fields::conj<fields::SHpp>::type, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*(Conj(UHpp(gt1,0))*UHpp(gt2,0) + Conj(UHpp(gt1,1))*UHpp(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::SHpp, typename fields::conj<fields::SHpp>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.1*(5*g2*Cos(ThetaW)*Cos(ThetaWp) - 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) - 3.1622776601683795*gN*Sin(ThetaWp))*(Conj(UHpp(gt1,0))*UHpp(gt2,0) + Conj(UHpp(gt1,1))*UHpp(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::SHpp, typename fields::conj<fields::SHpp>::type, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.1*(3.1622776601683795*gN*Cos(ThetaWp) + (5*g2*Cos(ThetaW) - 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*(Conj(UHpp(gt1,0))*UHpp(gt2,0) + Conj(UHpp(gt1,1))*UHpp(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::SSI0, typename fields::conj<fields::SSI0>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -0.7905694150420949*gN*KroneckerDelta(gt1,gt2)*Sin(ThetaWp);

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::SSI0, typename fields::conj<fields::SSI0>::type, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -0.7905694150420949*gN*Cos(ThetaWp)*KroneckerDelta(gt1,gt2);

   return {result, minuend_index, subtrahend_index};
}

ScalarVertex VertexImpl<fields::Ah, fields::Su, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
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
   const auto gN = MODELPARAMETER(gN);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.016666666666666666*(-((30*g2*Cos(ThetaW)*Cos(ThetaWp) - 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1))) + (30.983866769659336*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.016666666666666666*((-9.486832980505138*gN*Cos(ThetaWp) + 30*g2*Cos(ThetaW)*Sin(ThetaWp) - 7.745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1)) + 2.23606797749979*(4.242640687119286*gN*Cos(ThetaWp) - 13.856406460551018*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1)));

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
   const auto gN = MODELPARAMETER(gN);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -0.1*KroneckerDelta(gt1,gt2)*(5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::Sv, typename fields::conj<fields::Sv>::type, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.1*KroneckerDelta(gt1,gt2)*(3.1622776601683795*gN*Cos(ThetaWp) - (5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp));

   return {result, minuend_index, subtrahend_index};
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

   const std::complex<double> result = std::complex<double>(0,-0.5)*g2*(ZA(gt1,0)*ZP(gt2,0) + ZA(gt1,1)*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VP>::evaluate(
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

InverseMetricVertex VertexImpl<fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VZ>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*g2*(vd*(7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp))*ZP(gt1,0) - 2*vu*(3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*ZP(gt1,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VZp>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.05*g2*(vd*(9.486832980505138*gN*Cos(ThetaWp) + 7.745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt1,0) + 2*vu*(3.1622776601683795*gN*Cos(ThetaWp) - 3.872983346207417*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt1,1));

   return {result};
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
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -(g2*Cos(ThetaW)*Cos(ThetaWp));

   return {result, TripleVectorVertex::odd_permutation{}};
}

TripleVectorVertex VertexImpl<typename fields::conj<fields::VWm>::type, fields::VWm, fields::VZp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = g2*Cos(ThetaW)*Sin(ThetaWp);

   return {result, TripleVectorVertex::odd_permutation{}};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha>::type, fields::Fe, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
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

ChiralVertex VertexImpl<fields::Chi, fields::Fv, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZV = MODELPARAMETER(ZV);

   const std::complex<double> left = IF(gt2 < 3,0.5477225575051661*g1*Conj(ZN(gt1,0))*ZV(gt3,gt2),0) + IF(gt2 < 3,-0.7071067811865475*g2*Conj(ZN(gt1,1))*ZV(gt3,gt2),0) + IF(gt2 < 3,-0.4472135954999579*gN*Conj(ZN(gt1,5))*ZV(gt3,gt2),0);

   const std::complex<double> right = 0;

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fv>::type, fields::Chi, fields::Sv>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZN = MODELPARAMETER(ZN);

   const std::complex<double> left = 0;

   const std::complex<double> right = IF(gt1 < 3,0.5477225575051661*g1*Conj(ZV(gt3,gt1))*ZN(gt2,0),0) + IF(gt1 < 3,-0.7071067811865475*g2*Conj(ZV(gt3,gt1))*ZN(gt2,1),0) + IF(gt1 < 3,-0.4472135954999579*gN*Conj(ZV(gt3,gt1))*ZN(gt2,5),0);

   return {left, right};
}

ScalarVertex VertexImpl<fields::Se, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto TYe = MODELPARAMETER(TYe);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto vs = MODELPARAMETER(vs);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*Conj(TYe(j1,j1))*ZV(gt3,j1))*ZP(gt2,0) - 0.35355339059327373*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZV(gt3,j1))*(vd*ZP(gt2,0) + vu*ZP(gt2,1)) + 0.7071067811865475*(vd*SUM(j1,0,2,AbsSqr(Ye(j1,j1))*Conj(ZE(gt1,j1))*ZV(gt3,j1))*ZP(gt2,0) + vs*Lambdax*SUM(j1,0,2,Conj(Ye(j1,j1))*Conj(ZE(gt1,3 + j1))*ZV(gt3,j1))*ZP(gt2,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::Sv, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto TYe = MODELPARAMETER(TYe);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto vs = MODELPARAMETER(vs);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-1.4142135623730951*Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZE(gt3,j1))*(vd*ZP(gt1,0) + vu*ZP(gt1,1)) + 2*(1.4142135623730951*vd*SUM(j1,0,2,AbsSqr(Ye(j1,j1))*Conj(ZV(gt2,j1))*ZE(gt3,j1))*ZP(gt1,0) + 2*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZE(gt3,3 + j1)*TYe(j1,j1))*ZP(gt1,0) + 1.4142135623730951*vs*Conj(Lambdax)*SUM(j1,0,2,Conj(ZV(gt2,j1))*Ye(j1,j1)*ZE(gt3,3 + j1))*ZP(gt1,1)));

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

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Sv, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = -0.05*KroneckerDelta(gt3,gt4)*((3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZA(gt1,0)*ZA(gt2,0) - (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZA(gt1,1)*ZA(gt2,1) + 5*Sqr(gN)*ZA(gt1,2)*ZA(gt2,2));

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
   const auto gN = MODELPARAMETER(gN);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -(SUM(j1,0,2,AbsSqr(Ye(j1,j1))*Conj(ZV(gt2,j1))*ZV(gt4,j1))*ZP(gt1,0)*ZP(gt3,0)) + 0.05*KroneckerDelta(gt2,gt4)*((-3*Sqr(g1) + 5*Sqr(g2) + 3*Sqr(gN))*ZP(gt1,0)*ZP(gt3,0) + (3*Sqr(g1) - 5*Sqr(g2) + 2*Sqr(gN))*ZP(gt1,1)*ZP(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::SDX, fields::Sv, typename fields::conj<fields::SDX>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZDX = MODELPARAMETER(ZDX);

   const std::complex<double> result = 0.05*KroneckerDelta(gt2,gt4)*(-2*(Sqr(g1) - Sqr(gN))*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt3,j1)) + (2*Sqr(g1) + 3*Sqr(gN))*SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt3,3 + j1)));

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
   const auto gN = MODELPARAMETER(gN);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);

   const std::complex<double> result = 0.05*(-(KroneckerDelta(gt2,gt4)*((3*Sqr(g1) - 5*Sqr(g2) + 2*Sqr(gN))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1)) + (-6*Sqr(g1) + Sqr(gN))*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1)))) - 5*(Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZV(gt4,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZE(gt3,j2)) + Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZV(gt4,j2)) + 4*SUM(j1,0,2,Conj(ZV(gt2,j1))*Ye(j1,j1)*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(Ye(j2,j2))*Conj(ZE(gt1,3 + j2))*ZV(gt4,j2))));

   return {result};
}

ScalarVertex VertexImpl<fields::SHI0, fields::Sv, typename fields::conj<fields::SHI0>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto UHI0 = MODELPARAMETER(UHI0);

   const std::complex<double> result = -0.05*KroneckerDelta(gt2,gt4)*((3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt3,j1)) - (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt3,2 + j1)));

   return {result};
}

ScalarVertex VertexImpl<fields::SHIp, fields::Sv, typename fields::conj<fields::SHIp>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto UHIp = MODELPARAMETER(UHIp);

   const std::complex<double> result = -0.05*KroneckerDelta(gt2,gt4)*((3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt3,j1)) + (-3*Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN))*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt3,2 + j1)));

   return {result};
}

ScalarVertex VertexImpl<fields::SHp0, fields::Sv, typename fields::conj<fields::SHp0>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto UHp0 = MODELPARAMETER(UHp0);

   const std::complex<double> result = -0.05*KroneckerDelta(gt2,gt4)*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*(Conj(UHp0(gt1,0))*UHp0(gt3,0) - Conj(UHp0(gt1,1))*UHp0(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::SHpp, fields::Sv, typename fields::conj<fields::SHpp>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto UHpp = MODELPARAMETER(UHpp);

   const std::complex<double> result = -0.05*KroneckerDelta(gt2,gt4)*(3*Sqr(g1) - 5*Sqr(g2) + 2*Sqr(gN))*(Conj(UHpp(gt1,0))*UHpp(gt3,0) - Conj(UHpp(gt1,1))*UHpp(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::SSI0, fields::Sv, typename fields::conj<fields::SSI0>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto gN = MODELPARAMETER(gN);

   const std::complex<double> result = -0.25*KroneckerDelta(gt1,gt3)*KroneckerDelta(gt2,gt4)*Sqr(gN);

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
   const auto gN = MODELPARAMETER(gN);
   const auto ZU = MODELPARAMETER(ZU);

   const std::complex<double> result = 0.05*KroneckerDelta(gt2,gt4)*((Sqr(g1) - 5*Sqr(g2) - Sqr(gN))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt3,j1)) - (4*Sqr(g1) + Sqr(gN))*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt3,3 + j1)));

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
   const auto gN = MODELPARAMETER(gN);

   const std::complex<double> result = -0.05*(KroneckerDelta(gt1,gt4)*KroneckerDelta(gt2,gt3) + KroneckerDelta(gt1,gt3)*KroneckerDelta(gt2,gt4))*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN));

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
   const auto gN = MODELPARAMETER(gN);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.1*KroneckerDelta(gt1,gt2)*(3.1622776601683795*g2*gN*Cos(ThetaW)*Sin(2*ThetaWp) + g1*Sin(ThetaW)*(7.745966692414834*g2*Cos(ThetaW) + 3*g1*Sin(ThetaW))*Sqr(Cos(ThetaWp)) + 5*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Cos(ThetaWp)) + gN*(2.449489742783178*g1*Sin(ThetaW)*Sin(2*ThetaWp) + 2*gN*Sqr(Sin(ThetaWp))));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sv, typename fields::conj<fields::Sv>::type, fields::VZp, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.1*KroneckerDelta(gt1,gt2)*((g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*Sin(ThetaWp)*(-6.324555320336759*gN*Cos(ThetaWp) + (5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp)) + 2*Sqr(gN)*Sqr(Cos(ThetaWp)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sv, typename fields::conj<fields::Sv>::type, fields::VZ, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.0125*KroneckerDelta(gt1,gt2)*(25.298221281347036*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) - 15.491933384829668*g1*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) - 9*Sin(2*ThetaWp)*Sqr(g1) + 3*Cos(2*ThetaW)*Sin(2*ThetaWp)*Sqr(g1) - 5*Sin(2*ThetaWp)*Sqr(g2) - 5*Cos(2*ThetaW)*Sin(2*ThetaWp)*Sqr(g2) + 8*Sin(2*ThetaWp)*Sqr(gN) + 2*Sin(2*ThetaWp)*(3*Sqr(g1) - 5*Sqr(g2))*Sqr(Cos(ThetaW)) + 19.595917942265423*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) - 19.595917942265423*g1*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp)));

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
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -0.35355339059327373*(Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt3,j1))*ZE(gt4,j1))*(ZH(gt1,0)*ZP(gt2,0) + ZH(gt1,1)*ZP(gt2,1)) - 2*(SUM(j1,0,2,AbsSqr(Ye(j1,j1))*Conj(ZV(gt3,j1))*ZE(gt4,j1))*ZH(gt1,0)*ZP(gt2,0) + Conj(Lambdax)*SUM(j1,0,2,Conj(ZV(gt3,j1))*Ye(j1,j1)*ZE(gt4,3 + j1))*ZH(gt1,2)*ZP(gt2,1)));

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
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -0.35355339059327373*(Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZV(gt4,j1))*(ZH(gt1,0)*ZP(gt3,0) + ZH(gt1,1)*ZP(gt3,1)) - 2*(SUM(j1,0,2,AbsSqr(Ye(j1,j1))*Conj(ZE(gt2,j1))*ZV(gt4,j1))*ZH(gt1,0)*ZP(gt3,0) + Lambdax*SUM(j1,0,2,Conj(Ye(j1,j1))*Conj(ZE(gt2,3 + j1))*ZV(gt4,j1))*ZH(gt1,2)*ZP(gt3,1)));

   return {result};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha>::type, fields::Fd, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto UP = MODELPARAMETER(UP);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZDR = MODELPARAMETER(ZDR);
   const auto UM = MODELPARAMETER(UM);

   const std::complex<double> left = -(g2*Conj(UP(gt1,0))*SUM(j1,0,2,Conj(ZDL(gt2,j1))*ZU(gt3,j1))) + Conj(UP(gt1,1))*SUM(j1,0,2,Conj(ZDL(gt2,j1))*Yu(j1,j1)*ZU(gt3,3 + j1));

   const std::complex<double> right = SUM(j1,0,2,Conj(Yd(j1,j1))*ZDR(gt2,j1)*ZU(gt3,j1))*UM(gt1,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Cha, fields::Su>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Yu = MODELPARAMETER(Yu);
   const auto UM = MODELPARAMETER(UM);
   const auto ZDR = MODELPARAMETER(ZDR);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto UP = MODELPARAMETER(UP);

   const std::complex<double> left = Conj(UM(gt2,1))*SUM(j1,0,2,Conj(ZDR(gt1,j1))*Conj(ZU(gt3,j1))*Yd(j1,j1));

   const std::complex<double> right = -(g2*SUM(j1,0,2,Conj(ZU(gt3,j1))*ZDL(gt1,j1))*UP(gt2,0)) + SUM(j1,0,2,Conj(Yu(j1,j1))*Conj(ZU(gt3,3 + j1))*ZDL(gt1,j1))*UP(gt2,1);

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
   const auto gN = MODELPARAMETER(gN);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZUR = MODELPARAMETER(ZUR);

   const std::complex<double> left = -0.18257418583505536*g1*Conj(ZN(gt1,0))*SUM(j1,0,2,Conj(ZUL(gt2,j1))*ZU(gt3,j1)) - 0.7071067811865475*g2*Conj(ZN(gt1,1))*SUM(j1,0,2,Conj(ZUL(gt2,j1))*ZU(gt3,j1)) - 0.22360679774997896*gN*Conj(ZN(gt1,5))*SUM(j1,0,2,Conj(ZUL(gt2,j1))*ZU(gt3,j1)) - Conj(ZN(gt1,3))*SUM(j1,0,2,Conj(ZUL(gt2,j1))*Yu(j1,j1)*ZU(gt3,3 + j1));

   const std::complex<double> right = -(SUM(j1,0,2,Conj(Yu(j1,j1))*ZU(gt3,j1)*ZUR(gt2,j1))*ZN(gt1,3)) + 0.07453559924999298*SUM(j1,0,2,ZU(gt3,3 + j1)*ZUR(gt2,j1))*(9.797958971132712*g1*ZN(gt1,0) - 3*gN*ZN(gt1,5));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, fields::Chi, fields::Su>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZUR = MODELPARAMETER(ZUR);
   const auto ZUL = MODELPARAMETER(ZUL);

   const std::complex<double> left = 0.7302967433402214*g1*Conj(ZN(gt2,0))*SUM(j1,0,2,Conj(ZU(gt3,3 + j1))*Conj(ZUR(gt1,j1))) - 0.22360679774997896*gN*Conj(ZN(gt2,5))*SUM(j1,0,2,Conj(ZU(gt3,3 + j1))*Conj(ZUR(gt1,j1))) - Conj(ZN(gt2,3))*SUM(j1,0,2,Conj(ZU(gt3,j1))*Conj(ZUR(gt1,j1))*Yu(j1,j1));

   const std::complex<double> right = -(SUM(j1,0,2,Conj(Yu(j1,j1))*Conj(ZU(gt3,3 + j1))*ZUL(gt1,j1))*ZN(gt2,3)) - 0.03333333333333333*SUM(j1,0,2,Conj(ZU(gt3,j1))*ZUL(gt1,j1))*(5.477225575051661*g1*ZN(gt2,0) + 21.213203435596427*g2*ZN(gt2,1) + 6.708203932499369*gN*ZN(gt2,5));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Glu, fields::Fu, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt3 = indices[1];
   const auto g3 = MODELPARAMETER(g3);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZUR = MODELPARAMETER(ZUR);
   const auto PhaseGlu = PHASE(PhaseGlu);

   const std::complex<double> left = -0.7071067811865475*g3*PhaseGlu*SUM(j1,0,2,Conj(ZUL(gt2,j1))*ZU(gt3,j1));

   const std::complex<double> right = 0.7071067811865475*g3*Conj(PhaseGlu)*SUM(j1,0,2,ZU(gt3,3 + j1)*ZUR(gt2,j1));

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
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto PhaseGlu = PHASE(PhaseGlu);

   const std::complex<double> left = 0.7071067811865475*g3*PhaseGlu*SUM(j1,0,2,Conj(ZU(gt3,3 + j1))*Conj(ZUR(gt1,j1)));

   const std::complex<double> right = -0.7071067811865475*g3*Conj(PhaseGlu)*SUM(j1,0,2,Conj(ZU(gt3,j1))*ZUL(gt1,j1));

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
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = 0.025*(SUM(j1,0,2,Conj(ZU(gt3,3 + j1))*ZU(gt4,3 + j1))*((-8*Sqr(g1) + 3*Sqr(gN))*ZA(gt1,0)*ZA(gt2,0) + 2*(4*Sqr(g1) + Sqr(gN))*ZA(gt1,1)*ZA(gt2,1) - 5*Sqr(gN)*ZA(gt1,2)*ZA(gt2,2)) + SUM(j1,0,2,Conj(ZU(gt3,j1))*ZU(gt4,j1))*((2*Sqr(g1) - 10*Sqr(g2) + 3*Sqr(gN))*ZA(gt1,0)*ZA(gt2,0) + 2*(-Sqr(g1) + 5*Sqr(g2) + Sqr(gN))*ZA(gt1,1)*ZA(gt2,1) - 5*Sqr(gN)*ZA(gt1,2)*ZA(gt2,2)) - 20*(2*(SUM(j1,0,2,AbsSqr(Yu(j1,j1))*Conj(ZU(gt3,j1))*ZU(gt4,j1)) + SUM(j1,0,2,AbsSqr(Yu(j1,j1))*Conj(ZU(gt3,3 + j1))*ZU(gt4,3 + j1)))*ZA(gt1,1)*ZA(gt2,1) + Lambdax*SUM(j1,0,2,Conj(Yu(j1,j1))*Conj(ZU(gt3,3 + j1))*ZU(gt4,j1))*(ZA(gt1,2)*ZA(gt2,0) + ZA(gt1,0)*ZA(gt2,2)) + Conj(Lambdax)*SUM(j1,0,2,Conj(ZU(gt3,j1))*Yu(j1,j1)*ZU(gt4,3 + j1))*(ZA(gt1,2)*ZA(gt2,0) + ZA(gt1,0)*ZA(gt2,2))));

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
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.025*(-40*SUM(j1,0,2,AbsSqr(Yd(j1,j1))*Conj(ZU(gt2,j1))*ZU(gt4,j1))*ZP(gt1,0)*ZP(gt3,0) - 8*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt4,3 + j1))*ZP(gt1,0)*ZP(gt3,0) + 3*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt4,3 + j1))*ZP(gt1,0)*ZP(gt3,0) + 8*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt4,3 + j1))*ZP(gt1,1)*ZP(gt3,1) + 2*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt4,3 + j1))*ZP(gt1,1)*ZP(gt3,1) - 40*SUM(j1,0,2,AbsSqr(Yu(j1,j1))*Conj(ZU(gt2,3 + j1))*ZU(gt4,3 + j1))*ZP(gt1,1)*ZP(gt3,1) + SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt4,j1))*((2*Sqr(g1) + 10*Sqr(g2) + 3*Sqr(gN))*ZP(gt1,0)*ZP(gt3,0) - 2*(Sqr(g1) + 5*Sqr(g2) - Sqr(gN))*ZP(gt1,1)*ZP(gt3,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::SDX, fields::Su, typename fields::conj<fields::SDX>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g3 = MODELPARAMETER(g3);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZDX = MODELPARAMETER(ZDX);

   const std::complex<double> result = -0.25*Sqr(g3)*(SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt4,j1))*(SUM(j2,0,2,Conj(ZDX(gt1,j2))*ZDX(gt3,j2)) - SUM(j2,0,2,Conj(ZDX(gt1,3 + j2))*ZDX(gt3,3 + j2))) + SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt4,3 + j1))*(-SUM(j2,0,2,Conj(ZDX(gt1,j2))*ZDX(gt3,j2)) + SUM(j2,0,2,Conj(ZDX(gt1,3 + j2))*ZDX(gt3,3 + j2))) + (SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt3,j1)) - SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt3,3 + j1)))*(SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) - SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2))));

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
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZE = MODELPARAMETER(ZE);

   const std::complex<double> result = 0.0125*(-(SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt4,3 + j1))*(2*(4*Sqr(g1) + Sqr(gN))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt3,j2)) + (-16*Sqr(g1) + Sqr(gN))*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt3,3 + j2)))) + SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt4,j1))*(2*(Sqr(g1) + 5*Sqr(g2) - Sqr(gN))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt3,j2)) - (4*Sqr(g1) + Sqr(gN))*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt3,3 + j2))) + 2*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) + 10*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) - 2*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) - Sqr(gN)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) - 8*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)) - 2*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)) + 16*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)) - Sqr(gN)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)));

   return {result};
}

ScalarVertex VertexImpl<fields::SHI0, fields::Su, typename fields::conj<fields::SHI0>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZU = MODELPARAMETER(ZU);
   const auto UHI0 = MODELPARAMETER(UHI0);

   const std::complex<double> result = 0.0125*(SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt4,3 + j1))*((-8*Sqr(g1) + 3*Sqr(gN))*SUM(j2,0,1,Conj(UHI0(gt1,j2))*UHI0(gt3,j2)) + 2*(4*Sqr(g1) + Sqr(gN))*SUM(j2,0,1,Conj(UHI0(gt1,2 + j2))*UHI0(gt3,2 + j2))) + SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt4,j1))*((2*Sqr(g1) - 10*Sqr(g2) + 3*Sqr(gN))*SUM(j2,0,1,Conj(UHI0(gt1,j2))*UHI0(gt3,j2)) + 2*(-Sqr(g1) + 5*Sqr(g2) + Sqr(gN))*SUM(j2,0,1,Conj(UHI0(gt1,2 + j2))*UHI0(gt3,2 + j2))) + 2*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) - 10*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) + 3*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) - 2*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt3,2 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) + 10*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt3,2 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) + 2*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt3,2 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) - 8*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)) + 3*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)) + 8*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt3,2 + j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)) + 2*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt3,2 + j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)));

   return {result};
}

ScalarVertex VertexImpl<fields::SHIp, fields::Su, typename fields::conj<fields::SHIp>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto ZU = MODELPARAMETER(ZU);
   const auto UHIp = MODELPARAMETER(UHIp);

   const std::complex<double> result = 0.0125*(SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt4,j1))*((2*Sqr(g1) + 10*Sqr(g2) + 3*Sqr(gN))*SUM(j2,0,1,Conj(UHIp(gt1,j2))*UHIp(gt3,j2)) - 2*(Sqr(g1) + 5*Sqr(g2) - Sqr(gN))*SUM(j2,0,1,Conj(UHIp(gt1,2 + j2))*UHIp(gt3,2 + j2))) + SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt4,3 + j1))*((-8*Sqr(g1) + 3*Sqr(gN))*SUM(j2,0,1,Conj(UHIp(gt1,j2))*UHIp(gt3,j2)) + 2*(4*Sqr(g1) + Sqr(gN))*SUM(j2,0,1,Conj(UHIp(gt1,2 + j2))*UHIp(gt3,2 + j2))) + 2*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) + 10*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) + 3*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) - 2*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt3,2 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) - 10*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt3,2 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) + 2*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt3,2 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) - 8*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)) + 3*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)) + 8*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt3,2 + j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)) + 2*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt3,2 + j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)));

   return {result};
}

ScalarVertex VertexImpl<fields::SHp0, fields::Su, typename fields::conj<fields::SHp0>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto ZU = MODELPARAMETER(ZU);
   const auto UHp0 = MODELPARAMETER(UHp0);

   const std::complex<double> result = 0.05*((Sqr(g1) - 5*Sqr(g2) - Sqr(gN))*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt4,j1)) - (4*Sqr(g1) + Sqr(gN))*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt4,3 + j1)))*(Conj(UHp0(gt1,0))*UHp0(gt3,0) - Conj(UHp0(gt1,1))*UHp0(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::SHpp, fields::Su, typename fields::conj<fields::SHpp>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto ZU = MODELPARAMETER(ZU);
   const auto UHpp = MODELPARAMETER(UHpp);

   const std::complex<double> result = 0.05*((Sqr(g1) + 5*Sqr(g2) - Sqr(gN))*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt4,j1)) - (4*Sqr(g1) + Sqr(gN))*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt4,3 + j1)))*(Conj(UHpp(gt1,0))*UHpp(gt3,0) - Conj(UHpp(gt1,1))*UHpp(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::SSI0, fields::Su, typename fields::conj<fields::SSI0>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto gN = MODELPARAMETER(gN);

   const std::complex<double> result = -0.125*KroneckerDelta(gt1,gt3)*KroneckerDelta(gt2,gt4)*Sqr(gN);

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
   const auto gN = MODELPARAMETER(gN);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);

   const std::complex<double> result = 0.004166666666666667*(-2*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt4,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt3,j2)) - 30*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt4,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt3,j2)) + 20*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt4,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt3,j2)) - 3*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt4,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt3,j2)) + 8*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt4,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt3,j2)) - 20*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt4,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt3,j2)) - 3*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt4,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt3,j2)) - 240*SUM(j1,0,2,Conj(ZU(gt1,j1))*Yu(j1,j1)*ZU(gt4,3 + j1))*SUM(j2,0,2,Conj(Yu(j2,j2))*Conj(ZU(gt2,3 + j2))*ZU(gt3,j2)) - 60*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt4,j1))*(SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt3,j2)) - SUM(j2,0,2,Conj(ZU(gt1,3 + j2))*ZU(gt3,3 + j2))) + 60*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt4,3 + j1))*(SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt3,j2)) - SUM(j2,0,2,Conj(ZU(gt1,3 + j2))*ZU(gt3,3 + j2))) + 8*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt4,j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt3,3 + j2)) - 20*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt4,j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt3,3 + j2)) - 3*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt4,j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt3,3 + j2)) - 32*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt4,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt3,3 + j2)) + 20*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt4,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt3,3 + j2)) - 3*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt4,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt3,3 + j2)) - 2*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt4,j2)) - 30*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt4,j2)) + 20*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt4,j2)) - 3*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt4,j2)) + 8*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt4,j2)) - 20*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt4,j2)) - 3*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt4,j2)) - 240*SUM(j1,0,2,Conj(ZU(gt2,j1))*Yu(j1,j1)*ZU(gt3,3 + j1))*SUM(j2,0,2,Conj(Yu(j2,j2))*Conj(ZU(gt1,3 + j2))*ZU(gt4,j2)) - 60*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) + 60*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) + 8*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt1,3 + j2))*ZU(gt4,3 + j2)) - 20*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt1,3 + j2))*ZU(gt4,3 + j2)) - 3*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt1,3 + j2))*ZU(gt4,3 + j2)) - 32*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt1,3 + j2))*ZU(gt4,3 + j2)) + 20*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt1,3 + j2))*ZU(gt4,3 + j2)) - 3*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt1,3 + j2))*ZU(gt4,3 + j2)) + 60*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)) - 60*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)));

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
   const auto gN = MODELPARAMETER(gN);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.016666666666666666*((-2.449489742783178*g1*gN*Sin(ThetaW)*Sin(2*ThetaWp) + 30*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Cos(ThetaWp)) + Cos(ThetaW)*(9.486832980505138*g2*gN*Sin(2*ThetaWp) - 15.491933384829668*g1*g2*Sin(ThetaW)*Sqr(Cos(ThetaWp))) + 2*Sqr(g1)*Sqr(Cos(ThetaWp))*Sqr(Sin(ThetaW)) + 3*Sqr(gN)*Sqr(Sin(ThetaWp)))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1)) + (19.595917942265423*g1*gN*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) + 32*Sqr(g1)*Sqr(Cos(ThetaWp))*Sqr(Sin(ThetaW)) + 3*Sqr(gN)*Sqr(Sin(ThetaWp)))*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VZp, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.016666666666666666*((3*Sqr(gN)*Sqr(Cos(ThetaWp)) + 30*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Sin(ThetaWp)) + g1*Sin(ThetaW)*(2.449489742783178*gN*Sin(2*ThetaWp) + 2*g1*Sin(ThetaW)*Sqr(Sin(ThetaWp))) - Cos(ThetaW)*(9.486832980505138*g2*gN*Sin(2*ThetaWp) + 15.491933384829668*g1*g2*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1)) + (-19.595917942265423*g1*gN*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) + 3*Sqr(gN)*Sqr(Cos(ThetaWp)) + 32*Sqr(g1)*Sqr(Sin(ThetaW))*Sqr(Sin(ThetaWp)))*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VZ, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.016666666666666666*((9.486832980505138*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) + 3*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN) - 30*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) - 2.449489742783178*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + g1*(3.872983346207417*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) - g1*Sin(2*ThetaWp)*Sqr(Sin(ThetaW)) + 2.449489742783178*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1)) + (3*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN) + 9.797958971132712*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) - 4*g1*Sin(ThetaW)*(4*g1*Sin(ThetaW)*Sin(2*ThetaWp) + 2.449489742783178*gN*Sqr(Sin(ThetaWp))))*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1)));

   return {result};
}

ChiralVertex VertexImpl<fields::Chi, fields::Fe, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Chi, fields::Se>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
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

ChiralVertex VertexImpl<fields::Cha, fields::Fv, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto UM = MODELPARAMETER(UM);
   const auto ZE = MODELPARAMETER(ZE);

   const std::complex<double> left = IF(gt2 < 3,-(g2*Conj(UM(gt1,0))*ZE(gt3,gt2)),0) + IF(gt2 < 3,Conj(UM(gt1,1))*Ye(gt2,gt2)*ZE(gt3,3 + gt2),0);

   const std::complex<double> right = 0;

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha>::type, typename fields::bar<fields::Fv>::type, fields::Se>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto UM = MODELPARAMETER(UM);

   const std::complex<double> left = 0;

   const std::complex<double> right = IF(gt2 < 3,-(g2*Conj(ZE(gt3,gt2))*UM(gt1,0)),0) + IF(gt2 < 3,Conj(Ye(gt2,gt2))*Conj(ZE(gt3,3 + gt2))*UM(gt1,1),0);

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
   const auto gN = MODELPARAMETER(gN);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = 0.025*(-40*SUM(j1,0,2,AbsSqr(Ye(j1,j1))*Conj(ZE(gt3,j1))*ZE(gt4,j1))*ZA(gt1,0)*ZA(gt2,0) + 12*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*ZE(gt4,3 + j1))*ZA(gt1,0)*ZA(gt2,0) + 3*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*ZE(gt4,3 + j1))*ZA(gt1,0)*ZA(gt2,0) - 40*SUM(j1,0,2,AbsSqr(Ye(j1,j1))*Conj(ZE(gt3,3 + j1))*ZE(gt4,3 + j1))*ZA(gt1,0)*ZA(gt2,0) - 12*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*ZE(gt4,3 + j1))*ZA(gt1,1)*ZA(gt2,1) + 2*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*ZE(gt4,3 + j1))*ZA(gt1,1)*ZA(gt2,1) - 20*Lambdax*SUM(j1,0,2,Conj(Ye(j1,j1))*Conj(ZE(gt3,3 + j1))*ZE(gt4,j1))*ZA(gt1,2)*ZA(gt2,1) - 20*Conj(Lambdax)*SUM(j1,0,2,Conj(ZE(gt3,j1))*Ye(j1,j1)*ZE(gt4,3 + j1))*ZA(gt1,2)*ZA(gt2,1) - 20*Lambdax*SUM(j1,0,2,Conj(Ye(j1,j1))*Conj(ZE(gt3,3 + j1))*ZE(gt4,j1))*ZA(gt1,1)*ZA(gt2,2) - 20*Conj(Lambdax)*SUM(j1,0,2,Conj(ZE(gt3,j1))*Ye(j1,j1)*ZE(gt4,3 + j1))*ZA(gt1,1)*ZA(gt2,2) - 5*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*ZE(gt4,3 + j1))*ZA(gt1,2)*ZA(gt2,2) - 2*SUM(j1,0,2,Conj(ZE(gt3,j1))*ZE(gt4,j1))*((3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZA(gt1,0)*ZA(gt2,0) + (-3*Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN))*ZA(gt1,1)*ZA(gt2,1) + 5*Sqr(gN)*ZA(gt1,2)*ZA(gt2,2)));

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
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.025*(-40*SUM(j1,0,2,AbsSqr(Ye(j1,j1))*Conj(ZE(gt2,3 + j1))*ZE(gt4,3 + j1))*ZP(gt1,0)*ZP(gt3,0) + SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt4,3 + j1))*(3*(4*Sqr(g1) + Sqr(gN))*ZP(gt1,0)*ZP(gt3,0) + 2*(-6*Sqr(g1) + Sqr(gN))*ZP(gt1,1)*ZP(gt3,1)) + SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt4,j1))*(-2*(3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZP(gt1,0)*ZP(gt3,0) + 2*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZP(gt1,1)*ZP(gt3,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::SDX, fields::Se, typename fields::conj<fields::SDX>::type, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZDX = MODELPARAMETER(ZDX);

   const std::complex<double> result = 0.0125*(SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt4,3 + j1))*(2*(4*Sqr(g1) + Sqr(gN))*SUM(j2,0,2,Conj(ZDX(gt1,j2))*ZDX(gt3,j2)) + (-8*Sqr(g1) + 3*Sqr(gN))*SUM(j2,0,2,Conj(ZDX(gt1,3 + j2))*ZDX(gt3,3 + j2))) + SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt4,j1))*(-4*(Sqr(g1) - Sqr(gN))*SUM(j2,0,2,Conj(ZDX(gt1,j2))*ZDX(gt3,j2)) + 2*(2*Sqr(g1) + 3*Sqr(gN))*SUM(j2,0,2,Conj(ZDX(gt1,3 + j2))*ZDX(gt3,3 + j2))) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) + 4*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) + 4*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) + 6*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) + 8*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt4,3 + j2)) + 2*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt4,3 + j2)) - 8*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt4,3 + j2)) + 3*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt4,3 + j2)));

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
   const auto gN = MODELPARAMETER(gN);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);

   const std::complex<double> result = 0.0125*(-80*SUM(j1,0,2,Conj(ZE(gt2,j1))*Ye(j1,j1)*ZE(gt4,3 + j1))*SUM(j2,0,2,Conj(Ye(j2,j2))*Conj(ZE(gt1,3 + j2))*ZE(gt3,j2)) - 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt4,j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt3,j2)) - 10*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt4,j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt3,j2)) - 4*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt4,j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt3,j2)) + 12*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt4,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt3,j2)) - 2*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt4,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt3,j2)) - 80*SUM(j1,0,2,Conj(ZE(gt1,j1))*Ye(j1,j1)*ZE(gt4,3 + j1))*SUM(j2,0,2,Conj(Ye(j2,j2))*Conj(ZE(gt2,3 + j2))*ZE(gt3,j2)) - 2*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt4,j1))*((3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt3,j2)) + (-6*Sqr(g1) + Sqr(gN))*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt3,3 + j2))) - SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt4,3 + j1))*(2*(-6*Sqr(g1) + Sqr(gN))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt3,j2)) + (24*Sqr(g1) + Sqr(gN))*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt3,3 + j2))) + 12*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt4,j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt3,3 + j2)) - 2*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt4,j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt3,3 + j2)) - 24*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt4,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt3,3 + j2)) - Sqr(gN)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt4,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt3,3 + j2)) - 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt4,j2)) - 10*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt4,j2)) - 4*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt4,j2)) + 12*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt4,j2)) - 2*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt4,j2)) - 80*SUM(j1,0,2,Conj(ZE(gt2,j1))*Ye(j1,j1)*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(Ye(j2,j2))*Conj(ZE(gt1,3 + j2))*ZE(gt4,j2)) - 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) - 10*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) - 4*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) + 12*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) - 2*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) - 80*SUM(j1,0,2,Conj(ZE(gt1,j1))*Ye(j1,j1)*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(Ye(j2,j2))*Conj(ZE(gt2,3 + j2))*ZE(gt4,j2)) + 12*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt4,3 + j2)) - 2*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt4,3 + j2)) - 24*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt4,3 + j2)) - Sqr(gN)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt4,3 + j2)) + 12*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt4,3 + j2)) - 2*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt4,3 + j2)) - 24*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt4,3 + j2)) - Sqr(gN)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt4,3 + j2)));

   return {result};
}

ScalarVertex VertexImpl<fields::Se, fields::SHI0, typename fields::conj<fields::Se>::type, typename fields::conj<fields::SHI0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto ZE = MODELPARAMETER(ZE);
   const auto UHI0 = MODELPARAMETER(UHI0);

   const std::complex<double> result = 0.0125*(SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1))*(3*(4*Sqr(g1) + Sqr(gN))*SUM(j2,0,1,Conj(UHI0(gt2,j2))*UHI0(gt4,j2)) + 2*(-6*Sqr(g1) + Sqr(gN))*SUM(j2,0,1,Conj(UHI0(gt2,2 + j2))*UHI0(gt4,2 + j2))) + SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*((-6*Sqr(g1) + 10*Sqr(g2) + 6*Sqr(gN))*SUM(j2,0,1,Conj(UHI0(gt2,j2))*UHI0(gt4,j2)) + 2*(3*Sqr(g1) - 5*Sqr(g2) + 2*Sqr(gN))*SUM(j2,0,1,Conj(UHI0(gt2,2 + j2))*UHI0(gt4,2 + j2))) - 6*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt4,j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt3,j2)) + 10*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt4,j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt3,j2)) + 6*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt4,j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt3,j2)) + 6*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*UHI0(gt4,2 + j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt3,j2)) - 10*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*UHI0(gt4,2 + j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt3,j2)) + 4*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*UHI0(gt4,2 + j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt3,j2)) + 12*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt4,j1))*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt3,3 + j2)) + 3*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt4,j1))*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt3,3 + j2)) - 12*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*UHI0(gt4,2 + j1))*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt3,3 + j2)) + 2*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*UHI0(gt4,2 + j1))*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt3,3 + j2)));

   return {result};
}

ScalarVertex VertexImpl<fields::Se, fields::SHIp, typename fields::conj<fields::Se>::type, typename fields::conj<fields::SHIp>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZE = MODELPARAMETER(ZE);
   const auto UHIp = MODELPARAMETER(UHIp);

   const std::complex<double> result = 0.0125*(SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1))*(3*(4*Sqr(g1) + Sqr(gN))*SUM(j2,0,1,Conj(UHIp(gt2,j2))*UHIp(gt4,j2)) + 2*(-6*Sqr(g1) + Sqr(gN))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*UHIp(gt4,2 + j2))) + SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*(-2*(3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*SUM(j2,0,1,Conj(UHIp(gt2,j2))*UHIp(gt4,j2)) + 2*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*UHIp(gt4,2 + j2))) - 6*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt4,j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt3,j2)) - 10*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt4,j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt3,j2)) + 6*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt4,j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt3,j2)) + 6*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*UHIp(gt4,2 + j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt3,j2)) + 10*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*UHIp(gt4,2 + j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt3,j2)) + 4*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*UHIp(gt4,2 + j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt3,j2)) + 12*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt4,j1))*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt3,3 + j2)) + 3*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt4,j1))*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt3,3 + j2)) - 12*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*UHIp(gt4,2 + j1))*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt3,3 + j2)) + 2*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*UHIp(gt4,2 + j1))*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt3,3 + j2)));

   return {result};
}

ScalarVertex VertexImpl<fields::Se, fields::SHp0, typename fields::conj<fields::Se>::type, typename fields::conj<fields::SHp0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto ZE = MODELPARAMETER(ZE);
   const auto UHp0 = MODELPARAMETER(UHp0);

   const std::complex<double> result = -0.05*((3*Sqr(g1) - 5*Sqr(g2) + 2*Sqr(gN))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1)) + (-6*Sqr(g1) + Sqr(gN))*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1)))*(Conj(UHp0(gt2,0))*UHp0(gt4,0) - Conj(UHp0(gt2,1))*UHp0(gt4,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Se, fields::SHpp, typename fields::conj<fields::Se>::type, typename fields::conj<fields::SHpp>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto ZE = MODELPARAMETER(ZE);
   const auto UHpp = MODELPARAMETER(UHpp);

   const std::complex<double> result = -0.05*((3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1)) + (-6*Sqr(g1) + Sqr(gN))*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1)))*(Conj(UHpp(gt2,0))*UHpp(gt4,0) - Conj(UHpp(gt2,1))*UHpp(gt4,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Se, fields::SSI0, typename fields::conj<fields::Se>::type, typename fields::conj<fields::SSI0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto gN = MODELPARAMETER(gN);
   const auto ZE = MODELPARAMETER(ZE);

   const std::complex<double> result = -0.125*KroneckerDelta(gt2,gt4)*Sqr(gN)*(2*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1)) + SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1)));

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
   const auto gN = MODELPARAMETER(gN);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05*(2*(2.449489742783178*g1*gN*Sin(ThetaW)*Sin(2*ThetaWp) + 5*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Cos(ThetaWp)) - 2*Cos(ThetaW)*(3.1622776601683795*g2*gN*Cos(ThetaWp)*Sin(ThetaWp) + 3.872983346207417*g1*g2*Sin(ThetaW)*Sqr(Cos(ThetaWp))) + 3*Sqr(g1)*Sqr(Cos(ThetaWp))*Sqr(Sin(ThetaW)) + 2*Sqr(gN)*Sqr(Sin(ThetaWp)))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) + (-9.797958971132712*g1*gN*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) + 24*Sqr(g1)*Sqr(Cos(ThetaWp))*Sqr(Sin(ThetaW)) + Sqr(gN)*Sqr(Sin(ThetaWp)))*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Se, typename fields::conj<fields::Se>::type, fields::VZp, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*(2*(2*Sqr(gN)*Sqr(Cos(ThetaWp)) + 5*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Sin(ThetaWp)) + g1*Sin(ThetaW)*(-2.449489742783178*gN*Sin(2*ThetaWp) + 3*g1*Sin(ThetaW)*Sqr(Sin(ThetaWp))) + Cos(ThetaW)*(3.1622776601683795*g2*gN*Sin(2*ThetaWp) - 7.745966692414834*g1*g2*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) + (9.797958971132712*g1*gN*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) + Sqr(gN)*Sqr(Cos(ThetaWp)) + 24*Sqr(g1)*Sqr(Sin(ThetaW))*Sqr(Sin(ThetaWp)))*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Se, typename fields::conj<fields::Se>::type, fields::VZ, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05*((-6.324555320336759*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) + 4*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN) - 5*Sin(2*ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 4.898979485566356*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + g1*(3.872983346207417*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) - 3*g1*Sin(2*ThetaWp)*Sqr(Sin(ThetaW)) - 4.898979485566356*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) + (-4.898979485566356*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + Cos(ThetaWp)*Sin(ThetaWp)*(Sqr(gN) - 24*Sqr(g1)*Sqr(Sin(ThetaW))) + 4.898979485566356*g1*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp)))*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1)));

   return {result};
}

ChiralVertex VertexImpl<fields::Chi, fields::FDX, typename fields::conj<fields::SDX>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto Kappa = MODELPARAMETER(Kappa);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZDXL = MODELPARAMETER(ZDXL);
   const auto ZDX = MODELPARAMETER(ZDX);
   const auto ZDXR = MODELPARAMETER(ZDXR);

   const std::complex<double> left = 0.3651483716701107*g1*Conj(ZN(gt1,0))*SUM(j1,0,2,Conj(ZDXL(gt2,j1))*ZDX(gt3,j1)) + 0.4472135954999579*gN*Conj(ZN(gt1,5))*SUM(j1,0,2,Conj(ZDXL(gt2,j1))*ZDX(gt3,j1)) - Conj(ZN(gt1,4))*SUM(j1,0,2,Conj(ZDXL(gt2,j1))*ZDX(gt3,3 + j1)*Kappa(j1,j1));

   const std::complex<double> right = -(SUM(j1,0,2,Conj(Kappa(j1,j1))*ZDX(gt3,j1)*ZDXR(gt2,j1))*ZN(gt1,4)) - 0.07453559924999298*SUM(j1,0,2,ZDX(gt3,3 + j1)*ZDXR(gt2,j1))*(4.898979485566356*g1*ZN(gt1,0) - 9*gN*ZN(gt1,5));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::FDX>::type, fields::Chi, fields::SDX>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto Kappa = MODELPARAMETER(Kappa);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZDX = MODELPARAMETER(ZDX);
   const auto ZDXR = MODELPARAMETER(ZDXR);
   const auto ZDXL = MODELPARAMETER(ZDXL);

   const std::complex<double> left = -0.3651483716701107*g1*Conj(ZN(gt2,0))*SUM(j1,0,2,Conj(ZDX(gt3,3 + j1))*Conj(ZDXR(gt1,j1))) + 0.6708203932499369*gN*Conj(ZN(gt2,5))*SUM(j1,0,2,Conj(ZDX(gt3,3 + j1))*Conj(ZDXR(gt1,j1))) - Conj(ZN(gt2,4))*SUM(j1,0,2,Conj(ZDX(gt3,j1))*Conj(ZDXR(gt1,j1))*Kappa(j1,j1));

   const std::complex<double> right = -(SUM(j1,0,2,Conj(ZDX(gt3,3 + j1))*Conj(Kappa(j1,j1))*ZDXL(gt1,j1))*ZN(gt2,4)) + SUM(j1,0,2,Conj(ZDX(gt3,j1))*ZDXL(gt1,j1))*(0.3651483716701107*g1*ZN(gt2,0) + 0.4472135954999579*gN*ZN(gt2,5));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Glu, fields::FDX, typename fields::conj<fields::SDX>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt3 = indices[1];
   const auto g3 = MODELPARAMETER(g3);
   const auto ZDXL = MODELPARAMETER(ZDXL);
   const auto ZDX = MODELPARAMETER(ZDX);
   const auto ZDXR = MODELPARAMETER(ZDXR);
   const auto PhaseGlu = PHASE(PhaseGlu);

   const std::complex<double> left = -0.7071067811865475*g3*PhaseGlu*SUM(j1,0,2,Conj(ZDXL(gt2,j1))*ZDX(gt3,j1));

   const std::complex<double> right = 0.7071067811865475*g3*Conj(PhaseGlu)*SUM(j1,0,2,ZDX(gt3,3 + j1)*ZDXR(gt2,j1));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::FDX>::type, fields::Glu, fields::SDX>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const auto g3 = MODELPARAMETER(g3);
   const auto ZDX = MODELPARAMETER(ZDX);
   const auto ZDXR = MODELPARAMETER(ZDXR);
   const auto ZDXL = MODELPARAMETER(ZDXL);
   const auto PhaseGlu = PHASE(PhaseGlu);

   const std::complex<double> left = 0.7071067811865475*g3*PhaseGlu*SUM(j1,0,2,Conj(ZDX(gt3,3 + j1))*Conj(ZDXR(gt1,j1)));

   const std::complex<double> right = -0.7071067811865475*g3*Conj(PhaseGlu)*SUM(j1,0,2,Conj(ZDX(gt3,j1))*ZDXL(gt1,j1));

   return {left, right};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::SDX, typename fields::conj<fields::SDX>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Kappa = MODELPARAMETER(Kappa);
   const auto ZDX = MODELPARAMETER(ZDX);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = 0.025*(-2*SUM(j1,0,2,Conj(ZDX(gt3,j1))*ZDX(gt4,j1))*((2*Sqr(g1) + 3*Sqr(gN))*ZA(gt1,0)*ZA(gt2,0) + 2*(-Sqr(g1) + Sqr(gN))*ZA(gt1,1)*ZA(gt2,1) - 5*Sqr(gN)*ZA(gt1,2)*ZA(gt2,2)) + SUM(j1,0,2,Conj(ZDX(gt3,3 + j1))*ZDX(gt4,3 + j1))*((4*Sqr(g1) - 9*Sqr(gN))*ZA(gt1,0)*ZA(gt2,0) - 2*(2*Sqr(g1) + 3*Sqr(gN))*ZA(gt1,1)*ZA(gt2,1) + 15*Sqr(gN)*ZA(gt1,2)*ZA(gt2,2)) - 20*(Lambdax*SUM(j1,0,2,Conj(ZDX(gt3,3 + j1))*Conj(Kappa(j1,j1))*ZDX(gt4,j1))*(ZA(gt1,1)*ZA(gt2,0) + ZA(gt1,0)*ZA(gt2,1)) + Conj(Lambdax)*SUM(j1,0,2,Conj(ZDX(gt3,j1))*ZDX(gt4,3 + j1)*Kappa(j1,j1))*(ZA(gt1,1)*ZA(gt2,0) + ZA(gt1,0)*ZA(gt2,1)) + 2*(SUM(j1,0,2,AbsSqr(Kappa(j1,j1))*Conj(ZDX(gt3,j1))*ZDX(gt4,j1)) + SUM(j1,0,2,AbsSqr(Kappa(j1,j1))*Conj(ZDX(gt3,3 + j1))*ZDX(gt4,3 + j1)))*ZA(gt1,2)*ZA(gt2,2)));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::SDX, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::SDX>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto Kappa = MODELPARAMETER(Kappa);
   const auto ZDX = MODELPARAMETER(ZDX);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.025*(-40*(Conj(Lambdax)*SUM(j1,0,2,Conj(ZDX(gt2,j1))*ZDX(gt4,3 + j1)*Kappa(j1,j1))*ZP(gt1,1)*ZP(gt3,0) + Lambdax*SUM(j1,0,2,Conj(ZDX(gt2,3 + j1))*Conj(Kappa(j1,j1))*ZDX(gt4,j1))*ZP(gt1,0)*ZP(gt3,1)) - 2*SUM(j1,0,2,Conj(ZDX(gt2,j1))*ZDX(gt4,j1))*((2*Sqr(g1) + 3*Sqr(gN))*ZP(gt1,0)*ZP(gt3,0) + 2*(-Sqr(g1) + Sqr(gN))*ZP(gt1,1)*ZP(gt3,1)) + SUM(j1,0,2,Conj(ZDX(gt2,3 + j1))*ZDX(gt4,3 + j1))*((4*Sqr(g1) - 9*Sqr(gN))*ZP(gt1,0)*ZP(gt3,0) - 2*(2*Sqr(g1) + 3*Sqr(gN))*ZP(gt1,1)*ZP(gt3,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::SDX, fields::SDX, typename fields::conj<fields::SDX>::type, typename fields::conj<fields::SDX>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g3 = MODELPARAMETER(g3);
   const auto gN = MODELPARAMETER(gN);
   const auto Kappa = MODELPARAMETER(Kappa);
   const auto ZDX = MODELPARAMETER(ZDX);

   const std::complex<double> result = 0.004166666666666667*(-8*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt4,j1))*SUM(j2,0,2,Conj(ZDX(gt2,j2))*ZDX(gt3,j2)) + 20*Sqr(g3)*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt4,j1))*SUM(j2,0,2,Conj(ZDX(gt2,j2))*ZDX(gt3,j2)) - 12*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt4,j1))*SUM(j2,0,2,Conj(ZDX(gt2,j2))*ZDX(gt3,j2)) + 8*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt4,3 + j1))*SUM(j2,0,2,Conj(ZDX(gt2,j2))*ZDX(gt3,j2)) - 20*Sqr(g3)*SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt4,3 + j1))*SUM(j2,0,2,Conj(ZDX(gt2,j2))*ZDX(gt3,j2)) - 18*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt4,3 + j1))*SUM(j2,0,2,Conj(ZDX(gt2,j2))*ZDX(gt3,j2)) - 240*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt4,3 + j1)*Kappa(j1,j1))*SUM(j2,0,2,Conj(ZDX(gt2,3 + j2))*Conj(Kappa(j2,j2))*ZDX(gt3,j2)) - 60*Sqr(g3)*SUM(j1,0,2,Conj(ZDX(gt2,j1))*ZDX(gt4,j1))*(SUM(j2,0,2,Conj(ZDX(gt1,j2))*ZDX(gt3,j2)) - SUM(j2,0,2,Conj(ZDX(gt1,3 + j2))*ZDX(gt3,3 + j2))) + 60*Sqr(g3)*SUM(j1,0,2,Conj(ZDX(gt2,3 + j1))*ZDX(gt4,3 + j1))*(SUM(j2,0,2,Conj(ZDX(gt1,j2))*ZDX(gt3,j2)) - SUM(j2,0,2,Conj(ZDX(gt1,3 + j2))*ZDX(gt3,3 + j2))) + 8*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt4,j1))*SUM(j2,0,2,Conj(ZDX(gt2,3 + j2))*ZDX(gt3,3 + j2)) - 20*Sqr(g3)*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt4,j1))*SUM(j2,0,2,Conj(ZDX(gt2,3 + j2))*ZDX(gt3,3 + j2)) - 18*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt4,j1))*SUM(j2,0,2,Conj(ZDX(gt2,3 + j2))*ZDX(gt3,3 + j2)) - 8*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt4,3 + j1))*SUM(j2,0,2,Conj(ZDX(gt2,3 + j2))*ZDX(gt3,3 + j2)) + 20*Sqr(g3)*SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt4,3 + j1))*SUM(j2,0,2,Conj(ZDX(gt2,3 + j2))*ZDX(gt3,3 + j2)) - 27*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt4,3 + j1))*SUM(j2,0,2,Conj(ZDX(gt2,3 + j2))*ZDX(gt3,3 + j2)) - 8*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gt2,j1))*ZDX(gt3,j1))*SUM(j2,0,2,Conj(ZDX(gt1,j2))*ZDX(gt4,j2)) + 20*Sqr(g3)*SUM(j1,0,2,Conj(ZDX(gt2,j1))*ZDX(gt3,j1))*SUM(j2,0,2,Conj(ZDX(gt1,j2))*ZDX(gt4,j2)) - 12*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gt2,j1))*ZDX(gt3,j1))*SUM(j2,0,2,Conj(ZDX(gt1,j2))*ZDX(gt4,j2)) + 8*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gt2,3 + j1))*ZDX(gt3,3 + j1))*SUM(j2,0,2,Conj(ZDX(gt1,j2))*ZDX(gt4,j2)) - 20*Sqr(g3)*SUM(j1,0,2,Conj(ZDX(gt2,3 + j1))*ZDX(gt3,3 + j1))*SUM(j2,0,2,Conj(ZDX(gt1,j2))*ZDX(gt4,j2)) - 18*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gt2,3 + j1))*ZDX(gt3,3 + j1))*SUM(j2,0,2,Conj(ZDX(gt1,j2))*ZDX(gt4,j2)) - 60*Sqr(g3)*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt3,j1))*SUM(j2,0,2,Conj(ZDX(gt2,j2))*ZDX(gt4,j2)) + 60*Sqr(g3)*SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt3,3 + j1))*SUM(j2,0,2,Conj(ZDX(gt2,j2))*ZDX(gt4,j2)) - 240*SUM(j1,0,2,Conj(ZDX(gt2,j1))*ZDX(gt3,3 + j1)*Kappa(j1,j1))*SUM(j2,0,2,Conj(ZDX(gt1,3 + j2))*Conj(Kappa(j2,j2))*ZDX(gt4,j2)) + 8*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gt2,j1))*ZDX(gt3,j1))*SUM(j2,0,2,Conj(ZDX(gt1,3 + j2))*ZDX(gt4,3 + j2)) - 20*Sqr(g3)*SUM(j1,0,2,Conj(ZDX(gt2,j1))*ZDX(gt3,j1))*SUM(j2,0,2,Conj(ZDX(gt1,3 + j2))*ZDX(gt4,3 + j2)) - 18*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gt2,j1))*ZDX(gt3,j1))*SUM(j2,0,2,Conj(ZDX(gt1,3 + j2))*ZDX(gt4,3 + j2)) - 8*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gt2,3 + j1))*ZDX(gt3,3 + j1))*SUM(j2,0,2,Conj(ZDX(gt1,3 + j2))*ZDX(gt4,3 + j2)) + 20*Sqr(g3)*SUM(j1,0,2,Conj(ZDX(gt2,3 + j1))*ZDX(gt3,3 + j1))*SUM(j2,0,2,Conj(ZDX(gt1,3 + j2))*ZDX(gt4,3 + j2)) - 27*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gt2,3 + j1))*ZDX(gt3,3 + j1))*SUM(j2,0,2,Conj(ZDX(gt1,3 + j2))*ZDX(gt4,3 + j2)) + 60*Sqr(g3)*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt3,j1))*SUM(j2,0,2,Conj(ZDX(gt2,3 + j2))*ZDX(gt4,3 + j2)) - 60*Sqr(g3)*SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt3,3 + j1))*SUM(j2,0,2,Conj(ZDX(gt2,3 + j2))*ZDX(gt4,3 + j2)));

   return {result};
}

ScalarVertex VertexImpl<fields::SDX, fields::SHI0, typename fields::conj<fields::SDX>::type, typename fields::conj<fields::SHI0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto Kappa = MODELPARAMETER(Kappa);
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto ZDX = MODELPARAMETER(ZDX);
   const auto UHI0 = MODELPARAMETER(UHI0);

   const std::complex<double> result = 0.0125*(80*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt3,3 + j1)*Kappa(j1,j1))*SUM(j2,0,1,Conj(UHI0(gt2,2 + j2))*Conj(Lambda12(j2,j2))*UHI0(gt4,j2)) - 2*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt3,j1))*((2*Sqr(g1) + 3*Sqr(gN))*SUM(j2,0,1,Conj(UHI0(gt2,j2))*UHI0(gt4,j2)) + 2*(-Sqr(g1) + Sqr(gN))*SUM(j2,0,1,Conj(UHI0(gt2,2 + j2))*UHI0(gt4,2 + j2))) + SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt3,3 + j1))*((4*Sqr(g1) - 9*Sqr(gN))*SUM(j2,0,1,Conj(UHI0(gt2,j2))*UHI0(gt4,j2)) - 2*(2*Sqr(g1) + 3*Sqr(gN))*SUM(j2,0,1,Conj(UHI0(gt2,2 + j2))*UHI0(gt4,2 + j2))) - 4*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt4,j1))*SUM(j2,0,2,Conj(ZDX(gt1,j2))*ZDX(gt3,j2)) - 6*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt4,j1))*SUM(j2,0,2,Conj(ZDX(gt1,j2))*ZDX(gt3,j2)) + 4*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*UHI0(gt4,2 + j1))*SUM(j2,0,2,Conj(ZDX(gt1,j2))*ZDX(gt3,j2)) - 4*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*UHI0(gt4,2 + j1))*SUM(j2,0,2,Conj(ZDX(gt1,j2))*ZDX(gt3,j2)) + 80*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt4,2 + j1)*Lambda12(j1,j1))*SUM(j2,0,2,Conj(ZDX(gt1,3 + j2))*Conj(Kappa(j2,j2))*ZDX(gt3,j2)) + 4*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt4,j1))*SUM(j2,0,2,Conj(ZDX(gt1,3 + j2))*ZDX(gt3,3 + j2)) - 9*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt4,j1))*SUM(j2,0,2,Conj(ZDX(gt1,3 + j2))*ZDX(gt3,3 + j2)) - 4*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*UHI0(gt4,2 + j1))*SUM(j2,0,2,Conj(ZDX(gt1,3 + j2))*ZDX(gt3,3 + j2)) - 6*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*UHI0(gt4,2 + j1))*SUM(j2,0,2,Conj(ZDX(gt1,3 + j2))*ZDX(gt3,3 + j2)));

   return {result};
}

ScalarVertex VertexImpl<fields::SDX, fields::SHIp, typename fields::conj<fields::SDX>::type, typename fields::conj<fields::SHIp>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto Kappa = MODELPARAMETER(Kappa);
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto ZDX = MODELPARAMETER(ZDX);
   const auto UHIp = MODELPARAMETER(UHIp);

   const std::complex<double> result = 0.0125*(-80*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt3,3 + j1)*Kappa(j1,j1))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*Conj(Lambda12(j2,j2))*UHIp(gt4,j2)) - 2*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt3,j1))*((2*Sqr(g1) + 3*Sqr(gN))*SUM(j2,0,1,Conj(UHIp(gt2,j2))*UHIp(gt4,j2)) + 2*(-Sqr(g1) + Sqr(gN))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*UHIp(gt4,2 + j2))) + SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt3,3 + j1))*((4*Sqr(g1) - 9*Sqr(gN))*SUM(j2,0,1,Conj(UHIp(gt2,j2))*UHIp(gt4,j2)) - 2*(2*Sqr(g1) + 3*Sqr(gN))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*UHIp(gt4,2 + j2))) - 4*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt4,j1))*SUM(j2,0,2,Conj(ZDX(gt1,j2))*ZDX(gt3,j2)) - 6*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt4,j1))*SUM(j2,0,2,Conj(ZDX(gt1,j2))*ZDX(gt3,j2)) + 4*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*UHIp(gt4,2 + j1))*SUM(j2,0,2,Conj(ZDX(gt1,j2))*ZDX(gt3,j2)) - 4*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*UHIp(gt4,2 + j1))*SUM(j2,0,2,Conj(ZDX(gt1,j2))*ZDX(gt3,j2)) - 80*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt4,2 + j1)*Lambda12(j1,j1))*SUM(j2,0,2,Conj(ZDX(gt1,3 + j2))*Conj(Kappa(j2,j2))*ZDX(gt3,j2)) + 4*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt4,j1))*SUM(j2,0,2,Conj(ZDX(gt1,3 + j2))*ZDX(gt3,3 + j2)) - 9*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt4,j1))*SUM(j2,0,2,Conj(ZDX(gt1,3 + j2))*ZDX(gt3,3 + j2)) - 4*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*UHIp(gt4,2 + j1))*SUM(j2,0,2,Conj(ZDX(gt1,3 + j2))*ZDX(gt3,3 + j2)) - 6*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*UHIp(gt4,2 + j1))*SUM(j2,0,2,Conj(ZDX(gt1,3 + j2))*ZDX(gt3,3 + j2)));

   return {result};
}

ScalarVertex VertexImpl<fields::SDX, fields::SHp0, typename fields::conj<fields::SDX>::type, typename fields::conj<fields::SHp0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZDX = MODELPARAMETER(ZDX);
   const auto UHp0 = MODELPARAMETER(UHp0);

   const std::complex<double> result = -0.05*(2*(Sqr(g1) - Sqr(gN))*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt3,j1)) - (2*Sqr(g1) + 3*Sqr(gN))*SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt3,3 + j1)))*(Conj(UHp0(gt2,0))*UHp0(gt4,0) - Conj(UHp0(gt2,1))*UHp0(gt4,1));

   return {result};
}

ScalarVertex VertexImpl<fields::SDX, fields::SHpp, typename fields::conj<fields::SDX>::type, typename fields::conj<fields::SHpp>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZDX = MODELPARAMETER(ZDX);
   const auto UHpp = MODELPARAMETER(UHpp);

   const std::complex<double> result = -0.05*(2*(Sqr(g1) - Sqr(gN))*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt3,j1)) - (2*Sqr(g1) + 3*Sqr(gN))*SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt3,3 + j1)))*(Conj(UHpp(gt2,0))*UHpp(gt4,0) - Conj(UHpp(gt2,1))*UHpp(gt4,1));

   return {result};
}

ScalarVertex VertexImpl<fields::SDX, fields::SSI0, typename fields::conj<fields::SDX>::type, typename fields::conj<fields::SSI0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto gN = MODELPARAMETER(gN);
   const auto ZDX = MODELPARAMETER(ZDX);

   const std::complex<double> result = 0.125*KroneckerDelta(gt2,gt4)*Sqr(gN)*(2*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt3,j1)) + 3*SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt3,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SDX, typename fields::conj<fields::SDX>::type, fields::VG, fields::VG>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g3 = MODELPARAMETER(g3);

   const std::complex<double> result = 0.25*KroneckerDelta(gt1,gt2)*Sqr(g3);

   return {result};
}

InverseMetricVertex VertexImpl<fields::SDX, typename fields::conj<fields::SDX>::type, fields::VP, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.13333333333333333*KroneckerDelta(gt1,gt2)*Sqr(g1)*Sqr(Cos(ThetaW));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SDX, typename fields::conj<fields::SDX>::type, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZDX = MODELPARAMETER(ZDX);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.016666666666666666*(4*(-4.898979485566356*g1*gN*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) + 2*Sqr(g1)*Sqr(Cos(ThetaWp))*Sqr(Sin(ThetaW)) + 3*Sqr(gN)*Sqr(Sin(ThetaWp)))*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt2,j1)) + (8*Sqr(g1)*Sqr(Cos(ThetaWp))*Sqr(Sin(ThetaW)) + 3*gN*(4.898979485566356*g1*Sin(ThetaW)*Sin(2*ThetaWp) + 9*gN*Sqr(Sin(ThetaWp))))*SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SDX, typename fields::conj<fields::SDX>::type, fields::VZp, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZDX = MODELPARAMETER(ZDX);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.016666666666666666*(4*(3*Sqr(gN)*Sqr(Cos(ThetaWp)) + g1*Sin(ThetaW)*(2.449489742783178*gN*Sin(2*ThetaWp) + 2*g1*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt2,j1)) + (27*Sqr(gN)*Sqr(Cos(ThetaWp)) + 2*g1*Sin(ThetaW)*(-7.348469228349534*gN*Sin(2*ThetaWp) + 4*g1*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SDX, typename fields::conj<fields::SDX>::type, fields::VZ, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZDX = MODELPARAMETER(ZDX);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.016666666666666666*((6*Sin(2*ThetaWp)*Sqr(gN) - 9.797958971132712*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) - 4*Sin(2*ThetaWp)*Sqr(g1)*Sqr(Sin(ThetaW)) + 9.797958971132712*g1*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp)))*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt2,j1)) + (14.696938456699067*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + Cos(ThetaWp)*Sin(ThetaWp)*(27*Sqr(gN) - 8*Sqr(g1)*Sqr(Sin(ThetaW))) - 14.696938456699067*g1*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp)))*SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::hh, fields::VZ, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05*(-((9.486832980505138*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) - 9*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN) + 5*Sin(2*ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 7.348469228349534*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + g1*(3.872983346207417*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) + 3*g1*Sin(2*ThetaWp)*Sqr(Sin(ThetaW)) - 7.348469228349534*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*ZH(gt1,0)*ZH(gt2,0)) + (6.324555320336759*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) + 4*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN) - 5*Sin(2*ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 4.898979485566356*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) - g1*(3.872983346207417*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) + 3*g1*Sin(2*ThetaWp)*Sqr(Sin(ThetaW)) + 4.898979485566356*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*ZH(gt1,1)*ZH(gt2,1) + 25*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN)*ZH(gt1,2)*ZH(gt2,2));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Ah, fields::Ah>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto gN = MODELPARAMETER(gN);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = 0.025*(-(ZA(gt1,0)*(-2*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZA(gt2,1)*(ZA(gt3,1)*ZA(gt4,0) + ZA(gt3,0)*ZA(gt4,1)) - 5*(-8*AbsSqr(Lambdax) + 3*Sqr(gN))*ZA(gt2,2)*(ZA(gt3,2)*ZA(gt4,0) + ZA(gt3,0)*ZA(gt4,2)) + ZA(gt2,0)*(3*(6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*ZA(gt3,0)*ZA(gt4,0) + 2*(20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2) + 3*Sqr(gN))*ZA(gt3,1)*ZA(gt4,1) + 5*(8*AbsSqr(Lambdax) - 3*Sqr(gN))*ZA(gt3,2)*ZA(gt4,2)))) + 5*ZA(gt1,2)*((-8*AbsSqr(Lambdax) + 3*Sqr(gN))*ZA(gt2,0)*(ZA(gt3,2)*ZA(gt4,0) + ZA(gt3,0)*ZA(gt4,2)) + 2*(-4*AbsSqr(Lambdax) + Sqr(gN))*ZA(gt2,1)*(ZA(gt3,2)*ZA(gt4,1) + ZA(gt3,1)*ZA(gt4,2)) + ZA(gt2,2)*((-8*AbsSqr(Lambdax) + 3*Sqr(gN))*ZA(gt3,0)*ZA(gt4,0) + 2*(-4*AbsSqr(Lambdax) + Sqr(gN))*ZA(gt3,1)*ZA(gt4,1) - 15*Sqr(gN)*ZA(gt3,2)*ZA(gt4,2))) + 2*ZA(gt1,1)*((-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZA(gt2,0)*(ZA(gt3,1)*ZA(gt4,0) + ZA(gt3,0)*ZA(gt4,1)) + 5*(-4*AbsSqr(Lambdax) + Sqr(gN))*ZA(gt2,2)*(ZA(gt3,2)*ZA(gt4,1) + ZA(gt3,1)*ZA(gt4,2)) + ZA(gt2,1)*((-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZA(gt3,0)*ZA(gt4,0) - 3*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZA(gt3,1)*ZA(gt4,1) + 5*(-4*AbsSqr(Lambdax) + Sqr(gN))*ZA(gt3,2)*ZA(gt4,2))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Hpm, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto gN = MODELPARAMETER(gN);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.025*(5*ZA(gt1,2)*ZA(gt2,2)*((-8*AbsSqr(Lambdax) + 3*Sqr(gN))*ZP(gt3,0)*ZP(gt4,0) + 2*(-4*AbsSqr(Lambdax) + Sqr(gN))*ZP(gt3,1)*ZP(gt4,1)) + 2*ZA(gt1,1)*(5*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZA(gt2,0)*(ZP(gt3,1)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,1)) + ZA(gt2,1)*((3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZP(gt3,0)*ZP(gt4,0) - (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZP(gt3,1)*ZP(gt4,1))) - ZA(gt1,0)*(-10*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZA(gt2,1)*(ZP(gt3,1)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,1)) + ZA(gt2,0)*((6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*ZP(gt3,0)*ZP(gt4,0) + 2*(-3*Sqr(g1) + 5*Sqr(g2) + 3*Sqr(gN))*ZP(gt3,1)*ZP(gt4,1))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::SHI0, typename fields::conj<fields::SHI0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = 0.025*(2*SUM(j1,0,1,Conj(UHI0(gt3,2 + j1))*UHI0(gt4,2 + j1))*((3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZA(gt1,0)*ZA(gt2,0) - (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZA(gt1,1)*ZA(gt2,1) + 5*Sqr(gN)*ZA(gt1,2)*ZA(gt2,2)) + SUM(j1,0,1,Conj(UHI0(gt3,j1))*UHI0(gt4,j1))*(-((6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*ZA(gt1,0)*ZA(gt2,0)) + 2*(3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZA(gt1,1)*ZA(gt2,1) + 15*Sqr(gN)*ZA(gt1,2)*ZA(gt2,2)) + 20*(Lambdax*SUM(j1,0,1,Conj(UHI0(gt3,2 + j1))*Conj(Lambda12(j1,j1))*UHI0(gt4,j1))*(ZA(gt1,1)*ZA(gt2,0) + ZA(gt1,0)*ZA(gt2,1)) + Conj(Lambdax)*SUM(j1,0,1,Conj(UHI0(gt3,j1))*UHI0(gt4,2 + j1)*Lambda12(j1,j1))*(ZA(gt1,1)*ZA(gt2,0) + ZA(gt1,0)*ZA(gt2,1)) - 2*(SUM(j1,0,1,AbsSqr(Lambda12(j1,j1))*Conj(UHI0(gt3,j1))*UHI0(gt4,j1)) + SUM(j1,0,1,AbsSqr(Lambda12(j1,j1))*Conj(UHI0(gt3,2 + j1))*UHI0(gt4,2 + j1)))*ZA(gt1,2)*ZA(gt2,2)));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::SHIp, typename fields::conj<fields::SHIp>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = 0.025*(2*SUM(j1,0,1,Conj(UHIp(gt3,2 + j1))*UHIp(gt4,2 + j1))*((3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZA(gt1,0)*ZA(gt2,0) + (-3*Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN))*ZA(gt1,1)*ZA(gt2,1) + 5*Sqr(gN)*ZA(gt1,2)*ZA(gt2,2)) + SUM(j1,0,1,Conj(UHIp(gt3,j1))*UHIp(gt4,j1))*((-6*Sqr(g1) + 10*Sqr(g2) - 9*Sqr(gN))*ZA(gt1,0)*ZA(gt2,0) + 2*(3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZA(gt1,1)*ZA(gt2,1) + 15*Sqr(gN)*ZA(gt1,2)*ZA(gt2,2)) - 20*(Lambdax*SUM(j1,0,1,Conj(UHIp(gt3,2 + j1))*Conj(Lambda12(j1,j1))*UHIp(gt4,j1))*(ZA(gt1,1)*ZA(gt2,0) + ZA(gt1,0)*ZA(gt2,1)) + Conj(Lambdax)*SUM(j1,0,1,Conj(UHIp(gt3,j1))*UHIp(gt4,2 + j1)*Lambda12(j1,j1))*(ZA(gt1,1)*ZA(gt2,0) + ZA(gt1,0)*ZA(gt2,1)) + 2*(SUM(j1,0,1,AbsSqr(Lambda12(j1,j1))*Conj(UHIp(gt3,j1))*UHIp(gt4,j1)) + SUM(j1,0,1,AbsSqr(Lambda12(j1,j1))*Conj(UHIp(gt3,2 + j1))*UHIp(gt4,2 + j1)))*ZA(gt1,2)*ZA(gt2,2)));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::SHp0, typename fields::conj<fields::SHp0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = -0.05*(Conj(UHp0(gt3,0))*UHp0(gt4,0) - Conj(UHp0(gt3,1))*UHp0(gt4,1))*((3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZA(gt1,0)*ZA(gt2,0) - (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZA(gt1,1)*ZA(gt2,1) + 5*Sqr(gN)*ZA(gt1,2)*ZA(gt2,2));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::SHpp, typename fields::conj<fields::SHpp>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = -0.05*(Conj(UHpp(gt3,0))*UHpp(gt4,0) - Conj(UHpp(gt3,1))*UHpp(gt4,1))*((3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZA(gt1,0)*ZA(gt2,0) + (-3*Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN))*ZA(gt1,1)*ZA(gt2,1) + 5*Sqr(gN)*ZA(gt1,2)*ZA(gt2,2));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::SSI0, typename fields::conj<fields::SSI0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto gN = MODELPARAMETER(gN);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = 0.125*KroneckerDelta(gt3,gt4)*Sqr(gN)*(3*ZA(gt1,0)*ZA(gt2,0) + 2*ZA(gt1,1)*ZA(gt2,1) - 5*ZA(gt1,2)*ZA(gt2,2));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Ah, fields::Ah, typename fields::conj<fields::VWm>::type, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = 0.5*Sqr(g2)*(ZA(gt1,0)*ZA(gt2,0) + ZA(gt1,1)*ZA(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Ah, fields::Ah, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05*((-14.696938456699067*g1*gN*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) + 10*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Cos(ThetaWp)) + Cos(ThetaW)*(-18.973665961010276*g2*gN*Cos(ThetaWp)*Sin(ThetaWp) + 15.491933384829668*g1*g2*Sin(ThetaW)*Sqr(Cos(ThetaWp))) + 6*Sqr(g1)*Sqr(Cos(ThetaWp))*Sqr(Sin(ThetaW)) + 9*Sqr(gN)*Sqr(Sin(ThetaWp)))*ZA(gt1,0)*ZA(gt2,0) + 2*(3.1622776601683795*g2*gN*Cos(ThetaW)*Sin(2*ThetaWp) + g1*Sin(ThetaW)*(7.745966692414834*g2*Cos(ThetaW) + 3*g1*Sin(ThetaW))*Sqr(Cos(ThetaWp)) + 5*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Cos(ThetaWp)) + gN*(2.449489742783178*g1*Sin(ThetaW)*Sin(2*ThetaWp) + 2*gN*Sqr(Sin(ThetaWp))))*ZA(gt1,1)*ZA(gt2,1) + 25*Sqr(gN)*Sqr(Sin(ThetaWp))*ZA(gt1,2)*ZA(gt2,2));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Ah, fields::Ah, fields::VZp, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*((9*Sqr(gN)*Sqr(Cos(ThetaWp)) + (g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(9.486832980505138*gN*Sin(2*ThetaWp) + 10*g2*Cos(ThetaW)*Sqr(Sin(ThetaWp)) + 7.745966692414834*g1*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*ZA(gt1,0)*ZA(gt2,0) + 2*(-3.1622776601683795*gN*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*Sin(2*ThetaWp) + 2*Sqr(gN)*Sqr(Cos(ThetaWp)) + 5*Sqr(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*Sqr(Sin(ThetaWp)))*ZA(gt1,1)*ZA(gt2,1) + 25*Sqr(gN)*Sqr(Cos(ThetaWp))*ZA(gt1,2)*ZA(gt2,2));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Ah, fields::Ah, fields::VZ, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05*(-((9.486832980505138*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) - 9*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN) + 5*Sin(2*ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 7.348469228349534*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + g1*(3.872983346207417*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) + 3*g1*Sin(2*ThetaWp)*Sqr(Sin(ThetaW)) - 7.348469228349534*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*ZA(gt1,0)*ZA(gt2,0)) + (6.324555320336759*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) + 4*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN) - 5*Sin(2*ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 4.898979485566356*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) - g1*(3.872983346207417*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) + 3*g1*Sin(2*ThetaWp)*Sqr(Sin(ThetaW)) + 4.898979485566356*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*ZA(gt1,1)*ZA(gt2,1) + 25*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN)*ZA(gt1,2)*ZA(gt2,2));

   return {result};
}

ChiralVertex VertexImpl<fields::Chi, fields::Cha, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto UM = MODELPARAMETER(UM);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZP = MODELPARAMETER(ZP);
   const auto UP = MODELPARAMETER(UP);

   const std::complex<double> left = -(g2*Conj(UM(gt2,0))*Conj(ZN(gt1,2))*ZP(gt3,0)) + 0.1*Conj(UM(gt2,1))*(5.477225575051661*g1*Conj(ZN(gt1,0))*ZP(gt3,0) + 7.0710678118654755*g2*Conj(ZN(gt1,1))*ZP(gt3,0) + 6.708203932499369*gN*Conj(ZN(gt1,5))*ZP(gt3,0) - 10*Conj(ZN(gt1,4))*Lambdax*ZP(gt3,1));

   const std::complex<double> right = -(Conj(Lambdax)*UP(gt2,1)*ZN(gt1,4)*ZP(gt3,0)) - 0.1*(10*g2*UP(gt2,0)*ZN(gt1,3) + UP(gt2,1)*(5.477225575051661*g1*ZN(gt1,0) + 7.0710678118654755*g2*ZN(gt1,1) - 4.47213595499958*gN*ZN(gt1,5)))*ZP(gt3,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha>::type, fields::Chi, fields::Hpm>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto UP = MODELPARAMETER(UP);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZP = MODELPARAMETER(ZP);
   const auto UM = MODELPARAMETER(UM);

   const std::complex<double> left = -(g2*Conj(UP(gt1,0))*Conj(ZN(gt2,3))*ZP(gt3,1)) - 0.1*Conj(UP(gt1,1))*(10*Conj(ZN(gt2,4))*Lambdax*ZP(gt3,0) + (5.477225575051661*g1*Conj(ZN(gt2,0)) + 7.0710678118654755*g2*Conj(ZN(gt2,1)) - 4.47213595499958*gN*Conj(ZN(gt2,5)))*ZP(gt3,1));

   const std::complex<double> right = -(g2*UM(gt1,0)*ZN(gt2,2)*ZP(gt3,0)) + 0.1*UM(gt1,1)*(5.477225575051661*g1*ZN(gt2,0)*ZP(gt3,0) + 7.0710678118654755*g2*ZN(gt2,1)*ZP(gt3,0) + 6.708203932499369*gN*ZN(gt2,5)*ZP(gt3,0) - 10*Conj(Lambdax)*ZN(gt2,4)*ZP(gt3,1));

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

   const std::complex<double> left = SUM(j1,0,2,Conj(ZDL(gt2,j1))*Conj(ZUR(gt1,j1))*Yu(j1,j1))*ZP(gt3,1);

   const std::complex<double> right = SUM(j1,0,2,Conj(Yd(j1,j1))*ZDR(gt2,j1)*ZUL(gt1,j1))*ZP(gt3,0);

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
   const auto ZDR = MODELPARAMETER(ZDR);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto ZUR = MODELPARAMETER(ZUR);

   const std::complex<double> left = SUM(j1,0,2,Conj(ZDR(gt1,j1))*Conj(ZUL(gt2,j1))*Yd(j1,j1))*ZP(gt3,0);

   const std::complex<double> right = SUM(j1,0,2,Conj(Yu(j1,j1))*ZDL(gt1,j1)*ZUR(gt2,j1))*ZP(gt3,1);

   return {left, right};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gZ, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.025*g2*(vd*(10*g2*Cos(ThetaW)*Cos(ThetaWp) - 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp))*ZP(gt3,0) + 2*vu*(-5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*ZP(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gWmC, fields::Hpm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -0.025*g2*(vd*(10*g2*Cos(ThetaW)*Cos(ThetaWp) + 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp))*ZP(gt3,0) - 2*vu*(5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*ZP(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gZp, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.025*(g2*vd*(9.486832980505138*gN*Cos(ThetaWp) + 2*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*ZP(gt3,0) + 2*g2*vu*(3.1622776601683795*gN*Cos(ThetaWp) + (5*g2*Cos(ThetaW) - 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*ZP(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gZp>::type, fields::gWmC, fields::Hpm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.025*g2*(vd*(9.486832980505138*gN*Cos(ThetaWp) + 2*(5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*ZP(gt3,0) + 2*vu*(3.1622776601683795*gN*Cos(ThetaWp) - (5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*ZP(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gWm, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -0.025*g2*(vd*(10*g2*Cos(ThetaW)*Cos(ThetaWp) + 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp))*ZP(gt3,0) - 2*vu*(5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*ZP(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gWm>::type, fields::gZ, fields::Hpm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.025*g2*(vd*(10*g2*Cos(ThetaW)*Cos(ThetaWp) - 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp))*ZP(gt3,0) + 2*vu*(-5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*ZP(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gWm>::type, fields::gZp, fields::Hpm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.025*(g2*vd*(9.486832980505138*gN*Cos(ThetaWp) + 2*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*ZP(gt3,0) + 2*g2*vu*(3.1622776601683795*gN*Cos(ThetaWp) + (5*g2*Cos(ThetaW) - 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*ZP(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gZp>::type, fields::gWm, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.025*g2*(vd*(9.486832980505138*gN*Cos(ThetaWp) + 2*(5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*ZP(gt3,0) + 2*vu*(3.1622776601683795*gN*Cos(ThetaWp) - (5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*ZP(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::SHIp, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::SHI0>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -0.35355339059327373*Sqr(g2)*SUM(j1,0,3,Conj(UHIp(gt1,j1))*UHI0(gt3,j1))*(vd*ZP(gt2,0) + vu*ZP(gt2,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::SHI0, typename fields::conj<fields::SHIp>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -0.35355339059327373*Sqr(g2)*SUM(j1,0,3,Conj(UHI0(gt2,j1))*UHIp(gt3,j1))*(vd*ZP(gt1,0) + vu*ZP(gt1,1));

   return {result};
}

ScalarVertex VertexImpl<fields::SHpp, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::SHp0>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -0.35355339059327373*Sqr(g2)*(Conj(UHpp(gt1,0))*UHp0(gt3,0) + Conj(UHpp(gt1,1))*UHp0(gt3,1))*(vd*ZP(gt2,0) + vu*ZP(gt2,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::SHp0, typename fields::conj<fields::SHpp>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -0.35355339059327373*Sqr(g2)*(Conj(UHp0(gt2,0))*UHpp(gt3,0) + Conj(UHp0(gt2,1))*UHpp(gt3,1))*(vd*ZP(gt1,0) + vu*ZP(gt1,1));

   return {result};
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

   const std::complex<double> right = IF(gt1 < 3,Conj(Ye(gt1,gt1))*ZER(gt2,gt1)*ZP(gt3,0),0);

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

   const std::complex<double> left = IF(gt2 < 3,Conj(ZER(gt1,gt2))*Ye(gt2,gt2)*ZP(gt3,0),0);

   const std::complex<double> right = 0;

   return {left, right};
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
   const auto gN = MODELPARAMETER(gN);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.05*(-((6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*ZP(gt1,0)*ZP(gt2,0)*ZP(gt3,0)*ZP(gt4,0)) - 2*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZP(gt1,1)*ZP(gt2,1)*ZP(gt3,1)*ZP(gt4,1) + (-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZP(gt1,1)*ZP(gt2,0)*(ZP(gt3,1)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,1)) + (-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZP(gt1,0)*ZP(gt2,1)*(ZP(gt3,1)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::SHI0, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::SHI0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = Conj(Lambdax)*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt4,2 + j1)*Lambda12(j1,j1))*ZP(gt1,1)*ZP(gt3,0) + Lambdax*SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*Conj(Lambda12(j1,j1))*UHI0(gt4,j1))*ZP(gt1,0)*ZP(gt3,1) + 0.05*SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*UHI0(gt4,2 + j1))*((3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZP(gt1,0)*ZP(gt3,0) + (-3*Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN))*ZP(gt1,1)*ZP(gt3,1)) - 0.025*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt4,j1))*((6*Sqr(g1) - 10*Sqr(g2) + 9*Sqr(gN))*ZP(gt1,0)*ZP(gt3,0) + 2*(-3*Sqr(g1) + 5*Sqr(g2) + 3*Sqr(gN))*ZP(gt1,1)*ZP(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::SHIp, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::SHIp>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.025*(-40*(Conj(Lambdax)*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt4,2 + j1)*Lambda12(j1,j1))*ZP(gt1,1)*ZP(gt3,0) + Lambdax*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*Conj(Lambda12(j1,j1))*UHIp(gt4,j1))*ZP(gt1,0)*ZP(gt3,1)) - SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt4,j1))*((6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*ZP(gt1,0)*ZP(gt3,0) - 2*(3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZP(gt1,1)*ZP(gt3,1)) + 2*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*UHIp(gt4,2 + j1))*((3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZP(gt1,0)*ZP(gt3,0) - (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZP(gt1,1)*ZP(gt3,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::SHp0, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::SHp0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -0.05*(Conj(UHp0(gt2,0))*UHp0(gt4,0) - Conj(UHp0(gt2,1))*UHp0(gt4,1))*((3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZP(gt1,0)*ZP(gt3,0) + (-3*Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN))*ZP(gt1,1)*ZP(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::SHpp, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::SHpp>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -0.05*(Conj(UHpp(gt2,0))*UHpp(gt4,0) - Conj(UHpp(gt2,1))*UHpp(gt4,1))*((3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZP(gt1,0)*ZP(gt3,0) - (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZP(gt1,1)*ZP(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::SSI0, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::SSI0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto gN = MODELPARAMETER(gN);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.125*KroneckerDelta(gt2,gt4)*Sqr(gN)*(3*ZP(gt1,0)*ZP(gt3,0) + 2*ZP(gt1,1)*ZP(gt3,1));

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

   const std::complex<double> result = 0.5*Sqr(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*(ZP(gt1,0)*ZP(gt2,0) + ZP(gt1,1)*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::VWm>::type, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.5*Sqr(g2)*(ZP(gt1,0)*ZP(gt2,0) + ZP(gt1,1)*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05*((-14.696938456699067*g1*gN*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) + 2*g2*Cos(ThetaW)*Cos(ThetaWp)*(-7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp)) + 10*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Cos(ThetaWp)) + 6*Sqr(g1)*Sqr(Cos(ThetaWp))*Sqr(Sin(ThetaW)) + 9*Sqr(gN)*Sqr(Sin(ThetaWp)))*ZP(gt1,0)*ZP(gt2,0) + 2*(5*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Cos(ThetaWp)) - 2*Cos(ThetaW)*(3.1622776601683795*g2*gN*Cos(ThetaWp)*Sin(ThetaWp) + 3.872983346207417*g1*g2*Sin(ThetaW)*Sqr(Cos(ThetaWp))) + 3*Sqr(g1)*Sqr(Cos(ThetaWp))*Sqr(Sin(ThetaW)) + gN*(2.449489742783178*g1*Sin(ThetaW)*Sin(2*ThetaWp) + 2*gN*Sqr(Sin(ThetaWp))))*ZP(gt1,1)*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, fields::VZp, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*((9*Sqr(gN)*Sqr(Cos(ThetaWp)) + (g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*(-9.486832980505138*gN*Sin(2*ThetaWp) + 10*g2*Cos(ThetaW)*Sqr(Sin(ThetaWp)) - 7.745966692414834*g1*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*ZP(gt1,0)*ZP(gt2,0) + 2*(2*Sqr(gN)*Sqr(Cos(ThetaWp)) + 5*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Sin(ThetaWp)) + g1*Sin(ThetaW)*(-2.449489742783178*gN*Sin(2*ThetaWp) + 3*g1*Sin(ThetaW)*Sqr(Sin(ThetaWp))) + Cos(ThetaW)*(3.1622776601683795*g2*gN*Sin(2*ThetaWp) - 7.745966692414834*g1*g2*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*ZP(gt1,1)*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, fields::VZ, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05*((9.486832980505138*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) + 9*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN) - 5*Sin(2*ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) - 7.348469228349534*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + g1*(3.872983346207417*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) - 3*g1*Sin(2*ThetaWp)*Sqr(Sin(ThetaW)) + 7.348469228349534*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*ZP(gt1,0)*ZP(gt2,0) + (-6.324555320336759*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) + 4*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN) - 5*Sin(2*ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 4.898979485566356*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + g1*(3.872983346207417*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) - 3*g1*Sin(2*ThetaWp)*Sqr(Sin(ThetaW)) - 4.898979485566356*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*ZP(gt1,1)*ZP(gt2,1));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Hpm, fields::SHI0, typename fields::conj<fields::SHIp>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -0.35355339059327373*Sqr(g2)*SUM(j1,0,3,Conj(UHI0(gt3,j1))*UHIp(gt4,j1))*(ZH(gt1,0)*ZP(gt2,0) + ZH(gt1,1)*ZP(gt2,1));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Hpm, fields::SHp0, typename fields::conj<fields::SHpp>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -0.35355339059327373*Sqr(g2)*(Conj(UHp0(gt3,0))*UHpp(gt4,0) + Conj(UHp0(gt3,1))*UHpp(gt4,1))*(ZH(gt1,0)*ZP(gt2,0) + ZH(gt1,1)*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.3872983346207417*g1*g2*Cos(ThetaW)*(ZH(gt1,0)*ZP(gt2,0) - ZH(gt1,1)*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*g2*((7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp))*ZH(gt1,0)*ZP(gt2,0) - 2*(3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*ZH(gt1,1)*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.05*g2*((9.486832980505138*gN*Cos(ThetaWp) + 7.745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZH(gt1,0)*ZP(gt2,0) + 2*(3.1622776601683795*gN*Cos(ThetaWp) - 3.872983346207417*g1*Sin(ThetaW)*Sin(ThetaWp))*ZH(gt1,1)*ZP(gt2,1));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::SHIp, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::SHI0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -0.35355339059327373*Sqr(g2)*SUM(j1,0,3,Conj(UHIp(gt2,j1))*UHI0(gt4,j1))*(ZH(gt1,0)*ZP(gt3,0) + ZH(gt1,1)*ZP(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::SHpp, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::SHp0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -0.35355339059327373*Sqr(g2)*(Conj(UHpp(gt2,0))*UHp0(gt4,0) + Conj(UHpp(gt2,1))*UHp0(gt4,1))*(ZH(gt1,0)*ZP(gt3,0) + ZH(gt1,1)*ZP(gt3,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, typename fields::conj<fields::Hpm>::type, fields::VP, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.3872983346207417*g1*g2*Cos(ThetaW)*(ZH(gt1,0)*ZP(gt2,0) - ZH(gt1,1)*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, typename fields::conj<fields::Hpm>::type, fields::VWm, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*g2*((7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp))*ZH(gt1,0)*ZP(gt2,0) - 2*(3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*ZH(gt1,1)*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, typename fields::conj<fields::Hpm>::type, fields::VWm, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.05*g2*((9.486832980505138*gN*Cos(ThetaWp) + 7.745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZH(gt1,0)*ZP(gt2,0) + 2*(3.1622776601683795*gN*Cos(ThetaWp) - 3.872983346207417*g1*Sin(ThetaW)*Sin(ThetaWp))*ZH(gt1,1)*ZP(gt2,1));

   return {result};
}

ChiralVertex VertexImpl<fields::Chi, fields::Cha, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto UM = MODELPARAMETER(UM);
   const auto ZN = MODELPARAMETER(ZN);
   const auto UP = MODELPARAMETER(UP);

   const std::complex<double> left = -0.5*g2*(2*Conj(UM(gt2,0))*ZN(gt1,1) + 1.4142135623730951*Conj(UM(gt2,1))*ZN(gt1,2));

   const std::complex<double> right = -(g2*Conj(ZN(gt1,1))*UP(gt2,0)) + 0.7071067811865475*g2*Conj(ZN(gt1,3))*UP(gt2,1);

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
   const auto g1 = MODELPARAMETER(g1);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.25*g2*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*(vd*ZP(gt3,0) - vu*ZP(gt3,1));

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
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -(g2*Cos(ThetaW)*Cos(ThetaWp));

   return {result, 1};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gZp, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = g2*Cos(ThetaW)*Sin(ThetaWp);

   return {result, 1};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gWm, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = g2*Cos(ThetaW)*Cos(ThetaWp);

   return {result, 1};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gZp>::type, fields::gWm, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -(g2*Cos(ThetaW)*Sin(ThetaWp));

   return {result, 1};
}

MomentumDifferenceVertex VertexImpl<fields::SHIp, typename fields::conj<fields::SHI0>::type, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto UHI0 = MODELPARAMETER(UHI0);

   const std::complex<double> result = 0.7071067811865475*g2*(-SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHI0(gt2,j1)) + SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHI0(gt2,2 + j1)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::SHpp, typename fields::conj<fields::SHp0>::type, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto UHp0 = MODELPARAMETER(UHp0);

   const std::complex<double> result = 0.7071067811865475*(-(g2*Conj(UHpp(gt1,0))*UHp0(gt2,0)) + g2*Conj(UHpp(gt1,1))*UHp0(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha>::type, fields::Chi, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZN = MODELPARAMETER(ZN);
   const auto UM = MODELPARAMETER(UM);
   const auto UP = MODELPARAMETER(UP);

   const std::complex<double> left = -0.5*g2*(2*Conj(ZN(gt2,1))*UM(gt1,0) + 1.4142135623730951*Conj(ZN(gt2,2))*UM(gt1,1));

   const std::complex<double> right = -(g2*Conj(UP(gt1,0))*ZN(gt2,1)) + 0.7071067811865475*g2*Conj(UP(gt1,1))*ZN(gt2,3);

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
   const auto g1 = MODELPARAMETER(g1);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.25*g2*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*(vd*ZP(gt3,0) - vu*ZP(gt3,1));

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
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = g2*Cos(ThetaW)*Cos(ThetaWp);

   return {result, 1};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWm>::type, fields::gZp, fields::VWm>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -(g2*Cos(ThetaW)*Sin(ThetaWp));

   return {result, 1};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gWmC, fields::VWm>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -(g2*Cos(ThetaW)*Cos(ThetaWp));

   return {result, 1};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gZp>::type, fields::gWmC, fields::VWm>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = g2*Cos(ThetaW)*Sin(ThetaWp);

   return {result, 1};
}

MomentumDifferenceVertex VertexImpl<fields::SHI0, typename fields::conj<fields::SHIp>::type, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto UHIp = MODELPARAMETER(UHIp);

   const std::complex<double> result = 0.7071067811865475*g2*(-SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHIp(gt2,j1)) + SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHIp(gt2,2 + j1)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::SHp0, typename fields::conj<fields::SHpp>::type, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto UHpp = MODELPARAMETER(UHpp);

   const std::complex<double> result = 0.7071067811865475*(-(g2*Conj(UHp0(gt1,0))*UHpp(gt2,0)) + g2*Conj(UHp0(gt1,1))*UHpp(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

ChiralVertex VertexImpl<fields::Chi, fields::ChiI, typename fields::conj<fields::SHI0>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZNI = MODELPARAMETER(ZNI);
   const auto UHI0 = MODELPARAMETER(UHI0);

   const std::complex<double> left = 0.5477225575051661*g1*Conj(ZN(gt1,0))*SUM(j1,0,1,Conj(ZNI(gt2,j1))*UHI0(gt3,j1)) - 0.7071067811865475*g2*Conj(ZN(gt1,1))*SUM(j1,0,1,Conj(ZNI(gt2,j1))*UHI0(gt3,j1)) + 0.6708203932499369*gN*Conj(ZN(gt1,5))*SUM(j1,0,1,Conj(ZNI(gt2,j1))*UHI0(gt3,j1)) + Conj(ZN(gt1,4))*SUM(j1,0,1,Conj(ZNI(gt2,j1))*UHI0(gt3,2 + j1)*Lambda12(j1,j1));

   const std::complex<double> right = SUM(j1,0,1,Conj(Lambda12(j1,j1))*UHI0(gt3,j1)*ZNI(gt2,2 + j1))*ZN(gt1,4) + SUM(j1,0,1,UHI0(gt3,2 + j1)*ZNI(gt2,2 + j1))*(-0.5477225575051661*g1*ZN(gt1,0) + 0.7071067811865475*g2*ZN(gt1,1) + 0.4472135954999579*gN*ZN(gt1,5));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha>::type, fields::ChaI, typename fields::conj<fields::SHI0>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto UP = MODELPARAMETER(UP);
   const auto ZMI = MODELPARAMETER(ZMI);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto ZPI = MODELPARAMETER(ZPI);
   const auto UM = MODELPARAMETER(UM);

   const std::complex<double> left = -(g2*Conj(UP(gt1,0))*SUM(j1,0,1,Conj(ZMI(gt2,j1))*UHI0(gt3,j1)));

   const std::complex<double> right = -(g2*SUM(j1,0,1,UHI0(gt3,2 + j1)*ZPI(gt2,j1))*UM(gt1,0));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::ChaI>::type, fields::Cha, fields::SHI0>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto UM = MODELPARAMETER(UM);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto ZPI = MODELPARAMETER(ZPI);
   const auto ZMI = MODELPARAMETER(ZMI);
   const auto UP = MODELPARAMETER(UP);

   const std::complex<double> left = -(g2*Conj(UM(gt2,0))*SUM(j1,0,1,Conj(UHI0(gt3,2 + j1))*Conj(ZPI(gt1,j1))));

   const std::complex<double> right = -(g2*SUM(j1,0,1,Conj(UHI0(gt3,j1))*ZMI(gt1,j1))*UP(gt2,0));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Chi, fields::ChiI, fields::SHI0>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto ZN = MODELPARAMETER(ZN);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto ZNI = MODELPARAMETER(ZNI);

   const std::complex<double> left = -0.5477225575051661*g1*Conj(ZN(gt1,0))*SUM(j1,0,1,Conj(UHI0(gt3,2 + j1))*Conj(ZNI(gt2,2 + j1))) + 0.7071067811865475*g2*Conj(ZN(gt1,1))*SUM(j1,0,1,Conj(UHI0(gt3,2 + j1))*Conj(ZNI(gt2,2 + j1))) + 0.4472135954999579*gN*Conj(ZN(gt1,5))*SUM(j1,0,1,Conj(UHI0(gt3,2 + j1))*Conj(ZNI(gt2,2 + j1))) + Conj(ZN(gt1,4))*SUM(j1,0,1,Conj(UHI0(gt3,j1))*Conj(ZNI(gt2,2 + j1))*Lambda12(j1,j1));

   const std::complex<double> right = SUM(j1,0,1,Conj(UHI0(gt3,2 + j1))*Conj(Lambda12(j1,j1))*ZNI(gt2,j1))*ZN(gt1,4) + 0.1*SUM(j1,0,1,Conj(UHI0(gt3,j1))*ZNI(gt2,j1))*(5.477225575051661*g1*ZN(gt1,0) - 7.0710678118654755*g2*ZN(gt1,1) + 6.708203932499369*gN*ZN(gt1,5));

   return {left, right};
}

ScalarVertex VertexImpl<fields::SHI0, fields::SHI0, typename fields::conj<fields::SHI0>::type, typename fields::conj<fields::SHI0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto UHI0 = MODELPARAMETER(UHI0);

   const std::complex<double> result = 0.0125*(-6*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt4,j1))*SUM(j2,0,1,Conj(UHI0(gt2,j2))*UHI0(gt3,j2)) - 10*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt4,j1))*SUM(j2,0,1,Conj(UHI0(gt2,j2))*UHI0(gt3,j2)) - 9*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt4,j1))*SUM(j2,0,1,Conj(UHI0(gt2,j2))*UHI0(gt3,j2)) + 6*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt4,2 + j1))*SUM(j2,0,1,Conj(UHI0(gt2,j2))*UHI0(gt3,j2)) + 10*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt4,2 + j1))*SUM(j2,0,1,Conj(UHI0(gt2,j2))*UHI0(gt3,j2)) - 6*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt4,2 + j1))*SUM(j2,0,1,Conj(UHI0(gt2,j2))*UHI0(gt3,j2)) - 80*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt4,2 + j1)*Lambda12(j1,j1))*SUM(j2,0,1,Conj(UHI0(gt1,2 + j2))*Conj(Lambda12(j2,j2))*UHI0(gt3,j2)) - 80*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt4,2 + j1)*Lambda12(j1,j1))*SUM(j2,0,1,Conj(UHI0(gt2,2 + j2))*Conj(Lambda12(j2,j2))*UHI0(gt3,j2)) - SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt4,j1))*((6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*SUM(j2,0,1,Conj(UHI0(gt1,j2))*UHI0(gt3,j2)) - 2*(3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*SUM(j2,0,1,Conj(UHI0(gt1,2 + j2))*UHI0(gt3,2 + j2))) + 2*SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*UHI0(gt4,2 + j1))*((3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*SUM(j2,0,1,Conj(UHI0(gt1,j2))*UHI0(gt3,j2)) - (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*SUM(j2,0,1,Conj(UHI0(gt1,2 + j2))*UHI0(gt3,2 + j2))) + 6*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt4,j1))*SUM(j2,0,1,Conj(UHI0(gt2,2 + j2))*UHI0(gt3,2 + j2)) + 10*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt4,j1))*SUM(j2,0,1,Conj(UHI0(gt2,2 + j2))*UHI0(gt3,2 + j2)) - 6*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt4,j1))*SUM(j2,0,1,Conj(UHI0(gt2,2 + j2))*UHI0(gt3,2 + j2)) - 6*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt4,2 + j1))*SUM(j2,0,1,Conj(UHI0(gt2,2 + j2))*UHI0(gt3,2 + j2)) - 10*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt4,2 + j1))*SUM(j2,0,1,Conj(UHI0(gt2,2 + j2))*UHI0(gt3,2 + j2)) - 4*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt4,2 + j1))*SUM(j2,0,1,Conj(UHI0(gt2,2 + j2))*UHI0(gt3,2 + j2)) - 6*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt3,j1))*SUM(j2,0,1,Conj(UHI0(gt1,j2))*UHI0(gt4,j2)) - 10*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt3,j1))*SUM(j2,0,1,Conj(UHI0(gt1,j2))*UHI0(gt4,j2)) - 9*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt3,j1))*SUM(j2,0,1,Conj(UHI0(gt1,j2))*UHI0(gt4,j2)) + 6*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*UHI0(gt3,2 + j1))*SUM(j2,0,1,Conj(UHI0(gt1,j2))*UHI0(gt4,j2)) + 10*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*UHI0(gt3,2 + j1))*SUM(j2,0,1,Conj(UHI0(gt1,j2))*UHI0(gt4,j2)) - 6*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*UHI0(gt3,2 + j1))*SUM(j2,0,1,Conj(UHI0(gt1,j2))*UHI0(gt4,j2)) - 6*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt3,j1))*SUM(j2,0,1,Conj(UHI0(gt2,j2))*UHI0(gt4,j2)) - 10*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt3,j1))*SUM(j2,0,1,Conj(UHI0(gt2,j2))*UHI0(gt4,j2)) - 9*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt3,j1))*SUM(j2,0,1,Conj(UHI0(gt2,j2))*UHI0(gt4,j2)) + 6*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt3,2 + j1))*SUM(j2,0,1,Conj(UHI0(gt2,j2))*UHI0(gt4,j2)) + 10*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt3,2 + j1))*SUM(j2,0,1,Conj(UHI0(gt2,j2))*UHI0(gt4,j2)) - 6*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt3,2 + j1))*SUM(j2,0,1,Conj(UHI0(gt2,j2))*UHI0(gt4,j2)) - 80*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt3,2 + j1)*Lambda12(j1,j1))*SUM(j2,0,1,Conj(UHI0(gt1,2 + j2))*Conj(Lambda12(j2,j2))*UHI0(gt4,j2)) - 80*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt3,2 + j1)*Lambda12(j1,j1))*SUM(j2,0,1,Conj(UHI0(gt2,2 + j2))*Conj(Lambda12(j2,j2))*UHI0(gt4,j2)) + 6*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt3,j1))*SUM(j2,0,1,Conj(UHI0(gt1,2 + j2))*UHI0(gt4,2 + j2)) + 10*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt3,j1))*SUM(j2,0,1,Conj(UHI0(gt1,2 + j2))*UHI0(gt4,2 + j2)) - 6*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt2,j1))*UHI0(gt3,j1))*SUM(j2,0,1,Conj(UHI0(gt1,2 + j2))*UHI0(gt4,2 + j2)) - 6*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*UHI0(gt3,2 + j1))*SUM(j2,0,1,Conj(UHI0(gt1,2 + j2))*UHI0(gt4,2 + j2)) - 10*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*UHI0(gt3,2 + j1))*SUM(j2,0,1,Conj(UHI0(gt1,2 + j2))*UHI0(gt4,2 + j2)) - 4*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt2,2 + j1))*UHI0(gt3,2 + j1))*SUM(j2,0,1,Conj(UHI0(gt1,2 + j2))*UHI0(gt4,2 + j2)) + 6*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt3,j1))*SUM(j2,0,1,Conj(UHI0(gt2,2 + j2))*UHI0(gt4,2 + j2)) + 10*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt3,j1))*SUM(j2,0,1,Conj(UHI0(gt2,2 + j2))*UHI0(gt4,2 + j2)) - 6*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt3,j1))*SUM(j2,0,1,Conj(UHI0(gt2,2 + j2))*UHI0(gt4,2 + j2)) - 6*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt3,2 + j1))*SUM(j2,0,1,Conj(UHI0(gt2,2 + j2))*UHI0(gt4,2 + j2)) - 10*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt3,2 + j1))*SUM(j2,0,1,Conj(UHI0(gt2,2 + j2))*UHI0(gt4,2 + j2)) - 4*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt3,2 + j1))*SUM(j2,0,1,Conj(UHI0(gt2,2 + j2))*UHI0(gt4,2 + j2)));

   return {result};
}

ScalarVertex VertexImpl<fields::SHI0, fields::SHIp, typename fields::conj<fields::SHI0>::type, typename fields::conj<fields::SHIp>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto UHIp = MODELPARAMETER(UHIp);

   const std::complex<double> result = 0.0125*(-20*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHIp(gt4,j1))*SUM(j2,0,1,Conj(UHIp(gt2,j2))*UHI0(gt3,j2)) - 20*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHIp(gt4,2 + j1))*SUM(j2,0,1,Conj(UHIp(gt2,j2))*UHI0(gt3,j2)) + 80*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt4,2 + j1)*Lambda12(j1,j1))*SUM(j2,0,1,Conj(UHI0(gt1,2 + j2))*Conj(Lambda12(j2,j2))*UHI0(gt3,j2)) + 2*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*UHIp(gt4,2 + j1))*((3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*SUM(j2,0,1,Conj(UHI0(gt1,j2))*UHI0(gt3,j2)) + (-3*Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN))*SUM(j2,0,1,Conj(UHI0(gt1,2 + j2))*UHI0(gt3,2 + j2))) - SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt4,j1))*((6*Sqr(g1) - 10*Sqr(g2) + 9*Sqr(gN))*SUM(j2,0,1,Conj(UHI0(gt1,j2))*UHI0(gt3,j2)) + 2*(-3*Sqr(g1) + 5*Sqr(g2) + 3*Sqr(gN))*SUM(j2,0,1,Conj(UHI0(gt1,2 + j2))*UHI0(gt3,2 + j2))) - 20*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHIp(gt4,j1))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*UHI0(gt3,2 + j2)) - 20*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHIp(gt4,2 + j1))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*UHI0(gt3,2 + j2)) - 20*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHI0(gt3,j1))*SUM(j2,0,1,Conj(UHI0(gt1,j2))*UHIp(gt4,j2)) - 20*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*UHI0(gt3,2 + j1))*SUM(j2,0,1,Conj(UHI0(gt1,j2))*UHIp(gt4,j2)) - 6*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt3,j1))*SUM(j2,0,1,Conj(UHIp(gt2,j2))*UHIp(gt4,j2)) + 10*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt3,j1))*SUM(j2,0,1,Conj(UHIp(gt2,j2))*UHIp(gt4,j2)) - 9*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt3,j1))*SUM(j2,0,1,Conj(UHIp(gt2,j2))*UHIp(gt4,j2)) + 6*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt3,2 + j1))*SUM(j2,0,1,Conj(UHIp(gt2,j2))*UHIp(gt4,j2)) - 10*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt3,2 + j1))*SUM(j2,0,1,Conj(UHIp(gt2,j2))*UHIp(gt4,j2)) - 6*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt3,2 + j1))*SUM(j2,0,1,Conj(UHIp(gt2,j2))*UHIp(gt4,j2)) + 80*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt3,2 + j1)*Lambda12(j1,j1))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*Conj(Lambda12(j2,j2))*UHIp(gt4,j2)) - 20*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHI0(gt3,j1))*SUM(j2,0,1,Conj(UHI0(gt1,2 + j2))*UHIp(gt4,2 + j2)) - 20*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*UHI0(gt3,2 + j1))*SUM(j2,0,1,Conj(UHI0(gt1,2 + j2))*UHIp(gt4,2 + j2)) + 6*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt3,j1))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*UHIp(gt4,2 + j2)) - 10*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt3,j1))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*UHIp(gt4,2 + j2)) - 6*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt3,j1))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*UHIp(gt4,2 + j2)) - 6*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt3,2 + j1))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*UHIp(gt4,2 + j2)) + 10*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt3,2 + j1))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*UHIp(gt4,2 + j2)) - 4*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt3,2 + j1))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*UHIp(gt4,2 + j2)));

   return {result};
}

ScalarVertex VertexImpl<fields::SHI0, fields::SHp0, typename fields::conj<fields::SHI0>::type, typename fields::conj<fields::SHp0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto UHp0 = MODELPARAMETER(UHp0);

   const std::complex<double> result = -0.05*((3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt3,j1)) - (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt3,2 + j1)))*(Conj(UHp0(gt2,0))*UHp0(gt4,0) - Conj(UHp0(gt2,1))*UHp0(gt4,1));

   return {result};
}

ScalarVertex VertexImpl<fields::SHI0, fields::SHpp, typename fields::conj<fields::SHI0>::type, typename fields::conj<fields::SHpp>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto UHpp = MODELPARAMETER(UHpp);

   const std::complex<double> result = -0.05*((3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt3,j1)) + (-3*Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN))*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt3,2 + j1)))*(Conj(UHpp(gt2,0))*UHpp(gt4,0) - Conj(UHpp(gt2,1))*UHpp(gt4,1));

   return {result};
}

ScalarVertex VertexImpl<fields::SHI0, fields::SSI0, typename fields::conj<fields::SHI0>::type, typename fields::conj<fields::SSI0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto gN = MODELPARAMETER(gN);
   const auto UHI0 = MODELPARAMETER(UHI0);

   const std::complex<double> result = 0.125*KroneckerDelta(gt2,gt4)*Sqr(gN)*(3*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt3,j1)) + 2*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt3,2 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHI0, typename fields::conj<fields::SHI0>::type, typename fields::conj<fields::VWm>::type, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> result = 0.5*KroneckerDelta(gt1,gt2)*Sqr(g2);

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHI0, typename fields::conj<fields::SHI0>::type, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05*((-14.696938456699067*g1*gN*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) + 10*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Cos(ThetaWp)) + Cos(ThetaW)*(-9.486832980505138*g2*gN*Sin(2*ThetaWp) + 15.491933384829668*g1*g2*Sin(ThetaW)*Sqr(Cos(ThetaWp))) + 6*Sqr(g1)*Sqr(Cos(ThetaWp))*Sqr(Sin(ThetaW)) + 9*Sqr(gN)*Sqr(Sin(ThetaWp)))*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt2,j1)) + 2*(3.1622776601683795*g2*gN*Cos(ThetaW)*Sin(2*ThetaWp) + g1*Sin(ThetaW)*(7.745966692414834*g2*Cos(ThetaW) + 3*g1*Sin(ThetaW))*Sqr(Cos(ThetaWp)) + 5*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Cos(ThetaWp)) + gN*(2.449489742783178*g1*Sin(ThetaW)*Sin(2*ThetaWp) + 2*gN*Sqr(Sin(ThetaWp))))*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt2,2 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHI0, typename fields::conj<fields::SHI0>::type, fields::VZp, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*((9*Sqr(gN)*Sqr(Cos(ThetaWp)) + (g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(9.486832980505138*gN*Sin(2*ThetaWp) + 10*g2*Cos(ThetaW)*Sqr(Sin(ThetaWp)) + 7.745966692414834*g1*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt2,j1)) + 2*(-3.1622776601683795*gN*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*Sin(2*ThetaWp) + 2*Sqr(gN)*Sqr(Cos(ThetaWp)) + 5*Sqr(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*Sqr(Sin(ThetaWp)))*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt2,2 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHI0, typename fields::conj<fields::SHI0>::type, fields::VZ, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05*(-((9.486832980505138*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) - 9*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN) + 5*Sin(2*ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 7.348469228349534*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + g1*(3.872983346207417*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) + 3*g1*Sin(2*ThetaWp)*Sqr(Sin(ThetaW)) - 7.348469228349534*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt2,j1))) + (6.324555320336759*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) + 4*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN) - 5*Sin(2*ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 4.898979485566356*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) - g1*(3.872983346207417*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) + 3*g1*Sin(2*ThetaWp)*Sqr(Sin(ThetaW)) + 4.898979485566356*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt2,2 + j1)));

   return {result};
}

ChiralVertex VertexImpl<fields::Chi, fields::ChaI, typename fields::conj<fields::SHIp>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZMI = MODELPARAMETER(ZMI);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto ZPI = MODELPARAMETER(ZPI);

   const std::complex<double> left = 0.5477225575051661*g1*Conj(ZN(gt1,0))*SUM(j1,0,1,Conj(ZMI(gt2,j1))*UHIp(gt3,j1)) + 0.7071067811865475*g2*Conj(ZN(gt1,1))*SUM(j1,0,1,Conj(ZMI(gt2,j1))*UHIp(gt3,j1)) + 0.6708203932499369*gN*Conj(ZN(gt1,5))*SUM(j1,0,1,Conj(ZMI(gt2,j1))*UHIp(gt3,j1)) - Conj(ZN(gt1,4))*SUM(j1,0,1,Conj(ZMI(gt2,j1))*UHIp(gt3,2 + j1)*Lambda12(j1,j1));

   const std::complex<double> right = -(SUM(j1,0,1,Conj(Lambda12(j1,j1))*UHIp(gt3,j1)*ZPI(gt2,j1))*ZN(gt1,4)) - 0.1*SUM(j1,0,1,UHIp(gt3,2 + j1)*ZPI(gt2,j1))*(5.477225575051661*g1*ZN(gt1,0) + 7.0710678118654755*g2*ZN(gt1,1) - 4.47213595499958*gN*ZN(gt1,5));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::ChaI>::type, fields::Chi, fields::SHIp>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto ZN = MODELPARAMETER(ZN);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto ZPI = MODELPARAMETER(ZPI);
   const auto ZMI = MODELPARAMETER(ZMI);

   const std::complex<double> left = -0.5477225575051661*g1*Conj(ZN(gt2,0))*SUM(j1,0,1,Conj(UHIp(gt3,2 + j1))*Conj(ZPI(gt1,j1))) - 0.7071067811865475*g2*Conj(ZN(gt2,1))*SUM(j1,0,1,Conj(UHIp(gt3,2 + j1))*Conj(ZPI(gt1,j1))) + 0.4472135954999579*gN*Conj(ZN(gt2,5))*SUM(j1,0,1,Conj(UHIp(gt3,2 + j1))*Conj(ZPI(gt1,j1))) - Conj(ZN(gt2,4))*SUM(j1,0,1,Conj(UHIp(gt3,j1))*Conj(ZPI(gt1,j1))*Lambda12(j1,j1));

   const std::complex<double> right = -(SUM(j1,0,1,Conj(UHIp(gt3,2 + j1))*Conj(Lambda12(j1,j1))*ZMI(gt1,j1))*ZN(gt2,4)) + 0.1*SUM(j1,0,1,Conj(UHIp(gt3,j1))*ZMI(gt1,j1))*(5.477225575051661*g1*ZN(gt2,0) + 7.0710678118654755*g2*ZN(gt2,1) + 6.708203932499369*gN*ZN(gt2,5));

   return {left, right};
}

ChiralVertex VertexImpl<fields::ChiI, fields::Cha, typename fields::conj<fields::SHIp>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto UM = MODELPARAMETER(UM);
   const auto ZNI = MODELPARAMETER(ZNI);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto UP = MODELPARAMETER(UP);

   const std::complex<double> left = -(g2*Conj(UM(gt2,0))*SUM(j1,0,1,Conj(ZNI(gt1,j1))*UHIp(gt3,j1)));

   const std::complex<double> right = -(g2*SUM(j1,0,1,UHIp(gt3,2 + j1)*ZNI(gt1,2 + j1))*UP(gt2,0));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha>::type, fields::ChiI, fields::SHIp>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto UP = MODELPARAMETER(UP);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto ZNI = MODELPARAMETER(ZNI);
   const auto UM = MODELPARAMETER(UM);

   const std::complex<double> left = -(g2*Conj(UP(gt1,0))*SUM(j1,0,1,Conj(UHIp(gt3,2 + j1))*Conj(ZNI(gt2,2 + j1))));

   const std::complex<double> right = -(g2*SUM(j1,0,1,Conj(UHIp(gt3,j1))*ZNI(gt2,j1))*UM(gt1,0));

   return {left, right};
}

ScalarVertex VertexImpl<fields::SHIp, fields::SHIp, typename fields::conj<fields::SHIp>::type, typename fields::conj<fields::SHIp>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto Lambda12 = MODELPARAMETER(Lambda12);
   const auto UHIp = MODELPARAMETER(UHIp);

   const std::complex<double> result = 0.0125*(-6*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt4,j1))*SUM(j2,0,1,Conj(UHIp(gt2,j2))*UHIp(gt3,j2)) - 10*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt4,j1))*SUM(j2,0,1,Conj(UHIp(gt2,j2))*UHIp(gt3,j2)) - 9*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt4,j1))*SUM(j2,0,1,Conj(UHIp(gt2,j2))*UHIp(gt3,j2)) + 6*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt4,2 + j1))*SUM(j2,0,1,Conj(UHIp(gt2,j2))*UHIp(gt3,j2)) + 10*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt4,2 + j1))*SUM(j2,0,1,Conj(UHIp(gt2,j2))*UHIp(gt3,j2)) - 6*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt4,2 + j1))*SUM(j2,0,1,Conj(UHIp(gt2,j2))*UHIp(gt3,j2)) - 80*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt4,2 + j1)*Lambda12(j1,j1))*SUM(j2,0,1,Conj(UHIp(gt1,2 + j2))*Conj(Lambda12(j2,j2))*UHIp(gt3,j2)) - 80*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt4,2 + j1)*Lambda12(j1,j1))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*Conj(Lambda12(j2,j2))*UHIp(gt3,j2)) - SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt4,j1))*((6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*SUM(j2,0,1,Conj(UHIp(gt1,j2))*UHIp(gt3,j2)) - 2*(3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*SUM(j2,0,1,Conj(UHIp(gt1,2 + j2))*UHIp(gt3,2 + j2))) + 2*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*UHIp(gt4,2 + j1))*((3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*SUM(j2,0,1,Conj(UHIp(gt1,j2))*UHIp(gt3,j2)) - (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*SUM(j2,0,1,Conj(UHIp(gt1,2 + j2))*UHIp(gt3,2 + j2))) + 6*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt4,j1))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*UHIp(gt3,2 + j2)) + 10*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt4,j1))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*UHIp(gt3,2 + j2)) - 6*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt4,j1))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*UHIp(gt3,2 + j2)) - 6*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt4,2 + j1))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*UHIp(gt3,2 + j2)) - 10*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt4,2 + j1))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*UHIp(gt3,2 + j2)) - 4*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt4,2 + j1))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*UHIp(gt3,2 + j2)) - 6*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt3,j1))*SUM(j2,0,1,Conj(UHIp(gt1,j2))*UHIp(gt4,j2)) - 10*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt3,j1))*SUM(j2,0,1,Conj(UHIp(gt1,j2))*UHIp(gt4,j2)) - 9*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt3,j1))*SUM(j2,0,1,Conj(UHIp(gt1,j2))*UHIp(gt4,j2)) + 6*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*UHIp(gt3,2 + j1))*SUM(j2,0,1,Conj(UHIp(gt1,j2))*UHIp(gt4,j2)) + 10*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*UHIp(gt3,2 + j1))*SUM(j2,0,1,Conj(UHIp(gt1,j2))*UHIp(gt4,j2)) - 6*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*UHIp(gt3,2 + j1))*SUM(j2,0,1,Conj(UHIp(gt1,j2))*UHIp(gt4,j2)) - 6*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt3,j1))*SUM(j2,0,1,Conj(UHIp(gt2,j2))*UHIp(gt4,j2)) - 10*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt3,j1))*SUM(j2,0,1,Conj(UHIp(gt2,j2))*UHIp(gt4,j2)) - 9*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt3,j1))*SUM(j2,0,1,Conj(UHIp(gt2,j2))*UHIp(gt4,j2)) + 6*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt3,2 + j1))*SUM(j2,0,1,Conj(UHIp(gt2,j2))*UHIp(gt4,j2)) + 10*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt3,2 + j1))*SUM(j2,0,1,Conj(UHIp(gt2,j2))*UHIp(gt4,j2)) - 6*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt3,2 + j1))*SUM(j2,0,1,Conj(UHIp(gt2,j2))*UHIp(gt4,j2)) - 80*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt3,2 + j1)*Lambda12(j1,j1))*SUM(j2,0,1,Conj(UHIp(gt1,2 + j2))*Conj(Lambda12(j2,j2))*UHIp(gt4,j2)) - 80*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt3,2 + j1)*Lambda12(j1,j1))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*Conj(Lambda12(j2,j2))*UHIp(gt4,j2)) + 6*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt3,j1))*SUM(j2,0,1,Conj(UHIp(gt1,2 + j2))*UHIp(gt4,2 + j2)) + 10*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt3,j1))*SUM(j2,0,1,Conj(UHIp(gt1,2 + j2))*UHIp(gt4,2 + j2)) - 6*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt2,j1))*UHIp(gt3,j1))*SUM(j2,0,1,Conj(UHIp(gt1,2 + j2))*UHIp(gt4,2 + j2)) - 6*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*UHIp(gt3,2 + j1))*SUM(j2,0,1,Conj(UHIp(gt1,2 + j2))*UHIp(gt4,2 + j2)) - 10*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*UHIp(gt3,2 + j1))*SUM(j2,0,1,Conj(UHIp(gt1,2 + j2))*UHIp(gt4,2 + j2)) - 4*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt2,2 + j1))*UHIp(gt3,2 + j1))*SUM(j2,0,1,Conj(UHIp(gt1,2 + j2))*UHIp(gt4,2 + j2)) + 6*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt3,j1))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*UHIp(gt4,2 + j2)) + 10*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt3,j1))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*UHIp(gt4,2 + j2)) - 6*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt3,j1))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*UHIp(gt4,2 + j2)) - 6*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt3,2 + j1))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*UHIp(gt4,2 + j2)) - 10*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt3,2 + j1))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*UHIp(gt4,2 + j2)) - 4*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt3,2 + j1))*SUM(j2,0,1,Conj(UHIp(gt2,2 + j2))*UHIp(gt4,2 + j2)));

   return {result};
}

ScalarVertex VertexImpl<fields::SHIp, fields::SHp0, typename fields::conj<fields::SHIp>::type, typename fields::conj<fields::SHp0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto UHp0 = MODELPARAMETER(UHp0);

   const std::complex<double> result = -0.05*((3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt3,j1)) + (-3*Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN))*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt3,2 + j1)))*(Conj(UHp0(gt2,0))*UHp0(gt4,0) - Conj(UHp0(gt2,1))*UHp0(gt4,1));

   return {result};
}

ScalarVertex VertexImpl<fields::SHIp, fields::SHpp, typename fields::conj<fields::SHIp>::type, typename fields::conj<fields::SHpp>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto UHpp = MODELPARAMETER(UHpp);

   const std::complex<double> result = -0.05*((3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt3,j1)) - (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt3,2 + j1)))*(Conj(UHpp(gt2,0))*UHpp(gt4,0) - Conj(UHpp(gt2,1))*UHpp(gt4,1));

   return {result};
}

ScalarVertex VertexImpl<fields::SHIp, fields::SSI0, typename fields::conj<fields::SHIp>::type, typename fields::conj<fields::SSI0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto gN = MODELPARAMETER(gN);
   const auto UHIp = MODELPARAMETER(UHIp);

   const std::complex<double> result = 0.125*KroneckerDelta(gt2,gt4)*Sqr(gN)*(3*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt3,j1)) + 2*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt3,2 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHIp, typename fields::conj<fields::SHIp>::type, fields::VP, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*KroneckerDelta(gt1,gt2)*Sqr(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHIp, typename fields::conj<fields::SHIp>::type, typename fields::conj<fields::VWm>::type, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> result = 0.5*KroneckerDelta(gt1,gt2)*Sqr(g2);

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHIp, typename fields::conj<fields::SHIp>::type, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05*((-14.696938456699067*g1*gN*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) + 2*g2*Cos(ThetaW)*Cos(ThetaWp)*(-7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp)) + 10*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Cos(ThetaWp)) + 6*Sqr(g1)*Sqr(Cos(ThetaWp))*Sqr(Sin(ThetaW)) + 9*Sqr(gN)*Sqr(Sin(ThetaWp)))*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt2,j1)) + 2*(5*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Cos(ThetaWp)) - 2*Cos(ThetaW)*(3.1622776601683795*g2*gN*Cos(ThetaWp)*Sin(ThetaWp) + 3.872983346207417*g1*g2*Sin(ThetaW)*Sqr(Cos(ThetaWp))) + 3*Sqr(g1)*Sqr(Cos(ThetaWp))*Sqr(Sin(ThetaW)) + gN*(2.449489742783178*g1*Sin(ThetaW)*Sin(2*ThetaWp) + 2*gN*Sqr(Sin(ThetaWp))))*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt2,2 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHIp, typename fields::conj<fields::SHIp>::type, fields::VZp, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*((9*Sqr(gN)*Sqr(Cos(ThetaWp)) + (g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*(-9.486832980505138*gN*Sin(2*ThetaWp) + 10*g2*Cos(ThetaW)*Sqr(Sin(ThetaWp)) - 7.745966692414834*g1*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt2,j1)) + 2*(-4.898979485566356*g1*gN*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) + 2*Sqr(gN)*Sqr(Cos(ThetaWp)) + 5*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Sin(ThetaWp)) + 3*Sqr(g1)*Sqr(Sin(ThetaW))*Sqr(Sin(ThetaWp)) + Cos(ThetaW)*(3.1622776601683795*g2*gN*Sin(2*ThetaWp) - 7.745966692414834*g1*g2*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt2,2 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHIp, typename fields::conj<fields::SHIp>::type, fields::VZ, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05*((9.486832980505138*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) + 9*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN) - 5*Sin(2*ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) - 7.348469228349534*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + g1*(3.872983346207417*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) - 3*g1*Sin(2*ThetaWp)*Sqr(Sin(ThetaW)) + 7.348469228349534*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt2,j1)) + (-6.324555320336759*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) + 4*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN) - 5*Sin(2*ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 4.898979485566356*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + g1*(3.872983346207417*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) - 3*g1*Sin(2*ThetaWp)*Sqr(Sin(ThetaW)) - 4.898979485566356*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt2,2 + j1)));

   return {result};
}

ChiralVertex VertexImpl<fields::Chi, fields::FSI, typename fields::conj<fields::SSI0>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto gN = MODELPARAMETER(gN);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZFSI = MODELPARAMETER(ZFSI);
   const auto ZSSI = MODELPARAMETER(ZSSI);

   const std::complex<double> left = -1.118033988749895*gN*Conj(ZN(gt1,5))*SUM(j1,0,1,Conj(ZFSI(gt2,j1))*ZSSI(gt3,j1));

   const std::complex<double> right = 0;

   return {left, right};
}

ChiralVertex VertexImpl<fields::Chi, fields::FSI, fields::SSI0>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto gN = MODELPARAMETER(gN);
   const auto ZSSI = MODELPARAMETER(ZSSI);
   const auto ZFSI = MODELPARAMETER(ZFSI);
   const auto ZN = MODELPARAMETER(ZN);

   const std::complex<double> left = 0;

   const std::complex<double> right = -1.118033988749895*gN*SUM(j1,0,1,Conj(ZSSI(gt3,j1))*ZFSI(gt2,j1))*ZN(gt1,5);

   return {left, right};
}

ScalarVertex VertexImpl<fields::SHp0, fields::SSI0, typename fields::conj<fields::SHp0>::type, typename fields::conj<fields::SSI0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto gN = MODELPARAMETER(gN);
   const auto UHp0 = MODELPARAMETER(UHp0);

   const std::complex<double> result = -0.25*KroneckerDelta(gt2,gt4)*Sqr(gN)*(Conj(UHp0(gt1,0))*UHp0(gt3,0) - Conj(UHp0(gt1,1))*UHp0(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::SHpp, fields::SSI0, typename fields::conj<fields::SHpp>::type, typename fields::conj<fields::SSI0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto gN = MODELPARAMETER(gN);
   const auto UHpp = MODELPARAMETER(UHpp);

   const std::complex<double> result = -0.25*KroneckerDelta(gt2,gt4)*Sqr(gN)*(Conj(UHpp(gt1,0))*UHpp(gt3,0) - Conj(UHpp(gt1,1))*UHpp(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::SSI0, fields::SSI0, typename fields::conj<fields::SSI0>::type, typename fields::conj<fields::SSI0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto gN = MODELPARAMETER(gN);

   const std::complex<double> result = -0.625*(KroneckerDelta(gt1,gt4)*KroneckerDelta(gt2,gt3) + KroneckerDelta(gt1,gt3)*KroneckerDelta(gt2,gt4))*Sqr(gN);

   return {result};
}

InverseMetricVertex VertexImpl<fields::SSI0, typename fields::conj<fields::SSI0>::type, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 1.25*KroneckerDelta(gt1,gt2)*Sqr(gN)*Sqr(Sin(ThetaWp));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SSI0, typename fields::conj<fields::SSI0>::type, fields::VZp, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 1.25*KroneckerDelta(gt1,gt2)*Sqr(gN)*Sqr(Cos(ThetaWp));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SSI0, typename fields::conj<fields::SSI0>::type, fields::VZ, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 1.25*Cos(ThetaWp)*KroneckerDelta(gt1,gt2)*Sin(ThetaWp)*Sqr(gN);

   return {result};
}

ChiralVertex VertexImpl<fields::Chi, fields::ChiP, typename fields::conj<fields::SHp0>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZNp = MODELPARAMETER(ZNp);
   const auto UHp0 = MODELPARAMETER(UHp0);

   const std::complex<double> left = 0.1*(5.477225575051661*g1*Conj(ZN(gt1,0)) - 7.0710678118654755*g2*Conj(ZN(gt1,1)) - 4.47213595499958*gN*Conj(ZN(gt1,5)))*Conj(ZNp(gt2,0))*UHp0(gt3,0);

   const std::complex<double> right = -0.1*UHp0(gt3,1)*(5.477225575051661*g1*ZN(gt1,0) - 7.0710678118654755*g2*ZN(gt1,1) - 4.47213595499958*gN*ZN(gt1,5))*ZNp(gt2,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha>::type, fields::ChaP, typename fields::conj<fields::SHp0>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto UP = MODELPARAMETER(UP);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto UM = MODELPARAMETER(UM);

   const std::complex<double> left = -(g2*Conj(UP(gt1,0))*UHp0(gt3,0));

   const std::complex<double> right = -(g2*UHp0(gt3,1)*UM(gt1,0));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::ChaP>::type, fields::Cha, fields::SHp0>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt3 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto UM = MODELPARAMETER(UM);
   const auto UP = MODELPARAMETER(UP);

   const std::complex<double> left = -(g2*Conj(UHp0(gt3,1))*Conj(UM(gt2,0)));

   const std::complex<double> right = -(g2*Conj(UHp0(gt3,0))*UP(gt2,0));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Chi, fields::ChiP, fields::SHp0>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZNp = MODELPARAMETER(ZNp);

   const std::complex<double> left = -0.1*Conj(UHp0(gt3,1))*(5.477225575051661*g1*Conj(ZN(gt1,0)) - 7.0710678118654755*g2*Conj(ZN(gt1,1)) - 4.47213595499958*gN*Conj(ZN(gt1,5)))*Conj(ZNp(gt2,1));

   const std::complex<double> right = 0.1*Conj(UHp0(gt3,0))*(5.477225575051661*g1*ZN(gt1,0) - 7.0710678118654755*g2*ZN(gt1,1) - 4.47213595499958*gN*ZN(gt1,5))*ZNp(gt2,0);

   return {left, right};
}

ScalarVertex VertexImpl<fields::SHp0, fields::SHp0, typename fields::conj<fields::SHp0>::type, typename fields::conj<fields::SHp0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto UHp0 = MODELPARAMETER(UHp0);

   const std::complex<double> result = -0.05*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*(-(Conj(UHp0(gt1,1))*(-2*Conj(UHp0(gt2,1))*UHp0(gt3,1)*UHp0(gt4,1) + Conj(UHp0(gt2,0))*(UHp0(gt3,1)*UHp0(gt4,0) + UHp0(gt3,0)*UHp0(gt4,1)))) + Conj(UHp0(gt1,0))*(2*Conj(UHp0(gt2,0))*UHp0(gt3,0)*UHp0(gt4,0) - Conj(UHp0(gt2,1))*(UHp0(gt3,1)*UHp0(gt4,0) + UHp0(gt3,0)*UHp0(gt4,1))));

   return {result};
}

ScalarVertex VertexImpl<fields::SHp0, fields::SHpp, typename fields::conj<fields::SHp0>::type, typename fields::conj<fields::SHpp>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto UHpp = MODELPARAMETER(UHpp);

   const std::complex<double> result = 0.05*(Conj(UHp0(gt1,1))*(-(Conj(UHpp(gt2,1))*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*UHp0(gt3,1)*UHpp(gt4,1)) + Conj(UHpp(gt2,0))*((3*Sqr(g1) - 5*Sqr(g2) + 2*Sqr(gN))*UHp0(gt3,1)*UHpp(gt4,0) - 10*Sqr(g2)*UHp0(gt3,0)*UHpp(gt4,1))) - Conj(UHp0(gt1,0))*(Conj(UHpp(gt2,0))*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*UHp0(gt3,0)*UHpp(gt4,0) + Conj(UHpp(gt2,1))*(10*Sqr(g2)*UHp0(gt3,1)*UHpp(gt4,0) + (-3*Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN))*UHp0(gt3,0)*UHpp(gt4,1))));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHp0, typename fields::conj<fields::SHp0>::type, typename fields::conj<fields::VWm>::type, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto UHp0 = MODELPARAMETER(UHp0);

   const std::complex<double> result = 0.5*Sqr(g2)*(Conj(UHp0(gt1,0))*UHp0(gt2,0) + Conj(UHp0(gt1,1))*UHp0(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHp0, typename fields::conj<fields::SHp0>::type, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -0.0125*(-9.797958971132712*g1*gN*Cos(ThetaW - 2*ThetaWp) + 9.797958971132712*g1*gN*Cos(ThetaW + 2*ThetaWp) - 15.491933384829668*g1*g2*Sin(2*ThetaW) + 12.649110640673518*g2*gN*Sin(ThetaW - 2*ThetaWp) - 7.745966692414834*g1*g2*Sin(2*(ThetaW - ThetaWp)) - 7.745966692414834*g1*g2*Sin(2*(ThetaW + ThetaWp)) - 12.649110640673518*g2*gN*Sin(ThetaW + 2*ThetaWp) - 6*Sqr(g1) + 3*Cos(2*(ThetaW - ThetaWp))*Sqr(g1) - 6*Cos(2*ThetaWp)*Sqr(g1) + 3*Cos(2*(ThetaW + ThetaWp))*Sqr(g1) + 2*Cos(2*ThetaW)*(3*Sqr(g1) - 5*Sqr(g2)) - 10*Sqr(g2) - 5*Cos(2*(ThetaW - ThetaWp))*Sqr(g2) - 10*Cos(2*ThetaWp)*Sqr(g2) - 5*Cos(2*(ThetaW + ThetaWp))*Sqr(g2) - 8*Sqr(gN) + 8*Cos(2*ThetaWp)*Sqr(gN))*(Conj(UHp0(gt1,0))*UHp0(gt2,0) + Conj(UHp0(gt1,1))*UHp0(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHp0, typename fields::conj<fields::SHp0>::type, fields::VZp, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -0.0125*(9.797958971132712*g1*gN*Cos(ThetaW - 2*ThetaWp) - 9.797958971132712*g1*gN*Cos(ThetaW + 2*ThetaWp) - 15.491933384829668*g1*g2*Sin(2*ThetaW) - 12.649110640673518*g2*gN*Sin(ThetaW - 2*ThetaWp) + 7.745966692414834*g1*g2*Sin(2*(ThetaW - ThetaWp)) + 7.745966692414834*g1*g2*Sin(2*(ThetaW + ThetaWp)) + 12.649110640673518*g2*gN*Sin(ThetaW + 2*ThetaWp) - 6*Sqr(g1) - 3*Cos(2*(ThetaW - ThetaWp))*Sqr(g1) + 6*Cos(2*ThetaWp)*Sqr(g1) - 3*Cos(2*(ThetaW + ThetaWp))*Sqr(g1) + 2*Cos(2*ThetaW)*(3*Sqr(g1) - 5*Sqr(g2)) - 10*Sqr(g2) + 5*Cos(2*(ThetaW - ThetaWp))*Sqr(g2) + 10*Cos(2*ThetaWp)*Sqr(g2) + 5*Cos(2*(ThetaW + ThetaWp))*Sqr(g2) - 8*Sqr(gN) - 8*Cos(2*ThetaWp)*Sqr(gN))*(Conj(UHp0(gt1,0))*UHp0(gt2,0) + Conj(UHp0(gt1,1))*UHp0(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHp0, typename fields::conj<fields::SHp0>::type, fields::VZ, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.0125*(25.298221281347036*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) - 15.491933384829668*g1*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) + Sin(2*ThetaWp)*(-9*Sqr(g1) + Cos(2*ThetaW)*(3*Sqr(g1) - 5*Sqr(g2)) - 5*Sqr(g2)) + 8*Sin(2*ThetaWp)*Sqr(gN) + 2*Sin(2*ThetaWp)*(3*Sqr(g1) - 5*Sqr(g2))*Sqr(Cos(ThetaW)) + 19.595917942265423*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) - 19.595917942265423*g1*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp)))*(Conj(UHp0(gt1,0))*UHp0(gt2,0) + Conj(UHp0(gt1,1))*UHp0(gt2,1));

   return {result};
}

ChiralVertex VertexImpl<fields::Chi, fields::ChaP, typename fields::conj<fields::SHpp>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto ZN = MODELPARAMETER(ZN);
   const auto UHpp = MODELPARAMETER(UHpp);

   const std::complex<double> left = 0.1*(5.477225575051661*g1*Conj(ZN(gt1,0)) + 7.0710678118654755*g2*Conj(ZN(gt1,1)) - 4.47213595499958*gN*Conj(ZN(gt1,5)))*UHpp(gt3,0);

   const std::complex<double> right = -0.1*UHpp(gt3,1)*(5.477225575051661*g1*ZN(gt1,0) + 7.0710678118654755*g2*ZN(gt1,1) - 4.47213595499958*gN*ZN(gt1,5));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::ChaP>::type, fields::Chi, fields::SHpp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt3 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto ZN = MODELPARAMETER(ZN);

   const std::complex<double> left = -0.1*Conj(UHpp(gt3,1))*(5.477225575051661*g1*Conj(ZN(gt2,0)) + 7.0710678118654755*g2*Conj(ZN(gt2,1)) - 4.47213595499958*gN*Conj(ZN(gt2,5)));

   const std::complex<double> right = 0.1*Conj(UHpp(gt3,0))*(5.477225575051661*g1*ZN(gt2,0) + 7.0710678118654755*g2*ZN(gt2,1) - 4.47213595499958*gN*ZN(gt2,5));

   return {left, right};
}

ChiralVertex VertexImpl<fields::ChiP, fields::Cha, typename fields::conj<fields::SHpp>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto UM = MODELPARAMETER(UM);
   const auto ZNp = MODELPARAMETER(ZNp);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto UP = MODELPARAMETER(UP);

   const std::complex<double> left = -(g2*Conj(UM(gt2,0))*Conj(ZNp(gt1,0))*UHpp(gt3,0));

   const std::complex<double> right = -(g2*UHpp(gt3,1)*UP(gt2,0)*ZNp(gt1,1));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha>::type, fields::ChiP, fields::SHpp>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto UP = MODELPARAMETER(UP);
   const auto ZNp = MODELPARAMETER(ZNp);
   const auto UM = MODELPARAMETER(UM);

   const std::complex<double> left = -(g2*Conj(UHpp(gt3,1))*Conj(UP(gt1,0))*Conj(ZNp(gt2,1)));

   const std::complex<double> right = -(g2*Conj(UHpp(gt3,0))*UM(gt1,0)*ZNp(gt2,0));

   return {left, right};
}

ScalarVertex VertexImpl<fields::SHpp, fields::SHpp, typename fields::conj<fields::SHpp>::type, typename fields::conj<fields::SHpp>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto UHpp = MODELPARAMETER(UHpp);

   const std::complex<double> result = -0.05*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*(-(Conj(UHpp(gt1,1))*(-2*Conj(UHpp(gt2,1))*UHpp(gt3,1)*UHpp(gt4,1) + Conj(UHpp(gt2,0))*(UHpp(gt3,1)*UHpp(gt4,0) + UHpp(gt3,0)*UHpp(gt4,1)))) + Conj(UHpp(gt1,0))*(2*Conj(UHpp(gt2,0))*UHpp(gt3,0)*UHpp(gt4,0) - Conj(UHpp(gt2,1))*(UHpp(gt3,1)*UHpp(gt4,0) + UHpp(gt3,0)*UHpp(gt4,1))));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHpp, typename fields::conj<fields::SHpp>::type, fields::VP, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*Sqr(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*(Conj(UHpp(gt1,0))*UHpp(gt2,0) + Conj(UHpp(gt1,1))*UHpp(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHpp, typename fields::conj<fields::SHpp>::type, typename fields::conj<fields::VWm>::type, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto UHpp = MODELPARAMETER(UHpp);

   const std::complex<double> result = 0.5*Sqr(g2)*(Conj(UHpp(gt1,0))*UHpp(gt2,0) + Conj(UHpp(gt1,1))*UHpp(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHpp, typename fields::conj<fields::SHpp>::type, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -0.0125*(-9.797958971132712*g1*gN*Cos(ThetaW - 2*ThetaWp) + 9.797958971132712*g1*gN*Cos(ThetaW + 2*ThetaWp) + 15.491933384829668*g1*g2*Sin(2*ThetaW) - 12.649110640673518*g2*gN*Sin(ThetaW - 2*ThetaWp) + 7.745966692414834*g1*g2*Sin(2*(ThetaW - ThetaWp)) + 7.745966692414834*g1*g2*Sin(2*(ThetaW + ThetaWp)) + 12.649110640673518*g2*gN*Sin(ThetaW + 2*ThetaWp) - 6*Sqr(g1) + 3*Cos(2*(ThetaW - ThetaWp))*Sqr(g1) - 6*Cos(2*ThetaWp)*Sqr(g1) + 3*Cos(2*(ThetaW + ThetaWp))*Sqr(g1) + 2*Cos(2*ThetaW)*(3*Sqr(g1) - 5*Sqr(g2)) - 10*Sqr(g2) - 5*Cos(2*(ThetaW - ThetaWp))*Sqr(g2) - 10*Cos(2*ThetaWp)*Sqr(g2) - 5*Cos(2*(ThetaW + ThetaWp))*Sqr(g2) - 8*Sqr(gN) + 8*Cos(2*ThetaWp)*Sqr(gN))*(Conj(UHpp(gt1,0))*UHpp(gt2,0) + Conj(UHpp(gt1,1))*UHpp(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHpp, typename fields::conj<fields::SHpp>::type, fields::VZp, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -0.0125*(9.797958971132712*g1*gN*Cos(ThetaW - 2*ThetaWp) - 9.797958971132712*g1*gN*Cos(ThetaW + 2*ThetaWp) + 15.491933384829668*g1*g2*Sin(2*ThetaW) + 12.649110640673518*g2*gN*Sin(ThetaW - 2*ThetaWp) - 7.745966692414834*g1*g2*Sin(2*(ThetaW - ThetaWp)) - 7.745966692414834*g1*g2*Sin(2*(ThetaW + ThetaWp)) - 12.649110640673518*g2*gN*Sin(ThetaW + 2*ThetaWp) - 6*Sqr(g1) - 3*Cos(2*(ThetaW - ThetaWp))*Sqr(g1) + 6*Cos(2*ThetaWp)*Sqr(g1) - 3*Cos(2*(ThetaW + ThetaWp))*Sqr(g1) + 2*Cos(2*ThetaW)*(3*Sqr(g1) - 5*Sqr(g2)) - 10*Sqr(g2) + 5*Cos(2*(ThetaW - ThetaWp))*Sqr(g2) + 10*Cos(2*ThetaWp)*Sqr(g2) + 5*Cos(2*(ThetaW + ThetaWp))*Sqr(g2) - 8*Sqr(gN) - 8*Cos(2*ThetaWp)*Sqr(gN))*(Conj(UHpp(gt1,0))*UHpp(gt2,0) + Conj(UHpp(gt1,1))*UHpp(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHpp, typename fields::conj<fields::SHpp>::type, fields::VZ, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -0.0125*(25.298221281347036*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) - 15.491933384829668*g1*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) - Sin(2*ThetaWp)*(-9*Sqr(g1) + Cos(2*ThetaW)*(3*Sqr(g1) - 5*Sqr(g2)) - 5*Sqr(g2) + 8*Sqr(gN)) - 2*Sin(2*ThetaWp)*(3*Sqr(g1) - 5*Sqr(g2))*Sqr(Cos(ThetaW)) - 19.595917942265423*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + 19.595917942265423*g1*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp)))*(Conj(UHpp(gt1,0))*UHpp(gt2,0) + Conj(UHpp(gt1,1))*UHpp(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, fields::VP, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -0.05*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*((-10*g2*Cos(ThetaW)*Cos(ThetaWp) + 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp))*ZP(gt1,0)*ZP(gt2,0) + 2*(-5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*ZP(gt1,1)*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VP, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.005555555555555556*(-((0.7745966692414834*g1*Cos(ThetaW) - 3*g2*Sin(ThetaW))*(30*g2*Cos(ThetaW)*Cos(ThetaWp) + 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1))) - 3.0983866769659336*g1*Cos(ThetaW)*(7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SDX, typename fields::conj<fields::SDX>::type, fields::VP, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZDX = MODELPARAMETER(ZDX);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -0.008606629658238702*g1*Cos(ThetaW)*((15.491933384829668*g1*Cos(ThetaWp)*Sin(ThetaW) - 18.973665961010276*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt2,j1)) + (15.491933384829668*g1*Cos(ThetaWp)*Sin(ThetaW) + 28.460498941515414*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Se, typename fields::conj<fields::Se>::type, fields::VP, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.1*(-((0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*(-5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1))) - 0.7745966692414834*g1*Cos(ThetaW)*(15.491933384829668*g1*Cos(ThetaWp)*Sin(ThetaW) - 3.1622776601683795*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHIp, typename fields::conj<fields::SHIp>::type, fields::VP, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -0.05*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*((-10*g2*Cos(ThetaW)*Cos(ThetaWp) + 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp))*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt2,j1)) + 2*(-5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt2,2 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHpp, typename fields::conj<fields::SHpp>::type, fields::VP, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -0.02*(3.872983346207417*g1*Cos(ThetaW) + 5*g2*Sin(ThetaW))*(-5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*(Conj(UHpp(gt1,0))*UHpp(gt2,0) + Conj(UHpp(gt1,1))*UHpp(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VP, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.005555555555555556*(-((0.7745966692414834*g1*Cos(ThetaW) + 3*g2*Sin(ThetaW))*(-30*g2*Cos(ThetaW)*Cos(ThetaWp) + 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1))) - 3.0983866769659336*g1*(15.491933384829668*g1*Cos(ThetaWp)*Sin(2*ThetaW) + 9.486832980505138*gN*Cos(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1)));

   return {result};
}

QuadrupleVectorVertex VertexImpl<typename fields::conj<fields::VWm>::type, fields::VP, fields::VWm, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> part1 = Cos(ThetaW)*Cos(ThetaWp)*Sin(ThetaW)*Sqr(g2);

   const std::complex<double> part2 = -(Cos(ThetaWp)*Sin(2*ThetaW)*Sqr(g2));

   const std::complex<double> part3 = Cos(ThetaW)*Cos(ThetaWp)*Sin(ThetaW)*Sqr(g2);

   return {part1, part2, part3};
}

InverseMetricVertex VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, fields::VP, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.01*(3.872983346207417*g1*Cos(ThetaW) + 5*g2*Sin(ThetaW))*((9.486832980505138*gN*Cos(ThetaWp) + 2*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*ZP(gt1,0)*ZP(gt2,0) - 2*(3.1622776601683795*gN*Cos(ThetaWp) + (5*g2*Cos(ThetaW) - 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*ZP(gt1,1)*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VP, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.005555555555555556*((0.7745966692414834*g1*Cos(ThetaW) - 3*g2*Sin(ThetaW))*(9.486832980505138*gN*Cos(ThetaWp) + 2*(15*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1)) + 3.0983866769659336*g1*Cos(ThetaW)*(9.486832980505138*gN*Cos(ThetaWp) + 7.745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SDX, typename fields::conj<fields::SDX>::type, fields::VP, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZDX = MODELPARAMETER(ZDX);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.008606629658238702*g1*Cos(ThetaW)*((18.973665961010276*gN*Cos(ThetaWp) + 15.491933384829668*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt2,j1)) + (-28.460498941515414*gN*Cos(ThetaWp) + 15.491933384829668*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Se, typename fields::conj<fields::Se>::type, fields::VP, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.1*(-((0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*(3.1622776601683795*gN*Cos(ThetaWp) + (5*g2*Cos(ThetaW) - 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1))) + 0.7745966692414834*g1*Cos(ThetaW)*(3.1622776601683795*gN*Cos(ThetaWp) + 15.491933384829668*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHIp, typename fields::conj<fields::SHIp>::type, fields::VP, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.01*(3.872983346207417*g1*Cos(ThetaW) + 5*g2*Sin(ThetaW))*((9.486832980505138*gN*Cos(ThetaWp) + 2*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt2,j1)) - 2*(3.1622776601683795*gN*Cos(ThetaWp) + (5*g2*Cos(ThetaW) - 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt2,2 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHpp, typename fields::conj<fields::SHpp>::type, fields::VP, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.02*(3.872983346207417*g1*Cos(ThetaW) + 5*g2*Sin(ThetaW))*(-3.1622776601683795*gN*Cos(ThetaWp) + (-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*(Conj(UHpp(gt1,0))*UHpp(gt2,0) + Conj(UHpp(gt1,1))*UHpp(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VP, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.005555555555555556*((0.7745966692414834*g1*Cos(ThetaW) + 3*g2*Sin(ThetaW))*(9.486832980505138*gN*Cos(ThetaWp) + 2*(-15*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1)) + 3.0983866769659336*g1*(-9.486832980505138*gN*Cos(ThetaW)*Cos(ThetaWp) + 15.491933384829668*g1*Sin(2*ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1)));

   return {result};
}

QuadrupleVectorVertex VertexImpl<typename fields::conj<fields::VWm>::type, fields::VP, fields::VWm, fields::VZp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> part1 = -(Cos(ThetaW)*Sin(ThetaW)*Sin(ThetaWp)*Sqr(g2));

   const std::complex<double> part2 = Sin(2*ThetaW)*Sin(ThetaWp)*Sqr(g2);

   const std::complex<double> part3 = -(Cos(ThetaW)*Sin(ThetaW)*Sin(ThetaWp)*Sqr(g2));

   return {part1, part2, part3};
}

InverseMetricVertex VertexImpl<fields::Ah, typename fields::conj<fields::Hpm>::type, fields::VP, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = std::complex<double>(0.,-0.3872983346207417)*g1*g2*Cos(ThetaW)*(ZA(gt1,0)*ZP(gt2,0) + ZA(gt1,1)*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Ah, fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = std::complex<double>(0.,0.3872983346207417)*g1*g2*Cos(ThetaW)*(ZA(gt1,0)*ZP(gt2,0) + ZA(gt1,1)*ZP(gt2,1));

   return {result};
}

ChiralVertex VertexImpl<typename fields::bar<fields::ChaP>::type, fields::ChaP, fields::VP>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.5*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   const std::complex<double> right = 0.5*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::ChaP>::type, fields::ChaP, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = 0.1*(5*g2*Cos(ThetaW)*Cos(ThetaWp) - 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) - 3.1622776601683795*gN*Sin(ThetaWp));

   const std::complex<double> right = 0.1*(5*g2*Cos(ThetaW)*Cos(ThetaWp) - 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) - 3.1622776601683795*gN*Sin(ThetaWp));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::ChaP>::type, fields::ChaP, fields::VZp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.1*(-3.1622776601683795*gN*Cos(ThetaWp) - 5*g2*Cos(ThetaW)*Sin(ThetaWp) + 3.872983346207417*g1*Sin(ThetaW)*Sin(ThetaWp));

   const std::complex<double> right = 0.1*(-3.1622776601683795*gN*Cos(ThetaWp) - 5*g2*Cos(ThetaW)*Sin(ThetaWp) + 3.872983346207417*g1*Sin(ThetaW)*Sin(ThetaWp));

   return {left, right};
}

QuadrupleVectorVertex VertexImpl<typename fields::conj<fields::VWm>::type, fields::VWm, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> part1 = -2*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Cos(ThetaWp));

   const std::complex<double> part2 = Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Cos(ThetaWp));

   const std::complex<double> part3 = Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Cos(ThetaWp));

   return {part1, part2, part3};
}

QuadrupleVectorVertex VertexImpl<typename fields::conj<fields::VWm>::type, fields::VWm, fields::VZ, fields::VZp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> part1 = Sin(2*ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW));

   const std::complex<double> part2 = -(Cos(ThetaWp)*Sin(ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)));

   const std::complex<double> part3 = -(Cos(ThetaWp)*Sin(ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)));

   return {part1, part2, part3};
}

InverseMetricVertex VertexImpl<fields::Ah, typename fields::conj<fields::Hpm>::type, fields::VWm, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = std::complex<double>(0,0.05)*g2*((7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp))*ZA(gt1,0)*ZP(gt2,0) + 2*(3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*ZA(gt1,1)*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Ah, fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = std::complex<double>(0,-0.05)*g2*((7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp))*ZA(gt1,0)*ZP(gt2,0) + 2*(3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*ZA(gt1,1)*ZP(gt2,1));

   return {result};
}

ChiralVertex VertexImpl<fields::ChiP, fields::ChiP, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZNp = MODELPARAMETER(ZNp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = -0.1*(5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*(Conj(ZNp(gt2,0))*ZNp(gt1,0) - Conj(ZNp(gt2,1))*ZNp(gt1,1));

   const std::complex<double> right = 0.1*(5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*(Conj(ZNp(gt1,0))*ZNp(gt2,0) - Conj(ZNp(gt1,1))*ZNp(gt2,1));

   return {left, right};
}

ChiralVertex VertexImpl<fields::FSI, fields::FSI, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = -0.7905694150420949*gN*KroneckerDelta(gt1,gt2)*Sin(ThetaWp);

   const std::complex<double> right = 0.7905694150420949*gN*KroneckerDelta(gt1,gt2)*Sin(ThetaWp);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fv>::type, fields::Fv, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = -0.1*KroneckerDelta(gt1,gt2)*(5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp));

   const std::complex<double> right = 0;

   return {left, right};
}

ChiralVertex VertexImpl<fields::ChiP, fields::ChiP, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZNp = MODELPARAMETER(ZNp);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.1*(3.1622776601683795*gN*Cos(ThetaWp) - (5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*(Conj(ZNp(gt2,0))*ZNp(gt1,0) - Conj(ZNp(gt2,1))*ZNp(gt1,1));

   const std::complex<double> right = 0.1*(3.1622776601683795*gN*Cos(ThetaWp) - (5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*(Conj(ZNp(gt1,0))*ZNp(gt2,0) - Conj(ZNp(gt1,1))*ZNp(gt2,1));

   return {left, right};
}

ChiralVertex VertexImpl<fields::FSI, fields::FSI, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = -0.7905694150420949*gN*Cos(ThetaWp)*KroneckerDelta(gt1,gt2);

   const std::complex<double> right = 0.7905694150420949*gN*Cos(ThetaWp)*KroneckerDelta(gt1,gt2);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fv>::type, fields::Fv, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.1*KroneckerDelta(gt1,gt2)*(3.1622776601683795*gN*Cos(ThetaWp) - (5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp));

   const std::complex<double> right = 0;

   return {left, right};
}

QuadrupleVectorVertex VertexImpl<typename fields::conj<fields::VWm>::type, fields::VWm, fields::VZp, fields::VZp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> part1 = -2*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Sin(ThetaWp));

   const std::complex<double> part2 = Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Sin(ThetaWp));

   const std::complex<double> part3 = Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Sin(ThetaWp));

   return {part1, part2, part3};
}

InverseMetricVertex VertexImpl<fields::Ah, typename fields::conj<fields::Hpm>::type, fields::VWm, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = std::complex<double>(0,-0.05)*g2*((9.486832980505138*gN*Cos(ThetaWp) + 7.745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZA(gt1,0)*ZP(gt2,0) - 2*(3.1622776601683795*gN*Cos(ThetaWp) - 3.872983346207417*g1*Sin(ThetaW)*Sin(ThetaWp))*ZA(gt1,1)*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Ah, fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = std::complex<double>(0,0.05)*g2*((9.486832980505138*gN*Cos(ThetaWp) + 7.745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZA(gt1,0)*ZP(gt2,0) - 2*(3.1622776601683795*gN*Cos(ThetaWp) - 3.872983346207417*g1*Sin(ThetaW)*Sin(ThetaWp))*ZA(gt1,1)*ZP(gt2,1));

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

ChiralVertex VertexImpl<fields::ChiI, fields::ChaI, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZMI = MODELPARAMETER(ZMI);
   const auto ZNI = MODELPARAMETER(ZNI);
   const auto ZPI = MODELPARAMETER(ZPI);

   const std::complex<double> left = -0.7071067811865475*g2*SUM(j1,0,1,Conj(ZMI(gt2,j1))*ZNI(gt1,j1));

   const std::complex<double> right = 0.7071067811865475*g2*SUM(j1,0,1,Conj(ZNI(gt1,2 + j1))*ZPI(gt2,j1));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::ChaI>::type, fields::ChiI, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZNI = MODELPARAMETER(ZNI);
   const auto ZMI = MODELPARAMETER(ZMI);
   const auto ZPI = MODELPARAMETER(ZPI);

   const std::complex<double> left = -0.7071067811865475*g2*SUM(j1,0,1,Conj(ZNI(gt2,j1))*ZMI(gt1,j1));

   const std::complex<double> right = 0.7071067811865475*g2*SUM(j1,0,1,Conj(ZPI(gt1,j1))*ZNI(gt2,2 + j1));

   return {left, right};
}

ChiralVertex VertexImpl<fields::ChiP, fields::ChaP, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZNp = MODELPARAMETER(ZNp);

   const std::complex<double> left = -0.7071067811865475*g2*ZNp(gt1,0);

   const std::complex<double> right = 0.7071067811865475*g2*Conj(ZNp(gt1,1));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::ChaP>::type, fields::ChiP, fields::VWm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZNp = MODELPARAMETER(ZNp);

   const std::complex<double> left = -0.7071067811865475*g2*Conj(ZNp(gt2,0));

   const std::complex<double> right = 0.7071067811865475*g2*ZNp(gt2,1);

   return {left, right};
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

ScalarVertex VertexImpl<fields::Ah, fields::Hpm, fields::Su, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0.,0.35355339059327373)*(Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt3,j1))*ZD(gt4,j1))*(ZA(gt1,0)*ZP(gt2,0) - ZA(gt1,1)*ZP(gt2,1)) + 2*(-(SUM(j1,0,2,AbsSqr(Yd(j1,j1))*Conj(ZU(gt3,j1))*ZD(gt4,j1))*ZA(gt1,0)*ZP(gt2,0)) + Lambdax*SUM(j1,0,2,Conj(Yu(j1,j1))*Conj(ZU(gt3,3 + j1))*ZD(gt4,j1))*ZA(gt1,2)*ZP(gt2,0) + SUM(j1,0,2,AbsSqr(Yu(j1,j1))*Conj(ZU(gt3,j1))*ZD(gt4,j1))*ZA(gt1,1)*ZP(gt2,1) - Conj(Lambdax)*SUM(j1,0,2,Conj(ZU(gt3,j1))*Yd(j1,j1)*ZD(gt4,3 + j1))*ZA(gt1,2)*ZP(gt2,1) + SUM(j1,0,2,Conj(Yu(j1,j1))*Conj(ZU(gt3,3 + j1))*Yd(j1,j1)*ZD(gt4,3 + j1))*(-(ZA(gt1,1)*ZP(gt2,0)) + ZA(gt1,0)*ZP(gt2,1))));

   return {result};
}

ScalarVertex VertexImpl<fields::SHIp, fields::Su, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::SHI0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = -0.25*Sqr(g2)*(SUM(j1,0,3,Conj(UHIp(gt1,j1))*UHI0(gt4,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZD(gt3,j2)) + SUM(j1,0,2,Conj(ZU(gt2,j1))*ZD(gt3,j1))*SUM(j2,0,3,Conj(UHIp(gt1,j2))*UHI0(gt4,j2)));

   return {result};
}

ScalarVertex VertexImpl<fields::SHpp, fields::Su, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::SHp0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto UHp0 = MODELPARAMETER(UHp0);

   const std::complex<double> result = -0.5*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZD(gt3,j1))*(Conj(UHpp(gt1,0))*UHp0(gt4,0) + Conj(UHpp(gt1,1))*UHp0(gt4,1));

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

   const std::complex<double> result = 0.25*(-(Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZV(gt4,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZD(gt3,j2))) - Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZV(gt4,j2)) - 4*SUM(j1,0,2,Conj(ZU(gt2,j1))*Yd(j1,j1)*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(Ye(j2,j2))*Conj(ZE(gt1,3 + j2))*ZV(gt4,j2)));

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
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.07453559924999298*g2*(2.449489742783178*g1*Cos(ThetaWp)*Sin(ThetaW) - 3*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZD(gt2,j1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Su, typename fields::conj<fields::Sd>::type, fields::VWm, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.07453559924999298*g2*(3*gN*Cos(ThetaWp) + 2.449489742783178*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZD(gt2,j1));

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
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0.,-0.35355339059327373)*(Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt3,j1))*ZE(gt4,j1))*(-(ZA(gt1,0)*ZP(gt2,0)) + ZA(gt1,1)*ZP(gt2,1)) + 2*(SUM(j1,0,2,AbsSqr(Ye(j1,j1))*Conj(ZV(gt3,j1))*ZE(gt4,j1))*ZA(gt1,0)*ZP(gt2,0) + Conj(Lambdax)*SUM(j1,0,2,Conj(ZV(gt3,j1))*Ye(j1,j1)*ZE(gt4,3 + j1))*ZA(gt1,2)*ZP(gt2,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::SHIp, fields::Sv, typename fields::conj<fields::Se>::type, typename fields::conj<fields::SHI0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);

   const std::complex<double> result = -0.25*Sqr(g2)*(SUM(j1,0,3,Conj(UHIp(gt1,j1))*UHI0(gt4,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZE(gt3,j2)) + SUM(j1,0,2,Conj(ZV(gt2,j1))*ZE(gt3,j1))*SUM(j2,0,3,Conj(UHIp(gt1,j2))*UHI0(gt4,j2)));

   return {result};
}

ScalarVertex VertexImpl<fields::SHpp, fields::Sv, typename fields::conj<fields::Se>::type, typename fields::conj<fields::SHp0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto UHp0 = MODELPARAMETER(UHp0);

   const std::complex<double> result = -0.5*Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZE(gt3,j1))*(Conj(UHpp(gt1,0))*UHp0(gt4,0) + Conj(UHpp(gt1,1))*UHp0(gt4,1));

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

   const std::complex<double> result = 0.25*(-(Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZU(gt4,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZE(gt3,j2))) - Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZU(gt4,j2)) - 4*SUM(j1,0,2,Conj(ZV(gt2,j1))*Ye(j1,j1)*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(Yd(j2,j2))*Conj(ZD(gt1,3 + j2))*ZU(gt4,j2)));

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
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.22360679774997896*g2*(2.449489742783178*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZE(gt2,j1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sv, typename fields::conj<fields::Se>::type, fields::VWm, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.22360679774997896*g2*(2*gN*Cos(ThetaWp) - 2.449489742783178*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZE(gt2,j1));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Hpm, fields::SHI0, typename fields::conj<fields::SHIp>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0.,0.35355339059327373)*Sqr(g2)*SUM(j1,0,3,Conj(UHI0(gt3,j1))*UHIp(gt4,j1))*(ZA(gt1,0)*ZP(gt2,0) - ZA(gt1,1)*ZP(gt2,1));

   return {result};
}

ScalarVertex VertexImpl<fields::SHI0, fields::SHpp, typename fields::conj<fields::SHIp>::type, typename fields::conj<fields::SHp0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto UHp0 = MODELPARAMETER(UHp0);

   const std::complex<double> result = -0.5*Sqr(g2)*SUM(j1,0,3,Conj(UHI0(gt1,j1))*UHIp(gt3,j1))*(Conj(UHpp(gt2,0))*UHp0(gt4,0) + Conj(UHpp(gt2,1))*UHp0(gt4,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Sd, fields::SHI0, typename fields::conj<fields::SHIp>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);

   const std::complex<double> result = -0.25*Sqr(g2)*(SUM(j1,0,3,Conj(UHI0(gt2,j1))*UHIp(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZU(gt4,j2)) + SUM(j1,0,2,Conj(ZD(gt1,j1))*ZU(gt4,j1))*SUM(j2,0,3,Conj(UHI0(gt2,j2))*UHIp(gt3,j2)));

   return {result};
}

ScalarVertex VertexImpl<fields::Se, fields::SHI0, typename fields::conj<fields::SHIp>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);

   const std::complex<double> result = -0.25*Sqr(g2)*(SUM(j1,0,3,Conj(UHI0(gt2,j1))*UHIp(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZV(gt4,j2)) + SUM(j1,0,2,Conj(ZE(gt1,j1))*ZV(gt4,j1))*SUM(j2,0,3,Conj(UHI0(gt2,j2))*UHIp(gt3,j2)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHI0, typename fields::conj<fields::SHIp>::type, fields::VP, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.5477225575051661*g1*g2*Cos(ThetaW)*(SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHIp(gt2,j1)) - SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHIp(gt2,2 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHI0, typename fields::conj<fields::SHIp>::type, fields::VWm, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.22360679774997896*g2*((2.449489742783178*g1*Cos(ThetaWp)*Sin(ThetaW) - 3*gN*Sin(ThetaWp))*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHIp(gt2,j1)) - (2.449489742783178*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gN*Sin(ThetaWp))*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHIp(gt2,2 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHI0, typename fields::conj<fields::SHIp>::type, fields::VWm, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.22360679774997896*g2*((3*gN*Cos(ThetaWp) + 2.449489742783178*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHIp(gt2,j1)) + (2*gN*Cos(ThetaWp) - 2.449489742783178*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHIp(gt2,2 + j1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Hpm, fields::SHp0, typename fields::conj<fields::SHpp>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0.,0.35355339059327373)*Sqr(g2)*(Conj(UHp0(gt3,0))*UHpp(gt4,0) + Conj(UHp0(gt3,1))*UHpp(gt4,1))*(ZA(gt1,0)*ZP(gt2,0) - ZA(gt1,1)*ZP(gt2,1));

   return {result};
}

ScalarVertex VertexImpl<fields::SHIp, fields::SHp0, typename fields::conj<fields::SHI0>::type, typename fields::conj<fields::SHpp>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto UHpp = MODELPARAMETER(UHpp);

   const std::complex<double> result = -0.5*Sqr(g2)*SUM(j1,0,3,Conj(UHIp(gt1,j1))*UHI0(gt3,j1))*(Conj(UHp0(gt2,0))*UHpp(gt4,0) + Conj(UHp0(gt2,1))*UHpp(gt4,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Sd, fields::SHp0, typename fields::conj<fields::SHpp>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto UHpp = MODELPARAMETER(UHpp);

   const std::complex<double> result = -0.5*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZU(gt4,j1))*(Conj(UHp0(gt2,0))*UHpp(gt3,0) + Conj(UHp0(gt2,1))*UHpp(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Se, fields::SHp0, typename fields::conj<fields::SHpp>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto UHpp = MODELPARAMETER(UHpp);

   const std::complex<double> result = -0.5*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZV(gt4,j1))*(Conj(UHp0(gt2,0))*UHpp(gt3,0) + Conj(UHp0(gt2,1))*UHpp(gt3,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHp0, typename fields::conj<fields::SHpp>::type, fields::VP, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.5477225575051661*g1*g2*Cos(ThetaW)*(Conj(UHp0(gt1,0))*UHpp(gt2,0) - Conj(UHp0(gt1,1))*UHpp(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHp0, typename fields::conj<fields::SHpp>::type, fields::VWm, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.22360679774997896*g2*(2.449489742783178*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gN*Sin(ThetaWp))*(Conj(UHp0(gt1,0))*UHpp(gt2,0) - Conj(UHp0(gt1,1))*UHpp(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHp0, typename fields::conj<fields::SHpp>::type, fields::VWm, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.22360679774997896*g2*(2*gN*Cos(ThetaWp) - 2.449489742783178*g1*Sin(ThetaW)*Sin(ThetaWp))*(Conj(UHp0(gt1,0))*UHpp(gt2,0) - Conj(UHp0(gt1,1))*UHpp(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHIp, typename fields::conj<fields::SHI0>::type, typename fields::conj<fields::VWm>::type, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.5477225575051661*g1*g2*Cos(ThetaW)*(SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHI0(gt2,j1)) - SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHI0(gt2,2 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHpp, typename fields::conj<fields::SHp0>::type, typename fields::conj<fields::VWm>::type, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.5477225575051661*g1*g2*Cos(ThetaW)*(Conj(UHpp(gt1,0))*UHp0(gt2,0) - Conj(UHpp(gt1,1))*UHp0(gt2,1));

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

InverseMetricVertex VertexImpl<fields::SHIp, typename fields::conj<fields::SHI0>::type, typename fields::conj<fields::VWm>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.22360679774997896*g2*((2.449489742783178*g1*Cos(ThetaWp)*Sin(ThetaW) - 3*gN*Sin(ThetaWp))*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHI0(gt2,j1)) - (2.449489742783178*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gN*Sin(ThetaWp))*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHI0(gt2,2 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHpp, typename fields::conj<fields::SHp0>::type, typename fields::conj<fields::VWm>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.22360679774997896*g2*(2.449489742783178*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gN*Sin(ThetaWp))*(Conj(UHpp(gt1,0))*UHp0(gt2,0) - Conj(UHpp(gt1,1))*UHp0(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sd, typename fields::conj<fields::Su>::type, typename fields::conj<fields::VWm>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.07453559924999298*g2*(2.449489742783178*g1*Cos(ThetaWp)*Sin(ThetaW) - 3*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZU(gt2,j1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Se, typename fields::conj<fields::Sv>::type, typename fields::conj<fields::VWm>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.22360679774997896*g2*(2.449489742783178*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZV(gt2,j1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHIp, typename fields::conj<fields::SHI0>::type, typename fields::conj<fields::VWm>::type, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.22360679774997896*g2*((3*gN*Cos(ThetaWp) + 2.449489742783178*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHI0(gt2,j1)) + (2*gN*Cos(ThetaWp) - 2.449489742783178*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHI0(gt2,2 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::SHpp, typename fields::conj<fields::SHp0>::type, typename fields::conj<fields::VWm>::type, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.22360679774997896*g2*(2*gN*Cos(ThetaWp) - 2.449489742783178*g1*Sin(ThetaW)*Sin(ThetaWp))*(Conj(UHpp(gt1,0))*UHp0(gt2,0) - Conj(UHpp(gt1,1))*UHp0(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sd, typename fields::conj<fields::Su>::type, typename fields::conj<fields::VWm>::type, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.07453559924999298*g2*(3*gN*Cos(ThetaWp) + 2.449489742783178*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZU(gt2,j1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Se, typename fields::conj<fields::Sv>::type, typename fields::conj<fields::VWm>::type, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.22360679774997896*g2*(2*gN*Cos(ThetaWp) - 2.449489742783178*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZV(gt2,j1));

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
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0.,-0.35355339059327373)*(Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZU(gt4,j1))*(ZA(gt1,0)*ZP(gt3,0) - ZA(gt1,1)*ZP(gt3,1)) + 2*(-(SUM(j1,0,2,AbsSqr(Yd(j1,j1))*Conj(ZD(gt2,j1))*ZU(gt4,j1))*ZA(gt1,0)*ZP(gt3,0)) + Conj(Lambdax)*SUM(j1,0,2,Conj(ZD(gt2,j1))*Yu(j1,j1)*ZU(gt4,3 + j1))*ZA(gt1,2)*ZP(gt3,0) + SUM(j1,0,2,AbsSqr(Yu(j1,j1))*Conj(ZD(gt2,j1))*ZU(gt4,j1))*ZA(gt1,1)*ZP(gt3,1) - Lambdax*SUM(j1,0,2,Conj(Yd(j1,j1))*Conj(ZD(gt2,3 + j1))*ZU(gt4,j1))*ZA(gt1,2)*ZP(gt3,1) + SUM(j1,0,2,Conj(Yd(j1,j1))*Conj(ZD(gt2,3 + j1))*Yu(j1,j1)*ZU(gt4,3 + j1))*(-(ZA(gt1,1)*ZP(gt3,0)) + ZA(gt1,0)*ZP(gt3,1))));

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
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0.,-0.35355339059327373)*(Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZV(gt4,j1))*(ZA(gt1,0)*ZP(gt3,0) - ZA(gt1,1)*ZP(gt3,1)) - 2*(SUM(j1,0,2,AbsSqr(Ye(j1,j1))*Conj(ZE(gt2,j1))*ZV(gt4,j1))*ZA(gt1,0)*ZP(gt3,0) + Lambdax*SUM(j1,0,2,Conj(Ye(j1,j1))*Conj(ZE(gt2,3 + j1))*ZV(gt4,j1))*ZA(gt1,2)*ZP(gt3,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::SHIp, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::SHI0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0.,-0.35355339059327373)*Sqr(g2)*SUM(j1,0,3,Conj(UHIp(gt2,j1))*UHI0(gt4,j1))*(ZA(gt1,0)*ZP(gt3,0) - ZA(gt1,1)*ZP(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::SHpp, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::SHp0>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0.,-0.35355339059327373)*Sqr(g2)*(Conj(UHpp(gt2,0))*UHp0(gt4,0) + Conj(UHpp(gt2,1))*UHp0(gt4,1))*(ZA(gt1,0)*ZP(gt3,0) - ZA(gt1,1)*ZP(gt3,1));

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
   const auto gN = MODELPARAMETER(gN);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = std::complex<double>(0,-0.05)*((10*g2*Cos(ThetaW)*Cos(ThetaWp) + 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp))*ZA(gt1,0)*ZH(gt2,0) - 2*(5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*ZA(gt1,1)*ZH(gt2,1) + 15.811388300841898*gN*Sin(ThetaWp)*ZA(gt1,2)*ZH(gt2,2));

   return {result, minuend_index, subtrahend_index};
}

ChiralVertex VertexImpl<fields::VZ, fields::Cha, typename fields::bar<fields::Cha>::type>::evaluate(
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

ChiralVertex VertexImpl<fields::VZ, fields::ChaI, typename fields::bar<fields::ChaI>::type>::evaluate(
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

ChiralVertex VertexImpl<fields::VZ, fields::ChaP, typename fields::bar<fields::ChaP>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = 0.1*(-5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp));

   const std::complex<double> right = 0.1*(-5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp));

   return {left, right};
}

ChiralVertex VertexImpl<fields::VZ, fields::Chi, fields::Chi>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = 0.05*(Conj(ZN(gt2,2))*(-10*g2*Cos(ThetaW)*Cos(ThetaWp) - 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp))*ZN(gt1,2) + 2*Conj(ZN(gt2,3))*(5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*ZN(gt1,3) - 15.811388300841898*gN*Conj(ZN(gt2,4))*Sin(ThetaWp)*ZN(gt1,4));

   const std::complex<double> right = 0.05*(Conj(ZN(gt1,2))*(10*g2*Cos(ThetaW)*Cos(ThetaWp) + 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp))*ZN(gt2,2) - 2*Conj(ZN(gt1,3))*(5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*ZN(gt2,3) + 15.811388300841898*gN*Conj(ZN(gt1,4))*Sin(ThetaWp)*ZN(gt2,4));

   return {left, right};
}

ChiralVertex VertexImpl<fields::VZ, fields::ChiI, fields::ChiI>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZNI = MODELPARAMETER(ZNI);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = 0.05*((-10*g2*Cos(ThetaW)*Cos(ThetaWp) - 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp))*SUM(j1,0,1,Conj(ZNI(gt2,j1))*ZNI(gt1,j1)) + 2*(5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*SUM(j1,0,1,Conj(ZNI(gt2,2 + j1))*ZNI(gt1,2 + j1)));

   const std::complex<double> right = 0.05*((10*g2*Cos(ThetaW)*Cos(ThetaWp) + 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp))*SUM(j1,0,1,Conj(ZNI(gt1,j1))*ZNI(gt2,j1)) - 2*(5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*SUM(j1,0,1,Conj(ZNI(gt1,2 + j1))*ZNI(gt2,2 + j1)));

   return {left, right};
}

ChiralVertex VertexImpl<fields::VZ, fields::ChiP, fields::ChiP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZNp = MODELPARAMETER(ZNp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = -0.1*(5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*(Conj(ZNp(gt2,0))*ZNp(gt1,0) - Conj(ZNp(gt2,1))*ZNp(gt1,1));

   const std::complex<double> right = 0.1*(5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*(Conj(ZNp(gt1,0))*ZNp(gt2,0) - Conj(ZNp(gt1,1))*ZNp(gt2,1));

   return {left, right};
}

ChiralVertex VertexImpl<fields::VZ, fields::FDX, typename fields::bar<fields::FDX>::type>::evaluate(
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

ChiralVertex VertexImpl<fields::VZ, fields::FSI, fields::FSI>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = -0.7905694150420949*gN*KroneckerDelta(gt1,gt2)*Sin(ThetaWp);

   const std::complex<double> right = 0.7905694150420949*gN*KroneckerDelta(gt1,gt2)*Sin(ThetaWp);

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
   const auto gN = MODELPARAMETER(gN);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = std::complex<double>(0,-0.05)*((10*g2*Cos(ThetaW)*Cos(ThetaWp) + 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp))*ZA(gt1,0)*ZH(gt2,0) - 2*(5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*ZA(gt1,1)*ZH(gt2,1) + 15.811388300841898*gN*Sin(ThetaWp)*ZA(gt1,2)*ZH(gt2,2));

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::VZ, fields::hh, fields::VZp>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto vd = MODELPARAMETER(vd);
   const auto g2 = MODELPARAMETER(g2);
   const auto gN = MODELPARAMETER(gN);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto vs = MODELPARAMETER(vs);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05*(-(vd*(9.486832980505138*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) - 9*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN) + 5*Sin(2*ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 7.348469228349534*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + g1*(3.872983346207417*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) + 3*g1*Sin(2*ThetaWp)*Sqr(Sin(ThetaW)) - 7.348469228349534*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*ZH(gt1,0)) + vu*(6.324555320336759*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) + 4*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN) - 5*Sin(2*ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 4.898979485566356*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) - g1*(3.872983346207417*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) + 3*g1*Sin(2*ThetaWp)*Sqr(Sin(ThetaW)) + 4.898979485566356*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*ZH(gt1,1) + 25*vs*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN)*ZH(gt1,2));

   return {result};
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
   const auto gN = MODELPARAMETER(gN);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05*((10*g2*Cos(ThetaW)*Cos(ThetaWp) - 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp))*ZP(gt1,0)*ZP(gt2,0) + 2*(5*g2*Cos(ThetaW)*Cos(ThetaWp) - 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) - 3.1622776601683795*gN*Sin(ThetaWp))*ZP(gt1,1)*ZP(gt2,1));

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
   const auto gN = MODELPARAMETER(gN);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.016666666666666666*((30*g2*Cos(ThetaW)*Cos(ThetaWp) + 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1)) + 2*(-7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, fields::SDX, typename fields::conj<fields::SDX>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZDX = MODELPARAMETER(ZDX);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.016666666666666666*(-2*(7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt2,j1)) - (15.491933384829668*g1*Cos(ThetaWp)*Sin(ThetaW) + 28.460498941515414*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt2,3 + j1)));

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
   const auto gN = MODELPARAMETER(gN);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05*(2*(5*g2*Cos(ThetaW)*Cos(ThetaWp) - 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) - 3.1622776601683795*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) + (-15.491933384829668*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1)));

   return {result, minuend_index, subtrahend_index};
}

} // namespace detail
} // namespace CE6SSM_cxx_diagrams
} // namespace flexiblesusy
