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

// File generated at Wed 16 Oct 2019 22:02:29

/**
 * @file cxx_qft/E6SSM_vertices.cpp
 *
 * This file was generated at Wed 16 Oct 2019 22:02:29 with FlexibleSUSY
 * 2.4.1 and SARAH 4.14.3 .
 */

#include "E6SSM_context_base.hpp"
#include "E6SSM_vertices.hpp"

#include "concatenate.hpp"
#include "wrappers.hpp"

#define INPUTPARAMETER(p) context.model.get_input().p
#define MODELPARAMETER(p) context.model.get_##p()
#define DERIVEDPARAMETER(p) context.model.p()
#define PHASE(p) context.model.get_##p()

namespace flexiblesusy {
namespace E6SSM_cxx_diagrams {
namespace detail {

ChiralVertex VertexImpl<fields::Ah, typename bar<fields::Fe>::type, fields::Fe>::evaluate(
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

ChiralVertex VertexImpl<fields::Chi, typename conj<fields::Se>::type, fields::Fe>::evaluate(
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

   const std::complex<double> right = 0.1*(-10*SUM(j1,0,2,Conj(Ye(j1,j1))*ZE(gt3,j1)*ZER(gt2,j1))*ZN(gt1,2) - 2.23606797749979*SUM(j1,0,2,ZE(gt3,3 + j1)*ZER(gt2,j1))*(4.898979485566356*g1*ZN(gt1,0) + gN*ZN(gt1,5)));

   return {left, right};
}

ChiralVertex VertexImpl<fields::hh, typename bar<fields::Fe>::type, fields::Fe>::evaluate(
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

ChiralVertex VertexImpl<typename bar<fields::Cha>::type, fields::Cha, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UM = MODELPARAMETER(UM);
   const auto UP = MODELPARAMETER(UP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.5*(2*g2*Conj(UM(gt2,0))*Sin(ThetaW)*UM(gt1,0) + Conj(UM(gt2,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UM(gt1,1));

   const std::complex<double> right = 0.5*(2*g2*Conj(UP(gt1,0))*Sin(ThetaW)*UP(gt2,0) + Conj(UP(gt1,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UP(gt2,1));

   return {left, right};
}

ChiralVertex VertexImpl<typename bar<fields::Fe>::type, fields::Cha, fields::Sv>::evaluate(
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

ChiralVertex VertexImpl<typename bar<fields::Fe>::type, fields::Fe, fields::Ah>::evaluate(
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

ChiralVertex VertexImpl<typename bar<fields::Fe>::type, fields::Fe, fields::hh>::evaluate(
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

ChiralVertex VertexImpl<typename bar<fields::Fe>::type, fields::Fe, fields::VP>::evaluate(
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

ChiralVertex VertexImpl<typename bar<fields::Fe>::type, fields::Hpm, fields::Fv>::evaluate(
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

ChiralVertex VertexImpl<typename bar<fields::Fe>::type, fields::Se, fields::Chi>::evaluate(
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

   const std::complex<double> right = 0.1*(-10*SUM(j1,0,2,Conj(Ye(j1,j1))*Conj(ZE(gt3,3 + j1))*ZEL(gt1,j1))*ZN(gt2,2) + SUM(j1,0,2,Conj(ZE(gt3,j1))*ZEL(gt1,j1))*(5.477225575051661*g1*ZN(gt2,0) + 7.0710678118654755*g2*ZN(gt2,1) - 4.47213595499958*gN*ZN(gt2,5)));

   return {left, right};
}

ChiralVertex VertexImpl<typename bar<fields::Fv>::type, typename conj<fields::Hpm>::type, fields::Fe>::evaluate(
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

MomentumDifferenceVertex VertexImpl<typename conj<fields::Hpm>::type, fields::Hpm, fields::VP>::evaluate(
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

MomentumDifferenceVertex VertexImpl<typename conj<fields::Se>::type, fields::Se, fields::VP>::evaluate(
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

   const std::complex<double> result = 0.5*((0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) + 1.5491933384829668*g1*Cos(ThetaW)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1)));

   return {result, minuend_index, subtrahend_index};
}

ChiralVertex VertexImpl<typename conj<fields::Sv>::type, typename bar<fields::Cha>::type, fields::Fe>::evaluate(
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

} // namespace detail
} // namespace E6SSM_cxx_diagrams
} // namespace flexiblesusy
