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
 * @file cxx_qft/MSSMRHN_vertices.cpp
 *
 * This file was generated with FlexibleSUSY 2.6.1 and SARAH 4.14.5 .
 */

#include "MSSMRHN_context_base.hpp"
#include "MSSMRHN_input_parameters.hpp"
#include "MSSMRHN_vertices.hpp"

#include "concatenate.hpp"
#include "wrappers.hpp"

#define INPUTPARAMETER(p) context.model.get_input_parameters().p
#define MODELPARAMETER(p) context.model.get_##p()
#define DERIVEDPARAMETER(p) context.model.p()
#define PHASE(p) context.model.get_##p()

namespace flexiblesusy {
namespace MSSMRHN_cxx_diagrams {
namespace detail {

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Ah, fields::Ah>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = 0.25*(0.6*Sqr(g1) + Sqr(g2))*(ZA(gt1,1)*(ZA(gt2,0)*(ZA(gt3,1)*ZA(gt4,0) + ZA(gt3,0)*ZA(gt4,1)) + ZA(gt2,1)*(ZA(gt3,0)*ZA(gt4,0) - 3*ZA(gt3,1)*ZA(gt4,1))) + ZA(gt1,0)*(ZA(gt2,1)*(ZA(gt3,1)*ZA(gt4,0) + ZA(gt3,0)*ZA(gt4,1)) + ZA(gt2,0)*(-3*ZA(gt3,0)*ZA(gt4,0) + ZA(gt3,1)*ZA(gt4,1))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::hh, fields::hh>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = -0.25*(0.6*Sqr(g1) + Sqr(g2))*(ZA(gt1,0)*ZA(gt2,0) - ZA(gt1,1)*ZA(gt2,1))*(ZH(gt3,0)*ZH(gt4,0) - ZH(gt3,1)*ZH(gt4,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = -0.25*(0.6*Sqr(g1) + Sqr(g2))*(ZA(gt1,0)*ZA(gt2,0) - ZA(gt1,1)*ZA(gt2,1))*(vd*ZH(gt3,0) - vu*ZH(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Hpm, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(ZA(gt1,0)*(Sqr(g2)*ZA(gt2,1)*(ZP(gt3,1)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,1)) - ZA(gt2,0)*((0.6*Sqr(g1) + Sqr(g2))*ZP(gt3,0)*ZP(gt4,0) + (-0.6*Sqr(g1) + Sqr(g2))*ZP(gt3,1)*ZP(gt4,1))) + ZA(gt1,1)*(Sqr(g2)*ZA(gt2,0)*(ZP(gt3,1)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,1)) + ZA(gt2,1)*((0.6*Sqr(g1) - Sqr(g2))*ZP(gt3,0)*ZP(gt4,0) - (0.6*Sqr(g1) + Sqr(g2))*ZP(gt3,1)*ZP(gt4,1))));

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
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = 0.05*(-20*(SUM(j3,0,2,Conj(ZD(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,j1))*ZD(gt4,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt3,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gt4,j3)))*ZA(gt1,0)*ZA(gt2,0) + (Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZD(gt3,j1))*ZD(gt4,j1))*(ZA(gt1,0)*ZA(gt2,0) - ZA(gt1,1)*ZA(gt2,1)) + 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt3,3 + j1))*ZD(gt4,3 + j1))*(ZA(gt1,0)*ZA(gt2,0) - ZA(gt1,1)*ZA(gt2,1)));

   return {result};
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

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Sv, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = -((SUM(j3,0,2,Conj(ZV(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j3,j1))*Yv(j2,j1))*ZV(gt4,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gt3,j2))*SUM(j1,0,2,Conj(Yv(j1,j3))*Yv(j1,j2)))*ZV(gt4,j3)))*ZA(gt1,1)*ZA(gt2,1)) - 0.25*(0.6*Sqr(g1) + Sqr(g2))*SUM(j1,0,2,Conj(ZV(gt3,j1))*ZV(gt4,j1))*(ZA(gt1,0)*ZA(gt2,0) - ZA(gt1,1)*ZA(gt2,1));

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

ScalarVertex VertexImpl<fields::Ah, fields::hh, fields::Hpm, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,-0.25)*Sqr(g2)*(ZA(gt1,1)*ZH(gt2,0) + ZA(gt1,0)*ZH(gt2,1))*(ZP(gt3,1)*ZP(gt4,0) - ZP(gt3,0)*ZP(gt4,1));

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

   const std::complex<double> result = std::complex<double>(0.,-0.35355339059327373)*(2*SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gt4,j3))*ZA(gt1,0)*ZP(gt2,0) - 2*SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZD(gt4,j3))*ZA(gt1,1)*ZP(gt2,1) + 2*SUM(j3,0,2,Conj(ZU(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yd(j2,j1))*ZD(gt4,3 + j2)))*(ZA(gt1,1)*ZP(gt2,0) - ZA(gt1,0)*ZP(gt2,1)) + Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt3,j1))*ZD(gt4,j1))*(-(ZA(gt1,0)*ZP(gt2,0)) + ZA(gt1,1)*ZP(gt2,1)));

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
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0.,-0.35355339059327373)*(2*SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gt3,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZE(gt4,j3))*ZA(gt1,0)*ZP(gt2,0) - 2*SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gt3,j2))*SUM(j1,0,2,Conj(Yv(j1,j3))*Yv(j1,j2)))*ZE(gt4,j3))*ZA(gt1,1)*ZP(gt2,1) + 2*SUM(j3,0,2,Conj(ZV(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j3,j1))*Ye(j2,j1))*ZE(gt4,3 + j2)))*(ZA(gt1,1)*ZP(gt2,0) - ZA(gt1,0)*ZP(gt2,1)) + Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt3,j1))*ZE(gt4,j1))*(-(ZA(gt1,0)*ZP(gt2,0)) + ZA(gt1,1)*ZP(gt2,1)));

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
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,-0.25)*Sqr(g2)*(vu*ZA(gt1,0) + vd*ZA(gt1,1))*(ZP(gt2,1)*ZP(gt3,0) - ZP(gt2,0)*ZP(gt3,1));

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

InverseMetricVertex VertexImpl<fields::Ah, fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = std::complex<double>(0.,-0.3872983346207417)*g1*g2*Sin(ThetaW)*(ZA(gt1,0)*ZP(gt2,0) + ZA(gt1,1)*ZP(gt2,1));

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

   const std::complex<double> result = std::complex<double>(0,-0.5)*g2*(ZA(gt1,0)*ZP(gt2,0) + ZA(gt1,1)*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
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

   const std::complex<double> result = std::complex<double>(0.,0.35355339059327373)*(2*SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZU(gt4,j3))*ZA(gt1,0)*ZP(gt3,0) - 2*SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gt4,j3))*ZA(gt1,1)*ZP(gt3,1) + 2*SUM(j3,0,2,Conj(ZD(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yu(j2,j1))*ZU(gt4,3 + j2)))*(ZA(gt1,1)*ZP(gt3,0) - ZA(gt1,0)*ZP(gt3,1)) + Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZU(gt4,j1))*(-(ZA(gt1,0)*ZP(gt3,0)) + ZA(gt1,1)*ZP(gt3,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Sd, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto TYd = MODELPARAMETER(TYd);
   const auto Mu = MODELPARAMETER(Mu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475)*(SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,ZD(gt3,3 + j1)*TYd(j1,j2)))*ZA(gt1,0) - SUM(j2,0,2,SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j2)))*ZD(gt3,j2))*ZA(gt1,0) + (Conj(Mu)*SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1))) - Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1)))*ZD(gt3,j2)))*ZA(gt1,1));

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
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0.,0.35355339059327373)*(2*SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZV(gt4,j3))*ZA(gt1,0)*ZP(gt3,0) - 2*SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Conj(Yv(j1,j3))*Yv(j1,j2)))*ZV(gt4,j3))*ZA(gt1,1)*ZP(gt3,1) + 2*SUM(j3,0,2,Conj(ZE(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Yv(j2,j1))*ZV(gt4,3 + j2)))*(ZA(gt1,1)*ZP(gt3,0) - ZA(gt1,0)*ZP(gt3,1)) + Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZV(gt4,j1))*(-(ZA(gt1,0)*ZP(gt3,0)) + ZA(gt1,1)*ZP(gt3,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Se, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto TYe = MODELPARAMETER(TYe);
   const auto Mu = MODELPARAMETER(Mu);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475)*(SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,ZE(gt3,3 + j1)*TYe(j1,j2)))*ZA(gt1,0) - SUM(j2,0,2,SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j2)))*ZE(gt3,j2))*ZA(gt1,0) + (Conj(Mu)*SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1))) - Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1)))*ZE(gt3,j2)))*ZA(gt1,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Su, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto TYu = MODELPARAMETER(TYu);
   const auto Mu = MODELPARAMETER(Mu);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475)*(Conj(Mu)*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt3,3 + j1)))*ZA(gt1,0) - Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1)))*ZU(gt3,j2))*ZA(gt1,0) + (SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,ZU(gt3,3 + j1)*TYu(j1,j2))) - SUM(j2,0,2,SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j2)))*ZU(gt3,j2)))*ZA(gt1,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Sv, fields::Sv>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Mv = MODELPARAMETER(Mv);
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0.,-0.35355339059327373)*(SUM(j3,0,2,Conj(ZV(gt3,3 + j3))*SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Conj(Mv(j1,j3))*Yv(j1,j2)))) + SUM(j3,0,2,Conj(ZV(gt2,3 + j3))*SUM(j2,0,2,Conj(ZV(gt3,j2))*SUM(j1,0,2,Conj(Mv(j1,j3))*Yv(j1,j2)))) + SUM(j3,0,2,Conj(ZV(gt3,3 + j3))*SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Conj(Mv(j3,j1))*Yv(j1,j2)))) + SUM(j3,0,2,Conj(ZV(gt2,3 + j3))*SUM(j2,0,2,Conj(ZV(gt3,j2))*SUM(j1,0,2,Conj(Mv(j3,j1))*Yv(j1,j2)))))*ZA(gt1,1);

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Sv, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto TYv = MODELPARAMETER(TYv);
   const auto Mu = MODELPARAMETER(Mu);
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475)*(Conj(Mu)*SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Yv(j1,j2)*ZV(gt3,3 + j1)))*ZA(gt1,0) - Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j2))*Conj(ZV(gt2,3 + j1)))*ZV(gt3,j2))*ZA(gt1,0) + (SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,ZV(gt3,3 + j1)*TYv(j1,j2))) - SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(gt2,3 + j1))*Conj(TYv(j1,j2)))*ZV(gt3,j2)))*ZA(gt1,1));

   return {result};
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

InverseMetricVertex VertexImpl<fields::Ah, typename fields::conj<fields::Hpm>::type, fields::VWm, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = std::complex<double>(0.,0.3872983346207417)*g1*g2*Sin(ThetaW)*(ZA(gt1,0)*ZP(gt2,0) + ZA(gt1,1)*ZP(gt2,1));

   return {result};
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

ScalarVertex VertexImpl<fields::Ah, typename fields::conj<fields::Sv>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Yv = MODELPARAMETER(Yv);
   const auto Mv = MODELPARAMETER(Mv);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0.,0.35355339059327373)*(SUM(j3,0,2,SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j3))*Mv(j1,j2))*ZV(gt3,3 + j2))*ZV(gt2,j3)) + SUM(j3,0,2,SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j3))*Mv(j2,j1))*ZV(gt3,3 + j2))*ZV(gt2,j3)) + SUM(j3,0,2,SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j3))*Mv(j1,j2))*ZV(gt2,3 + j2))*ZV(gt3,j3)) + SUM(j3,0,2,SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j3))*Mv(j2,j1))*ZV(gt2,3 + j2))*ZV(gt3,j3)))*ZA(gt1,1);

   return {result};
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

   const std::complex<double> left = -(g2*Conj(UM(gt1,0))*SUM(j1,0,2,Conj(ZUL(gt2,j1))*ZD(gt3,j1))) + Conj(UM(gt1,1))*SUM(j2,0,2,Conj(ZUL(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1)));

   const std::complex<double> right = SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*ZUR(gt2,j1))*ZD(gt3,j2))*UP(gt1,1);

   return {left, right};
}

ChiralVertex VertexImpl<fields::Chi, fields::Cha, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto UM = MODELPARAMETER(UM);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZP = MODELPARAMETER(ZP);
   const auto UP = MODELPARAMETER(UP);

   const std::complex<double> left = 0.5*(Conj(UM(gt2,1))*(1.0954451150103321*g1*Conj(ZN(gt1,0)) + 1.4142135623730951*g2*Conj(ZN(gt1,1))) - 2*g2*Conj(UM(gt2,0))*Conj(ZN(gt1,2)))*ZP(gt3,0);

   const std::complex<double> right = -0.5*(UP(gt2,1)*(1.0954451150103321*g1*ZN(gt1,0) + 1.4142135623730951*g2*ZN(gt1,1)) + 2*g2*UP(gt2,0)*ZN(gt1,3))*ZP(gt3,1);

   return {left, right};
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

ChiralVertex VertexImpl<fields::Chi, fields::Chi, fields::Ah>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0,-0.5)*(Conj(ZN(gt1,2))*(0.7745966692414834*g1*Conj(ZN(gt2,0)) - g2*Conj(ZN(gt2,1)))*ZA(gt3,0) - g2*Conj(ZN(gt1,1))*Conj(ZN(gt2,2))*ZA(gt3,0) - 0.7745966692414834*g1*Conj(ZN(gt1,3))*Conj(ZN(gt2,0))*ZA(gt3,1) + g2*Conj(ZN(gt1,3))*Conj(ZN(gt2,1))*ZA(gt3,1) + g2*Conj(ZN(gt1,1))*Conj(ZN(gt2,3))*ZA(gt3,1) + 0.7745966692414834*g1*Conj(ZN(gt1,0))*(Conj(ZN(gt2,2))*ZA(gt3,0) - Conj(ZN(gt2,3))*ZA(gt3,1)));

   const std::complex<double> right = std::complex<double>(0,0.1)*(ZA(gt3,0)*(ZN(gt1,2)*(3.872983346207417*g1*ZN(gt2,0) - 5*g2*ZN(gt2,1)) + (3.872983346207417*g1*ZN(gt1,0) - 5*g2*ZN(gt1,1))*ZN(gt2,2)) - ZA(gt3,1)*(ZN(gt1,3)*(3.872983346207417*g1*ZN(gt2,0) - 5*g2*ZN(gt2,1)) + (3.872983346207417*g1*ZN(gt1,0) - 5*g2*ZN(gt1,1))*ZN(gt2,3)));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Chi, fields::Chi, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = 0.5*(Conj(ZN(gt1,2))*(0.7745966692414834*g1*Conj(ZN(gt2,0)) - g2*Conj(ZN(gt2,1)))*ZH(gt3,0) - g2*Conj(ZN(gt1,1))*Conj(ZN(gt2,2))*ZH(gt3,0) - 0.7745966692414834*g1*Conj(ZN(gt1,3))*Conj(ZN(gt2,0))*ZH(gt3,1) + g2*Conj(ZN(gt1,3))*Conj(ZN(gt2,1))*ZH(gt3,1) + g2*Conj(ZN(gt1,1))*Conj(ZN(gt2,3))*ZH(gt3,1) + 0.7745966692414834*g1*Conj(ZN(gt1,0))*(Conj(ZN(gt2,2))*ZH(gt3,0) - Conj(ZN(gt2,3))*ZH(gt3,1)));

   const std::complex<double> right = 0.1*(ZH(gt3,0)*(ZN(gt1,2)*(3.872983346207417*g1*ZN(gt2,0) - 5*g2*ZN(gt2,1)) + (3.872983346207417*g1*ZN(gt1,0) - 5*g2*ZN(gt1,1))*ZN(gt2,2)) - ZH(gt3,1)*(ZN(gt1,3)*(3.872983346207417*g1*ZN(gt2,0) - 5*g2*ZN(gt2,1)) + (3.872983346207417*g1*ZN(gt1,0) - 5*g2*ZN(gt1,1))*ZN(gt2,3)));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Chi, fields::Chi, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.5*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(Conj(ZN(gt2,2))*ZN(gt1,2) - Conj(ZN(gt2,3))*ZN(gt1,3));

   const std::complex<double> right = 0.5*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(Conj(ZN(gt1,2))*ZN(gt2,2) - Conj(ZN(gt1,3))*ZN(gt2,3));

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
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZDR = MODELPARAMETER(ZDR);

   const std::complex<double> left = -0.18257418583505536*g1*Conj(ZN(gt1,0))*SUM(j1,0,2,Conj(ZDL(gt2,j1))*ZD(gt3,j1)) + 0.7071067811865475*g2*Conj(ZN(gt1,1))*SUM(j1,0,2,Conj(ZDL(gt2,j1))*ZD(gt3,j1)) - Conj(ZN(gt1,2))*SUM(j2,0,2,Conj(ZDL(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1)));

   const std::complex<double> right = -0.3651483716701107*g1*SUM(j1,0,2,ZD(gt3,3 + j1)*ZDR(gt2,j1))*ZN(gt1,0) - SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*ZDR(gt2,j1))*ZD(gt3,j2))*ZN(gt1,2);

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
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZER = MODELPARAMETER(ZER);

   const std::complex<double> left = 0.5477225575051661*g1*Conj(ZN(gt1,0))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZE(gt3,j1)) + 0.7071067811865475*g2*Conj(ZN(gt1,1))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZE(gt3,j1)) - Conj(ZN(gt1,2))*SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)));

   const std::complex<double> right = -1.0954451150103321*g1*SUM(j1,0,2,ZE(gt3,3 + j1)*ZER(gt2,j1))*ZN(gt1,0) - SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZE(gt3,j2))*ZN(gt1,2);

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
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZUR = MODELPARAMETER(ZUR);

   const std::complex<double> left = -0.18257418583505536*g1*Conj(ZN(gt1,0))*SUM(j1,0,2,Conj(ZUL(gt2,j1))*ZU(gt3,j1)) - 0.7071067811865475*g2*Conj(ZN(gt1,1))*SUM(j1,0,2,Conj(ZUL(gt2,j1))*ZU(gt3,j1)) - Conj(ZN(gt1,3))*SUM(j2,0,2,Conj(ZUL(gt2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt3,3 + j1)));

   const std::complex<double> right = 0.7302967433402214*g1*SUM(j1,0,2,ZU(gt3,3 + j1)*ZUR(gt2,j1))*ZN(gt1,0) - SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*ZUR(gt2,j1))*ZU(gt3,j2))*ZN(gt1,3);

   return {left, right};
}

ChiralVertex VertexImpl<fields::Chi, fields::Fv, fields::Sv>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZV = MODELPARAMETER(ZV);
   const auto UV = MODELPARAMETER(UV);

   const std::complex<double> left = -(Conj(ZN(gt1,3))*SUM(j2,0,2,Conj(ZV(gt3,j2))*SUM(j1,0,2,Conj(UV(gt2,3 + j1))*Yv(j1,j2))));

   const std::complex<double> right = 0.7071067811865475*SUM(j1,0,2,Conj(ZV(gt3,j1))*UV(gt2,j1))*(0.7745966692414834*g1*ZN(gt1,0) - g2*ZN(gt1,1)) - SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j2))*Conj(ZV(gt3,3 + j1)))*UV(gt2,j2))*ZN(gt1,3);

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
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZN = MODELPARAMETER(ZN);
   const auto UV = MODELPARAMETER(UV);
   const auto ZV = MODELPARAMETER(ZV);

   const std::complex<double> left = 0.5477225575051661*g1*Conj(ZN(gt1,0))*SUM(j1,0,2,Conj(UV(gt2,j1))*ZV(gt3,j1)) - 0.7071067811865475*g2*Conj(ZN(gt1,1))*SUM(j1,0,2,Conj(UV(gt2,j1))*ZV(gt3,j1)) - Conj(ZN(gt1,3))*SUM(j2,0,2,Conj(UV(gt2,j2))*SUM(j1,0,2,Yv(j1,j2)*ZV(gt3,3 + j1)));

   const std::complex<double> right = -(SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j2))*UV(gt2,3 + j1))*ZV(gt3,j2))*ZN(gt1,3));

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
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZER = MODELPARAMETER(ZER);

   const std::complex<double> left = 0.5477225575051661*g1*Conj(ZN(gt1,0))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZE(gt3,j1)) + 0.7071067811865475*g2*Conj(ZN(gt1,1))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZE(gt3,j1)) - Conj(ZN(gt1,2))*SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)));

   const std::complex<double> right = -1.0954451150103321*g1*SUM(j1,0,2,ZE(gt3,3 + j1)*ZER(gt2,j1))*ZN(gt1,0) - SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZE(gt3,j2))*ZN(gt1,2);

   return {left, right};
}

ChiralVertex VertexImpl<fields::Fv, fields::Cha, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Yv = MODELPARAMETER(Yv);
   const auto UM = MODELPARAMETER(UM);
   const auto UV = MODELPARAMETER(UV);
   const auto ZE = MODELPARAMETER(ZE);
   const auto UP = MODELPARAMETER(UP);

   const std::complex<double> left = -(g2*Conj(UM(gt2,0))*SUM(j1,0,2,Conj(UV(gt1,j1))*ZE(gt3,j1))) + Conj(UM(gt2,1))*SUM(j2,0,2,Conj(UV(gt1,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)));

   const std::complex<double> right = SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j2))*UV(gt1,3 + j1))*ZE(gt3,j2))*UP(gt2,1);

   return {left, right};
}

ChiralVertex VertexImpl<fields::Fv, fields::Fe, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Yv = MODELPARAMETER(Yv);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto UV = MODELPARAMETER(UV);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ZER = MODELPARAMETER(ZER);

   const std::complex<double> left = SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Conj(UV(gt1,3 + j1))*Yv(j1,j2)))*ZP(gt3,1);

   const std::complex<double> right = SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*UV(gt1,j2))*ZP(gt3,0);

   return {left, right};
}

ChiralVertex VertexImpl<fields::Fv, fields::Fe, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto UV = MODELPARAMETER(UV);

   const std::complex<double> left = -0.7071067811865475*g2*SUM(j1,0,2,Conj(ZEL(gt2,j1))*UV(gt1,j1));

   const std::complex<double> right = 0;

   return {left, right};
}

ChiralVertex VertexImpl<fields::Fv, fields::Fv, fields::Ah>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Yv = MODELPARAMETER(Yv);
   const auto UV = MODELPARAMETER(UV);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*(SUM(j2,0,2,Conj(UV(gt2,j2))*SUM(j1,0,2,Conj(UV(gt1,3 + j1))*Yv(j1,j2))) + SUM(j2,0,2,Conj(UV(gt1,j2))*SUM(j1,0,2,Conj(UV(gt2,3 + j1))*Yv(j1,j2))))*ZA(gt3,1);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*(SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j2))*UV(gt2,3 + j1))*UV(gt1,j2)) + SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j2))*UV(gt1,3 + j1))*UV(gt2,j2)))*ZA(gt3,1);

   return {left, right};
}

ChiralVertex VertexImpl<fields::Fv, fields::Fv, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Yv = MODELPARAMETER(Yv);
   const auto UV = MODELPARAMETER(UV);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*(SUM(j2,0,2,Conj(UV(gt2,j2))*SUM(j1,0,2,Conj(UV(gt1,3 + j1))*Yv(j1,j2))) + SUM(j2,0,2,Conj(UV(gt1,j2))*SUM(j1,0,2,Conj(UV(gt2,3 + j1))*Yv(j1,j2))))*ZH(gt3,1);

   const std::complex<double> right = -0.7071067811865475*(SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j2))*UV(gt2,3 + j1))*UV(gt1,j2)) + SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j2))*UV(gt1,3 + j1))*UV(gt2,j2)))*ZH(gt3,1);

   return {left, right};
}

ChiralVertex VertexImpl<fields::Fv, fields::Fv, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UV = MODELPARAMETER(UV);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.5*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*SUM(j1,0,2,Conj(UV(gt2,j1))*UV(gt1,j1));

   const std::complex<double> right = 0.5*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*SUM(j1,0,2,Conj(UV(gt1,j1))*UV(gt2,j1));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Fv, typename fields::conj<fields::Hpm>::type, fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Yv = MODELPARAMETER(Yv);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto UV = MODELPARAMETER(UV);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ZER = MODELPARAMETER(ZER);

   const std::complex<double> left = SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Conj(UV(gt1,3 + j1))*Yv(j1,j2)))*ZP(gt3,1);

   const std::complex<double> right = SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*UV(gt1,j2))*ZP(gt3,0);

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

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::hh, fields::hh>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.25*(0.6*Sqr(g1) + Sqr(g2))*(ZH(gt1,1)*(ZH(gt2,0)*(ZH(gt3,1)*ZH(gt4,0) + ZH(gt3,0)*ZH(gt4,1)) + ZH(gt2,1)*(ZH(gt3,0)*ZH(gt4,0) - 3*ZH(gt3,1)*ZH(gt4,1))) + ZH(gt1,0)*(ZH(gt2,1)*(ZH(gt3,1)*ZH(gt4,0) + ZH(gt3,0)*ZH(gt4,1)) + ZH(gt2,0)*(-3*ZH(gt3,0)*ZH(gt4,0) + ZH(gt3,1)*ZH(gt4,1))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vu = MODELPARAMETER(vu);
   const auto vd = MODELPARAMETER(vd);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.25*(0.6*Sqr(g1) + Sqr(g2))*(ZH(gt1,1)*(ZH(gt2,0)*(vu*ZH(gt3,0) + vd*ZH(gt3,1)) + ZH(gt2,1)*(vd*ZH(gt3,0) - 3*vu*ZH(gt3,1))) + ZH(gt1,0)*(ZH(gt2,1)*(vu*ZH(gt3,0) + vd*ZH(gt3,1)) + ZH(gt2,0)*(-3*vd*ZH(gt3,0) + vu*ZH(gt3,1))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::Hpm, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-(ZH(gt1,0)*(Sqr(g2)*ZH(gt2,1)*(ZP(gt3,1)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,1)) + ZH(gt2,0)*((0.6*Sqr(g1) + Sqr(g2))*ZP(gt3,0)*ZP(gt4,0) + (-0.6*Sqr(g1) + Sqr(g2))*ZP(gt3,1)*ZP(gt4,1)))) - ZH(gt1,1)*(Sqr(g2)*ZH(gt2,0)*(ZP(gt3,1)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,1)) + ZH(gt2,1)*((-0.6*Sqr(g1) + Sqr(g2))*ZP(gt3,0)*ZP(gt4,0) + (0.6*Sqr(g1) + Sqr(g2))*ZP(gt3,1)*ZP(gt4,1))));

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
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = -((SUM(j3,0,2,Conj(ZV(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j3,j1))*Yv(j2,j1))*ZV(gt4,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gt3,j2))*SUM(j1,0,2,Conj(Yv(j1,j3))*Yv(j1,j2)))*ZV(gt4,j3)))*ZH(gt1,1)*ZH(gt2,1)) - 0.25*(0.6*Sqr(g1) + Sqr(g2))*SUM(j1,0,2,Conj(ZV(gt3,j1))*ZV(gt4,j1))*(ZH(gt1,0)*ZH(gt2,0) - ZH(gt1,1)*ZH(gt2,1));

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

   const std::complex<double> result = -0.35355339059327373*(Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt3,j1))*ZD(gt4,j1))*(ZH(gt1,0)*ZP(gt2,0) + ZH(gt1,1)*ZP(gt2,1)) - 2*(SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gt4,j3))*ZH(gt1,0)*ZP(gt2,0) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZD(gt4,j3))*ZH(gt1,1)*ZP(gt2,1) + SUM(j3,0,2,Conj(ZU(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yd(j2,j1))*ZD(gt4,3 + j2)))*(ZH(gt1,1)*ZP(gt2,0) + ZH(gt1,0)*ZP(gt2,1))));

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
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -0.35355339059327373*(Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt3,j1))*ZE(gt4,j1))*(ZH(gt1,0)*ZP(gt2,0) + ZH(gt1,1)*ZP(gt2,1)) - 2*(SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gt3,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZE(gt4,j3))*ZH(gt1,0)*ZP(gt2,0) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gt3,j2))*SUM(j1,0,2,Conj(Yv(j1,j3))*Yv(j1,j2)))*ZE(gt4,j3))*ZH(gt1,1)*ZP(gt2,1) + SUM(j3,0,2,Conj(ZV(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j3,j1))*Ye(j2,j1))*ZE(gt4,3 + j2)))*(ZH(gt1,1)*ZP(gt2,0) + ZH(gt1,0)*ZP(gt2,1))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Hpm, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto vu = MODELPARAMETER(vu);
   const auto g1 = MODELPARAMETER(g1);
   const auto vd = MODELPARAMETER(vd);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-(ZH(gt1,0)*(ZP(gt2,0)*(vd*(0.6*Sqr(g1) + Sqr(g2))*ZP(gt3,0) + vu*Sqr(g2)*ZP(gt3,1)) + ZP(gt2,1)*(vu*Sqr(g2)*ZP(gt3,0) + vd*(-0.6*Sqr(g1) + Sqr(g2))*ZP(gt3,1)))) - ZH(gt1,1)*(ZP(gt2,0)*(vu*(-0.6*Sqr(g1) + Sqr(g2))*ZP(gt3,0) + vd*Sqr(g2)*ZP(gt3,1)) + ZP(gt2,1)*(vd*Sqr(g2)*ZP(gt3,0) + vu*(0.6*Sqr(g1) + Sqr(g2))*ZP(gt3,1))));

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
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.3872983346207417*g1*g2*Sin(ThetaW)*(ZH(gt1,0)*ZP(gt2,0) - ZH(gt1,1)*ZP(gt2,1));

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

   const std::complex<double> result = 0.5*g2*(ZH(gt1,0)*ZP(gt2,0) - ZH(gt1,1)*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
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

   const std::complex<double> result = -0.35355339059327373*(Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZU(gt4,j1))*(ZH(gt1,0)*ZP(gt3,0) + ZH(gt1,1)*ZP(gt3,1)) - 2*(SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZU(gt4,j3))*ZH(gt1,0)*ZP(gt3,0) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gt4,j3))*ZH(gt1,1)*ZP(gt3,1) + SUM(j3,0,2,Conj(ZD(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yu(j2,j1))*ZU(gt4,3 + j2)))*(ZH(gt1,1)*ZP(gt3,0) + ZH(gt1,0)*ZP(gt3,1))));

   return {result};
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
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Mu = MODELPARAMETER(Mu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.05*((Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*(vd*ZH(gt1,0) - vu*ZH(gt1,1)) + 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*(vd*ZH(gt1,0) - vu*ZH(gt1,1)) - 10*(1.4142135623730951*SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,ZD(gt3,3 + j1)*TYd(j1,j2)))*ZH(gt1,0) + 1.4142135623730951*SUM(j2,0,2,SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j2)))*ZD(gt3,j2))*ZH(gt1,0) + 2*vd*SUM(j3,0,2,Conj(ZD(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,j1))*ZD(gt3,3 + j2)))*ZH(gt1,0) + 2*vd*SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gt3,j3))*ZH(gt1,0) - 1.4142135623730951*Conj(Mu)*SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1)))*ZH(gt1,1) - 1.4142135623730951*Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1)))*ZD(gt3,j2))*ZH(gt1,1)));

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
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -0.35355339059327373*(Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZV(gt4,j1))*(ZH(gt1,0)*ZP(gt3,0) + ZH(gt1,1)*ZP(gt3,1)) - 2*(SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZV(gt4,j3))*ZH(gt1,0)*ZP(gt3,0) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Conj(Yv(j1,j3))*Yv(j1,j2)))*ZV(gt4,j3))*ZH(gt1,1)*ZP(gt3,1) + SUM(j3,0,2,Conj(ZE(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Yv(j2,j1))*ZV(gt4,3 + j2)))*(ZH(gt1,1)*ZP(gt3,0) + ZH(gt1,0)*ZP(gt3,1))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Se, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto TYe = MODELPARAMETER(TYe);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Mu = MODELPARAMETER(Mu);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.05*(-((3*Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*(vd*ZH(gt1,0) - vu*ZH(gt1,1))) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*(vd*ZH(gt1,0) - vu*ZH(gt1,1)) - 10*(1.4142135623730951*SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,ZE(gt3,3 + j1)*TYe(j1,j2)))*ZH(gt1,0) + 1.4142135623730951*SUM(j2,0,2,SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j2)))*ZE(gt3,j2))*ZH(gt1,0) + 2*vd*SUM(j3,0,2,Conj(ZE(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Ye(j2,j1))*ZE(gt3,3 + j2)))*ZH(gt1,0) + 2*vd*SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZE(gt3,j3))*ZH(gt1,0) - 1.4142135623730951*Conj(Mu)*SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)))*ZH(gt1,1) - 1.4142135623730951*Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1)))*ZE(gt3,j2))*ZH(gt1,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Su, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto TYu = MODELPARAMETER(TYu);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Mu = MODELPARAMETER(Mu);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.05*((Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*(vd*ZH(gt1,0) - vu*ZH(gt1,1)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*(vd*ZH(gt1,0) - vu*ZH(gt1,1)) + 10*(1.4142135623730951*Conj(Mu)*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt3,3 + j1)))*ZH(gt1,0) + 1.4142135623730951*Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1)))*ZU(gt3,j2))*ZH(gt1,0) - (1.4142135623730951*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,ZU(gt3,3 + j1)*TYu(j1,j2))) + 1.4142135623730951*SUM(j2,0,2,SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j2)))*ZU(gt3,j2)) + 2*vu*(SUM(j3,0,2,Conj(ZU(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))*ZU(gt3,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gt3,j3))))*ZH(gt1,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Sv, fields::Sv>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Mv = MODELPARAMETER(Mv);
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = -0.35355339059327373*(SUM(j3,0,2,Conj(ZV(gt3,3 + j3))*SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Conj(Mv(j1,j3))*Yv(j1,j2)))) + SUM(j3,0,2,Conj(ZV(gt2,3 + j3))*SUM(j2,0,2,Conj(ZV(gt3,j2))*SUM(j1,0,2,Conj(Mv(j1,j3))*Yv(j1,j2)))) + SUM(j3,0,2,Conj(ZV(gt3,3 + j3))*SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Conj(Mv(j3,j1))*Yv(j1,j2)))) + SUM(j3,0,2,Conj(ZV(gt2,3 + j3))*SUM(j2,0,2,Conj(ZV(gt3,j2))*SUM(j1,0,2,Conj(Mv(j3,j1))*Yv(j1,j2)))))*ZH(gt1,1);

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Sv, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto TYv = MODELPARAMETER(TYv);
   const auto Mu = MODELPARAMETER(Mu);
   const auto vu = MODELPARAMETER(vu);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.7071067811865475*Conj(Mu)*SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Yv(j1,j2)*ZV(gt3,3 + j1)))*ZH(gt1,0) + 0.7071067811865475*Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j2))*Conj(ZV(gt2,3 + j1)))*ZV(gt3,j2))*ZH(gt1,0) - 0.5*(1.4142135623730951*SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,ZV(gt3,3 + j1)*TYv(j1,j2))) + 1.4142135623730951*SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(gt2,3 + j1))*Conj(TYv(j1,j2)))*ZV(gt3,j2)) + 2*vu*(SUM(j3,0,2,Conj(ZV(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j3,j1))*Yv(j2,j1))*ZV(gt3,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Conj(Yv(j1,j3))*Yv(j1,j2)))*ZV(gt3,j3))))*ZH(gt1,1) - 0.05*(3*Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZV(gt3,j1))*(vd*ZH(gt1,0) - vu*ZH(gt1,1));

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
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.3872983346207417*g1*g2*Sin(ThetaW)*(ZH(gt1,0)*ZP(gt2,0) - ZH(gt1,1)*ZP(gt2,1));

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

ScalarVertex VertexImpl<fields::hh, typename fields::conj<fields::Sv>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Yv = MODELPARAMETER(Yv);
   const auto Mv = MODELPARAMETER(Mv);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = -0.35355339059327373*(SUM(j3,0,2,SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j3))*Mv(j1,j2))*ZV(gt3,3 + j2))*ZV(gt2,j3)) + SUM(j3,0,2,SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j3))*Mv(j2,j1))*ZV(gt3,3 + j2))*ZV(gt2,j3)) + SUM(j3,0,2,SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j3))*Mv(j1,j2))*ZV(gt2,3 + j2))*ZV(gt3,j3)) + SUM(j3,0,2,SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j3))*Mv(j2,j1))*ZV(gt2,3 + j2))*ZV(gt3,j3)))*ZH(gt1,1);

   return {result};
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

ScalarVertex VertexImpl<fields::Hpm, fields::Hpm, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -0.25*(0.6*Sqr(g1) + Sqr(g2))*(-(ZP(gt1,1)*(-2*ZP(gt2,1)*ZP(gt3,1)*ZP(gt4,1) + ZP(gt2,0)*(ZP(gt3,1)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,1)))) + ZP(gt1,0)*(2*ZP(gt2,0)*ZP(gt3,0)*ZP(gt4,0) - ZP(gt2,1)*(ZP(gt3,1)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,1))));

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

   const std::complex<double> result = 0.05*((Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt4,j1))*(ZP(gt1,0)*ZP(gt3,0) - ZP(gt1,1)*ZP(gt3,1)) + 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt4,3 + j1))*(ZP(gt1,0)*ZP(gt3,0) - ZP(gt1,1)*ZP(gt3,1)) - 20*(SUM(j3,0,2,Conj(ZD(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,j1))*ZD(gt4,3 + j2)))*ZP(gt1,0)*ZP(gt3,0) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZD(gt4,j3))*ZP(gt1,1)*ZP(gt3,1)));

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
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.05*(-((3*Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt4,j1))*(ZP(gt1,0)*ZP(gt3,0) - ZP(gt1,1)*ZP(gt3,1))) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt4,3 + j1))*(ZP(gt1,0)*ZP(gt3,0) - ZP(gt1,1)*ZP(gt3,1)) - 20*(SUM(j3,0,2,Conj(ZE(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Ye(j2,j1))*ZE(gt4,3 + j2)))*ZP(gt1,0)*ZP(gt3,0) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Conj(Yv(j1,j3))*Yv(j1,j2)))*ZE(gt4,j3))*ZP(gt1,1)*ZP(gt3,1)));

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

   const std::complex<double> result = 0.05*((Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt4,j1))*(ZP(gt1,0)*ZP(gt3,0) - ZP(gt1,1)*ZP(gt3,1)) - 4*(Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt4,3 + j1))*(ZP(gt1,0)*ZP(gt3,0) - ZP(gt1,1)*ZP(gt3,1)) + 5*(SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZU(gt4,j3))*ZP(gt1,0)*ZP(gt3,0) + SUM(j3,0,2,Conj(ZU(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))*ZU(gt4,3 + j2)))*ZP(gt1,1)*ZP(gt3,1))));

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
   const auto Mu = MODELPARAMETER(Mu);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-1.4142135623730951*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZD(gt3,j1))*(vd*ZP(gt1,0) + vu*ZP(gt1,1)) + 2*(2*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,ZD(gt3,3 + j1)*TYd(j1,j2)))*ZP(gt1,0) + 2*Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1)))*ZD(gt3,j2))*ZP(gt1,0) + 1.4142135623730951*vu*SUM(j3,0,2,Conj(ZU(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yd(j2,j1))*ZD(gt3,3 + j2)))*ZP(gt1,0) + 1.4142135623730951*vd*SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gt3,j3))*ZP(gt1,0) + 2*Conj(Mu)*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1)))*ZP(gt1,1) + 2*SUM(j2,0,2,SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j2)))*ZD(gt3,j2))*ZP(gt1,1) + 1.4142135623730951*vd*SUM(j3,0,2,Conj(ZU(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yd(j2,j1))*ZD(gt3,3 + j2)))*ZP(gt1,1) + 1.4142135623730951*vu*SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZD(gt3,j3))*ZP(gt1,1)));

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
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -(SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZV(gt4,j3))*ZP(gt1,0)*ZP(gt3,0)) - SUM(j3,0,2,Conj(ZV(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j3,j1))*Yv(j2,j1))*ZV(gt4,3 + j2)))*ZP(gt1,1)*ZP(gt3,1) - 0.05*(3*Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZV(gt4,j1))*(ZP(gt1,0)*ZP(gt3,0) - ZP(gt1,1)*ZP(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::Sv, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto TYe = MODELPARAMETER(TYe);
   const auto TYv = MODELPARAMETER(TYv);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Mu = MODELPARAMETER(Mu);
   const auto Yv = MODELPARAMETER(Yv);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-1.4142135623730951*Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZE(gt3,j1))*(vd*ZP(gt1,0) + vu*ZP(gt1,1)) + 2*(2*SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,ZE(gt3,3 + j1)*TYe(j1,j2)))*ZP(gt1,0) + 2*Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j2))*Conj(ZV(gt2,3 + j1)))*ZE(gt3,j2))*ZP(gt1,0) + 1.4142135623730951*vu*SUM(j3,0,2,Conj(ZV(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j3,j1))*Ye(j2,j1))*ZE(gt3,3 + j2)))*ZP(gt1,0) + 1.4142135623730951*vd*SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZE(gt3,j3))*ZP(gt1,0) + 2*Conj(Mu)*SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)))*ZP(gt1,1) + 2*SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(gt2,3 + j1))*Conj(TYv(j1,j2)))*ZE(gt3,j2))*ZP(gt1,1) + 1.4142135623730951*vd*SUM(j3,0,2,Conj(ZV(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j3,j1))*Ye(j2,j1))*ZE(gt3,3 + j2)))*ZP(gt1,1) + 1.4142135623730951*vu*SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Conj(Yv(j1,j3))*Yv(j1,j2)))*ZE(gt3,j3))*ZP(gt1,1)));

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

InverseMetricVertex VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, fields::VP, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*(7.745966692414834*g1*g2*Cos(2*ThetaW) + Sin(2*ThetaW)*(-3*Sqr(g1) + 5*Sqr(g2)))*(ZP(gt1,0)*ZP(gt2,0) + ZP(gt1,1)*ZP(gt2,1));

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

InverseMetricVertex VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*Sqr(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*(ZP(gt1,0)*ZP(gt2,0) + ZP(gt1,1)*ZP(gt2,1));

   return {result};
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

   const std::complex<double> result = 0.5*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*(ZP(gt1,0)*ZP(gt2,0) + ZP(gt1,1)*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
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

ScalarVertex VertexImpl<fields::Hpm, typename fields::conj<fields::Se>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Yv = MODELPARAMETER(Yv);
   const auto Mv = MODELPARAMETER(Mv);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.5*(SUM(j3,0,2,SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j3))*Mv(j1,j2))*ZV(gt3,3 + j2))*ZE(gt2,j3)) + SUM(j3,0,2,SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j3))*Mv(j2,j1))*ZV(gt3,3 + j2))*ZE(gt2,j3)))*ZP(gt1,1);

   return {result};
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
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.3872983346207417*g1*g2*Sin(ThetaW)*(vd*ZP(gt1,0) - vu*ZP(gt1,1));

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
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = 0.025*(SUM(j1,0,2,Conj(ZV(gt2,j1))*ZV(gt4,j1))*((Sqr(g1) + 5*Sqr(g2))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) + 2*Sqr(g1)*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2))) + ((Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1)) + 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1)))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt4,j2)));

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
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yv = MODELPARAMETER(Yv);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);

   const std::complex<double> result = 0.25*(-(Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZU(gt4,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZE(gt3,j2))) - Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZU(gt4,j2)) - 4*(SUM(j2,0,2,Conj(ZD(gt1,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt4,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yv(j3,j4))*Conj(ZV(gt2,3 + j3)))*ZE(gt3,j4)) + SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yd(j3,j4))*Conj(ZD(gt1,3 + j3)))*ZU(gt4,j4))));

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
   const auto Mu = MODELPARAMETER(Mu);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-1.4142135623730951*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZU(gt3,j1))*(vd*ZP(gt2,0) + vu*ZP(gt2,1)) + 2*(2*Conj(Mu)*SUM(j2,0,2,Conj(ZD(gt1,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt3,3 + j1)))*ZP(gt2,0) + 2*SUM(j2,0,2,SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*Conj(TYd(j1,j2)))*ZU(gt3,j2))*ZP(gt2,0) + 1.4142135623730951*vu*SUM(j3,0,2,Conj(ZD(gt1,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yu(j2,j1))*ZU(gt3,3 + j2)))*ZP(gt2,0) + 1.4142135623730951*vd*SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt1,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZU(gt3,j3))*ZP(gt2,0) + 2*SUM(j2,0,2,Conj(ZD(gt1,j2))*SUM(j1,0,2,ZU(gt3,3 + j1)*TYu(j1,j2)))*ZP(gt2,1) + 2*Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt1,3 + j1)))*ZU(gt3,j2))*ZP(gt2,1) + 1.4142135623730951*vd*SUM(j3,0,2,Conj(ZD(gt1,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yu(j2,j1))*ZU(gt3,3 + j2)))*ZP(gt2,1) + 1.4142135623730951*vu*SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt1,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gt3,j3))*ZP(gt2,1)));

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

ScalarVertex VertexImpl<fields::Se, fields::Su, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto Yv = MODELPARAMETER(Yv);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = 0.25*(-(Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZV(gt4,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZD(gt3,j2))) - Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZV(gt4,j2)) - 4*(SUM(j2,0,2,Conj(ZE(gt1,j2))*SUM(j1,0,2,Yv(j1,j2)*ZV(gt4,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yu(j3,j4))*Conj(ZU(gt2,3 + j3)))*ZD(gt3,j4)) + SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Ye(j3,j4))*Conj(ZE(gt1,3 + j3)))*ZV(gt4,j4))));

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

ScalarVertex VertexImpl<fields::Se, fields::Sv, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Mv = MODELPARAMETER(Mv);
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.5*(SUM(j3,0,2,Conj(ZV(gt2,3 + j3))*SUM(j2,0,2,Conj(ZE(gt1,j2))*SUM(j1,0,2,Conj(Mv(j1,j3))*Yv(j1,j2)))) + SUM(j3,0,2,Conj(ZV(gt2,3 + j3))*SUM(j2,0,2,Conj(ZE(gt1,j2))*SUM(j1,0,2,Conj(Mv(j3,j1))*Yv(j1,j2)))))*ZP(gt3,1);

   return {result};
}

ScalarVertex VertexImpl<fields::Se, fields::Sv, typename fields::conj<fields::Se>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto Yv = MODELPARAMETER(Yv);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);

   const std::complex<double> result = 0.025*(-10*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZV(gt4,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZE(gt3,j2)) + SUM(j1,0,2,Conj(ZV(gt2,j1))*ZV(gt4,j1))*((-3*Sqr(g1) + 5*Sqr(g2))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt3,j2)) + 6*Sqr(g1)*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt3,3 + j2))) - 10*Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZV(gt4,j2)) - 3*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt4,j2)) + 5*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt4,j2)) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt4,j2)) - 40*SUM(j2,0,2,Conj(ZE(gt1,j2))*SUM(j1,0,2,Yv(j1,j2)*ZV(gt4,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yv(j3,j4))*Conj(ZV(gt2,3 + j3)))*ZE(gt3,j4)) - 40*SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Ye(j3,j4))*Conj(ZE(gt1,3 + j3)))*ZV(gt4,j4)));

   return {result};
}

ScalarVertex VertexImpl<fields::Se, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto TYe = MODELPARAMETER(TYe);
   const auto TYv = MODELPARAMETER(TYv);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Mu = MODELPARAMETER(Mu);
   const auto Yv = MODELPARAMETER(Yv);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-1.4142135623730951*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZV(gt3,j1))*(vd*ZP(gt2,0) + vu*ZP(gt2,1)) + 2*(2*Conj(Mu)*SUM(j2,0,2,Conj(ZE(gt1,j2))*SUM(j1,0,2,Yv(j1,j2)*ZV(gt3,3 + j1)))*ZP(gt2,0) + 2*SUM(j2,0,2,SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*Conj(TYe(j1,j2)))*ZV(gt3,j2))*ZP(gt2,0) + 1.4142135623730951*vu*SUM(j3,0,2,Conj(ZE(gt1,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Yv(j2,j1))*ZV(gt3,3 + j2)))*ZP(gt2,0) + 1.4142135623730951*vd*SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gt1,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZV(gt3,j3))*ZP(gt2,0) + 2*SUM(j2,0,2,Conj(ZE(gt1,j2))*SUM(j1,0,2,ZV(gt3,3 + j1)*TYv(j1,j2)))*ZP(gt2,1) + 2*Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt1,3 + j1)))*ZV(gt3,j2))*ZP(gt2,1) + 1.4142135623730951*vd*SUM(j3,0,2,Conj(ZE(gt1,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Yv(j2,j1))*ZV(gt3,3 + j2)))*ZP(gt2,1) + 1.4142135623730951*vu*SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gt1,j2))*SUM(j1,0,2,Conj(Yv(j1,j3))*Yv(j1,j2)))*ZV(gt3,j3))*ZP(gt2,1)));

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

ScalarVertex VertexImpl<fields::Su, fields::Sv, typename fields::conj<fields::Su>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yv = MODELPARAMETER(Yv);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZU = MODELPARAMETER(ZU);

   const std::complex<double> result = 0.025*(SUM(j1,0,2,Conj(ZV(gt2,j1))*ZV(gt4,j1))*((Sqr(g1) - 5*Sqr(g2))*SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt3,j2)) - 4*Sqr(g1)*SUM(j2,0,2,Conj(ZU(gt1,3 + j2))*ZU(gt3,3 + j2))) + (Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt4,j2)) - 4*(Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt3,3 + j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt4,j2)) + 10*SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Yv(j1,j2)*ZV(gt4,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yu(j3,j4))*Conj(ZU(gt1,3 + j3)))*ZU(gt3,j4)) + 10*SUM(j2,0,2,Conj(ZU(gt1,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt3,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yv(j3,j4))*Conj(ZV(gt2,3 + j3)))*ZV(gt4,j4))));

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

InverseMetricVertex VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VG, fields::VG>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g3 = MODELPARAMETER(g3);

   const std::complex<double> result = 0.25*KroneckerDelta(gt1,gt2)*Sqr(g3);

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

ScalarVertex VertexImpl<fields::Sv, fields::Sv, typename fields::conj<fields::Sv>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZV = MODELPARAMETER(ZV);

   const std::complex<double> result = 0.125*(-((0.6*Sqr(g1) + Sqr(g2))*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZV(gt4,j1))*SUM(j2,0,2,Conj(ZV(gt1,j2))*ZV(gt3,j2))) - (0.6*Sqr(g1) + Sqr(g2))*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZV(gt4,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt3,j2)) - 0.6*Sqr(g1)*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZV(gt3,j1))*SUM(j2,0,2,Conj(ZV(gt1,j2))*ZV(gt4,j2)) - Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZV(gt3,j1))*SUM(j2,0,2,Conj(ZV(gt1,j2))*ZV(gt4,j2)) - 0.6*Sqr(g1)*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZV(gt3,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt4,j2)) - Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZV(gt3,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt4,j2)) - 8*SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Yv(j1,j2)*ZV(gt4,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yv(j3,j4))*Conj(ZV(gt1,3 + j3)))*ZV(gt3,j4)) - 8*SUM(j2,0,2,Conj(ZV(gt1,j2))*SUM(j1,0,2,Yv(j1,j2)*ZV(gt4,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yv(j3,j4))*Conj(ZV(gt2,3 + j3)))*ZV(gt3,j4)) - 8*SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Yv(j1,j2)*ZV(gt3,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yv(j3,j4))*Conj(ZV(gt1,3 + j3)))*ZV(gt4,j4)) - 8*SUM(j2,0,2,Conj(ZV(gt1,j2))*SUM(j1,0,2,Yv(j1,j2)*ZV(gt3,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yv(j3,j4))*Conj(ZV(gt2,3 + j3)))*ZV(gt4,j4)));

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

InverseMetricVertex VertexImpl<fields::Sv, typename fields::conj<fields::Sv>::type, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*Sqr(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZV(gt2,j1));

   return {result};
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
   const auto ZV = MODELPARAMETER(ZV);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.5*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZV(gt2,j1));

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::Sv, typename fields::conj<fields::Sv>::type, typename fields::conj<fields::VWm>::type, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZV = MODELPARAMETER(ZV);

   const std::complex<double> result = 0.5*Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZV(gt2,j1));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::VWm, fields::Ah, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,-0.5)*g2*(ZA(gt1,0)*ZP(gt2,0) + ZA(gt1,1)*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

ChiralVertex VertexImpl<fields::VWm, fields::Chi, typename fields::bar<fields::Cha>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto UP = MODELPARAMETER(UP);
   const auto ZN = MODELPARAMETER(ZN);
   const auto UM = MODELPARAMETER(UM);

   const std::complex<double> left = g2*Conj(UP(gt1,0))*ZN(gt2,1) - 0.7071067811865475*g2*Conj(UP(gt1,1))*ZN(gt2,3);

   const std::complex<double> right = g2*Conj(ZN(gt2,1))*UM(gt1,0) + 0.7071067811865475*g2*Conj(ZN(gt2,2))*UM(gt1,1);

   return {left, right};
}

MomentumDifferenceVertex VertexImpl<fields::VWm, fields::hh, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -0.5*g2*(ZH(gt1,0)*ZP(gt2,0) - ZH(gt1,1)*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::VWm, fields::hh, typename fields::conj<fields::VWm>::type>::evaluate(
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

MomentumDifferenceVertex VertexImpl<fields::VWm, fields::Su, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = -0.7071067811865475*g2*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZD(gt2,j1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VWm, fields::Sv, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);

   const std::complex<double> result = -0.7071067811865475*g2*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZE(gt2,j1));

   return {result, minuend_index, subtrahend_index};
}

ChiralVertex VertexImpl<fields::VWm, typename fields::bar<fields::Cha>::type, fields::Chi>::evaluate(
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

MomentumDifferenceVertex VertexImpl<fields::VWm, typename fields::conj<fields::Hpm>::type, fields::Ah>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,-0.5)*g2*(ZA(gt1,0)*ZP(gt2,0) + ZA(gt1,1)*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VWm, typename fields::conj<fields::Hpm>::type, fields::hh>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -0.5*g2*(ZH(gt1,0)*ZP(gt2,0) - ZH(gt1,1)*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::VWm, typename fields::conj<fields::Hpm>::type, fields::VP>::evaluate(
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

InverseMetricVertex VertexImpl<fields::VWm, typename fields::conj<fields::Hpm>::type, fields::VZ>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.3872983346207417*g1*g2*Sin(ThetaW)*(vd*ZP(gt1,0) - vu*ZP(gt1,1));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::VWm, typename fields::conj<fields::Sd>::type, fields::Su>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = -0.7071067811865475*g2*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZD(gt2,j1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VWm, typename fields::conj<fields::Se>::type, fields::Sv>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);

   const std::complex<double> result = -0.7071067811865475*g2*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZE(gt2,j1));

   return {result, minuend_index, subtrahend_index};
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

ChiralVertex VertexImpl<fields::VZ, fields::Cha, typename fields::bar<fields::Cha>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UP = MODELPARAMETER(UP);
   const auto UM = MODELPARAMETER(UM);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -(g2*Conj(UP(gt1,0))*Cos(ThetaW)*UP(gt2,0)) + 0.1*Conj(UP(gt1,1))*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*UP(gt2,1);

   const std::complex<double> right = -(g2*Conj(UM(gt2,0))*Cos(ThetaW)*UM(gt1,0)) + 0.1*Conj(UM(gt2,1))*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*UM(gt1,1);

   return {left, right};
}

ChiralVertex VertexImpl<fields::VZ, fields::Chi, fields::Chi>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.5*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(Conj(ZN(gt2,2))*ZN(gt1,2) - Conj(ZN(gt2,3))*ZN(gt1,3));

   const std::complex<double> right = 0.5*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(Conj(ZN(gt1,2))*ZN(gt2,2) - Conj(ZN(gt1,3))*ZN(gt2,3));

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

   const std::complex<double> result = 0.5*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*(ZP(gt1,0)*ZP(gt2,0) + ZP(gt1,1)*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::VZ, fields::Hpm, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.3872983346207417*g1*g2*Sin(ThetaW)*(vd*ZP(gt1,0) - vu*ZP(gt1,1));

   return {result};
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
   const auto ZV = MODELPARAMETER(ZV);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.5*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZV(gt2,j1));

   return {result, minuend_index, subtrahend_index};
}

ChiralVertex VertexImpl<fields::VZ, typename fields::bar<fields::Cha>::type, fields::Cha>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UM = MODELPARAMETER(UM);
   const auto UP = MODELPARAMETER(UP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = g2*Conj(UM(gt2,0))*Cos(ThetaW)*UM(gt1,0) + 0.5*Conj(UM(gt2,1))*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*UM(gt1,1);

   const std::complex<double> right = g2*Conj(UP(gt1,0))*Cos(ThetaW)*UP(gt2,0) + 0.5*Conj(UP(gt1,1))*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*UP(gt2,1);

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

   const std::complex<double> result = 0.5*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*(ZP(gt1,0)*ZP(gt2,0) + ZP(gt1,1)*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::VZ, typename fields::conj<fields::Hpm>::type, fields::VWm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.3872983346207417*g1*g2*Sin(ThetaW)*(vd*ZP(gt1,0) - vu*ZP(gt1,1));

   return {result};
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
   const auto ZV = MODELPARAMETER(ZV);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.5*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZV(gt2,j1));

   return {result, minuend_index, subtrahend_index};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha>::type, fields::Cha, fields::Ah>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto UM = MODELPARAMETER(UM);
   const auto UP = MODELPARAMETER(UP);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,0.7071067811865475)*g2*(Conj(UM(gt2,1))*Conj(UP(gt1,0))*ZA(gt3,0) + Conj(UM(gt2,0))*Conj(UP(gt1,1))*ZA(gt3,1));

   const std::complex<double> right = std::complex<double>(0.,-0.7071067811865475)*g2*(UM(gt1,1)*UP(gt2,0)*ZA(gt3,0) + UM(gt1,0)*UP(gt2,1)*ZA(gt3,1));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha>::type, fields::Cha, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto UM = MODELPARAMETER(UM);
   const auto UP = MODELPARAMETER(UP);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*g2*(Conj(UM(gt2,1))*Conj(UP(gt1,0))*ZH(gt3,0) + Conj(UM(gt2,0))*Conj(UP(gt1,1))*ZH(gt3,1));

   const std::complex<double> right = -0.7071067811865475*g2*(UM(gt1,1)*UP(gt2,0)*ZH(gt3,0) + UM(gt1,0)*UP(gt2,1)*ZH(gt3,1));

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

ChiralVertex VertexImpl<typename fields::bar<fields::Cha>::type, fields::Cha, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UM = MODELPARAMETER(UM);
   const auto UP = MODELPARAMETER(UP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = g2*Conj(UM(gt2,0))*Cos(ThetaW)*UM(gt1,0) + 0.5*Conj(UM(gt2,1))*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*UM(gt1,1);

   const std::complex<double> right = g2*Conj(UP(gt1,0))*Cos(ThetaW)*UP(gt2,0) + 0.5*Conj(UP(gt1,1))*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*UP(gt2,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha>::type, fields::Chi, fields::Hpm>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto UP = MODELPARAMETER(UP);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZP = MODELPARAMETER(ZP);
   const auto UM = MODELPARAMETER(UM);

   const std::complex<double> left = -0.5*(Conj(UP(gt1,1))*(1.0954451150103321*g1*Conj(ZN(gt2,0)) + 1.4142135623730951*g2*Conj(ZN(gt2,1))) + 2*g2*Conj(UP(gt1,0))*Conj(ZN(gt2,3)))*ZP(gt3,1);

   const std::complex<double> right = 0.5*(UM(gt1,1)*(1.0954451150103321*g1*ZN(gt2,0) + 1.4142135623730951*g2*ZN(gt2,1)) - 2*g2*UM(gt1,0)*ZN(gt2,2))*ZP(gt3,0);

   return {left, right};
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

   const std::complex<double> left = -(g2*Conj(UP(gt1,0))*SUM(j1,0,2,Conj(ZDL(gt2,j1))*ZU(gt3,j1))) + Conj(UP(gt1,1))*SUM(j2,0,2,Conj(ZDL(gt2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt3,3 + j1)));

   const std::complex<double> right = SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*ZDR(gt2,j1))*ZU(gt3,j2))*UM(gt1,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha>::type, fields::Fe, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Yv = MODELPARAMETER(Yv);
   const auto Ye = MODELPARAMETER(Ye);
   const auto UP = MODELPARAMETER(UP);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZER = MODELPARAMETER(ZER);
   const auto UM = MODELPARAMETER(UM);

   const std::complex<double> left = -(g2*Conj(UP(gt1,0))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZV(gt3,j1))) + Conj(UP(gt1,1))*SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Yv(j1,j2)*ZV(gt3,3 + j1)));

   const std::complex<double> right = SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZV(gt3,j2))*UM(gt1,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha>::type, fields::Fv, fields::Se>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Yv = MODELPARAMETER(Yv);
   const auto Ye = MODELPARAMETER(Ye);
   const auto UP = MODELPARAMETER(UP);
   const auto ZE = MODELPARAMETER(ZE);
   const auto UV = MODELPARAMETER(UV);
   const auto UM = MODELPARAMETER(UM);

   const std::complex<double> left = Conj(UP(gt1,1))*SUM(j2,0,2,Conj(ZE(gt3,j2))*SUM(j1,0,2,Conj(UV(gt2,3 + j1))*Yv(j1,j2)));

   const std::complex<double> right = -(g2*SUM(j1,0,2,Conj(ZE(gt3,j1))*UV(gt2,j1))*UM(gt1,0)) + SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt3,3 + j1)))*UV(gt2,j2))*UM(gt1,1);

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

   const std::complex<double> left = Conj(UP(gt1,1))*SUM(j2,0,2,Conj(ZD(gt3,j2))*SUM(j1,0,2,Conj(ZUR(gt2,j1))*Yu(j1,j2)));

   const std::complex<double> right = -(g2*SUM(j1,0,2,Conj(ZD(gt3,j1))*ZUL(gt2,j1))*UM(gt1,0)) + SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt3,3 + j1)))*ZUL(gt2,j2))*UM(gt1,1);

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
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZDR = MODELPARAMETER(ZDR);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto UP = MODELPARAMETER(UP);

   const std::complex<double> left = Conj(UM(gt2,1))*SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Conj(ZDR(gt1,j1))*Yd(j1,j2)));

   const std::complex<double> right = -(g2*SUM(j1,0,2,Conj(ZU(gt3,j1))*ZDL(gt1,j1))*UP(gt2,0)) + SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt3,3 + j1)))*ZDL(gt1,j2))*UP(gt2,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Chi, fields::Sd>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZDR = MODELPARAMETER(ZDR);
   const auto ZDL = MODELPARAMETER(ZDL);

   const std::complex<double> left = -0.3651483716701107*g1*Conj(ZN(gt2,0))*SUM(j1,0,2,Conj(ZD(gt3,3 + j1))*Conj(ZDR(gt1,j1))) - Conj(ZN(gt2,2))*SUM(j2,0,2,Conj(ZD(gt3,j2))*SUM(j1,0,2,Conj(ZDR(gt1,j1))*Yd(j1,j2)));

   const std::complex<double> right = SUM(j1,0,2,Conj(ZD(gt3,j1))*ZDL(gt1,j1))*(-0.18257418583505536*g1*ZN(gt2,0) + 0.7071067811865475*g2*ZN(gt2,1)) - SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt3,3 + j1)))*ZDL(gt1,j2))*ZN(gt2,2);

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
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.16666666666666666*KroneckerDelta(gt1,gt2)*(3*g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW));

   const std::complex<double> right = -0.2581988897471611*g1*KroneckerDelta(gt1,gt2)*Sin(ThetaW);

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

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Cha, fields::Sv>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Yv = MODELPARAMETER(Yv);
   const auto UM = MODELPARAMETER(UM);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto UP = MODELPARAMETER(UP);

   const std::complex<double> left = Conj(UM(gt2,1))*SUM(j2,0,2,Conj(ZV(gt3,j2))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,j2)));

   const std::complex<double> right = -(g2*SUM(j1,0,2,Conj(ZV(gt3,j1))*ZEL(gt1,j1))*UP(gt2,0)) + SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j2))*Conj(ZV(gt3,3 + j1)))*ZEL(gt1,j2))*UP(gt2,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Chi, fields::Se>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZEL = MODELPARAMETER(ZEL);

   const std::complex<double> left = -1.0954451150103321*g1*Conj(ZN(gt2,0))*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*Conj(ZER(gt1,j1))) - Conj(ZN(gt2,2))*SUM(j2,0,2,Conj(ZE(gt3,j2))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,j2)));

   const std::complex<double> right = 0.7071067811865475*SUM(j1,0,2,Conj(ZE(gt3,j1))*ZEL(gt1,j1))*(0.7745966692414834*g1*ZN(gt2,0) + g2*ZN(gt2,1)) - SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt3,3 + j1)))*ZEL(gt1,j2))*ZN(gt2,2);

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

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fv, fields::Hpm>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Yv = MODELPARAMETER(Yv);
   const auto UV = MODELPARAMETER(UV);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ZEL = MODELPARAMETER(ZEL);

   const std::complex<double> left = SUM(j2,0,2,Conj(UV(gt2,j2))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,j2)))*ZP(gt3,0);

   const std::complex<double> right = SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j2))*UV(gt2,3 + j1))*ZEL(gt1,j2))*ZP(gt3,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fv, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto UV = MODELPARAMETER(UV);
   const auto ZEL = MODELPARAMETER(ZEL);

   const std::complex<double> left = -0.7071067811865475*g2*SUM(j1,0,2,Conj(UV(gt2,j1))*ZEL(gt1,j1));

   const std::complex<double> right = 0;

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Hpm, fields::Fv>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Yv = MODELPARAMETER(Yv);
   const auto UV = MODELPARAMETER(UV);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ZEL = MODELPARAMETER(ZEL);

   const std::complex<double> left = SUM(j2,0,2,Conj(UV(gt2,j2))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,j2)))*ZP(gt3,0);

   const std::complex<double> right = SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j2))*UV(gt2,3 + j1))*ZEL(gt1,j2))*ZP(gt3,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Se, fields::Chi>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZEL = MODELPARAMETER(ZEL);

   const std::complex<double> left = -1.0954451150103321*g1*Conj(ZN(gt2,0))*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*Conj(ZER(gt1,j1))) - Conj(ZN(gt2,2))*SUM(j2,0,2,Conj(ZE(gt3,j2))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,j2)));

   const std::complex<double> right = 0.7071067811865475*SUM(j1,0,2,Conj(ZE(gt3,j1))*ZEL(gt1,j1))*(0.7745966692414834*g1*ZN(gt2,0) + g2*ZN(gt2,1)) - SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt3,3 + j1)))*ZEL(gt1,j2))*ZN(gt2,2);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, fields::Chi, fields::Su>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZUR = MODELPARAMETER(ZUR);
   const auto ZUL = MODELPARAMETER(ZUL);

   const std::complex<double> left = 0.7302967433402214*g1*Conj(ZN(gt2,0))*SUM(j1,0,2,Conj(ZU(gt3,3 + j1))*Conj(ZUR(gt1,j1))) - Conj(ZN(gt2,3))*SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Conj(ZUR(gt1,j1))*Yu(j1,j2)));

   const std::complex<double> right = SUM(j1,0,2,Conj(ZU(gt3,j1))*ZUL(gt1,j1))*(-0.18257418583505536*g1*ZN(gt2,0) - 0.7071067811865475*g2*ZN(gt2,1)) - SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt3,3 + j1)))*ZUL(gt1,j2))*ZN(gt2,3);

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
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.03333333333333333*KroneckerDelta(gt1,gt2)*(-15*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW));

   const std::complex<double> right = 0.5163977794943222*g1*KroneckerDelta(gt1,gt2)*Sin(ThetaW);

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

MomentumVertex VertexImpl<typename fields::bar<fields::gP>::type, fields::gWmC, fields::VWm>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Sin(ThetaW));

   return {result, 1};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gP>::type, fields::gWm, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Sin(ThetaW);

   return {result, 1};
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

MomentumVertex VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gP, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Sin(ThetaW));

   return {result, 1};
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

ScalarVertex VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gZ, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.25*g2*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*(vd*ZP(gt3,0) - vu*ZP(gt3,1));

   return {result};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gZ, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Cos(ThetaW));

   return {result, 1};
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

MomentumVertex VertexImpl<typename fields::bar<fields::gWm>::type, fields::gP, fields::VWm>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Sin(ThetaW);

   return {result, 1};
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

ScalarVertex VertexImpl<typename fields::bar<fields::gWm>::type, fields::gZ, fields::Hpm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.25*g2*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*(vd*ZP(gt3,0) - vu*ZP(gt3,1));

   return {result};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWm>::type, fields::gZ, fields::VWm>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Cos(ThetaW);

   return {result, 1};
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

MomentumVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gWmC, fields::VWm>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Cos(ThetaW));

   return {result, 1};
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

MomentumVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gWm, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Cos(ThetaW);

   return {result, 1};
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
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.3872983346207417*g1*g2*Sin(ThetaW)*(vd*ZP(gt1,0) - vu*ZP(gt1,1));

   return {result};
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

ChiralVertex VertexImpl<typename fields::conj<fields::Sv>::type, typename fields::bar<fields::Cha>::type, fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Yv = MODELPARAMETER(Yv);
   const auto Ye = MODELPARAMETER(Ye);
   const auto UP = MODELPARAMETER(UP);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZER = MODELPARAMETER(ZER);
   const auto UM = MODELPARAMETER(UM);

   const std::complex<double> left = -(g2*Conj(UP(gt1,0))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZV(gt3,j1))) + Conj(UP(gt1,1))*SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Yv(j1,j2)*ZV(gt3,3 + j1)));

   const std::complex<double> right = SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZV(gt3,j2))*UM(gt1,1);

   return {left, right};
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

TripleVectorVertex VertexImpl<typename fields::conj<fields::VWm>::type, fields::VP, fields::VWm>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Sin(ThetaW);

   return {result, TripleVectorVertex::odd_permutation{}};
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

TripleVectorVertex VertexImpl<typename fields::conj<fields::VWm>::type, fields::VWm, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Cos(ThetaW));

   return {result, TripleVectorVertex::odd_permutation{}};
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

} // namespace detail
} // namespace MSSMRHN_cxx_diagrams
} // namespace flexiblesusy
