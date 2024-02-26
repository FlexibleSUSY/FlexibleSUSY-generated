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
 * @file cxx_qft/THDMII_vertices.cpp
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#include "THDMII_context_base.hpp"
#include "THDMII_input_parameters.hpp"
#include "THDMII_vertices.hpp"

#include "concatenate.hpp"
#include "wrappers.hpp"

#define INPUTPARAMETER(p) context.model.get_input_parameters().p
#define MODELPARAMETER(p) context.model.get_##p()
#define DERIVEDPARAMETER(p) context.model.p()
#define PHASE(p) context.model.get_##p()

namespace flexiblesusy::THDMII_cxx_diagrams::detail {

cxx_diagrams::ChiralVertex VertexImpl<THDMII_cxx_diagrams::fields::Ah, typename THDMII_cxx_diagrams::fields::bar<THDMII_cxx_diagrams::fields::Fd>::type, THDMII_cxx_diagrams::fields::Fd>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Yd = MODELPARAMETER(Yd);
   const auto Vd = MODELPARAMETER(Vd);
   const auto Ud = MODELPARAMETER(Ud);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,Conj(Vd(gt2,j2))*SUM(j1,0,2,Conj(Ud(gt1,j1))*Yd(j1,j2)))*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Ud(gt2,j1))*Vd(gt1,j2))*ZA(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<THDMII_cxx_diagrams::fields::Ah, typename THDMII_cxx_diagrams::fields::bar<THDMII_cxx_diagrams::fields::Fe>::type, THDMII_cxx_diagrams::fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ve = MODELPARAMETER(Ve);
   const auto Ue = MODELPARAMETER(Ue);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,Conj(Ve(gt2,j2))*SUM(j1,0,2,Conj(Ue(gt1,j1))*Ye(j1,j2)))*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gt2,j1))*Ve(gt1,j2))*ZA(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<THDMII_cxx_diagrams::fields::Ah, typename THDMII_cxx_diagrams::fields::bar<THDMII_cxx_diagrams::fields::Fu>::type, THDMII_cxx_diagrams::fields::Fu>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Yu = MODELPARAMETER(Yu);
   const auto Vu = MODELPARAMETER(Vu);
   const auto Uu = MODELPARAMETER(Uu);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,Conj(Vu(gt2,j2))*SUM(j1,0,2,Conj(Uu(gt1,j1))*Yu(j1,j2)))*ZA(gt3,1);

   const std::complex<double> right = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Uu(gt2,j1))*Vu(gt1,j2))*ZA(gt3,1);

   return {left, right};
}

cxx_diagrams::ScalarVertex VertexImpl<THDMII_cxx_diagrams::fields::Ah, typename THDMII_cxx_diagrams::fields::conj<THDMII_cxx_diagrams::fields::Hm>::type, THDMII_cxx_diagrams::fields::Hm>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto v2 = MODELPARAMETER(v2);
   const auto v1 = MODELPARAMETER(v1);
   const auto Lambda6 = MODELPARAMETER(Lambda6);
   const auto Lambda4 = MODELPARAMETER(Lambda4);
   const auto Lambda5 = MODELPARAMETER(Lambda5);
   const auto Lambda7 = MODELPARAMETER(Lambda7);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,0.5)*(v2*ZA(gt1,0) - v1*ZA(gt1,1))*(ZP(gt2,0)*((Lambda6 - Conj(Lambda6))*ZP(gt3,0) + (-Lambda4 + Conj(Lambda5))*ZP(gt3,1)) + ZP(gt2,1)*((Lambda4 - Lambda5)*ZP(gt3,0) + (Lambda7 - Conj(Lambda7))*ZP(gt3,1)));

   return {result};
}

cxx_diagrams::ChiralVertex VertexImpl<THDMII_cxx_diagrams::fields::Fe, typename THDMII_cxx_diagrams::fields::bar<THDMII_cxx_diagrams::fields::Fe>::type, THDMII_cxx_diagrams::fields::VP>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<THDMII_cxx_diagrams::fields::Fe, typename THDMII_cxx_diagrams::fields::bar<THDMII_cxx_diagrams::fields::Fe>::type, THDMII_cxx_diagrams::fields::VZ>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<THDMII_cxx_diagrams::fields::hh, typename THDMII_cxx_diagrams::fields::bar<THDMII_cxx_diagrams::fields::Fd>::type, THDMII_cxx_diagrams::fields::Fd>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Yd = MODELPARAMETER(Yd);
   const auto Vd = MODELPARAMETER(Vd);
   const auto Ud = MODELPARAMETER(Ud);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(Vd(gt2,j2))*SUM(j1,0,2,Conj(Ud(gt1,j1))*Yd(j1,j2)))*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Ud(gt2,j1))*Vd(gt1,j2))*ZH(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<THDMII_cxx_diagrams::fields::hh, typename THDMII_cxx_diagrams::fields::bar<THDMII_cxx_diagrams::fields::Fe>::type, THDMII_cxx_diagrams::fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ve = MODELPARAMETER(Ve);
   const auto Ue = MODELPARAMETER(Ue);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(Ve(gt2,j2))*SUM(j1,0,2,Conj(Ue(gt1,j1))*Ye(j1,j2)))*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gt2,j1))*Ve(gt1,j2))*ZH(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<THDMII_cxx_diagrams::fields::hh, typename THDMII_cxx_diagrams::fields::bar<THDMII_cxx_diagrams::fields::Fu>::type, THDMII_cxx_diagrams::fields::Fu>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Yu = MODELPARAMETER(Yu);
   const auto Vu = MODELPARAMETER(Vu);
   const auto Uu = MODELPARAMETER(Uu);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = 0.7071067811865475*SUM(j2,0,2,Conj(Vu(gt2,j2))*SUM(j1,0,2,Conj(Uu(gt1,j1))*Yu(j1,j2)))*ZH(gt3,1);

   const std::complex<double> right = 0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Uu(gt2,j1))*Vu(gt1,j2))*ZH(gt3,1);

   return {left, right};
}

cxx_diagrams::ScalarVertex VertexImpl<THDMII_cxx_diagrams::fields::hh, typename THDMII_cxx_diagrams::fields::conj<THDMII_cxx_diagrams::fields::Hm>::type, THDMII_cxx_diagrams::fields::Hm>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Lambda6 = MODELPARAMETER(Lambda6);
   const auto v1 = MODELPARAMETER(v1);
   const auto Lambda3 = MODELPARAMETER(Lambda3);
   const auto v2 = MODELPARAMETER(v2);
   const auto Lambda4 = MODELPARAMETER(Lambda4);
   const auto Lambda7 = MODELPARAMETER(Lambda7);
   const auto Lambda5 = MODELPARAMETER(Lambda5);
   const auto Lambda2 = MODELPARAMETER(Lambda2);
   const auto Lambda1 = MODELPARAMETER(Lambda1);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.5*(-(ZH(gt1,1)*(ZP(gt2,0)*((Lambda6*v1 + 2*Lambda3*v2 + v1*Conj(Lambda6))*ZP(gt3,0) + (Lambda4*v1 + 2*Lambda7*v2 + v1*Conj(Lambda5))*ZP(gt3,1)) + ZP(gt2,1)*(((Lambda4 + Lambda5)*v1 + 2*v2*Conj(Lambda7))*ZP(gt3,0) + (Lambda7*v1 + 4*Lambda2*v2 + v1*Conj(Lambda7))*ZP(gt3,1)))) - ZH(gt1,0)*(ZP(gt2,0)*((4*Lambda1*v1 + Lambda6*v2 + v2*Conj(Lambda6))*ZP(gt3,0) + (2*Lambda6*v1 + Lambda4*v2 + v2*Conj(Lambda5))*ZP(gt3,1)) + ZP(gt2,1)*(((Lambda4 + Lambda5)*v2 + 2*v1*Conj(Lambda6))*ZP(gt3,0) + (2*Lambda3*v1 + Lambda7*v2 + v2*Conj(Lambda7))*ZP(gt3,1))));

   return {result};
}

cxx_diagrams::InverseMetricVertex VertexImpl<THDMII_cxx_diagrams::fields::hh, typename THDMII_cxx_diagrams::fields::conj<THDMII_cxx_diagrams::fields::VWm>::type, THDMII_cxx_diagrams::fields::VWm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.5*Sqr(g2)*(v1*ZH(gt1,0) + v2*ZH(gt1,1));

   return {result};
}

cxx_diagrams::ChiralVertex VertexImpl<THDMII_cxx_diagrams::fields::VP, THDMII_cxx_diagrams::fields::Fd, typename THDMII_cxx_diagrams::fields::bar<THDMII_cxx_diagrams::fields::Fd>::type>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<THDMII_cxx_diagrams::fields::VP, THDMII_cxx_diagrams::fields::Fe, typename THDMII_cxx_diagrams::fields::bar<THDMII_cxx_diagrams::fields::Fe>::type>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<THDMII_cxx_diagrams::fields::VP, THDMII_cxx_diagrams::fields::Fu, typename THDMII_cxx_diagrams::fields::bar<THDMII_cxx_diagrams::fields::Fu>::type>::evaluate(
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

cxx_diagrams::MomentumDifferenceVertex VertexImpl<THDMII_cxx_diagrams::fields::VP, THDMII_cxx_diagrams::fields::Hm, typename THDMII_cxx_diagrams::fields::conj<THDMII_cxx_diagrams::fields::Hm>::type>::evaluate(
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

cxx_diagrams::TripleVectorVertex VertexImpl<THDMII_cxx_diagrams::fields::VP, THDMII_cxx_diagrams::fields::VWm, typename THDMII_cxx_diagrams::fields::conj<THDMII_cxx_diagrams::fields::VWm>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Sin(ThetaW);

   return {result, cxx_diagrams::TripleVectorVertex::odd_permutation{}};
}

cxx_diagrams::ChiralVertex VertexImpl<THDMII_cxx_diagrams::fields::VP, typename THDMII_cxx_diagrams::fields::bar<THDMII_cxx_diagrams::fields::Fe>::type, THDMII_cxx_diagrams::fields::Fe>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<THDMII_cxx_diagrams::fields::VZ, THDMII_cxx_diagrams::fields::Fd, typename THDMII_cxx_diagrams::fields::bar<THDMII_cxx_diagrams::fields::Fd>::type>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<THDMII_cxx_diagrams::fields::VZ, THDMII_cxx_diagrams::fields::Fe, typename THDMII_cxx_diagrams::fields::bar<THDMII_cxx_diagrams::fields::Fe>::type>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<THDMII_cxx_diagrams::fields::VZ, THDMII_cxx_diagrams::fields::Fu, typename THDMII_cxx_diagrams::fields::bar<THDMII_cxx_diagrams::fields::Fu>::type>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<THDMII_cxx_diagrams::fields::VZ, typename THDMII_cxx_diagrams::fields::bar<THDMII_cxx_diagrams::fields::Fe>::type, THDMII_cxx_diagrams::fields::Fe>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<typename THDMII_cxx_diagrams::fields::bar<THDMII_cxx_diagrams::fields::Fe>::type, THDMII_cxx_diagrams::fields::Fe, THDMII_cxx_diagrams::fields::Ah>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ve = MODELPARAMETER(Ve);
   const auto Ue = MODELPARAMETER(Ue);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,Conj(Ve(gt2,j2))*SUM(j1,0,2,Conj(Ue(gt1,j1))*Ye(j1,j2)))*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gt2,j1))*Ve(gt1,j2))*ZA(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename THDMII_cxx_diagrams::fields::bar<THDMII_cxx_diagrams::fields::Fe>::type, THDMII_cxx_diagrams::fields::Fe, THDMII_cxx_diagrams::fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ve = MODELPARAMETER(Ve);
   const auto Ue = MODELPARAMETER(Ue);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(Ve(gt2,j2))*SUM(j1,0,2,Conj(Ue(gt1,j1))*Ye(j1,j2)))*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gt2,j1))*Ve(gt1,j2))*ZH(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename THDMII_cxx_diagrams::fields::bar<THDMII_cxx_diagrams::fields::Fe>::type, THDMII_cxx_diagrams::fields::Fe, THDMII_cxx_diagrams::fields::VP>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<typename THDMII_cxx_diagrams::fields::bar<THDMII_cxx_diagrams::fields::Fe>::type, THDMII_cxx_diagrams::fields::Fe, THDMII_cxx_diagrams::fields::VZ>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<typename THDMII_cxx_diagrams::fields::bar<THDMII_cxx_diagrams::fields::Fe>::type, THDMII_cxx_diagrams::fields::Hm, THDMII_cxx_diagrams::fields::Fv>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ue = MODELPARAMETER(Ue);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> left = -(SUM(j1,0,2,Conj(Ue(gt1,j1))*Ye(j1,gt2))*ZP(gt3,0));

   const std::complex<double> right = 0;

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename THDMII_cxx_diagrams::fields::bar<THDMII_cxx_diagrams::fields::Fe>::type, THDMII_cxx_diagrams::fields::VWm, THDMII_cxx_diagrams::fields::Fv>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ve = MODELPARAMETER(Ve);

   const std::complex<double> left = IF(gt2 < 3,-0.7071067811865475*g2*Ve(gt1,gt2),0);

   const std::complex<double> right = 0;

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename THDMII_cxx_diagrams::fields::bar<THDMII_cxx_diagrams::fields::Fv>::type, typename THDMII_cxx_diagrams::fields::conj<THDMII_cxx_diagrams::fields::Hm>::type, THDMII_cxx_diagrams::fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ue = MODELPARAMETER(Ue);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> left = 0;

   const std::complex<double> right = -(SUM(j1,0,2,Conj(Ye(j1,gt1))*Ue(gt2,j1))*ZP(gt3,0));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename THDMII_cxx_diagrams::fields::bar<THDMII_cxx_diagrams::fields::Fv>::type, typename THDMII_cxx_diagrams::fields::conj<THDMII_cxx_diagrams::fields::VWm>::type, THDMII_cxx_diagrams::fields::Fe>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ve = MODELPARAMETER(Ve);

   const std::complex<double> left = IF(gt1 < 3,-0.7071067811865475*g2*Conj(Ve(gt2,gt1)),0);

   const std::complex<double> right = 0;

   return {left, right};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<typename THDMII_cxx_diagrams::fields::conj<THDMII_cxx_diagrams::fields::Hm>::type, THDMII_cxx_diagrams::fields::Hm, THDMII_cxx_diagrams::fields::VP>::evaluate(
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

cxx_diagrams::InverseMetricVertex VertexImpl<typename THDMII_cxx_diagrams::fields::conj<THDMII_cxx_diagrams::fields::Hm>::type, THDMII_cxx_diagrams::fields::VWm, THDMII_cxx_diagrams::fields::VP>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.3872983346207417*g1*g2*Cos(ThetaW)*(v1*ZP(gt1,0) + v2*ZP(gt1,1));

   return {result};
}

cxx_diagrams::InverseMetricVertex VertexImpl<typename THDMII_cxx_diagrams::fields::conj<THDMII_cxx_diagrams::fields::VWm>::type, THDMII_cxx_diagrams::fields::Hm, THDMII_cxx_diagrams::fields::VP>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.3872983346207417*g1*g2*Cos(ThetaW)*(v1*ZP(gt1,0) + v2*ZP(gt1,1));

   return {result};
}

cxx_diagrams::TripleVectorVertex VertexImpl<typename THDMII_cxx_diagrams::fields::conj<THDMII_cxx_diagrams::fields::VWm>::type, THDMII_cxx_diagrams::fields::VWm, THDMII_cxx_diagrams::fields::VP>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Sin(ThetaW));

   return {result, cxx_diagrams::TripleVectorVertex::odd_permutation{}};
}

} // namespace flexiblesusy::THDMII_cxx_diagrams::detail
