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
 * This file was generated with FlexibleSUSY 2.7.1 and SARAH 4.14.5 .
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

namespace flexiblesusy {
namespace THDMII_cxx_diagrams {
namespace detail {

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Ah, fields::Ah>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambda7 = MODELPARAMETER(Lambda7);
   const auto Lambda2 = MODELPARAMETER(Lambda2);
   const auto Lambda3 = MODELPARAMETER(Lambda3);
   const auto Lambda4 = MODELPARAMETER(Lambda4);
   const auto Lambda5 = MODELPARAMETER(Lambda5);
   const auto Lambda6 = MODELPARAMETER(Lambda6);
   const auto Lambda1 = MODELPARAMETER(Lambda1);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = 0.5*(-(ZA(gt1,1)*(ZA(gt2,1)*(3*ZA(gt3,1)*((Lambda7 + Conj(Lambda7))*ZA(gt4,0) + 4*Lambda2*ZA(gt4,1)) + ZA(gt3,0)*((2*Lambda3 + 2*Lambda4 + Lambda5 + Conj(Lambda5))*ZA(gt4,0) + 3*(Lambda7 + Conj(Lambda7))*ZA(gt4,1))) + ZA(gt2,0)*(ZA(gt3,0)*(3*(Lambda6 + Conj(Lambda6))*ZA(gt4,0) + (2*Lambda3 + 2*Lambda4 + Lambda5 + Conj(Lambda5))*ZA(gt4,1)) + ZA(gt3,1)*((2*Lambda3 + 2*Lambda4 + Lambda5 + Conj(Lambda5))*ZA(gt4,0) + 3*(Lambda7 + Conj(Lambda7))*ZA(gt4,1))))) - ZA(gt1,0)*(ZA(gt2,0)*(ZA(gt3,1)*(3*(Lambda6 + Conj(Lambda6))*ZA(gt4,0) + (2*Lambda3 + 2*Lambda4 + Lambda5 + Conj(Lambda5))*ZA(gt4,1)) + 3*ZA(gt3,0)*(4*Lambda1*ZA(gt4,0) + (Lambda6 + Conj(Lambda6))*ZA(gt4,1))) + ZA(gt2,1)*(ZA(gt3,0)*(3*(Lambda6 + Conj(Lambda6))*ZA(gt4,0) + (2*Lambda3 + 2*Lambda4 + Lambda5 + Conj(Lambda5))*ZA(gt4,1)) + ZA(gt3,1)*((2*Lambda3 + 2*Lambda4 + Lambda5 + Conj(Lambda5))*ZA(gt4,0) + 3*(Lambda7 + Conj(Lambda7))*ZA(gt4,1)))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Ah, fields::hh>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambda7 = MODELPARAMETER(Lambda7);
   const auto Lambda5 = MODELPARAMETER(Lambda5);
   const auto Lambda6 = MODELPARAMETER(Lambda6);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = std::complex<double>(0,-0.5)*(ZA(gt1,1)*(ZA(gt2,1)*(3*(Lambda7 - Conj(Lambda7))*ZA(gt3,1)*ZH(gt4,0) + ZA(gt3,0)*((-Lambda5 + Conj(Lambda5))*ZH(gt4,0) + (-Lambda7 + Conj(Lambda7))*ZH(gt4,1))) + ZA(gt2,0)*(ZA(gt3,0)*((Lambda6 - Conj(Lambda6))*ZH(gt4,0) + (Lambda5 - Conj(Lambda5))*ZH(gt4,1)) + ZA(gt3,1)*((-Lambda5 + Conj(Lambda5))*ZH(gt4,0) + (-Lambda7 + Conj(Lambda7))*ZH(gt4,1)))) + ZA(gt1,0)*(ZA(gt2,0)*(3*(-Lambda6 + Conj(Lambda6))*ZA(gt3,0)*ZH(gt4,1) + ZA(gt3,1)*((Lambda6 - Conj(Lambda6))*ZH(gt4,0) + (Lambda5 - Conj(Lambda5))*ZH(gt4,1))) + ZA(gt2,1)*(ZA(gt3,0)*((Lambda6 - Conj(Lambda6))*ZH(gt4,0) + (Lambda5 - Conj(Lambda5))*ZH(gt4,1)) + ZA(gt3,1)*((-Lambda5 + Conj(Lambda5))*ZH(gt4,0) + (-Lambda7 + Conj(Lambda7))*ZH(gt4,1)))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Ah>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Lambda5 = MODELPARAMETER(Lambda5);
   const auto v1 = MODELPARAMETER(v1);
   const auto Lambda7 = MODELPARAMETER(Lambda7);
   const auto v2 = MODELPARAMETER(v2);
   const auto Lambda6 = MODELPARAMETER(Lambda6);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,0.5)*(ZA(gt1,1)*(ZA(gt2,1)*((Lambda5*v1 + Lambda7*v2 - v1*Conj(Lambda5) - v2*Conj(Lambda7))*ZA(gt3,0) + 3*v1*(-Lambda7 + Conj(Lambda7))*ZA(gt3,1)) + ZA(gt2,0)*((-(Lambda6*v1) - Lambda5*v2 + v2*Conj(Lambda5) + v1*Conj(Lambda6))*ZA(gt3,0) + (Lambda5*v1 + Lambda7*v2 - v1*Conj(Lambda5) - v2*Conj(Lambda7))*ZA(gt3,1))) + ZA(gt1,0)*(ZA(gt2,0)*(3*v2*(Lambda6 - Conj(Lambda6))*ZA(gt3,0) + (-(Lambda6*v1) - Lambda5*v2 + v2*Conj(Lambda5) + v1*Conj(Lambda6))*ZA(gt3,1)) + ZA(gt2,1)*((-(Lambda6*v1) - Lambda5*v2 + v2*Conj(Lambda5) + v1*Conj(Lambda6))*ZA(gt3,0) + (Lambda5*v1 + Lambda7*v2 - v1*Conj(Lambda5) - v2*Conj(Lambda7))*ZA(gt3,1))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::hh, fields::hh>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambda7 = MODELPARAMETER(Lambda7);
   const auto Lambda2 = MODELPARAMETER(Lambda2);
   const auto Lambda3 = MODELPARAMETER(Lambda3);
   const auto Lambda4 = MODELPARAMETER(Lambda4);
   const auto Lambda5 = MODELPARAMETER(Lambda5);
   const auto Lambda6 = MODELPARAMETER(Lambda6);
   const auto Lambda1 = MODELPARAMETER(Lambda1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.5*(-(ZA(gt1,1)*(ZA(gt2,1)*(ZH(gt3,1)*((Lambda7 + Conj(Lambda7))*ZH(gt4,0) + 4*Lambda2*ZH(gt4,1)) + ZH(gt3,0)*((2*Lambda3 + 2*Lambda4 - Lambda5 - Conj(Lambda5))*ZH(gt4,0) + (Lambda7 + Conj(Lambda7))*ZH(gt4,1))) + ZA(gt2,0)*(ZH(gt3,0)*((Lambda6 + Conj(Lambda6))*ZH(gt4,0) + (Lambda5 + Conj(Lambda5))*ZH(gt4,1)) + ZH(gt3,1)*((Lambda5 + Conj(Lambda5))*ZH(gt4,0) + (Lambda7 + Conj(Lambda7))*ZH(gt4,1))))) - ZA(gt1,0)*(ZA(gt2,0)*(ZH(gt3,1)*((Lambda6 + Conj(Lambda6))*ZH(gt4,0) + (2*Lambda3 + 2*Lambda4 - Lambda5 - Conj(Lambda5))*ZH(gt4,1)) + ZH(gt3,0)*(4*Lambda1*ZH(gt4,0) + (Lambda6 + Conj(Lambda6))*ZH(gt4,1))) + ZA(gt2,1)*(ZH(gt3,0)*((Lambda6 + Conj(Lambda6))*ZH(gt4,0) + (Lambda5 + Conj(Lambda5))*ZH(gt4,1)) + ZH(gt3,1)*((Lambda5 + Conj(Lambda5))*ZH(gt4,0) + (Lambda7 + Conj(Lambda7))*ZH(gt4,1)))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Lambda3 = MODELPARAMETER(Lambda3);
   const auto v1 = MODELPARAMETER(v1);
   const auto Lambda4 = MODELPARAMETER(Lambda4);
   const auto Lambda5 = MODELPARAMETER(Lambda5);
   const auto Lambda7 = MODELPARAMETER(Lambda7);
   const auto v2 = MODELPARAMETER(v2);
   const auto Lambda2 = MODELPARAMETER(Lambda2);
   const auto Lambda6 = MODELPARAMETER(Lambda6);
   const auto Lambda1 = MODELPARAMETER(Lambda1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.5*(-(ZA(gt1,1)*(ZA(gt2,1)*((2*Lambda3*v1 + 2*Lambda4*v1 - Lambda5*v1 + Lambda7*v2 - v1*Conj(Lambda5) + v2*Conj(Lambda7))*ZH(gt3,0) + (Lambda7*v1 + 4*Lambda2*v2 + v1*Conj(Lambda7))*ZH(gt3,1)) + ZA(gt2,0)*((Lambda6*v1 + Lambda5*v2 + v2*Conj(Lambda5) + v1*Conj(Lambda6))*ZH(gt3,0) + (Lambda5*v1 + Lambda7*v2 + v1*Conj(Lambda5) + v2*Conj(Lambda7))*ZH(gt3,1)))) - ZA(gt1,0)*(ZA(gt2,0)*((4*Lambda1*v1 + Lambda6*v2 + v2*Conj(Lambda6))*ZH(gt3,0) + (Lambda6*v1 + 2*Lambda3*v2 + 2*Lambda4*v2 - Lambda5*v2 - v2*Conj(Lambda5) + v1*Conj(Lambda6))*ZH(gt3,1)) + ZA(gt2,1)*((Lambda6*v1 + Lambda5*v2 + v2*Conj(Lambda5) + v1*Conj(Lambda6))*ZH(gt3,0) + (Lambda5*v1 + Lambda7*v2 + v1*Conj(Lambda5) + v2*Conj(Lambda7))*ZH(gt3,1))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Hm, typename fields::conj<fields::Hm>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambda7 = MODELPARAMETER(Lambda7);
   const auto Lambda2 = MODELPARAMETER(Lambda2);
   const auto Lambda3 = MODELPARAMETER(Lambda3);
   const auto Lambda6 = MODELPARAMETER(Lambda6);
   const auto Lambda4 = MODELPARAMETER(Lambda4);
   const auto Lambda5 = MODELPARAMETER(Lambda5);
   const auto Lambda1 = MODELPARAMETER(Lambda1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.5*(-(ZA(gt1,1)*(2*ZA(gt2,1)*(ZP(gt3,1)*(Conj(Lambda7)*ZP(gt4,0) + 2*Lambda2*ZP(gt4,1)) + ZP(gt3,0)*(Lambda3*ZP(gt4,0) + Lambda7*ZP(gt4,1))) + ZA(gt2,0)*(ZP(gt3,0)*((Lambda6 + Conj(Lambda6))*ZP(gt4,0) + (Lambda4 + Conj(Lambda5))*ZP(gt4,1)) + ZP(gt3,1)*((Lambda4 + Lambda5)*ZP(gt4,0) + (Lambda7 + Conj(Lambda7))*ZP(gt4,1))))) - ZA(gt1,0)*(2*ZA(gt2,0)*(ZP(gt3,1)*(Conj(Lambda6)*ZP(gt4,0) + Lambda3*ZP(gt4,1)) + ZP(gt3,0)*(2*Lambda1*ZP(gt4,0) + Lambda6*ZP(gt4,1))) + ZA(gt2,1)*(ZP(gt3,0)*((Lambda6 + Conj(Lambda6))*ZP(gt4,0) + (Lambda4 + Conj(Lambda5))*ZP(gt4,1)) + ZP(gt3,1)*((Lambda4 + Lambda5)*ZP(gt4,0) + (Lambda7 + Conj(Lambda7))*ZP(gt4,1)))));

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

ScalarVertex VertexImpl<fields::Ah, fields::hh, fields::hh, fields::hh>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambda7 = MODELPARAMETER(Lambda7);
   const auto Lambda5 = MODELPARAMETER(Lambda5);
   const auto Lambda6 = MODELPARAMETER(Lambda6);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = std::complex<double>(0,-0.5)*(ZA(gt1,1)*(ZH(gt2,1)*((Lambda7 - Conj(Lambda7))*ZH(gt3,1)*ZH(gt4,0) + ZH(gt3,0)*((-Lambda5 + Conj(Lambda5))*ZH(gt4,0) + (Lambda7 - Conj(Lambda7))*ZH(gt4,1))) + ZH(gt2,0)*(ZH(gt3,0)*(3*(Lambda6 - Conj(Lambda6))*ZH(gt4,0) + (-Lambda5 + Conj(Lambda5))*ZH(gt4,1)) + ZH(gt3,1)*((-Lambda5 + Conj(Lambda5))*ZH(gt4,0) + (Lambda7 - Conj(Lambda7))*ZH(gt4,1)))) + ZA(gt1,0)*(ZH(gt2,0)*((-Lambda6 + Conj(Lambda6))*ZH(gt3,0)*ZH(gt4,1) + ZH(gt3,1)*((-Lambda6 + Conj(Lambda6))*ZH(gt4,0) + (Lambda5 - Conj(Lambda5))*ZH(gt4,1))) + ZH(gt2,1)*(ZH(gt3,0)*((-Lambda6 + Conj(Lambda6))*ZH(gt4,0) + (Lambda5 - Conj(Lambda5))*ZH(gt4,1)) + ZH(gt3,1)*((Lambda5 - Conj(Lambda5))*ZH(gt4,0) + 3*(-Lambda7 + Conj(Lambda7))*ZH(gt4,1)))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::hh, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Lambda5 = MODELPARAMETER(Lambda5);
   const auto v1 = MODELPARAMETER(v1);
   const auto Lambda7 = MODELPARAMETER(Lambda7);
   const auto v2 = MODELPARAMETER(v2);
   const auto Lambda6 = MODELPARAMETER(Lambda6);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = std::complex<double>(0,0.5)*(ZA(gt1,0)*(ZH(gt2,0)*(v2*(Lambda6 - Conj(Lambda6))*ZH(gt3,0) + (Lambda6*v1 - Lambda5*v2 + v2*Conj(Lambda5) - v1*Conj(Lambda6))*ZH(gt3,1)) + ZH(gt2,1)*((Lambda6*v1 - Lambda5*v2 + v2*Conj(Lambda5) - v1*Conj(Lambda6))*ZH(gt3,0) + (-(Lambda5*v1) + 3*Lambda7*v2 + v1*Conj(Lambda5) - 3*v2*Conj(Lambda7))*ZH(gt3,1))) + ZA(gt1,1)*(ZH(gt2,1)*((Lambda5*v1 - Lambda7*v2 - v1*Conj(Lambda5) + v2*Conj(Lambda7))*ZH(gt3,0) + v1*(-Lambda7 + Conj(Lambda7))*ZH(gt3,1)) + ZH(gt2,0)*((-3*Lambda6*v1 + Lambda5*v2 - v2*Conj(Lambda5) + 3*v1*Conj(Lambda6))*ZH(gt3,0) + (Lambda5*v1 - Lambda7*v2 - v1*Conj(Lambda5) + v2*Conj(Lambda7))*ZH(gt3,1))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::hh, fields::Hm, typename fields::conj<fields::Hm>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambda6 = MODELPARAMETER(Lambda6);
   const auto Lambda4 = MODELPARAMETER(Lambda4);
   const auto Lambda5 = MODELPARAMETER(Lambda5);
   const auto Lambda7 = MODELPARAMETER(Lambda7);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,-0.5)*(ZA(gt1,1)*ZH(gt2,0) - ZA(gt1,0)*ZH(gt2,1))*(ZP(gt3,0)*((Lambda6 - Conj(Lambda6))*ZP(gt4,0) + (-Lambda4 + Conj(Lambda5))*ZP(gt4,1)) + ZP(gt3,1)*((Lambda4 - Lambda5)*ZP(gt4,0) + (Lambda7 - Conj(Lambda7))*ZP(gt4,1)));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::Ah, fields::hh, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt715 = indices[0];
   const int gt716 = indices[1];

   const std::complex<double> result = 0;

   return {result, minuend_index, subtrahend_index};
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

   const std::complex<double> result = std::complex<double>(0,0.5)*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(ZA(gt1,0)*ZH(gt2,0) + ZA(gt1,1)*ZH(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

ScalarVertex VertexImpl<fields::Ah, fields::Hm, typename fields::conj<fields::Hm>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
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

InverseMetricVertex VertexImpl<fields::Ah, fields::Hm, typename fields::conj<fields::VWm>::type, fields::VP>::evaluate(
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

InverseMetricVertex VertexImpl<fields::Ah, fields::Hm, typename fields::conj<fields::VWm>::type, fields::VZ>::evaluate(
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

MomentumDifferenceVertex VertexImpl<fields::Ah, fields::Hm, typename fields::conj<fields::VWm>::type>::evaluate(
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

ChiralVertex VertexImpl<fields::Ah, typename fields::bar<fields::Fd>::type, fields::Fd>::evaluate(
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

ChiralVertex VertexImpl<fields::Ah, typename fields::bar<fields::Fe>::type, fields::Fe>::evaluate(
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

ChiralVertex VertexImpl<fields::Ah, typename fields::bar<fields::Fu>::type, fields::Fu>::evaluate(
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

ScalarVertex VertexImpl<fields::Ah, typename fields::conj<fields::Hm>::type, fields::Hm>::evaluate(
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

InverseMetricVertex VertexImpl<fields::Ah, typename fields::conj<fields::Hm>::type, fields::VP, fields::VWm>::evaluate(
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

InverseMetricVertex VertexImpl<fields::Ah, typename fields::conj<fields::Hm>::type, fields::VWm, fields::VZ>::evaluate(
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

MomentumDifferenceVertex VertexImpl<fields::Ah, typename fields::conj<fields::Hm>::type, fields::VWm>::evaluate(
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

ChiralVertex VertexImpl<fields::Fe, fields::Ah, typename fields::bar<fields::Fe>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt3 = indices[1];
   const int gt1 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ve = MODELPARAMETER(Ve);
   const auto Ue = MODELPARAMETER(Ue);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,Conj(Ve(gt2,j2))*SUM(j1,0,2,Conj(Ue(gt1,j1))*Ye(j1,j2)))*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gt2,j1))*Ve(gt1,j2))*ZA(gt3,0);

   return {left, right};
}

ChiralVertex VertexImpl<fields::Fe, fields::hh, typename fields::bar<fields::Fe>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt3 = indices[1];
   const int gt1 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ve = MODELPARAMETER(Ve);
   const auto Ue = MODELPARAMETER(Ue);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(Ve(gt2,j2))*SUM(j1,0,2,Conj(Ue(gt1,j1))*Ye(j1,j2)))*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gt2,j1))*Ve(gt1,j2))*ZH(gt3,0);

   return {left, right};
}

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

ChiralVertex VertexImpl<fields::Fe, typename fields::bar<fields::Fv>::type, typename fields::conj<fields::Hm>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const int gt3 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ue = MODELPARAMETER(Ue);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> left = 0;

   const std::complex<double> right = -(SUM(j1,0,2,Conj(Ye(j1,gt1))*Ue(gt2,j1))*ZP(gt3,0));

   return {left, right};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::hh, fields::hh>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambda7 = MODELPARAMETER(Lambda7);
   const auto Lambda2 = MODELPARAMETER(Lambda2);
   const auto Lambda3 = MODELPARAMETER(Lambda3);
   const auto Lambda4 = MODELPARAMETER(Lambda4);
   const auto Lambda5 = MODELPARAMETER(Lambda5);
   const auto Lambda6 = MODELPARAMETER(Lambda6);
   const auto Lambda1 = MODELPARAMETER(Lambda1);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.5*(-(ZH(gt1,1)*(ZH(gt2,1)*(3*ZH(gt3,1)*((Lambda7 + Conj(Lambda7))*ZH(gt4,0) + 4*Lambda2*ZH(gt4,1)) + ZH(gt3,0)*((2*Lambda3 + 2*Lambda4 + Lambda5 + Conj(Lambda5))*ZH(gt4,0) + 3*(Lambda7 + Conj(Lambda7))*ZH(gt4,1))) + ZH(gt2,0)*(ZH(gt3,0)*(3*(Lambda6 + Conj(Lambda6))*ZH(gt4,0) + (2*Lambda3 + 2*Lambda4 + Lambda5 + Conj(Lambda5))*ZH(gt4,1)) + ZH(gt3,1)*((2*Lambda3 + 2*Lambda4 + Lambda5 + Conj(Lambda5))*ZH(gt4,0) + 3*(Lambda7 + Conj(Lambda7))*ZH(gt4,1))))) - ZH(gt1,0)*(ZH(gt2,0)*(ZH(gt3,1)*(3*(Lambda6 + Conj(Lambda6))*ZH(gt4,0) + (2*Lambda3 + 2*Lambda4 + Lambda5 + Conj(Lambda5))*ZH(gt4,1)) + 3*ZH(gt3,0)*(4*Lambda1*ZH(gt4,0) + (Lambda6 + Conj(Lambda6))*ZH(gt4,1))) + ZH(gt2,1)*(ZH(gt3,0)*(3*(Lambda6 + Conj(Lambda6))*ZH(gt4,0) + (2*Lambda3 + 2*Lambda4 + Lambda5 + Conj(Lambda5))*ZH(gt4,1)) + ZH(gt3,1)*((2*Lambda3 + 2*Lambda4 + Lambda5 + Conj(Lambda5))*ZH(gt4,0) + 3*(Lambda7 + Conj(Lambda7))*ZH(gt4,1)))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Lambda3 = MODELPARAMETER(Lambda3);
   const auto v1 = MODELPARAMETER(v1);
   const auto Lambda4 = MODELPARAMETER(Lambda4);
   const auto Lambda5 = MODELPARAMETER(Lambda5);
   const auto Lambda7 = MODELPARAMETER(Lambda7);
   const auto v2 = MODELPARAMETER(v2);
   const auto Lambda2 = MODELPARAMETER(Lambda2);
   const auto Lambda6 = MODELPARAMETER(Lambda6);
   const auto Lambda1 = MODELPARAMETER(Lambda1);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.5*(-(ZH(gt1,1)*(ZH(gt2,1)*((2*Lambda3*v1 + 2*Lambda4*v1 + Lambda5*v1 + 3*Lambda7*v2 + v1*Conj(Lambda5) + 3*v2*Conj(Lambda7))*ZH(gt3,0) + 3*(Lambda7*v1 + 4*Lambda2*v2 + v1*Conj(Lambda7))*ZH(gt3,1)) + ZH(gt2,0)*((3*Lambda6*v1 + 2*Lambda3*v2 + 2*Lambda4*v2 + Lambda5*v2 + v2*Conj(Lambda5) + 3*v1*Conj(Lambda6))*ZH(gt3,0) + (2*Lambda3*v1 + 2*Lambda4*v1 + Lambda5*v1 + 3*Lambda7*v2 + v1*Conj(Lambda5) + 3*v2*Conj(Lambda7))*ZH(gt3,1)))) - ZH(gt1,0)*(ZH(gt2,0)*(3*(4*Lambda1*v1 + Lambda6*v2 + v2*Conj(Lambda6))*ZH(gt3,0) + (3*Lambda6*v1 + 2*Lambda3*v2 + 2*Lambda4*v2 + Lambda5*v2 + v2*Conj(Lambda5) + 3*v1*Conj(Lambda6))*ZH(gt3,1)) + ZH(gt2,1)*((3*Lambda6*v1 + 2*Lambda3*v2 + 2*Lambda4*v2 + Lambda5*v2 + v2*Conj(Lambda5) + 3*v1*Conj(Lambda6))*ZH(gt3,0) + (2*Lambda3*v1 + 2*Lambda4*v1 + Lambda5*v1 + 3*Lambda7*v2 + v1*Conj(Lambda5) + 3*v2*Conj(Lambda7))*ZH(gt3,1))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::Hm, typename fields::conj<fields::Hm>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambda7 = MODELPARAMETER(Lambda7);
   const auto Lambda2 = MODELPARAMETER(Lambda2);
   const auto Lambda3 = MODELPARAMETER(Lambda3);
   const auto Lambda6 = MODELPARAMETER(Lambda6);
   const auto Lambda4 = MODELPARAMETER(Lambda4);
   const auto Lambda5 = MODELPARAMETER(Lambda5);
   const auto Lambda1 = MODELPARAMETER(Lambda1);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.5*(-(ZH(gt1,1)*(2*ZH(gt2,1)*(ZP(gt3,1)*(Conj(Lambda7)*ZP(gt4,0) + 2*Lambda2*ZP(gt4,1)) + ZP(gt3,0)*(Lambda3*ZP(gt4,0) + Lambda7*ZP(gt4,1))) + ZH(gt2,0)*(ZP(gt3,0)*((Lambda6 + Conj(Lambda6))*ZP(gt4,0) + (Lambda4 + Conj(Lambda5))*ZP(gt4,1)) + ZP(gt3,1)*((Lambda4 + Lambda5)*ZP(gt4,0) + (Lambda7 + Conj(Lambda7))*ZP(gt4,1))))) - ZH(gt1,0)*(2*ZH(gt2,0)*(ZP(gt3,1)*(Conj(Lambda6)*ZP(gt4,0) + Lambda3*ZP(gt4,1)) + ZP(gt3,0)*(2*Lambda1*ZP(gt4,0) + Lambda6*ZP(gt4,1))) + ZH(gt2,1)*(ZP(gt3,0)*((Lambda6 + Conj(Lambda6))*ZP(gt4,0) + (Lambda4 + Conj(Lambda5))*ZP(gt4,1)) + ZP(gt3,1)*((Lambda4 + Lambda5)*ZP(gt4,0) + (Lambda7 + Conj(Lambda7))*ZP(gt4,1)))));

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

ScalarVertex VertexImpl<fields::hh, fields::Hm, typename fields::conj<fields::Hm>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
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

InverseMetricVertex VertexImpl<fields::hh, fields::Hm, typename fields::conj<fields::VWm>::type, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.3872983346207417*g1*g2*Cos(ThetaW)*(ZH(gt1,0)*ZP(gt2,0) + ZH(gt1,1)*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::Hm, typename fields::conj<fields::VWm>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.3872983346207417*g1*g2*Sin(ThetaW)*(ZH(gt1,0)*ZP(gt2,0) + ZH(gt1,1)*ZP(gt2,1));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::hh, fields::Hm, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -0.5*g2*(ZH(gt1,0)*ZP(gt2,0) + ZH(gt1,1)*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::hh, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*Sqr(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(v1*ZH(gt1,0) + v2*ZH(gt1,1));

   return {result};
}

ChiralVertex VertexImpl<fields::hh, typename fields::bar<fields::Fd>::type, fields::Fd>::evaluate(
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

ChiralVertex VertexImpl<fields::hh, typename fields::bar<fields::Fe>::type, fields::Fe>::evaluate(
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

ChiralVertex VertexImpl<fields::hh, typename fields::bar<fields::Fu>::type, fields::Fu>::evaluate(
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

ScalarVertex VertexImpl<fields::hh, typename fields::conj<fields::Hm>::type, fields::Hm>::evaluate(
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

InverseMetricVertex VertexImpl<fields::hh, typename fields::conj<fields::Hm>::type, fields::VP, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.3872983346207417*g1*g2*Cos(ThetaW)*(ZH(gt1,0)*ZP(gt2,0) + ZH(gt1,1)*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, typename fields::conj<fields::Hm>::type, fields::VWm, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.3872983346207417*g1*g2*Sin(ThetaW)*(ZH(gt1,0)*ZP(gt2,0) + ZH(gt1,1)*ZP(gt2,1));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::hh, typename fields::conj<fields::Hm>::type, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.5*g2*(ZH(gt1,0)*ZP(gt2,0) + ZH(gt1,1)*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::hh, typename fields::conj<fields::VWm>::type, fields::VWm>::evaluate(
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

ScalarVertex VertexImpl<fields::Hm, fields::Hm, typename fields::conj<fields::Hm>::type, typename fields::conj<fields::Hm>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambda6 = MODELPARAMETER(Lambda6);
   const auto Lambda3 = MODELPARAMETER(Lambda3);
   const auto Lambda4 = MODELPARAMETER(Lambda4);
   const auto Lambda7 = MODELPARAMETER(Lambda7);
   const auto Lambda1 = MODELPARAMETER(Lambda1);
   const auto Lambda5 = MODELPARAMETER(Lambda5);
   const auto Lambda2 = MODELPARAMETER(Lambda2);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = -(ZP(gt1,0)*(ZP(gt2,1)*(2*Conj(Lambda6)*ZP(gt3,0)*ZP(gt4,0) + (Lambda3 + Lambda4)*ZP(gt3,0)*ZP(gt4,1) + ZP(gt3,1)*((Lambda3 + Lambda4)*ZP(gt4,0) + 2*Lambda7*ZP(gt4,1))) + 2*ZP(gt2,0)*(ZP(gt3,0)*(2*Lambda1*ZP(gt4,0) + Lambda6*ZP(gt4,1)) + ZP(gt3,1)*(Lambda6*ZP(gt4,0) + Conj(Lambda5)*ZP(gt4,1))))) - ZP(gt1,1)*(2*Conj(Lambda6)*ZP(gt2,0)*ZP(gt3,0)*ZP(gt4,0) + ZP(gt2,0)*((Lambda3 + Lambda4)*ZP(gt3,0)*ZP(gt4,1) + ZP(gt3,1)*((Lambda3 + Lambda4)*ZP(gt4,0) + 2*Lambda7*ZP(gt4,1))) + 2*ZP(gt2,1)*(ZP(gt3,1)*(Conj(Lambda7)*ZP(gt4,0) + 2*Lambda2*ZP(gt4,1)) + ZP(gt3,0)*(Lambda5*ZP(gt4,0) + Conj(Lambda7)*ZP(gt4,1))));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Hm, typename fields::conj<fields::Hm>::type, fields::VP, fields::VP>::evaluate(
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

InverseMetricVertex VertexImpl<fields::Hm, typename fields::conj<fields::Hm>::type, fields::VP, fields::VZ>::evaluate(
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

MomentumDifferenceVertex VertexImpl<fields::Hm, typename fields::conj<fields::Hm>::type, fields::VP>::evaluate(
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

InverseMetricVertex VertexImpl<fields::Hm, typename fields::conj<fields::Hm>::type, fields::VZ, fields::VZ>::evaluate(
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

MomentumDifferenceVertex VertexImpl<fields::Hm, typename fields::conj<fields::Hm>::type, fields::VZ>::evaluate(
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

InverseMetricVertex VertexImpl<fields::Hm, typename fields::conj<fields::Hm>::type, typename fields::conj<fields::VWm>::type, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.5*Sqr(g2)*(ZP(gt1,0)*ZP(gt2,0) + ZP(gt1,1)*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Hm, typename fields::conj<fields::VWm>::type, fields::VP>::evaluate(
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

InverseMetricVertex VertexImpl<fields::Hm, typename fields::conj<fields::VWm>::type, fields::VZ>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.3872983346207417*g1*g2*Sin(ThetaW)*(v1*ZP(gt1,0) + v2*ZP(gt1,1));

   return {result};
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

MomentumDifferenceVertex VertexImpl<fields::VP, fields::Hm, typename fields::conj<fields::Hm>::type>::evaluate(
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

TripleVectorVertex VertexImpl<fields::VP, fields::VWm, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Sin(ThetaW);

   return {result, TripleVectorVertex::odd_permutation{}};
}

MomentumDifferenceVertex VertexImpl<fields::VWm, fields::Ah, typename fields::conj<fields::Hm>::type>::evaluate(
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

MomentumDifferenceVertex VertexImpl<fields::VWm, fields::hh, typename fields::conj<fields::Hm>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.5*g2*(ZH(gt1,0)*ZP(gt2,0) + ZH(gt1,1)*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::VWm, fields::hh, typename fields::conj<fields::VWm>::type>::evaluate(
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

MomentumDifferenceVertex VertexImpl<fields::VWm, typename fields::conj<fields::Hm>::type, fields::Ah>::evaluate(
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

MomentumDifferenceVertex VertexImpl<fields::VWm, typename fields::conj<fields::Hm>::type, fields::hh>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.5*g2*(ZH(gt1,0)*ZP(gt2,0) + ZH(gt1,1)*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::VWm, typename fields::conj<fields::Hm>::type, fields::VP>::evaluate(
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

InverseMetricVertex VertexImpl<fields::VWm, typename fields::conj<fields::Hm>::type, fields::VZ>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.3872983346207417*g1*g2*Sin(ThetaW)*(v1*ZP(gt1,0) + v2*ZP(gt1,1));

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

   const std::complex<double> result = std::complex<double>(0,0.5)*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(ZA(gt1,0)*ZH(gt2,0) + ZA(gt1,1)*ZH(gt2,1));

   return {result, minuend_index, subtrahend_index};
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

   const std::complex<double> result = std::complex<double>(0,0.5)*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(ZA(gt1,0)*ZH(gt2,0) + ZA(gt1,1)*ZH(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::VZ, fields::hh, fields::VP>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt724 = indices[0];

   const std::complex<double> result = 0;

   return {result};
}

InverseMetricVertex VertexImpl<fields::VZ, fields::hh, fields::VZ>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*Sqr(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(v1*ZH(gt1,0) + v2*ZH(gt1,1));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, fields::Hm, typename fields::conj<fields::Hm>::type>::evaluate(
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

InverseMetricVertex VertexImpl<fields::VZ, fields::Hm, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.3872983346207417*g1*g2*Sin(ThetaW)*(v1*ZP(gt1,0) + v2*ZP(gt1,1));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, typename fields::conj<fields::Hm>::type, fields::Hm>::evaluate(
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

InverseMetricVertex VertexImpl<fields::VZ, typename fields::conj<fields::Hm>::type, fields::VWm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.3872983346207417*g1*g2*Sin(ThetaW)*(v1*ZP(gt1,0) + v2*ZP(gt1,1));

   return {result};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fd, fields::Ah>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Yd = MODELPARAMETER(Yd);
   const auto Vd = MODELPARAMETER(Vd);
   const auto Ud = MODELPARAMETER(Ud);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,Conj(Vd(gt2,j2))*SUM(j1,0,2,Conj(Ud(gt1,j1))*Yd(j1,j2)))*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Ud(gt2,j1))*Vd(gt1,j2))*ZA(gt3,0);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fd, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Yd = MODELPARAMETER(Yd);
   const auto Vd = MODELPARAMETER(Vd);
   const auto Ud = MODELPARAMETER(Ud);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(Vd(gt2,j2))*SUM(j1,0,2,Conj(Ud(gt1,j1))*Yd(j1,j2)))*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Ud(gt2,j1))*Vd(gt1,j2))*ZH(gt3,0);

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

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fu, fields::Hm>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Yd = MODELPARAMETER(Yd);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Vu = MODELPARAMETER(Vu);
   const auto Ud = MODELPARAMETER(Ud);
   const auto ZP = MODELPARAMETER(ZP);
   const auto Uu = MODELPARAMETER(Uu);
   const auto Vd = MODELPARAMETER(Vd);

   const std::complex<double> left = -(SUM(j2,0,2,Conj(Vu(gt2,j2))*SUM(j1,0,2,Conj(Ud(gt1,j1))*Yd(j1,j2)))*ZP(gt3,0));

   const std::complex<double> right = -(SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Uu(gt2,j1))*Vd(gt1,j2))*ZP(gt3,1));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fu, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto Vu = MODELPARAMETER(Vu);
   const auto Vd = MODELPARAMETER(Vd);

   const std::complex<double> left = -0.7071067811865475*g2*SUM(j1,0,2,Conj(Vu(gt2,j1))*Vd(gt1,j1));

   const std::complex<double> right = 0;

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Ah, fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ve = MODELPARAMETER(Ve);
   const auto Ue = MODELPARAMETER(Ue);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,Conj(Ve(gt2,j2))*SUM(j1,0,2,Conj(Ue(gt1,j1))*Ye(j1,j2)))*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gt2,j1))*Ve(gt1,j2))*ZA(gt3,0);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fe, fields::Ah>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fe, fields::hh>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fv, fields::Hm>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ue = MODELPARAMETER(Ue);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> left = -(SUM(j1,0,2,Conj(Ue(gt1,j1))*Ye(j1,gt2))*ZP(gt3,0));

   const std::complex<double> right = 0;

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fv, fields::VWm>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::hh, fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ve = MODELPARAMETER(Ve);
   const auto Ue = MODELPARAMETER(Ue);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(Ve(gt2,j2))*SUM(j1,0,2,Conj(Ue(gt1,j1))*Ye(j1,j2)))*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gt2,j1))*Ve(gt1,j2))*ZH(gt3,0);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Hm, fields::Fv>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fd, typename fields::conj<fields::Hm>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Vd = MODELPARAMETER(Vd);
   const auto Uu = MODELPARAMETER(Uu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto Ud = MODELPARAMETER(Ud);
   const auto Vu = MODELPARAMETER(Vu);

   const std::complex<double> left = -(SUM(j2,0,2,Conj(Vd(gt2,j2))*SUM(j1,0,2,Conj(Uu(gt1,j1))*Yu(j1,j2)))*ZP(gt3,1));

   const std::complex<double> right = -(SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Ud(gt2,j1))*Vu(gt1,j2))*ZP(gt3,0));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fd, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto Vd = MODELPARAMETER(Vd);
   const auto Vu = MODELPARAMETER(Vu);

   const std::complex<double> left = -0.7071067811865475*g2*SUM(j1,0,2,Conj(Vd(gt2,j1))*Vu(gt1,j1));

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
   const auto Vu = MODELPARAMETER(Vu);
   const auto Uu = MODELPARAMETER(Uu);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,Conj(Vu(gt2,j2))*SUM(j1,0,2,Conj(Uu(gt1,j1))*Yu(j1,j2)))*ZA(gt3,1);

   const std::complex<double> right = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Uu(gt2,j1))*Vu(gt1,j2))*ZA(gt3,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fu, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Yu = MODELPARAMETER(Yu);
   const auto Vu = MODELPARAMETER(Vu);
   const auto Uu = MODELPARAMETER(Uu);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = 0.7071067811865475*SUM(j2,0,2,Conj(Vu(gt2,j2))*SUM(j1,0,2,Conj(Uu(gt1,j1))*Yu(j1,j2)))*ZH(gt3,1);

   const std::complex<double> right = 0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Uu(gt2,j1))*Vu(gt1,j2))*ZH(gt3,1);

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

ChiralVertex VertexImpl<typename fields::bar<fields::Fv>::type, fields::Fe, typename fields::conj<fields::Hm>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ue = MODELPARAMETER(Ue);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> left = 0;

   const std::complex<double> right = -(SUM(j1,0,2,Conj(Ye(j1,gt1))*Ue(gt2,j1))*ZP(gt3,0));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fv>::type, fields::Fe, typename fields::conj<fields::VWm>::type>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fv>::type, typename fields::conj<fields::Hm>::type, fields::Fe>::evaluate(
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

ScalarVertex VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gP, typename fields::conj<fields::Hm>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.25*g2*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*(v1*ZP(gt3,0) + v2*ZP(gt3,1));

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
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,-0.25)*Sqr(g2)*(v1*ZA(gt3,0) + v2*ZA(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gWmC, fields::hh>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = -0.25*Sqr(g2)*(v1*ZH(gt3,0) + v2*ZH(gt3,1));

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

ScalarVertex VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gZ, typename fields::conj<fields::Hm>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.25*g2*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*(v1*ZP(gt3,0) + v2*ZP(gt3,1));

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

ScalarVertex VertexImpl<typename fields::bar<fields::gWm>::type, fields::gP, fields::Hm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.25*g2*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*(v1*ZP(gt3,0) + v2*ZP(gt3,1));

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
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,0.25)*Sqr(g2)*(v1*ZA(gt3,0) + v2*ZA(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gWm>::type, fields::gWm, fields::hh>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = -0.25*Sqr(g2)*(v1*ZH(gt3,0) + v2*ZH(gt3,1));

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

ScalarVertex VertexImpl<typename fields::bar<fields::gWm>::type, fields::gZ, fields::Hm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.25*g2*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*(v1*ZP(gt3,0) + v2*ZP(gt3,1));

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

ScalarVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gP, fields::hh>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.025*(7.745966692414834*g1*g2*Cos(2*ThetaW) + Sin(2*ThetaW)*(3*Sqr(g1) - 5*Sqr(g2)))*(v1*ZH(gt3,0) + v2*ZH(gt3,1));

   return {result};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gWmC, fields::Hm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.25*g2*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(v1*ZP(gt3,0) + v2*ZP(gt3,1));

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

ScalarVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gWm, typename fields::conj<fields::Hm>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.25*g2*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(v1*ZP(gt3,0) + v2*ZP(gt3,1));

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
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.25*Sqr(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*(v1*ZH(gt3,0) + v2*ZH(gt3,1));

   return {result};
}

MomentumDifferenceVertex VertexImpl<typename fields::conj<fields::Hm>::type, fields::Hm, fields::VP>::evaluate(
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

InverseMetricVertex VertexImpl<typename fields::conj<fields::Hm>::type, fields::VP, fields::VWm>::evaluate(
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

InverseMetricVertex VertexImpl<typename fields::conj<fields::Hm>::type, fields::VWm, fields::VZ>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto v1 = MODELPARAMETER(v1);
   const auto v2 = MODELPARAMETER(v2);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.3872983346207417*g1*g2*Sin(ThetaW)*(v1*ZP(gt1,0) + v2*ZP(gt1,1));

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
} // namespace THDMII_cxx_diagrams
} // namespace flexiblesusy
