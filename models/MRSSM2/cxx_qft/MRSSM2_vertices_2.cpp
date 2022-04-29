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

   const std::complex<double> result = std::complex<double>(0,-0.5)*g2*(ZA(gt1,0)*ZP(gt2,0) + ZA(gt1,1)*ZP(gt2,1) + 1.4142135623730951*ZA(gt1,3)*(ZP(gt2,2) - ZP(gt2,3)));

   return {result, minuend_index, subtrahend_index};
}

ChiralVertex VertexImpl<fields::VWm, fields::Cha1, typename fields::bar<fields::Chi>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZN2 = MODELPARAMETER(ZN2);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto ZN1 = MODELPARAMETER(ZN1);

   const std::complex<double> left = -0.5*g2*(2*Conj(ZN2(gt1,1))*UM1(gt2,0) + 1.4142135623730951*Conj(ZN2(gt1,2))*UM1(gt2,1));

   const std::complex<double> right = -(g2*Conj(UP1(gt2,0))*ZN1(gt1,1)) + 0.7071067811865475*g2*Conj(UP1(gt2,1))*ZN1(gt1,2);

   return {left, right};
}

ChiralVertex VertexImpl<fields::VWm, fields::Chi, typename fields::bar<fields::Cha2>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto UP2 = MODELPARAMETER(UP2);
   const auto ZN2 = MODELPARAMETER(ZN2);
   const auto ZN1 = MODELPARAMETER(ZN1);
   const auto UM2 = MODELPARAMETER(UM2);

   const std::complex<double> left = g2*Conj(UP2(gt1,0))*ZN2(gt2,1) - 0.7071067811865475*g2*Conj(UP2(gt1,1))*ZN2(gt2,3);

   const std::complex<double> right = g2*Conj(ZN1(gt2,1))*UM2(gt1,0) + 0.7071067811865475*g2*Conj(ZN1(gt2,3))*UM2(gt1,1);

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

   const std::complex<double> result = -0.5*g2*(ZH(gt1,0)*ZP(gt2,0) - ZH(gt1,1)*ZP(gt2,1) + 1.4142135623730951*ZH(gt1,3)*(ZP(gt2,2) + ZP(gt2,3)));

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::VWm, fields::hh, typename fields::conj<fields::VWm>::type>::evaluate(
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

MomentumDifferenceVertex VertexImpl<fields::VWm, fields::Rh, typename fields::conj<fields::SRum>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.7071067811865475*g2*ZHR(gt1,1);

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VWm, fields::SRdp, typename fields::conj<fields::Rh>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt2 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.7071067811865475*g2*ZHR(gt2,0);

   return {result, minuend_index, subtrahend_index};
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

ChiralVertex VertexImpl<fields::VWm, typename fields::bar<fields::Cha2>::type, fields::Chi>::evaluate(
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

ChiralVertex VertexImpl<fields::VWm, typename fields::bar<fields::Chi>::type, fields::Cha1>::evaluate(
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

   const std::complex<double> result = std::complex<double>(0,-0.5)*g2*(ZA(gt1,0)*ZP(gt2,0) + ZA(gt1,1)*ZP(gt2,1) + 1.4142135623730951*ZA(gt1,3)*(ZP(gt2,2) - ZP(gt2,3)));

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

   const std::complex<double> result = -0.5*g2*(ZH(gt1,0)*ZP(gt2,0) - ZH(gt1,1)*ZP(gt2,1) + 1.4142135623730951*ZH(gt1,3)*(ZP(gt2,2) + ZP(gt2,3)));

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::VWm, typename fields::conj<fields::Hpm>::type, fields::VP>::evaluate(
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

InverseMetricVertex VertexImpl<fields::VWm, typename fields::conj<fields::Hpm>::type, fields::VZ>::evaluate(
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

MomentumDifferenceVertex VertexImpl<fields::VWm, typename fields::conj<fields::Rh>::type, fields::SRdp>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.7071067811865475*g2*ZHR(gt2,0);

   return {result, minuend_index, subtrahend_index};
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

MomentumDifferenceVertex VertexImpl<fields::VWm, typename fields::conj<fields::SRum>::type, fields::Rh>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZHR = MODELPARAMETER(ZHR);

   const std::complex<double> result = -0.7071067811865475*g2*ZHR(gt1,1);

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::VZ, fields::Hpm, typename fields::conj<fields::VWm>::type>::evaluate(
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

InverseMetricVertex VertexImpl<fields::VZ, typename fields::conj<fields::Hpm>::type, fields::VWm>::evaluate(
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

} // namespace detail
} // namespace MRSSM2_cxx_diagrams
} // namespace flexiblesusy
