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
 * @file cxx_qft/HSSUSY_vertices.cpp
 *
 * This file was generated with FlexibleSUSY 2.6.2 and SARAH 4.14.5 .
 */

#include "HSSUSY_context_base.hpp"
#include "HSSUSY_input_parameters.hpp"
#include "HSSUSY_vertices.hpp"

#include "concatenate.hpp"
#include "wrappers.hpp"

#define INPUTPARAMETER(p) context.model.get_input_parameters().p
#define MODELPARAMETER(p) context.model.get_##p()
#define DERIVEDPARAMETER(p) context.model.p()
#define PHASE(p) context.model.get_##p()

namespace flexiblesusy {
namespace HSSUSY_cxx_diagrams {
namespace detail {

ChiralVertex VertexImpl<fields::Ah, typename fields::bar<fields::Fd>::type, fields::Fd>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Yd = MODELPARAMETER(Yd);
   const auto Vd = MODELPARAMETER(Vd);
   const auto Ud = MODELPARAMETER(Ud);

   const std::complex<double> left = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,Conj(Vd(gt2,j2))*SUM(j1,0,2,Conj(Ud(gt1,j1))*Yd(j1,j2)));

   const std::complex<double> right = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Ud(gt2,j1))*Vd(gt1,j2));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Ah, typename fields::bar<fields::Fe>::type, fields::Fe>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ve = MODELPARAMETER(Ve);
   const auto Ue = MODELPARAMETER(Ue);

   const std::complex<double> left = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,Conj(Ve(gt2,j2))*SUM(j1,0,2,Conj(Ue(gt1,j1))*Ye(j1,j2)));

   const std::complex<double> right = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gt2,j1))*Ve(gt1,j2));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Ah, typename fields::bar<fields::Fu>::type, fields::Fu>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Yu = MODELPARAMETER(Yu);
   const auto Vu = MODELPARAMETER(Vu);
   const auto Uu = MODELPARAMETER(Uu);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,Conj(Vu(gt2,j2))*SUM(j1,0,2,Conj(Uu(gt1,j1))*Yu(j1,j2)));

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Uu(gt2,j1))*Vu(gt1,j2));

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

ChiralVertex VertexImpl<fields::hh, typename fields::bar<fields::Fd>::type, fields::Fd>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Yd = MODELPARAMETER(Yd);
   const auto Vd = MODELPARAMETER(Vd);
   const auto Ud = MODELPARAMETER(Ud);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(Vd(gt2,j2))*SUM(j1,0,2,Conj(Ud(gt1,j1))*Yd(j1,j2)));

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Ud(gt2,j1))*Vd(gt1,j2));

   return {left, right};
}

ChiralVertex VertexImpl<fields::hh, typename fields::bar<fields::Fe>::type, fields::Fe>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ve = MODELPARAMETER(Ve);
   const auto Ue = MODELPARAMETER(Ue);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(Ve(gt2,j2))*SUM(j1,0,2,Conj(Ue(gt1,j1))*Ye(j1,j2)));

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gt2,j1))*Ve(gt1,j2));

   return {left, right};
}

ChiralVertex VertexImpl<fields::hh, typename fields::bar<fields::Fu>::type, fields::Fu>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Yu = MODELPARAMETER(Yu);
   const auto Vu = MODELPARAMETER(Vu);
   const auto Uu = MODELPARAMETER(Uu);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(Vu(gt2,j2))*SUM(j1,0,2,Conj(Uu(gt1,j1))*Yu(j1,j2)));

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Uu(gt2,j1))*Vu(gt1,j2));

   return {left, right};
}

ScalarVertex VertexImpl<fields::hh, typename fields::conj<fields::Hp>::type, fields::Hp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto v = MODELPARAMETER(v);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   const std::complex<double> result = -(v*Lambdax);

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, typename fields::conj<fields::VWp>::type, fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);

   const std::complex<double> result = 0.5*v*Sqr(g2);

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

MomentumDifferenceVertex VertexImpl<fields::VP, fields::Hp, typename fields::conj<fields::Hp>::type>::evaluate(
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

TripleVectorVertex VertexImpl<fields::VP, fields::VWp, typename fields::conj<fields::VWp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Sin(ThetaW));

   return {result, TripleVectorVertex::odd_permutation{}};
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fe, fields::Ah>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ve = MODELPARAMETER(Ve);
   const auto Ue = MODELPARAMETER(Ue);

   const std::complex<double> left = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,Conj(Ve(gt2,j2))*SUM(j1,0,2,Conj(Ue(gt1,j1))*Ye(j1,j2)));

   const std::complex<double> right = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gt2,j1))*Ve(gt1,j2));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fe, fields::hh>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ve = MODELPARAMETER(Ve);
   const auto Ue = MODELPARAMETER(Ue);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(Ve(gt2,j2))*SUM(j1,0,2,Conj(Ue(gt1,j1))*Ye(j1,j2)));

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gt2,j1))*Ve(gt1,j2));

   return {left, right};
}

} // namespace detail
} // namespace HSSUSY_cxx_diagrams
} // namespace flexiblesusy
