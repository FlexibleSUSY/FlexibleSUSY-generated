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

// File generated at Wed 16 Oct 2019 22:54:17

/**
 * @file cxx_qft/MSSMNoFV_vertices.cpp
 *
 * This file was generated at Wed 16 Oct 2019 22:54:17 with FlexibleSUSY
 * 2.4.1 and SARAH 4.14.3 .
 */

#include "MSSMNoFV_context_base.hpp"
#include "MSSMNoFV_vertices.hpp"

#include "concatenate.hpp"
#include "wrappers.hpp"

#define INPUTPARAMETER(p) context.model.get_input().p
#define MODELPARAMETER(p) context.model.get_##p()
#define DERIVEDPARAMETER(p) context.model.p()
#define PHASE(p) context.model.get_##p()

namespace flexiblesusy {
namespace MSSMNoFV_cxx_diagrams {
namespace detail {

ChiralVertex VertexImpl<fields::Ah, typename bar<fields::Fm>::type, fields::Fm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*Ye(1,1)*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*Conj(Ye(1,1))*ZA(gt3,0);

   return {left, right};
}

ChiralVertex VertexImpl<fields::Chi, typename conj<fields::Sm>::type, fields::Fm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZM = MODELPARAMETER(ZM);

   const std::complex<double> left = 0.5477225575051661*g1*Conj(ZN(gt1,0))*ZM(gt3,0) + 0.7071067811865475*g2*Conj(ZN(gt1,1))*ZM(gt3,0) - Conj(ZN(gt1,2))*Ye(1,1)*ZM(gt3,1);

   const std::complex<double> right = -1.0954451150103321*g1*ZM(gt3,1)*ZN(gt1,0) - Conj(Ye(1,1))*ZM(gt3,0)*ZN(gt1,2);

   return {left, right};
}

ChiralVertex VertexImpl<fields::hh, typename bar<fields::Fm>::type, fields::Fm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*Ye(1,1)*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*Conj(Ye(1,1))*ZH(gt3,0);

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

ChiralVertex VertexImpl<typename bar<fields::Fm>::type, fields::Cha, fields::SvmL>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto UM = MODELPARAMETER(UM);
   const auto UP = MODELPARAMETER(UP);

   const std::complex<double> left = Conj(UM(gt2,1))*Ye(1,1);

   const std::complex<double> right = -(g2*UP(gt2,0));

   return {left, right};
}

ChiralVertex VertexImpl<typename bar<fields::Fm>::type, fields::Fm, fields::Ah>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*Ye(1,1)*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*Conj(Ye(1,1))*ZA(gt3,0);

   return {left, right};
}

ChiralVertex VertexImpl<typename bar<fields::Fm>::type, fields::Fm, fields::hh>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*Ye(1,1)*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*Conj(Ye(1,1))*ZH(gt3,0);

   return {left, right};
}

ChiralVertex VertexImpl<typename bar<fields::Fm>::type, fields::Fm, fields::VP>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.5*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   const std::complex<double> right = 0.7745966692414834*g1*Cos(ThetaW);

   return {left, right};
}

ChiralVertex VertexImpl<typename bar<fields::Fm>::type, fields::Hpm, fields::Fvm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> left = Ye(1,1)*ZP(gt3,0);

   const std::complex<double> right = 0;

   return {left, right};
}

ChiralVertex VertexImpl<typename bar<fields::Fm>::type, fields::Sm, fields::Chi>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZM = MODELPARAMETER(ZM);
   const auto ZN = MODELPARAMETER(ZN);

   const std::complex<double> left = -1.0954451150103321*g1*Conj(ZM(gt3,1))*Conj(ZN(gt2,0)) - Conj(ZM(gt3,0))*Conj(ZN(gt2,2))*Ye(1,1);

   const std::complex<double> right = 0.5*(Conj(ZM(gt3,0))*(1.0954451150103321*g1*ZN(gt2,0) + 1.4142135623730951*g2*ZN(gt2,1)) - 2*Conj(Ye(1,1))*Conj(ZM(gt3,1))*ZN(gt2,2));

   return {left, right};
}

ChiralVertex VertexImpl<typename bar<fields::Fvm>::type, typename conj<fields::Hpm>::type, fields::Fm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> left = 0;

   const std::complex<double> right = Conj(Ye(1,1))*ZP(gt3,0);

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

MomentumDifferenceVertex VertexImpl<typename conj<fields::Sm>::type, fields::Sm, fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 0;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZM = MODELPARAMETER(ZM);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(Conj(ZM(gt1,0))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*ZM(gt2,0) + 1.5491933384829668*g1*Conj(ZM(gt1,1))*Cos(ThetaW)*ZM(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

ChiralVertex VertexImpl<typename conj<fields::SvmL>::type, typename bar<fields::Cha>::type, fields::Fm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto UP = MODELPARAMETER(UP);
   const auto UM = MODELPARAMETER(UM);

   const std::complex<double> left = -(g2*Conj(UP(gt1,0)));

   const std::complex<double> right = Conj(Ye(1,1))*UM(gt1,1);

   return {left, right};
}

} // namespace detail
} // namespace MSSMNoFV_cxx_diagrams
} // namespace flexiblesusy
