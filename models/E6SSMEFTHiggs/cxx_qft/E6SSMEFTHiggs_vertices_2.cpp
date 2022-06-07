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
 * @file cxx_qft/E6SSMEFTHiggs_vertices.cpp
 *
 * This file was generated with FlexibleSUSY 2.7.1 and SARAH 4.14.5 .
 */

#include "E6SSMEFTHiggs_context_base.hpp"
#include "E6SSMEFTHiggs_input_parameters.hpp"
#include "E6SSMEFTHiggs_vertices.hpp"

#include "concatenate.hpp"
#include "wrappers.hpp"

#define INPUTPARAMETER(p) context.model.get_input_parameters().p
#define MODELPARAMETER(p) context.model.get_##p()
#define DERIVEDPARAMETER(p) context.model.p()
#define PHASE(p) context.model.get_##p()

namespace flexiblesusy {
namespace E6SSMEFTHiggs_cxx_diagrams {
namespace detail {

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

ChiralVertex VertexImpl<fields::VWm, fields::ChiI, typename fields::bar<fields::ChaI>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZPI = MODELPARAMETER(ZPI);
   const auto ZNI = MODELPARAMETER(ZNI);
   const auto ZMI = MODELPARAMETER(ZMI);

   const std::complex<double> left = -0.7071067811865475*g2*SUM(j1,0,1,Conj(ZPI(gt1,j1))*ZNI(gt2,2 + j1));

   const std::complex<double> right = 0.7071067811865475*g2*SUM(j1,0,1,Conj(ZNI(gt2,j1))*ZMI(gt1,j1));

   return {left, right};
}

ChiralVertex VertexImpl<fields::VWm, fields::ChiP, typename fields::bar<fields::ChaP>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZNp = MODELPARAMETER(ZNp);

   const std::complex<double> left = -0.7071067811865475*g2*ZNp(gt2,1);

   const std::complex<double> right = 0.7071067811865475*g2*Conj(ZNp(gt2,0));

   return {left, right};
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

MomentumDifferenceVertex VertexImpl<fields::VWm, fields::SHI0, typename fields::conj<fields::SHIp>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto UHIp = MODELPARAMETER(UHIp);

   const std::complex<double> result = 0.7071067811865475*g2*(-SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHIp(gt2,j1)) + SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHIp(gt2,2 + j1)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VWm, fields::SHp0, typename fields::conj<fields::SHpp>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto UHpp = MODELPARAMETER(UHpp);

   const std::complex<double> result = 0.7071067811865475*(-(g2*Conj(UHp0(gt1,0))*UHpp(gt2,0)) + g2*Conj(UHp0(gt1,1))*UHpp(gt2,1));

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

InverseMetricVertex VertexImpl<fields::VWm, fields::VZp, typename fields::conj<fields::Hpm>::type>::evaluate(
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

TripleVectorVertex VertexImpl<fields::VWm, fields::VZp, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = g2*Cos(ThetaW)*Sin(ThetaWp);

   return {result, TripleVectorVertex::odd_permutation{}};
}

ChiralVertex VertexImpl<fields::VWm, typename fields::bar<fields::ChaI>::type, fields::ChiI>::evaluate(
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

ChiralVertex VertexImpl<fields::VWm, typename fields::bar<fields::ChaP>::type, fields::ChiP>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZNp = MODELPARAMETER(ZNp);

   const std::complex<double> left = -0.7071067811865475*g2*Conj(ZNp(gt2,0));

   const std::complex<double> right = 0.7071067811865475*g2*ZNp(gt2,1);

   return {left, right};
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

InverseMetricVertex VertexImpl<fields::VWm, typename fields::conj<fields::Hpm>::type, fields::VZp>::evaluate(
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

InverseMetricVertex VertexImpl<fields::VWm, typename fields::conj<fields::Hpm>::type, fields::VZ>::evaluate(
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

MomentumDifferenceVertex VertexImpl<fields::VWm, typename fields::conj<fields::SHIp>::type, fields::SHI0>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto UHIp = MODELPARAMETER(UHIp);

   const std::complex<double> result = 0.7071067811865475*g2*(-SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHIp(gt2,j1)) + SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHIp(gt2,2 + j1)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VWm, typename fields::conj<fields::SHpp>::type, fields::SHp0>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto UHpp = MODELPARAMETER(UHpp);

   const std::complex<double> result = 0.7071067811865475*(-(g2*Conj(UHp0(gt1,0))*UHpp(gt2,0)) + g2*Conj(UHp0(gt1,1))*UHpp(gt2,1));

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
   const auto gN = MODELPARAMETER(gN);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = std::complex<double>(0,-0.05)*((10*g2*Cos(ThetaW)*Cos(ThetaWp) + 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp))*ZA(gt1,0)*ZH(gt2,0) - 2*(5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*ZA(gt1,1)*ZH(gt2,1) + 15.811388300841898*gN*Sin(ThetaWp)*ZA(gt1,2)*ZH(gt2,2));

   return {result, minuend_index, subtrahend_index};
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

InverseMetricVertex VertexImpl<fields::VZ, fields::hh, fields::VP>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt4191 = indices[0];

   const std::complex<double> result = 0;

   return {result};
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

   const std::complex<double> result = 0.05*(-(vd*(9.486832980505138*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) - 9*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN) + 5*Sin(2*ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 7.348469228349534*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + g1*(3.872983346207417*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) + 3*g1*Sin(2*ThetaWp)*Sqr(Sin(ThetaW)) - 7.348469228349534*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*ZH(gt1,0)) + vu*(6.324555320336759*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) - 3.872983346207417*g1*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) + 2*Sin(2*ThetaWp)*Sqr(gN) - 5*Sin(2*ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 4.898979485566356*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) - 3*Sin(2*ThetaWp)*Sqr(g1)*Sqr(Sin(ThetaW)) - 4.898979485566356*g1*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp)))*ZH(gt1,1) + 25*vs*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN)*ZH(gt1,2));

   return {result};
}

InverseMetricVertex VertexImpl<fields::VZ, fields::hh, fields::VZ>::evaluate(
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

InverseMetricVertex VertexImpl<fields::VZ, fields::Hpm, typename fields::conj<fields::VWm>::type>::evaluate(
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

MomentumDifferenceVertex VertexImpl<fields::VZ, fields::SHI0, typename fields::conj<fields::SHI0>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

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

MomentumDifferenceVertex VertexImpl<fields::VZ, fields::SHIp, typename fields::conj<fields::SHIp>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

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

MomentumDifferenceVertex VertexImpl<fields::VZ, fields::SHp0, typename fields::conj<fields::SHp0>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

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

MomentumDifferenceVertex VertexImpl<fields::VZ, fields::SHpp, typename fields::conj<fields::SHpp>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

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

MomentumDifferenceVertex VertexImpl<fields::VZ, fields::SSI0, typename fields::conj<fields::SSI0>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -0.7905694150420949*gN*KroneckerDelta(gt1,gt2)*Sin(ThetaWp);

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
   const auto gN = MODELPARAMETER(gN);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.016666666666666666*(-((30*g2*Cos(ThetaW)*Cos(ThetaWp) - 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1))) + (30.983866769659336*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1)));

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
   const auto gN = MODELPARAMETER(gN);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -0.1*KroneckerDelta(gt1,gt2)*(5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp));

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::VZ, fields::VZp, fields::hh>::evaluate(
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

   const std::complex<double> result = 0.05*(-(vd*(9.486832980505138*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) - 9*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN) + 5*Sin(2*ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 7.348469228349534*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + g1*(3.872983346207417*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) + 3*g1*Sin(2*ThetaWp)*Sqr(Sin(ThetaW)) - 7.348469228349534*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp))))*ZH(gt1,0)) + vu*(6.324555320336759*g2*gN*Cos(ThetaW)*Cos(2*ThetaWp) - 3.872983346207417*g1*g2*Sin(2*ThetaW)*Sin(2*ThetaWp) + 2*Sin(2*ThetaWp)*Sqr(gN) - 5*Sin(2*ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 4.898979485566356*g1*gN*Sin(ThetaW)*Sqr(Cos(ThetaWp)) - 3*Sin(2*ThetaWp)*Sqr(g1)*Sqr(Sin(ThetaW)) - 4.898979485566356*g1*gN*Sin(ThetaW)*Sqr(Sin(ThetaWp)))*ZH(gt1,1) + 25*vs*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gN)*ZH(gt1,2));

   return {result};
}

ChiralVertex VertexImpl<fields::VZ, typename fields::bar<fields::ChaI>::type, fields::ChaI>::evaluate(
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

ChiralVertex VertexImpl<fields::VZ, typename fields::bar<fields::ChaP>::type, fields::ChaP>::evaluate(
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

ChiralVertex VertexImpl<fields::VZ, typename fields::bar<fields::Cha>::type, fields::Cha>::evaluate(
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

ChiralVertex VertexImpl<fields::VZ, typename fields::bar<fields::FDX>::type, fields::FDX>::evaluate(
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

MomentumDifferenceVertex VertexImpl<fields::VZ, typename fields::conj<fields::Hpm>::type, fields::Hpm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05*((10*g2*Cos(ThetaW)*Cos(ThetaWp) - 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp))*ZP(gt1,0)*ZP(gt2,0) + 2*(5*g2*Cos(ThetaW)*Cos(ThetaWp) - 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) - 3.1622776601683795*gN*Sin(ThetaWp))*ZP(gt1,1)*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::VZ, typename fields::conj<fields::Hpm>::type, fields::VWm>::evaluate(
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

MomentumDifferenceVertex VertexImpl<fields::VZ, typename fields::conj<fields::Sd>::type, fields::Sd>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.016666666666666666*((30*g2*Cos(ThetaW)*Cos(ThetaWp) + 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1)) + 2*(-7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, typename fields::conj<fields::SDX>::type, fields::SDX>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto ZDX = MODELPARAMETER(ZDX);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.016666666666666666*(-2*(7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 9.486832980505138*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZDX(gt1,j1))*ZDX(gt2,j1)) - (15.491933384829668*g1*Cos(ThetaWp)*Sin(ThetaW) + 28.460498941515414*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZDX(gt1,3 + j1))*ZDX(gt2,3 + j1)));

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
   const auto gN = MODELPARAMETER(gN);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05*(2*(5*g2*Cos(ThetaW)*Cos(ThetaWp) - 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) - 3.1622776601683795*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) + (-15.491933384829668*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, typename fields::conj<fields::SHI0>::type, fields::SHI0>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto UHI0 = MODELPARAMETER(UHI0);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05*((-10*g2*Cos(ThetaW)*Cos(ThetaWp) - 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp))*SUM(j1,0,1,Conj(UHI0(gt1,j1))*UHI0(gt2,j1)) - 2*(5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*SUM(j1,0,1,Conj(UHI0(gt1,2 + j1))*UHI0(gt2,2 + j1)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, typename fields::conj<fields::SHIp>::type, fields::SHIp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto UHIp = MODELPARAMETER(UHIp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05*((10*g2*Cos(ThetaW)*Cos(ThetaWp) - 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp))*SUM(j1,0,1,Conj(UHIp(gt1,j1))*UHIp(gt2,j1)) + 2*(5*g2*Cos(ThetaW)*Cos(ThetaWp) - 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) - 3.1622776601683795*gN*Sin(ThetaWp))*SUM(j1,0,1,Conj(UHIp(gt1,2 + j1))*UHIp(gt2,2 + j1)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, typename fields::conj<fields::SHp0>::type, fields::SHp0>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto UHp0 = MODELPARAMETER(UHp0);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -0.1*(5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp))*(Conj(UHp0(gt1,0))*UHp0(gt2,0) + Conj(UHp0(gt1,1))*UHp0(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, typename fields::conj<fields::SHpp>::type, fields::SHpp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gN = MODELPARAMETER(gN);
   const auto UHpp = MODELPARAMETER(UHpp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.1*(5*g2*Cos(ThetaW)*Cos(ThetaWp) - 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) - 3.1622776601683795*gN*Sin(ThetaWp))*(Conj(UHpp(gt1,0))*UHpp(gt2,0) + Conj(UHpp(gt1,1))*UHpp(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, typename fields::conj<fields::SSI0>::type, fields::SSI0>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto gN = MODELPARAMETER(gN);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -0.7905694150420949*gN*KroneckerDelta(gt1,gt2)*Sin(ThetaWp);

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
   const auto gN = MODELPARAMETER(gN);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.016666666666666666*(-((30*g2*Cos(ThetaW)*Cos(ThetaWp) - 7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1))) + (30.983866769659336*g1*Cos(ThetaWp)*Sin(ThetaW) + 9.486832980505138*gN*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1)));

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
   const auto gN = MODELPARAMETER(gN);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -0.1*KroneckerDelta(gt1,gt2)*(5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 3.1622776601683795*gN*Sin(ThetaWp));

   return {result, minuend_index, subtrahend_index};
}

} // namespace detail
} // namespace E6SSMEFTHiggs_cxx_diagrams
} // namespace flexiblesusy
