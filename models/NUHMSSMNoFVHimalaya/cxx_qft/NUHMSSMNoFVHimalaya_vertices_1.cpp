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
 * @file cxx_qft/NUHMSSMNoFVHimalaya_vertices.cpp
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#include "NUHMSSMNoFVHimalaya_context_base.hpp"
#include "NUHMSSMNoFVHimalaya_input_parameters.hpp"
#include "NUHMSSMNoFVHimalaya_vertices.hpp"

#include "concatenate.hpp"
#include "wrappers.hpp"

#define INPUTPARAMETER(p) context.model.get_input_parameters().p
#define MODELPARAMETER(p) context.model.get_##p()
#define DERIVEDPARAMETER(p) context.model.p()
#define PHASE(p) context.model.get_##p()

namespace flexiblesusy::NUHMSSMNoFVHimalaya_cxx_diagrams::detail {

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Cha>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Cha>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto UM = MODELPARAMETER(UM);
   const auto UP = MODELPARAMETER(UP);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,0.7071067811865475)*g2*(Conj(UM(gt2,1))*Conj(UP(gt1,0))*ZA(gt3,0) + Conj(UM(gt2,0))*Conj(UP(gt1,1))*ZA(gt3,1));

   const std::complex<double> right = std::complex<double>(0.,-0.7071067811865475)*g2*(UM(gt1,1)*UP(gt2,0)*ZA(gt3,0) + UM(gt1,0)*UP(gt2,1)*ZA(gt3,1));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fb>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fb>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*Yd(2,2)*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*Conj(Yd(2,2))*ZA(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fc>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fc>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*Yu(1,1)*ZA(gt3,1);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*Conj(Yu(1,1))*ZA(gt3,1);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fd>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fd>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*Yd(0,0)*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*Conj(Yd(0,0))*ZA(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fe>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fe>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*Ye(0,0)*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*Conj(Ye(0,0))*ZA(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*Ye(1,1)*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*Conj(Ye(1,1))*ZA(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fs>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fs>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*Yd(1,1)*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*Conj(Yd(1,1))*ZA(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ftau>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ftau>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*Ye(2,2)*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*Conj(Ye(2,2))*ZA(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ft>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ft>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*Yu(2,2)*ZA(gt3,1);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*Conj(Yu(2,2))*ZA(gt3,1);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fu>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fu>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*Yu(0,0)*ZA(gt3,1);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*Conj(Yu(0,0))*ZA(gt3,1);

   return {left, right};
}

cxx_diagrams::ScalarVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Hpm>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Hpm>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto vu = MODELPARAMETER(vu);
   const auto vd = MODELPARAMETER(vd);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,-0.25)*Sqr(g2)*(vu*ZA(gt1,0) + vd*ZA(gt1,1))*(ZP(gt2,1)*ZP(gt3,0) - ZP(gt2,0)*ZP(gt3,1));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sb>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sb>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TYd = MODELPARAMETER(TYd);
   const auto Mu = MODELPARAMETER(Mu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZB = MODELPARAMETER(ZB);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)*(Conj(ZB(gt2,1))*(Conj(TYd(2,2))*ZA(gt1,0) + Conj(Yd(2,2))*Mu*ZA(gt1,1))*ZB(gt3,0) - Conj(ZB(gt2,0))*ZB(gt3,1)*(Conj(Mu)*Yd(2,2)*ZA(gt1,1) + ZA(gt1,0)*TYd(2,2)));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sc>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sc>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TYu = MODELPARAMETER(TYu);
   const auto Mu = MODELPARAMETER(Mu);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZC = MODELPARAMETER(ZC);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)*(Conj(Yu(1,1))*Conj(ZC(gt2,1))*Mu*ZA(gt1,0)*ZC(gt3,0) + Conj(ZC(gt2,1))*Conj(TYu(1,1))*ZA(gt1,1)*ZC(gt3,0) - Conj(ZC(gt2,0))*ZC(gt3,1)*(Conj(Mu)*Yu(1,1)*ZA(gt1,0) + ZA(gt1,1)*TYu(1,1)));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sd>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sd>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TYd = MODELPARAMETER(TYd);
   const auto Mu = MODELPARAMETER(Mu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)*(Conj(ZD(gt2,1))*(Conj(TYd(0,0))*ZA(gt1,0) + Conj(Yd(0,0))*Mu*ZA(gt1,1))*ZD(gt3,0) - Conj(ZD(gt2,0))*ZD(gt3,1)*(Conj(Mu)*Yd(0,0)*ZA(gt1,1) + ZA(gt1,0)*TYd(0,0)));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Se>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Se>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TYe = MODELPARAMETER(TYe);
   const auto Mu = MODELPARAMETER(Mu);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)*(Conj(ZE(gt2,1))*(Conj(TYe(0,0))*ZA(gt1,0) + Conj(Ye(0,0))*Mu*ZA(gt1,1))*ZE(gt3,0) - Conj(ZE(gt2,0))*ZE(gt3,1)*(Conj(Mu)*Ye(0,0)*ZA(gt1,1) + ZA(gt1,0)*TYe(0,0)));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sm>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sm>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TYe = MODELPARAMETER(TYe);
   const auto Mu = MODELPARAMETER(Mu);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZM = MODELPARAMETER(ZM);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)*(Conj(ZM(gt2,1))*(Conj(TYe(1,1))*ZA(gt1,0) + Conj(Ye(1,1))*Mu*ZA(gt1,1))*ZM(gt3,0) - Conj(ZM(gt2,0))*ZM(gt3,1)*(Conj(Mu)*Ye(1,1)*ZA(gt1,1) + ZA(gt1,0)*TYe(1,1)));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ss>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ss>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TYd = MODELPARAMETER(TYd);
   const auto Mu = MODELPARAMETER(Mu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZS = MODELPARAMETER(ZS);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)*(Conj(ZS(gt2,1))*(Conj(TYd(1,1))*ZA(gt1,0) + Conj(Yd(1,1))*Mu*ZA(gt1,1))*ZS(gt3,0) - Conj(ZS(gt2,0))*ZS(gt3,1)*(Conj(Mu)*Yd(1,1)*ZA(gt1,1) + ZA(gt1,0)*TYd(1,1)));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Stau>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Stau>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TYe = MODELPARAMETER(TYe);
   const auto Mu = MODELPARAMETER(Mu);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZTau = MODELPARAMETER(ZTau);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)*(Conj(ZTau(gt2,1))*(Conj(TYe(2,2))*ZA(gt1,0) + Conj(Ye(2,2))*Mu*ZA(gt1,1))*ZTau(gt3,0) - Conj(ZTau(gt2,0))*ZTau(gt3,1)*(Conj(Mu)*Ye(2,2)*ZA(gt1,1) + ZA(gt1,0)*TYe(2,2)));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::St>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::St>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TYu = MODELPARAMETER(TYu);
   const auto Mu = MODELPARAMETER(Mu);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZT = MODELPARAMETER(ZT);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)*(Conj(Yu(2,2))*Conj(ZT(gt2,1))*Mu*ZA(gt1,0)*ZT(gt3,0) + Conj(ZT(gt2,1))*Conj(TYu(2,2))*ZA(gt1,1)*ZT(gt3,0) - Conj(ZT(gt2,0))*ZT(gt3,1)*(Conj(Mu)*Yu(2,2)*ZA(gt1,0) + ZA(gt1,1)*TYu(2,2)));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Su>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Su>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TYu = MODELPARAMETER(TYu);
   const auto Mu = MODELPARAMETER(Mu);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)*(Conj(Yu(0,0))*Conj(ZU(gt2,1))*Mu*ZA(gt1,0)*ZU(gt3,0) + Conj(ZU(gt2,1))*Conj(TYu(0,0))*ZA(gt1,1)*ZU(gt3,0) - Conj(ZU(gt2,0))*ZU(gt3,1)*(Conj(Mu)*Yu(0,0)*ZA(gt1,0) + ZA(gt1,1)*TYu(0,0)));

   return {result};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Chi, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sm>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.7745966692414834*g1*Cos(ThetaW);

   const std::complex<double> right = 0.1*(-3.872983346207417*g1*Cos(ThetaW) - 5*g2*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.7745966692414834*g1*Sin(ThetaW);

   const std::complex<double> right = 0.1*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Cha>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Cha>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto UM = MODELPARAMETER(UM);
   const auto UP = MODELPARAMETER(UP);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*g2*(Conj(UM(gt2,1))*Conj(UP(gt1,0))*ZH(gt3,0) + Conj(UM(gt2,0))*Conj(UP(gt1,1))*ZH(gt3,1));

   const std::complex<double> right = -0.7071067811865475*g2*(UM(gt1,1)*UP(gt2,0)*ZH(gt3,0) + UM(gt1,0)*UP(gt2,1)*ZH(gt3,1));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fb>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fb>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*Yd(2,2)*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*Conj(Yd(2,2))*ZH(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fc>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fc>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*Yu(1,1)*ZH(gt3,1);

   const std::complex<double> right = -0.7071067811865475*Conj(Yu(1,1))*ZH(gt3,1);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fd>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fd>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*Yd(0,0)*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*Conj(Yd(0,0))*ZH(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fe>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fe>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*Ye(0,0)*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*Conj(Ye(0,0))*ZH(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*Ye(1,1)*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*Conj(Ye(1,1))*ZH(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fs>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fs>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*Yd(1,1)*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*Conj(Yd(1,1))*ZH(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ftau>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ftau>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*Ye(2,2)*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*Conj(Ye(2,2))*ZH(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ft>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ft>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*Yu(2,2)*ZH(gt3,1);

   const std::complex<double> right = -0.7071067811865475*Conj(Yu(2,2))*ZH(gt3,1);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fu>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fu>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*Yu(0,0)*ZH(gt3,1);

   const std::complex<double> right = -0.7071067811865475*Conj(Yu(0,0))*ZH(gt3,1);

   return {left, right};
}

cxx_diagrams::ScalarVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Hpm>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Hpm>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto vu = MODELPARAMETER(vu);
   const auto g1 = MODELPARAMETER(g1);
   const auto vd = MODELPARAMETER(vd);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-(ZH(gt1,0)*(ZP(gt2,0)*(vd*(0.6*Sqr(g1) + Sqr(g2))*ZP(gt3,0) + vu*Sqr(g2)*ZP(gt3,1)) + ZP(gt2,1)*(vu*Sqr(g2)*ZP(gt3,0) + vd*(-0.6*Sqr(g1) + Sqr(g2))*ZP(gt3,1)))) - ZH(gt1,1)*(ZP(gt2,0)*(vu*(-0.6*Sqr(g1) + Sqr(g2))*ZP(gt3,0) + vd*Sqr(g2)*ZP(gt3,1)) + ZP(gt2,1)*(vd*Sqr(g2)*ZP(gt3,0) + vu*(0.6*Sqr(g1) + Sqr(g2))*ZP(gt3,1))));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sb>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sb>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TYd = MODELPARAMETER(TYd);
   const auto Mu = MODELPARAMETER(Mu);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZB = MODELPARAMETER(ZB);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.05*(-2*Conj(ZB(gt2,1))*(7.0710678118654755*Conj(TYd(2,2))*ZB(gt3,0)*ZH(gt1,0) - 7.0710678118654755*Conj(Yd(2,2))*Mu*ZB(gt3,0)*ZH(gt1,1) + ZB(gt3,1)*(-(vd*(-10*AbsSqr(Yd(2,2)) + Sqr(g1))*ZH(gt1,0)) + vu*Sqr(g1)*ZH(gt1,1))) + Conj(ZB(gt2,0))*(ZB(gt3,0)*(vd*(-20*AbsSqr(Yd(2,2)) + Sqr(g1) + 5*Sqr(g2))*ZH(gt1,0) - vu*(Sqr(g1) + 5*Sqr(g2))*ZH(gt1,1)) + 14.142135623730951*ZB(gt3,1)*(Conj(Mu)*Yd(2,2)*ZH(gt1,1) - ZH(gt1,0)*TYd(2,2))));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sc>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sc>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TYu = MODELPARAMETER(TYu);
   const auto Mu = MODELPARAMETER(Mu);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZC = MODELPARAMETER(ZC);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.05*(2*Conj(ZC(gt2,1))*(-7.0710678118654755*Conj(TYu(1,1))*ZC(gt3,0)*ZH(gt1,1) + 2*Sqr(g1)*ZC(gt3,1)*(-(vd*ZH(gt1,0)) + vu*ZH(gt1,1)) + 5*Conj(Yu(1,1))*(1.4142135623730951*Mu*ZC(gt3,0)*ZH(gt1,0) - 2*vu*Yu(1,1)*ZC(gt3,1)*ZH(gt1,1))) + Conj(ZC(gt2,0))*(ZC(gt3,0)*(vd*(Sqr(g1) - 5*Sqr(g2))*ZH(gt1,0) - vu*(20*AbsSqr(Yu(1,1)) + Sqr(g1) - 5*Sqr(g2))*ZH(gt1,1)) + 14.142135623730951*ZC(gt3,1)*(Conj(Mu)*Yu(1,1)*ZH(gt1,0) - ZH(gt1,1)*TYu(1,1))));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sd>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sd>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TYd = MODELPARAMETER(TYd);
   const auto Mu = MODELPARAMETER(Mu);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.05*(-2*Conj(ZD(gt2,1))*(7.0710678118654755*Conj(TYd(0,0))*ZD(gt3,0)*ZH(gt1,0) - 7.0710678118654755*Conj(Yd(0,0))*Mu*ZD(gt3,0)*ZH(gt1,1) + ZD(gt3,1)*(-(vd*(-10*AbsSqr(Yd(0,0)) + Sqr(g1))*ZH(gt1,0)) + vu*Sqr(g1)*ZH(gt1,1))) + Conj(ZD(gt2,0))*(ZD(gt3,0)*(vd*(-20*AbsSqr(Yd(0,0)) + Sqr(g1) + 5*Sqr(g2))*ZH(gt1,0) - vu*(Sqr(g1) + 5*Sqr(g2))*ZH(gt1,1)) + 14.142135623730951*ZD(gt3,1)*(Conj(Mu)*Yd(0,0)*ZH(gt1,1) - ZH(gt1,0)*TYd(0,0))));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Se>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Se>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TYe = MODELPARAMETER(TYe);
   const auto Mu = MODELPARAMETER(Mu);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.05*(-2*Conj(ZE(gt2,1))*(7.0710678118654755*Conj(TYe(0,0))*ZE(gt3,0)*ZH(gt1,0) - 7.0710678118654755*Conj(Ye(0,0))*Mu*ZE(gt3,0)*ZH(gt1,1) + ZE(gt3,1)*((10*vd*AbsSqr(Ye(0,0)) - 3*vd*Sqr(g1))*ZH(gt1,0) + 3*vu*Sqr(g1)*ZH(gt1,1))) + Conj(ZE(gt2,0))*(ZE(gt3,0)*(vd*(-20*AbsSqr(Ye(0,0)) - 3*Sqr(g1) + 5*Sqr(g2))*ZH(gt1,0) + vu*(3*Sqr(g1) - 5*Sqr(g2))*ZH(gt1,1)) + 14.142135623730951*ZE(gt3,1)*(Conj(Mu)*Ye(0,0)*ZH(gt1,1) - ZH(gt1,0)*TYe(0,0))));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sm>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sm>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TYe = MODELPARAMETER(TYe);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vu = MODELPARAMETER(vu);
   const auto Mu = MODELPARAMETER(Mu);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZM = MODELPARAMETER(ZM);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.25*(-2*Conj(ZM(gt2,1))*(1.4142135623730951*Conj(TYe(1,1))*ZH(gt1,0)*ZM(gt3,0) + 0.6*Sqr(g1)*(-(vd*ZH(gt1,0)) + vu*ZH(gt1,1))*ZM(gt3,1) + Conj(Ye(1,1))*(-1.4142135623730951*Mu*ZH(gt1,1)*ZM(gt3,0) + 2*vd*Ye(1,1)*ZH(gt1,0)*ZM(gt3,1))) - Conj(ZM(gt2,0))*(ZH(gt1,1)*(vu*(-0.6*Sqr(g1) + Sqr(g2))*ZM(gt3,0) - 2.8284271247461903*Conj(Mu)*Ye(1,1)*ZM(gt3,1)) + ZH(gt1,0)*(vd*(4*AbsSqr(Ye(1,1)) + 0.6*Sqr(g1) - Sqr(g2))*ZM(gt3,0) + 2.8284271247461903*ZM(gt3,1)*TYe(1,1))));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ss>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ss>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TYd = MODELPARAMETER(TYd);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vu = MODELPARAMETER(vu);
   const auto Mu = MODELPARAMETER(Mu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZS = MODELPARAMETER(ZS);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.05*(-2*Conj(ZS(gt2,1))*(7.0710678118654755*Conj(TYd(1,1))*ZH(gt1,0)*ZS(gt3,0) + Sqr(g1)*(-(vd*ZH(gt1,0)) + vu*ZH(gt1,1))*ZS(gt3,1) - 5*Conj(Yd(1,1))*(1.4142135623730951*Mu*ZH(gt1,1)*ZS(gt3,0) - 2*vd*Yd(1,1)*ZH(gt1,0)*ZS(gt3,1))) + Conj(ZS(gt2,0))*(-(ZH(gt1,1)*(vu*(Sqr(g1) + 5*Sqr(g2))*ZS(gt3,0) - 14.142135623730951*Conj(Mu)*Yd(1,1)*ZS(gt3,1))) + ZH(gt1,0)*(vd*(-20*AbsSqr(Yd(1,1)) + Sqr(g1) + 5*Sqr(g2))*ZS(gt3,0) - 14.142135623730951*ZS(gt3,1)*TYd(1,1))));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Stau>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Stau>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TYe = MODELPARAMETER(TYe);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vu = MODELPARAMETER(vu);
   const auto Mu = MODELPARAMETER(Mu);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZTau = MODELPARAMETER(ZTau);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.25*(-2*Conj(ZTau(gt2,1))*(1.4142135623730951*Conj(TYe(2,2))*ZH(gt1,0)*ZTau(gt3,0) + 0.6*Sqr(g1)*(-(vd*ZH(gt1,0)) + vu*ZH(gt1,1))*ZTau(gt3,1) + Conj(Ye(2,2))*(-1.4142135623730951*Mu*ZH(gt1,1)*ZTau(gt3,0) + 2*vd*Ye(2,2)*ZH(gt1,0)*ZTau(gt3,1))) - Conj(ZTau(gt2,0))*(ZH(gt1,1)*(vu*(-0.6*Sqr(g1) + Sqr(g2))*ZTau(gt3,0) - 2.8284271247461903*Conj(Mu)*Ye(2,2)*ZTau(gt3,1)) + ZH(gt1,0)*(vd*(4*AbsSqr(Ye(2,2)) + 0.6*Sqr(g1) - Sqr(g2))*ZTau(gt3,0) + 2.8284271247461903*ZTau(gt3,1)*TYe(2,2))));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::St>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::St>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TYu = MODELPARAMETER(TYu);
   const auto vu = MODELPARAMETER(vu);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto Mu = MODELPARAMETER(Mu);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZT = MODELPARAMETER(ZT);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.05*(2*Conj(ZT(gt2,1))*(-7.0710678118654755*Conj(TYu(2,2))*ZH(gt1,1)*ZT(gt3,0) + 2*Sqr(g1)*(-(vd*ZH(gt1,0)) + vu*ZH(gt1,1))*ZT(gt3,1) + 5*Conj(Yu(2,2))*(1.4142135623730951*Mu*ZH(gt1,0)*ZT(gt3,0) - 2*vu*Yu(2,2)*ZH(gt1,1)*ZT(gt3,1))) + Conj(ZT(gt2,0))*(ZH(gt1,0)*(vd*(Sqr(g1) - 5*Sqr(g2))*ZT(gt3,0) + 14.142135623730951*Conj(Mu)*Yu(2,2)*ZT(gt3,1)) - ZH(gt1,1)*(vu*(20*AbsSqr(Yu(2,2)) + Sqr(g1) - 5*Sqr(g2))*ZT(gt3,0) + 14.142135623730951*ZT(gt3,1)*TYu(2,2))));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Su>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Su>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto TYu = MODELPARAMETER(TYu);
   const auto vu = MODELPARAMETER(vu);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto Mu = MODELPARAMETER(Mu);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.05*(2*Conj(ZU(gt2,1))*(-7.0710678118654755*Conj(TYu(0,0))*ZH(gt1,1)*ZU(gt3,0) + 2*Sqr(g1)*(-(vd*ZH(gt1,0)) + vu*ZH(gt1,1))*ZU(gt3,1) + 5*Conj(Yu(0,0))*(1.4142135623730951*Mu*ZH(gt1,0)*ZU(gt3,0) - 2*vu*Yu(0,0)*ZH(gt1,1)*ZU(gt3,1))) + Conj(ZU(gt2,0))*(ZH(gt1,0)*(vd*(Sqr(g1) - 5*Sqr(g2))*ZU(gt3,0) + 14.142135623730951*Conj(Mu)*Yu(0,0)*ZU(gt3,1)) - ZH(gt1,1)*(vu*(20*AbsSqr(Yu(0,0)) + Sqr(g1) - 5*Sqr(g2))*ZU(gt3,0) + 14.142135623730951*ZU(gt3,1)*TYu(0,0))));

   return {result};
}

cxx_diagrams::InverseMetricVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VWm>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VWm>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Cha, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Cha>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UP = MODELPARAMETER(UP);
   const auto UM = MODELPARAMETER(UM);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -(g2*Conj(UP(gt1,0))*Sin(ThetaW)*UP(gt2,0)) - 0.5*Conj(UP(gt1,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UP(gt2,1);

   const std::complex<double> right = -(g2*Conj(UM(gt2,0))*Sin(ThetaW)*UM(gt1,0)) - 0.5*Conj(UM(gt2,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UM(gt1,1);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fb, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fb>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.2581988897471611*g1*Cos(ThetaW);

   const std::complex<double> right = 0.16666666666666666*(0.7745966692414834*g1*Cos(ThetaW) - 3*g2*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fc, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fc>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.5163977794943222*g1*Cos(ThetaW);

   const std::complex<double> right = 0.16666666666666666*(0.7745966692414834*g1*Cos(ThetaW) + 3*g2*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fd, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fd>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.2581988897471611*g1*Cos(ThetaW);

   const std::complex<double> right = 0.16666666666666666*(0.7745966692414834*g1*Cos(ThetaW) - 3*g2*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fe, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fe>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.7745966692414834*g1*Cos(ThetaW);

   const std::complex<double> right = 0.1*(-3.872983346207417*g1*Cos(ThetaW) - 5*g2*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.7745966692414834*g1*Cos(ThetaW);

   const std::complex<double> right = 0.1*(-3.872983346207417*g1*Cos(ThetaW) - 5*g2*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fs, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fs>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.2581988897471611*g1*Cos(ThetaW);

   const std::complex<double> right = 0.16666666666666666*(0.7745966692414834*g1*Cos(ThetaW) - 3*g2*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ftau, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ftau>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.7745966692414834*g1*Cos(ThetaW);

   const std::complex<double> right = 0.1*(-3.872983346207417*g1*Cos(ThetaW) - 5*g2*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ft, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ft>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.5163977794943222*g1*Cos(ThetaW);

   const std::complex<double> right = 0.16666666666666666*(0.7745966692414834*g1*Cos(ThetaW) + 3*g2*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fu, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fu>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.5163977794943222*g1*Cos(ThetaW);

   const std::complex<double> right = 0.16666666666666666*(0.7745966692414834*g1*Cos(ThetaW) + 3*g2*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Hpm, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Hpm>::type>::evaluate(
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

cxx_diagrams::MomentumDifferenceVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sb, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sb>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZB = MODELPARAMETER(ZB);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.16666666666666666*Conj(ZB(gt1,0))*(0.7745966692414834*g1*Cos(ThetaW) - 3*g2*Sin(ThetaW))*ZB(gt2,0) + 0.2581988897471611*g1*Conj(ZB(gt1,1))*Cos(ThetaW)*ZB(gt2,1);

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sc, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sc>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZC = MODELPARAMETER(ZC);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.16666666666666666*Conj(ZC(gt1,0))*(0.7745966692414834*g1*Cos(ThetaW) + 3*g2*Sin(ThetaW))*ZC(gt2,0) - 0.5163977794943222*g1*Conj(ZC(gt1,1))*Cos(ThetaW)*ZC(gt2,1);

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sd, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sd>::type>::evaluate(
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

   const std::complex<double> result = -0.16666666666666666*Conj(ZD(gt1,0))*(0.7745966692414834*g1*Cos(ThetaW) - 3*g2*Sin(ThetaW))*ZD(gt2,0) + 0.2581988897471611*g1*Conj(ZD(gt1,1))*Cos(ThetaW)*ZD(gt2,1);

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Se, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Se>::type>::evaluate(
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

   const std::complex<double> result = 0.5*Conj(ZE(gt1,0))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*ZE(gt2,0) + 0.7745966692414834*g1*Conj(ZE(gt1,1))*Cos(ThetaW)*ZE(gt2,1);

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sm, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sm>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZM = MODELPARAMETER(ZM);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*Conj(ZM(gt1,0))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*ZM(gt2,0) + 0.7745966692414834*g1*Conj(ZM(gt1,1))*Cos(ThetaW)*ZM(gt2,1);

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ss, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ss>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZS = MODELPARAMETER(ZS);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.16666666666666666*Conj(ZS(gt1,0))*(0.7745966692414834*g1*Cos(ThetaW) - 3*g2*Sin(ThetaW))*ZS(gt2,0) + 0.2581988897471611*g1*Conj(ZS(gt1,1))*Cos(ThetaW)*ZS(gt2,1);

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Stau, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Stau>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZTau = MODELPARAMETER(ZTau);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*Conj(ZTau(gt1,0))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*ZTau(gt2,0) + 0.7745966692414834*g1*Conj(ZTau(gt1,1))*Cos(ThetaW)*ZTau(gt2,1);

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::St, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::St>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZT = MODELPARAMETER(ZT);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.16666666666666666*Conj(ZT(gt1,0))*(0.7745966692414834*g1*Cos(ThetaW) + 3*g2*Sin(ThetaW))*ZT(gt2,0) - 0.5163977794943222*g1*Conj(ZT(gt1,1))*Cos(ThetaW)*ZT(gt2,1);

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Su, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Su>::type>::evaluate(
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

   const std::complex<double> result = -0.16666666666666666*Conj(ZU(gt1,0))*(0.7745966692414834*g1*Cos(ThetaW) + 3*g2*Sin(ThetaW))*ZU(gt2,0) - 0.5163977794943222*g1*Conj(ZU(gt1,1))*Cos(ThetaW)*ZU(gt2,1);

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::TripleVectorVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VWm, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VWm>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Sin(ThetaW);

   return {result, cxx_diagrams::TripleVectorVertex::odd_permutation{}};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.5*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   const std::complex<double> right = 0.7745966692414834*g1*Cos(ThetaW);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VZ, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Cha, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Cha>::type>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VZ, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fb, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fb>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.2581988897471611*g1*Sin(ThetaW);

   const std::complex<double> right = 0.16666666666666666*(-3*g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VZ, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fc, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fc>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.5163977794943222*g1*Sin(ThetaW);

   const std::complex<double> right = 0.16666666666666666*(3*g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VZ, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fd, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fd>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.2581988897471611*g1*Sin(ThetaW);

   const std::complex<double> right = 0.16666666666666666*(-3*g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VZ, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fe, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fe>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.7745966692414834*g1*Sin(ThetaW);

   const std::complex<double> right = 0.1*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VZ, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.7745966692414834*g1*Sin(ThetaW);

   const std::complex<double> right = 0.1*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VZ, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fs, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fs>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.2581988897471611*g1*Sin(ThetaW);

   const std::complex<double> right = 0.16666666666666666*(-3*g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VZ, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ftau, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ftau>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.7745966692414834*g1*Sin(ThetaW);

   const std::complex<double> right = 0.1*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VZ, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ft, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ft>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.5163977794943222*g1*Sin(ThetaW);

   const std::complex<double> right = 0.16666666666666666*(3*g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VZ, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fu, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fu>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.5163977794943222*g1*Sin(ThetaW);

   const std::complex<double> right = 0.16666666666666666*(3*g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VZ, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.5*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW));

   const std::complex<double> right = -0.7745966692414834*g1*Sin(ThetaW);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Cha>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Cha, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Cha, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::SvmL>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Ah>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*Ye(1,1)*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*Conj(Ye(1,1))*ZA(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::hh>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*Ye(1,1)*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*Conj(Ye(1,1))*ZH(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.5*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   const std::complex<double> right = 0.7745966692414834*g1*Cos(ThetaW);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.5*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW));

   const std::complex<double> right = -0.7745966692414834*g1*Sin(ThetaW);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Hpm, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fvm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> left = Ye(1,1)*ZP(gt3,0);

   const std::complex<double> right = 0;

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sm, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Chi>::evaluate(
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

   const std::complex<double> right = 0.7071067811865475*Conj(ZM(gt3,0))*(0.7745966692414834*g1*ZN(gt2,0) + g2*ZN(gt2,1)) - Conj(Ye(1,1))*Conj(ZM(gt3,1))*ZN(gt2,2);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VWm, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fvm>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> left = -0.7071067811865475*g2;

   const std::complex<double> right = 0;

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fvm>::type, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Hpm>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> left = 0;

   const std::complex<double> right = Conj(Ye(1,1))*ZP(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fvm>::type, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VWm>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> left = -0.7071067811865475*g2;

   const std::complex<double> right = 0;

   return {left, right};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Hpm>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Hpm, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP>::evaluate(
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

cxx_diagrams::InverseMetricVertex VertexImpl<typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Hpm>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VWm, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP>::evaluate(
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

cxx_diagrams::MomentumDifferenceVertex VertexImpl<typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sm>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Sm, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP>::evaluate(
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

   const std::complex<double> result = 0.5*Conj(ZM(gt1,0))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*ZM(gt2,0) + 0.7745966692414834*g1*Conj(ZM(gt1,1))*Cos(ThetaW)*ZM(gt2,1);

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::ChiralVertex VertexImpl<typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::SvmL>::type, typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::bar<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Cha>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::evaluate(
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

cxx_diagrams::InverseMetricVertex VertexImpl<typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VWm>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Hpm, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP>::evaluate(
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

cxx_diagrams::TripleVectorVertex VertexImpl<typename NUHMSSMNoFVHimalaya_cxx_diagrams::fields::conj<NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VWm>::type, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VWm, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Sin(ThetaW));

   return {result, cxx_diagrams::TripleVectorVertex::odd_permutation{}};
}

} // namespace flexiblesusy::NUHMSSMNoFVHimalaya_cxx_diagrams::detail
