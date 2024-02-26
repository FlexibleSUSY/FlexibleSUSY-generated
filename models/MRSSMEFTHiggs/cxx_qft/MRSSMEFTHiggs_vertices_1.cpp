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
 * @file cxx_qft/MRSSMEFTHiggs_vertices.cpp
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#include "MRSSMEFTHiggs_context_base.hpp"
#include "MRSSMEFTHiggs_input_parameters.hpp"
#include "MRSSMEFTHiggs_vertices.hpp"

#include "concatenate.hpp"
#include "wrappers.hpp"

#define INPUTPARAMETER(p) context.model.get_input_parameters().p
#define MODELPARAMETER(p) context.model.get_##p()
#define DERIVEDPARAMETER(p) context.model.p()
#define PHASE(p) context.model.get_##p()

namespace flexiblesusy::MRSSMEFTHiggs_cxx_diagrams::detail {

cxx_diagrams::ChiralVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha1>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Cha1>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0,0.5)*(Conj(UM1(gt1,0))*(-1.4142135623730951*LamTD*Conj(UP1(gt2,1))*ZA(gt3,0) + 2*g2*Conj(UP1(gt2,0))*ZA(gt3,3)) + Conj(UM1(gt1,1))*(1.4142135623730951*g2*Conj(UP1(gt2,0))*ZA(gt3,0) + Conj(UP1(gt2,1))*(-1.4142135623730951*LamSD*ZA(gt3,2) + LamTD*ZA(gt3,3))));

   const std::complex<double> right = std::complex<double>(0,-0.5)*(UM1(gt2,0)*(-1.4142135623730951*Conj(LamTD)*UP1(gt1,1)*ZA(gt3,0) + 2*g2*UP1(gt1,0)*ZA(gt3,3)) + UM1(gt2,1)*(1.4142135623730951*g2*UP1(gt1,0)*ZA(gt3,0) + UP1(gt1,1)*(-1.4142135623730951*Conj(LamSD)*ZA(gt3,2) + Conj(LamTD)*ZA(gt3,3))));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha2>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Cha2>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto UM2 = MODELPARAMETER(UM2);
   const auto UP2 = MODELPARAMETER(UP2);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0,0.5)*(g2*Conj(UM2(gt2,0))*(1.4142135623730951*Conj(UP2(gt1,1))*ZA(gt3,1) - 2*Conj(UP2(gt1,0))*ZA(gt3,3)) + Conj(UM2(gt2,1))*(1.4142135623730951*LamTU*Conj(UP2(gt1,0))*ZA(gt3,1) + Conj(UP2(gt1,1))*(1.4142135623730951*LamSU*ZA(gt3,2) + LamTU*ZA(gt3,3))));

   const std::complex<double> right = std::complex<double>(0,-0.5)*(1.4142135623730951*Conj(LamSU)*UM2(gt1,1)*UP2(gt2,1)*ZA(gt3,2) + g2*UM2(gt1,0)*(1.4142135623730951*UP2(gt2,1)*ZA(gt3,1) - 2*UP2(gt2,0)*ZA(gt3,3)) + Conj(LamTU)*UM2(gt1,1)*(1.4142135623730951*UP2(gt2,0)*ZA(gt3,1) + UP2(gt2,1)*ZA(gt3,3)));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fd>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Fd>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto ZDR = MODELPARAMETER(ZDR);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,Conj(ZDL(gt2,j2))*SUM(j1,0,2,Conj(ZDR(gt1,j1))*Yd(j1,j2)))*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*ZDR(gt2,j1))*ZDL(gt1,j2))*ZA(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fu>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Fu>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto ZUR = MODELPARAMETER(ZUR);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,Conj(ZUL(gt2,j2))*SUM(j1,0,2,Conj(ZUR(gt1,j1))*Yu(j1,j2)))*ZA(gt3,1);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*ZUR(gt2,j1))*ZUL(gt1,j2))*ZA(gt3,1);

   return {left, right};
}

cxx_diagrams::ScalarVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto vu = MODELPARAMETER(vu);
   const auto vd = MODELPARAMETER(vd);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto vT = MODELPARAMETER(vT);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto vS = MODELPARAMETER(vS);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto MuD = MODELPARAMETER(MuD);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto MuU = MODELPARAMETER(MuU);
   const auto g1 = MODELPARAMETER(g1);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,0.05)*(ZA(gt1,2)*(10*LamTD*vd*Conj(LamSD)*ZP(gt2,2)*ZP(gt3,0) - 10*LamSD*vd*Conj(LamTD)*ZP(gt2,3)*ZP(gt3,0) - 7.745966692414834*g1*MDBS*ZP(gt2,1)*ZP(gt3,1) + 14.142135623730951*MuU*Conj(LamSU)*ZP(gt2,1)*ZP(gt3,1) + 7.0710678118654755*LamTU*vT*Conj(LamSU)*ZP(gt2,1)*ZP(gt3,1) - 7.0710678118654755*LamSU*vT*Conj(LamTU)*ZP(gt2,1)*ZP(gt3,1) + 7.745966692414834*g1*Conj(MDBS)*ZP(gt2,1)*ZP(gt3,1) - 14.142135623730951*LamSU*Conj(MuU)*ZP(gt2,1)*ZP(gt3,1) + 10*LamTU*vu*Conj(LamSU)*ZP(gt2,2)*ZP(gt3,1) - 10*LamSU*vu*Conj(LamTU)*ZP(gt2,3)*ZP(gt3,1) - 10*LamSU*vu*Conj(LamTU)*ZP(gt2,1)*ZP(gt3,2) + 10*LamTU*vu*Conj(LamSU)*ZP(gt2,1)*ZP(gt3,3) + ZP(gt2,0)*((7.745966692414834*g1*MDBS + 7.0710678118654755*(2*MuD - LamTD*vT)*Conj(LamSD) + 7.0710678118654755*LamSD*vT*Conj(LamTD) - 7.745966692414834*g1*Conj(MDBS) - 14.142135623730951*LamSD*Conj(MuD))*ZP(gt3,0) + 10*vd*(-(LamSD*Conj(LamTD)*ZP(gt3,2)) + LamTD*Conj(LamSD)*ZP(gt3,3)))) - 5*(ZA(gt1,0)*(vu*Sqr(g2)*ZP(gt2,1)*ZP(gt3,0) + (2*LamTD*vS*Conj(LamSD) + 1.4142135623730951*(2*g2*MDWBT - vT*AbsSqr(LamTD) + 2*LamTD*Conj(MuD) + vT*Sqr(g2)))*ZP(gt2,2)*ZP(gt3,0) + 1.4142135623730951*vT*AbsSqr(LamTD)*ZP(gt2,3)*ZP(gt3,0) + 2.8284271247461903*MuD*Conj(LamTD)*ZP(gt2,3)*ZP(gt3,0) + 2*LamSD*vS*Conj(LamTD)*ZP(gt2,3)*ZP(gt3,0) + 2.8284271247461903*g2*Conj(MDWBT)*ZP(gt2,3)*ZP(gt3,0) - 1.4142135623730951*vT*Sqr(g2)*ZP(gt2,3)*ZP(gt3,0) - vu*Sqr(g2)*ZP(gt2,0)*ZP(gt3,1) + 1.4142135623730951*vT*AbsSqr(LamTD)*ZP(gt2,0)*ZP(gt3,2) - 2.8284271247461903*MuD*Conj(LamTD)*ZP(gt2,0)*ZP(gt3,2) - 2*LamSD*vS*Conj(LamTD)*ZP(gt2,0)*ZP(gt3,2) - 2.8284271247461903*g2*Conj(MDWBT)*ZP(gt2,0)*ZP(gt3,2) - 1.4142135623730951*vT*Sqr(g2)*ZP(gt2,0)*ZP(gt3,2) - 2.8284271247461903*g2*MDWBT*ZP(gt2,0)*ZP(gt3,3) - 1.4142135623730951*vT*AbsSqr(LamTD)*ZP(gt2,0)*ZP(gt3,3) - 2*LamTD*vS*Conj(LamSD)*ZP(gt2,0)*ZP(gt3,3) - 2.8284271247461903*LamTD*Conj(MuD)*ZP(gt2,0)*ZP(gt3,3) + 1.4142135623730951*vT*Sqr(g2)*ZP(gt2,0)*ZP(gt3,3)) + ZA(gt1,1)*(-((vd*Sqr(g2)*ZP(gt2,0) + (2*LamTU*vS*Conj(LamSU) + 1.4142135623730951*(2*g2*MDWBT - vT*AbsSqr(LamTU) + 2*LamTU*Conj(MuU) + vT*Sqr(g2)))*ZP(gt2,2) + ((2.8284271247461903*MuU + 2*LamSU*vS + 1.4142135623730951*LamTU*vT)*Conj(LamTU) + 1.4142135623730951*g2*(-(g2*vT) + 2*Conj(MDWBT)))*ZP(gt2,3))*ZP(gt3,1)) + ZP(gt2,1)*(vd*Sqr(g2)*ZP(gt3,0) + ((2.8284271247461903*MuU + 2*LamSU*vS - 1.4142135623730951*LamTU*vT)*Conj(LamTU) + 1.4142135623730951*g2*(g2*vT + 2*Conj(MDWBT)))*ZP(gt3,2) + (2*LamTU*vS*Conj(LamSU) + 1.4142135623730951*(2*g2*MDWBT + vT*AbsSqr(LamTU) + 2*LamTU*Conj(MuU) - vT*Sqr(g2)))*ZP(gt3,3))) + ZA(gt1,3)*(1.4142135623730951*vd*AbsSqr(LamTD)*ZP(gt2,3)*ZP(gt3,0) - 1.4142135623730951*vd*Sqr(g2)*ZP(gt2,3)*ZP(gt3,0) + 2*g2*MDWBT*ZP(gt2,1)*ZP(gt3,1) + 1.4142135623730951*LamTU*vS*Conj(LamSU)*ZP(gt2,1)*ZP(gt3,1) - 2*MuU*Conj(LamTU)*ZP(gt2,1)*ZP(gt3,1) - 1.4142135623730951*LamSU*vS*Conj(LamTU)*ZP(gt2,1)*ZP(gt3,1) - 2*g2*Conj(MDWBT)*ZP(gt2,1)*ZP(gt3,1) + 2*LamTU*Conj(MuU)*ZP(gt2,1)*ZP(gt3,1) + 1.4142135623730951*vu*AbsSqr(LamTU)*ZP(gt2,3)*ZP(gt3,1) - 1.4142135623730951*vu*Sqr(g2)*ZP(gt2,3)*ZP(gt3,1) - 1.4142135623730951*vu*AbsSqr(LamTU)*ZP(gt2,1)*ZP(gt3,2) + 1.4142135623730951*vu*Sqr(g2)*ZP(gt2,1)*ZP(gt3,2) - 4*vT*Sqr(g2)*ZP(gt2,3)*ZP(gt3,2) - 1.4142135623730951*vu*AbsSqr(LamTU)*ZP(gt2,1)*ZP(gt3,3) + 1.4142135623730951*vu*Sqr(g2)*ZP(gt2,1)*ZP(gt3,3) + ZP(gt2,2)*(1.4142135623730951*vd*(AbsSqr(LamTD) - Sqr(g2))*ZP(gt3,0) + 1.4142135623730951*vu*(AbsSqr(LamTU) - Sqr(g2))*ZP(gt3,1) + 4*vT*Sqr(g2)*ZP(gt3,3)) + ZP(gt2,0)*((-2*g2*MDWBT - 1.4142135623730951*LamTD*vS*Conj(LamSD) + 2*MuD*Conj(LamTD) + 1.4142135623730951*LamSD*vS*Conj(LamTD) + 2*g2*Conj(MDWBT) - 2*LamTD*Conj(MuD))*ZP(gt3,0) + 1.4142135623730951*vd*(-AbsSqr(LamTD) + Sqr(g2))*(ZP(gt3,2) + ZP(gt3,3))))));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Sd>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Sd>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Mu = MODELPARAMETER(Mu);
   const auto g1 = MODELPARAMETER(g1);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,-0.16666666666666666)*(4.242640687119286*Conj(Mu)*SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1)))*ZA(gt1,1) - 4.242640687119286*Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1)))*ZD(gt3,j2))*ZA(gt1,1) + 0.7745966692414834*g1*MDBS*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*ZA(gt1,2) - 0.7745966692414834*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*ZA(gt1,2) + 1.5491933384829668*g1*MDBS*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*ZA(gt1,2) - 1.5491933384829668*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*ZA(gt1,2) - 3*g2*MDWBT*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*ZA(gt1,3) + 3*g2*Conj(MDWBT)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*ZA(gt1,3));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Se>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Se>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Mu = MODELPARAMETER(Mu);
   const auto g1 = MODELPARAMETER(g1);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,-0.1)*(7.0710678118654755*Conj(Mu)*SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)))*ZA(gt1,1) - 7.0710678118654755*Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1)))*ZE(gt3,j2))*ZA(gt1,1) - 3.872983346207417*g1*MDBS*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*ZA(gt1,2) + 3.872983346207417*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*ZA(gt1,2) + 7.745966692414834*g1*MDBS*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*ZA(gt1,2) - 7.745966692414834*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*ZA(gt1,2) - 5*g2*MDWBT*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*ZA(gt1,3) + 5*g2*Conj(MDWBT)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*ZA(gt1,3));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type, MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MuD = MODELPARAMETER(MuD);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto vT = MODELPARAMETER(vT);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto vS = MODELPARAMETER(vS);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,-0.25)*((1.5491933384829668*g1*MDBS + 1.4142135623730951*(-2*MuD + LamTD*vT)*Conj(LamSD) - 1.4142135623730951*LamSD*vT*Conj(LamTD) - 1.5491933384829668*g1*Conj(MDBS) + 2.8284271247461903*LamSD*Conj(MuD))*ZA(gt1,2) + (2*g2*MDWBT - 1.4142135623730951*LamTD*vS*Conj(LamSD) + (2*MuD + 1.4142135623730951*LamSD*vS)*Conj(LamTD) - 2*g2*Conj(MDWBT) - 2*LamTD*Conj(MuD))*ZA(gt1,3));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRum>::type, MRSSMEFTHiggs_cxx_diagrams::fields::SRum>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MuU = MODELPARAMETER(MuU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto vT = MODELPARAMETER(vT);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto vS = MODELPARAMETER(vS);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,-0.25)*((-1.5491933384829668*g1*MDBS - 1.4142135623730951*(2*MuU + LamTU*vT)*Conj(LamSU) + 1.4142135623730951*LamSU*vT*Conj(LamTU) + 1.5491933384829668*g1*Conj(MDBS) + 2.8284271247461903*LamSU*Conj(MuU))*ZA(gt1,2) + (-2*g2*MDWBT + 1.4142135623730951*LamTU*vS*Conj(LamSU) - (2*MuU + 1.4142135623730951*LamSU*vS)*Conj(LamTU) + 2*g2*Conj(MDWBT) + 2*LamTU*Conj(MuU))*ZA(gt1,3));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Su>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Su>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Mu = MODELPARAMETER(Mu);
   const auto g1 = MODELPARAMETER(g1);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0,-0.16666666666666666)*(4.242640687119286*Conj(Mu)*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt3,3 + j1)))*ZA(gt1,0) - 4.242640687119286*Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1)))*ZU(gt3,j2))*ZA(gt1,0) + 0.7745966692414834*g1*MDBS*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*ZA(gt1,2) - 0.7745966692414834*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*ZA(gt1,2) - 3.0983866769659336*g1*MDBS*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*ZA(gt1,2) + 3.0983866769659336*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*ZA(gt1,2) + 3*g2*MDWBT*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*ZA(gt1,3) - 3*g2*Conj(MDWBT)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*ZA(gt1,3));

   return {result};
}

cxx_diagrams::ChiralVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::Cha1, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha1>::type, MRSSMEFTHiggs_cxx_diagrams::fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = g2*Conj(UM1(gt1,0))*Sin(ThetaW)*UM1(gt2,0) + 0.5*Conj(UM1(gt1,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UM1(gt2,1);

   const std::complex<double> right = g2*Conj(UP1(gt2,0))*Sin(ThetaW)*UP1(gt1,0) + 0.5*Conj(UP1(gt2,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UP1(gt1,1);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::Chi, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Se>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZN1 = MODELPARAMETER(ZN1);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZN2 = MODELPARAMETER(ZN2);

   const std::complex<double> left = 0.7071067811865475*(0.7745966692414834*g1*Conj(ZN1(gt1,0)) + g2*Conj(ZN1(gt1,1)))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZE(gt3,j1));

   const std::complex<double> right = -(SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZE(gt3,j2))*ZN2(gt1,2));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::Fe, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::type, MRSSMEFTHiggs_cxx_diagrams::fields::VP>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::Fe, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::type, MRSSMEFTHiggs_cxx_diagrams::fields::VZ>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha1>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Cha1>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = 0.5*(-(Conj(UM1(gt1,0))*(1.4142135623730951*LamTD*Conj(UP1(gt2,1))*ZH(gt3,0) + 2*g2*Conj(UP1(gt2,0))*ZH(gt3,3))) - Conj(UM1(gt1,1))*(1.4142135623730951*g2*Conj(UP1(gt2,0))*ZH(gt3,0) + Conj(UP1(gt2,1))*(1.4142135623730951*LamSD*ZH(gt3,2) - LamTD*ZH(gt3,3))));

   const std::complex<double> right = 0.5*(-(UM1(gt2,0)*(1.4142135623730951*Conj(LamTD)*UP1(gt1,1)*ZH(gt3,0) + 2*g2*UP1(gt1,0)*ZH(gt3,3))) - UM1(gt2,1)*(1.4142135623730951*g2*UP1(gt1,0)*ZH(gt3,0) + UP1(gt1,1)*(1.4142135623730951*Conj(LamSD)*ZH(gt3,2) - Conj(LamTD)*ZH(gt3,3))));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha2>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Cha2>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto UM2 = MODELPARAMETER(UM2);
   const auto UP2 = MODELPARAMETER(UP2);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = 0.5*(g2*Conj(UM2(gt2,0))*(-1.4142135623730951*Conj(UP2(gt1,1))*ZH(gt3,1) + 2*Conj(UP2(gt1,0))*ZH(gt3,3)) + Conj(UM2(gt2,1))*(1.4142135623730951*LamTU*Conj(UP2(gt1,0))*ZH(gt3,1) + Conj(UP2(gt1,1))*(1.4142135623730951*LamSU*ZH(gt3,2) + LamTU*ZH(gt3,3))));

   const std::complex<double> right = 0.5*(1.4142135623730951*Conj(LamSU)*UM2(gt1,1)*UP2(gt2,1)*ZH(gt3,2) + UM2(gt1,0)*(-1.4142135623730951*g2*UP2(gt2,1)*ZH(gt3,1) + 2*g2*UP2(gt2,0)*ZH(gt3,3)) + Conj(LamTU)*UM2(gt1,1)*(1.4142135623730951*UP2(gt2,0)*ZH(gt3,1) + UP2(gt2,1)*ZH(gt3,3)));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fd>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Fd>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto ZDR = MODELPARAMETER(ZDR);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(ZDL(gt2,j2))*SUM(j1,0,2,Conj(ZDR(gt1,j1))*Yd(j1,j2)))*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*ZDR(gt2,j1))*ZDL(gt1,j2))*ZH(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fu>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Fu>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto ZUR = MODELPARAMETER(ZUR);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(ZUL(gt2,j2))*SUM(j1,0,2,Conj(ZUR(gt1,j1))*Yu(j1,j2)))*ZH(gt3,1);

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*ZUR(gt2,j1))*ZUL(gt1,j2))*ZH(gt3,1);

   return {left, right};
}

cxx_diagrams::ScalarVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MuD = MODELPARAMETER(MuD);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto vS = MODELPARAMETER(vS);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto vT = MODELPARAMETER(vT);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto vd = MODELPARAMETER(vd);
   const auto MuU = MODELPARAMETER(MuU);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto vu = MODELPARAMETER(vu);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.05*(7.745966692414834*g1*MDBS*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) - 20*vS*AbsSqr(LamSD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) - 14.142135623730951*MuD*Conj(LamSD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) + 7.0710678118654755*LamTD*vT*Conj(LamSD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) + 7.0710678118654755*LamSD*vT*Conj(LamTD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) + 7.745966692414834*g1*Conj(MDBS)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) - 14.142135623730951*LamSD*Conj(MuD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,0) + 10*g2*MDWBT*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) - 10*vT*AbsSqr(LamTD)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) + 7.0710678118654755*LamTD*vS*Conj(LamSD)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) + 10*MuD*Conj(LamTD)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) + 7.0710678118654755*LamSD*vS*Conj(LamTD)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) + 10*g2*Conj(MDWBT)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) + 10*LamTD*Conj(MuD)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,0) - 10*LamTD*vd*Conj(LamSD)*ZH(gt1,2)*ZP(gt2,2)*ZP(gt3,0) + 7.0710678118654755*vd*AbsSqr(LamTD)*ZH(gt1,3)*ZP(gt2,2)*ZP(gt3,0) - 7.0710678118654755*vd*Sqr(g2)*ZH(gt1,3)*ZP(gt2,2)*ZP(gt3,0) - 10*LamSD*vd*Conj(LamTD)*ZH(gt1,2)*ZP(gt2,3)*ZP(gt3,0) - 7.0710678118654755*vd*AbsSqr(LamTD)*ZH(gt1,3)*ZP(gt2,3)*ZP(gt3,0) + 7.0710678118654755*vd*Sqr(g2)*ZH(gt1,3)*ZP(gt2,3)*ZP(gt3,0) - 7.745966692414834*g1*MDBS*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) - 20*vS*AbsSqr(LamSU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) - 14.142135623730951*MuU*Conj(LamSU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) - 7.0710678118654755*LamTU*vT*Conj(LamSU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) - 7.0710678118654755*LamSU*vT*Conj(LamTU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) - 7.745966692414834*g1*Conj(MDBS)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) - 14.142135623730951*LamSU*Conj(MuU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1) - 10*g2*MDWBT*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,1) - 10*vT*AbsSqr(LamTU)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,1) - 7.0710678118654755*LamTU*vS*Conj(LamSU)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,1) - 10*MuU*Conj(LamTU)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,1) - 7.0710678118654755*LamSU*vS*Conj(LamTU)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,1) - 10*g2*Conj(MDWBT)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,1) - 10*LamTU*Conj(MuU)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,1) - 10*LamTU*vu*Conj(LamSU)*ZH(gt1,2)*ZP(gt2,2)*ZP(gt3,1) + 7.0710678118654755*vu*AbsSqr(LamTU)*ZH(gt1,3)*ZP(gt2,2)*ZP(gt3,1) - 7.0710678118654755*vu*Sqr(g2)*ZH(gt1,3)*ZP(gt2,2)*ZP(gt3,1) - 10*LamSU*vu*Conj(LamTU)*ZH(gt1,2)*ZP(gt2,3)*ZP(gt3,1) - 7.0710678118654755*vu*AbsSqr(LamTU)*ZH(gt1,3)*ZP(gt2,3)*ZP(gt3,1) + 7.0710678118654755*vu*Sqr(g2)*ZH(gt1,3)*ZP(gt2,3)*ZP(gt3,1) - 10*LamSD*vd*Conj(LamTD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,2) + 7.0710678118654755*vd*AbsSqr(LamTD)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,2) - 7.0710678118654755*vd*Sqr(g2)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,2) - 10*LamSU*vu*Conj(LamTU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,2) + 7.0710678118654755*vu*AbsSqr(LamTU)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,2) - 7.0710678118654755*vu*Sqr(g2)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,2) - 20*vT*Sqr(g2)*ZH(gt1,3)*ZP(gt2,2)*ZP(gt3,2) + 20*vT*Sqr(g2)*ZH(gt1,3)*ZP(gt2,3)*ZP(gt3,2) - 10*LamTD*vd*Conj(LamSD)*ZH(gt1,2)*ZP(gt2,0)*ZP(gt3,3) - 7.0710678118654755*vd*AbsSqr(LamTD)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,3) + 7.0710678118654755*vd*Sqr(g2)*ZH(gt1,3)*ZP(gt2,0)*ZP(gt3,3) - 10*LamTU*vu*Conj(LamSU)*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,3) - 7.0710678118654755*vu*AbsSqr(LamTU)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,3) + 7.0710678118654755*vu*Sqr(g2)*ZH(gt1,3)*ZP(gt2,1)*ZP(gt3,3) + 20*vT*Sqr(g2)*ZH(gt1,3)*ZP(gt2,2)*ZP(gt3,3) - 20*vT*Sqr(g2)*ZH(gt1,3)*ZP(gt2,3)*ZP(gt3,3) - ZH(gt1,1)*(ZP(gt2,0)*((-3*vu*Sqr(g1) + 5*vu*Sqr(g2))*ZP(gt3,0) + 5*vd*Sqr(g2)*ZP(gt3,1)) + 5*ZP(gt2,2)*((2*LamTU*vS*Conj(LamSU) + 1.4142135623730951*(2*g2*MDWBT - vT*AbsSqr(LamTU) + 2*LamTU*Conj(MuU) + vT*Sqr(g2)))*ZP(gt3,1) + 2*vu*Sqr(g2)*ZP(gt3,2)) + 5*ZP(gt2,3)*(((2.8284271247461903*MuU + 2*LamSU*vS + 1.4142135623730951*LamTU*vT)*Conj(LamTU) + 1.4142135623730951*g2*(-(g2*vT) + 2*Conj(MDWBT)))*ZP(gt3,1) - 2*vu*(-2*AbsSqr(LamTU) + Sqr(g2))*ZP(gt3,3)) + ZP(gt2,1)*(5*vd*Sqr(g2)*ZP(gt3,0) + vu*(3*Sqr(g1) + 5*Sqr(g2))*ZP(gt3,1) + 5*((2.8284271247461903*MuU + 2*LamSU*vS - 1.4142135623730951*LamTU*vT)*Conj(LamTU) + 1.4142135623730951*g2*(g2*vT + 2*Conj(MDWBT)))*ZP(gt3,2) + 5*(2*LamTU*vS*Conj(LamSU) + 1.4142135623730951*(2*g2*MDWBT + vT*AbsSqr(LamTU) + 2*LamTU*Conj(MuU) - vT*Sqr(g2)))*ZP(gt3,3))) - ZH(gt1,0)*(ZP(gt2,1)*(5*vu*Sqr(g2)*ZP(gt3,0) + vd*(-3*Sqr(g1) + 5*Sqr(g2))*ZP(gt3,1)) + 5*(ZP(gt2,2)*((2*LamTD*vS*Conj(LamSD) + 1.4142135623730951*(2*g2*MDWBT - vT*AbsSqr(LamTD) + 2*LamTD*Conj(MuD) + vT*Sqr(g2)))*ZP(gt3,0) - 2*vd*(-2*AbsSqr(LamTD) + Sqr(g2))*ZP(gt3,2)) + ZP(gt2,3)*(((2.8284271247461903*MuD + 2*LamSD*vS + 1.4142135623730951*LamTD*vT)*Conj(LamTD) + 1.4142135623730951*g2*(-(g2*vT) + 2*Conj(MDWBT)))*ZP(gt3,0) + 2*vd*Sqr(g2)*ZP(gt3,3))) + ZP(gt2,0)*(vd*(3*Sqr(g1) + 5*Sqr(g2))*ZP(gt3,0) + 5*(vu*Sqr(g2)*ZP(gt3,1) + ((2.8284271247461903*MuD + 2*LamSD*vS - 1.4142135623730951*LamTD*vT)*Conj(LamTD) + 1.4142135623730951*g2*(g2*vT + 2*Conj(MDWBT)))*ZP(gt3,2) + (2*LamTD*vS*Conj(LamSD) + 1.4142135623730951*(2*g2*MDWBT + vT*AbsSqr(LamTD) + 2*LamTD*Conj(MuD) - vT*Sqr(g2)))*ZP(gt3,3)))));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Sd>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Sd>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto vd = MODELPARAMETER(vd);
   const auto Mu = MODELPARAMETER(Mu);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.016666666666666666*(30*(-2*vd*SUM(j3,0,2,Conj(ZD(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,j1))*ZD(gt3,3 + j2)))*ZH(gt1,0) - 2*vd*SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gt3,j3))*ZH(gt1,0) + 1.4142135623730951*(Conj(Mu)*SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1))) + Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1)))*ZD(gt3,j2)))*ZH(gt1,1)) + 2*g1*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*(3*g1*vd*ZH(gt1,0) - 3*g1*vu*ZH(gt1,1) - 7.745966692414834*(MDBS + Conj(MDBS))*ZH(gt1,2)) + SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*(3*vd*(Sqr(g1) + 5*Sqr(g2))*ZH(gt1,0) - 3*vu*(Sqr(g1) + 5*Sqr(g2))*ZH(gt1,1) - 7.745966692414834*g1*MDBS*ZH(gt1,2) - 7.745966692414834*g1*Conj(MDBS)*ZH(gt1,2) + 30*g2*MDWBT*ZH(gt1,3) + 30*g2*Conj(MDWBT)*ZH(gt1,3)));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Se>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Se>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto vd = MODELPARAMETER(vd);
   const auto Mu = MODELPARAMETER(Mu);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.05*(2*(5*(-2*vd*SUM(j3,0,2,Conj(ZE(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Ye(j2,j1))*ZE(gt3,3 + j2)))*ZH(gt1,0) - 2*vd*SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZE(gt3,j3))*ZH(gt1,0) + 1.4142135623730951*(Conj(Mu)*SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1))) + Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1)))*ZE(gt3,j2)))*ZH(gt1,1)) + g1*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*(3*g1*vd*ZH(gt1,0) - 3*g1*vu*ZH(gt1,1) - 7.745966692414834*(MDBS + Conj(MDBS))*ZH(gt1,2))) + SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*((-3*vd*Sqr(g1) + 5*vd*Sqr(g2))*ZH(gt1,0) + vu*(3*Sqr(g1) - 5*Sqr(g2))*ZH(gt1,1) + 7.745966692414834*g1*(MDBS + Conj(MDBS))*ZH(gt1,2) + 10*g2*(MDWBT + Conj(MDWBT))*ZH(gt1,3)));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type, MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto LamTD = MODELPARAMETER(LamTD);
   const auto vu = MODELPARAMETER(vu);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MuD = MODELPARAMETER(MuD);
   const auto LamSD = MODELPARAMETER(LamSD);
   const auto vS = MODELPARAMETER(vS);
   const auto vT = MODELPARAMETER(vT);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.25*(vd*(-4*AbsSqr(LamTD) + 0.6*Sqr(g1) - Sqr(g2))*ZH(gt1,0) + vu*(-0.6*Sqr(g1) + Sqr(g2))*ZH(gt1,1) - 1.5491933384829668*g1*MDBS*ZH(gt1,2) - 4*vS*AbsSqr(LamSD)*ZH(gt1,2) - 2.8284271247461903*MuD*Conj(LamSD)*ZH(gt1,2) + 1.4142135623730951*LamTD*vT*Conj(LamSD)*ZH(gt1,2) + 1.4142135623730951*LamSD*vT*Conj(LamTD)*ZH(gt1,2) - 1.5491933384829668*g1*Conj(MDBS)*ZH(gt1,2) - 2.8284271247461903*LamSD*Conj(MuD)*ZH(gt1,2) - 2*g2*MDWBT*ZH(gt1,3) - 2*vT*AbsSqr(LamTD)*ZH(gt1,3) + 1.4142135623730951*LamTD*vS*Conj(LamSD)*ZH(gt1,3) + 2*MuD*Conj(LamTD)*ZH(gt1,3) + 1.4142135623730951*LamSD*vS*Conj(LamTD)*ZH(gt1,3) - 2*g2*Conj(MDWBT)*ZH(gt1,3) + 2*LamTD*Conj(MuD)*ZH(gt1,3));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRum>::type, MRSSMEFTHiggs_cxx_diagrams::fields::SRum>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto LamTU = MODELPARAMETER(LamTU);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto MuU = MODELPARAMETER(MuU);
   const auto LamSU = MODELPARAMETER(LamSU);
   const auto vS = MODELPARAMETER(vS);
   const auto vT = MODELPARAMETER(vT);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.25*(-0.2*vd*(3*Sqr(g1) - 5*Sqr(g2))*ZH(gt1,0) - vu*(4*AbsSqr(LamTU) - 0.6*Sqr(g1) + Sqr(g2))*ZH(gt1,1) + 1.5491933384829668*g1*MDBS*ZH(gt1,2) - 4*vS*AbsSqr(LamSU)*ZH(gt1,2) - 2.8284271247461903*MuU*Conj(LamSU)*ZH(gt1,2) - 1.4142135623730951*LamTU*vT*Conj(LamSU)*ZH(gt1,2) - 1.4142135623730951*LamSU*vT*Conj(LamTU)*ZH(gt1,2) + 1.5491933384829668*g1*Conj(MDBS)*ZH(gt1,2) - 2.8284271247461903*LamSU*Conj(MuU)*ZH(gt1,2) + 2*g2*MDWBT*ZH(gt1,3) - 2*vT*AbsSqr(LamTU)*ZH(gt1,3) - 1.4142135623730951*LamTU*vS*Conj(LamSU)*ZH(gt1,3) - 2*MuU*Conj(LamTU)*ZH(gt1,3) - 1.4142135623730951*LamSU*vS*Conj(LamTU)*ZH(gt1,3) + 2*g2*Conj(MDWBT)*ZH(gt1,3) - 2*LamTU*Conj(MuU)*ZH(gt1,3));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Su>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Su>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Mu = MODELPARAMETER(Mu);
   const auto vu = MODELPARAMETER(vu);
   const auto g1 = MODELPARAMETER(g1);
   const auto vd = MODELPARAMETER(vd);
   const auto MDBS = MODELPARAMETER(MDBS);
   const auto g2 = MODELPARAMETER(g2);
   const auto MDWBT = MODELPARAMETER(MDWBT);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.016666666666666666*(30*(1.4142135623730951*Conj(Mu)*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt3,3 + j1)))*ZH(gt1,0) + 1.4142135623730951*Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1)))*ZU(gt3,j2))*ZH(gt1,0) - 2*vu*(SUM(j3,0,2,Conj(ZU(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))*ZU(gt3,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gt3,j3)))*ZH(gt1,1)) + 4*g1*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*(-3*g1*vd*ZH(gt1,0) + 3*g1*vu*ZH(gt1,1) + 7.745966692414834*(MDBS + Conj(MDBS))*ZH(gt1,2)) + SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*(3*vd*(Sqr(g1) - 5*Sqr(g2))*ZH(gt1,0) - 3*vu*(Sqr(g1) - 5*Sqr(g2))*ZH(gt1,1) - 7.745966692414834*g1*(MDBS + Conj(MDBS))*ZH(gt1,2) - 30*g2*(MDWBT + Conj(MDWBT))*ZH(gt1,3)));

   return {result};
}

cxx_diagrams::InverseMetricVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::VWm>::type, MRSSMEFTHiggs_cxx_diagrams::fields::VWm>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::VP, MRSSMEFTHiggs_cxx_diagrams::fields::Cha1, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha1>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = g2*Conj(UM1(gt1,0))*Sin(ThetaW)*UM1(gt2,0) + 0.5*Conj(UM1(gt1,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UM1(gt2,1);

   const std::complex<double> right = g2*Conj(UP1(gt2,0))*Sin(ThetaW)*UP1(gt1,0) + 0.5*Conj(UP1(gt2,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UP1(gt1,1);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::VP, MRSSMEFTHiggs_cxx_diagrams::fields::Cha2, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha2>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UP2 = MODELPARAMETER(UP2);
   const auto UM2 = MODELPARAMETER(UM2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -(g2*Conj(UP2(gt1,0))*Sin(ThetaW)*UP2(gt2,0)) - 0.5*Conj(UP2(gt1,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UP2(gt2,1);

   const std::complex<double> right = -(g2*Conj(UM2(gt2,0))*Sin(ThetaW)*UM2(gt1,0)) - 0.5*Conj(UM2(gt2,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UM2(gt1,1);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::VP, MRSSMEFTHiggs_cxx_diagrams::fields::Fd, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fd>::type>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::VP, MRSSMEFTHiggs_cxx_diagrams::fields::Fe, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::type>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::VP, MRSSMEFTHiggs_cxx_diagrams::fields::Fu, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fu>::type>::evaluate(
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

cxx_diagrams::MomentumDifferenceVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::VP, MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type>::evaluate(
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

   const std::complex<double> result = 0.5*((0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*ZP(gt1,0)*ZP(gt2,0) + (0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*ZP(gt1,1)*ZP(gt2,1) + 2*g2*Sin(ThetaW)*(ZP(gt1,2)*ZP(gt2,2) + ZP(gt1,3)*ZP(gt2,3)));

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::VP, MRSSMEFTHiggs_cxx_diagrams::fields::Sd, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Sd>::type>::evaluate(
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

   const std::complex<double> result = -0.03333333333333333*(3.872983346207417*g1*Cos(ThetaW) - 15*g2*Sin(ThetaW))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1)) + 0.2581988897471611*g1*Cos(ThetaW)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1));

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::VP, MRSSMEFTHiggs_cxx_diagrams::fields::Se, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Se>::type>::evaluate(
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

   const std::complex<double> result = 0.5*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) + 0.7745966692414834*g1*Cos(ThetaW)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1));

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::VP, MRSSMEFTHiggs_cxx_diagrams::fields::SRdp, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type>::evaluate(
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

cxx_diagrams::MomentumDifferenceVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::VP, MRSSMEFTHiggs_cxx_diagrams::fields::SRum, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRum>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::VP, MRSSMEFTHiggs_cxx_diagrams::fields::Su, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Su>::type>::evaluate(
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

   const std::complex<double> result = -0.16666666666666666*(0.7745966692414834*g1*Cos(ThetaW) + 3*g2*Sin(ThetaW))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1)) - 0.5163977794943222*g1*Cos(ThetaW)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1));

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::TripleVectorVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::VP, MRSSMEFTHiggs_cxx_diagrams::fields::VWm, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::VWm>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Sin(ThetaW);

   return {result, cxx_diagrams::TripleVectorVertex::odd_permutation{}};
}

cxx_diagrams::ChiralVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::VP, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::VZ, MRSSMEFTHiggs_cxx_diagrams::fields::Cha1, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha1>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = g2*Conj(UM1(gt1,0))*Cos(ThetaW)*UM1(gt2,0) + 0.5*Conj(UM1(gt1,1))*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*UM1(gt2,1);

   const std::complex<double> right = g2*Conj(UP1(gt2,0))*Cos(ThetaW)*UP1(gt1,0) + 0.5*Conj(UP1(gt2,1))*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW))*UP1(gt1,1);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::VZ, MRSSMEFTHiggs_cxx_diagrams::fields::Cha2, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha2>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UP2 = MODELPARAMETER(UP2);
   const auto UM2 = MODELPARAMETER(UM2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -(g2*Conj(UP2(gt1,0))*Cos(ThetaW)*UP2(gt2,0)) + 0.1*Conj(UP2(gt1,1))*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*UP2(gt2,1);

   const std::complex<double> right = -(g2*Conj(UM2(gt2,0))*Cos(ThetaW)*UM2(gt1,0)) + 0.1*Conj(UM2(gt2,1))*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*UM2(gt1,1);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::VZ, MRSSMEFTHiggs_cxx_diagrams::fields::Fd, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fd>::type>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::VZ, MRSSMEFTHiggs_cxx_diagrams::fields::Fe, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::type>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::VZ, MRSSMEFTHiggs_cxx_diagrams::fields::Fu, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fu>::type>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<MRSSMEFTHiggs_cxx_diagrams::fields::VZ, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Chi>::type, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Se>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZN2 = MODELPARAMETER(ZN2);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZN1 = MODELPARAMETER(ZN1);

   const std::complex<double> left = -(Conj(ZN2(gt1,2))*SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1))));

   const std::complex<double> right = -1.0954451150103321*g1*SUM(j1,0,2,ZE(gt3,3 + j1)*ZER(gt2,j1))*ZN1(gt1,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Fe, MRSSMEFTHiggs_cxx_diagrams::fields::Ah>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Fe, MRSSMEFTHiggs_cxx_diagrams::fields::hh>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Fe, MRSSMEFTHiggs_cxx_diagrams::fields::VP>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Fe, MRSSMEFTHiggs_cxx_diagrams::fields::VZ>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::Fv>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> left = SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,gt2))*ZP(gt3,0);

   const std::complex<double> right = 0;

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Se, MRSSMEFTHiggs_cxx_diagrams::fields::Chi>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZN1 = MODELPARAMETER(ZN1);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZN2 = MODELPARAMETER(ZN2);

   const std::complex<double> left = -1.0954451150103321*g1*Conj(ZN1(gt2,0))*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*Conj(ZER(gt1,j1)));

   const std::complex<double> right = -(SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt3,3 + j1)))*ZEL(gt1,j2))*ZN2(gt2,2));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Se, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Chi>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt3 = indices[1];
   const int gt1 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZN2 = MODELPARAMETER(ZN2);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZN1 = MODELPARAMETER(ZN1);

   const std::complex<double> left = -(Conj(ZN2(gt1,2))*SUM(j2,0,2,Conj(ZE(gt3,j2))*SUM(j1,0,2,Conj(ZER(gt2,j1))*Ye(j1,j2))));

   const std::complex<double> right = 0.7071067811865475*SUM(j1,0,2,Conj(ZE(gt3,j1))*ZEL(gt2,j1))*(0.7745966692414834*g1*ZN1(gt1,0) + g2*ZN1(gt1,1));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::type, MRSSMEFTHiggs_cxx_diagrams::fields::VWm, MRSSMEFTHiggs_cxx_diagrams::fields::Fv>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::type, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha1>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Sv>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto UM1 = MODELPARAMETER(UM1);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto UP1 = MODELPARAMETER(UP1);

   const std::complex<double> left = Conj(UM1(gt1,1))*SUM(j2,0,2,Conj(ZV(gt3,j2))*SUM(j1,0,2,Conj(ZER(gt2,j1))*Ye(j1,j2)));

   const std::complex<double> right = -(g2*SUM(j1,0,2,Conj(ZV(gt3,j1))*ZEL(gt2,j1))*UP1(gt1,0));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fv>::type, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> left = 0;

   const std::complex<double> right = SUM(j1,0,2,Conj(Ye(j1,gt1))*ZER(gt2,j1))*ZP(gt3,0);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fv>::type, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::VWm>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::evaluate(
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

cxx_diagrams::MomentumDifferenceVertex VertexImpl<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::VP>::evaluate(
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

   const std::complex<double> result = 0.5*((0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*ZP(gt1,0)*ZP(gt2,0) + (0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*ZP(gt1,1)*ZP(gt2,1) + 2*g2*Sin(ThetaW)*(ZP(gt1,2)*ZP(gt2,2) + ZP(gt1,3)*ZP(gt2,3)));

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::InverseMetricVertex VertexImpl<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type, MRSSMEFTHiggs_cxx_diagrams::fields::VWm, MRSSMEFTHiggs_cxx_diagrams::fields::VP>::evaluate(
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

cxx_diagrams::MomentumDifferenceVertex VertexImpl<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Se>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Se, MRSSMEFTHiggs_cxx_diagrams::fields::VP>::evaluate(
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

cxx_diagrams::ChiralVertex VertexImpl<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Sv>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Cha1, MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto UP1 = MODELPARAMETER(UP1);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZER = MODELPARAMETER(ZER);
   const auto UM1 = MODELPARAMETER(UM1);

   const std::complex<double> left = -(g2*Conj(UP1(gt1,0))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZV(gt3,j1)));

   const std::complex<double> right = SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZV(gt3,j2))*UM1(gt1,1);

   return {left, right};
}

cxx_diagrams::InverseMetricVertex VertexImpl<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::VWm>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::VP>::evaluate(
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

cxx_diagrams::TripleVectorVertex VertexImpl<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::VWm>::type, MRSSMEFTHiggs_cxx_diagrams::fields::VWm, MRSSMEFTHiggs_cxx_diagrams::fields::VP>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Sin(ThetaW));

   return {result, cxx_diagrams::TripleVectorVertex::odd_permutation{}};
}

} // namespace flexiblesusy::MRSSMEFTHiggs_cxx_diagrams::detail
