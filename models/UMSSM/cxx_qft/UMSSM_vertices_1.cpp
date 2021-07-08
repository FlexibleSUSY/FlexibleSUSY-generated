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
 * @file cxx_qft/UMSSM_vertices.cpp
 *
 * This file was generated with FlexibleSUSY 2.6.1 and SARAH 4.14.5 .
 */

#include "UMSSM_context_base.hpp"
#include "UMSSM_input_parameters.hpp"
#include "UMSSM_vertices.hpp"

#include "concatenate.hpp"
#include "wrappers.hpp"

#define INPUTPARAMETER(p) context.model.get_input_parameters().p
#define MODELPARAMETER(p) context.model.get_##p()
#define DERIVEDPARAMETER(p) context.model.p()
#define PHASE(p) context.model.get_##p()

namespace flexiblesusy {
namespace UMSSM_cxx_diagrams {
namespace detail {

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Ah, fields::Ah>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = 0.25*(-(Conj(ZA(gt1,0))*(-(Conj(ZA(gt2,1))*(Conj(ZA(gt3,1))*Conj(ZA(gt4,0)) + Conj(ZA(gt3,0))*Conj(ZA(gt4,1)))*(-4*AbsSqr(Lambdax) + 0.6*Sqr(g1) + Sqr(g2) - 4*QHd*QHu*Sqr(gp))) + 4*Conj(ZA(gt2,2))*(Conj(ZA(gt3,2))*Conj(ZA(gt4,0)) + Conj(ZA(gt3,0))*Conj(ZA(gt4,2)))*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + Conj(ZA(gt2,0))*(-(Conj(ZA(gt3,1))*Conj(ZA(gt4,1))*(-4*AbsSqr(Lambdax) + 0.6*Sqr(g1) + Sqr(g2) - 4*QHd*QHu*Sqr(gp))) + 4*Conj(ZA(gt3,2))*Conj(ZA(gt4,2))*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + 3*Conj(ZA(gt3,0))*Conj(ZA(gt4,0))*(0.6*Sqr(g1) + Sqr(g2) + 4*Sqr(gp)*Sqr(QHd))))) - Conj(ZA(gt1,1))*(-(Conj(ZA(gt2,0))*(Conj(ZA(gt3,1))*Conj(ZA(gt4,0)) + Conj(ZA(gt3,0))*Conj(ZA(gt4,1)))*(-4*AbsSqr(Lambdax) + 0.6*Sqr(g1) + Sqr(g2) - 4*QHd*QHu*Sqr(gp))) + 4*Conj(ZA(gt2,2))*(Conj(ZA(gt3,2))*Conj(ZA(gt4,1)) + Conj(ZA(gt3,1))*Conj(ZA(gt4,2)))*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) + Conj(ZA(gt2,1))*(-(Conj(ZA(gt3,0))*Conj(ZA(gt4,0))*(-4*AbsSqr(Lambdax) + 0.6*Sqr(g1) + Sqr(g2) - 4*QHd*QHu*Sqr(gp))) + 4*Conj(ZA(gt3,2))*Conj(ZA(gt4,2))*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) + 3*Conj(ZA(gt3,1))*Conj(ZA(gt4,1))*(0.6*Sqr(g1) + Sqr(g2) + 4*Sqr(gp)*Sqr(QHu)))) - 4*Conj(ZA(gt1,2))*(Conj(ZA(gt2,0))*(Conj(ZA(gt3,2))*Conj(ZA(gt4,0)) + Conj(ZA(gt3,0))*Conj(ZA(gt4,2)))*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + Conj(ZA(gt2,1))*(Conj(ZA(gt3,2))*Conj(ZA(gt4,1)) + Conj(ZA(gt3,1))*Conj(ZA(gt4,2)))*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) + Conj(ZA(gt2,2))*(Conj(ZA(gt3,0))*Conj(ZA(gt4,0))*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + Conj(ZA(gt3,1))*Conj(ZA(gt4,1))*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) + 3*Conj(ZA(gt3,2))*Conj(ZA(gt4,2))*Sqr(gp)*Sqr(Qs))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Ah>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto TLambdax = MODELPARAMETER(TLambdax);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = std::complex<double>(0.,0.35355339059327373)*(Conj(ZA(gt1,2))*(Conj(ZA(gt2,1))*Conj(ZA(gt3,0)) + Conj(ZA(gt2,0))*Conj(ZA(gt3,1))) + Conj(ZA(gt1,1))*(Conj(ZA(gt2,2))*Conj(ZA(gt3,0)) + Conj(ZA(gt2,0))*Conj(ZA(gt3,2))) + Conj(ZA(gt1,0))*(Conj(ZA(gt2,2))*Conj(ZA(gt3,1)) + Conj(ZA(gt2,1))*Conj(ZA(gt3,2))))*(Conj(TLambdax) - TLambdax);

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::hh, fields::hh>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.25*(-(Conj(ZA(gt1,0))*Conj(ZA(gt2,0))*(-(Conj(ZH(gt3,1))*Conj(ZH(gt4,1))*(-4*AbsSqr(Lambdax) + 0.6*Sqr(g1) + Sqr(g2) - 4*QHd*QHu*Sqr(gp))) + 4*Conj(ZH(gt3,2))*Conj(ZH(gt4,2))*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + Conj(ZH(gt3,0))*Conj(ZH(gt4,0))*(0.6*Sqr(g1) + Sqr(g2) + 4*Sqr(gp)*Sqr(QHd)))) - Conj(ZA(gt1,1))*Conj(ZA(gt2,1))*(-(Conj(ZH(gt3,0))*Conj(ZH(gt4,0))*(-4*AbsSqr(Lambdax) + 0.6*Sqr(g1) + Sqr(g2) - 4*QHd*QHu*Sqr(gp))) + 4*Conj(ZH(gt3,2))*Conj(ZH(gt4,2))*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) + Conj(ZH(gt3,1))*Conj(ZH(gt4,1))*(0.6*Sqr(g1) + Sqr(g2) + 4*Sqr(gp)*Sqr(QHu))) - 4*Conj(ZA(gt1,2))*Conj(ZA(gt2,2))*(Conj(ZH(gt3,0))*Conj(ZH(gt4,0))*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + Conj(ZH(gt3,1))*Conj(ZH(gt4,1))*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) + Conj(ZH(gt3,2))*Conj(ZH(gt4,2))*Sqr(gp)*Sqr(Qs)));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto TLambdax = MODELPARAMETER(TLambdax);
   const auto vd = MODELPARAMETER(vd);
   const auto gp = MODELPARAMETER(gp);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto vu = MODELPARAMETER(vu);
   const auto vS = MODELPARAMETER(vS);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.25*(-(Conj(ZA(gt1,2))*(1.4142135623730951*Conj(TLambdax)*(Conj(ZA(gt2,1))*Conj(ZH(gt3,0)) + Conj(ZA(gt2,0))*Conj(ZH(gt3,1))) + 4*Conj(ZA(gt2,2))*(vd*Conj(ZH(gt3,0))*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + vu*Conj(ZH(gt3,1))*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) + vS*Conj(ZH(gt3,2))*Sqr(gp)*Sqr(Qs)) + 1.4142135623730951*(Conj(ZA(gt2,1))*Conj(ZH(gt3,0)) + Conj(ZA(gt2,0))*Conj(ZH(gt3,1)))*TLambdax)) - Conj(ZA(gt1,1))*(Conj(ZA(gt2,1))*(-(vd*Conj(ZH(gt3,0))*(-4*AbsSqr(Lambdax) + 0.6*Sqr(g1) + Sqr(g2) - 4*QHd*QHu*Sqr(gp))) + 4*vS*Conj(ZH(gt3,2))*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) + vu*Conj(ZH(gt3,1))*(0.6*Sqr(g1) + Sqr(g2) + 4*Sqr(gp)*Sqr(QHu))) + 1.4142135623730951*(Conj(ZA(gt2,2))*Conj(ZH(gt3,0)) + Conj(ZA(gt2,0))*Conj(ZH(gt3,2)))*(Conj(TLambdax) + TLambdax)) - Conj(ZA(gt1,0))*(Conj(ZA(gt2,0))*(-(vu*Conj(ZH(gt3,1))*(-4*AbsSqr(Lambdax) + 0.6*Sqr(g1) + Sqr(g2) - 4*QHd*QHu*Sqr(gp))) + 4*vS*Conj(ZH(gt3,2))*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + vd*Conj(ZH(gt3,0))*(0.6*Sqr(g1) + Sqr(g2) + 4*Sqr(gp)*Sqr(QHd))) + 1.4142135623730951*(Conj(ZA(gt2,2))*Conj(ZH(gt3,1)) + Conj(ZA(gt2,1))*Conj(ZH(gt3,2)))*(Conj(TLambdax) + TLambdax)));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Hpm, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-4*Conj(ZA(gt1,2))*Conj(ZA(gt2,2))*((AbsSqr(Lambdax) + QHd*Qs*Sqr(gp))*ZP(gt3,0)*ZP(gt4,0) + (AbsSqr(Lambdax) + QHu*Qs*Sqr(gp))*ZP(gt3,1)*ZP(gt4,1)) - Conj(ZA(gt1,0))*(-(Conj(ZA(gt2,1))*(-2*AbsSqr(Lambdax) + Sqr(g2))*(ZP(gt3,1)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,1))) + Conj(ZA(gt2,0))*((0.6*Sqr(g1) + Sqr(g2) + 4*Sqr(gp)*Sqr(QHd))*ZP(gt3,0)*ZP(gt4,0) + (-0.6*Sqr(g1) + Sqr(g2) + 4*QHd*QHu*Sqr(gp))*ZP(gt3,1)*ZP(gt4,1))) - Conj(ZA(gt1,1))*(-(Conj(ZA(gt2,0))*(-2*AbsSqr(Lambdax) + Sqr(g2))*(ZP(gt3,1)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,1))) + Conj(ZA(gt2,1))*((-0.6*Sqr(g1) + Sqr(g2) + 4*QHd*QHu*Sqr(gp))*ZP(gt3,0)*ZP(gt4,0) + (0.6*Sqr(g1) + Sqr(g2) + 4*Sqr(gp)*Sqr(QHu))*ZP(gt3,1)*ZP(gt4,1))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Sd, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qd = INPUTPARAMETER(Qd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto gp = MODELPARAMETER(gp);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = 0.05*(-10*Conj(ZA(gt1,2))*(2*Qs*Conj(ZA(gt2,2))*Sqr(gp)*(Qq*SUM(j1,0,2,Conj(ZD(gt3,j1))*ZD(gt4,j1)) + Qd*SUM(j1,0,2,Conj(ZD(gt3,3 + j1))*ZD(gt4,3 + j1))) + Conj(ZA(gt2,1))*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZD(gt3,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt4,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt3,3 + j1)))*ZD(gt4,j2)))) - Conj(ZA(gt1,1))*(Conj(ZA(gt2,1))*((Sqr(g1) + 5*Sqr(g2) + 20*QHu*Qq*Sqr(gp))*SUM(j1,0,2,Conj(ZD(gt3,j1))*ZD(gt4,j1)) + 2*(Sqr(g1) + 10*Qd*QHu*Sqr(gp))*SUM(j1,0,2,Conj(ZD(gt3,3 + j1))*ZD(gt4,3 + j1))) + 10*Conj(ZA(gt2,2))*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZD(gt3,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt4,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt3,3 + j1)))*ZD(gt4,j2)))) + Conj(ZA(gt1,0))*Conj(ZA(gt2,0))*((Sqr(g1) + 5*(Sqr(g2) - 4*QHd*Qq*Sqr(gp)))*SUM(j1,0,2,Conj(ZD(gt3,j1))*ZD(gt4,j1)) + 2*(Sqr(g1) - 10*Qd*QHd*Sqr(gp))*SUM(j1,0,2,Conj(ZD(gt3,3 + j1))*ZD(gt4,3 + j1)) - 20*(SUM(j3,0,2,Conj(ZD(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,j1))*ZD(gt4,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt3,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gt4,j3)))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Se, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qe = INPUTPARAMETER(Qe);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto gp = MODELPARAMETER(gp);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZE = MODELPARAMETER(ZE);

   const std::complex<double> result = 0.25*(-2*Conj(ZA(gt1,2))*(2*Qs*Conj(ZA(gt2,2))*Sqr(gp)*(Ql*SUM(j1,0,2,Conj(ZE(gt3,j1))*ZE(gt4,j1)) + Qe*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*ZE(gt4,3 + j1))) + Conj(ZA(gt2,1))*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZE(gt3,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt4,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt3,3 + j1)))*ZE(gt4,j2)))) - Conj(ZA(gt1,1))*(Conj(ZA(gt2,1))*((-0.6*Sqr(g1) + Sqr(g2) + 4*QHu*Ql*Sqr(gp))*SUM(j1,0,2,Conj(ZE(gt3,j1))*ZE(gt4,j1)) + 2*(0.6*Sqr(g1) + 2*Qe*QHu*Sqr(gp))*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*ZE(gt4,3 + j1))) + 2*Conj(ZA(gt2,2))*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZE(gt3,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt4,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt3,3 + j1)))*ZE(gt4,j2)))) - 0.2*Conj(ZA(gt1,0))*Conj(ZA(gt2,0))*((3*Sqr(g1) - 5*Sqr(g2) + 20*QHd*Ql*Sqr(gp))*SUM(j1,0,2,Conj(ZE(gt3,j1))*ZE(gt4,j1)) + (-6*Sqr(g1) + 20*Qe*QHd*Sqr(gp))*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*ZE(gt4,3 + j1)) + 20*(SUM(j3,0,2,Conj(ZE(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Ye(j2,j1))*ZE(gt4,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gt3,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZE(gt4,j3)))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Su, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qu = INPUTPARAMETER(Qu);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZU = MODELPARAMETER(ZU);

   const std::complex<double> result = 0.05*(-10*Conj(ZA(gt1,2))*(2*Qs*Conj(ZA(gt2,2))*Sqr(gp)*(Qq*SUM(j1,0,2,Conj(ZU(gt3,j1))*ZU(gt4,j1)) + Qu*SUM(j1,0,2,Conj(ZU(gt3,3 + j1))*ZU(gt4,3 + j1))) + Conj(ZA(gt2,0))*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt4,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt3,3 + j1)))*ZU(gt4,j2)))) + Conj(ZA(gt1,0))*(Conj(ZA(gt2,0))*((Sqr(g1) - 5*(Sqr(g2) + 4*QHd*Qq*Sqr(gp)))*SUM(j1,0,2,Conj(ZU(gt3,j1))*ZU(gt4,j1)) - 4*(Sqr(g1) + 5*QHd*Qu*Sqr(gp))*SUM(j1,0,2,Conj(ZU(gt3,3 + j1))*ZU(gt4,3 + j1))) - 10*Conj(ZA(gt2,2))*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt4,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt3,3 + j1)))*ZU(gt4,j2)))) - Conj(ZA(gt1,1))*Conj(ZA(gt2,1))*((Sqr(g1) - 5*Sqr(g2) + 20*QHu*Qq*Sqr(gp))*SUM(j1,0,2,Conj(ZU(gt3,j1))*ZU(gt4,j1)) - 4*(Sqr(g1) - 5*QHu*Qu*Sqr(gp))*SUM(j1,0,2,Conj(ZU(gt3,3 + j1))*ZU(gt4,3 + j1)) + 20*(SUM(j3,0,2,Conj(ZU(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))*ZU(gt4,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gt4,j3)))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::Sv, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qv = INPUTPARAMETER(Qv);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZV = MODELPARAMETER(ZV);

   const std::complex<double> result = 0.25*(-2*Conj(ZA(gt1,2))*(2*Qs*Conj(ZA(gt2,2))*Sqr(gp)*(Ql*SUM(j1,0,2,Conj(ZV(gt3,j1))*ZV(gt4,j1)) + Qv*SUM(j1,0,2,Conj(ZV(gt3,3 + j1))*ZV(gt4,3 + j1))) + Conj(ZA(gt2,0))*(Lambdax*SUM(j2,0,2,Conj(ZV(gt3,3 + j2))*SUM(j1,0,2,Conj(Yv(j1,j2))*ZV(gt4,j1))) + Conj(Lambdax)*SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(gt3,j1))*Yv(j1,j2))*ZV(gt4,3 + j2)))) - Conj(ZA(gt1,0))*(Conj(ZA(gt2,0))*((0.6*Sqr(g1) + Sqr(g2) + 4*QHd*Ql*Sqr(gp))*SUM(j1,0,2,Conj(ZV(gt3,j1))*ZV(gt4,j1)) + 4*QHd*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gt3,3 + j1))*ZV(gt4,3 + j1))) + 2*Conj(ZA(gt2,2))*(Lambdax*SUM(j2,0,2,Conj(ZV(gt3,3 + j2))*SUM(j1,0,2,Conj(Yv(j1,j2))*ZV(gt4,j1))) + Conj(Lambdax)*SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(gt3,j1))*Yv(j1,j2))*ZV(gt4,3 + j2)))) - Conj(ZA(gt1,1))*Conj(ZA(gt2,1))*(-((0.6*Sqr(g1) + Sqr(g2) - 4*QHu*Ql*Sqr(gp))*SUM(j1,0,2,Conj(ZV(gt3,j1))*ZV(gt4,j1))) + 4*(QHu*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gt3,3 + j1))*ZV(gt4,3 + j1)) + SUM(j3,0,2,Conj(ZV(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j3))*Yv(j1,j2))*ZV(gt4,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gt3,j2))*SUM(j1,0,2,Conj(Yv(j3,j1))*Yv(j2,j1)))*ZV(gt4,j3)))));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Ah, fields::Ah, fields::VZ, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.1*(20*Conj(ZA(gt1,2))*Conj(ZA(gt2,2))*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gp)*Sqr(Qs) - Conj(ZA(gt1,0))*Conj(ZA(gt2,0))*(5*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) - 7.745966692414834*g1*gp*QHd*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + Cos(ThetaWp)*Sin(ThetaWp)*(-20*Sqr(gp)*Sqr(QHd) + 3*Sqr(g1)*Sqr(Sin(ThetaW))) + 7.745966692414834*g1*gp*QHd*Sin(ThetaW)*Sqr(Sin(ThetaWp)) + 2*g2*Cos(ThetaW)*(3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) - 5*gp*QHd*Sqr(Cos(ThetaWp)) + 5*gp*QHd*Sqr(Sin(ThetaWp)))) - Conj(ZA(gt1,1))*Conj(ZA(gt2,1))*(5*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 7.745966692414834*g1*gp*QHu*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + Cos(ThetaWp)*Sin(ThetaWp)*(-20*Sqr(gp)*Sqr(QHu) + 3*Sqr(g1)*Sqr(Sin(ThetaW))) - 7.745966692414834*g1*gp*QHu*Sin(ThetaW)*Sqr(Sin(ThetaWp)) + 2*g2*Cos(ThetaW)*(3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) + 5*gp*QHu*Sqr(Cos(ThetaWp)) - 5*gp*QHu*Sqr(Sin(ThetaWp)))));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Ah, fields::Ah, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(4*Conj(ZA(gt1,2))*Conj(ZA(gt2,2))*Sqr(gp)*Sqr(Qs)*Sqr(Sin(ThetaWp)) + Conj(ZA(gt1,0))*Conj(ZA(gt2,0))*Sqr(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHd*Sin(ThetaWp)) + Conj(ZA(gt1,1))*Conj(ZA(gt2,1))*Sqr(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*QHu*Sin(ThetaWp)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Ah, fields::Ah, fields::VZp, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(4*Conj(ZA(gt1,2))*Conj(ZA(gt2,2))*Sqr(gp)*Sqr(Qs)*Sqr(Cos(ThetaWp)) + Conj(ZA(gt1,0))*Conj(ZA(gt2,0))*Sqr(-2*gp*QHd*Cos(ThetaWp) + g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp)) + Conj(ZA(gt1,1))*Conj(ZA(gt2,1))*Sqr(2*gp*QHu*Cos(ThetaWp) + g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Ah, fields::Ah, typename fields::conj<fields::VWm>::type, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> result = 0.5*(Conj(ZA(gt1,0))*Conj(ZA(gt2,0)) + Conj(ZA(gt1,1))*Conj(ZA(gt2,1)))*Sqr(g2);

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::hh, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto TLambdax = MODELPARAMETER(TLambdax);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = std::complex<double>(0.,-0.35355339059327373)*(Conj(ZA(gt1,2))*(Conj(ZH(gt2,1))*Conj(ZH(gt3,0)) + Conj(ZH(gt2,0))*Conj(ZH(gt3,1))) + Conj(ZA(gt1,1))*(Conj(ZH(gt2,2))*Conj(ZH(gt3,0)) + Conj(ZH(gt2,0))*Conj(ZH(gt3,2))) + Conj(ZA(gt1,0))*(Conj(ZH(gt2,2))*Conj(ZH(gt3,1)) + Conj(ZH(gt2,1))*Conj(ZH(gt3,2))))*(Conj(TLambdax) - TLambdax);

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
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,-0.25)*(Conj(ZA(gt1,1))*Conj(ZH(gt2,0)) + Conj(ZA(gt1,0))*Conj(ZH(gt2,1)))*(-2*AbsSqr(Lambdax) + Sqr(g2))*(ZP(gt3,1)*ZP(gt4,0) - ZP(gt3,0)*ZP(gt4,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::hh, fields::Sd, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = std::complex<double>(0,-0.5)*(Conj(ZA(gt1,2))*Conj(ZH(gt2,1)) + Conj(ZA(gt1,1))*Conj(ZH(gt2,2)))*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZD(gt3,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt4,3 + j1))) - Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt3,3 + j1)))*ZD(gt4,j2)));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::hh, fields::Se, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZE = MODELPARAMETER(ZE);

   const std::complex<double> result = std::complex<double>(0,-0.5)*(Conj(ZA(gt1,2))*Conj(ZH(gt2,1)) + Conj(ZA(gt1,1))*Conj(ZH(gt2,2)))*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZE(gt3,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt4,3 + j1))) - Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt3,3 + j1)))*ZE(gt4,j2)));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::hh, fields::Su, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZU = MODELPARAMETER(ZU);

   const std::complex<double> result = std::complex<double>(0,-0.5)*(Conj(ZA(gt1,2))*Conj(ZH(gt2,0)) + Conj(ZA(gt1,0))*Conj(ZH(gt2,2)))*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt4,3 + j1))) - Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt3,3 + j1)))*ZU(gt4,j2)));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::hh, fields::Sv, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZV = MODELPARAMETER(ZV);

   const std::complex<double> result = std::complex<double>(0,0.5)*(Conj(ZA(gt1,2))*Conj(ZH(gt2,0)) + Conj(ZA(gt1,0))*Conj(ZH(gt2,2)))*(Lambdax*SUM(j2,0,2,Conj(ZV(gt3,3 + j2))*SUM(j1,0,2,Conj(Yv(j1,j2))*ZV(gt4,j1))) - Conj(Lambdax)*SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(gt3,j1))*Yv(j1,j2))*ZV(gt4,3 + j2)));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::Ah, fields::hh, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = std::complex<double>(0,-0.5)*(2*gp*Qs*Conj(ZA(gt1,2))*Conj(ZH(gt2,2))*Cos(ThetaWp) + Conj(ZA(gt1,1))*Conj(ZH(gt2,1))*(2*gp*QHu*Cos(ThetaWp) + g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp)) + Conj(ZA(gt1,0))*Conj(ZH(gt2,0))*(2*gp*QHd*Cos(ThetaWp) - (g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*Sin(ThetaWp)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::Ah, fields::hh, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = std::complex<double>(0,-0.5)*(2*gp*Qs*Conj(ZA(gt1,2))*Conj(ZH(gt2,2))*Sin(ThetaWp) + Conj(ZA(gt1,0))*Conj(ZH(gt2,0))*(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHd*Sin(ThetaWp)) - Conj(ZA(gt1,1))*Conj(ZH(gt2,1))*(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*QHu*Sin(ThetaWp)));

   return {result, minuend_index, subtrahend_index};
}

ScalarVertex VertexImpl<fields::Ah, fields::Hpm, fields::Su, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0.,0.35355339059327373)*(2*Conj(ZA(gt1,2))*(Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt3,3 + j1)))*ZD(gt4,j2))*ZP(gt2,0) - Conj(Lambdax)*SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt4,3 + j1)))*ZP(gt2,1)) + Conj(ZA(gt1,0))*(Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt3,j1))*ZD(gt4,j1))*ZP(gt2,0) - 2*SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gt4,j3))*ZP(gt2,0) + 2*SUM(j3,0,2,Conj(ZU(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yd(j2,j1))*ZD(gt4,3 + j2)))*ZP(gt2,1)) - Conj(ZA(gt1,1))*(2*SUM(j3,0,2,Conj(ZU(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yd(j2,j1))*ZD(gt4,3 + j2)))*ZP(gt2,0) + (Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt3,j1))*ZD(gt4,j1)) - 2*SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZD(gt4,j3)))*ZP(gt2,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Hpm, fields::Sv, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yv = MODELPARAMETER(Yv);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0.,0.35355339059327373)*(2*Conj(ZA(gt1,2))*(Lambdax*SUM(j2,0,2,Conj(ZV(gt3,3 + j2))*SUM(j1,0,2,Conj(Yv(j1,j2))*ZE(gt4,j1)))*ZP(gt2,0) - Conj(Lambdax)*SUM(j2,0,2,Conj(ZV(gt3,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt4,3 + j1)))*ZP(gt2,1)) + Conj(ZA(gt1,0))*(Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt3,j1))*ZE(gt4,j1))*ZP(gt2,0) - 2*SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gt3,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZE(gt4,j3))*ZP(gt2,0) + 2*SUM(j3,0,2,Conj(ZV(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j3))*Ye(j2,j1))*ZE(gt4,3 + j2)))*ZP(gt2,1)) - Conj(ZA(gt1,1))*(2*SUM(j3,0,2,Conj(ZV(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j3))*Ye(j2,j1))*ZE(gt4,3 + j2)))*ZP(gt2,0) + (Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt3,j1))*ZE(gt4,j1)) - 2*SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gt3,j2))*SUM(j1,0,2,Conj(Yv(j3,j1))*Yv(j2,j1)))*ZE(gt4,j3)))*ZP(gt2,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Hpm, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto TLambdax = MODELPARAMETER(TLambdax);
   const auto vu = MODELPARAMETER(vu);
   const auto g2 = MODELPARAMETER(g2);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto vd = MODELPARAMETER(vd);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0,-0.25)*(vu*Conj(ZA(gt1,0))*(-2*AbsSqr(Lambdax) + Sqr(g2))*(ZP(gt2,1)*ZP(gt3,0) - ZP(gt2,0)*ZP(gt3,1)) + vd*Conj(ZA(gt1,1))*(-2*AbsSqr(Lambdax) + Sqr(g2))*(ZP(gt2,1)*ZP(gt3,0) - ZP(gt2,0)*ZP(gt3,1)) + 2.8284271247461903*Conj(ZA(gt1,2))*(-(Conj(TLambdax)*ZP(gt2,1)*ZP(gt3,0)) + TLambdax*ZP(gt2,0)*ZP(gt3,1)));

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

   const std::complex<double> result = std::complex<double>(0.,0.3872983346207417)*g1*g2*Cos(ThetaW)*(Conj(ZA(gt1,0))*ZP(gt2,0) + Conj(ZA(gt1,1))*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Ah, fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = std::complex<double>(0,-0.5)*g2*(Conj(ZA(gt1,0))*(2*gp*QHd*Cos(ThetaWp) - 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt2,0) - Conj(ZA(gt1,1))*(2*gp*QHu*Cos(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Ah, fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = std::complex<double>(0,-0.5)*g2*(Conj(ZA(gt1,0))*(0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHd*Sin(ThetaWp))*ZP(gt2,0) + Conj(ZA(gt1,1))*(0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*QHu*Sin(ThetaWp))*ZP(gt2,1));

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

   const std::complex<double> result = std::complex<double>(0,-0.5)*g2*(Conj(ZA(gt1,0))*ZP(gt2,0) + Conj(ZA(gt1,1))*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

ScalarVertex VertexImpl<fields::Ah, fields::Sd, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0.,-0.35355339059327373)*(2*Conj(Lambdax)*Conj(ZA(gt1,2))*SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt4,3 + j1)))*ZP(gt3,0) - 2*Conj(ZA(gt1,1))*SUM(j3,0,2,Conj(ZD(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yu(j2,j1))*ZU(gt4,3 + j2)))*ZP(gt3,0) - Conj(ZA(gt1,1))*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZU(gt4,j1))*ZP(gt3,1) - 2*Conj(ZA(gt1,2))*Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1)))*ZU(gt4,j2))*ZP(gt3,1) + 2*Conj(ZA(gt1,1))*SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gt4,j3))*ZP(gt3,1) + Conj(ZA(gt1,0))*(Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZU(gt4,j1))*ZP(gt3,0) - 2*SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZU(gt4,j3))*ZP(gt3,0) + 2*SUM(j3,0,2,Conj(ZD(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yu(j2,j1))*ZU(gt4,3 + j2)))*ZP(gt3,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Sd, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto TYd = MODELPARAMETER(TYd);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto vS = MODELPARAMETER(vS);
   const auto vu = MODELPARAMETER(vu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = std::complex<double>(0,-0.5)*(Conj(Lambdax)*(vS*Conj(ZA(gt1,1)) + vu*Conj(ZA(gt1,2)))*SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1))) - (vS*Conj(ZA(gt1,1)) + vu*Conj(ZA(gt1,2)))*Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1)))*ZD(gt3,j2)) + 1.4142135623730951*Conj(ZA(gt1,0))*(SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,ZD(gt3,3 + j1)*TYd(j1,j2))) - SUM(j2,0,2,SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j2)))*ZD(gt3,j2))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Se, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yv = MODELPARAMETER(Yv);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = std::complex<double>(0.,-0.35355339059327373)*(2*Conj(Lambdax)*Conj(ZA(gt1,2))*SUM(j2,0,2,SUM(j1,0,2,Conj(ZE(gt2,j1))*Yv(j1,j2))*ZV(gt4,3 + j2))*ZP(gt3,0) - 2*Conj(ZA(gt1,1))*SUM(j3,0,2,Conj(ZE(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Yv(j1,j2))*ZV(gt4,3 + j2)))*ZP(gt3,0) - Conj(ZA(gt1,1))*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZV(gt4,j1))*ZP(gt3,1) - 2*Conj(ZA(gt1,2))*Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1)))*ZV(gt4,j2))*ZP(gt3,1) + 2*Conj(ZA(gt1,1))*SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Conj(Yv(j3,j1))*Yv(j2,j1)))*ZV(gt4,j3))*ZP(gt3,1) + Conj(ZA(gt1,0))*(Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZV(gt4,j1))*ZP(gt3,0) - 2*SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZV(gt4,j3))*ZP(gt3,0) + 2*SUM(j3,0,2,Conj(ZE(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Yv(j1,j2))*ZV(gt4,3 + j2)))*ZP(gt3,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Se, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto TYe = MODELPARAMETER(TYe);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto vS = MODELPARAMETER(vS);
   const auto vu = MODELPARAMETER(vu);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZE = MODELPARAMETER(ZE);

   const std::complex<double> result = std::complex<double>(0,-0.5)*(Conj(Lambdax)*(vS*Conj(ZA(gt1,1)) + vu*Conj(ZA(gt1,2)))*SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1))) - (vS*Conj(ZA(gt1,1)) + vu*Conj(ZA(gt1,2)))*Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1)))*ZE(gt3,j2)) + 1.4142135623730951*Conj(ZA(gt1,0))*(SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,ZE(gt3,3 + j1)*TYe(j1,j2))) - SUM(j2,0,2,SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j2)))*ZE(gt3,j2))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Su, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto TYu = MODELPARAMETER(TYu);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto vS = MODELPARAMETER(vS);
   const auto vd = MODELPARAMETER(vd);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZU = MODELPARAMETER(ZU);

   const std::complex<double> result = std::complex<double>(0,-0.5)*(Conj(Lambdax)*(vS*Conj(ZA(gt1,0)) + vd*Conj(ZA(gt1,2)))*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt3,3 + j1))) - (vS*Conj(ZA(gt1,0)) + vd*Conj(ZA(gt1,2)))*Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1)))*ZU(gt3,j2)) + 1.4142135623730951*Conj(ZA(gt1,1))*(SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,ZU(gt3,3 + j1)*TYu(j1,j2))) - SUM(j2,0,2,SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j2)))*ZU(gt3,j2))));

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Sv, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto TYv = MODELPARAMETER(TYv);
   const auto vd = MODELPARAMETER(vd);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto vS = MODELPARAMETER(vS);
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZV = MODELPARAMETER(ZV);

   const std::complex<double> result = std::complex<double>(0,0.5)*(vS*Conj(ZA(gt1,0))*(Lambdax*SUM(j2,0,2,Conj(ZV(gt2,3 + j2))*SUM(j1,0,2,Conj(Yv(j1,j2))*ZV(gt3,j1))) - Conj(Lambdax)*SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(gt2,j1))*Yv(j1,j2))*ZV(gt3,3 + j2))) + vd*Conj(ZA(gt1,2))*(Lambdax*SUM(j2,0,2,Conj(ZV(gt2,3 + j2))*SUM(j1,0,2,Conj(Yv(j1,j2))*ZV(gt3,j1))) - Conj(Lambdax)*SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(gt2,j1))*Yv(j1,j2))*ZV(gt3,3 + j2))) + 1.4142135623730951*Conj(ZA(gt1,1))*(SUM(j2,0,2,Conj(ZV(gt2,3 + j2))*SUM(j1,0,2,Conj(TYv(j1,j2))*ZV(gt3,j1))) - SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(gt2,j1))*TYv(j1,j2))*ZV(gt3,3 + j2))));

   return {result};
}

ChiralVertex VertexImpl<fields::Ah, typename fields::bar<fields::Fe>::type, fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZER = MODELPARAMETER(ZER);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gt3,0))*SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,j2)));

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*Conj(ZA(gt3,0))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZEL(gt1,j2));

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

   const std::complex<double> result = std::complex<double>(0.,-0.3872983346207417)*g1*g2*Cos(ThetaW)*(Conj(ZA(gt1,0))*ZP(gt2,0) + Conj(ZA(gt1,1))*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Ah, typename fields::conj<fields::Hpm>::type, fields::VWm, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = std::complex<double>(0,-0.5)*g2*(Conj(ZA(gt1,0))*(-2*gp*QHd*Cos(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt2,0) + Conj(ZA(gt1,1))*(2*gp*QHu*Cos(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Ah, typename fields::conj<fields::Hpm>::type, fields::VWm, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = std::complex<double>(0,0.5)*g2*(Conj(ZA(gt1,0))*(0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHd*Sin(ThetaWp))*ZP(gt2,0) + Conj(ZA(gt1,1))*(0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*QHu*Sin(ThetaWp))*ZP(gt2,1));

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

   const std::complex<double> result = std::complex<double>(0,-0.5)*g2*(Conj(ZA(gt1,0))*ZP(gt2,0) + Conj(ZA(gt1,1))*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
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

ChiralVertex VertexImpl<fields::Cha, fields::Fv, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Yv = MODELPARAMETER(Yv);
   const auto UM = MODELPARAMETER(UM);
   const auto ZVL = MODELPARAMETER(ZVL);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZVR = MODELPARAMETER(ZVR);
   const auto UP = MODELPARAMETER(UP);

   const std::complex<double> left = -(g2*Conj(UM(gt1,0))*SUM(j1,0,2,Conj(ZVL(gt2,j1))*ZE(gt3,j1))) + Conj(UM(gt1,1))*SUM(j2,0,2,Conj(ZVL(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)));

   const std::complex<double> right = SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j2))*ZE(gt3,j1))*ZVR(gt2,j2))*UP(gt1,1);

   return {left, right};
}

ChiralVertex VertexImpl<fields::Chi, fields::Cha, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto UM = MODELPARAMETER(UM);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZP = MODELPARAMETER(ZP);
   const auto UP = MODELPARAMETER(UP);

   const std::complex<double> left = -(g2*Conj(UM(gt2,0))*Conj(ZN(gt1,3))*ZP(gt3,0)) + Conj(UM(gt2,1))*(-1.4142135623730951*gp*QHd*Conj(ZN(gt1,0))*ZP(gt3,0) + 0.5477225575051661*g1*Conj(ZN(gt1,1))*ZP(gt3,0) + 0.7071067811865475*g2*Conj(ZN(gt1,2))*ZP(gt3,0) - Conj(ZN(gt1,5))*Lambdax*ZP(gt3,1));

   const std::complex<double> right = -(Conj(Lambdax)*UP(gt2,1)*ZN(gt1,5)*ZP(gt3,0)) - 0.5*(1.4142135623730951*UP(gt2,1)*(2*gp*QHu*ZN(gt1,0) + 0.7745966692414834*g1*ZN(gt1,1) + g2*ZN(gt1,2)) + 2*g2*UP(gt2,0)*ZN(gt1,4))*ZP(gt3,1);

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

   const std::complex<double> left = -0.5*g2*(2*Conj(UM(gt2,0))*ZN(gt1,2) + 1.4142135623730951*Conj(UM(gt2,1))*ZN(gt1,3));

   const std::complex<double> right = -(g2*Conj(ZN(gt1,2))*UP(gt2,0)) + 0.7071067811865475*g2*Conj(ZN(gt1,4))*UP(gt2,1);

   return {left, right};
}

ChiralVertex VertexImpl<fields::Chi, fields::Chi, fields::Ah>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto gp = MODELPARAMETER(gp);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZN = MODELPARAMETER(ZN);

   const std::complex<double> left = std::complex<double>(0,-0.5)*(-(Conj(ZA(gt3,2))*(2*gp*Qs*Conj(ZN(gt1,5))*Conj(ZN(gt2,0)) + 2*gp*Qs*Conj(ZN(gt1,0))*Conj(ZN(gt2,5)) + 1.4142135623730951*Conj(ZN(gt1,4))*Conj(ZN(gt2,3))*Lambdax + 1.4142135623730951*Conj(ZN(gt1,3))*Conj(ZN(gt2,4))*Lambdax)) - Conj(ZA(gt3,1))*(Conj(ZN(gt1,4))*(2*gp*QHu*Conj(ZN(gt2,0)) + 0.7745966692414834*g1*Conj(ZN(gt2,1)) - g2*Conj(ZN(gt2,2))) + 2*gp*QHu*Conj(ZN(gt1,0))*Conj(ZN(gt2,4)) + 0.7745966692414834*g1*Conj(ZN(gt1,1))*Conj(ZN(gt2,4)) - g2*Conj(ZN(gt1,2))*Conj(ZN(gt2,4)) + 1.4142135623730951*Conj(ZN(gt1,5))*Conj(ZN(gt2,3))*Lambdax + 1.4142135623730951*Conj(ZN(gt1,3))*Conj(ZN(gt2,5))*Lambdax) - Conj(ZA(gt3,0))*(Conj(ZN(gt1,3))*(2*gp*QHd*Conj(ZN(gt2,0)) - 0.7745966692414834*g1*Conj(ZN(gt2,1)) + g2*Conj(ZN(gt2,2))) + 2*gp*QHd*Conj(ZN(gt1,0))*Conj(ZN(gt2,3)) - 0.7745966692414834*g1*Conj(ZN(gt1,1))*Conj(ZN(gt2,3)) + g2*Conj(ZN(gt1,2))*Conj(ZN(gt2,3)) + 1.4142135623730951*Conj(ZN(gt1,5))*Conj(ZN(gt2,4))*Lambdax + 1.4142135623730951*Conj(ZN(gt1,4))*Conj(ZN(gt2,5))*Lambdax));

   const std::complex<double> right = std::complex<double>(0,-0.5)*(Conj(ZA(gt3,2))*(2*gp*Qs*ZN(gt1,5)*ZN(gt2,0) + 1.4142135623730951*Conj(Lambdax)*(ZN(gt1,4)*ZN(gt2,3) + ZN(gt1,3)*ZN(gt2,4)) + 2*gp*Qs*ZN(gt1,0)*ZN(gt2,5)) + Conj(ZA(gt3,0))*(ZN(gt1,3)*(2*gp*QHd*ZN(gt2,0) - 0.7745966692414834*g1*ZN(gt2,1) + g2*ZN(gt2,2)) + 2*gp*QHd*ZN(gt1,0)*ZN(gt2,3) - 0.7745966692414834*g1*ZN(gt1,1)*ZN(gt2,3) + g2*ZN(gt1,2)*ZN(gt2,3) + 1.4142135623730951*Conj(Lambdax)*ZN(gt1,5)*ZN(gt2,4) + 1.4142135623730951*Conj(Lambdax)*ZN(gt1,4)*ZN(gt2,5)) + Conj(ZA(gt3,1))*(ZN(gt1,4)*(2*gp*QHu*ZN(gt2,0) + 0.7745966692414834*g1*ZN(gt2,1) - g2*ZN(gt2,2)) + (2*gp*QHu*ZN(gt1,0) + 0.7745966692414834*g1*ZN(gt1,1) - g2*ZN(gt1,2))*ZN(gt2,4) + 1.4142135623730951*Conj(Lambdax)*(ZN(gt1,5)*ZN(gt2,3) + ZN(gt1,3)*ZN(gt2,5))));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Chi, fields::Chi, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto gp = MODELPARAMETER(gp);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZN = MODELPARAMETER(ZN);

   const std::complex<double> left = 0.1*(5*Conj(ZH(gt3,2))*(-2*gp*Qs*Conj(ZN(gt1,5))*Conj(ZN(gt2,0)) - 2*gp*Qs*Conj(ZN(gt1,0))*Conj(ZN(gt2,5)) + 1.4142135623730951*Conj(ZN(gt1,4))*Conj(ZN(gt2,3))*Lambdax + 1.4142135623730951*Conj(ZN(gt1,3))*Conj(ZN(gt2,4))*Lambdax) - Conj(ZH(gt3,1))*(Conj(ZN(gt1,4))*(10*gp*QHu*Conj(ZN(gt2,0)) + 3.872983346207417*g1*Conj(ZN(gt2,1)) - 5*g2*Conj(ZN(gt2,2))) + 10*gp*QHu*Conj(ZN(gt1,0))*Conj(ZN(gt2,4)) + 3.872983346207417*g1*Conj(ZN(gt1,1))*Conj(ZN(gt2,4)) - 5*g2*Conj(ZN(gt1,2))*Conj(ZN(gt2,4)) - 7.0710678118654755*Conj(ZN(gt1,5))*Conj(ZN(gt2,3))*Lambdax - 7.0710678118654755*Conj(ZN(gt1,3))*Conj(ZN(gt2,5))*Lambdax) + Conj(ZH(gt3,0))*(Conj(ZN(gt1,3))*(-10*gp*QHd*Conj(ZN(gt2,0)) + 3.872983346207417*g1*Conj(ZN(gt2,1)) - 5*g2*Conj(ZN(gt2,2))) - 10*gp*QHd*Conj(ZN(gt1,0))*Conj(ZN(gt2,3)) + 3.872983346207417*g1*Conj(ZN(gt1,1))*Conj(ZN(gt2,3)) - 5*g2*Conj(ZN(gt1,2))*Conj(ZN(gt2,3)) + 7.0710678118654755*Conj(ZN(gt1,5))*Conj(ZN(gt2,4))*Lambdax + 7.0710678118654755*Conj(ZN(gt1,4))*Conj(ZN(gt2,5))*Lambdax));

   const std::complex<double> right = 0.1*(5*Conj(ZH(gt3,2))*(-2*gp*Qs*ZN(gt1,5)*ZN(gt2,0) + 1.4142135623730951*Conj(Lambdax)*(ZN(gt1,4)*ZN(gt2,3) + ZN(gt1,3)*ZN(gt2,4)) - 2*gp*Qs*ZN(gt1,0)*ZN(gt2,5)) + Conj(ZH(gt3,0))*(ZN(gt1,3)*(-10*gp*QHd*ZN(gt2,0) + 3.872983346207417*g1*ZN(gt2,1) - 5*g2*ZN(gt2,2)) - 10*gp*QHd*ZN(gt1,0)*ZN(gt2,3) + 3.872983346207417*g1*ZN(gt1,1)*ZN(gt2,3) - 5*g2*ZN(gt1,2)*ZN(gt2,3) + 7.0710678118654755*Conj(Lambdax)*ZN(gt1,5)*ZN(gt2,4) + 7.0710678118654755*Conj(Lambdax)*ZN(gt1,4)*ZN(gt2,5)) - Conj(ZH(gt3,1))*(ZN(gt1,4)*(10*gp*QHu*ZN(gt2,0) + 3.872983346207417*g1*ZN(gt2,1) - 5*g2*ZN(gt2,2)) + (10*gp*QHu*ZN(gt1,0) + 3.872983346207417*g1*ZN(gt1,1) - 5*g2*ZN(gt1,2))*ZN(gt2,4) - 7.0710678118654755*Conj(Lambdax)*(ZN(gt1,5)*ZN(gt2,3) + ZN(gt1,3)*ZN(gt2,5))));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Chi, fields::Chi, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto Qs = INPUTPARAMETER(Qs);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.5*(-(Conj(ZN(gt2,3))*(2*gp*QHd*Cos(ThetaWp) - (g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*Sin(ThetaWp))*ZN(gt1,3)) - Conj(ZN(gt2,4))*(2*gp*QHu*Cos(ThetaWp) + g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZN(gt1,4) - 2*gp*Qs*Conj(ZN(gt2,5))*Cos(ThetaWp)*ZN(gt1,5));

   const std::complex<double> right = 0.5*(Conj(ZN(gt1,3))*(2*gp*QHd*Cos(ThetaWp) - (g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*Sin(ThetaWp))*ZN(gt2,3) + Conj(ZN(gt1,4))*(2*gp*QHu*Cos(ThetaWp) + g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZN(gt2,4) + 2*gp*Qs*Conj(ZN(gt1,5))*Cos(ThetaWp)*ZN(gt2,5));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Chi, fields::Chi, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto Qs = INPUTPARAMETER(Qs);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = 0.5*(-(Conj(ZN(gt2,3))*(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHd*Sin(ThetaWp))*ZN(gt1,3)) + Conj(ZN(gt2,4))*(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*QHu*Sin(ThetaWp))*ZN(gt1,4) - 2*gp*Qs*Conj(ZN(gt2,5))*Sin(ThetaWp)*ZN(gt1,5));

   const std::complex<double> right = 0.5*(Conj(ZN(gt1,3))*(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHd*Sin(ThetaWp))*ZN(gt2,3) - Conj(ZN(gt1,4))*(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*QHu*Sin(ThetaWp))*ZN(gt2,4) + 2*gp*Qs*Conj(ZN(gt1,5))*Sin(ThetaWp)*ZN(gt2,5));

   return {left, right};
}

ChiralVertex VertexImpl<fields::Chi, fields::Fd, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qd = INPUTPARAMETER(Qd);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZDR = MODELPARAMETER(ZDR);

   const std::complex<double> left = -1.4142135623730951*gp*Qq*Conj(ZN(gt1,0))*SUM(j1,0,2,Conj(ZDL(gt2,j1))*ZD(gt3,j1)) - 0.18257418583505536*g1*Conj(ZN(gt1,1))*SUM(j1,0,2,Conj(ZDL(gt2,j1))*ZD(gt3,j1)) + 0.7071067811865475*g2*Conj(ZN(gt1,2))*SUM(j1,0,2,Conj(ZDL(gt2,j1))*ZD(gt3,j1)) - Conj(ZN(gt1,3))*SUM(j2,0,2,Conj(ZDL(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1)));

   const std::complex<double> right = SUM(j1,0,2,ZD(gt3,3 + j1)*ZDR(gt2,j1))*(-1.4142135623730951*gp*Qd*ZN(gt1,0) - 0.3651483716701107*g1*ZN(gt1,1)) - SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*ZDR(gt2,j1))*ZD(gt3,j2))*ZN(gt1,3);

   return {left, right};
}

ChiralVertex VertexImpl<fields::Chi, fields::Fe, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qe = INPUTPARAMETER(Qe);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZER = MODELPARAMETER(ZER);

   const std::complex<double> left = -1.4142135623730951*gp*Ql*Conj(ZN(gt1,0))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZE(gt3,j1)) + 0.5477225575051661*g1*Conj(ZN(gt1,1))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZE(gt3,j1)) + 0.7071067811865475*g2*Conj(ZN(gt1,2))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZE(gt3,j1)) - Conj(ZN(gt1,3))*SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)));

   const std::complex<double> right = -1.4142135623730951*SUM(j1,0,2,ZE(gt3,3 + j1)*ZER(gt2,j1))*(gp*Qe*ZN(gt1,0) + 0.7745966692414834*g1*ZN(gt1,1)) - SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZE(gt3,j2))*ZN(gt1,3);

   return {left, right};
}

ChiralVertex VertexImpl<fields::Chi, fields::Fu, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qu = INPUTPARAMETER(Qu);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZUR = MODELPARAMETER(ZUR);

   const std::complex<double> left = -1.4142135623730951*gp*Qq*Conj(ZN(gt1,0))*SUM(j1,0,2,Conj(ZUL(gt2,j1))*ZU(gt3,j1)) - 0.18257418583505536*g1*Conj(ZN(gt1,1))*SUM(j1,0,2,Conj(ZUL(gt2,j1))*ZU(gt3,j1)) - 0.7071067811865475*g2*Conj(ZN(gt1,2))*SUM(j1,0,2,Conj(ZUL(gt2,j1))*ZU(gt3,j1)) - Conj(ZN(gt1,4))*SUM(j2,0,2,Conj(ZUL(gt2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt3,3 + j1)));

   const std::complex<double> right = SUM(j1,0,2,ZU(gt3,3 + j1)*ZUR(gt2,j1))*(-1.4142135623730951*gp*Qu*ZN(gt1,0) + 0.7302967433402214*g1*ZN(gt1,1)) - SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*ZUR(gt2,j1))*ZU(gt3,j2))*ZN(gt1,4);

   return {left, right};
}

ChiralVertex VertexImpl<fields::Chi, fields::Fv, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qv = INPUTPARAMETER(Qv);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZVL = MODELPARAMETER(ZVL);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZVR = MODELPARAMETER(ZVR);

   const std::complex<double> left = -1.4142135623730951*gp*Ql*Conj(ZN(gt1,0))*SUM(j1,0,2,Conj(ZVL(gt2,j1))*ZV(gt3,j1)) + 0.5477225575051661*g1*Conj(ZN(gt1,1))*SUM(j1,0,2,Conj(ZVL(gt2,j1))*ZV(gt3,j1)) - 0.7071067811865475*g2*Conj(ZN(gt1,2))*SUM(j1,0,2,Conj(ZVL(gt2,j1))*ZV(gt3,j1)) - Conj(ZN(gt1,4))*SUM(j2,0,2,SUM(j1,0,2,Conj(ZVL(gt2,j1))*Yv(j1,j2))*ZV(gt3,3 + j2));

   const std::complex<double> right = -1.4142135623730951*gp*Qv*SUM(j1,0,2,ZV(gt3,3 + j1)*ZVR(gt2,j1))*ZN(gt1,0) - SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j2))*ZV(gt3,j1))*ZVR(gt2,j2))*ZN(gt1,4);

   return {left, right};
}

ChiralVertex VertexImpl<fields::Chi, typename fields::conj<fields::Se>::type, fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qe = INPUTPARAMETER(Qe);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZER = MODELPARAMETER(ZER);

   const std::complex<double> left = -1.4142135623730951*gp*Ql*Conj(ZN(gt1,0))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZE(gt3,j1)) + 0.5477225575051661*g1*Conj(ZN(gt1,1))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZE(gt3,j1)) + 0.7071067811865475*g2*Conj(ZN(gt1,2))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZE(gt3,j1)) - Conj(ZN(gt1,3))*SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)));

   const std::complex<double> right = -1.4142135623730951*SUM(j1,0,2,ZE(gt3,3 + j1)*ZER(gt2,j1))*(gp*Qe*ZN(gt1,0) + 0.7745966692414834*g1*ZN(gt1,1)) - SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZE(gt3,j2))*ZN(gt1,3);

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
   const auto QHd = INPUTPARAMETER(QHd);
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.25*(-(Conj(ZH(gt1,0))*(-(Conj(ZH(gt2,1))*(Conj(ZH(gt3,1))*Conj(ZH(gt4,0)) + Conj(ZH(gt3,0))*Conj(ZH(gt4,1)))*(-4*AbsSqr(Lambdax) + 0.6*Sqr(g1) + Sqr(g2) - 4*QHd*QHu*Sqr(gp))) + 4*Conj(ZH(gt2,2))*(Conj(ZH(gt3,2))*Conj(ZH(gt4,0)) + Conj(ZH(gt3,0))*Conj(ZH(gt4,2)))*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + Conj(ZH(gt2,0))*(-(Conj(ZH(gt3,1))*Conj(ZH(gt4,1))*(-4*AbsSqr(Lambdax) + 0.6*Sqr(g1) + Sqr(g2) - 4*QHd*QHu*Sqr(gp))) + 4*Conj(ZH(gt3,2))*Conj(ZH(gt4,2))*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + 3*Conj(ZH(gt3,0))*Conj(ZH(gt4,0))*(0.6*Sqr(g1) + Sqr(g2) + 4*Sqr(gp)*Sqr(QHd))))) - Conj(ZH(gt1,1))*(-(Conj(ZH(gt2,0))*(Conj(ZH(gt3,1))*Conj(ZH(gt4,0)) + Conj(ZH(gt3,0))*Conj(ZH(gt4,1)))*(-4*AbsSqr(Lambdax) + 0.6*Sqr(g1) + Sqr(g2) - 4*QHd*QHu*Sqr(gp))) + 4*Conj(ZH(gt2,2))*(Conj(ZH(gt3,2))*Conj(ZH(gt4,1)) + Conj(ZH(gt3,1))*Conj(ZH(gt4,2)))*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) + Conj(ZH(gt2,1))*(-(Conj(ZH(gt3,0))*Conj(ZH(gt4,0))*(-4*AbsSqr(Lambdax) + 0.6*Sqr(g1) + Sqr(g2) - 4*QHd*QHu*Sqr(gp))) + 4*Conj(ZH(gt3,2))*Conj(ZH(gt4,2))*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) + 3*Conj(ZH(gt3,1))*Conj(ZH(gt4,1))*(0.6*Sqr(g1) + Sqr(g2) + 4*Sqr(gp)*Sqr(QHu)))) - 4*Conj(ZH(gt1,2))*(Conj(ZH(gt2,0))*(Conj(ZH(gt3,2))*Conj(ZH(gt4,0)) + Conj(ZH(gt3,0))*Conj(ZH(gt4,2)))*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + Conj(ZH(gt2,1))*(Conj(ZH(gt3,2))*Conj(ZH(gt4,1)) + Conj(ZH(gt3,1))*Conj(ZH(gt4,2)))*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) + Conj(ZH(gt2,2))*(Conj(ZH(gt3,0))*Conj(ZH(gt4,0))*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + Conj(ZH(gt3,1))*Conj(ZH(gt4,1))*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) + 3*Conj(ZH(gt3,2))*Conj(ZH(gt4,2))*Sqr(gp)*Sqr(Qs))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto TLambdax = MODELPARAMETER(TLambdax);
   const auto gp = MODELPARAMETER(gp);
   const auto vd = MODELPARAMETER(vd);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto vS = MODELPARAMETER(vS);
   const auto vu = MODELPARAMETER(vu);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.25*(-(Conj(ZH(gt1,2))*(-1.4142135623730951*Conj(TLambdax)*Conj(ZH(gt2,1))*Conj(ZH(gt3,0)) + 4*vd*AbsSqr(Lambdax)*Conj(ZH(gt2,2))*Conj(ZH(gt3,0)) + 4*vS*AbsSqr(Lambdax)*Conj(ZH(gt2,1))*Conj(ZH(gt3,1)) + 4*vu*AbsSqr(Lambdax)*Conj(ZH(gt2,2))*Conj(ZH(gt3,1)) + 4*vu*AbsSqr(Lambdax)*Conj(ZH(gt2,1))*Conj(ZH(gt3,2)) + 4*QHd*Qs*vd*Conj(ZH(gt2,2))*Conj(ZH(gt3,0))*Sqr(gp) + 4*QHu*Qs*vS*Conj(ZH(gt2,1))*Conj(ZH(gt3,1))*Sqr(gp) + 4*QHu*Qs*vu*Conj(ZH(gt2,2))*Conj(ZH(gt3,1))*Sqr(gp) + 4*QHu*Qs*vu*Conj(ZH(gt2,1))*Conj(ZH(gt3,2))*Sqr(gp) + 12*vS*Conj(ZH(gt2,2))*Conj(ZH(gt3,2))*Sqr(gp)*Sqr(Qs) - 1.4142135623730951*Conj(ZH(gt2,1))*Conj(ZH(gt3,0))*TLambdax + Conj(ZH(gt2,0))*(-1.4142135623730951*Conj(TLambdax)*Conj(ZH(gt3,1)) + 4*vd*AbsSqr(Lambdax)*Conj(ZH(gt3,2)) + 4*QHd*Qs*vd*Conj(ZH(gt3,2))*Sqr(gp) + 4*vS*Conj(ZH(gt3,0))*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) - 1.4142135623730951*Conj(ZH(gt3,1))*TLambdax))) + Conj(ZH(gt1,1))*(Conj(ZH(gt2,1))*(vd*Conj(ZH(gt3,0))*(-4*AbsSqr(Lambdax) + 0.6*Sqr(g1) + Sqr(g2) - 4*QHd*QHu*Sqr(gp)) - 4*vS*Conj(ZH(gt3,2))*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) - 3*vu*Conj(ZH(gt3,1))*(0.6*Sqr(g1) + Sqr(g2) + 4*Sqr(gp)*Sqr(QHu))) + Conj(ZH(gt2,2))*(1.4142135623730951*Conj(TLambdax)*Conj(ZH(gt3,0)) - 4*vu*AbsSqr(Lambdax)*Conj(ZH(gt3,2)) - 4*QHu*Qs*vu*Conj(ZH(gt3,2))*Sqr(gp) - 4*vS*Conj(ZH(gt3,1))*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) + 1.4142135623730951*Conj(ZH(gt3,0))*TLambdax) + Conj(ZH(gt2,0))*(vu*Conj(ZH(gt3,0))*(-4*AbsSqr(Lambdax) + 0.6*Sqr(g1) + Sqr(g2) - 4*QHd*QHu*Sqr(gp)) + vd*Conj(ZH(gt3,1))*(-4*AbsSqr(Lambdax) + 0.6*Sqr(g1) + Sqr(g2) - 4*QHd*QHu*Sqr(gp)) + 1.4142135623730951*Conj(ZH(gt3,2))*(Conj(TLambdax) + TLambdax))) + Conj(ZH(gt1,0))*(-(Conj(ZH(gt2,0))*(-(vu*Conj(ZH(gt3,1))*(-4*AbsSqr(Lambdax) + 0.6*Sqr(g1) + Sqr(g2) - 4*QHd*QHu*Sqr(gp))) + 4*vS*Conj(ZH(gt3,2))*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + 3*vd*Conj(ZH(gt3,0))*(0.6*Sqr(g1) + Sqr(g2) + 4*Sqr(gp)*Sqr(QHd)))) + Conj(ZH(gt2,2))*(1.4142135623730951*Conj(TLambdax)*Conj(ZH(gt3,1)) - 4*vd*AbsSqr(Lambdax)*Conj(ZH(gt3,2)) - 4*QHd*Qs*vd*Conj(ZH(gt3,2))*Sqr(gp) - 4*vS*Conj(ZH(gt3,0))*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + 1.4142135623730951*Conj(ZH(gt3,1))*TLambdax) + Conj(ZH(gt2,1))*(vu*Conj(ZH(gt3,0))*(-4*AbsSqr(Lambdax) + 0.6*Sqr(g1) + Sqr(g2) - 4*QHd*QHu*Sqr(gp)) + vd*Conj(ZH(gt3,1))*(-4*AbsSqr(Lambdax) + 0.6*Sqr(g1) + Sqr(g2) - 4*QHd*QHu*Sqr(gp)) + 1.4142135623730951*Conj(ZH(gt3,2))*(Conj(TLambdax) + TLambdax))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::Hpm, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-4*Conj(ZH(gt1,2))*Conj(ZH(gt2,2))*((AbsSqr(Lambdax) + QHd*Qs*Sqr(gp))*ZP(gt3,0)*ZP(gt4,0) + (AbsSqr(Lambdax) + QHu*Qs*Sqr(gp))*ZP(gt3,1)*ZP(gt4,1)) - Conj(ZH(gt1,0))*(Conj(ZH(gt2,1))*(-2*AbsSqr(Lambdax) + Sqr(g2))*(ZP(gt3,1)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,1)) + Conj(ZH(gt2,0))*((0.6*Sqr(g1) + Sqr(g2) + 4*Sqr(gp)*Sqr(QHd))*ZP(gt3,0)*ZP(gt4,0) + (-0.6*Sqr(g1) + Sqr(g2) + 4*QHd*QHu*Sqr(gp))*ZP(gt3,1)*ZP(gt4,1))) - Conj(ZH(gt1,1))*(Conj(ZH(gt2,0))*(-2*AbsSqr(Lambdax) + Sqr(g2))*(ZP(gt3,1)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,1)) + Conj(ZH(gt2,1))*((-0.6*Sqr(g1) + Sqr(g2) + 4*QHd*QHu*Sqr(gp))*ZP(gt3,0)*ZP(gt4,0) + (0.6*Sqr(g1) + Sqr(g2) + 4*Sqr(gp)*Sqr(QHu))*ZP(gt3,1)*ZP(gt4,1))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::Sd, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qd = INPUTPARAMETER(Qd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto gp = MODELPARAMETER(gp);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = 0.05*(10*Conj(ZH(gt1,2))*(-2*Qs*Conj(ZH(gt2,2))*Sqr(gp)*(Qq*SUM(j1,0,2,Conj(ZD(gt3,j1))*ZD(gt4,j1)) + Qd*SUM(j1,0,2,Conj(ZD(gt3,3 + j1))*ZD(gt4,3 + j1))) + Conj(ZH(gt2,1))*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZD(gt3,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt4,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt3,3 + j1)))*ZD(gt4,j2)))) - Conj(ZH(gt1,1))*(Conj(ZH(gt2,1))*((Sqr(g1) + 5*Sqr(g2) + 20*QHu*Qq*Sqr(gp))*SUM(j1,0,2,Conj(ZD(gt3,j1))*ZD(gt4,j1)) + 2*(Sqr(g1) + 10*Qd*QHu*Sqr(gp))*SUM(j1,0,2,Conj(ZD(gt3,3 + j1))*ZD(gt4,3 + j1))) - 10*Conj(ZH(gt2,2))*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZD(gt3,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt4,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt3,3 + j1)))*ZD(gt4,j2)))) + Conj(ZH(gt1,0))*Conj(ZH(gt2,0))*((Sqr(g1) + 5*(Sqr(g2) - 4*QHd*Qq*Sqr(gp)))*SUM(j1,0,2,Conj(ZD(gt3,j1))*ZD(gt4,j1)) + 2*(Sqr(g1) - 10*Qd*QHd*Sqr(gp))*SUM(j1,0,2,Conj(ZD(gt3,3 + j1))*ZD(gt4,3 + j1)) - 20*(SUM(j3,0,2,Conj(ZD(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,j1))*ZD(gt4,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt3,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gt4,j3)))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::Se, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qe = INPUTPARAMETER(Qe);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto gp = MODELPARAMETER(gp);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZE = MODELPARAMETER(ZE);

   const std::complex<double> result = 0.05*(10*Conj(ZH(gt1,2))*(-2*Qs*Conj(ZH(gt2,2))*Sqr(gp)*(Ql*SUM(j1,0,2,Conj(ZE(gt3,j1))*ZE(gt4,j1)) + Qe*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*ZE(gt4,3 + j1))) + Conj(ZH(gt2,1))*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZE(gt3,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt4,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt3,3 + j1)))*ZE(gt4,j2)))) + Conj(ZH(gt1,1))*(Conj(ZH(gt2,1))*((3*Sqr(g1) - 5*(Sqr(g2) + 4*QHu*Ql*Sqr(gp)))*SUM(j1,0,2,Conj(ZE(gt3,j1))*ZE(gt4,j1)) - 2*(3*Sqr(g1) + 10*Qe*QHu*Sqr(gp))*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*ZE(gt4,3 + j1))) + 10*Conj(ZH(gt2,2))*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZE(gt3,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt4,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt3,3 + j1)))*ZE(gt4,j2)))) - Conj(ZH(gt1,0))*Conj(ZH(gt2,0))*((3*Sqr(g1) - 5*Sqr(g2) + 20*QHd*Ql*Sqr(gp))*SUM(j1,0,2,Conj(ZE(gt3,j1))*ZE(gt4,j1)) + (-6*Sqr(g1) + 20*Qe*QHd*Sqr(gp))*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*ZE(gt4,3 + j1)) + 20*(SUM(j3,0,2,Conj(ZE(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Ye(j2,j1))*ZE(gt4,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gt3,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZE(gt4,j3)))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::Su, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qu = INPUTPARAMETER(Qu);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZU = MODELPARAMETER(ZU);

   const std::complex<double> result = 0.05*(10*Conj(ZH(gt1,2))*(-2*Qs*Conj(ZH(gt2,2))*Sqr(gp)*(Qq*SUM(j1,0,2,Conj(ZU(gt3,j1))*ZU(gt4,j1)) + Qu*SUM(j1,0,2,Conj(ZU(gt3,3 + j1))*ZU(gt4,3 + j1))) + Conj(ZH(gt2,0))*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt4,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt3,3 + j1)))*ZU(gt4,j2)))) + Conj(ZH(gt1,0))*(Conj(ZH(gt2,0))*((Sqr(g1) - 5*(Sqr(g2) + 4*QHd*Qq*Sqr(gp)))*SUM(j1,0,2,Conj(ZU(gt3,j1))*ZU(gt4,j1)) - 4*(Sqr(g1) + 5*QHd*Qu*Sqr(gp))*SUM(j1,0,2,Conj(ZU(gt3,3 + j1))*ZU(gt4,3 + j1))) + 10*Conj(ZH(gt2,2))*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt4,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt3,3 + j1)))*ZU(gt4,j2)))) - Conj(ZH(gt1,1))*Conj(ZH(gt2,1))*((Sqr(g1) - 5*Sqr(g2) + 20*QHu*Qq*Sqr(gp))*SUM(j1,0,2,Conj(ZU(gt3,j1))*ZU(gt4,j1)) - 4*(Sqr(g1) - 5*QHu*Qu*Sqr(gp))*SUM(j1,0,2,Conj(ZU(gt3,3 + j1))*ZU(gt4,3 + j1)) + 20*(SUM(j3,0,2,Conj(ZU(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))*ZU(gt4,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gt4,j3)))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::Sv, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qv = INPUTPARAMETER(Qv);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZV = MODELPARAMETER(ZV);

   const std::complex<double> result = 0.25*(-2*Conj(ZH(gt1,2))*(2*Qs*Conj(ZH(gt2,2))*Sqr(gp)*(Ql*SUM(j1,0,2,Conj(ZV(gt3,j1))*ZV(gt4,j1)) + Qv*SUM(j1,0,2,Conj(ZV(gt3,3 + j1))*ZV(gt4,3 + j1))) - Conj(ZH(gt2,0))*(Lambdax*SUM(j2,0,2,Conj(ZV(gt3,3 + j2))*SUM(j1,0,2,Conj(Yv(j1,j2))*ZV(gt4,j1))) + Conj(Lambdax)*SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(gt3,j1))*Yv(j1,j2))*ZV(gt4,3 + j2)))) - Conj(ZH(gt1,0))*(Conj(ZH(gt2,0))*((0.6*Sqr(g1) + Sqr(g2) + 4*QHd*Ql*Sqr(gp))*SUM(j1,0,2,Conj(ZV(gt3,j1))*ZV(gt4,j1)) + 4*QHd*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gt3,3 + j1))*ZV(gt4,3 + j1))) - 2*Conj(ZH(gt2,2))*(Lambdax*SUM(j2,0,2,Conj(ZV(gt3,3 + j2))*SUM(j1,0,2,Conj(Yv(j1,j2))*ZV(gt4,j1))) + Conj(Lambdax)*SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(gt3,j1))*Yv(j1,j2))*ZV(gt4,3 + j2)))) - Conj(ZH(gt1,1))*Conj(ZH(gt2,1))*(-((0.6*Sqr(g1) + Sqr(g2) - 4*QHu*Ql*Sqr(gp))*SUM(j1,0,2,Conj(ZV(gt3,j1))*ZV(gt4,j1))) + 4*(QHu*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gt3,3 + j1))*ZV(gt4,3 + j1)) + SUM(j3,0,2,Conj(ZV(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j3))*Yv(j1,j2))*ZV(gt4,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gt3,j2))*SUM(j1,0,2,Conj(Yv(j3,j1))*Yv(j2,j1)))*ZV(gt4,j3)))));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::hh, fields::VZ, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.1*(20*Conj(ZH(gt1,2))*Conj(ZH(gt2,2))*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gp)*Sqr(Qs) - Conj(ZH(gt1,0))*Conj(ZH(gt2,0))*(5*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) - 7.745966692414834*g1*gp*QHd*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + Cos(ThetaWp)*Sin(ThetaWp)*(-20*Sqr(gp)*Sqr(QHd) + 3*Sqr(g1)*Sqr(Sin(ThetaW))) + 7.745966692414834*g1*gp*QHd*Sin(ThetaW)*Sqr(Sin(ThetaWp)) + 2*g2*Cos(ThetaW)*(3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) - 5*gp*QHd*Sqr(Cos(ThetaWp)) + 5*gp*QHd*Sqr(Sin(ThetaWp)))) - Conj(ZH(gt1,1))*Conj(ZH(gt2,1))*(5*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 7.745966692414834*g1*gp*QHu*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + Cos(ThetaWp)*Sin(ThetaWp)*(-20*Sqr(gp)*Sqr(QHu) + 3*Sqr(g1)*Sqr(Sin(ThetaW))) - 7.745966692414834*g1*gp*QHu*Sin(ThetaW)*Sqr(Sin(ThetaWp)) + 2*g2*Cos(ThetaW)*(3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) + 5*gp*QHu*Sqr(Cos(ThetaWp)) - 5*gp*QHu*Sqr(Sin(ThetaWp)))));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::hh, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(4*Conj(ZH(gt1,2))*Conj(ZH(gt2,2))*Sqr(gp)*Sqr(Qs)*Sqr(Sin(ThetaWp)) + Conj(ZH(gt1,0))*Conj(ZH(gt2,0))*Sqr(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHd*Sin(ThetaWp)) + Conj(ZH(gt1,1))*Conj(ZH(gt2,1))*Sqr(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*QHu*Sin(ThetaWp)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::hh, fields::VZp, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(4*Conj(ZH(gt1,2))*Conj(ZH(gt2,2))*Sqr(gp)*Sqr(Qs)*Sqr(Cos(ThetaWp)) + Conj(ZH(gt1,0))*Conj(ZH(gt2,0))*Sqr(-2*gp*QHd*Cos(ThetaWp) + g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp)) + Conj(ZH(gt1,1))*Conj(ZH(gt2,1))*Sqr(2*gp*QHu*Cos(ThetaWp) + g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::hh, typename fields::conj<fields::VWm>::type, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.5*(Conj(ZH(gt1,0))*Conj(ZH(gt2,0)) + Conj(ZH(gt1,1))*Conj(ZH(gt2,1)))*Sqr(g2);

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Hpm, fields::Su, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.35355339059327373*(2*Conj(ZH(gt1,2))*(Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt3,3 + j1)))*ZD(gt4,j2))*ZP(gt2,0) + Conj(Lambdax)*SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt4,3 + j1)))*ZP(gt2,1)) + Conj(ZH(gt1,0))*(-(Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt3,j1))*ZD(gt4,j1))*ZP(gt2,0)) + 2*SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gt4,j3))*ZP(gt2,0) + 2*SUM(j3,0,2,Conj(ZU(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yd(j2,j1))*ZD(gt4,3 + j2)))*ZP(gt2,1)) + Conj(ZH(gt1,1))*(2*SUM(j3,0,2,Conj(ZU(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yd(j2,j1))*ZD(gt4,3 + j2)))*ZP(gt2,0) + (-(Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt3,j1))*ZD(gt4,j1))) + 2*SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZD(gt4,j3)))*ZP(gt2,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Hpm, fields::Sv, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yv = MODELPARAMETER(Yv);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.35355339059327373*(2*Conj(ZH(gt1,2))*(Lambdax*SUM(j2,0,2,Conj(ZV(gt3,3 + j2))*SUM(j1,0,2,Conj(Yv(j1,j2))*ZE(gt4,j1)))*ZP(gt2,0) + Conj(Lambdax)*SUM(j2,0,2,Conj(ZV(gt3,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt4,3 + j1)))*ZP(gt2,1)) + Conj(ZH(gt1,0))*(-(Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt3,j1))*ZE(gt4,j1))*ZP(gt2,0)) + 2*SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gt3,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZE(gt4,j3))*ZP(gt2,0) + 2*SUM(j3,0,2,Conj(ZV(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j3))*Ye(j2,j1))*ZE(gt4,3 + j2)))*ZP(gt2,1)) + Conj(ZH(gt1,1))*(2*SUM(j3,0,2,Conj(ZV(gt3,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j3))*Ye(j2,j1))*ZE(gt4,3 + j2)))*ZP(gt2,0) + (-(Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt3,j1))*ZE(gt4,j1))) + 2*SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gt3,j2))*SUM(j1,0,2,Conj(Yv(j3,j1))*Yv(j2,j1)))*ZE(gt4,j3)))*ZP(gt2,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Hpm, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto QHu = INPUTPARAMETER(QHu);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto Qs = INPUTPARAMETER(Qs);
   const auto TLambdax = MODELPARAMETER(TLambdax);
   const auto vd = MODELPARAMETER(vd);
   const auto g2 = MODELPARAMETER(g2);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto vu = MODELPARAMETER(vu);
   const auto vS = MODELPARAMETER(vS);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-(Conj(ZH(gt1,0))*(ZP(gt2,0)*(vd*(0.6*Sqr(g1) + Sqr(g2) + 4*Sqr(gp)*Sqr(QHd))*ZP(gt3,0) + vu*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gt3,1)) + ZP(gt2,1)*(vu*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gt3,0) + vd*(-0.6*Sqr(g1) + Sqr(g2) + 4*QHd*QHu*Sqr(gp))*ZP(gt3,1)))) - Conj(ZH(gt1,1))*(ZP(gt2,0)*(vu*(-0.6*Sqr(g1) + Sqr(g2) + 4*QHd*QHu*Sqr(gp))*ZP(gt3,0) + vd*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gt3,1)) + ZP(gt2,1)*(vd*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gt3,0) + vu*(0.6*Sqr(g1) + Sqr(g2) + 4*Sqr(gp)*Sqr(QHu))*ZP(gt3,1))) - 2*Conj(ZH(gt1,2))*(ZP(gt2,1)*(1.4142135623730951*Conj(TLambdax)*ZP(gt3,0) + 2*vS*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp))*ZP(gt3,1)) + ZP(gt2,0)*(2*vS*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp))*ZP(gt3,0) + 1.4142135623730951*TLambdax*ZP(gt3,1))));

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

   const std::complex<double> result = -0.3872983346207417*g1*g2*Cos(ThetaW)*(Conj(ZH(gt1,0))*ZP(gt2,0) - Conj(ZH(gt1,1))*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*g2*(Conj(ZH(gt1,0))*(2*gp*QHd*Cos(ThetaWp) - 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt2,0) + Conj(ZH(gt1,1))*(2*gp*QHu*Cos(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*g2*(Conj(ZH(gt1,0))*(0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHd*Sin(ThetaWp))*ZP(gt2,0) + Conj(ZH(gt1,1))*(-0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHu*Sin(ThetaWp))*ZP(gt2,1));

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

   const std::complex<double> result = 0.5*g2*(Conj(ZH(gt1,0))*ZP(gt2,0) - Conj(ZH(gt1,1))*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

ScalarVertex VertexImpl<fields::hh, fields::Sd, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.35355339059327373*(2*Conj(Lambdax)*Conj(ZH(gt1,2))*SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt4,3 + j1)))*ZP(gt3,0) + 2*Conj(ZH(gt1,1))*SUM(j3,0,2,Conj(ZD(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yu(j2,j1))*ZU(gt4,3 + j2)))*ZP(gt3,0) - Conj(ZH(gt1,1))*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZU(gt4,j1))*ZP(gt3,1) + 2*Conj(ZH(gt1,2))*Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1)))*ZU(gt4,j2))*ZP(gt3,1) + 2*Conj(ZH(gt1,1))*SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gt4,j3))*ZP(gt3,1) + Conj(ZH(gt1,0))*(-(Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZU(gt4,j1))*ZP(gt3,0)) + 2*SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZU(gt4,j3))*ZP(gt3,0) + 2*SUM(j3,0,2,Conj(ZD(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yu(j2,j1))*ZU(gt4,3 + j2)))*ZP(gt3,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Sd, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qs = INPUTPARAMETER(Qs);
   const auto Qd = INPUTPARAMETER(Qd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto TYd = MODELPARAMETER(TYd);
   const auto gp = MODELPARAMETER(gp);
   const auto vS = MODELPARAMETER(vS);
   const auto vu = MODELPARAMETER(vu);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = 0.05*(10*Conj(ZH(gt1,2))*(-2*Qq*Qs*vS*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1)) - 2*Qd*Qs*vS*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1)) + vu*Conj(Lambdax)*SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1))) + vu*Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1)))*ZD(gt3,j2))) - Conj(ZH(gt1,1))*(vu*(Sqr(g1) + 5*(Sqr(g2) + 4*QHu*Qq*Sqr(gp)))*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1)) + 2*(vu*(Sqr(g1) + 10*Qd*QHu*Sqr(gp))*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1)) - 5*vS*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt2,3 + j1)))*ZD(gt3,j2))))) + Conj(ZH(gt1,0))*(vd*(Sqr(g1) + 5*(Sqr(g2) - 4*QHd*Qq*Sqr(gp)))*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1)) + 2*(vd*(Sqr(g1) - 10*Qd*QHd*Sqr(gp))*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1)) - 5*(1.4142135623730951*SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,ZD(gt3,3 + j1)*TYd(j1,j2))) + 1.4142135623730951*SUM(j2,0,2,SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*Conj(TYd(j1,j2)))*ZD(gt3,j2)) + 2*vd*(SUM(j3,0,2,Conj(ZD(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,j1))*ZD(gt3,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gt3,j3)))))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Se, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yv = MODELPARAMETER(Yv);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.35355339059327373*(2*Conj(Lambdax)*Conj(ZH(gt1,2))*SUM(j2,0,2,SUM(j1,0,2,Conj(ZE(gt2,j1))*Yv(j1,j2))*ZV(gt4,3 + j2))*ZP(gt3,0) + 2*Conj(ZH(gt1,1))*SUM(j3,0,2,Conj(ZE(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Yv(j1,j2))*ZV(gt4,3 + j2)))*ZP(gt3,0) - Conj(ZH(gt1,1))*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZV(gt4,j1))*ZP(gt3,1) + 2*Conj(ZH(gt1,2))*Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1)))*ZV(gt4,j2))*ZP(gt3,1) + 2*Conj(ZH(gt1,1))*SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Conj(Yv(j3,j1))*Yv(j2,j1)))*ZV(gt4,j3))*ZP(gt3,1) + Conj(ZH(gt1,0))*(-(Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZV(gt4,j1))*ZP(gt3,0)) + 2*SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZV(gt4,j3))*ZP(gt3,0) + 2*SUM(j3,0,2,Conj(ZE(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Yv(j1,j2))*ZV(gt4,3 + j2)))*ZP(gt3,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Se, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto QHu = INPUTPARAMETER(QHu);
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qe = INPUTPARAMETER(Qe);
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto TYe = MODELPARAMETER(TYe);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto vu = MODELPARAMETER(vu);
   const auto vS = MODELPARAMETER(vS);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto vd = MODELPARAMETER(vd);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZE = MODELPARAMETER(ZE);

   const std::complex<double> result = 0.05*(10*Conj(ZH(gt1,2))*(-2*Ql*Qs*vS*Sqr(gp)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1)) - 2*Qe*Qs*vS*Sqr(gp)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1)) + vu*Conj(Lambdax)*SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1))) + vu*Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1)))*ZE(gt3,j2))) + Conj(ZH(gt1,1))*(vu*(3*Sqr(g1) - 5*(Sqr(g2) + 4*QHu*Ql*Sqr(gp)))*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1)) - 2*vu*(3*Sqr(g1) + 10*Qe*QHu*Sqr(gp))*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1)) + 10*vS*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt2,3 + j1)))*ZE(gt3,j2)))) - Conj(ZH(gt1,0))*(vd*(3*Sqr(g1) - 5*Sqr(g2) + 20*QHd*Ql*Sqr(gp))*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1)) + (-6*vd*Sqr(g1) + 20*Qe*QHd*vd*Sqr(gp))*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1)) + 10*(1.4142135623730951*SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,ZE(gt3,3 + j1)*TYe(j1,j2))) + 1.4142135623730951*SUM(j2,0,2,SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*Conj(TYe(j1,j2)))*ZE(gt3,j2)) + 2*vd*(SUM(j3,0,2,Conj(ZE(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Ye(j2,j1))*ZE(gt3,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZE(gt3,j3))))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Su, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qs = INPUTPARAMETER(Qs);
   const auto Qu = INPUTPARAMETER(Qu);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto TYu = MODELPARAMETER(TYu);
   const auto gp = MODELPARAMETER(gp);
   const auto vS = MODELPARAMETER(vS);
   const auto vd = MODELPARAMETER(vd);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vu = MODELPARAMETER(vu);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZU = MODELPARAMETER(ZU);

   const std::complex<double> result = 0.05*(10*Conj(ZH(gt1,2))*(-2*Qq*Qs*vS*Sqr(gp)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1)) - 2*Qs*Qu*vS*Sqr(gp)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1)) + vd*Conj(Lambdax)*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt3,3 + j1))) + vd*Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1)))*ZU(gt3,j2))) + Conj(ZH(gt1,0))*(vd*(Sqr(g1) - 5*(Sqr(g2) + 4*QHd*Qq*Sqr(gp)))*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1)) - 4*vd*(Sqr(g1) + 5*QHd*Qu*Sqr(gp))*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1)) + 10*vS*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt3,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1)))*ZU(gt3,j2)))) - Conj(ZH(gt1,1))*(vu*(Sqr(g1) - 5*Sqr(g2) + 20*QHu*Qq*Sqr(gp))*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1)) - 4*vu*(Sqr(g1) - 5*QHu*Qu*Sqr(gp))*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1)) + 10*(1.4142135623730951*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,ZU(gt3,3 + j1)*TYu(j1,j2))) + 1.4142135623730951*SUM(j2,0,2,SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j2)))*ZU(gt3,j2)) + 2*vu*(SUM(j3,0,2,Conj(ZU(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))*ZU(gt3,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gt3,j3))))));

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Sv, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qs = INPUTPARAMETER(Qs);
   const auto Qv = INPUTPARAMETER(Qv);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto TYv = MODELPARAMETER(TYv);
   const auto gp = MODELPARAMETER(gp);
   const auto vS = MODELPARAMETER(vS);
   const auto vd = MODELPARAMETER(vd);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vu = MODELPARAMETER(vu);
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZV = MODELPARAMETER(ZV);

   const std::complex<double> result = 0.25*(-2*Conj(ZH(gt1,2))*(2*Ql*Qs*vS*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZV(gt3,j1)) + 2*Qs*Qv*vS*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gt2,3 + j1))*ZV(gt3,3 + j1)) - vd*(Lambdax*SUM(j2,0,2,Conj(ZV(gt2,3 + j2))*SUM(j1,0,2,Conj(Yv(j1,j2))*ZV(gt3,j1))) + Conj(Lambdax)*SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(gt2,j1))*Yv(j1,j2))*ZV(gt3,3 + j2)))) - Conj(ZH(gt1,0))*(vd*(0.6*Sqr(g1) + Sqr(g2) + 4*QHd*Ql*Sqr(gp))*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZV(gt3,j1)) + 4*QHd*Qv*vd*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gt2,3 + j1))*ZV(gt3,3 + j1)) - 2*vS*(Lambdax*SUM(j2,0,2,Conj(ZV(gt2,3 + j2))*SUM(j1,0,2,Conj(Yv(j1,j2))*ZV(gt3,j1))) + Conj(Lambdax)*SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(gt2,j1))*Yv(j1,j2))*ZV(gt3,3 + j2)))) - Conj(ZH(gt1,1))*(-(vu*(0.6*Sqr(g1) + Sqr(g2) - 4*QHu*Ql*Sqr(gp))*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZV(gt3,j1))) + 2*(2*QHu*Qv*vu*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gt2,3 + j1))*ZV(gt3,3 + j1)) + 1.4142135623730951*SUM(j2,0,2,Conj(ZV(gt2,3 + j2))*SUM(j1,0,2,Conj(TYv(j1,j2))*ZV(gt3,j1))) + 1.4142135623730951*SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(gt2,j1))*TYv(j1,j2))*ZV(gt3,3 + j2)) + 2*vu*SUM(j3,0,2,Conj(ZV(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j3))*Yv(j1,j2))*ZV(gt3,3 + j2))) + 2*vu*SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Conj(Yv(j3,j1))*Yv(j2,j1)))*ZV(gt3,j3)))));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::VZ, fields::VZp>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto vS = MODELPARAMETER(vS);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vu = MODELPARAMETER(vu);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.1*(20*vS*Conj(ZH(gt1,2))*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gp)*Sqr(Qs) - vd*Conj(ZH(gt1,0))*(5*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) - 7.745966692414834*g1*gp*QHd*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + Cos(ThetaWp)*Sin(ThetaWp)*(-20*Sqr(gp)*Sqr(QHd) + 3*Sqr(g1)*Sqr(Sin(ThetaW))) + 7.745966692414834*g1*gp*QHd*Sin(ThetaW)*Sqr(Sin(ThetaWp)) + 2*g2*Cos(ThetaW)*(3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) - 5*gp*QHd*Sqr(Cos(ThetaWp)) + 5*gp*QHd*Sqr(Sin(ThetaWp)))) - vu*Conj(ZH(gt1,1))*(5*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 7.745966692414834*g1*gp*QHu*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + Cos(ThetaWp)*Sin(ThetaWp)*(-20*Sqr(gp)*Sqr(QHu) + 3*Sqr(g1)*Sqr(Sin(ThetaW))) - 7.745966692414834*g1*gp*QHu*Sin(ThetaW)*Sqr(Sin(ThetaWp)) + 2*g2*Cos(ThetaW)*(3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) + 5*gp*QHu*Sqr(Cos(ThetaWp)) - 5*gp*QHu*Sqr(Sin(ThetaWp)))));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto vS = MODELPARAMETER(vS);
   const auto vd = MODELPARAMETER(vd);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(4*vS*Conj(ZH(gt1,2))*Sqr(gp)*Sqr(Qs)*Sqr(Sin(ThetaWp)) + vd*Conj(ZH(gt1,0))*Sqr(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHd*Sin(ThetaWp)) + vu*Conj(ZH(gt1,1))*Sqr(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*QHu*Sin(ThetaWp)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::VZp, fields::VZp>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto vS = MODELPARAMETER(vS);
   const auto vd = MODELPARAMETER(vd);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(4*vS*Conj(ZH(gt1,2))*Sqr(gp)*Sqr(Qs)*Sqr(Cos(ThetaWp)) + vd*Conj(ZH(gt1,0))*Sqr(-2*gp*QHd*Cos(ThetaWp) + g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp)) + vu*Conj(ZH(gt1,1))*Sqr(2*gp*QHu*Cos(ThetaWp) + g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp)));

   return {result};
}

ChiralVertex VertexImpl<fields::hh, typename fields::bar<fields::Fe>::type, fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const int gt1 = indices[1];
   const int gt2 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZER = MODELPARAMETER(ZER);

   const std::complex<double> left = -0.7071067811865475*Conj(ZH(gt3,0))*SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,j2)));

   const std::complex<double> right = -0.7071067811865475*Conj(ZH(gt3,0))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZEL(gt1,j2));

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

   const std::complex<double> result = -0.3872983346207417*g1*g2*Cos(ThetaW)*(Conj(ZH(gt1,0))*ZP(gt2,0) - Conj(ZH(gt1,1))*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, typename fields::conj<fields::Hpm>::type, fields::VWm, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*g2*(Conj(ZH(gt1,0))*(2*gp*QHd*Cos(ThetaWp) - 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt2,0) + Conj(ZH(gt1,1))*(2*gp*QHu*Cos(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, typename fields::conj<fields::Hpm>::type, fields::VWm, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*g2*(Conj(ZH(gt1,0))*(0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHd*Sin(ThetaWp))*ZP(gt2,0) + Conj(ZH(gt1,1))*(-0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHu*Sin(ThetaWp))*ZP(gt2,1));

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

   const std::complex<double> result = -0.5*g2*(Conj(ZH(gt1,0))*ZP(gt2,0) - Conj(ZH(gt1,1))*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::hh, typename fields::conj<fields::VWm>::type, fields::VWm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> result = 0.5*(vd*Conj(ZH(gt1,0)) + vu*Conj(ZH(gt1,1)))*Sqr(g2);

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::Hpm, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto QHu = INPUTPARAMETER(QHu);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-2*(0.6*Sqr(g1) + Sqr(g2) + 4*Sqr(gp)*Sqr(QHd))*ZP(gt1,0)*ZP(gt2,0)*ZP(gt3,0)*ZP(gt4,0) - 2*(0.6*Sqr(g1) + Sqr(g2) + 4*Sqr(gp)*Sqr(QHu))*ZP(gt1,1)*ZP(gt2,1)*ZP(gt3,1)*ZP(gt4,1) + (-4*AbsSqr(Lambdax) + 0.6*Sqr(g1) + Sqr(g2) - 4*QHd*QHu*Sqr(gp))*ZP(gt1,1)*ZP(gt2,0)*(ZP(gt3,1)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,1)) + (-4*AbsSqr(Lambdax) + 0.6*Sqr(g1) + Sqr(g2) - 4*QHd*QHu*Sqr(gp))*ZP(gt1,0)*ZP(gt2,1)*(ZP(gt3,1)*ZP(gt4,0) + ZP(gt3,0)*ZP(gt4,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::Sd, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto Qq = INPUTPARAMETER(Qq);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto Qd = INPUTPARAMETER(Qd);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.05*(2*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt4,3 + j1))*((Sqr(g1) - 10*Qd*QHd*Sqr(gp))*ZP(gt1,0)*ZP(gt3,0) - (Sqr(g1) + 10*Qd*QHu*Sqr(gp))*ZP(gt1,1)*ZP(gt3,1)) + SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt4,j1))*((Sqr(g1) - 5*(Sqr(g2) + 4*QHd*Qq*Sqr(gp)))*ZP(gt1,0)*ZP(gt3,0) - (Sqr(g1) - 5*Sqr(g2) + 20*QHu*Qq*Sqr(gp))*ZP(gt1,1)*ZP(gt3,1)) - 20*(SUM(j3,0,2,Conj(ZD(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,j1))*ZD(gt4,3 + j2)))*ZP(gt1,0)*ZP(gt3,0) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZD(gt4,j3))*ZP(gt1,1)*ZP(gt3,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::Se, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Qe = INPUTPARAMETER(Qe);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto Ql = INPUTPARAMETER(Ql);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.05*(SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt4,3 + j1))*(2*(3*Sqr(g1) - 10*Qe*QHd*Sqr(gp))*ZP(gt1,0)*ZP(gt3,0) - 2*(3*Sqr(g1) + 10*Qe*QHu*Sqr(gp))*ZP(gt1,1)*ZP(gt3,1)) - SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt4,j1))*((3*Sqr(g1) + 5*Sqr(g2) + 20*QHd*Ql*Sqr(gp))*ZP(gt1,0)*ZP(gt3,0) + (-3*Sqr(g1) - 5*Sqr(g2) + 20*QHu*Ql*Sqr(gp))*ZP(gt1,1)*ZP(gt3,1)) - 20*(SUM(j3,0,2,Conj(ZE(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Ye(j2,j1))*ZE(gt4,3 + j2)))*ZP(gt1,0)*ZP(gt3,0) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Conj(Yv(j3,j1))*Yv(j2,j1)))*ZE(gt4,j3))*ZP(gt1,1)*ZP(gt3,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::Su, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto Qq = INPUTPARAMETER(Qq);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto Qu = INPUTPARAMETER(Qu);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.05*(SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt4,j1))*((Sqr(g1) + 5*Sqr(g2) - 20*QHd*Qq*Sqr(gp))*ZP(gt1,0)*ZP(gt3,0) - (Sqr(g1) + 5*Sqr(g2) + 20*QHu*Qq*Sqr(gp))*ZP(gt1,1)*ZP(gt3,1)) - 4*(SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt4,3 + j1))*((Sqr(g1) + 5*QHd*Qu*Sqr(gp))*ZP(gt1,0)*ZP(gt3,0) - (Sqr(g1) - 5*QHu*Qu*Sqr(gp))*ZP(gt1,1)*ZP(gt3,1)) + 5*(SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZU(gt4,j3))*ZP(gt1,0)*ZP(gt3,0) + SUM(j3,0,2,Conj(ZU(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))*ZU(gt4,3 + j2)))*ZP(gt1,1)*ZP(gt3,1))));

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
   const auto vS = MODELPARAMETER(vS);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-1.4142135623730951*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZD(gt3,j1))*(vd*ZP(gt1,0) + vu*ZP(gt1,1)) + 2*(2*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,ZD(gt3,3 + j1)*TYd(j1,j2)))*ZP(gt1,0) + 1.4142135623730951*vS*Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt2,3 + j1)))*ZD(gt3,j2))*ZP(gt1,0) + 1.4142135623730951*vu*SUM(j3,0,2,Conj(ZU(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yd(j2,j1))*ZD(gt3,3 + j2)))*ZP(gt1,0) + 1.4142135623730951*vd*SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gt3,j3))*ZP(gt1,0) + 1.4142135623730951*vS*Conj(Lambdax)*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1)))*ZP(gt1,1) + 2*SUM(j2,0,2,SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*Conj(TYu(j1,j2)))*ZD(gt3,j2))*ZP(gt1,1) + 1.4142135623730951*vd*SUM(j3,0,2,Conj(ZU(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yd(j2,j1))*ZD(gt3,3 + j2)))*ZP(gt1,1) + 1.4142135623730951*vu*SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZD(gt3,j3))*ZP(gt1,1)));

   return {result};
}

ScalarVertex VertexImpl<fields::Hpm, fields::Sv, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto Ql = INPUTPARAMETER(Ql);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto Qv = INPUTPARAMETER(Qv);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-(SUM(j1,0,2,Conj(ZV(gt2,j1))*ZV(gt4,j1))*((0.6*Sqr(g1) - Sqr(g2) + 4*QHd*Ql*Sqr(gp))*ZP(gt1,0)*ZP(gt3,0) + (-0.6*Sqr(g1) + Sqr(g2) + 4*QHu*Ql*Sqr(gp))*ZP(gt1,1)*ZP(gt3,1))) - 4*(SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZV(gt4,j3))*ZP(gt1,0)*ZP(gt3,0) + SUM(j3,0,2,Conj(ZV(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j3))*Yv(j1,j2))*ZV(gt4,3 + j2)))*ZP(gt1,1)*ZP(gt3,1) + Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gt2,3 + j1))*ZV(gt4,3 + j1))*(QHd*ZP(gt1,0)*ZP(gt3,0) + QHu*ZP(gt1,1)*ZP(gt3,1))));

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
   const auto vS = MODELPARAMETER(vS);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Yv = MODELPARAMETER(Yv);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-1.4142135623730951*Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZE(gt3,j1))*(vd*ZP(gt1,0) + vu*ZP(gt1,1)) + 2*(1.4142135623730951*vS*Lambdax*SUM(j2,0,2,Conj(ZV(gt2,3 + j2))*SUM(j1,0,2,Conj(Yv(j1,j2))*ZE(gt3,j1)))*ZP(gt1,0) + 2*SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,ZE(gt3,3 + j1)*TYe(j1,j2)))*ZP(gt1,0) + 1.4142135623730951*vu*SUM(j3,0,2,Conj(ZV(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j3))*Ye(j2,j1))*ZE(gt3,3 + j2)))*ZP(gt1,0) + 1.4142135623730951*vd*SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZE(gt3,j3))*ZP(gt1,0) + 2*SUM(j2,0,2,Conj(ZV(gt2,3 + j2))*SUM(j1,0,2,Conj(TYv(j1,j2))*ZE(gt3,j1)))*ZP(gt1,1) + 1.4142135623730951*vS*Conj(Lambdax)*SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)))*ZP(gt1,1) + 1.4142135623730951*vd*SUM(j3,0,2,Conj(ZV(gt2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j3))*Ye(j2,j1))*ZE(gt3,3 + j2)))*ZP(gt1,1) + 1.4142135623730951*vu*SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Conj(Yv(j3,j1))*Yv(j2,j1)))*ZE(gt3,j3))*ZP(gt1,1)));

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

InverseMetricVertex VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, fields::VP, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.5*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*((-2*gp*QHd*Cos(ThetaWp) - g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt1,0)*ZP(gt2,0) + (2*gp*QHu*Cos(ThetaWp) - g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt1,1)*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, fields::VP, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -0.5*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*((-(g2*Cos(ThetaW)*Cos(ThetaWp)) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHd*Sin(ThetaWp))*ZP(gt1,0)*ZP(gt2,0) + (-(g2*Cos(ThetaW)*Cos(ThetaWp)) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*QHu*Sin(ThetaWp))*ZP(gt1,1)*ZP(gt2,1));

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

InverseMetricVertex VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, fields::VZ, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.1*((-5*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 7.745966692414834*g1*gp*QHd*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + Cos(ThetaWp)*Sin(ThetaWp)*(20*Sqr(gp)*Sqr(QHd) - 3*Sqr(g1)*Sqr(Sin(ThetaW))) - 7.745966692414834*g1*gp*QHd*Sin(ThetaW)*Sqr(Sin(ThetaWp)) + 2*g2*Cos(ThetaW)*(3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) - 5*gp*QHd*Sqr(Cos(ThetaWp)) + 5*gp*QHd*Sqr(Sin(ThetaWp))))*ZP(gt1,0)*ZP(gt2,0) + (-5*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) - 7.745966692414834*g1*gp*QHu*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + Cos(ThetaWp)*Sin(ThetaWp)*(20*Sqr(gp)*Sqr(QHu) - 3*Sqr(g1)*Sqr(Sin(ThetaW))) + 7.745966692414834*g1*gp*QHu*Sin(ThetaW)*Sqr(Sin(ThetaWp)) + 2*g2*Cos(ThetaW)*(3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) + 5*gp*QHu*Sqr(Cos(ThetaWp)) - 5*gp*QHu*Sqr(Sin(ThetaWp))))*ZP(gt1,1)*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.5*(Sqr(-(g2*Cos(ThetaW)*Cos(ThetaWp)) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHd*Sin(ThetaWp))*ZP(gt1,0)*ZP(gt2,0) + Sqr(g2*Cos(ThetaW)*Cos(ThetaWp) - 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHu*Sin(ThetaWp))*ZP(gt1,1)*ZP(gt2,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, fields::VZp, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(Sqr(2*gp*QHd*Cos(ThetaWp) + g2*Cos(ThetaW)*Sin(ThetaWp) - 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt1,0)*ZP(gt2,0) + Sqr(2*gp*QHu*Cos(ThetaWp) - g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt1,1)*ZP(gt2,1));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.1*((-10*gp*QHd*Cos(ThetaWp) + (-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*ZP(gt1,0)*ZP(gt2,0) + (10*gp*QHu*Cos(ThetaWp) + (-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*ZP(gt1,1)*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.5*((g2*Cos(ThetaW)*Cos(ThetaWp) - 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*QHd*Sin(ThetaWp))*ZP(gt1,0)*ZP(gt2,0) + (g2*Cos(ThetaW)*Cos(ThetaWp) - 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHu*Sin(ThetaWp))*ZP(gt1,1)*ZP(gt2,1));

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

InverseMetricVertex VertexImpl<fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VZp>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*g2*(vd*(2*gp*QHd*Cos(ThetaWp) - 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt1,0) + vu*(2*gp*QHu*Cos(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt1,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VZ>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*g2*(vd*(0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHd*Sin(ThetaWp))*ZP(gt1,0) + vu*(-0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHu*Sin(ThetaWp))*ZP(gt1,1));

   return {result};
}

ScalarVertex VertexImpl<fields::Sd, fields::Sd, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qd = INPUTPARAMETER(Qd);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto g3 = MODELPARAMETER(g3);
   const auto gp = MODELPARAMETER(gp);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = 0.008333333333333333*(-(Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt2,j2))*ZD(gt3,j2))) - 15*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt2,j2))*ZD(gt3,j2)) + 10*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt2,j2))*ZD(gt3,j2)) - 60*Sqr(gp)*Sqr(Qq)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt2,j2))*ZD(gt3,j2)) - 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt4,3 + j1))*SUM(j2,0,2,Conj(ZD(gt2,j2))*ZD(gt3,j2)) - 10*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt4,3 + j1))*SUM(j2,0,2,Conj(ZD(gt2,j2))*ZD(gt3,j2)) - 60*Qd*Qq*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt4,3 + j1))*SUM(j2,0,2,Conj(ZD(gt2,j2))*ZD(gt3,j2)) - 30*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt4,j1))*(SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) - SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2))) + 30*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt4,3 + j1))*(SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) - SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2))) - 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt2,3 + j2))*ZD(gt3,3 + j2)) - 10*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt2,3 + j2))*ZD(gt3,3 + j2)) - 60*Qd*Qq*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt4,j1))*SUM(j2,0,2,Conj(ZD(gt2,3 + j2))*ZD(gt3,3 + j2)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt4,3 + j1))*SUM(j2,0,2,Conj(ZD(gt2,3 + j2))*ZD(gt3,3 + j2)) + 10*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt4,3 + j1))*SUM(j2,0,2,Conj(ZD(gt2,3 + j2))*ZD(gt3,3 + j2)) - 60*Sqr(gp)*Sqr(Qd)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt4,3 + j1))*SUM(j2,0,2,Conj(ZD(gt2,3 + j2))*ZD(gt3,3 + j2)) - Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt4,j2)) - 15*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt4,j2)) + 10*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt4,j2)) - 60*Sqr(gp)*Sqr(Qq)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt4,j2)) - 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt4,j2)) - 10*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt4,j2)) - 60*Qd*Qq*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt4,j2)) - 30*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt2,j2))*ZD(gt4,j2)) + 30*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZD(gt2,j2))*ZD(gt4,j2)) - 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt4,3 + j2)) - 10*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt4,3 + j2)) - 60*Qd*Qq*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gt2,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt4,3 + j2)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt4,3 + j2)) + 10*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt4,3 + j2)) - 60*Sqr(gp)*Sqr(Qd)*SUM(j1,0,2,Conj(ZD(gt2,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt4,3 + j2)) + 30*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt2,3 + j2))*ZD(gt4,3 + j2)) - 30*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZD(gt2,3 + j2))*ZD(gt4,3 + j2)) - 120*SUM(j2,0,2,Conj(ZD(gt1,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt4,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yd(j3,j4))*Conj(ZD(gt2,3 + j3)))*ZD(gt3,j4)) - 120*SUM(j2,0,2,Conj(ZD(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yd(j3,j4))*Conj(ZD(gt1,3 + j3)))*ZD(gt4,j4)));

   return {result};
}

ScalarVertex VertexImpl<fields::Sd, fields::Se, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Qe = INPUTPARAMETER(Qe);
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qd = INPUTPARAMETER(Qd);
   const auto Ql = INPUTPARAMETER(Ql);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = 0.025*(-2*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt4,3 + j1))*((Sqr(g1) + 10*Qe*Qq*Sqr(gp))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) + 2*(Sqr(g1) + 5*Qd*Qe*Sqr(gp))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2))) + SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt4,j1))*((Sqr(g1) - 5*(Sqr(g2) + 4*Ql*Qq*Sqr(gp)))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) + 2*(Sqr(g1) - 10*Qd*Ql*Sqr(gp))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2))) + Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) - 5*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) - 20*Ql*Qq*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) + 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) - 20*Qd*Ql*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) - 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt4,3 + j2)) - 20*Qe*Qq*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt4,3 + j2)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt4,3 + j2)) - 20*Qd*Qe*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt4,3 + j2)) - 40*SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt4,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yd(j3,j4))*Conj(ZD(gt1,3 + j3)))*ZD(gt3,j4)) - 40*SUM(j2,0,2,Conj(ZD(gt1,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Ye(j3,j4))*Conj(ZE(gt2,3 + j3)))*ZE(gt4,j4)));

   return {result};
}

ScalarVertex VertexImpl<fields::Sd, fields::Su, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qd = INPUTPARAMETER(Qd);
   const auto Qu = INPUTPARAMETER(Qu);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto g3 = MODELPARAMETER(g3);
   const auto gp = MODELPARAMETER(gp);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = 0.008333333333333333*(-(SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt4,j1))*((Sqr(g1) - 15*Sqr(g2) - 10*Sqr(g3) + 60*Sqr(gp)*Sqr(Qq))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) + 2*(Sqr(g1) + 5*Sqr(g3) + 30*Qd*Qq*Sqr(gp))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2)))) + 2*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt4,3 + j1))*((2*Sqr(g1) - 5*(Sqr(g3) + 6*Qq*Qu*Sqr(gp)))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) + (4*Sqr(g1) + 5*(Sqr(g3) - 6*Qd*Qu*Sqr(gp)))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2))) - Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) + 15*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) + 10*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) - 60*Sqr(gp)*Sqr(Qq)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) - 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) - 10*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) - 60*Qd*Qq*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) + 4*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)) - 10*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)) - 60*Qq*Qu*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)) + 8*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)) + 10*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)) - 60*Qd*Qu*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)));

   return {result};
}

ScalarVertex VertexImpl<fields::Sd, fields::Sv, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Qv = INPUTPARAMETER(Qv);
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qd = INPUTPARAMETER(Qd);
   const auto Ql = INPUTPARAMETER(Ql);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZD = MODELPARAMETER(ZD);

   const std::complex<double> result = 0.025*(-20*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gt2,3 + j1))*ZV(gt4,3 + j1))*(Qq*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) + Qd*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2))) + SUM(j1,0,2,Conj(ZV(gt2,j1))*ZV(gt4,j1))*((Sqr(g1) + 5*Sqr(g2) - 20*Ql*Qq*Sqr(gp))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZD(gt3,j2)) + 2*(Sqr(g1) - 10*Qd*Ql*Sqr(gp))*SUM(j2,0,2,Conj(ZD(gt1,3 + j2))*ZD(gt3,3 + j2))) + Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt4,j2)) + 5*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt4,j2)) - 20*Ql*Qq*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt4,j2)) + 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt4,j2)) - 20*Qd*Ql*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt4,j2)) - 20*Qq*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZV(gt2,3 + j2))*ZV(gt4,3 + j2)) - 20*Qd*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt3,3 + j1))*SUM(j2,0,2,Conj(ZV(gt2,3 + j2))*ZV(gt4,3 + j2)));

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

   const std::complex<double> result = 0.25*(-(Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZU(gt4,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZE(gt3,j2))) - Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZD(gt1,j2))*ZU(gt4,j2)) - 4*(SUM(j2,0,2,Conj(ZD(gt1,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt4,3 + j1)))*SUM(j4,0,2,Conj(ZV(gt2,3 + j4))*SUM(j3,0,2,Conj(Yv(j3,j4))*ZE(gt3,j3))) + SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yd(j3,j4))*Conj(ZD(gt1,3 + j3)))*ZU(gt4,j4))));

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
   const auto vS = MODELPARAMETER(vS);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-1.4142135623730951*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZU(gt3,j1))*(vd*ZP(gt2,0) + vu*ZP(gt2,1)) + 2*(1.4142135623730951*vS*Conj(Lambdax)*SUM(j2,0,2,Conj(ZD(gt1,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt3,3 + j1)))*ZP(gt2,0) + 2*SUM(j2,0,2,SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*Conj(TYd(j1,j2)))*ZU(gt3,j2))*ZP(gt2,0) + 1.4142135623730951*vu*SUM(j3,0,2,Conj(ZD(gt1,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yu(j2,j1))*ZU(gt3,3 + j2)))*ZP(gt2,0) + 1.4142135623730951*vd*SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt1,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZU(gt3,j3))*ZP(gt2,0) + 2*SUM(j2,0,2,Conj(ZD(gt1,j2))*SUM(j1,0,2,ZU(gt3,3 + j1)*TYu(j1,j2)))*ZP(gt2,1) + 1.4142135623730951*vS*Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt1,3 + j1)))*ZU(gt3,j2))*ZP(gt2,1) + 1.4142135623730951*vd*SUM(j3,0,2,Conj(ZD(gt1,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yu(j2,j1))*ZU(gt3,3 + j2)))*ZP(gt2,1) + 1.4142135623730951*vu*SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gt1,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gt3,j3))*ZP(gt2,1)));

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

InverseMetricVertex VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VP, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qd = INPUTPARAMETER(Qd);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05555555555555555*((0.7745966692414834*g1*Cos(ThetaW) - 3*g2*Sin(ThetaW))*(6*gp*Qq*Cos(ThetaWp) + 3*g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1)) + 3.0983866769659336*g1*Cos(ThetaW)*(3*gp*Qd*Cos(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VP, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qd = INPUTPARAMETER(Qd);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05555555555555555*(-((0.7745966692414834*g1*Cos(ThetaW) - 3*g2*Sin(ThetaW))*(3*g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 6*gp*Qq*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1))) - 3.0983866769659336*g1*Cos(ThetaW)*(0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 3*gp*Qd*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1)));

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

InverseMetricVertex VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VZ, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qd = INPUTPARAMETER(Qd);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.03333333333333333*(-((15*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 7.745966692414834*g1*gp*Qq*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + Cos(ThetaWp)*Sin(ThetaWp)*(-60*Sqr(gp)*Sqr(Qq) + Sqr(g1)*Sqr(Sin(ThetaW))) - 7.745966692414834*g1*gp*Qq*Sin(ThetaW)*Sqr(Sin(ThetaWp)) + 2*g2*Cos(ThetaW)*(3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) + 15*gp*Qq*Sqr(Cos(ThetaWp)) - 15*gp*Qq*Sqr(Sin(ThetaWp))))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1))) - 2*(-15*Sin(2*ThetaWp)*Sqr(gp)*Sqr(Qd) + 7.745966692414834*g1*gp*Qd*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + Sin(2*ThetaWp)*Sqr(g1)*Sqr(Sin(ThetaW)) - 7.745966692414834*g1*gp*Qd*Sin(ThetaW)*Sqr(Sin(ThetaWp)))*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qd = INPUTPARAMETER(Qd);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05555555555555555*(Sqr(3*g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 6*gp*Qq*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1)) + 4*Sqr(0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 3*gp*Qd*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VZp, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qd = INPUTPARAMETER(Qd);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05555555555555555*(Sqr(6*gp*Qq*Cos(ThetaWp) + 3*g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1)) + 4*Sqr(3*gp*Qd*Cos(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1)));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qd = INPUTPARAMETER(Qd);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.16666666666666666*(-((6*gp*Qq*Cos(ThetaWp) + 3*g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1))) + 2*(3*gp*Qd*Cos(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qd = INPUTPARAMETER(Qd);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.16666666666666666*((3*g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 6*gp*Qq*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1)) - 2*(0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 3*gp*Qd*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1)));

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

InverseMetricVertex VertexImpl<fields::Sd, typename fields::conj<fields::Su>::type, typename fields::conj<fields::VWm>::type, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.2357022603955158*g2*(6*gp*Qq*Cos(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZU(gt2,j1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sd, typename fields::conj<fields::Su>::type, typename fields::conj<fields::VWm>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.2357022603955158*g2*(0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 6*gp*Qq*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZU(gt2,j1));

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
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qe = INPUTPARAMETER(Qe);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);

   const std::complex<double> result = 0.025*(-3*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt4,j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt3,j2)) - 5*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt4,j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt3,j2)) - 20*Sqr(gp)*Sqr(Ql)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt4,j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt3,j2)) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt4,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt3,j2)) - 20*Qe*Ql*Sqr(gp)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt4,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt3,j2)) - SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt4,j1))*((3*Sqr(g1) + 5*Sqr(g2) + 20*Sqr(gp)*Sqr(Ql))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt3,j2)) + 2*(-3*Sqr(g1) + 10*Qe*Ql*Sqr(gp))*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt3,3 + j2))) + 2*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt4,3 + j1))*((3*Sqr(g1) - 10*Qe*Ql*Sqr(gp))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt3,j2)) - 2*(3*Sqr(g1) + 5*Sqr(gp)*Sqr(Qe))*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt3,3 + j2))) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt4,j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt3,3 + j2)) - 20*Qe*Ql*Sqr(gp)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt4,j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt3,3 + j2)) - 12*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt4,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt3,3 + j2)) - 20*Sqr(gp)*Sqr(Qe)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt4,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt3,3 + j2)) - 3*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt4,j2)) - 5*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt4,j2)) - 20*Sqr(gp)*Sqr(Ql)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt4,j2)) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt4,j2)) - 20*Qe*Ql*Sqr(gp)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt4,j2)) - 3*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) - 5*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) - 20*Sqr(gp)*Sqr(Ql)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) - 20*Qe*Ql*Sqr(gp)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,j2))*ZE(gt4,j2)) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt4,3 + j2)) - 20*Qe*Ql*Sqr(gp)*SUM(j1,0,2,Conj(ZE(gt2,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt4,3 + j2)) - 12*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt4,3 + j2)) - 20*Sqr(gp)*Sqr(Qe)*SUM(j1,0,2,Conj(ZE(gt2,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt4,3 + j2)) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt4,3 + j2)) - 20*Qe*Ql*Sqr(gp)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt4,3 + j2)) - 12*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt4,3 + j2)) - 20*Sqr(gp)*Sqr(Qe)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZE(gt2,3 + j2))*ZE(gt4,3 + j2)) - 40*SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt4,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Ye(j3,j4))*Conj(ZE(gt1,3 + j3)))*ZE(gt3,j4)) - 40*SUM(j2,0,2,Conj(ZE(gt1,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt4,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Ye(j3,j4))*Conj(ZE(gt2,3 + j3)))*ZE(gt3,j4)) - 40*SUM(j2,0,2,Conj(ZE(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Ye(j3,j4))*Conj(ZE(gt1,3 + j3)))*ZE(gt4,j4)) - 40*SUM(j2,0,2,Conj(ZE(gt1,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Ye(j3,j4))*Conj(ZE(gt2,3 + j3)))*ZE(gt4,j4)));

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

   const std::complex<double> result = 0.25*(-(Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZV(gt4,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZD(gt3,j2))) - Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZD(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZV(gt4,j2)) - 4*(SUM(j2,0,2,SUM(j1,0,2,Conj(ZE(gt1,j1))*Yv(j1,j2))*ZV(gt4,3 + j2))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yu(j3,j4))*Conj(ZU(gt2,3 + j3)))*ZD(gt3,j4)) + SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gt3,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Ye(j3,j4))*Conj(ZE(gt1,3 + j3)))*ZV(gt4,j4))));

   return {result};
}

ScalarVertex VertexImpl<fields::Se, fields::Su, typename fields::conj<fields::Se>::type, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qe = INPUTPARAMETER(Qe);
   const auto Qu = INPUTPARAMETER(Qu);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZE = MODELPARAMETER(ZE);

   const std::complex<double> result = 0.025*(SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt4,j1))*((Sqr(g1) + 5*Sqr(g2) - 20*Ql*Qq*Sqr(gp))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt3,j2)) - 2*(Sqr(g1) + 10*Qe*Qq*Sqr(gp))*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt3,3 + j2))) - 4*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt4,3 + j1))*((Sqr(g1) + 5*Ql*Qu*Sqr(gp))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt3,j2)) + (-2*Sqr(g1) + 5*Qe*Qu*Sqr(gp))*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt3,3 + j2))) + Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) + 5*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) - 20*Ql*Qq*Sqr(gp)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) - 2*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) - 20*Qe*Qq*Sqr(gp)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)) - 20*Ql*Qu*Sqr(gp)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)) + 8*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)) - 20*Qe*Qu*Sqr(gp)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)));

   return {result};
}

ScalarVertex VertexImpl<fields::Se, fields::Sv, typename fields::conj<fields::Se>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Qv = INPUTPARAMETER(Qv);
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qe = INPUTPARAMETER(Qe);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto Yv = MODELPARAMETER(Yv);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);

   const std::complex<double> result = 0.125*(-2*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZV(gt4,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZE(gt3,j2)) - 4*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gt2,3 + j1))*ZV(gt4,3 + j1))*(Ql*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt3,j2)) + Qe*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt3,3 + j2))) - 0.2*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZV(gt4,j1))*((3*Sqr(g1) - 5*Sqr(g2) + 20*Sqr(gp)*Sqr(Ql))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZE(gt3,j2)) + 2*(-3*Sqr(g1) + 10*Qe*Ql*Sqr(gp))*SUM(j2,0,2,Conj(ZE(gt1,3 + j2))*ZE(gt3,3 + j2))) - 2*Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZE(gt1,j2))*ZV(gt4,j2)) - 0.6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt4,j2)) + Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt4,j2)) - 4*Sqr(gp)*Sqr(Ql)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt4,j2)) + 1.2*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt4,j2)) - 4*Qe*Ql*Sqr(gp)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt4,j2)) - 4*Ql*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt3,j1))*SUM(j2,0,2,Conj(ZV(gt2,3 + j2))*ZV(gt4,3 + j2)) - 4*Qe*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt3,3 + j1))*SUM(j2,0,2,Conj(ZV(gt2,3 + j2))*ZV(gt4,3 + j2)) - 8*SUM(j2,0,2,SUM(j1,0,2,Conj(ZE(gt1,j1))*Yv(j1,j2))*ZV(gt4,3 + j2))*SUM(j4,0,2,Conj(ZV(gt2,3 + j4))*SUM(j3,0,2,Conj(Yv(j3,j4))*ZE(gt3,j3))) - 8*SUM(j2,0,2,Conj(ZV(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Ye(j3,j4))*Conj(ZE(gt1,3 + j3)))*ZV(gt4,j4)));

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
   const auto vS = MODELPARAMETER(vS);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Yv = MODELPARAMETER(Yv);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> result = 0.25*(-1.4142135623730951*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZV(gt3,j1))*(vd*ZP(gt2,0) + vu*ZP(gt2,1)) + 2*(2*SUM(j2,0,2,SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*Conj(TYe(j1,j2)))*ZV(gt3,j2))*ZP(gt2,0) + 1.4142135623730951*vS*Conj(Lambdax)*SUM(j2,0,2,SUM(j1,0,2,Conj(ZE(gt1,j1))*Yv(j1,j2))*ZV(gt3,3 + j2))*ZP(gt2,0) + 1.4142135623730951*vu*SUM(j3,0,2,Conj(ZE(gt1,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Yv(j1,j2))*ZV(gt3,3 + j2)))*ZP(gt2,0) + 1.4142135623730951*vd*SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gt1,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZV(gt3,j3))*ZP(gt2,0) + 1.4142135623730951*vS*Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt1,3 + j1)))*ZV(gt3,j2))*ZP(gt2,1) + 2*SUM(j2,0,2,SUM(j1,0,2,Conj(ZE(gt1,j1))*TYv(j1,j2))*ZV(gt3,3 + j2))*ZP(gt2,1) + 1.4142135623730951*vd*SUM(j3,0,2,Conj(ZE(gt1,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Yv(j1,j2))*ZV(gt3,3 + j2)))*ZP(gt2,1) + 1.4142135623730951*vu*SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gt1,j2))*SUM(j1,0,2,Conj(Yv(j3,j1))*Yv(j2,j1)))*ZV(gt3,j3))*ZP(gt2,1)));

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

InverseMetricVertex VertexImpl<fields::Se, typename fields::conj<fields::Se>::type, fields::VP, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qe = INPUTPARAMETER(Qe);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.5*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*(-2*gp*Ql*Cos(ThetaWp) - g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) + 1.5491933384829668*g1*Cos(ThetaW)*(gp*Qe*Cos(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Se, typename fields::conj<fields::Se>::type, fields::VP, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qe = INPUTPARAMETER(Qe);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -0.5*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*(-(g2*Cos(ThetaW)*Cos(ThetaWp)) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*Ql*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) - 1.5491933384829668*g1*Cos(ThetaW)*(0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - gp*Qe*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1));

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

InverseMetricVertex VertexImpl<fields::Se, typename fields::conj<fields::Se>::type, fields::VZ, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qe = INPUTPARAMETER(Qe);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.1*((-5*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 7.745966692414834*g1*gp*Ql*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + Cos(ThetaWp)*Sin(ThetaWp)*(20*Sqr(gp)*Sqr(Ql) - 3*Sqr(g1)*Sqr(Sin(ThetaW))) - 7.745966692414834*g1*gp*Ql*Sin(ThetaW)*Sqr(Sin(ThetaWp)) + 2*g2*Cos(ThetaW)*(3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) - 5*gp*Ql*Sqr(Cos(ThetaWp)) + 5*gp*Ql*Sqr(Sin(ThetaWp))))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) + 2*(5*Sin(2*ThetaWp)*Sqr(gp)*Sqr(Qe) - 7.745966692414834*g1*gp*Qe*Sin(ThetaW)*Sqr(Cos(ThetaWp)) - 3*Sin(2*ThetaWp)*Sqr(g1)*Sqr(Sin(ThetaW)) + 7.745966692414834*g1*gp*Qe*Sin(ThetaW)*Sqr(Sin(ThetaWp)))*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Se, typename fields::conj<fields::Se>::type, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qe = INPUTPARAMETER(Qe);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.5*(Sqr(-(g2*Cos(ThetaW)*Cos(ThetaWp)) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*Ql*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) + 4*Sqr(0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - gp*Qe*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Se, typename fields::conj<fields::Se>::type, fields::VZp, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qe = INPUTPARAMETER(Qe);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(Sqr(2*gp*Ql*Cos(ThetaWp) + g2*Cos(ThetaW)*Sin(ThetaWp) - 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) + 4*Sqr(gp*Qe*Cos(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1)));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::Se, typename fields::conj<fields::Se>::type, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qe = INPUTPARAMETER(Qe);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.5*(2*gp*Ql*Cos(ThetaWp) + g2*Cos(ThetaW)*Sin(ThetaWp) - 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) + (gp*Qe*Cos(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::Se, typename fields::conj<fields::Se>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qe = INPUTPARAMETER(Qe);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.5*(g2*Cos(ThetaW)*Cos(ThetaWp) - 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*Ql*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) + (-0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + gp*Qe*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1));

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

InverseMetricVertex VertexImpl<fields::Se, typename fields::conj<fields::Sv>::type, typename fields::conj<fields::VWm>::type, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.7071067811865475*g2*(2*gp*Ql*Cos(ThetaWp) - 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZV(gt2,j1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Se, typename fields::conj<fields::Sv>::type, typename fields::conj<fields::VWm>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.7071067811865475*g2*(0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*Ql*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZV(gt2,j1));

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
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qu = INPUTPARAMETER(Qu);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto g3 = MODELPARAMETER(g3);
   const auto gp = MODELPARAMETER(gp);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZU = MODELPARAMETER(ZU);

   const std::complex<double> result = 0.008333333333333333*(-(Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt4,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt3,j2))) - 15*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt4,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt3,j2)) + 10*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt4,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt3,j2)) - 60*Sqr(gp)*Sqr(Qq)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt4,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt3,j2)) + 4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt4,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt3,j2)) - 10*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt4,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt3,j2)) - 60*Qq*Qu*Sqr(gp)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt4,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt3,j2)) - 30*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt4,j1))*(SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt3,j2)) - SUM(j2,0,2,Conj(ZU(gt1,3 + j2))*ZU(gt3,3 + j2))) + 30*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt4,3 + j1))*(SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt3,j2)) - SUM(j2,0,2,Conj(ZU(gt1,3 + j2))*ZU(gt3,3 + j2))) + 4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt4,j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt3,3 + j2)) - 10*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt4,j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt3,3 + j2)) - 60*Qq*Qu*Sqr(gp)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt4,j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt3,3 + j2)) - 16*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt4,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt3,3 + j2)) + 10*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt4,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt3,3 + j2)) - 60*Sqr(gp)*Sqr(Qu)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt4,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt3,3 + j2)) - Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt4,j2)) - 15*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt4,j2)) + 10*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt4,j2)) - 60*Sqr(gp)*Sqr(Qq)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt4,j2)) + 4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt4,j2)) - 10*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt4,j2)) - 60*Qq*Qu*Sqr(gp)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt4,j2)) - 30*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) + 30*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,j2))*ZU(gt4,j2)) + 4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt1,3 + j2))*ZU(gt4,3 + j2)) - 10*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt1,3 + j2))*ZU(gt4,3 + j2)) - 60*Qq*Qu*Sqr(gp)*SUM(j1,0,2,Conj(ZU(gt2,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt1,3 + j2))*ZU(gt4,3 + j2)) - 16*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt1,3 + j2))*ZU(gt4,3 + j2)) + 10*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt1,3 + j2))*ZU(gt4,3 + j2)) - 60*Sqr(gp)*Sqr(Qu)*SUM(j1,0,2,Conj(ZU(gt2,3 + j1))*ZU(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt1,3 + j2))*ZU(gt4,3 + j2)) + 30*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)) - 30*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt3,3 + j1))*SUM(j2,0,2,Conj(ZU(gt2,3 + j2))*ZU(gt4,3 + j2)) - 120*SUM(j2,0,2,Conj(ZU(gt1,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt4,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yu(j3,j4))*Conj(ZU(gt2,3 + j3)))*ZU(gt3,j4)) - 120*SUM(j2,0,2,Conj(ZU(gt2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt3,3 + j1)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yu(j3,j4))*Conj(ZU(gt1,3 + j3)))*ZU(gt4,j4)));

   return {result};
}

ScalarVertex VertexImpl<fields::Su, fields::Sv, typename fields::conj<fields::Su>::type, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 4>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const int gt4 = indices[3];
   const auto Qv = INPUTPARAMETER(Qv);
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qu = INPUTPARAMETER(Qu);
   const auto Ql = INPUTPARAMETER(Ql);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZU = MODELPARAMETER(ZU);

   const std::complex<double> result = 0.025*(-20*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gt2,3 + j1))*ZV(gt4,3 + j1))*(Qq*SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt3,j2)) + Qu*SUM(j2,0,2,Conj(ZU(gt1,3 + j2))*ZU(gt3,3 + j2))) + SUM(j1,0,2,Conj(ZV(gt2,j1))*ZV(gt4,j1))*((Sqr(g1) - 5*(Sqr(g2) + 4*Ql*Qq*Sqr(gp)))*SUM(j2,0,2,Conj(ZU(gt1,j2))*ZU(gt3,j2)) - 4*(Sqr(g1) + 5*Ql*Qu*Sqr(gp))*SUM(j2,0,2,Conj(ZU(gt1,3 + j2))*ZU(gt3,3 + j2))) + Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt4,j2)) - 5*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt4,j2)) - 20*Ql*Qq*Sqr(gp)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt4,j2)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt3,3 + j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt4,j2)) - 20*Ql*Qu*Sqr(gp)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt3,3 + j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt4,j2)) - 20*Qq*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt3,j1))*SUM(j2,0,2,Conj(ZV(gt2,3 + j2))*ZV(gt4,3 + j2)) - 20*Qu*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt3,3 + j1))*SUM(j2,0,2,Conj(ZV(gt2,3 + j2))*ZV(gt4,3 + j2)) - 40*SUM(j2,0,2,Conj(ZU(gt1,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gt3,3 + j1)))*SUM(j4,0,2,Conj(ZV(gt2,3 + j4))*SUM(j3,0,2,Conj(Yv(j3,j4))*ZV(gt4,j3))) - 40*SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(gt2,j1))*Yv(j1,j2))*ZV(gt4,3 + j2))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yu(j3,j4))*Conj(ZU(gt1,3 + j3)))*ZU(gt3,j4)));

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

InverseMetricVertex VertexImpl<fields::Su, typename fields::conj<fields::Sd>::type, fields::VWm, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.2357022603955158*g2*(6*gp*Qq*Cos(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZD(gt2,j1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Su, typename fields::conj<fields::Sd>::type, fields::VWm, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.2357022603955158*g2*(0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 6*gp*Qq*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZD(gt2,j1));

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

InverseMetricVertex VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VP, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qu = INPUTPARAMETER(Qu);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05555555555555555*((0.7745966692414834*g1*Cos(ThetaW) + 3*g2*Sin(ThetaW))*(6*gp*Qq*Cos(ThetaWp) - 3*g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1)) + 6.196773353931867*g1*(-3*gp*Qu*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Sin(2*ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VP, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qu = INPUTPARAMETER(Qu);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05555555555555555*(-((0.7745966692414834*g1*Cos(ThetaW) + 3*g2*Sin(ThetaW))*(-3*g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 6*gp*Qq*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1))) - 6.196773353931867*g1*(0.7745966692414834*g1*Cos(ThetaWp)*Sin(2*ThetaW) + 3*gp*Qu*Cos(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1)));

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

InverseMetricVertex VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VZ, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qu = INPUTPARAMETER(Qu);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.03333333333333333*((-15*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) - 7.745966692414834*g1*gp*Qq*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + Cos(ThetaWp)*Sin(ThetaWp)*(60*Sqr(gp)*Sqr(Qq) - Sqr(g1)*Sqr(Sin(ThetaW))) + 7.745966692414834*g1*gp*Qq*Sin(ThetaW)*Sqr(Sin(ThetaWp)) + 2*g2*Cos(ThetaW)*(3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) + 15*gp*Qq*Sqr(Cos(ThetaWp)) - 15*gp*Qq*Sqr(Sin(ThetaWp))))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1)) + 2*(7.745966692414834*g1*gp*Qu*Sin(ThetaW - 2*ThetaWp) + 7.745966692414834*g1*gp*Qu*Sin(ThetaW + 2*ThetaWp) - Sin(2*(ThetaW - ThetaWp))*Sqr(g1) - 2*Sin(2*ThetaWp)*Sqr(g1) + Sin(2*(ThetaW + ThetaWp))*Sqr(g1) + 15*Sin(2*ThetaWp)*Sqr(gp)*Sqr(Qu))*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qu = INPUTPARAMETER(Qu);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.05555555555555555*(Sqr(3*g2*Cos(ThetaW)*Cos(ThetaWp) - 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 6*gp*Qq*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1)) + 4*Sqr(1.5491933384829668*g1*Cos(ThetaWp)*Sin(ThetaW) + 3*gp*Qu*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VZp, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qu = INPUTPARAMETER(Qu);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05555555555555555*(Sqr(6*gp*Qq*Cos(ThetaWp) - 3*g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1)) + 4*Sqr(3*gp*Qu*Cos(ThetaWp) - 1.5491933384829668*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1)));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qu = INPUTPARAMETER(Qu);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.16666666666666666*(-((6*gp*Qq*Cos(ThetaWp) - 3*g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1))) - 2*(-3*gp*Qu*Cos(ThetaWp) + 1.5491933384829668*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qu = INPUTPARAMETER(Qu);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.16666666666666666*((-3*g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 6*gp*Qq*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1)) + 0.4*(7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 15*gp*Qu*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1)));

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
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qv = INPUTPARAMETER(Qv);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZV = MODELPARAMETER(ZV);

   const std::complex<double> result = 0.125*(-0.6*Sqr(g1)*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZV(gt4,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt3,j2)) - Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZV(gt4,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt3,j2)) - 4*Sqr(gp)*Sqr(Ql)*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZV(gt4,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt3,j2)) - 4*Ql*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gt1,3 + j1))*ZV(gt4,3 + j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt3,j2)) - 4*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gt2,3 + j1))*ZV(gt4,3 + j1))*(Ql*SUM(j2,0,2,Conj(ZV(gt1,j2))*ZV(gt3,j2)) + Qv*SUM(j2,0,2,Conj(ZV(gt1,3 + j2))*ZV(gt3,3 + j2))) - SUM(j1,0,2,Conj(ZV(gt2,j1))*ZV(gt4,j1))*((0.6*Sqr(g1) + Sqr(g2) + 4*Sqr(gp)*Sqr(Ql))*SUM(j2,0,2,Conj(ZV(gt1,j2))*ZV(gt3,j2)) + 4*Ql*Qv*Sqr(gp)*SUM(j2,0,2,Conj(ZV(gt1,3 + j2))*ZV(gt3,3 + j2))) - 4*Ql*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZV(gt4,j1))*SUM(j2,0,2,Conj(ZV(gt2,3 + j2))*ZV(gt3,3 + j2)) - 4*Sqr(gp)*Sqr(Qv)*SUM(j1,0,2,Conj(ZV(gt1,3 + j1))*ZV(gt4,3 + j1))*SUM(j2,0,2,Conj(ZV(gt2,3 + j2))*ZV(gt3,3 + j2)) - 0.6*Sqr(g1)*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZV(gt3,j1))*SUM(j2,0,2,Conj(ZV(gt1,j2))*ZV(gt4,j2)) - Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZV(gt3,j1))*SUM(j2,0,2,Conj(ZV(gt1,j2))*ZV(gt4,j2)) - 4*Sqr(gp)*Sqr(Ql)*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZV(gt3,j1))*SUM(j2,0,2,Conj(ZV(gt1,j2))*ZV(gt4,j2)) - 4*Ql*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gt2,3 + j1))*ZV(gt3,3 + j1))*SUM(j2,0,2,Conj(ZV(gt1,j2))*ZV(gt4,j2)) - 0.6*Sqr(g1)*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZV(gt3,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt4,j2)) - Sqr(g2)*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZV(gt3,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt4,j2)) - 4*Sqr(gp)*Sqr(Ql)*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZV(gt3,j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt4,j2)) - 4*Ql*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gt1,3 + j1))*ZV(gt3,3 + j1))*SUM(j2,0,2,Conj(ZV(gt2,j2))*ZV(gt4,j2)) - 4*Ql*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gt2,j1))*ZV(gt3,j1))*SUM(j2,0,2,Conj(ZV(gt1,3 + j2))*ZV(gt4,3 + j2)) - 4*Sqr(gp)*Sqr(Qv)*SUM(j1,0,2,Conj(ZV(gt2,3 + j1))*ZV(gt3,3 + j1))*SUM(j2,0,2,Conj(ZV(gt1,3 + j2))*ZV(gt4,3 + j2)) - 4*Ql*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZV(gt3,j1))*SUM(j2,0,2,Conj(ZV(gt2,3 + j2))*ZV(gt4,3 + j2)) - 4*Sqr(gp)*Sqr(Qv)*SUM(j1,0,2,Conj(ZV(gt1,3 + j1))*ZV(gt3,3 + j1))*SUM(j2,0,2,Conj(ZV(gt2,3 + j2))*ZV(gt4,3 + j2)) - 8*SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(gt2,j1))*Yv(j1,j2))*ZV(gt4,3 + j2))*SUM(j4,0,2,Conj(ZV(gt1,3 + j4))*SUM(j3,0,2,Conj(Yv(j3,j4))*ZV(gt3,j3))) - 8*SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(gt1,j1))*Yv(j1,j2))*ZV(gt4,3 + j2))*SUM(j4,0,2,Conj(ZV(gt2,3 + j4))*SUM(j3,0,2,Conj(Yv(j3,j4))*ZV(gt3,j3))) - 8*SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(gt2,j1))*Yv(j1,j2))*ZV(gt3,3 + j2))*SUM(j4,0,2,Conj(ZV(gt1,3 + j4))*SUM(j3,0,2,Conj(Yv(j3,j4))*ZV(gt4,j3))) - 8*SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(gt1,j1))*Yv(j1,j2))*ZV(gt3,3 + j2))*SUM(j4,0,2,Conj(ZV(gt2,3 + j4))*SUM(j3,0,2,Conj(Yv(j3,j4))*ZV(gt4,j3))));

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

InverseMetricVertex VertexImpl<fields::Sv, typename fields::conj<fields::Se>::type, fields::VWm, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.7071067811865475*g2*(2*gp*Ql*Cos(ThetaWp) - 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZE(gt2,j1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sv, typename fields::conj<fields::Se>::type, fields::VWm, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.7071067811865475*g2*(0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*Ql*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZE(gt2,j1));

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

InverseMetricVertex VertexImpl<fields::Sv, typename fields::conj<fields::Sv>::type, fields::VZ, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qv = INPUTPARAMETER(Qv);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.1*(5*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) - 7.745966692414834*g1*gp*Ql*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + Cos(ThetaWp)*Sin(ThetaWp)*(-20*Sqr(gp)*Sqr(Ql) + 3*Sqr(g1)*Sqr(Sin(ThetaW))) + 7.745966692414834*g1*gp*Ql*Sin(ThetaW)*Sqr(Sin(ThetaWp)) + 2*g2*Cos(ThetaW)*(3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) - 5*gp*Ql*Sqr(Cos(ThetaWp)) + 5*gp*Ql*Sqr(Sin(ThetaWp))))*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZV(gt2,j1)) + Sin(2*ThetaWp)*Sqr(gp)*Sqr(Qv)*SUM(j1,0,2,Conj(ZV(gt1,3 + j1))*ZV(gt2,3 + j1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sv, typename fields::conj<fields::Sv>::type, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qv = INPUTPARAMETER(Qv);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.5*(Sqr(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*Ql*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZV(gt2,j1)) + 4*Sqr(gp)*Sqr(Qv)*Sqr(Sin(ThetaWp))*SUM(j1,0,2,Conj(ZV(gt1,3 + j1))*ZV(gt2,3 + j1)));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Sv, typename fields::conj<fields::Sv>::type, fields::VZp, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qv = INPUTPARAMETER(Qv);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(Sqr(-2*gp*Ql*Cos(ThetaWp) + g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZV(gt2,j1)) + 4*Sqr(gp)*Sqr(Qv)*Sqr(Cos(ThetaWp))*SUM(j1,0,2,Conj(ZV(gt1,3 + j1))*ZV(gt2,3 + j1)));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::Sv, typename fields::conj<fields::Sv>::type, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qv = INPUTPARAMETER(Qv);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.1*(-10*gp*Ql*Cos(ThetaWp) + (5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZV(gt2,j1)) + gp*Qv*Cos(ThetaWp)*SUM(j1,0,2,Conj(ZV(gt1,3 + j1))*ZV(gt2,3 + j1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::Sv, typename fields::conj<fields::Sv>::type, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qv = INPUTPARAMETER(Qv);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -0.5*(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*Ql*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZV(gt2,j1)) + gp*Qv*Sin(ThetaWp)*SUM(j1,0,2,Conj(ZV(gt1,3 + j1))*ZV(gt2,3 + j1));

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

   const std::complex<double> result = std::complex<double>(0,-0.5)*g2*(Conj(ZA(gt1,0))*ZP(gt2,0) + Conj(ZA(gt1,1))*ZP(gt2,1));

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

   const std::complex<double> left = g2*Conj(UP(gt1,0))*ZN(gt2,2) - 0.7071067811865475*g2*Conj(UP(gt1,1))*ZN(gt2,4);

   const std::complex<double> right = g2*Conj(ZN(gt2,2))*UM(gt1,0) + 0.7071067811865475*g2*Conj(ZN(gt2,3))*UM(gt1,1);

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

   const std::complex<double> result = -0.5*g2*(Conj(ZH(gt1,0))*ZP(gt2,0) - Conj(ZH(gt1,1))*ZP(gt2,1));

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

   const std::complex<double> result = 0.5*(vd*Conj(ZH(gt1,0)) + vu*Conj(ZH(gt1,1)))*Sqr(g2);

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

InverseMetricVertex VertexImpl<fields::VWm, fields::VZp, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*g2*(vd*(2*gp*QHd*Cos(ThetaWp) - 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt1,0) + vu*(2*gp*QHu*Cos(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt1,1));

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

ChiralVertex VertexImpl<fields::VWm, typename fields::bar<fields::Cha>::type, fields::Chi>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZN = MODELPARAMETER(ZN);
   const auto UM = MODELPARAMETER(UM);
   const auto UP = MODELPARAMETER(UP);

   const std::complex<double> left = -0.5*g2*(2*Conj(ZN(gt2,2))*UM(gt1,0) + 1.4142135623730951*Conj(ZN(gt2,3))*UM(gt1,1));

   const std::complex<double> right = -(g2*Conj(UP(gt1,0))*ZN(gt2,2)) + 0.7071067811865475*g2*Conj(UP(gt1,1))*ZN(gt2,4);

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

   const std::complex<double> result = std::complex<double>(0,-0.5)*g2*(Conj(ZA(gt1,0))*ZP(gt2,0) + Conj(ZA(gt1,1))*ZP(gt2,1));

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

   const std::complex<double> result = -0.5*g2*(Conj(ZH(gt1,0))*ZP(gt2,0) - Conj(ZH(gt1,1))*ZP(gt2,1));

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
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*g2*(vd*(2*gp*QHd*Cos(ThetaWp) - 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt1,0) + vu*(2*gp*QHu*Cos(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt1,1));

   return {result};
}

InverseMetricVertex VertexImpl<fields::VWm, typename fields::conj<fields::Hpm>::type, fields::VZ>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*g2*(vd*(0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHd*Sin(ThetaWp))*ZP(gt1,0) + vu*(-0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHu*Sin(ThetaWp))*ZP(gt1,1));

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
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = std::complex<double>(0,-0.5)*(2*gp*Qs*Conj(ZA(gt1,2))*Conj(ZH(gt2,2))*Sin(ThetaWp) + Conj(ZA(gt1,0))*Conj(ZH(gt2,0))*(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHd*Sin(ThetaWp)) - Conj(ZA(gt1,1))*Conj(ZH(gt2,1))*(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*QHu*Sin(ThetaWp)));

   return {result, minuend_index, subtrahend_index};
}

ChiralVertex VertexImpl<fields::VZ, fields::Cha, typename fields::bar<fields::Cha>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto QHu = INPUTPARAMETER(QHu);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto UP = MODELPARAMETER(UP);
   const auto UM = MODELPARAMETER(UM);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = -(g2*Conj(UP(gt1,0))*Cos(ThetaW)*Cos(ThetaWp)*UP(gt2,0)) - 0.5*Conj(UP(gt1,1))*(g2*Cos(ThetaW)*Cos(ThetaWp) - 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHu*Sin(ThetaWp))*UP(gt2,1);

   const std::complex<double> right = -(g2*Conj(UM(gt2,0))*Cos(ThetaW)*Cos(ThetaWp)*UM(gt1,0)) + 0.1*Conj(UM(gt2,1))*(-5*g2*Cos(ThetaW)*Cos(ThetaWp) + 3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW) + 10*gp*QHd*Sin(ThetaWp))*UM(gt1,1);

   return {left, right};
}

ChiralVertex VertexImpl<fields::VZ, fields::Chi, fields::Chi>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto Qs = INPUTPARAMETER(Qs);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = 0.5*(-(Conj(ZN(gt2,3))*(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHd*Sin(ThetaWp))*ZN(gt1,3)) + Conj(ZN(gt2,4))*(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*QHu*Sin(ThetaWp))*ZN(gt1,4) - 2*gp*Qs*Conj(ZN(gt2,5))*Sin(ThetaWp)*ZN(gt1,5));

   const std::complex<double> right = 0.5*(Conj(ZN(gt1,3))*(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHd*Sin(ThetaWp))*ZN(gt2,3) - Conj(ZN(gt1,4))*(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*QHu*Sin(ThetaWp))*ZN(gt2,4) + 2*gp*Qs*Conj(ZN(gt1,5))*Sin(ThetaWp)*ZN(gt2,5));

   return {left, right};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, fields::hh, fields::Ah>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = std::complex<double>(0,-0.5)*(2*gp*Qs*Conj(ZA(gt1,2))*Conj(ZH(gt2,2))*Sin(ThetaWp) + Conj(ZA(gt1,0))*Conj(ZH(gt2,0))*(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHd*Sin(ThetaWp)) - Conj(ZA(gt1,1))*Conj(ZH(gt2,1))*(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*QHu*Sin(ThetaWp)));

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::VZ, fields::hh, fields::VZp>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto vS = MODELPARAMETER(vS);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vu = MODELPARAMETER(vu);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.1*(20*vS*Conj(ZH(gt1,2))*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gp)*Sqr(Qs) - vd*Conj(ZH(gt1,0))*(5*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) - 7.745966692414834*g1*gp*QHd*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + Cos(ThetaWp)*Sin(ThetaWp)*(-20*Sqr(gp)*Sqr(QHd) + 3*Sqr(g1)*Sqr(Sin(ThetaW))) + 7.745966692414834*g1*gp*QHd*Sin(ThetaW)*Sqr(Sin(ThetaWp)) + 2*g2*Cos(ThetaW)*(3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) - 5*gp*QHd*Sqr(Cos(ThetaWp)) + 5*gp*QHd*Sqr(Sin(ThetaWp)))) - vu*Conj(ZH(gt1,1))*(5*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 7.745966692414834*g1*gp*QHu*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + Cos(ThetaWp)*Sin(ThetaWp)*(-20*Sqr(gp)*Sqr(QHu) + 3*Sqr(g1)*Sqr(Sin(ThetaW))) - 7.745966692414834*g1*gp*QHu*Sin(ThetaW)*Sqr(Sin(ThetaWp)) + 2*g2*Cos(ThetaW)*(3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) + 5*gp*QHu*Sqr(Cos(ThetaWp)) - 5*gp*QHu*Sqr(Sin(ThetaWp)))));

   return {result};
}

InverseMetricVertex VertexImpl<fields::VZ, fields::hh, fields::VZ>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto vS = MODELPARAMETER(vS);
   const auto vd = MODELPARAMETER(vd);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(4*vS*Conj(ZH(gt1,2))*Sqr(gp)*Sqr(Qs)*Sqr(Sin(ThetaWp)) + vd*Conj(ZH(gt1,0))*Sqr(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHd*Sin(ThetaWp)) + vu*Conj(ZH(gt1,1))*Sqr(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*QHu*Sin(ThetaWp)));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, fields::Hpm, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.5*((g2*Cos(ThetaW)*Cos(ThetaWp) - 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*QHd*Sin(ThetaWp))*ZP(gt1,0)*ZP(gt2,0) + (g2*Cos(ThetaW)*Cos(ThetaWp) - 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHu*Sin(ThetaWp))*ZP(gt1,1)*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::VZ, fields::Hpm, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*g2*(vd*(0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHd*Sin(ThetaWp))*ZP(gt1,0) + vu*(-0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHu*Sin(ThetaWp))*ZP(gt1,1));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, fields::Sd, typename fields::conj<fields::Sd>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qd = INPUTPARAMETER(Qd);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.16666666666666666*((3*g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 6*gp*Qq*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1)) - 2*(0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 3*gp*Qd*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, fields::Se, typename fields::conj<fields::Se>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qe = INPUTPARAMETER(Qe);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.5*(g2*Cos(ThetaW)*Cos(ThetaWp) - 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*Ql*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) + (-0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + gp*Qe*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, fields::Su, typename fields::conj<fields::Su>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qu = INPUTPARAMETER(Qu);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.16666666666666666*((-3*g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 6*gp*Qq*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1)) + 0.4*(7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 15*gp*Qu*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, fields::Sv, typename fields::conj<fields::Sv>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qv = INPUTPARAMETER(Qv);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -0.5*(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*Ql*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZV(gt2,j1)) + gp*Qv*Sin(ThetaWp)*SUM(j1,0,2,Conj(ZV(gt1,3 + j1))*ZV(gt2,3 + j1));

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::VZ, fields::VZp, fields::hh>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto vS = MODELPARAMETER(vS);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vu = MODELPARAMETER(vu);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.1*(20*vS*Conj(ZH(gt1,2))*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gp)*Sqr(Qs) - vd*Conj(ZH(gt1,0))*(5*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) - 7.745966692414834*g1*gp*QHd*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + Cos(ThetaWp)*Sin(ThetaWp)*(-20*Sqr(gp)*Sqr(QHd) + 3*Sqr(g1)*Sqr(Sin(ThetaW))) + 7.745966692414834*g1*gp*QHd*Sin(ThetaW)*Sqr(Sin(ThetaWp)) + 2*g2*Cos(ThetaW)*(3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) - 5*gp*QHd*Sqr(Cos(ThetaWp)) + 5*gp*QHd*Sqr(Sin(ThetaWp)))) - vu*Conj(ZH(gt1,1))*(5*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 7.745966692414834*g1*gp*QHu*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + Cos(ThetaWp)*Sin(ThetaWp)*(-20*Sqr(gp)*Sqr(QHu) + 3*Sqr(g1)*Sqr(Sin(ThetaW))) - 7.745966692414834*g1*gp*QHu*Sin(ThetaW)*Sqr(Sin(ThetaWp)) + 2*g2*Cos(ThetaW)*(3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) + 5*gp*QHu*Sqr(Cos(ThetaWp)) - 5*gp*QHu*Sqr(Sin(ThetaWp)))));

   return {result};
}

ChiralVertex VertexImpl<fields::VZ, typename fields::bar<fields::Cha>::type, fields::Cha>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto UM = MODELPARAMETER(UM);
   const auto UP = MODELPARAMETER(UP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = g2*Conj(UM(gt2,0))*Cos(ThetaW)*Cos(ThetaWp)*UM(gt1,0) + 0.5*Conj(UM(gt2,1))*(g2*Cos(ThetaW)*Cos(ThetaWp) - 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*QHd*Sin(ThetaWp))*UM(gt1,1);

   const std::complex<double> right = g2*Conj(UP(gt1,0))*Cos(ThetaW)*Cos(ThetaWp)*UP(gt2,0) + 0.5*Conj(UP(gt1,1))*(g2*Cos(ThetaW)*Cos(ThetaWp) - 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHu*Sin(ThetaWp))*UP(gt2,1);

   return {left, right};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, typename fields::conj<fields::Hpm>::type, fields::Hpm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.5*((g2*Cos(ThetaW)*Cos(ThetaWp) - 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*QHd*Sin(ThetaWp))*ZP(gt1,0)*ZP(gt2,0) + (g2*Cos(ThetaW)*Cos(ThetaWp) - 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHu*Sin(ThetaWp))*ZP(gt1,1)*ZP(gt2,1));

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::VZ, typename fields::conj<fields::Hpm>::type, fields::VWm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*g2*(vd*(0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHd*Sin(ThetaWp))*ZP(gt1,0) + vu*(-0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHu*Sin(ThetaWp))*ZP(gt1,1));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, typename fields::conj<fields::Sd>::type, fields::Sd>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qd = INPUTPARAMETER(Qd);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.16666666666666666*((3*g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 6*gp*Qq*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,j1))*ZD(gt2,j1)) - 2*(0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 3*gp*Qd*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZD(gt1,3 + j1))*ZD(gt2,3 + j1)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, typename fields::conj<fields::Se>::type, fields::Se>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qe = INPUTPARAMETER(Qe);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.5*(g2*Cos(ThetaW)*Cos(ThetaWp) - 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*Ql*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) + (-0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + gp*Qe*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, typename fields::conj<fields::Su>::type, fields::Su>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qu = INPUTPARAMETER(Qu);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.16666666666666666*((-3*g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 6*gp*Qq*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,j1))*ZU(gt2,j1)) + 0.4*(7.745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 15*gp*Qu*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZU(gt1,3 + j1))*ZU(gt2,3 + j1)));

   return {result, minuend_index, subtrahend_index};
}

MomentumDifferenceVertex VertexImpl<fields::VZ, typename fields::conj<fields::Sv>::type, fields::Sv>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qv = INPUTPARAMETER(Qv);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -0.5*(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*Ql*Sin(ThetaWp))*SUM(j1,0,2,Conj(ZV(gt1,j1))*ZV(gt2,j1)) + gp*Qv*Sin(ThetaWp)*SUM(j1,0,2,Conj(ZV(gt1,3 + j1))*ZV(gt2,3 + j1));

   return {result, minuend_index, subtrahend_index};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha>::type, fields::Cha, fields::Ah>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto UM = MODELPARAMETER(UM);
   const auto UP = MODELPARAMETER(UP);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,0.7071067811865475)*(g2*Conj(UM(gt2,0))*Conj(UP(gt1,1))*Conj(ZA(gt3,1)) + Conj(UM(gt2,1))*(g2*Conj(UP(gt1,0))*Conj(ZA(gt3,0)) - Conj(UP(gt1,1))*Conj(ZA(gt3,2))*Lambdax));

   const std::complex<double> right = std::complex<double>(0.,-0.7071067811865475)*(g2*Conj(ZA(gt3,0))*UM(gt1,1)*UP(gt2,0) + (g2*Conj(ZA(gt3,1))*UM(gt1,0) - Conj(Lambdax)*Conj(ZA(gt3,2))*UM(gt1,1))*UP(gt2,1));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha>::type, fields::Cha, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto UM = MODELPARAMETER(UM);
   const auto UP = MODELPARAMETER(UP);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*(g2*Conj(UM(gt2,0))*Conj(UP(gt1,1))*Conj(ZH(gt3,1)) + Conj(UM(gt2,1))*(g2*Conj(UP(gt1,0))*Conj(ZH(gt3,0)) + Conj(UP(gt1,1))*Conj(ZH(gt3,2))*Lambdax));

   const std::complex<double> right = -0.7071067811865475*(g2*Conj(ZH(gt3,0))*UM(gt1,1)*UP(gt2,0) + (g2*Conj(ZH(gt3,1))*UM(gt1,0) + Conj(Lambdax)*Conj(ZH(gt3,2))*UM(gt1,1))*UP(gt2,1));

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

ChiralVertex VertexImpl<typename fields::bar<fields::Cha>::type, fields::Cha, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto UM = MODELPARAMETER(UM);
   const auto UP = MODELPARAMETER(UP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = -(g2*Conj(UM(gt2,0))*Cos(ThetaW)*Sin(ThetaWp)*UM(gt1,0)) - 0.5*Conj(UM(gt2,1))*(2*gp*QHd*Cos(ThetaWp) + g2*Cos(ThetaW)*Sin(ThetaWp) - 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*UM(gt1,1);

   const std::complex<double> right = -(g2*Conj(UP(gt1,0))*Cos(ThetaW)*Sin(ThetaWp)*UP(gt2,0)) + 0.1*Conj(UP(gt1,1))*(10*gp*QHu*Cos(ThetaWp) + (-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW))*Sin(ThetaWp))*UP(gt2,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha>::type, fields::Cha, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto UM = MODELPARAMETER(UM);
   const auto UP = MODELPARAMETER(UP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = g2*Conj(UM(gt2,0))*Cos(ThetaW)*Cos(ThetaWp)*UM(gt1,0) + 0.5*Conj(UM(gt2,1))*(g2*Cos(ThetaW)*Cos(ThetaWp) - 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*QHd*Sin(ThetaWp))*UM(gt1,1);

   const std::complex<double> right = g2*Conj(UP(gt1,0))*Cos(ThetaW)*Cos(ThetaWp)*UP(gt2,0) + 0.5*Conj(UP(gt1,1))*(g2*Cos(ThetaW)*Cos(ThetaWp) - 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHu*Sin(ThetaWp))*UP(gt2,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Cha>::type, fields::Chi, fields::Hpm>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto QHu = INPUTPARAMETER(QHu);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto g2 = MODELPARAMETER(g2);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto UP = MODELPARAMETER(UP);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZP = MODELPARAMETER(ZP);
   const auto UM = MODELPARAMETER(UM);

   const std::complex<double> left = -(g2*Conj(UP(gt1,0))*Conj(ZN(gt2,4))*ZP(gt3,1)) - 0.5*Conj(UP(gt1,1))*(2*Conj(ZN(gt2,5))*Lambdax*ZP(gt3,0) + 1.4142135623730951*(2*gp*QHu*Conj(ZN(gt2,0)) + 0.7745966692414834*g1*Conj(ZN(gt2,1)) + g2*Conj(ZN(gt2,2)))*ZP(gt3,1));

   const std::complex<double> right = -(g2*UM(gt1,0)*ZN(gt2,3)*ZP(gt3,0)) + UM(gt1,1)*(-1.4142135623730951*gp*QHd*ZN(gt2,0)*ZP(gt3,0) + 0.5477225575051661*g1*ZN(gt2,1)*ZP(gt3,0) + 0.7071067811865475*g2*ZN(gt2,2)*ZP(gt3,0) - Conj(Lambdax)*ZN(gt2,5)*ZP(gt3,1));

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

   const std::complex<double> left = -0.5*g2*(2*Conj(ZN(gt2,2))*UM(gt1,0) + 1.4142135623730951*Conj(ZN(gt2,3))*UM(gt1,1));

   const std::complex<double> right = -(g2*Conj(UP(gt1,0))*ZN(gt2,2)) + 0.7071067811865475*g2*Conj(UP(gt1,1))*ZN(gt2,4);

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

   const std::complex<double> left = -(g2*Conj(UP(gt1,0))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZV(gt3,j1))) + Conj(UP(gt1,1))*SUM(j2,0,2,SUM(j1,0,2,Conj(ZEL(gt2,j1))*Yv(j1,j2))*ZV(gt3,3 + j2));

   const std::complex<double> right = SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZV(gt3,j2))*UM(gt1,1);

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

ChiralVertex VertexImpl<typename fields::bar<fields::Cha>::type, typename fields::bar<fields::Fv>::type, fields::Se>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Yv = MODELPARAMETER(Yv);
   const auto Ye = MODELPARAMETER(Ye);
   const auto UP = MODELPARAMETER(UP);
   const auto ZVR = MODELPARAMETER(ZVR);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZVL = MODELPARAMETER(ZVL);
   const auto UM = MODELPARAMETER(UM);

   const std::complex<double> left = Conj(UP(gt1,1))*SUM(j2,0,2,Conj(ZVR(gt2,j2))*SUM(j1,0,2,Conj(ZE(gt3,j1))*Yv(j1,j2)));

   const std::complex<double> right = -(g2*SUM(j1,0,2,Conj(ZE(gt3,j1))*ZVL(gt2,j1))*UM(gt1,0)) + SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt3,3 + j1)))*ZVL(gt2,j2))*UM(gt1,1);

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
   const auto Qd = INPUTPARAMETER(Qd);
   const auto Qq = INPUTPARAMETER(Qq);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZD = MODELPARAMETER(ZD);
   const auto ZDR = MODELPARAMETER(ZDR);
   const auto ZDL = MODELPARAMETER(ZDL);

   const std::complex<double> left = -1.4142135623730951*gp*Qd*Conj(ZN(gt2,0))*SUM(j1,0,2,Conj(ZD(gt3,3 + j1))*Conj(ZDR(gt1,j1))) - 0.3651483716701107*g1*Conj(ZN(gt2,1))*SUM(j1,0,2,Conj(ZD(gt3,3 + j1))*Conj(ZDR(gt1,j1))) - Conj(ZN(gt2,3))*SUM(j2,0,2,Conj(ZD(gt3,j2))*SUM(j1,0,2,Conj(ZDR(gt1,j1))*Yd(j1,j2)));

   const std::complex<double> right = -0.2357022603955158*SUM(j1,0,2,Conj(ZD(gt3,j1))*ZDL(gt1,j1))*(6*gp*Qq*ZN(gt2,0) + 0.7745966692414834*g1*ZN(gt2,1) - 3*g2*ZN(gt2,2)) - SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gt3,3 + j1)))*ZDL(gt1,j2))*ZN(gt2,3);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fd, fields::Ah>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto ZDR = MODELPARAMETER(ZDR);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gt3,0))*SUM(j2,0,2,Conj(ZDL(gt2,j2))*SUM(j1,0,2,Conj(ZDR(gt1,j1))*Yd(j1,j2)));

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*Conj(ZA(gt3,0))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*ZDR(gt2,j1))*ZDL(gt1,j2));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fd, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Yd = MODELPARAMETER(Yd);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZDL = MODELPARAMETER(ZDL);
   const auto ZDR = MODELPARAMETER(ZDR);

   const std::complex<double> left = -0.7071067811865475*Conj(ZH(gt3,0))*SUM(j2,0,2,Conj(ZDL(gt2,j2))*SUM(j1,0,2,Conj(ZDR(gt1,j1))*Yd(j1,j2)));

   const std::complex<double> right = -0.7071067811865475*Conj(ZH(gt3,0))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*ZDR(gt2,j1))*ZDL(gt1,j2));

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

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fd, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qd = INPUTPARAMETER(Qd);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.16666666666666666*KroneckerDelta(gt1,gt2)*(6*gp*Qq*Cos(ThetaWp) + 3*g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp));

   const std::complex<double> right = KroneckerDelta(gt1,gt2)*(gp*Qd*Cos(ThetaWp) + 0.2581988897471611*g1*Sin(ThetaW)*Sin(ThetaWp));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fd, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qd = INPUTPARAMETER(Qd);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = 0.16666666666666666*KroneckerDelta(gt1,gt2)*(3*g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 6*gp*Qq*Sin(ThetaWp));

   const std::complex<double> right = KroneckerDelta(gt1,gt2)*(-0.2581988897471611*g1*Cos(ThetaWp)*Sin(ThetaW) + gp*Qd*Sin(ThetaWp));

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

   const std::complex<double> right = -(g2*SUM(j1,0,2,Conj(ZV(gt3,j1))*ZEL(gt1,j1))*UP(gt2,0)) + SUM(j2,0,2,Conj(ZV(gt3,3 + j2))*SUM(j1,0,2,Conj(Yv(j1,j2))*ZEL(gt1,j1)))*UP(gt2,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Chi, fields::Se>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Qe = INPUTPARAMETER(Qe);
   const auto Ql = INPUTPARAMETER(Ql);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZEL = MODELPARAMETER(ZEL);

   const std::complex<double> left = -1.4142135623730951*gp*Qe*Conj(ZN(gt2,0))*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*Conj(ZER(gt1,j1))) - 1.0954451150103321*g1*Conj(ZN(gt2,1))*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*Conj(ZER(gt1,j1))) - Conj(ZN(gt2,3))*SUM(j2,0,2,Conj(ZE(gt3,j2))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,j2)));

   const std::complex<double> right = 0.1414213562373095*SUM(j1,0,2,Conj(ZE(gt3,j1))*ZEL(gt1,j1))*(-10*gp*Ql*ZN(gt2,0) + 3.872983346207417*g1*ZN(gt2,1) + 5*g2*ZN(gt2,2)) - SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt3,3 + j1)))*ZEL(gt1,j2))*ZN(gt2,3);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fe, fields::Ah>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZER = MODELPARAMETER(ZER);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gt3,0))*SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,j2)));

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*Conj(ZA(gt3,0))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZEL(gt1,j2));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fe, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZER = MODELPARAMETER(ZER);

   const std::complex<double> left = -0.7071067811865475*Conj(ZH(gt3,0))*SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,j2)));

   const std::complex<double> right = -0.7071067811865475*Conj(ZH(gt3,0))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZEL(gt1,j2));

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

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fe, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qe = INPUTPARAMETER(Qe);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.5*KroneckerDelta(gt1,gt2)*(2*gp*Ql*Cos(ThetaWp) + g2*Cos(ThetaW)*Sin(ThetaWp) - 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp));

   const std::complex<double> right = KroneckerDelta(gt1,gt2)*(gp*Qe*Cos(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fe, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qe = INPUTPARAMETER(Qe);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = 0.5*KroneckerDelta(gt1,gt2)*(g2*Cos(ThetaW)*Cos(ThetaWp) - 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*Ql*Sin(ThetaWp));

   const std::complex<double> right = -(KroneckerDelta(gt1,gt2)*(0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - gp*Qe*Sin(ThetaWp)));

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
   const auto ZVL = MODELPARAMETER(ZVL);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZVR = MODELPARAMETER(ZVR);

   const std::complex<double> left = SUM(j2,0,2,Conj(ZVL(gt2,j2))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,j2)))*ZP(gt3,0);

   const std::complex<double> right = SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j2))*ZEL(gt1,j1))*ZVR(gt2,j2))*ZP(gt3,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fv, fields::VWm>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZVL = MODELPARAMETER(ZVL);
   const auto ZEL = MODELPARAMETER(ZEL);

   const std::complex<double> left = -0.7071067811865475*g2*SUM(j1,0,2,Conj(ZVL(gt2,j1))*ZEL(gt1,j1));

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
   const auto ZVL = MODELPARAMETER(ZVL);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZVR = MODELPARAMETER(ZVR);

   const std::complex<double> left = SUM(j2,0,2,Conj(ZVL(gt2,j2))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,j2)))*ZP(gt3,0);

   const std::complex<double> right = SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j2))*ZEL(gt1,j1))*ZVR(gt2,j2))*ZP(gt3,1);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Se, fields::Chi>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Qe = INPUTPARAMETER(Qe);
   const auto Ql = INPUTPARAMETER(Ql);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZEL = MODELPARAMETER(ZEL);

   const std::complex<double> left = -1.4142135623730951*gp*Qe*Conj(ZN(gt2,0))*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*Conj(ZER(gt1,j1))) - 1.0954451150103321*g1*Conj(ZN(gt2,1))*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*Conj(ZER(gt1,j1))) - Conj(ZN(gt2,3))*SUM(j2,0,2,Conj(ZE(gt3,j2))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,j2)));

   const std::complex<double> right = 0.1414213562373095*SUM(j1,0,2,Conj(ZE(gt3,j1))*ZEL(gt1,j1))*(-10*gp*Ql*ZN(gt2,0) + 3.872983346207417*g1*ZN(gt2,1) + 5*g2*ZN(gt2,2)) - SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt3,3 + j1)))*ZEL(gt1,j2))*ZN(gt2,3);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, fields::Chi, fields::Su>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Qu = INPUTPARAMETER(Qu);
   const auto Qq = INPUTPARAMETER(Qq);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZU = MODELPARAMETER(ZU);
   const auto ZUR = MODELPARAMETER(ZUR);
   const auto ZUL = MODELPARAMETER(ZUL);

   const std::complex<double> left = -1.4142135623730951*gp*Qu*Conj(ZN(gt2,0))*SUM(j1,0,2,Conj(ZU(gt3,3 + j1))*Conj(ZUR(gt1,j1))) + 0.7302967433402214*g1*Conj(ZN(gt2,1))*SUM(j1,0,2,Conj(ZU(gt3,3 + j1))*Conj(ZUR(gt1,j1))) - Conj(ZN(gt2,4))*SUM(j2,0,2,Conj(ZU(gt3,j2))*SUM(j1,0,2,Conj(ZUR(gt1,j1))*Yu(j1,j2)));

   const std::complex<double> right = -0.2357022603955158*SUM(j1,0,2,Conj(ZU(gt3,j1))*ZUL(gt1,j1))*(6*gp*Qq*ZN(gt2,0) + 0.7745966692414834*g1*ZN(gt2,1) + 3*g2*ZN(gt2,2)) - SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gt3,3 + j1)))*ZUL(gt1,j2))*ZN(gt2,4);

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
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto ZUR = MODELPARAMETER(ZUR);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gt3,1))*SUM(j2,0,2,Conj(ZUL(gt2,j2))*SUM(j1,0,2,Conj(ZUR(gt1,j1))*Yu(j1,j2)));

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*Conj(ZA(gt3,1))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*ZUR(gt2,j1))*ZUL(gt1,j2));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fu, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Yu = MODELPARAMETER(Yu);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZUL = MODELPARAMETER(ZUL);
   const auto ZUR = MODELPARAMETER(ZUR);

   const std::complex<double> left = -0.7071067811865475*Conj(ZH(gt3,1))*SUM(j2,0,2,Conj(ZUL(gt2,j2))*SUM(j1,0,2,Conj(ZUR(gt1,j1))*Yu(j1,j2)));

   const std::complex<double> right = -0.7071067811865475*Conj(ZH(gt3,1))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*ZUR(gt2,j1))*ZUL(gt1,j2));

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

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fu, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qu = INPUTPARAMETER(Qu);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.16666666666666666*KroneckerDelta(gt1,gt2)*(6*gp*Qq*Cos(ThetaWp) - 3*g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp));

   const std::complex<double> right = KroneckerDelta(gt1,gt2)*(gp*Qu*Cos(ThetaWp) - 0.5163977794943222*g1*Sin(ThetaW)*Sin(ThetaWp));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fu, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Qq = INPUTPARAMETER(Qq);
   const auto Qu = INPUTPARAMETER(Qu);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = -0.16666666666666666*KroneckerDelta(gt1,gt2)*(3*g2*Cos(ThetaW)*Cos(ThetaWp) - 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 6*gp*Qq*Sin(ThetaWp));

   const std::complex<double> right = KroneckerDelta(gt1,gt2)*(0.5163977794943222*g1*Cos(ThetaWp)*Sin(ThetaW) + gp*Qu*Sin(ThetaWp));

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

ChiralVertex VertexImpl<typename fields::bar<fields::Fv>::type, fields::Chi, fields::Sv>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Qv = INPUTPARAMETER(Qv);
   const auto Ql = INPUTPARAMETER(Ql);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZVR = MODELPARAMETER(ZVR);
   const auto ZVL = MODELPARAMETER(ZVL);

   const std::complex<double> left = -1.4142135623730951*gp*Qv*Conj(ZN(gt2,0))*SUM(j1,0,2,Conj(ZV(gt3,3 + j1))*Conj(ZVR(gt1,j1))) - Conj(ZN(gt2,4))*SUM(j2,0,2,Conj(ZVR(gt1,j2))*SUM(j1,0,2,Conj(ZV(gt3,j1))*Yv(j1,j2)));

   const std::complex<double> right = -0.7071067811865475*SUM(j1,0,2,Conj(ZV(gt3,j1))*ZVL(gt1,j1))*(2*gp*Ql*ZN(gt2,0) - 0.7745966692414834*g1*ZN(gt2,1) + g2*ZN(gt2,2)) - SUM(j2,0,2,Conj(ZV(gt3,3 + j2))*SUM(j1,0,2,Conj(Yv(j1,j2))*ZVL(gt1,j1)))*ZN(gt2,4);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fv>::type, fields::Fe, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Yv = MODELPARAMETER(Yv);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZVR = MODELPARAMETER(ZVR);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZVL = MODELPARAMETER(ZVL);

   const std::complex<double> left = SUM(j2,0,2,Conj(ZVR(gt1,j2))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*Yv(j1,j2)))*ZP(gt3,1);

   const std::complex<double> right = SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZVL(gt1,j2))*ZP(gt3,0);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fv>::type, fields::Fe, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZVL = MODELPARAMETER(ZVL);

   const std::complex<double> left = -0.7071067811865475*g2*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZVL(gt1,j1));

   const std::complex<double> right = 0;

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fv>::type, fields::Fv, fields::Ah>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZA = MODELPARAMETER(ZA);
   const auto ZVR = MODELPARAMETER(ZVR);
   const auto ZVL = MODELPARAMETER(ZVL);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gt3,1))*SUM(j2,0,2,Conj(ZVR(gt1,j2))*SUM(j1,0,2,Conj(ZVL(gt2,j1))*Yv(j1,j2)));

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*Conj(ZA(gt3,1))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j2))*ZVL(gt1,j1))*ZVR(gt2,j2));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fv>::type, fields::Fv, fields::hh>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Yv = MODELPARAMETER(Yv);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ZVR = MODELPARAMETER(ZVR);
   const auto ZVL = MODELPARAMETER(ZVL);

   const std::complex<double> left = -0.7071067811865475*Conj(ZH(gt3,1))*SUM(j2,0,2,Conj(ZVR(gt1,j2))*SUM(j1,0,2,Conj(ZVL(gt2,j1))*Yv(j1,j2)));

   const std::complex<double> right = -0.7071067811865475*Conj(ZH(gt3,1))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j2))*ZVL(gt1,j1))*ZVR(gt2,j2));

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fv>::type, fields::Fv, fields::VZp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qv = INPUTPARAMETER(Qv);
   const auto gp = MODELPARAMETER(gp);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.5*KroneckerDelta(gt1,gt2)*(2*gp*Ql*Cos(ThetaWp) - (g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW))*Sin(ThetaWp));

   const std::complex<double> right = gp*Qv*Cos(ThetaWp)*KroneckerDelta(gt1,gt2);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fv>::type, fields::Fv, fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ql = INPUTPARAMETER(Ql);
   const auto Qv = INPUTPARAMETER(Qv);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> left = -0.5*KroneckerDelta(gt1,gt2)*(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*Ql*Sin(ThetaWp));

   const std::complex<double> right = gp*Qv*KroneckerDelta(gt1,gt2)*Sin(ThetaWp);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fv>::type, typename fields::conj<fields::Hpm>::type, fields::Fe>::evaluate(
   const std::array<int, 3>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Yv = MODELPARAMETER(Yv);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZVR = MODELPARAMETER(ZVR);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZVL = MODELPARAMETER(ZVL);

   const std::complex<double> left = SUM(j2,0,2,Conj(ZVR(gt1,j2))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*Yv(j1,j2)))*ZP(gt3,1);

   const std::complex<double> right = SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZVL(gt1,j2))*ZP(gt3,0);

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

   const std::complex<double> result = std::complex<double>(0,0.25)*(vd*Conj(ZA(gt3,0)) - vu*Conj(ZA(gt3,1)))*Sqr(g2);

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

   const std::complex<double> result = -0.25*(vd*Conj(ZH(gt3,0)) + vu*Conj(ZH(gt3,1)))*Sqr(g2);

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

MomentumVertex VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gWmC, fields::VZp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -(g2*Cos(ThetaW)*Sin(ThetaWp));

   return {result, 1};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gWmC, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = g2*Cos(ThetaW)*Cos(ThetaWp);

   return {result, 1};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gZp, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.25*g2*(vd*(2*gp*QHd*Cos(ThetaWp) + g2*Cos(ThetaW)*Sin(ThetaWp) - 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt3,0) + vu*(2*gp*QHu*Cos(ThetaWp) - g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt3,1));

   return {result};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gZp, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = g2*Cos(ThetaW)*Sin(ThetaWp);

   return {result, 1};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gZ, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.25*g2*(vd*(g2*Cos(ThetaW)*Cos(ThetaWp) - 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*QHd*Sin(ThetaWp))*ZP(gt3,0) + vu*(-(g2*Cos(ThetaW)*Cos(ThetaWp)) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*QHu*Sin(ThetaWp))*ZP(gt3,1));

   return {result};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gZ, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -(g2*Cos(ThetaW)*Cos(ThetaWp));

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

   const std::complex<double> result = std::complex<double>(0,-0.25)*(vd*Conj(ZA(gt3,0)) - vu*Conj(ZA(gt3,1)))*Sqr(g2);

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

   const std::complex<double> result = -0.25*(vd*Conj(ZH(gt3,0)) + vu*Conj(ZH(gt3,1)))*Sqr(g2);

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

MomentumVertex VertexImpl<typename fields::bar<fields::gWm>::type, fields::gWm, fields::VZp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = g2*Cos(ThetaW)*Sin(ThetaWp);

   return {result, 1};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWm>::type, fields::gWm, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -(g2*Cos(ThetaW)*Cos(ThetaWp));

   return {result, 1};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gWm>::type, fields::gZ, fields::Hpm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = 0.25*g2*(vd*(g2*Cos(ThetaW)*Cos(ThetaWp) - 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*QHd*Sin(ThetaWp))*ZP(gt3,0) + vu*(-(g2*Cos(ThetaW)*Cos(ThetaWp)) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*QHu*Sin(ThetaWp))*ZP(gt3,1));

   return {result};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWm>::type, fields::gZ, fields::VWm>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = g2*Cos(ThetaW)*Cos(ThetaWp);

   return {result, 1};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gWm>::type, fields::gZp, fields::Hpm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.25*g2*(vd*(2*gp*QHd*Cos(ThetaWp) + g2*Cos(ThetaW)*Sin(ThetaWp) - 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt3,0) + vu*(2*gp*QHu*Cos(ThetaWp) - g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt3,1));

   return {result};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWm>::type, fields::gZp, fields::VWm>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -(g2*Cos(ThetaW)*Sin(ThetaWp));

   return {result, 1};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gZp>::type, fields::gWmC, fields::Hpm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.25*g2*(vd*(-2*gp*QHd*Cos(ThetaWp) + g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt3,0) - vu*(2*gp*QHu*Cos(ThetaWp) + g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt3,1));

   return {result};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gZp>::type, fields::gWmC, fields::VWm>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = g2*Cos(ThetaW)*Sin(ThetaWp);

   return {result, 1};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gZp>::type, fields::gWm, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.25*g2*(vd*(-2*gp*QHd*Cos(ThetaWp) + g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt3,0) - vu*(2*gp*QHu*Cos(ThetaWp) + g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt3,1));

   return {result};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gZp>::type, fields::gWm, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -(g2*Cos(ThetaW)*Sin(ThetaWp));

   return {result, 1};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gZp>::type, fields::gZ, fields::hh>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto vS = MODELPARAMETER(vS);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vu = MODELPARAMETER(vu);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*(-20*vS*Conj(ZH(gt3,2))*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gp)*Sqr(Qs) + vd*Conj(ZH(gt3,0))*(5*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) - 7.745966692414834*g1*gp*QHd*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + Cos(ThetaWp)*Sin(ThetaWp)*(-20*Sqr(gp)*Sqr(QHd) + 3*Sqr(g1)*Sqr(Sin(ThetaW))) + 7.745966692414834*g1*gp*QHd*Sin(ThetaW)*Sqr(Sin(ThetaWp)) + 2*g2*Cos(ThetaW)*(3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) - 5*gp*QHd*Sqr(Cos(ThetaWp)) + 5*gp*QHd*Sqr(Sin(ThetaWp)))) + vu*Conj(ZH(gt3,1))*(5*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 7.745966692414834*g1*gp*QHu*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + Cos(ThetaWp)*Sin(ThetaWp)*(-20*Sqr(gp)*Sqr(QHu) + 3*Sqr(g1)*Sqr(Sin(ThetaW))) - 7.745966692414834*g1*gp*QHu*Sin(ThetaW)*Sqr(Sin(ThetaWp)) + 2*g2*Cos(ThetaW)*(3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) + 5*gp*QHu*Sqr(Cos(ThetaWp)) - 5*gp*QHu*Sqr(Sin(ThetaWp)))));

   return {result};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gZp>::type, fields::gZp, fields::hh>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto vS = MODELPARAMETER(vS);
   const auto vd = MODELPARAMETER(vd);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.25*(-4*vS*Conj(ZH(gt3,2))*Sqr(gp)*Sqr(Qs)*Sqr(Cos(ThetaWp)) - vd*Conj(ZH(gt3,0))*Sqr(-2*gp*QHd*Cos(ThetaWp) + g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp)) - vu*Conj(ZH(gt3,1))*Sqr(2*gp*QHu*Cos(ThetaWp) + g2*Cos(ThetaW)*Sin(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp)));

   return {result};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gWmC, fields::Hpm>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -0.25*g2*(vd*(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHd*Sin(ThetaWp))*ZP(gt3,0) - vu*(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*QHu*Sin(ThetaWp))*ZP(gt3,1));

   return {result};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gWmC, fields::VWm>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -(g2*Cos(ThetaW)*Cos(ThetaWp));

   return {result, 1};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gWm, typename fields::conj<fields::Hpm>::type>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -0.25*g2*(vd*(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHd*Sin(ThetaWp))*ZP(gt3,0) - vu*(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*QHu*Sin(ThetaWp))*ZP(gt3,1));

   return {result};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gWm, typename fields::conj<fields::VWm>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = g2*Cos(ThetaW)*Cos(ThetaWp);

   return {result, 1};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gZ, fields::hh>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto vS = MODELPARAMETER(vS);
   const auto vd = MODELPARAMETER(vd);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.25*(-4*vS*Conj(ZH(gt3,2))*Sqr(gp)*Sqr(Qs)*Sqr(Sin(ThetaWp)) - vd*Conj(ZH(gt3,0))*Sqr(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHd*Sin(ThetaWp)) - vu*Conj(ZH(gt3,1))*Sqr(g2*Cos(ThetaW)*Cos(ThetaWp) + 0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) - 2*gp*QHu*Sin(ThetaWp)));

   return {result};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gZp, fields::hh>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt3 = indices[0];
   const auto Qs = INPUTPARAMETER(Qs);
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto gp = MODELPARAMETER(gp);
   const auto vS = MODELPARAMETER(vS);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto vu = MODELPARAMETER(vu);
   const auto ZH = MODELPARAMETER(ZH);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*(-20*vS*Conj(ZH(gt3,2))*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(gp)*Sqr(Qs) + vd*Conj(ZH(gt3,0))*(5*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) - 7.745966692414834*g1*gp*QHd*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + Cos(ThetaWp)*Sin(ThetaWp)*(-20*Sqr(gp)*Sqr(QHd) + 3*Sqr(g1)*Sqr(Sin(ThetaW))) + 7.745966692414834*g1*gp*QHd*Sin(ThetaW)*Sqr(Sin(ThetaWp)) + 2*g2*Cos(ThetaW)*(3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) - 5*gp*QHd*Sqr(Cos(ThetaWp)) + 5*gp*QHd*Sqr(Sin(ThetaWp)))) + vu*Conj(ZH(gt3,1))*(5*Cos(ThetaWp)*Sin(ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)) + 7.745966692414834*g1*gp*QHu*Sin(ThetaW)*Sqr(Cos(ThetaWp)) + Cos(ThetaWp)*Sin(ThetaWp)*(-20*Sqr(gp)*Sqr(QHu) + 3*Sqr(g1)*Sqr(Sin(ThetaW))) - 7.745966692414834*g1*gp*QHu*Sin(ThetaW)*Sqr(Sin(ThetaWp)) + 2*g2*Cos(ThetaW)*(3.872983346207417*g1*Cos(ThetaWp)*Sin(ThetaW)*Sin(ThetaWp) + 5*gp*QHu*Sqr(Cos(ThetaWp)) - 5*gp*QHu*Sqr(Sin(ThetaWp)))));

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

InverseMetricVertex VertexImpl<typename fields::conj<fields::Hpm>::type, fields::VWm, fields::VZp>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto gp = MODELPARAMETER(gp);
   const auto g1 = MODELPARAMETER(g1);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*g2*(vd*(2*gp*QHd*Cos(ThetaWp) - 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt1,0) + vu*(2*gp*QHu*Cos(ThetaWp) + 0.7745966692414834*g1*Sin(ThetaW)*Sin(ThetaWp))*ZP(gt1,1));

   return {result};
}

InverseMetricVertex VertexImpl<typename fields::conj<fields::Hpm>::type, fields::VWm, fields::VZ>::evaluate(
   const std::array<int, 1>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const auto QHd = INPUTPARAMETER(QHd);
   const auto QHu = INPUTPARAMETER(QHu);
   const auto g2 = MODELPARAMETER(g2);
   const auto vd = MODELPARAMETER(vd);
   const auto g1 = MODELPARAMETER(g1);
   const auto gp = MODELPARAMETER(gp);
   const auto vu = MODELPARAMETER(vu);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*g2*(vd*(0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHd*Sin(ThetaWp))*ZP(gt1,0) + vu*(-0.7745966692414834*g1*Cos(ThetaWp)*Sin(ThetaW) + 2*gp*QHu*Sin(ThetaWp))*ZP(gt1,1));

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

   const std::complex<double> left = -(g2*Conj(UP(gt1,0))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZV(gt3,j1))) + Conj(UP(gt1,1))*SUM(j2,0,2,SUM(j1,0,2,Conj(ZEL(gt2,j1))*Yv(j1,j2))*ZV(gt3,3 + j2));

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

QuadrupleVectorVertex VertexImpl<typename fields::conj<fields::VWm>::type, fields::VP, fields::VWm, fields::VZp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> part1 = -(Cos(ThetaW)*Sin(ThetaW)*Sin(ThetaWp)*Sqr(g2));

   const std::complex<double> part2 = Sin(2*ThetaW)*Sin(ThetaWp)*Sqr(g2);

   const std::complex<double> part3 = -(Cos(ThetaW)*Sin(ThetaW)*Sin(ThetaWp)*Sqr(g2));

   return {part1, part2, part3};
}

QuadrupleVectorVertex VertexImpl<typename fields::conj<fields::VWm>::type, fields::VP, fields::VWm, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> part1 = Cos(ThetaW)*Cos(ThetaWp)*Sin(ThetaW)*Sqr(g2);

   const std::complex<double> part2 = -(Cos(ThetaWp)*Sin(2*ThetaW)*Sqr(g2));

   const std::complex<double> part3 = Cos(ThetaW)*Cos(ThetaWp)*Sin(ThetaW)*Sqr(g2);

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

QuadrupleVectorVertex VertexImpl<typename fields::conj<fields::VWm>::type, fields::VWm, fields::VZ, fields::VZp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> part1 = Sin(2*ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW));

   const std::complex<double> part2 = -(Cos(ThetaWp)*Sin(ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)));

   const std::complex<double> part3 = -(Cos(ThetaWp)*Sin(ThetaWp)*Sqr(g2)*Sqr(Cos(ThetaW)));

   return {part1, part2, part3};
}

QuadrupleVectorVertex VertexImpl<typename fields::conj<fields::VWm>::type, fields::VWm, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> part1 = -2*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Cos(ThetaWp));

   const std::complex<double> part2 = Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Cos(ThetaWp));

   const std::complex<double> part3 = Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Cos(ThetaWp));

   return {part1, part2, part3};
}

QuadrupleVectorVertex VertexImpl<typename fields::conj<fields::VWm>::type, fields::VWm, fields::VZp, fields::VZp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> part1 = -2*Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Sin(ThetaWp));

   const std::complex<double> part2 = Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Sin(ThetaWp));

   const std::complex<double> part3 = Sqr(g2)*Sqr(Cos(ThetaW))*Sqr(Sin(ThetaWp));

   return {part1, part2, part3};
}

TripleVectorVertex VertexImpl<typename fields::conj<fields::VWm>::type, fields::VWm, fields::VZp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = g2*Cos(ThetaW)*Sin(ThetaWp);

   return {result, TripleVectorVertex::odd_permutation{}};
}

TripleVectorVertex VertexImpl<typename fields::conj<fields::VWm>::type, fields::VWm, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);
   const auto ThetaWp = DERIVEDPARAMETER(ThetaWp);

   const std::complex<double> result = -(g2*Cos(ThetaW)*Cos(ThetaWp));

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
} // namespace UMSSM_cxx_diagrams
} // namespace flexiblesusy
