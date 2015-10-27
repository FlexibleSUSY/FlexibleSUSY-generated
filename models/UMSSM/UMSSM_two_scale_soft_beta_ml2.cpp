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

// File generated at Tue 27 Oct 2015 15:14:56

#include "UMSSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of ml2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_ml2_one_loop(const Soft_traces& soft_traces) const
{
   const auto Ql = INPUT(Ql);
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_ml2;

   beta_ml2 = (oneOver16PiSqr*(2*mHd2*(Ye.adjoint()*Ye) + 2*((TYe)
      .adjoint()*TYe) + ml2*Ye.adjoint()*Ye + 2*(Ye.adjoint()*me2*Ye) +
      Ye.adjoint()*Ye*ml2 - 0.7745966692414834*g1*Tr11*UNITMATRIX(3) + 2*gp*Ql*
      Tr14*UNITMATRIX(3) - 1.2*AbsSqr(MassB)*Sqr(g1)*UNITMATRIX(3) - 6*AbsSqr(
      MassWB)*Sqr(g2)*UNITMATRIX(3) - 8*AbsSqr(MassU)*Sqr(gp)*Sqr(Ql)*
      UNITMATRIX(3))).real();


   return beta_ml2;
}

/**
 * Calculates the two-loop beta function of ml2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_ml2_two_loop(const Soft_traces& soft_traces) const
{
   const auto Ql = INPUT(Ql);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto Qd = INPUT(Qd);
   const auto QHu = INPUT(QHu);
   const auto Qq = INPUT(Qq);
   const auto Qu = INPUT(Qu);
   const auto Qs = INPUT(Qs);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_ml2;

   beta_ml2 = (twoLoop*(-6*traceconjTYdTpTYd*(Ye.adjoint()*Ye) - 2*
      traceconjTYeTpTYe*(Ye.adjoint()*Ye) - 6*tracemd2YdAdjYd*(Ye.adjoint()*Ye)
      - 2*traceme2YeAdjYe*(Ye.adjoint()*Ye) - 2*traceml2AdjYeYe*(Ye.adjoint()*
      Ye) - 6*tracemq2AdjYdYd*(Ye.adjoint()*Ye) - 12*mHd2*traceYdAdjYd*(
      Ye.adjoint()*Ye) - 4*mHd2*traceYeAdjYe*(Ye.adjoint()*Ye) - 4*mHd2*AbsSqr(
      Lambdax)*(Ye.adjoint()*Ye) - 2*mHu2*AbsSqr(Lambdax)*(Ye.adjoint()*Ye) - 2
      *ms2*AbsSqr(Lambdax)*(Ye.adjoint()*Ye) - 2*AbsSqr(TLambdax)*(Ye.adjoint()
      *Ye) + 2.4*mHd2*Sqr(g1)*(Ye.adjoint()*Ye) + 4*mHd2*Sqr(gp)*Sqr(Qe)*(
      Ye.adjoint()*Ye) + 4*mHd2*Sqr(gp)*Sqr(QHd)*(Ye.adjoint()*Ye) - 4*mHd2*Sqr
      (gp)*Sqr(Ql)*(Ye.adjoint()*Ye) - 6*traceconjTYdTpYd*(Ye.adjoint()*TYe) -
      2*traceconjTYeTpYe*(Ye.adjoint()*TYe) - 2*Conj(TLambdax)*Lambdax*(
      Ye.adjoint()*TYe) - 6*traceAdjYdTYd*((TYe).adjoint()*Ye) - 2*
      traceAdjYeTYe*((TYe).adjoint()*Ye) - 2.4*MassB*Sqr(g1)*((TYe).adjoint()*
      Ye) - 4*MassU*Sqr(gp)*Sqr(Qe)*((TYe).adjoint()*Ye) - 4*MassU*Sqr(gp)*Sqr(
      QHd)*((TYe).adjoint()*Ye) + 4*MassU*Sqr(gp)*Sqr(Ql)*((TYe).adjoint()*Ye)
      - 2*Conj(Lambdax)*TLambdax*((TYe).adjoint()*Ye) - 6*traceYdAdjYd*((TYe)
      .adjoint()*TYe) - 2*traceYeAdjYe*((TYe).adjoint()*TYe) - 2*AbsSqr(Lambdax
      )*((TYe).adjoint()*TYe) + 2.4*Sqr(g1)*((TYe).adjoint()*TYe) + 4*Sqr(gp)*
      Sqr(Qe)*((TYe).adjoint()*TYe) + 4*Sqr(gp)*Sqr(QHd)*((TYe).adjoint()*TYe)
      - 4*Sqr(gp)*Sqr(Ql)*((TYe).adjoint()*TYe) - 3*traceYdAdjYd*(ml2*
      Ye.adjoint()*Ye) - traceYeAdjYe*(ml2*Ye.adjoint()*Ye) - AbsSqr(Lambdax)*(
      ml2*Ye.adjoint()*Ye) + 1.2*Sqr(g1)*(ml2*Ye.adjoint()*Ye) + 2*Sqr(gp)*Sqr(
      Qe)*(ml2*Ye.adjoint()*Ye) + 2*Sqr(gp)*Sqr(QHd)*(ml2*Ye.adjoint()*Ye) - 2*
      Sqr(gp)*Sqr(Ql)*(ml2*Ye.adjoint()*Ye) - 6*traceYdAdjYd*(Ye.adjoint()*me2*
      Ye) - 2*traceYeAdjYe*(Ye.adjoint()*me2*Ye) - 2*AbsSqr(Lambdax)*(
      Ye.adjoint()*me2*Ye) + 2.4*Sqr(g1)*(Ye.adjoint()*me2*Ye) + 4*Sqr(gp)*Sqr(
      Qe)*(Ye.adjoint()*me2*Ye) + 4*Sqr(gp)*Sqr(QHd)*(Ye.adjoint()*me2*Ye) - 4*
      Sqr(gp)*Sqr(Ql)*(Ye.adjoint()*me2*Ye) - 3*traceYdAdjYd*(Ye.adjoint()*Ye*
      ml2) - traceYeAdjYe*(Ye.adjoint()*Ye*ml2) - AbsSqr(Lambdax)*(Ye.adjoint()
      *Ye*ml2) + 1.2*Sqr(g1)*(Ye.adjoint()*Ye*ml2) + 2*Sqr(gp)*Sqr(Qe)*(
      Ye.adjoint()*Ye*ml2) + 2*Sqr(gp)*Sqr(QHd)*(Ye.adjoint()*Ye*ml2) - 2*Sqr(
      gp)*Sqr(Ql)*(Ye.adjoint()*Ye*ml2) - 8*mHd2*(Ye.adjoint()*Ye*Ye.adjoint()*
      Ye) - 4*(Ye.adjoint()*Ye*(TYe).adjoint()*TYe) - 4*(Ye.adjoint()*TYe*(TYe)
      .adjoint()*Ye) - 4*((TYe).adjoint()*Ye*Ye.adjoint()*TYe) - 4*((TYe)
      .adjoint()*TYe*Ye.adjoint()*Ye) - 2*(ml2*Ye.adjoint()*Ye*Ye.adjoint()*Ye)
      - 4*(Ye.adjoint()*me2*Ye*Ye.adjoint()*Ye) - 4*(Ye.adjoint()*Ye*ml2*
      Ye.adjoint()*Ye) - 4*(Ye.adjoint()*Ye*Ye.adjoint()*me2*Ye) - 2*(
      Ye.adjoint()*Ye*Ye.adjoint()*Ye*ml2) + 6*Power(g2,4)*Tr22*UNITMATRIX(3) -
      3.0983866769659336*g1*gp*Ql*Tr2U114*UNITMATRIX(3) - 3.0983866769659336*
      g1*gp*Ql*Tr2U141*UNITMATRIX(3) - 3.0983866769659336*g1*Tr31*UNITMATRIX(3)
      + 8*gp*Ql*Tr34*UNITMATRIX(3) + 33*Power(g2,4)*AbsSqr(MassWB)*UNITMATRIX(
      3) + 1.2*Tr2U111*Sqr(g1)*UNITMATRIX(3) + 3.6*AbsSqr(MassWB)*Sqr(g1)*Sqr(
      g2)*UNITMATRIX(3) + 1.8*MassB*Conj(MassWB)*Sqr(g1)*Sqr(g2)*UNITMATRIX(3)
      + 8*Tr2U144*Sqr(gp)*Sqr(Ql)*UNITMATRIX(3) + 24*AbsSqr(MassWB)*Sqr(g2)*Sqr
      (gp)*Sqr(Ql)*UNITMATRIX(3) + 12*MassU*Conj(MassWB)*Sqr(g2)*Sqr(gp)*Sqr(Ql
      )*UNITMATRIX(3) + 0.12*Conj(MassB)*Sqr(g1)*(40*MassB*(Ye.adjoint()*Ye) -
      20*(Ye.adjoint()*TYe) + (207*MassB*Sqr(g1) + 5*(3*(2*MassB + MassWB)*Sqr(
      g2) + 4*(2*MassB + MassU)*Ql*(-3*Qd - 3*Qe + QHd - QHu + 4*Ql - 3*Qq + 6*
      Qu)*Sqr(gp)))*UNITMATRIX(3)) + 0.8*Conj(MassU)*Sqr(gp)*(5*(Sqr(Qe) + Sqr(
      QHd) - Sqr(Ql))*(2*MassU*(Ye.adjoint()*Ye) - Ye.adjoint()*TYe) + 3*Ql*(-(
      (MassB + 2*MassU)*(3*Qd + 3*Qe - QHd + QHu - 4*Ql + 3*Qq - 6*Qu)*Sqr(g1))
      + 5*Ql*((2*MassU + MassWB)*Sqr(g2) + 2*MassU*Sqr(gp)*(9*Sqr(Qd) + 3*Sqr(
      Qe) + 2*Sqr(QHd) + 2*Sqr(QHu) + 8*Sqr(Ql) + 18*Sqr(Qq) + Sqr(Qs) + 9*Sqr(
      Qu))))*UNITMATRIX(3)))).real();


   return beta_ml2;
}

/**
 * Calculates the three-loop beta function of ml2.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_ml2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_ml2;

   beta_ml2 = ZEROMATRIX(3,3);


   return beta_ml2;
}

} // namespace flexiblesusy
