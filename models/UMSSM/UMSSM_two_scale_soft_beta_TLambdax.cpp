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

// File generated at Tue 27 Oct 2015 15:14:49

#include "UMSSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TLambdax.
 *
 * @return one-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_TLambdax_one_loop(const Soft_traces& soft_traces) const
{
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Qs = INPUT(Qs);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;


   double beta_TLambdax;

   beta_TLambdax = Re(oneOver16PiSqr*(6*traceAdjYdTYd*Lambdax + 2*
      traceAdjYeTYe*Lambdax + 6*traceAdjYuTYu*Lambdax + 1.2*MassB*Lambdax*Sqr(
      g1) + 6*MassWB*Lambdax*Sqr(g2) + 4*MassU*Lambdax*Sqr(gp)*Sqr(QHd) + 4*
      MassU*Lambdax*Sqr(gp)*Sqr(QHu) + 4*MassU*Lambdax*Sqr(gp)*Sqr(Qs) + (3*
      traceYdAdjYd + traceYeAdjYe + 3*traceYuAdjYu + 12*AbsSqr(Lambdax) - 0.6*
      Sqr(g1) - 3*Sqr(g2) - 2*Sqr(gp)*Sqr(QHd) - 2*Sqr(gp)*Sqr(QHu) - 2*Sqr(gp)
      *Sqr(Qs))*TLambdax));


   return beta_TLambdax;
}

/**
 * Calculates the two-loop beta function of TLambdax.
 *
 * @return two-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_TLambdax_two_loop(const Soft_traces& soft_traces) const
{
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qs = INPUT(Qs);
   const auto Qu = INPUT(Qu);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYuAdjYuTYuAdjYu = TRACE_STRUCT.traceYuAdjYuTYuAdjYu;


   double beta_TLambdax;

   beta_TLambdax = Re(twoLoop*(-0.08*Lambdax*(207*Power(g1,4)*MassB + 375
      *Power(g2,4)*MassWB + 400*Power(gp,4)*MassU*Power(QHd,4) + 400*Power(gp,4
      )*MassU*Power(QHu,4) + 300*Power(gp,4)*MassU*Power(Qs,4) + 450*
      traceYdAdjYdTYdAdjYd + 150*traceYdAdjYuTYuAdjYd + 150*
      traceYeAdjYeTYeAdjYe + 150*traceYuAdjYdTYdAdjYu + 450*
      traceYuAdjYuTYuAdjYu + 10*traceAdjYdTYd*Sqr(g1) - 30*traceAdjYeTYe*Sqr(g1
      ) - 20*traceAdjYuTYu*Sqr(g1) + 20*MassB*traceYuAdjYu*Sqr(g1) + 45*MassB*
      Sqr(g1)*Sqr(g2) + 45*MassWB*Sqr(g1)*Sqr(g2) - 400*traceAdjYdTYd*Sqr(g3) -
      400*traceAdjYuTYu*Sqr(g3) + 400*MassG*traceYuAdjYu*Sqr(g3) - 90*MassB*Qd
      *QHd*Sqr(g1)*Sqr(gp) - 90*MassU*Qd*QHd*Sqr(g1)*Sqr(gp) - 90*MassB*Qe*QHd*
      Sqr(g1)*Sqr(gp) - 90*MassU*Qe*QHd*Sqr(g1)*Sqr(gp) + 90*MassB*Qd*QHu*Sqr(
      g1)*Sqr(gp) + 90*MassU*Qd*QHu*Sqr(g1)*Sqr(gp) + 90*MassB*Qe*QHu*Sqr(g1)*
      Sqr(gp) + 90*MassU*Qe*QHu*Sqr(g1)*Sqr(gp) - 60*MassB*QHd*QHu*Sqr(g1)*Sqr(
      gp) - 60*MassU*QHd*QHu*Sqr(g1)*Sqr(gp) + 90*MassB*QHd*Ql*Sqr(g1)*Sqr(gp)
      + 90*MassU*QHd*Ql*Sqr(g1)*Sqr(gp) - 90*MassB*QHu*Ql*Sqr(g1)*Sqr(gp) - 90*
      MassU*QHu*Ql*Sqr(g1)*Sqr(gp) - 90*MassB*QHd*Qq*Sqr(g1)*Sqr(gp) - 90*MassU
      *QHd*Qq*Sqr(g1)*Sqr(gp) + 90*MassB*QHu*Qq*Sqr(g1)*Sqr(gp) + 90*MassU*QHu*
      Qq*Sqr(g1)*Sqr(gp) + 180*MassB*QHd*Qu*Sqr(g1)*Sqr(gp) + 180*MassU*QHd*Qu*
      Sqr(g1)*Sqr(gp) - 180*MassB*QHu*Qu*Sqr(g1)*Sqr(gp) - 180*MassU*QHu*Qu*Sqr
      (g1)*Sqr(gp) - 150*traceAdjYdTYd*Sqr(gp)*Sqr(Qd) - 50*traceAdjYeTYe*Sqr(
      gp)*Sqr(Qe) + 150*traceAdjYdTYd*Sqr(gp)*Sqr(QHd) + 50*traceAdjYeTYe*Sqr(
      gp)*Sqr(QHd) + 60*MassB*Sqr(g1)*Sqr(gp)*Sqr(QHd) + 60*MassU*Sqr(g1)*Sqr(
      gp)*Sqr(QHd) + 150*MassU*Sqr(g2)*Sqr(gp)*Sqr(QHd) + 150*MassWB*Sqr(g2)*
      Sqr(gp)*Sqr(QHd) + 900*Power(gp,4)*MassU*Sqr(Qd)*Sqr(QHd) + 300*Power(gp,
      4)*MassU*Sqr(Qe)*Sqr(QHd) + 150*traceAdjYuTYu*Sqr(gp)*Sqr(QHu) - 150*
      MassU*traceYuAdjYu*Sqr(gp)*Sqr(QHu) + 60*MassB*Sqr(g1)*Sqr(gp)*Sqr(QHu) +
      60*MassU*Sqr(g1)*Sqr(gp)*Sqr(QHu) + 150*MassU*Sqr(g2)*Sqr(gp)*Sqr(QHu) +
      150*MassWB*Sqr(g2)*Sqr(gp)*Sqr(QHu) + 900*Power(gp,4)*MassU*Sqr(Qd)*Sqr(
      QHu) + 300*Power(gp,4)*MassU*Sqr(Qe)*Sqr(QHu) + 400*Power(gp,4)*MassU*Sqr
      (QHd)*Sqr(QHu) - 50*traceAdjYeTYe*Sqr(gp)*Sqr(Ql) + 600*Power(gp,4)*MassU
      *Sqr(QHd)*Sqr(Ql) + 600*Power(gp,4)*MassU*Sqr(QHu)*Sqr(Ql) + 10*
      traceYeAdjYe*(3*MassB*Sqr(g1) + 5*MassU*Sqr(gp)*(Sqr(Qe) - Sqr(QHd) + Sqr
      (Ql))) - 150*traceAdjYdTYd*Sqr(gp)*Sqr(Qq) - 150*traceAdjYuTYu*Sqr(gp)*
      Sqr(Qq) + 150*MassU*traceYuAdjYu*Sqr(gp)*Sqr(Qq) + 1800*Power(gp,4)*MassU
      *Sqr(QHd)*Sqr(Qq) + 1800*Power(gp,4)*MassU*Sqr(QHu)*Sqr(Qq) - 10*
      traceYdAdjYd*(MassB*Sqr(g1) - 5*(8*MassG*Sqr(g3) + 3*MassU*Sqr(gp)*(Sqr(
      Qd) - Sqr(QHd) + Sqr(Qq)))) + 900*Power(gp,4)*MassU*Sqr(Qd)*Sqr(Qs) + 300
      *Power(gp,4)*MassU*Sqr(Qe)*Sqr(Qs) + 300*Power(gp,4)*MassU*Sqr(QHd)*Sqr(
      Qs) + 300*Power(gp,4)*MassU*Sqr(QHu)*Sqr(Qs) + 600*Power(gp,4)*MassU*Sqr(
      Ql)*Sqr(Qs) + 1800*Power(gp,4)*MassU*Sqr(Qq)*Sqr(Qs) - 150*traceAdjYuTYu*
      Sqr(gp)*Sqr(Qu) + 150*MassU*traceYuAdjYu*Sqr(gp)*Sqr(Qu) + 900*Power(gp,4
      )*MassU*Sqr(QHd)*Sqr(Qu) + 900*Power(gp,4)*MassU*Sqr(QHu)*Sqr(Qu) + 900*
      Power(gp,4)*MassU*Sqr(Qs)*Sqr(Qu)) + (4.14*Power(g1,4) + 7.5*Power(g2,4)
      + 8*Power(gp,4)*Power(QHd,4) + 8*Power(gp,4)*Power(QHu,4) + 6*Power(gp,4)
      *Power(Qs,4) - 9*traceYdAdjYdYdAdjYd - 6*traceYdAdjYuYuAdjYd - 3*
      traceYeAdjYeYeAdjYe - 9*traceYuAdjYuYuAdjYu + 0.8*traceYuAdjYu*Sqr(g1) +
      1.8*Sqr(g1)*Sqr(g2) + 16*traceYuAdjYu*Sqr(g3) - 3.6*Qd*QHd*Sqr(g1)*Sqr(gp
      ) - 3.6*Qe*QHd*Sqr(g1)*Sqr(gp) + 3.6*Qd*QHu*Sqr(g1)*Sqr(gp) + 3.6*Qe*QHu*
      Sqr(g1)*Sqr(gp) - 2.4*QHd*QHu*Sqr(g1)*Sqr(gp) + 3.6*QHd*Ql*Sqr(g1)*Sqr(gp
      ) - 3.6*QHu*Ql*Sqr(g1)*Sqr(gp) - 3.6*QHd*Qq*Sqr(g1)*Sqr(gp) + 3.6*QHu*Qq*
      Sqr(g1)*Sqr(gp) + 7.2*QHd*Qu*Sqr(g1)*Sqr(gp) - 7.2*QHu*Qu*Sqr(g1)*Sqr(gp)
      + 2.4*Sqr(g1)*Sqr(gp)*Sqr(QHd) + 6*Sqr(g2)*Sqr(gp)*Sqr(QHd) + 18*Power(
      gp,4)*Sqr(Qd)*Sqr(QHd) + 6*Power(gp,4)*Sqr(Qe)*Sqr(QHd) - 6*traceYuAdjYu*
      Sqr(gp)*Sqr(QHu) + 2.4*Sqr(g1)*Sqr(gp)*Sqr(QHu) + 6*Sqr(g2)*Sqr(gp)*Sqr(
      QHu) + 18*Power(gp,4)*Sqr(Qd)*Sqr(QHu) + 6*Power(gp,4)*Sqr(Qe)*Sqr(QHu) +
      8*Power(gp,4)*Sqr(QHd)*Sqr(QHu) + 12*Power(gp,4)*Sqr(QHd)*Sqr(Ql) + 12*
      Power(gp,4)*Sqr(QHu)*Sqr(Ql) + 0.4*traceYeAdjYe*(3*Sqr(g1) + 5*Sqr(gp)*(
      Sqr(Qe) - Sqr(QHd) + Sqr(Ql))) + 6*traceYuAdjYu*Sqr(gp)*Sqr(Qq) + 36*
      Power(gp,4)*Sqr(QHd)*Sqr(Qq) + 36*Power(gp,4)*Sqr(QHu)*Sqr(Qq) - 0.4*
      traceYdAdjYd*(Sqr(g1) - 5*(8*Sqr(g3) + 3*Sqr(gp)*(Sqr(Qd) - Sqr(QHd) +
      Sqr(Qq)))) + 18*Power(gp,4)*Sqr(Qd)*Sqr(Qs) + 6*Power(gp,4)*Sqr(Qe)*Sqr(
      Qs) + 6*Power(gp,4)*Sqr(QHd)*Sqr(Qs) + 6*Power(gp,4)*Sqr(QHu)*Sqr(Qs) +
      12*Power(gp,4)*Sqr(Ql)*Sqr(Qs) + 36*Power(gp,4)*Sqr(Qq)*Sqr(Qs) + 6*
      traceYuAdjYu*Sqr(gp)*Sqr(Qu) + 18*Power(gp,4)*Sqr(QHd)*Sqr(Qu) + 18*Power
      (gp,4)*Sqr(QHu)*Sqr(Qu) + 18*Power(gp,4)*Sqr(Qs)*Sqr(Qu))*TLambdax - 50*
      Sqr(Conj(Lambdax))*Sqr(Lambdax)*TLambdax - 0.2*AbsSqr(Lambdax)*(2*Lambdax
      *(45*traceAdjYdTYd + 15*traceAdjYeTYe + 45*traceAdjYuTYu + 6*MassB*Sqr(g1
      ) + 30*MassWB*Sqr(g2) + 20*MassU*Sqr(gp)*Sqr(QHd) + 20*MassU*Sqr(gp)*Sqr(
      QHu)) - 3*(-45*traceYdAdjYd - 15*traceYeAdjYe - 45*traceYuAdjYu + 6*Sqr(
      g1) + 30*Sqr(g2) + 20*Sqr(gp)*Sqr(QHd) + 20*Sqr(gp)*Sqr(QHu))*TLambdax)))
      ;


   return beta_TLambdax;
}

/**
 * Calculates the three-loop beta function of TLambdax.
 *
 * @return three-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_TLambdax_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_TLambdax;

   beta_TLambdax = 0;


   return beta_TLambdax;
}

} // namespace flexiblesusy
