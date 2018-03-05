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

// File generated at Mon 5 Mar 2018 18:29:53

#include "UMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambdax.
 *
 * @return 1-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_Lambdax_1_loop(const Susy_traces& susy_traces) const
{
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Qs = INPUT(Qs);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   double beta_Lambdax;

   beta_Lambdax = Re(oneOver16PiSqr*(-0.6*Lambdax*Sqr(g1) + Lambdax*(3*
      traceYdAdjYd + traceYeAdjYe + 3*traceYuAdjYu + traceYvAdjYv - 3*Sqr(g2) -
      2*Sqr(gp)*(Sqr(QHd) + Sqr(QHu) + Sqr(Qs))) + 4*Conj(Lambdax)*Sqr(Lambdax
      )));


   return beta_Lambdax;
}

/**
 * Calculates the 2-loop beta function of Lambdax.
 *
 * @return 2-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_Lambdax_2_loop(const Susy_traces& susy_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto QHd = INPUT(QHd);
   const auto Qe = INPUT(Qe);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qs = INPUT(Qs);
   const auto Qu = INPUT(Qu);
   const auto Qv = INPUT(Qv);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYvAdjYvYvAdjYv = TRACE_STRUCT.traceYvAdjYvYvAdjYv;
   const double traceYvAdjYvTpYeconjYe =
      TRACE_STRUCT.traceYvAdjYvTpYeconjYe;


   double beta_Lambdax;

   beta_Lambdax = Re(-0.02*twoLoop*Lambdax*(-207*Quad(g1) - 10*Sqr(g1)*(9
      *Sqr(g2) - 2*(traceYdAdjYd - 3*traceYeAdjYe - 2*traceYuAdjYu + 3*Sqr(gp)*
      (3*Qd*QHd + 3*Qe*QHd - 3*Qd*QHu - 3*Qe*QHu + 2*QHd*QHu - 3*QHd*Ql + 3*QHu
      *Ql + 3*QHd*Qq - 3*QHu*Qq - 6*QHd*Qu + 6*QHu*Qu - 2*Sqr(QHd) - 2*Sqr(QHu)
      ))) - 10*AbsSqr(Lambdax)*(6*Sqr(g1) + 5*(-3*(3*traceYdAdjYd +
      traceYeAdjYe + 3*traceYuAdjYu + traceYvAdjYv) + 6*Sqr(g2) + 4*Sqr(gp)*(
      Sqr(QHd) + Sqr(QHu)))) - 25*(15*Quad(g2) + 12*Sqr(g2)*Sqr(gp)*(Sqr(QHd) +
      Sqr(QHu)) + 2*(-9*traceYdAdjYdYdAdjYd - 6*traceYdAdjYuYuAdjYd - 3*
      traceYeAdjYeYeAdjYe - 9*traceYuAdjYuYuAdjYu - 2*traceYvAdjYvTpYeconjYe -
      3*traceYvAdjYvYvAdjYv + 16*traceYdAdjYd*Sqr(g3) + 16*traceYuAdjYu*Sqr(g3)
      + 2*Sqr(gp)*(3*traceYdAdjYd*Sqr(Qd) + traceYeAdjYe*Sqr(Qe) - (3*
      traceYdAdjYd + traceYeAdjYe)*Sqr(QHd) - 3*traceYuAdjYu*Sqr(QHu) -
      traceYvAdjYv*Sqr(QHu) + traceYeAdjYe*Sqr(Ql) + traceYvAdjYv*Sqr(Ql) + 3*
      traceYdAdjYd*Sqr(Qq) + 3*traceYuAdjYu*Sqr(Qq) + 3*traceYuAdjYu*Sqr(Qu) +
      traceYvAdjYv*Sqr(Qv)) + 2*Quad(gp)*(4*Quad(QHd) + 4*Quad(QHu) + 3*Quad(Qs
      ) + 4*Sqr(QHd)*Sqr(QHu) + 6*Sqr(QHd)*Sqr(Ql) + 6*Sqr(QHu)*Sqr(Ql) + 18*
      Sqr(QHd)*Sqr(Qq) + 18*Sqr(QHu)*Sqr(Qq) + 3*Sqr(QHd)*Sqr(Qs) + 3*Sqr(QHu)*
      Sqr(Qs) + 6*Sqr(Ql)*Sqr(Qs) + 18*Sqr(Qq)*Sqr(Qs) + 9*Sqr(Qd)*(Sqr(QHd) +
      Sqr(QHu) + Sqr(Qs)) + 3*Sqr(Qe)*(Sqr(QHd) + Sqr(QHu) + Sqr(Qs)) + 9*Sqr(
      QHd)*Sqr(Qu) + 9*Sqr(QHu)*Sqr(Qu) + 9*Sqr(Qs)*Sqr(Qu) + 3*Sqr(QHd)*Sqr(Qv
      ) + 3*Sqr(QHu)*Sqr(Qv) + 3*Sqr(Qs)*Sqr(Qv)))) + 500*Sqr(Conj(Lambdax))*
      Sqr(Lambdax)));


   return beta_Lambdax;
}

/**
 * Calculates the 3-loop beta function of Lambdax.
 *
 * @return 3-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_Lambdax_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambdax;

   beta_Lambdax = 0;


   return beta_Lambdax;
}

/**
 * Calculates the 4-loop beta function of Lambdax.
 *
 * @return 4-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_Lambdax_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambdax;

   beta_Lambdax = 0;


   return beta_Lambdax;
}

} // namespace flexiblesusy
