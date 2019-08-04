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

// File generated at Sun 4 Aug 2019 19:33:40

#include "UMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of gp.
 *
 * @return 1-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_gp_1_loop(const Susy_traces& susy_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qs = INPUT(Qs);
   const auto Qu = INPUT(Qu);
   const auto Qv = INPUT(Qv);


   double beta_gp;

   beta_gp = Re(oneOver16PiSqr*Cube(gp)*(9*Sqr(Qd) + 3*Sqr(Qe) + 2*Sqr(QHd) + 2
      *Sqr(QHu) + 6*Sqr(Ql) + 18*Sqr(Qq) + Sqr(Qs) + 9*Sqr(Qu) + 3*Sqr(Qv)));


   return beta_gp;
}

/**
 * Calculates the 2-loop beta function of gp.
 *
 * @return 2-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_gp_2_loop(const Susy_traces& susy_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
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


   double beta_gp;

   beta_gp = Re(0.4*twoLoop*Cube(gp)*(90*Quad(Qd)*Sqr(gp) + 30*Quad(Qe)*Sqr(gp)
      + 20*Quad(QHd)*Sqr(gp) + 20*Quad(QHu)*Sqr(gp) + 60*Quad(Ql)*Sqr(gp) + 180
      *Quad(Qq)*Sqr(gp) + 10*Quad(Qs)*Sqr(gp) + 90*Quad(Qu)*Sqr(gp) + 30*Quad(
      Qv)*Sqr(gp) - 30*traceYdAdjYd*Sqr(Qd) + 6*Sqr(g1)*Sqr(Qd) + 120*Sqr(g3)*
      Sqr(Qd) - 10*traceYeAdjYe*Sqr(Qe) + 18*Sqr(g1)*Sqr(Qe) - 30*traceYdAdjYd*
      Sqr(QHd) - 10*traceYeAdjYe*Sqr(QHd) - 10*AbsSqr(Lambdax)*Sqr(QHd) + 3*Sqr
      (g1)*Sqr(QHd) + 15*Sqr(g2)*Sqr(QHd) - 30*traceYuAdjYu*Sqr(QHu) - 10*
      traceYvAdjYv*Sqr(QHu) - 10*AbsSqr(Lambdax)*Sqr(QHu) + 3*Sqr(g1)*Sqr(QHu)
      + 15*Sqr(g2)*Sqr(QHu) - 10*traceYeAdjYe*Sqr(Ql) - 10*traceYvAdjYv*Sqr(Ql)
      + 9*Sqr(g1)*Sqr(Ql) + 45*Sqr(g2)*Sqr(Ql) - 30*traceYdAdjYd*Sqr(Qq) - 30*
      traceYuAdjYu*Sqr(Qq) + 3*Sqr(g1)*Sqr(Qq) + 135*Sqr(g2)*Sqr(Qq) + 240*Sqr(
      g3)*Sqr(Qq) - 10*AbsSqr(Lambdax)*Sqr(Qs) - 30*traceYuAdjYu*Sqr(Qu) + 24*
      Sqr(g1)*Sqr(Qu) + 120*Sqr(g3)*Sqr(Qu) - 10*traceYvAdjYv*Sqr(Qv)));


   return beta_gp;
}

/**
 * Calculates the 3-loop beta function of gp.
 *
 * @return 3-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_gp_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_gp;

   beta_gp = 0;


   return beta_gp;
}

/**
 * Calculates the 4-loop beta function of gp.
 *
 * @return 4-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_gp_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_gp;

   beta_gp = 0;


   return beta_gp;
}

/**
 * Calculates the 5-loop beta function of gp.
 *
 * @return 5-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_gp_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_gp;

   beta_gp = 0;


   return beta_gp;
}

} // namespace flexiblesusy
