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

// File generated at Tue 22 Jan 2019 17:28:28

#include "UMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of vd.
 *
 * @return 1-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_vd_1_loop(const Susy_traces& susy_traces) const
{
   const auto QHd = INPUT(QHd);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_vd;

   beta_vd = Re(0.1*oneOver16PiSqr*vd*(-30*traceYdAdjYd - 10*traceYeAdjYe - 10*
      AbsSqr(Lambdax) + 3*Sqr(g1) + 15*Sqr(g2) + 20*Sqr(gp)*Sqr(QHd)));


   return beta_vd;
}

/**
 * Calculates the 2-loop beta function of vd.
 *
 * @return 2-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_vd_2_loop(const Susy_traces& susy_traces) const
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
   const double traceYvAdjYvTpYeconjYe = TRACE_STRUCT.traceYvAdjYvTpYeconjYe;


   double beta_vd;

   beta_vd = Re(-0.005*twoLoop*vd*(-1800*traceYdAdjYdYdAdjYd - 600*
      traceYdAdjYuYuAdjYd - 600*traceYeAdjYeYeAdjYe - 200*
      traceYvAdjYvTpYeconjYe - 600*traceYuAdjYu*AbsSqr(Lambdax) - 200*
      traceYvAdjYv*AbsSqr(Lambdax) + 207*Quad(g1) + 275*Quad(g2) + 800*Quad(gp)
      *Quad(QHd) + 100*traceYdAdjYd*Sqr(g1) + 300*traceYeAdjYe*Sqr(g1) + 60*
      AbsSqr(Lambdax)*Sqr(g1) + 900*traceYdAdjYd*Sqr(g2) + 300*traceYeAdjYe*Sqr
      (g2) + 300*AbsSqr(Lambdax)*Sqr(g2) + 90*Sqr(g1)*Sqr(g2) + 3200*
      traceYdAdjYd*Sqr(g3) - 360*Qd*QHd*Sqr(g1)*Sqr(gp) - 360*Qe*QHd*Sqr(g1)*
      Sqr(gp) - 120*QHd*QHu*Sqr(g1)*Sqr(gp) + 360*QHd*Ql*Sqr(g1)*Sqr(gp) - 360*
      QHd*Qq*Sqr(g1)*Sqr(gp) + 720*QHd*Qu*Sqr(g1)*Sqr(gp) + 1200*traceYdAdjYd*
      Sqr(gp)*Sqr(Qd) + 400*traceYeAdjYe*Sqr(gp)*Sqr(Qe) + 240*Sqr(g1)*Sqr(gp)*
      Sqr(QHd) + 600*Sqr(g2)*Sqr(gp)*Sqr(QHd) + 1800*Quad(gp)*Sqr(Qd)*Sqr(QHd)
      + 600*Quad(gp)*Sqr(Qe)*Sqr(QHd) + 400*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHu) +
      400*Quad(gp)*Sqr(QHd)*Sqr(QHu) + 400*traceYeAdjYe*Sqr(gp)*Sqr(Ql) + 1200*
      Quad(gp)*Sqr(QHd)*Sqr(Ql) + 1200*traceYdAdjYd*Sqr(gp)*Sqr(Qq) + 3600*Quad
      (gp)*Sqr(QHd)*Sqr(Qq) + 400*AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs) + 200*Quad(gp
      )*Sqr(QHd)*Sqr(Qs) + 1800*Quad(gp)*Sqr(QHd)*Sqr(Qu) + 600*Quad(gp)*Sqr(
      QHd)*Sqr(Qv) - 600*Sqr(Conj(Lambdax))*Sqr(Lambdax)));


   return beta_vd;
}

/**
 * Calculates the 3-loop beta function of vd.
 *
 * @return 3-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_vd_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vd;

   beta_vd = 0;


   return beta_vd;
}

/**
 * Calculates the 4-loop beta function of vd.
 *
 * @return 4-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_vd_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vd;

   beta_vd = 0;


   return beta_vd;
}

/**
 * Calculates the 5-loop beta function of vd.
 *
 * @return 5-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_vd_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vd;

   beta_vd = 0;


   return beta_vd;
}

} // namespace flexiblesusy
