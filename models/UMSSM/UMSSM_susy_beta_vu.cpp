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


#include "UMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of vu.
 *
 * @return 1-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_vu_1_loop(const Susy_traces& susy_traces) const
{
   const auto QHu = INPUT(QHu);
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   double beta_vu;

   beta_vu = Re(0.1*vu*(-30*traceYuAdjYu - 10*traceYvAdjYv - 10*AbsSqr(Lambdax)
      + 3*Sqr(g1) + 15*Sqr(g2) + 20*Sqr(gp)*Sqr(QHu)));


   return oneLoop * beta_vu;
}

/**
 * Calculates the 2-loop beta function of vu.
 *
 * @return 2-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_vu_2_loop(const Susy_traces& susy_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto QHu = INPUT(QHu);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qs = INPUT(Qs);
   const auto Qu = INPUT(Qu);
   const auto Qv = INPUT(Qv);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYvAdjYvYvAdjYv = TRACE_STRUCT.traceYvAdjYvYvAdjYv;
   const double traceYvAdjYvTpYeconjYe = TRACE_STRUCT.traceYvAdjYvTpYeconjYe;


   double beta_vu;

   beta_vu = Re(-0.005*vu*(-600*traceYdAdjYuYuAdjYd - 1800*traceYuAdjYuYuAdjYu
      - 200*traceYvAdjYvTpYeconjYe - 600*traceYvAdjYvYvAdjYv - 600*traceYdAdjYd
      *AbsSqr(Lambdax) - 200*traceYeAdjYe*AbsSqr(Lambdax) + 207*Quad(g1) + 275*
      Quad(g2) + 800*Quad(gp)*Quad(QHu) + 340*traceYuAdjYu*Sqr(g1) + 60*
      traceYvAdjYv*Sqr(g1) + 60*AbsSqr(Lambdax)*Sqr(g1) + 900*traceYuAdjYu*Sqr(
      g2) + 300*traceYvAdjYv*Sqr(g2) + 300*AbsSqr(Lambdax)*Sqr(g2) + 90*Sqr(g1)
      *Sqr(g2) + 3200*traceYuAdjYu*Sqr(g3) + 360*Qd*QHu*Sqr(g1)*Sqr(gp) + 360*
      Qe*QHu*Sqr(g1)*Sqr(gp) - 120*QHd*QHu*Sqr(g1)*Sqr(gp) - 360*QHu*Ql*Sqr(g1)
      *Sqr(gp) + 360*QHu*Qq*Sqr(g1)*Sqr(gp) - 720*QHu*Qu*Sqr(g1)*Sqr(gp) + 400*
      AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHd) + 240*Sqr(g1)*Sqr(gp)*Sqr(QHu) + 600*Sqr
      (g2)*Sqr(gp)*Sqr(QHu) + 1800*Quad(gp)*Sqr(Qd)*Sqr(QHu) + 600*Quad(gp)*Sqr
      (Qe)*Sqr(QHu) + 400*Quad(gp)*Sqr(QHd)*Sqr(QHu) + 400*traceYvAdjYv*Sqr(gp)
      *Sqr(Ql) + 1200*Quad(gp)*Sqr(QHu)*Sqr(Ql) + 1200*traceYuAdjYu*Sqr(gp)*Sqr
      (Qq) + 3600*Quad(gp)*Sqr(QHu)*Sqr(Qq) + 400*AbsSqr(Lambdax)*Sqr(gp)*Sqr(
      Qs) + 200*Quad(gp)*Sqr(QHu)*Sqr(Qs) + 1200*traceYuAdjYu*Sqr(gp)*Sqr(Qu) +
      1800*Quad(gp)*Sqr(QHu)*Sqr(Qu) + 400*traceYvAdjYv*Sqr(gp)*Sqr(Qv) + 600*
      Quad(gp)*Sqr(QHu)*Sqr(Qv) - 600*Sqr(Conj(Lambdax))*Sqr(Lambdax)));


   return twoLoop * beta_vu;
}

/**
 * Calculates the 3-loop beta function of vu.
 *
 * @return 3-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_vu_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vu;

   beta_vu = 0;


   return threeLoop * beta_vu;
}

/**
 * Calculates the 4-loop beta function of vu.
 *
 * @return 4-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_vu_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vu;

   beta_vu = 0;


   return fourLoop * beta_vu;
}

/**
 * Calculates the 5-loop beta function of vu.
 *
 * @return 5-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_vu_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vu;

   beta_vu = 0;


   return fiveLoop * beta_vu;
}

} // namespace flexiblesusy
