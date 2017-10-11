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

// File generated at Tue 10 Oct 2017 22:13:32

#include "UMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of vS.
 *
 * @return 1-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_vS_1_loop(const Susy_traces& susy_traces) const
{
   const auto Qs = INPUT(Qs);


   double beta_vS;

   beta_vS = Re(2*oneOver16PiSqr*vS*(-AbsSqr(Lambdax) + Sqr(gp)*Sqr(Qs)))
      ;


   return beta_vS;
}

/**
 * Calculates the 2-loop beta function of vS.
 *
 * @return 2-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_vS_2_loop(const Susy_traces& susy_traces) const
{
   const auto Qs = INPUT(Qs);
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qu = INPUT(Qu);
   const auto Qv = INPUT(Qv);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   double beta_vS;

   beta_vS = Re(-0.2*twoLoop*vS*(2*AbsSqr(Lambdax)*(3*Sqr(g1) + 5*(-3*
      traceYdAdjYd - traceYeAdjYe - 3*traceYuAdjYu - traceYvAdjYv + 3*Sqr(g2) +
      2*Sqr(gp)*(Sqr(QHd) + Sqr(QHu)))) + 5*Quad(gp)*Sqr(Qs)*(9*Sqr(Qd) + 3*
      Sqr(Qe) + 2*Sqr(QHd) + 2*Sqr(QHu) + 6*Sqr(Ql) + 18*Sqr(Qq) + 3*Sqr(Qs) +
      9*Sqr(Qu) + 3*Sqr(Qv)) - 20*Sqr(Conj(Lambdax))*Sqr(Lambdax)));


   return beta_vS;
}

/**
 * Calculates the 3-loop beta function of vS.
 *
 * @return 3-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_vS_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vS;

   beta_vS = 0;


   return beta_vS;
}

} // namespace flexiblesusy
