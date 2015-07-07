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

// File generated at Tue 7 Jul 2015 12:46:27

#include "UMSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of vS.
 *
 * @return one-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_vS_one_loop(const Susy_traces& susy_traces) const
{
   const auto Qs = INPUT(Qs);


   double beta_vS;

   beta_vS = Re(2*oneOver16PiSqr*vS*(-AbsSqr(Lambdax) + Sqr(gp)*Sqr(Qs)))
      ;


   return beta_vS;
}

/**
 * Calculates the two-loop beta function of vS.
 *
 * @return two-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_vS_two_loop(const Susy_traces& susy_traces) const
{
   const auto Qs = INPUT(Qs);
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qu = INPUT(Qu);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_vS;

   beta_vS = Re(-0.2*twoLoop*vS*(2*AbsSqr(Lambdax)*(-15*traceYdAdjYd - 5*
      traceYeAdjYe - 15*traceYuAdjYu + 3*Sqr(g1) + 15*Sqr(g2) + 10*Sqr(gp)*Sqr(
      QHd) + 10*Sqr(gp)*Sqr(QHu)) + 5*Power(gp,4)*Sqr(Qs)*(9*Sqr(Qd) + 3*Sqr(Qe
      ) + 2*Sqr(QHd) + 2*Sqr(QHu) + 6*Sqr(Ql) + 18*Sqr(Qq) + 3*Sqr(Qs) + 9*Sqr(
      Qu)) - 20*Sqr(Conj(Lambdax))*Sqr(Lambdax)));


   return beta_vS;
}

/**
 * Calculates the three-loop beta function of vS.
 *
 * @return three-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_vS_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vS;

   beta_vS = 0;


   return beta_vS;
}

} // namespace flexiblesusy
