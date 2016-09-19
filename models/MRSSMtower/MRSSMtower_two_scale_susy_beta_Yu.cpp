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

// File generated at Mon 19 Sep 2016 09:37:05

#include "MRSSMtower_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Yu.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSMtower_susy_parameters::calc_beta_Yu_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (oneOver16PiSqr*(Yu*(3*traceYuAdjYu + AbsSqr(LamSU) + 1.5*
      AbsSqr(LamTU) - 0.8666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3)) + Yu*Yd.adjoint()*Yd + 3*(Yu*Yu.adjoint()*Yu))
      ).real();


   return beta_Yu;
}

/**
 * Calculates the two-loop beta function of Yu.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSMtower_susy_parameters::calc_beta_Yu_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (twoLoop*(Yu*(6.615555555555556*Power(g1,4) + 16.5*Power(g2,
      4) + 14.222222222222221*Power(g3,4) - 3*traceYdAdjYuYuAdjYd - 9*
      traceYuAdjYuYuAdjYu - 2*AbsSqr(LamSD)*AbsSqr(LamSU) - 3*AbsSqr(LamSU)*
      AbsSqr(LamTU) - 1.5*AbsSqr(LamTD)*AbsSqr(LamTU) + 0.8*traceYuAdjYu*Sqr(g1
      ) + 6*AbsSqr(LamTU)*Sqr(g2) + Sqr(g1)*Sqr(g2) + 16*traceYuAdjYu*Sqr(g3) +
      3.022222222222222*Sqr(g1)*Sqr(g3) + 8*Sqr(g2)*Sqr(g3) - 3*Sqr(LamSU)*Sqr
      (Conj(LamSU)) - 3.75*Sqr(LamTU)*Sqr(Conj(LamTU))) + (-3*traceYdAdjYd -
      traceYeAdjYe - AbsSqr(LamSD) - 1.5*AbsSqr(LamTD) + 0.4*Sqr(g1))*(Yu*
      Yd.adjoint()*Yd) + (-9*traceYuAdjYu - 3*AbsSqr(LamSU) - 4.5*AbsSqr(LamTU)
      + 0.4*Sqr(g1) + 6*Sqr(g2))*(Yu*Yu.adjoint()*Yu) - 2*(Yu*Yd.adjoint()*Yd*
      Yd.adjoint()*Yd) - 2*(Yu*Yd.adjoint()*Yd*Yu.adjoint()*Yu) - 4*(Yu*
      Yu.adjoint()*Yu*Yu.adjoint()*Yu))).real();


   return beta_Yu;
}

/**
 * Calculates the three-loop beta function of Yu.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSMtower_susy_parameters::calc_beta_Yu_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = ZEROMATRIX(3,3);


   return beta_Yu;
}

} // namespace flexiblesusy