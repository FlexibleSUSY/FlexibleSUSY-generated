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

// File generated at Thu 15 Dec 2016 12:42:51

#include "SM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Yd.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> SM_susy_parameters::calc_beta_Yd_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (oneOver16PiSqr*(-0.25*Yd*(Sqr(g1) + 9*Sqr(g2) + 4*(-3*
      traceYdAdjYd - traceYeAdjYe - 3*traceYuAdjYu + 8*Sqr(g3))) + 1.5*(Yd*
      Yd.adjoint()*Yd) - 1.5*(Yd*Yu.adjoint()*Yu))).real();


   return beta_Yd;
}

/**
 * Calculates the two-loop beta function of Yd.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> SM_susy_parameters::calc_beta_Yd_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (twoLoop*(-0.0016666666666666668*Yd*(127*Power(g1,4) + 5*Sqr
      (g1)*(-15*(5*traceYdAdjYd + 15*traceYeAdjYe + 17*traceYuAdjYu) + 162*Sqr(
      g2) - 248*Sqr(g3)) + 75*(46*Power(g2,4) - 3*Sqr(g2)*(5*(3*traceYdAdjYd +
      traceYeAdjYe + 3*traceYuAdjYu) + 24*Sqr(g3)) + 2*(432*Power(g3,4) - 80*(
      traceYdAdjYd + traceYuAdjYu)*Sqr(g3) + 3*(9*traceYdAdjYdYdAdjYd - 2*
      traceYdAdjYuYuAdjYd + 3*traceYeAdjYeYeAdjYe + 9*traceYuAdjYuYuAdjYu - 2*
      Sqr(Lambdax))))) + 0.0125*(187*Sqr(g1) + 675*Sqr(g2) + 20*(-3*(9*
      traceYdAdjYd + 3*traceYeAdjYe + 9*traceYuAdjYu + 8*Lambdax) + 64*Sqr(g3))
      )*(Yd*Yd.adjoint()*Yd) + 0.0125*(-79*Sqr(g1) + 5*(20*(3*traceYdAdjYd +
      traceYeAdjYe + 3*traceYuAdjYu) + 9*Sqr(g2) - 256*Sqr(g3)))*(Yd*Yu.adjoint
      ()*Yu) + 1.5*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - Yd*Yd.adjoint()*Yd*
      Yu.adjoint()*Yu - 0.25*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*Yd) + 2.75*(Yd*
      Yu.adjoint()*Yu*Yu.adjoint()*Yu))).real();


   return beta_Yd;
}

/**
 * Calculates the three-loop beta function of Yd.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> SM_susy_parameters::calc_beta_Yd_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (ZEROMATRIX(3,3)).real();


   return beta_Yd;
}

} // namespace flexiblesusy
