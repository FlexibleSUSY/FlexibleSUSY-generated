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

// File generated at Fri 20 Oct 2017 08:34:00

#include "MRSSMEFTHiggs_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Yd.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSMEFTHiggs_susy_parameters::calc_beta_Yd_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (oneOver16PiSqr*(Yd*(3*traceYdAdjYd + traceYeAdjYe + AbsSqr(
      LamSD) + 1.5*AbsSqr(LamTD) - 0.4666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3)) + 3*(Yd*Yd.adjoint()*Yd) + Yd*Yu.adjoint()*Yu)
      ).real();


   return beta_Yd;
}

/**
 * Calculates the 2-loop beta function of Yd.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSMEFTHiggs_susy_parameters::calc_beta_Yd_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (twoLoop*(Yd*(-9*traceYdAdjYdYdAdjYd - 3*traceYdAdjYuYuAdjYd
      - 3*traceYeAdjYeYeAdjYe - AbsSqr(LamSD)*(2*AbsSqr(LamSU) + 3*AbsSqr(
      LamTD)) + 3.468888888888889*Quad(g1) + 16.5*Quad(g2) + 14.222222222222221
      *Quad(g3) - 0.4*traceYdAdjYd*Sqr(g1) + 1.2*traceYeAdjYe*Sqr(g1) + Sqr(g1)
      *Sqr(g2) + Conj(LamTD)*(-1.5*LamTD*AbsSqr(LamTU) + 6*LamTD*Sqr(g2)) + 16*
      traceYdAdjYd*Sqr(g3) + 0.8888888888888888*Sqr(g1)*Sqr(g3) + 8*Sqr(g2)*Sqr
      (g3) - 3*Sqr(LamSD)*Sqr(Conj(LamSD)) - 3.75*Sqr(LamTD)*Sqr(Conj(LamTD)))
      + (-9*traceYdAdjYd - 3*traceYeAdjYe - 3*AbsSqr(LamSD) - 4.5*AbsSqr(LamTD)
      + 0.8*Sqr(g1) + 6*Sqr(g2))*(Yd*Yd.adjoint()*Yd) + (-3*traceYuAdjYu -
      AbsSqr(LamSU) - 1.5*AbsSqr(LamTU) + 0.8*Sqr(g1))*(Yd*Yu.adjoint()*Yu) - 4
      *(Yd*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 2*(Yd*Yu.adjoint()*Yu*Yd.adjoint(
      )*Yd) - 2*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*Yu))).real();


   return beta_Yd;
}

/**
 * Calculates the 3-loop beta function of Yd.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSMEFTHiggs_susy_parameters::calc_beta_Yd_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = ZEROMATRIX(3,3);


   return beta_Yd;
}

} // namespace flexiblesusy
