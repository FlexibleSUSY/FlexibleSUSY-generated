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

// File generated at Sun 26 Aug 2018 14:06:34

#include "HGTHDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Ye.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> HGTHDMIIMSSMBC_susy_parameters::calc_beta_Ye_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (oneOver16PiSqr*(0.25*Ye*(12*traceYdAdjYd + 4*traceYeAdjYe - 9*Sqr
      (g1) + 6*Sqr(g1d) + 2*Sqr(g1dp) - 9*Sqr(g2)) + 1.5*(Ye*Ye.adjoint()*Ye)))
      .real();


   return beta_Ye;
}

/**
 * Calculates the 2-loop beta function of Ye.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> HGTHDMIIMSSMBC_susy_parameters::calc_beta_Ye_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (twoLoop*(Ye*(Lambda3*Lambda4 - 6.75*traceYdAdjYdYdAdjYd - 2.25*
      traceYdAdjYuYuAdjYd - 2.25*traceYeAdjYeYeAdjYe + 1.5*AbsSqr(Lambda5) +
      4.5*AbsSqr(Lambda6) + 1.5*AbsSqr(Lambda7) + 7.905*Quad(g1) - 2.8125*Quad(
      g1d) - 0.5625*Quad(g1dp) - 3.75*Quad(g2) + 0.625*traceYdAdjYd*Sqr(g1) +
      1.875*traceYeAdjYe*Sqr(g1) + 0.5625*Sqr(g1)*Sqr(g1d) + 0.1875*Sqr(g1)*Sqr
      (g1dp) - 1.125*Sqr(g1d)*Sqr(g1dp) + 5.625*traceYdAdjYd*Sqr(g2) + 1.875*
      traceYeAdjYe*Sqr(g2) + 1.35*Sqr(g1)*Sqr(g2) + 10.3125*Sqr(g1d)*Sqr(g2) +
      0.9375*Sqr(g1dp)*Sqr(g2) - 1.125*Sqr(g1d)*Sqr(g2u) - 0.375*Sqr(g1dp)*Sqr(
      g2up) + 20*traceYdAdjYd*Sqr(g3) + 6*Sqr(Lambda1) + Sqr(Lambda3) + Sqr(
      Lambda4)) + 0.0375*(129*Sqr(g1) - 5*(64*Lambda1 + 36*traceYdAdjYd + 12*
      traceYeAdjYe + 18*Sqr(g1d) + 6*Sqr(g1dp) - 45*Sqr(g2)))*(Ye*Ye.adjoint()*
      Ye) + 1.5*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*Ye))).real();


   return beta_Ye;
}

/**
 * Calculates the 3-loop beta function of Ye.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> HGTHDMIIMSSMBC_susy_parameters::calc_beta_Ye_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = ZEROMATRIX(3,3);


   return beta_Ye;
}

/**
 * Calculates the 4-loop beta function of Ye.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> HGTHDMIIMSSMBC_susy_parameters::calc_beta_Ye_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = ZEROMATRIX(3,3);


   return beta_Ye;
}

} // namespace flexiblesusy
