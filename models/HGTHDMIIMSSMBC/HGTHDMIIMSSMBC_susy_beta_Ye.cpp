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

// File generated at Tue 10 Oct 2017 20:53:50

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

   beta_Ye = (oneOver16PiSqr*(0.25*Ye*(12*traceYdAdjYd + 4*traceYeAdjYe -
      9*Sqr(g1) + 6*Sqr(g1d) + 2*Sqr(g1dp) - 9*Sqr(g2)) + 1.5*(Ye*Ye.adjoint()
      *Ye))).real();


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

   beta_Ye = (twoLoop*(0.0025*Ye*(3162*Quad(g1) + 5*Sqr(g1)*(50*
      traceYdAdjYd + 150*traceYeAdjYe + 45*Sqr(g1d) + 15*Sqr(g1dp) + 108*Sqr(g2
      )) - 25*(45*Quad(g1d) + 9*Quad(g1dp) + 3*Sqr(g1d)*(6*Sqr(g1dp) - 55*Sqr(
      g2) + 6*Sqr(g2u)) + Sqr(g1dp)*(-15*Sqr(g2) + 6*Sqr(g2up)) + 2*(30*Quad(g2
      ) - 15*(3*traceYdAdjYd + traceYeAdjYe)*Sqr(g2) - 2*(4*Lambda3*Lambda4 -
      27*traceYdAdjYdYdAdjYd - 9*traceYdAdjYuYuAdjYd - 9*traceYeAdjYeYeAdjYe +
      80*traceYdAdjYd*Sqr(g3) + 24*Sqr(Lambda1) + 4*Sqr(Lambda3) + 4*Sqr(
      Lambda4) + 6*Sqr(Lambda5) + 18*Sqr(Lambda6) + 6*Sqr(Lambda7))))) + 0.0375
      *(129*Sqr(g1) - 5*(64*Lambda1 + 36*traceYdAdjYd + 12*traceYeAdjYe + 18*
      Sqr(g1d) + 6*Sqr(g1dp) - 45*Sqr(g2)))*(Ye*Ye.adjoint()*Ye) + 1.5*(Ye*
      Ye.adjoint()*Ye*Ye.adjoint()*Ye))).real();


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

} // namespace flexiblesusy
