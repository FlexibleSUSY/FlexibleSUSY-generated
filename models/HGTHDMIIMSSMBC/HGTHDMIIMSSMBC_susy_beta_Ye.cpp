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

// File generated at Tue 22 Jan 2019 16:22:19

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

   beta_Ye = (oneOver16PiSqr*(-0.25*Ye*(-12*traceYdAdjYd - 4*traceYeAdjYe + 9*
      Sqr(g1) - 6*Sqr(g1d) - 2*Sqr(g1dp) + 9*Sqr(g2)) + 1.5*(Ye*Ye.adjoint()*Ye
      ))).real();


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

   beta_Ye = (twoLoop*(0.0025*Ye*(400*Lambda3*Lambda4 - 2700*
      traceYdAdjYdYdAdjYd - 900*traceYdAdjYuYuAdjYd - 900*traceYeAdjYeYeAdjYe +
      600*AbsSqr(Lambda5) + 1800*AbsSqr(Lambda6) + 600*AbsSqr(Lambda7) + 3162*
      Quad(g1) - 1125*Quad(g1d) - 225*Quad(g1dp) - 1500*Quad(g2) + 250*
      traceYdAdjYd*Sqr(g1) + 750*traceYeAdjYe*Sqr(g1) + 225*Sqr(g1)*Sqr(g1d) +
      75*Sqr(g1)*Sqr(g1dp) - 450*Sqr(g1d)*Sqr(g1dp) + 2250*traceYdAdjYd*Sqr(g2)
      + 750*traceYeAdjYe*Sqr(g2) + 540*Sqr(g1)*Sqr(g2) + 4125*Sqr(g1d)*Sqr(g2)
      + 375*Sqr(g1dp)*Sqr(g2) - 450*Sqr(g1d)*Sqr(g2u) - 150*Sqr(g1dp)*Sqr(g2up)
      + 8000*traceYdAdjYd*Sqr(g3) + 2400*Sqr(Lambda1) + 400*Sqr(Lambda3) + 400*
      Sqr(Lambda4)) + 0.0375*(-320*Lambda1 - 180*traceYdAdjYd - 60*traceYeAdjYe
       + 129*Sqr(g1) - 90*Sqr(g1d) - 30*Sqr(g1dp) + 225*Sqr(g2))*(Ye*Ye.adjoint
      ()*Ye) + 1.5*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*Ye))).real();


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

/**
 * Calculates the 5-loop beta function of Ye.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> HGTHDMIIMSSMBC_susy_parameters::calc_beta_Ye_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = ZEROMATRIX(3,3);


   return beta_Ye;
}

} // namespace flexiblesusy
