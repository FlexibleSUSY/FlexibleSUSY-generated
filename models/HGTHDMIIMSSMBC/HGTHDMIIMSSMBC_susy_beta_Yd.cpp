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
 * Calculates the 1-loop beta function of Yd.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> HGTHDMIIMSSMBC_susy_parameters::calc_beta_Yd_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (oneOver16PiSqr*(-0.25*Yd*(-12*traceYdAdjYd - 4*traceYeAdjYe + Sqr
      (g1) - 6*Sqr(g1d) - 2*Sqr(g1dp) + 9*Sqr(g2) + 32*Sqr(g3)) + 1.5*(Yd*Yd.
      adjoint()*Yd) + 0.5*(Yd*Yu.adjoint()*Yu))).real();


   return beta_Yd;
}

/**
 * Calculates the 2-loop beta function of Yd.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> HGTHDMIIMSSMBC_susy_parameters::calc_beta_Yd_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (twoLoop*(-0.0008333333333333334*Yd*(-1200*Lambda3*Lambda4 + 8100*
      traceYdAdjYdYdAdjYd + 2700*traceYdAdjYuYuAdjYd + 2700*traceYeAdjYeYeAdjYe
       - 1800*AbsSqr(Lambda5) - 5400*AbsSqr(Lambda6) - 1800*AbsSqr(Lambda7) +
      234*Quad(g1) + 3375*Quad(g1d) + 675*Quad(g1dp) + 4500*Quad(g2) + 113600*
      Quad(g3) - 750*traceYdAdjYd*Sqr(g1) - 2250*traceYeAdjYe*Sqr(g1) - 675*Sqr
      (g1)*Sqr(g1d) - 225*Sqr(g1)*Sqr(g1dp) + 1350*Sqr(g1d)*Sqr(g1dp) - 6750*
      traceYdAdjYd*Sqr(g2) - 2250*traceYeAdjYe*Sqr(g2) + 1620*Sqr(g1)*Sqr(g2) -
      12375*Sqr(g1d)*Sqr(g2) - 1125*Sqr(g1dp)*Sqr(g2) + 1350*Sqr(g1d)*Sqr(g2u)
      + 450*Sqr(g1dp)*Sqr(g2up) - 24000*traceYdAdjYd*Sqr(g3) - 2480*Sqr(g1)*Sqr
      (g3) - 10800*Sqr(g2)*Sqr(g3) - 7200*Sqr(Lambda1) - 1200*Sqr(Lambda3) -
      1200*Sqr(Lambda4)) + 0.0125*(-960*Lambda1 - 540*traceYdAdjYd - 180*
      traceYeAdjYe + 187*Sqr(g1) - 270*Sqr(g1d) - 90*Sqr(g1dp) + 675*Sqr(g2) +
      1280*Sqr(g3))*(Yd*Yd.adjoint()*Yd) + 0.004166666666666667*(-480*Lambda3 +
      480*Lambda4 - 540*traceYuAdjYu - 53*Sqr(g1) + 495*Sqr(g2) - 270*Sqr(g2u)
      - 90*Sqr(g2up) + 1280*Sqr(g3))*(Yd*Yu.adjoint()*Yu) + 1.5*(Yd*Yd.adjoint(
      )*Yd*Yd.adjoint()*Yd) - 0.25*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*Yd) - 0.25*
      (Yd*Yu.adjoint()*Yu*Yu.adjoint()*Yu))).real();


   return beta_Yd;
}

/**
 * Calculates the 3-loop beta function of Yd.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> HGTHDMIIMSSMBC_susy_parameters::calc_beta_Yd_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = ZEROMATRIX(3,3);


   return beta_Yd;
}

/**
 * Calculates the 4-loop beta function of Yd.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> HGTHDMIIMSSMBC_susy_parameters::calc_beta_Yd_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = ZEROMATRIX(3,3);


   return beta_Yd;
}

/**
 * Calculates the 5-loop beta function of Yd.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> HGTHDMIIMSSMBC_susy_parameters::calc_beta_Yd_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = ZEROMATRIX(3,3);


   return beta_Yd;
}

} // namespace flexiblesusy
