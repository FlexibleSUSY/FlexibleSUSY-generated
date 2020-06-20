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


#include "HGTHDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Yu.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> HGTHDMIIMSSMBC_susy_parameters::calc_beta_Yu_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (-0.05*Yu*(-60*traceYuAdjYu + 17*Sqr(g1) + 45*Sqr(g2) - 30*Sqr(g2u
      ) - 10*Sqr(g2up) + 160*Sqr(g3)) + 0.5*(Yu*Yd.adjoint()*Yd) + 1.5*(Yu*Yu.
      adjoint()*Yu)).real();


   return oneLoop * beta_Yu;
}

/**
 * Calculates the 2-loop beta function of Yu.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> HGTHDMIIMSSMBC_susy_parameters::calc_beta_Yu_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (0.0008333333333333334*Yu*(1200*Lambda3*Lambda4 - 2700*
      traceYdAdjYuYuAdjYd - 8100*traceYuAdjYuYuAdjYu + 1800*AbsSqr(Lambda5) +
      1800*AbsSqr(Lambda6) + 5400*AbsSqr(Lambda7) + 2766*Quad(g1) - 4500*Quad(
      g2) - 3375*Quad(g2u) - 675*Quad(g2up) - 113600*Quad(g3) + 2550*
      traceYuAdjYu*Sqr(g1) + 6750*traceYuAdjYu*Sqr(g2) - 540*Sqr(g1)*Sqr(g2) +
      675*Sqr(g1)*Sqr(g2u) - 1350*Sqr(g1d)*Sqr(g2u) + 12375*Sqr(g2)*Sqr(g2u) +
      225*Sqr(g1)*Sqr(g2up) - 450*Sqr(g1dp)*Sqr(g2up) + 1125*Sqr(g2)*Sqr(g2up)
      - 1350*Sqr(g2u)*Sqr(g2up) + 24000*traceYuAdjYu*Sqr(g3) + 1520*Sqr(g1)*Sqr
      (g3) + 10800*Sqr(g2)*Sqr(g3) + 7200*Sqr(Lambda2) + 1200*Sqr(Lambda3) +
      1200*Sqr(Lambda4)) + 0.004166666666666667*(-480*Lambda3 + 480*Lambda4 -
      540*traceYdAdjYd - 180*traceYeAdjYe - 41*Sqr(g1) - 270*Sqr(g1d) - 90*Sqr(
      g1dp) + 495*Sqr(g2) + 1280*Sqr(g3))*(Yu*Yd.adjoint()*Yd) + 0.0125*(-960*
      Lambda2 - 540*traceYuAdjYu + 223*Sqr(g1) + 675*Sqr(g2) - 270*Sqr(g2u) -
      90*Sqr(g2up) + 1280*Sqr(g3))*(Yu*Yu.adjoint()*Yu) - 0.25*(Yu*Yd.adjoint()
      *Yd*Yd.adjoint()*Yd) - 0.25*(Yu*Yd.adjoint()*Yd*Yu.adjoint()*Yu) + 1.5*(
      Yu*Yu.adjoint()*Yu*Yu.adjoint()*Yu)).real();


   return twoLoop * beta_Yu;
}

/**
 * Calculates the 3-loop beta function of Yu.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> HGTHDMIIMSSMBC_susy_parameters::calc_beta_Yu_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = ZEROMATRIX(3,3);


   return threeLoop * beta_Yu;
}

/**
 * Calculates the 4-loop beta function of Yu.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> HGTHDMIIMSSMBC_susy_parameters::calc_beta_Yu_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = ZEROMATRIX(3,3);


   return fourLoop * beta_Yu;
}

/**
 * Calculates the 5-loop beta function of Yu.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> HGTHDMIIMSSMBC_susy_parameters::calc_beta_Yu_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = ZEROMATRIX(3,3);


   return fiveLoop * beta_Yu;
}

} // namespace flexiblesusy
