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

// File generated at Tue 22 Jan 2019 16:35:37

#include "HSSUSY_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Yu.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> HSSUSY_susy_parameters::calc_beta_Yu_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (oneOver16PiSqr*(-0.05*Yu*(-60*traceYdAdjYd - 20*traceYeAdjYe - 60
      *traceYuAdjYu + 17*Sqr(g1) + 45*Sqr(g2) + 160*Sqr(g3)) - 1.5*(Yu*Yd.
      adjoint()*Yd) + 1.5*(Yu*Yu.adjoint()*Yu))).real();


   return beta_Yu;
}

/**
 * Calculates the 2-loop beta function of Yu.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> HSSUSY_susy_parameters::calc_beta_Yu_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (twoLoop*(0.0016666666666666668*Yu*(-4050*traceYdAdjYdYdAdjYd +
      900*traceYdAdjYuYuAdjYd - 1350*traceYeAdjYeYeAdjYe - 4050*
      traceYuAdjYuYuAdjYu + 1187*Quad(g1) - 3450*Quad(g2) - 64800*Quad(g3) +
      375*traceYdAdjYd*Sqr(g1) + 1125*traceYeAdjYe*Sqr(g1) + 1275*traceYuAdjYu*
      Sqr(g1) + 3375*traceYdAdjYd*Sqr(g2) + 1125*traceYeAdjYe*Sqr(g2) + 3375*
      traceYuAdjYu*Sqr(g2) - 270*Sqr(g1)*Sqr(g2) + 12000*traceYdAdjYd*Sqr(g3) +
      12000*traceYuAdjYu*Sqr(g3) + 760*Sqr(g1)*Sqr(g3) + 5400*Sqr(g2)*Sqr(g3) +
      900*Sqr(Lambdax)) + 0.0125*(300*traceYdAdjYd + 100*traceYeAdjYe + 300*
      traceYuAdjYu - 43*Sqr(g1) + 45*Sqr(g2) - 1280*Sqr(g3))*(Yu*Yd.adjoint()*
      Yd) + 0.0125*(-540*traceYdAdjYd - 180*traceYeAdjYe - 540*traceYuAdjYu -
      480*Lambdax + 223*Sqr(g1) + 675*Sqr(g2) + 1280*Sqr(g3))*(Yu*Yu.adjoint()*
      Yu) + 2.75*(Yu*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 0.25*(Yu*Yd.adjoint()*
      Yd*Yu.adjoint()*Yu) - Yu*Yu.adjoint()*Yu*Yd.adjoint()*Yd + 1.5*(Yu*Yu.
      adjoint()*Yu*Yu.adjoint()*Yu))).real();


   return beta_Yu;
}

/**
 * Calculates the 3-loop beta function of Yu.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> HSSUSY_susy_parameters::calc_beta_Yu_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (0.00005*PROJECTOR*threeLoop*(-90000*Cube(Lambdax) + 321980*Power6
      (g1) + 3396580*Power6(g2) - 12387000*Power6(g3) + 1172056*Power6(Yu(2,2))
      - 27225*Lambdax*Quad(g1) - 106875*Lambdax*Quad(g2) + 1980000*Lambdax*Quad
      (Yu(2,2)) - 94860*Quad(g2)*Sqr(g1) - 301920*Quad(g3)*Sqr(g1) - 609250*
      Quad(Yu(2,2))*Sqr(g1) - 88840*Quad(g1)*Sqr(g2) + 1473080*Quad(g3)*Sqr(g2)
      - 1991250*Quad(Yu(2,2))*Sqr(g2) + 29250*Lambdax*Sqr(g1)*Sqr(g2) - 446380*
      Quad(g1)*Sqr(g3) - 421440*Quad(g2)*Sqr(g3) - 3140000*Quad(Yu(2,2))*Sqr(g3
      ) - 321000*Sqr(g1)*Sqr(g2)*Sqr(g3) + 45000*Sqr(g1)*Sqr(Lambdax) + 225000*
      Sqr(g2)*Sqr(Lambdax) - 488440*Quad(g1)*Sqr(Yu(2,2)) + 339800*Quad(g2)*Sqr
      (Yu(2,2)) + 7275280*Quad(g3)*Sqr(Yu(2,2)) - 127000*Lambdax*Sqr(g1)*Sqr(Yu
      (2,2)) - 675000*Lambdax*Sqr(g2)*Sqr(Yu(2,2)) + 696580*Sqr(g1)*Sqr(g2)*Sqr
      (Yu(2,2)) + 160000*Lambdax*Sqr(g3)*Sqr(Yu(2,2)) + 361480*Sqr(g1)*Sqr(g3)*
      Sqr(Yu(2,2)) + 967400*Sqr(g2)*Sqr(g3)*Sqr(Yu(2,2)) + 18750*Sqr(Lambdax)*
      Sqr(Yu(2,2)))*Yu(2,2)).real();


   return beta_Yu;
}

/**
 * Calculates the 4-loop beta function of Yu.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> HSSUSY_susy_parameters::calc_beta_Yu_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (2308.18*PROJECTOR*Power8(g3)*Quad(oneOver16PiSqr)*Yu(2,2)).real()
      ;


   return beta_Yu;
}

/**
 * Calculates the 5-loop beta function of Yu.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> HSSUSY_susy_parameters::calc_beta_Yu_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = ZEROMATRIX(3,3);


   return beta_Yu;
}

} // namespace flexiblesusy
