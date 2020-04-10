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

// File generated at Fri 10 Apr 2020 19:53:08

#include "SplitMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Yu.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> SplitMSSM_susy_parameters::calc_beta_Yu_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (oneOver16PiSqr*(-0.05*Yu*(-60*traceYdAdjYd - 20*traceYeAdjYe - 60
      *traceYuAdjYu + 17*Sqr(g1) + 45*Sqr(g2) - 30*Sqr(g2d) - 30*Sqr(g2u) + 160
      *Sqr(g3) - 10*Sqr(gYd) - 10*Sqr(gYu)) - 1.5*(Yu*Yd.adjoint()*Yd) + 1.5*(
      Yu*Yu.adjoint()*Yu))).real();


   return beta_Yu;
}

/**
 * Calculates the 2-loop beta function of Yu.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> SplitMSSM_susy_parameters::calc_beta_Yu_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (twoLoop*(0.0008333333333333334*Yu*(-3600*g2d*g2u*gYd*gYu - 8100*
      traceYdAdjYdYdAdjYd + 1800*traceYdAdjYuYuAdjYd - 2700*traceYeAdjYeYeAdjYe
       - 8100*traceYuAdjYuYuAdjYu + 2606*Quad(g1) - 5100*Quad(g2) - 3375*Quad(
      g2d) - 3375*Quad(g2u) - 113600*Quad(g3) - 675*Quad(gYd) - 675*Quad(gYu) +
      750*traceYdAdjYd*Sqr(g1) + 2250*traceYeAdjYe*Sqr(g1) + 2550*traceYuAdjYu*
      Sqr(g1) + 6750*traceYdAdjYd*Sqr(g2) + 2250*traceYeAdjYe*Sqr(g2) + 6750*
      traceYuAdjYu*Sqr(g2) - 540*Sqr(g1)*Sqr(g2) + 675*Sqr(g1)*Sqr(g2d) + 12375
      *Sqr(g2)*Sqr(g2d) + 675*Sqr(g1)*Sqr(g2u) + 12375*Sqr(g2)*Sqr(g2u) - 900*
      Sqr(g2d)*Sqr(g2u) + 24000*traceYdAdjYd*Sqr(g3) + 24000*traceYuAdjYu*Sqr(
      g3) + 1520*Sqr(g1)*Sqr(g3) + 10800*Sqr(g2)*Sqr(g3) + 225*Sqr(g1)*Sqr(gYd)
      + 1125*Sqr(g2)*Sqr(gYd) - 1350*Sqr(g2d)*Sqr(gYd) + 225*Sqr(g1)*Sqr(gYu) +
      1125*Sqr(g2)*Sqr(gYu) - 1350*Sqr(g2u)*Sqr(gYu) - 1500*Sqr(gYd)*Sqr(gYu) +
      1800*Sqr(Lambdax)) + 0.0125*(300*traceYdAdjYd + 100*traceYeAdjYe + 300*
      traceYuAdjYu - 43*Sqr(g1) + 45*Sqr(g2) + 150*Sqr(g2d) + 150*Sqr(g2u) -
      1280*Sqr(g3) + 50*Sqr(gYd) + 50*Sqr(gYu))*(Yu*Yd.adjoint()*Yd) + 0.0125*(
      -540*traceYdAdjYd - 180*traceYeAdjYe - 540*traceYuAdjYu - 480*Lambdax +
      223*Sqr(g1) + 675*Sqr(g2) - 270*Sqr(g2d) - 270*Sqr(g2u) + 1280*Sqr(g3) -
      90*Sqr(gYd) - 90*Sqr(gYu))*(Yu*Yu.adjoint()*Yu) + 2.75*(Yu*Yd.adjoint()*
      Yd*Yd.adjoint()*Yd) - 0.25*(Yu*Yd.adjoint()*Yd*Yu.adjoint()*Yu) - Yu*Yu.
      adjoint()*Yu*Yd.adjoint()*Yd + 1.5*(Yu*Yu.adjoint()*Yu*Yu.adjoint()*Yu)))
      .real();


   return beta_Yu;
}

/**
 * Calculates the 3-loop beta function of Yu.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> SplitMSSM_susy_parameters::calc_beta_Yu_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = ZEROMATRIX(3,3);


   return beta_Yu;
}

/**
 * Calculates the 4-loop beta function of Yu.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> SplitMSSM_susy_parameters::calc_beta_Yu_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = ZEROMATRIX(3,3);


   return beta_Yu;
}

/**
 * Calculates the 5-loop beta function of Yu.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> SplitMSSM_susy_parameters::calc_beta_Yu_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = ZEROMATRIX(3,3);


   return beta_Yu;
}

} // namespace flexiblesusy
