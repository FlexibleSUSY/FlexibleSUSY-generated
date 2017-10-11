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

// File generated at Tue 10 Oct 2017 20:53:49

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

   beta_Yu = (oneOver16PiSqr*(-0.05*Yu*(17*Sqr(g1) + 5*(9*Sqr(g2) - 2*(6*
      traceYuAdjYu + 3*Sqr(g2u) + Sqr(g2up) - 16*Sqr(g3)))) + 0.5*(Yu*
      Yd.adjoint()*Yd) + 1.5*(Yu*Yu.adjoint()*Yu))).real();


   return beta_Yu;
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

   beta_Yu = (twoLoop*(0.0008333333333333334*Yu*(2766*Quad(g1) + 5*Sqr(g1
      )*(510*traceYuAdjYu - 108*Sqr(g2) + 135*Sqr(g2u) + 45*Sqr(g2up) + 304*Sqr
      (g3)) - 25*(-48*Lambda3*Lambda4 + 108*traceYdAdjYuYuAdjYd + 324*
      traceYuAdjYuYuAdjYu + 180*Quad(g2) + 135*Quad(g2u) + 27*Quad(g2up) + 4544
      *Quad(g3) + 54*Sqr(g1d)*Sqr(g2u) + 18*Sqr(g1dp)*Sqr(g2up) + 54*Sqr(g2u)*
      Sqr(g2up) - 960*traceYuAdjYu*Sqr(g3) - 9*Sqr(g2)*(30*traceYuAdjYu + 55*
      Sqr(g2u) + 5*Sqr(g2up) + 48*Sqr(g3)) - 288*Sqr(Lambda2) - 48*Sqr(Lambda3)
      - 48*Sqr(Lambda4) - 72*Sqr(Lambda5) - 72*Sqr(Lambda6) - 216*Sqr(Lambda7)
      )) + (-2*Lambda3 + 2*Lambda4 - 2.25*traceYdAdjYd - 0.75*traceYeAdjYe -
      0.17083333333333334*Sqr(g1) - 1.125*Sqr(g1d) - 0.375*Sqr(g1dp) + 2.0625*
      Sqr(g2) + 5.333333333333333*Sqr(g3))*(Yu*Yd.adjoint()*Yd) + 0.0125*(223*
      Sqr(g1) + 5*(135*Sqr(g2) - 2*(96*Lambda2 + 54*traceYuAdjYu + 27*Sqr(g2u)
      + 9*Sqr(g2up) - 128*Sqr(g3))))*(Yu*Yu.adjoint()*Yu) - 0.25*(Yu*Yd.adjoint
      ()*Yd*Yd.adjoint()*Yd) - 0.25*(Yu*Yd.adjoint()*Yd*Yu.adjoint()*Yu) + 1.5*
      (Yu*Yu.adjoint()*Yu*Yu.adjoint()*Yu))).real();


   return beta_Yu;
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


   return beta_Yu;
}

} // namespace flexiblesusy
