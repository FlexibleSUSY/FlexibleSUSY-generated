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


#include "SplitMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Ye.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> SplitMSSM_susy_parameters::calc_beta_Ye_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (-0.25*Ye*(-12*traceYdAdjYd - 4*traceYeAdjYe - 12*traceYuAdjYu + 9
      *Sqr(g1) + 9*Sqr(g2) - 6*Sqr(g2d) - 6*Sqr(g2u) - 2*Sqr(gYd) - 2*Sqr(gYu))
      + 1.5*(Ye*Ye.adjoint()*Ye)).real();


   return oneLoop * beta_Ye;
}

/**
 * Calculates the 2-loop beta function of Ye.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> SplitMSSM_susy_parameters::calc_beta_Ye_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (0.0025*Ye*(-1200*g2d*g2u*gYd*gYu - 2700*traceYdAdjYdYdAdjYd + 600
      *traceYdAdjYuYuAdjYd - 900*traceYeAdjYeYeAdjYe - 2700*traceYuAdjYuYuAdjYu
       + 3006*Quad(g1) - 1700*Quad(g2) - 1125*Quad(g2d) - 1125*Quad(g2u) - 225*
      Quad(gYd) - 225*Quad(gYu) + 250*traceYdAdjYd*Sqr(g1) + 750*traceYeAdjYe*
      Sqr(g1) + 850*traceYuAdjYu*Sqr(g1) + 2250*traceYdAdjYd*Sqr(g2) + 750*
      traceYeAdjYe*Sqr(g2) + 2250*traceYuAdjYu*Sqr(g2) + 540*Sqr(g1)*Sqr(g2) +
      225*Sqr(g1)*Sqr(g2d) + 4125*Sqr(g2)*Sqr(g2d) + 225*Sqr(g1)*Sqr(g2u) +
      4125*Sqr(g2)*Sqr(g2u) - 300*Sqr(g2d)*Sqr(g2u) + 8000*traceYdAdjYd*Sqr(g3)
      + 8000*traceYuAdjYu*Sqr(g3) + 75*Sqr(g1)*Sqr(gYd) + 375*Sqr(g2)*Sqr(gYd)
      - 450*Sqr(g2d)*Sqr(gYd) + 75*Sqr(g1)*Sqr(gYu) + 375*Sqr(g2)*Sqr(gYu) -
      450*Sqr(g2u)*Sqr(gYu) - 500*Sqr(gYd)*Sqr(gYu) + 600*Sqr(Lambdax)) +
      0.0375*(-180*traceYdAdjYd - 60*traceYeAdjYe - 180*traceYuAdjYu - 160*
      Lambdax + 129*Sqr(g1) + 225*Sqr(g2) - 90*Sqr(g2d) - 90*Sqr(g2u) - 30*Sqr(
      gYd) - 30*Sqr(gYu))*(Ye*Ye.adjoint()*Ye) + 1.5*(Ye*Ye.adjoint()*Ye*Ye.
      adjoint()*Ye)).real();


   return twoLoop * beta_Ye;
}

/**
 * Calculates the 3-loop beta function of Ye.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> SplitMSSM_susy_parameters::calc_beta_Ye_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = ZEROMATRIX(3,3);


   return threeLoop * beta_Ye;
}

/**
 * Calculates the 4-loop beta function of Ye.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> SplitMSSM_susy_parameters::calc_beta_Ye_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = ZEROMATRIX(3,3);


   return fourLoop * beta_Ye;
}

/**
 * Calculates the 5-loop beta function of Ye.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> SplitMSSM_susy_parameters::calc_beta_Ye_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = ZEROMATRIX(3,3);


   return fiveLoop * beta_Ye;
}

} // namespace flexiblesusy
