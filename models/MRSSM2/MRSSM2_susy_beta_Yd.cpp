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


#include "MRSSM2_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Yd.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSM2_susy_parameters::calc_beta_Yd_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (-0.03333333333333333*Yd*(-90*traceYdAdjYd - 30*traceYeAdjYe - 30*
      AbsSqr(LamSD) - 45*AbsSqr(LamTD) + 14*Sqr(g1) + 90*Sqr(g2) + 160*Sqr(g3))
      + 3*(Yd*Yd.adjoint()*Yd) + Yd*Yu.adjoint()*Yu).real();


   return oneLoop * beta_Yd;
}

/**
 * Calculates the 2-loop beta function of Yd.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSM2_susy_parameters::calc_beta_Yd_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (0.0011111111111111111*Yd*(-8100*traceYdAdjYdYdAdjYd - 2700*
      traceYdAdjYuYuAdjYd - 2700*traceYeAdjYeYeAdjYe - 1800*AbsSqr(LamSD)*
      AbsSqr(LamSU) - 2700*AbsSqr(LamSD)*AbsSqr(LamTD) - 1350*AbsSqr(LamTD)*
      AbsSqr(LamTU) + 3122*Quad(g1) + 14850*Quad(g2) + 12800*Quad(g3) - 360*
      traceYdAdjYd*Sqr(g1) + 1080*traceYeAdjYe*Sqr(g1) + 5400*AbsSqr(LamTD)*Sqr
      (g2) + 900*Sqr(g1)*Sqr(g2) + 14400*traceYdAdjYd*Sqr(g3) + 800*Sqr(g1)*Sqr
      (g3) + 7200*Sqr(g2)*Sqr(g3) - 2700*Sqr(LamSD)*Sqr(Conj(LamSD)) - 3375*Sqr
      (LamTD)*Sqr(Conj(LamTD))) + 0.1*(-90*traceYdAdjYd - 30*traceYeAdjYe - 30*
      AbsSqr(LamSD) - 45*AbsSqr(LamTD) + 8*Sqr(g1) + 60*Sqr(g2))*(Yd*Yd.adjoint
      ()*Yd) + 0.1*(-30*traceYuAdjYu - 10*AbsSqr(LamSU) - 15*AbsSqr(LamTU) + 8*
      Sqr(g1))*(Yd*Yu.adjoint()*Yu) - 4*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*Yd) -
      2*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*Yd) - 2*(Yd*Yu.adjoint()*Yu*Yu.adjoint
      ()*Yu)).real();


   return twoLoop * beta_Yd;
}

/**
 * Calculates the 3-loop beta function of Yd.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSM2_susy_parameters::calc_beta_Yd_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = ZEROMATRIX(3,3);


   return threeLoop * beta_Yd;
}

/**
 * Calculates the 4-loop beta function of Yd.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSM2_susy_parameters::calc_beta_Yd_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = ZEROMATRIX(3,3);


   return fourLoop * beta_Yd;
}

/**
 * Calculates the 5-loop beta function of Yd.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSM2_susy_parameters::calc_beta_Yd_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = ZEROMATRIX(3,3);


   return fiveLoop * beta_Yd;
}

} // namespace flexiblesusy
