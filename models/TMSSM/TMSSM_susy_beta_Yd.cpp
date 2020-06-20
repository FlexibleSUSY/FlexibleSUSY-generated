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


#include "TMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Yd.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_susy_parameters::calc_beta_Yd_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (-0.03333333333333333*Yd*(-90*traceYdAdjYd - 30*traceYeAdjYe - 45*
      AbsSqr(Lambdax) + 14*Sqr(g1) + 90*Sqr(g2) + 160*Sqr(g3)) + 3*(Yd*Yd.
      adjoint()*Yd) + Yd*Yu.adjoint()*Yu).real();


   return oneLoop * beta_Yd;
}

/**
 * Calculates the 2-loop beta function of Yd.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_susy_parameters::calc_beta_Yd_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (0.005555555555555556*Yd*(-1620*traceYdAdjYdYdAdjYd - 540*
      traceYdAdjYuYuAdjYd - 540*traceYeAdjYeYeAdjYe - 810*traceYuAdjYu*AbsSqr(
      Lambdax) + 574*Quad(g1) + 2430*Quad(g2) - 320*Quad(g3) - 72*traceYdAdjYd*
      Sqr(g1) + 216*traceYeAdjYe*Sqr(g1) + 1080*AbsSqr(Lambdax)*Sqr(g2) + 180*
      Sqr(g1)*Sqr(g2) + 2880*traceYdAdjYd*Sqr(g3) + 160*Sqr(g1)*Sqr(g3) + 1440*
      Sqr(g2)*Sqr(g3) - 675*Sqr(Conj(Lambdax))*Sqr(Lambdax)) + 0.1*(-90*
      traceYdAdjYd - 30*traceYeAdjYe - 45*AbsSqr(Lambdax) + 8*Sqr(g1) + 60*Sqr(
      g2))*(Yd*Yd.adjoint()*Yd) + 0.1*(-30*traceYuAdjYu - 15*AbsSqr(Lambdax) +
      8*Sqr(g1))*(Yd*Yu.adjoint()*Yu) - 4*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*Yd)
      - 2*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*Yd) - 2*(Yd*Yu.adjoint()*Yu*Yu.
      adjoint()*Yu)).real();


   return twoLoop * beta_Yd;
}

/**
 * Calculates the 3-loop beta function of Yd.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_susy_parameters::calc_beta_Yd_3_loop(const Susy_traces& susy_traces) const
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
Eigen::Matrix<double,3,3> TMSSM_susy_parameters::calc_beta_Yd_4_loop(const Susy_traces& susy_traces) const
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
Eigen::Matrix<double,3,3> TMSSM_susy_parameters::calc_beta_Yd_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = ZEROMATRIX(3,3);


   return fiveLoop * beta_Yd;
}

} // namespace flexiblesusy
