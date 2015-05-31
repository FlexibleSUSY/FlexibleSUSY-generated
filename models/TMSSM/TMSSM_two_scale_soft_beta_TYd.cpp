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

// File generated at Sun 31 May 2015 12:23:40

#include "TMSSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TYd.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_TYd_one_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = (oneOver16PiSqr*(3*traceYdAdjYd*TYd + traceYeAdjYe*TYd +
      1.5*AbsSqr(Lambdax)*TYd - 0.4666666666666667*Sqr(g1)*TYd - 3*Sqr(g2)*TYd
      - 5.333333333333333*Sqr(g3)*TYd + Yd*(6*traceAdjYdTYd + 2*traceAdjYeTYe +
      0.9333333333333333*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 10.666666666666666
      *MassG*Sqr(g3) + 3*Conj(Lambdax)*TLambdax) + 4*(Yd*Yd.adjoint()*TYd) + 2*
      (Yd*Yu.adjoint()*TYu) + 5*(TYd*Yd.adjoint()*Yd) + TYd*Yu.adjoint()*Yu))
      .real();


   return beta_TYd;
}

/**
 * Calculates the two-loop beta function of TYd.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_TYd_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = (twoLoop*(3.188888888888889*Power(g1,4)*TYd + 13.5*Power(g2
      ,4)*TYd - 1.7777777777777777*Power(g3,4)*TYd - 9*traceYdAdjYdYdAdjYd*TYd
      - 3*traceYdAdjYuYuAdjYd*TYd - 3*traceYeAdjYeYeAdjYe*TYd - 4.5*
      traceYuAdjYu*AbsSqr(Lambdax)*TYd - 0.4*traceYdAdjYd*Sqr(g1)*TYd + 1.2*
      traceYeAdjYe*Sqr(g1)*TYd + 6*AbsSqr(Lambdax)*Sqr(g2)*TYd + Sqr(g1)*Sqr(g2
      )*TYd + 16*traceYdAdjYd*Sqr(g3)*TYd + 0.8888888888888888*Sqr(g1)*Sqr(g3)*
      TYd + 8*Sqr(g2)*Sqr(g3)*TYd - 3.75*Sqr(Conj(Lambdax))*Sqr(Lambdax)*TYd -
      0.022222222222222223*Yd*(2*(287*Power(g1,4)*MassB - 160*Power(g3,4)*MassG
      + 1215*Power(g2,4)*MassWB + 810*traceYdAdjYdTYdAdjYd + 135*
      traceYdAdjYuTYuAdjYd + 270*traceYeAdjYeTYeAdjYe + 135*
      traceYuAdjYdTYdAdjYu + 18*traceAdjYdTYd*Sqr(g1) - 54*traceAdjYeTYe*Sqr(g1
      ) + 54*MassB*traceYeAdjYe*Sqr(g1) + 45*MassB*Sqr(g1)*Sqr(g2) + 45*MassWB*
      Sqr(g1)*Sqr(g2) - 720*traceAdjYdTYd*Sqr(g3) + 40*MassB*Sqr(g1)*Sqr(g3) +
      40*MassG*Sqr(g1)*Sqr(g3) + 360*MassG*Sqr(g2)*Sqr(g3) + 360*MassWB*Sqr(g2)
      *Sqr(g3) - 18*traceYdAdjYd*(MassB*Sqr(g1) - 40*MassG*Sqr(g3))) + 675*
      Lambdax*Sqr(Conj(Lambdax))*TLambdax + 135*Conj(Lambdax)*(3*traceAdjYuTYu*
      Lambdax + 4*MassWB*Lambdax*Sqr(g2) + 3*traceYuAdjYu*TLambdax - 4*Sqr(g2)*
      TLambdax)) - 0.2*(90*traceAdjYdTYd + 30*traceAdjYeTYe + 8*MassB*Sqr(g1) +
      60*MassWB*Sqr(g2) + 45*Conj(Lambdax)*TLambdax)*(Yd*Yd.adjoint()*Yd) - 12
      *traceYdAdjYd*(Yd*Yd.adjoint()*TYd) - 4*traceYeAdjYe*(Yd*Yd.adjoint()*TYd
      ) - 6*AbsSqr(Lambdax)*(Yd*Yd.adjoint()*TYd) + 1.2*Sqr(g1)*(Yd*Yd.adjoint(
      )*TYd) + 6*Sqr(g2)*(Yd*Yd.adjoint()*TYd) - 6*traceAdjYuTYu*(Yd*Yu.adjoint
      ()*Yu) - 1.6*MassB*Sqr(g1)*(Yd*Yu.adjoint()*Yu) - 3*Conj(Lambdax)*
      TLambdax*(Yd*Yu.adjoint()*Yu) - 6*traceYuAdjYu*(Yd*Yu.adjoint()*TYu) - 3*
      AbsSqr(Lambdax)*(Yd*Yu.adjoint()*TYu) + 1.6*Sqr(g1)*(Yd*Yu.adjoint()*TYu)
      - 15*traceYdAdjYd*(TYd*Yd.adjoint()*Yd) - 5*traceYeAdjYe*(TYd*Yd.adjoint
      ()*Yd) - 7.5*AbsSqr(Lambdax)*(TYd*Yd.adjoint()*Yd) + 1.2*Sqr(g1)*(TYd*
      Yd.adjoint()*Yd) + 12*Sqr(g2)*(TYd*Yd.adjoint()*Yd) - 3*traceYuAdjYu*(TYd
      *Yu.adjoint()*Yu) - 1.5*AbsSqr(Lambdax)*(TYd*Yu.adjoint()*Yu) + 0.8*Sqr(
      g1)*(TYd*Yu.adjoint()*Yu) - 6*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*TYd) - 8*(
      Yd*Yd.adjoint()*TYd*Yd.adjoint()*Yd) - 2*(Yd*Yu.adjoint()*Yu*Yd.adjoint()
      *TYd) - 4*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*TYu) - 4*(Yd*Yu.adjoint()*TYu*
      Yd.adjoint()*Yd) - 4*(Yd*Yu.adjoint()*TYu*Yu.adjoint()*Yu) - 6*(TYd*
      Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 4*(TYd*Yu.adjoint()*Yu*Yd.adjoint()*Yd
      ) - 2*(TYd*Yu.adjoint()*Yu*Yu.adjoint()*Yu))).real();


   return beta_TYd;
}

} // namespace flexiblesusy
