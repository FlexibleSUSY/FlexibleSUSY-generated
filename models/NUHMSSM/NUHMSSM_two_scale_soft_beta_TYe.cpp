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

// File generated at Tue 8 Sep 2015 13:21:10

#include "NUHMSSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TYe.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> NUHMSSM_soft_parameters::calc_beta_TYe_one_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = (oneOver16PiSqr*(Ye*(6*traceAdjYdTYd + 2*traceAdjYeTYe +
      3.6*MassB*Sqr(g1) + 6*MassWB*Sqr(g2)) + 3*traceYdAdjYd*TYe + traceYeAdjYe
      *TYe - 1.8*Sqr(g1)*TYe - 3*Sqr(g2)*TYe + 4*(Ye*Ye.adjoint()*TYe) + 5*(TYe
      *Ye.adjoint()*Ye))).real();


   return beta_TYe;
}

/**
 * Calculates the two-loop beta function of TYe.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> NUHMSSM_soft_parameters::calc_beta_TYe_two_loop(const Soft_traces& soft_traces) const
{
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


   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = (0.1*twoLoop*(-4*Ye*(135*Power(g1,4)*MassB + 75*Power(g2,4)
      *MassWB + 90*traceYdAdjYdTYdAdjYd + 15*traceYdAdjYuTYuAdjYd + 30*
      traceYeAdjYeTYeAdjYe + 15*traceYuAdjYdTYdAdjYu + 2*traceAdjYdTYd*Sqr(g1)
      - 6*traceAdjYeTYe*Sqr(g1) + 6*MassB*traceYeAdjYe*Sqr(g1) + 9*MassB*Sqr(g1
      )*Sqr(g2) + 9*MassWB*Sqr(g1)*Sqr(g2) - 80*traceAdjYdTYd*Sqr(g3) +
      traceYdAdjYd*(-2*MassB*Sqr(g1) + 80*MassG*Sqr(g3))) + 135*Power(g1,4)*TYe
      + 75*Power(g2,4)*TYe - 90*traceYdAdjYdYdAdjYd*TYe - 30*
      traceYdAdjYuYuAdjYd*TYe - 30*traceYeAdjYeYeAdjYe*TYe - 4*traceYdAdjYd*Sqr
      (g1)*TYe + 12*traceYeAdjYe*Sqr(g1)*TYe + 18*Sqr(g1)*Sqr(g2)*TYe + 160*
      traceYdAdjYd*Sqr(g3)*TYe - 60*(3*traceAdjYdTYd + traceAdjYeTYe + 2*MassWB
      *Sqr(g2))*(Ye*Ye.adjoint()*Ye) - 120*traceYdAdjYd*(Ye*Ye.adjoint()*TYe) -
      40*traceYeAdjYe*(Ye*Ye.adjoint()*TYe) + 12*Sqr(g1)*(Ye*Ye.adjoint()*TYe)
      + 60*Sqr(g2)*(Ye*Ye.adjoint()*TYe) - 150*traceYdAdjYd*(TYe*Ye.adjoint()*
      Ye) - 50*traceYeAdjYe*(TYe*Ye.adjoint()*Ye) - 12*Sqr(g1)*(TYe*Ye.adjoint(
      )*Ye) + 120*Sqr(g2)*(TYe*Ye.adjoint()*Ye) - 60*(Ye*Ye.adjoint()*Ye*
      Ye.adjoint()*TYe) - 80*(Ye*Ye.adjoint()*TYe*Ye.adjoint()*Ye) - 60*(TYe*
      Ye.adjoint()*Ye*Ye.adjoint()*Ye))).real();


   return beta_TYe;
}

/**
 * Calculates the three-loop beta function of TYe.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> NUHMSSM_soft_parameters::calc_beta_TYe_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = ZEROMATRIX(3,3);


   return beta_TYe;
}

} // namespace flexiblesusy
