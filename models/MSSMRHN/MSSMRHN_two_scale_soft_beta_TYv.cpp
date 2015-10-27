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

// File generated at Tue 27 Oct 2015 15:25:47

#include "MSSMRHN_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TYv.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_TYv_one_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   Eigen::Matrix<double,3,3> beta_TYv;

   beta_TYv = (oneOver16PiSqr*(Yv*(6*traceAdjYuTYu + 2*traceAdjYvTYv +
      1.2*MassB*Sqr(g1) + 6*MassWB*Sqr(g2)) + 3*traceYuAdjYu*TYv + traceYvAdjYv
      *TYv - 0.6*Sqr(g1)*TYv - 3*Sqr(g2)*TYv + 2*(Yv*Ye.adjoint()*TYe) + 4*(Yv*
      Yv.adjoint()*TYv) + TYv*Ye.adjoint()*Ye + 5*(TYv*Yv.adjoint()*Yv))).real(
      );


   return beta_TYv;
}

/**
 * Calculates the two-loop beta function of TYv.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_TYv_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYvTYvAdjYe = TRACE_STRUCT.traceYeAdjYvTYvAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYuAdjYuTYuAdjYu = TRACE_STRUCT.traceYuAdjYuTYuAdjYu;
   const double traceYvAdjYeTYeAdjYv = TRACE_STRUCT.traceYvAdjYeTYeAdjYv;
   const double traceYvAdjYvTYvAdjYv = TRACE_STRUCT.traceYvAdjYvTYvAdjYv;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYvYvAdjYe = TRACE_STRUCT.traceYeAdjYvYvAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYvAdjYvYvAdjYv = TRACE_STRUCT.traceYvAdjYvYvAdjYv;


   Eigen::Matrix<double,3,3> beta_TYv;

   beta_TYv = (twoLoop*(-0.08*Yv*(207*Power(g1,4)*MassB + 375*Power(g2,4)
      *MassWB + 75*traceYdAdjYuTYuAdjYd + 25*traceYeAdjYvTYvAdjYe + 75*
      traceYuAdjYdTYdAdjYu + 450*traceYuAdjYuTYuAdjYu + 25*traceYvAdjYeTYeAdjYv
      + 150*traceYvAdjYvTYvAdjYv + 45*MassB*Sqr(g1)*Sqr(g2) + 45*MassWB*Sqr(g1
      )*Sqr(g2) - 20*traceAdjYuTYu*(Sqr(g1) + 20*Sqr(g3)) + 20*traceYuAdjYu*(
      MassB*Sqr(g1) + 20*MassG*Sqr(g3))) + 4.14*Power(g1,4)*TYv + 7.5*Power(g2,
      4)*TYv - 3*traceYdAdjYuYuAdjYd*TYv - traceYeAdjYvYvAdjYe*TYv - 9*
      traceYuAdjYuYuAdjYu*TYv - 3*traceYvAdjYvYvAdjYv*TYv + 0.8*traceYuAdjYu*
      Sqr(g1)*TYv + 1.8*Sqr(g1)*Sqr(g2)*TYv + 16*traceYuAdjYu*Sqr(g3)*TYv - 0.4
      *(15*traceAdjYdTYd + 5*traceAdjYeTYe + 6*MassB*Sqr(g1))*(Yv*Ye.adjoint()*
      Ye) - 6*traceYdAdjYd*(Yv*Ye.adjoint()*TYe) - 2*traceYeAdjYe*(Yv*
      Ye.adjoint()*TYe) + 2.4*Sqr(g1)*(Yv*Ye.adjoint()*TYe) - 18*traceAdjYuTYu*
      (Yv*Yv.adjoint()*Yv) - 6*traceAdjYvTYv*(Yv*Yv.adjoint()*Yv) - 2.4*MassB*
      Sqr(g1)*(Yv*Yv.adjoint()*Yv) - 12*MassWB*Sqr(g2)*(Yv*Yv.adjoint()*Yv) -
      12*traceYuAdjYu*(Yv*Yv.adjoint()*TYv) - 4*traceYvAdjYv*(Yv*Yv.adjoint()*
      TYv) + 1.2*Sqr(g1)*(Yv*Yv.adjoint()*TYv) + 6*Sqr(g2)*(Yv*Yv.adjoint()*TYv
      ) - 3*traceYdAdjYd*(TYv*Ye.adjoint()*Ye) - traceYeAdjYe*(TYv*Ye.adjoint()
      *Ye) + 1.2*Sqr(g1)*(TYv*Ye.adjoint()*Ye) - 15*traceYuAdjYu*(TYv*
      Yv.adjoint()*Yv) - 5*traceYvAdjYv*(TYv*Yv.adjoint()*Yv) + 2.4*Sqr(g1)*(
      TYv*Yv.adjoint()*Yv) + 12*Sqr(g2)*(TYv*Yv.adjoint()*Yv) - 4*(Yv*
      Ye.adjoint()*Ye*Ye.adjoint()*TYe) - 2*(Yv*Ye.adjoint()*Ye*Yv.adjoint()*
      TYv) - 4*(Yv*Ye.adjoint()*TYe*Ye.adjoint()*Ye) - 4*(Yv*Ye.adjoint()*TYe*
      Yv.adjoint()*Yv) - 6*(Yv*Yv.adjoint()*Yv*Yv.adjoint()*TYv) - 8*(Yv*
      Yv.adjoint()*TYv*Yv.adjoint()*Yv) - 2*(TYv*Ye.adjoint()*Ye*Ye.adjoint()*
      Ye) - 4*(TYv*Ye.adjoint()*Ye*Yv.adjoint()*Yv) - 6*(TYv*Yv.adjoint()*Yv*
      Yv.adjoint()*Yv))).real();


   return beta_TYv;
}

/**
 * Calculates the three-loop beta function of TYv.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_TYv_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYv;

   beta_TYv = ZEROMATRIX(3,3);


   return beta_TYv;
}

} // namespace flexiblesusy
