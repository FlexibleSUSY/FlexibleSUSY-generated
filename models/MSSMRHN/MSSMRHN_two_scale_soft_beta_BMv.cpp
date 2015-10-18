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

// File generated at Sun 18 Oct 2015 13:01:25

#include "MSSMRHN_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of BMv.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_BMv_one_loop(const Soft_traces& soft_traces) const
{


   Eigen::Matrix<double,3,3> beta_BMv;

   beta_BMv = (2*oneOver16PiSqr*(2*(Mv*Yv.conjugate()*(TYv).transpose())
      + Yv*Yv.adjoint()*BMv + BMv*Yv.conjugate()*Yv.transpose() + 2*(TYv*
      Yv.adjoint()*Mv))).real();


   return beta_BMv;
}

/**
 * Calculates the two-loop beta function of BMv.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_BMv_two_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   Eigen::Matrix<double,3,3> beta_BMv;

   beta_BMv = (-0.4*twoLoop*(2*(15*traceAdjYuTYu + 5*traceAdjYvTYv + 3*
      MassB*Sqr(g1) + 15*MassWB*Sqr(g2))*(Mv*Yv.conjugate()*Yv.transpose()) + (
      30*traceYuAdjYu + 10*traceYvAdjYv - 6*(Sqr(g1) + 5*Sqr(g2)))*(Mv*
      Yv.conjugate()*(TYv).transpose()) + 30*traceAdjYuTYu*(Yv*Yv.adjoint()*Mv)
      + 10*traceAdjYvTYv*(Yv*Yv.adjoint()*Mv) + 6*MassB*Sqr(g1)*(Yv*Yv.adjoint
      ()*Mv) + 30*MassWB*Sqr(g2)*(Yv*Yv.adjoint()*Mv) + 15*traceYuAdjYu*(Yv*
      Yv.adjoint()*BMv) + 5*traceYvAdjYv*(Yv*Yv.adjoint()*BMv) - 3*Sqr(g1)*(Yv*
      Yv.adjoint()*BMv) - 15*Sqr(g2)*(Yv*Yv.adjoint()*BMv) + 15*traceYuAdjYu*(
      BMv*Yv.conjugate()*Yv.transpose()) + 5*traceYvAdjYv*(BMv*Yv.conjugate()*
      Yv.transpose()) - 3*Sqr(g1)*(BMv*Yv.conjugate()*Yv.transpose()) - 15*Sqr(
      g2)*(BMv*Yv.conjugate()*Yv.transpose()) + 30*traceYuAdjYu*(TYv*Yv.adjoint
      ()*Mv) + 10*traceYvAdjYv*(TYv*Yv.adjoint()*Mv) - 6*Sqr(g1)*(TYv*
      Yv.adjoint()*Mv) - 30*Sqr(g2)*(TYv*Yv.adjoint()*Mv) + 10*(Mv*Yv.conjugate
      ()*Ye.transpose()*Ye.conjugate()*(TYv).transpose()) + 10*(Mv*Yv.conjugate
      ()*Yv.transpose()*Yv.conjugate()*(TYv).transpose()) + 10*(Mv*Yv.conjugate
      ()*(TYe).transpose()*Ye.conjugate()*Yv.transpose()) + 10*(Mv*Yv.conjugate
      ()*(TYv).transpose()*Yv.conjugate()*Yv.transpose()) + 5*(Yv*Ye.adjoint()*
      Ye*Yv.adjoint()*BMv) + 10*(Yv*Ye.adjoint()*TYe*Yv.adjoint()*Mv) + 5*(Yv*
      Yv.adjoint()*Yv*Yv.adjoint()*BMv) + 10*(Yv*Yv.adjoint()*TYv*Yv.adjoint()*
      Mv) + 5*(BMv*Yv.conjugate()*Ye.transpose()*Ye.conjugate()*Yv.transpose())
      + 5*(BMv*Yv.conjugate()*Yv.transpose()*Yv.conjugate()*Yv.transpose()) +
      10*(TYv*Ye.adjoint()*Ye*Yv.adjoint()*Mv) + 10*(TYv*Yv.adjoint()*Yv*
      Yv.adjoint()*Mv))).real();


   return beta_BMv;
}

/**
 * Calculates the three-loop beta function of BMv.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_BMv_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_BMv;

   beta_BMv = ZEROMATRIX(3,3);


   return beta_BMv;
}

} // namespace flexiblesusy
