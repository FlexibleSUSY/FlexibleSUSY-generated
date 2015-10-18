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

// File generated at Sun 18 Oct 2015 13:01:35

#include "MSSMRHN_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of mv2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_mv2_one_loop(const Soft_traces& soft_traces) const
{


   Eigen::Matrix<double,3,3> beta_mv2;

   beta_mv2 = (2*oneOver16PiSqr*(2*mHu2*(Yv*Yv.adjoint()) + 2*(TYv*(TYv)
      .adjoint()) + mv2*Yv*Yv.adjoint() + 2*(Yv*ml2*Yv.adjoint()) + Yv*
      Yv.adjoint()*mv2)).real();


   return beta_mv2;
}

/**
 * Calculates the two-loop beta function of mv2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_mv2_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double traceconjTYvTpTYv = TRACE_STRUCT.traceconjTYvTpTYv;
   const double traceml2AdjYvYv = TRACE_STRUCT.traceml2AdjYvYv;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double tracemv2YvAdjYv = TRACE_STRUCT.tracemv2YvAdjYv;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYvTpYv = TRACE_STRUCT.traceconjTYvTpYv;


   Eigen::Matrix<double,3,3> beta_mv2;

   beta_mv2 = (twoLoop*(0.8*(-15*traceconjTYuTpTYu - 5*traceconjTYvTpTYv
      - 5*traceml2AdjYvYv - 15*tracemq2AdjYuYu - 15*tracemu2YuAdjYu - 5*
      tracemv2YvAdjYv - 30*mHu2*traceYuAdjYu - 10*mHu2*traceYvAdjYv + 3*mHu2*
      Sqr(g1) + 6*AbsSqr(MassB)*Sqr(g1) + 15*mHu2*Sqr(g2) + 30*AbsSqr(MassWB)*
      Sqr(g2))*(Yv*Yv.adjoint()) - 0.4*(2*(15*traceAdjYuTYu + 5*traceAdjYvTYv +
      3*MassB*Sqr(g1) + 15*MassWB*Sqr(g2))*(Yv*(TYv).adjoint()) + 30*
      traceconjTYuTpYu*(TYv*Yv.adjoint()) + 10*traceconjTYvTpYv*(TYv*Yv.adjoint
      ()) + 6*Conj(MassB)*Sqr(g1)*(TYv*Yv.adjoint()) + 30*Conj(MassWB)*Sqr(g2)*
      (TYv*Yv.adjoint()) + 30*traceYuAdjYu*(TYv*(TYv).adjoint()) + 10*
      traceYvAdjYv*(TYv*(TYv).adjoint()) - 6*Sqr(g1)*(TYv*(TYv).adjoint()) - 30
      *Sqr(g2)*(TYv*(TYv).adjoint()) + 15*traceYuAdjYu*(mv2*Yv*Yv.adjoint()) +
      5*traceYvAdjYv*(mv2*Yv*Yv.adjoint()) - 3*Sqr(g1)*(mv2*Yv*Yv.adjoint()) -
      15*Sqr(g2)*(mv2*Yv*Yv.adjoint()) + 30*traceYuAdjYu*(Yv*ml2*Yv.adjoint())
      + 10*traceYvAdjYv*(Yv*ml2*Yv.adjoint()) - 6*Sqr(g1)*(Yv*ml2*Yv.adjoint())
      - 30*Sqr(g2)*(Yv*ml2*Yv.adjoint()) + 15*traceYuAdjYu*(Yv*Yv.adjoint()*
      mv2) + 5*traceYvAdjYv*(Yv*Yv.adjoint()*mv2) - 3*Sqr(g1)*(Yv*Yv.adjoint()*
      mv2) - 15*Sqr(g2)*(Yv*Yv.adjoint()*mv2) + 10*mHd2*(Yv*Ye.adjoint()*Ye*
      Yv.adjoint()) + 10*mHu2*(Yv*Ye.adjoint()*Ye*Yv.adjoint()) + 10*(Yv*
      Ye.adjoint()*TYe*(TYv).adjoint()) + 20*mHu2*(Yv*Yv.adjoint()*Yv*
      Yv.adjoint()) + 10*(Yv*Yv.adjoint()*TYv*(TYv).adjoint()) + 10*(Yv*(TYe)
      .adjoint()*TYe*Yv.adjoint()) + 10*(Yv*(TYv).adjoint()*TYv*Yv.adjoint()) +
      10*(TYv*Ye.adjoint()*Ye*(TYv).adjoint()) + 10*(TYv*Yv.adjoint()*Yv*(TYv)
      .adjoint()) + 10*(TYv*(TYe).adjoint()*Ye*Yv.adjoint()) + 10*(TYv*(TYv)
      .adjoint()*Yv*Yv.adjoint()) + 5*(mv2*Yv*Ye.adjoint()*Ye*Yv.adjoint()) + 5
      *(mv2*Yv*Yv.adjoint()*Yv*Yv.adjoint()) + 10*(Yv*ml2*Ye.adjoint()*Ye*
      Yv.adjoint()) + 10*(Yv*ml2*Yv.adjoint()*Yv*Yv.adjoint()) + 10*(Yv*
      Ye.adjoint()*me2*Ye*Yv.adjoint()) + 10*(Yv*Ye.adjoint()*Ye*ml2*Yv.adjoint
      ()) + 5*(Yv*Ye.adjoint()*Ye*Yv.adjoint()*mv2) + 10*(Yv*Yv.adjoint()*mv2*
      Yv*Yv.adjoint()) + 10*(Yv*Yv.adjoint()*Yv*ml2*Yv.adjoint()) + 5*(Yv*
      Yv.adjoint()*Yv*Yv.adjoint()*mv2)))).real();


   return beta_mv2;
}

/**
 * Calculates the three-loop beta function of mv2.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_mv2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mv2;

   beta_mv2 = ZEROMATRIX(3,3);


   return beta_mv2;
}

} // namespace flexiblesusy
