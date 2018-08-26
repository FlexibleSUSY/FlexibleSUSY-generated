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

// File generated at Sun 26 Aug 2018 14:56:14

#include "MSSMRHN_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

namespace {

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator+(const Eigen::MatrixBase<Derived>& m, double n)
{
   return m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator+(double n, const Eigen::MatrixBase<Derived>& m)
{
   return m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator-(const Eigen::MatrixBase<Derived>& m, double n)
{
   return m - Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator-(double n, const Eigen::MatrixBase<Derived>& m)
{
   return - m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

} // anonymous namespace

/**
 * Calculates the 1-loop beta function of mv2.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_mv2_1_loop(const Soft_traces& soft_traces) const
{


   Eigen::Matrix<double,3,3> beta_mv2;

   beta_mv2 = (oneOver16PiSqr*(4*mHu2*(Yv*Yv.adjoint()) + 4*(TYv*(TYv).adjoint(
      )) + 2*(mv2*Yv*Yv.adjoint()) + 4*(Yv*ml2*Yv.adjoint()) + 2*(Yv*Yv.adjoint
      ()*mv2))).real();


   return beta_mv2;
}

/**
 * Calculates the 2-loop beta function of mv2.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_mv2_2_loop(const Soft_traces& soft_traces) const
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

   beta_mv2 = (twoLoop*(0.8*(-15*traceconjTYuTpTYu - 5*traceconjTYvTpTYv - 5*
      traceml2AdjYvYv - 15*tracemq2AdjYuYu - 15*tracemu2YuAdjYu - 5*
      tracemv2YvAdjYv - 30*mHu2*traceYuAdjYu - 10*mHu2*traceYvAdjYv + 3*mHu2*
      Sqr(g1) + 6*AbsSqr(MassB)*Sqr(g1) + 15*mHu2*Sqr(g2) + 30*AbsSqr(MassWB)*
      Sqr(g2))*(Yv*Yv.adjoint()) - 0.8*(3*MassB*Sqr(g1) + 5*(3*traceAdjYuTYu +
      traceAdjYvTYv + 3*MassWB*Sqr(g2)))*(Yv*(TYv).adjoint()) - 0.8*(3*Conj(
      MassB)*Sqr(g1) + 5*(3*traceconjTYuTpYu + traceconjTYvTpYv + 3*Conj(MassWB
      )*Sqr(g2)))*(TYv*Yv.adjoint()) + (-12*traceYuAdjYu - 4*traceYvAdjYv + 2.4
      *Sqr(g1) + 12*Sqr(g2))*(TYv*(TYv).adjoint()) + (-2*(3*traceYuAdjYu +
      traceYvAdjYv) + 1.2*Sqr(g1) + 6*Sqr(g2))*(mv2*Yv*Yv.adjoint()) + (-12*
      traceYuAdjYu - 4*traceYvAdjYv + 2.4*Sqr(g1) + 12*Sqr(g2))*(Yv*ml2*Yv.
      adjoint()) + (-2*(3*traceYuAdjYu + traceYvAdjYv) + 1.2*Sqr(g1) + 6*Sqr(g2
      ))*(Yv*Yv.adjoint()*mv2) - 4*(mHd2 + mHu2)*(Yv*Ye.adjoint()*Ye*Yv.adjoint
      ()) - 4*(Yv*Ye.adjoint()*TYe*(TYv).adjoint()) - 8*mHu2*(Yv*Yv.adjoint()*
      Yv*Yv.adjoint()) - 4*(Yv*Yv.adjoint()*TYv*(TYv).adjoint()) - 4*(Yv*(TYe).
      adjoint()*TYe*Yv.adjoint()) - 4*(Yv*(TYv).adjoint()*TYv*Yv.adjoint()) - 4
      *(TYv*Ye.adjoint()*Ye*(TYv).adjoint()) - 4*(TYv*Yv.adjoint()*Yv*(TYv).
      adjoint()) - 4*(TYv*(TYe).adjoint()*Ye*Yv.adjoint()) - 4*(TYv*(TYv).
      adjoint()*Yv*Yv.adjoint()) - 2*(mv2*Yv*Ye.adjoint()*Ye*Yv.adjoint()) - 2*
      (mv2*Yv*Yv.adjoint()*Yv*Yv.adjoint()) - 4*(Yv*ml2*Ye.adjoint()*Ye*Yv.
      adjoint()) - 4*(Yv*ml2*Yv.adjoint()*Yv*Yv.adjoint()) - 4*(Yv*Ye.adjoint()
      *me2*Ye*Yv.adjoint()) - 4*(Yv*Ye.adjoint()*Ye*ml2*Yv.adjoint()) - 2*(Yv*
      Ye.adjoint()*Ye*Yv.adjoint()*mv2) - 4*(Yv*Yv.adjoint()*mv2*Yv*Yv.adjoint(
      )) - 4*(Yv*Yv.adjoint()*Yv*ml2*Yv.adjoint()) - 2*(Yv*Yv.adjoint()*Yv*Yv.
      adjoint()*mv2))).real();


   return beta_mv2;
}

/**
 * Calculates the 3-loop beta function of mv2.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_mv2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mv2;

   beta_mv2 = ZEROMATRIX(3,3);


   return beta_mv2;
}

/**
 * Calculates the 4-loop beta function of mv2.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_mv2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mv2;

   beta_mv2 = ZEROMATRIX(3,3);


   return beta_mv2;
}

} // namespace flexiblesusy
