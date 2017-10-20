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

// File generated at Fri 20 Oct 2017 09:03:43

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
 * Calculates the 1-loop beta function of TYv.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_TYv_1_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   Eigen::Matrix<double,3,3> beta_TYv;

   beta_TYv = (oneOver16PiSqr*(1.2*MassB*Yv*Sqr(g1) + 2*Yv*(3*
      traceAdjYuTYu + traceAdjYvTYv + 3*MassWB*Sqr(g2)) + (3*traceYuAdjYu +
      traceYvAdjYv - 0.6*Sqr(g1) - 3*Sqr(g2))*TYv + 2*(Yv*Ye.adjoint()*TYe) + 4
      *(Yv*Yv.adjoint()*TYv) + TYv*Ye.adjoint()*Ye + 5*(TYv*Yv.adjoint()*Yv)))
      .real();


   return beta_TYv;
}

/**
 * Calculates the 2-loop beta function of TYv.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_TYv_2_loop(const Soft_traces& soft_traces) const
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

   beta_TYv = (twoLoop*(0.02*(-4*Yv*(207*MassB*Quad(g1) + 5*Sqr(g1)*(-4*
      traceAdjYuTYu + 4*MassB*traceYuAdjYu + 9*(MassB + MassWB)*Sqr(g2)) + 25*(
      3*traceYdAdjYuTYuAdjYd + traceYeAdjYvTYvAdjYe + 3*traceYuAdjYdTYdAdjYu +
      18*traceYuAdjYuTYuAdjYu + traceYvAdjYeTYeAdjYv + 6*traceYvAdjYvTYvAdjYv +
      15*MassWB*Quad(g2) - 16*(traceAdjYuTYu - MassG*traceYuAdjYu)*Sqr(g3))) +
      (207*Quad(g1) + 10*Sqr(g1)*(4*traceYuAdjYu + 9*Sqr(g2)) + 25*(15*Quad(g2
      ) - 2*(3*traceYdAdjYuYuAdjYd + traceYeAdjYvYvAdjYe + 9*
      traceYuAdjYuYuAdjYu + 3*traceYvAdjYvYvAdjYv - 16*traceYuAdjYu*Sqr(g3))))*
      TYv) - 0.4*(5*(3*traceAdjYdTYd + traceAdjYeTYe) + 6*MassB*Sqr(g1))*(Yv*
      Ye.adjoint()*Ye) + (-2*(3*traceYdAdjYd + traceYeAdjYe) + 2.4*Sqr(g1))*(Yv
      *Ye.adjoint()*TYe) - 1.2*(2*MassB*Sqr(g1) + 5*(3*traceAdjYuTYu +
      traceAdjYvTYv + 2*MassWB*Sqr(g2)))*(Yv*Yv.adjoint()*Yv) + (-4*(3*
      traceYuAdjYu + traceYvAdjYv) + 1.2*Sqr(g1) + 6*Sqr(g2))*(Yv*Yv.adjoint()*
      TYv) + (-3*traceYdAdjYd - traceYeAdjYe + 1.2*Sqr(g1))*(TYv*Ye.adjoint()*
      Ye) + (-5*(3*traceYuAdjYu + traceYvAdjYv) + 2.4*Sqr(g1) + 12*Sqr(g2))*(
      TYv*Yv.adjoint()*Yv) - 4*(Yv*Ye.adjoint()*Ye*Ye.adjoint()*TYe) - 2*(Yv*
      Ye.adjoint()*Ye*Yv.adjoint()*TYv) - 4*(Yv*Ye.adjoint()*TYe*Ye.adjoint()*
      Ye) - 4*(Yv*Ye.adjoint()*TYe*Yv.adjoint()*Yv) - 6*(Yv*Yv.adjoint()*Yv*
      Yv.adjoint()*TYv) - 8*(Yv*Yv.adjoint()*TYv*Yv.adjoint()*Yv) - 2*(TYv*
      Ye.adjoint()*Ye*Ye.adjoint()*Ye) - 4*(TYv*Ye.adjoint()*Ye*Yv.adjoint()*Yv
      ) - 6*(TYv*Yv.adjoint()*Yv*Yv.adjoint()*Yv))).real();


   return beta_TYv;
}

/**
 * Calculates the 3-loop beta function of TYv.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_TYv_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYv;

   beta_TYv = ZEROMATRIX(3,3);


   return beta_TYv;
}

} // namespace flexiblesusy
