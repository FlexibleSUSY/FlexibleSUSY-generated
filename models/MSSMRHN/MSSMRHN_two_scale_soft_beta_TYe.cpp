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

// File generated at Fri 8 Jan 2016 15:19:51

#include "MSSMRHN_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of TYe.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_TYe_one_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = (oneOver16PiSqr*(Ye*(6*traceAdjYdTYd + 2*traceAdjYeTYe +
      3.6*MassB*Sqr(g1) + 6*MassWB*Sqr(g2)) + 3*traceYdAdjYd*TYe + traceYeAdjYe
      *TYe - 1.8*Sqr(g1)*TYe - 3*Sqr(g2)*TYe + 4*(Ye*Ye.adjoint()*TYe) + 2*(Ye*
      Yv.adjoint()*TYv) + 5*(TYe*Ye.adjoint()*Ye) + TYe*Yv.adjoint()*Yv)).real(
      );


   return beta_TYe;
}

/**
 * Calculates the two-loop beta function of TYe.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_TYe_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYeAdjYvTYvAdjYe = TRACE_STRUCT.traceYeAdjYvTYvAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYvAdjYeTYeAdjYv = TRACE_STRUCT.traceYvAdjYeTYeAdjYv;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYeAdjYvYvAdjYe = TRACE_STRUCT.traceYeAdjYvYvAdjYe;


   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = (twoLoop*(-0.4*Ye*(135*Power(g1,4)*MassB + 75*Power(g2,4)*
      MassWB + 90*traceYdAdjYdTYdAdjYd + 15*traceYdAdjYuTYuAdjYd + 30*
      traceYeAdjYeTYeAdjYe + 5*traceYeAdjYvTYvAdjYe + 15*traceYuAdjYdTYdAdjYu +
      5*traceYvAdjYeTYeAdjYv + 2*traceAdjYdTYd*Sqr(g1) - 6*traceAdjYeTYe*Sqr(
      g1) + 6*MassB*traceYeAdjYe*Sqr(g1) + 9*MassB*Sqr(g1)*Sqr(g2) + 9*MassWB*
      Sqr(g1)*Sqr(g2) - 80*traceAdjYdTYd*Sqr(g3) + traceYdAdjYd*(-2*MassB*Sqr(
      g1) + 80*MassG*Sqr(g3))) + 13.5*Power(g1,4)*TYe + 7.5*Power(g2,4)*TYe - 9
      *traceYdAdjYdYdAdjYd*TYe - 3*traceYdAdjYuYuAdjYd*TYe - 3*
      traceYeAdjYeYeAdjYe*TYe - traceYeAdjYvYvAdjYe*TYe - 0.4*traceYdAdjYd*Sqr(
      g1)*TYe + 1.2*traceYeAdjYe*Sqr(g1)*TYe + 1.8*Sqr(g1)*Sqr(g2)*TYe + 16*
      traceYdAdjYd*Sqr(g3)*TYe - 6*(3*traceAdjYdTYd + traceAdjYeTYe + 2*MassWB*
      Sqr(g2))*(Ye*Ye.adjoint()*Ye) - 12*traceYdAdjYd*(Ye*Ye.adjoint()*TYe) - 4
      *traceYeAdjYe*(Ye*Ye.adjoint()*TYe) + 1.2*Sqr(g1)*(Ye*Ye.adjoint()*TYe) +
      6*Sqr(g2)*(Ye*Ye.adjoint()*TYe) - 6*traceAdjYuTYu*(Ye*Yv.adjoint()*Yv) -
      2*traceAdjYvTYv*(Ye*Yv.adjoint()*Yv) - 6*traceYuAdjYu*(Ye*Yv.adjoint()*
      TYv) - 2*traceYvAdjYv*(Ye*Yv.adjoint()*TYv) - 15*traceYdAdjYd*(TYe*
      Ye.adjoint()*Ye) - 5*traceYeAdjYe*(TYe*Ye.adjoint()*Ye) - 1.2*Sqr(g1)*(
      TYe*Ye.adjoint()*Ye) + 12*Sqr(g2)*(TYe*Ye.adjoint()*Ye) - 3*traceYuAdjYu*
      (TYe*Yv.adjoint()*Yv) - traceYvAdjYv*(TYe*Yv.adjoint()*Yv) - 6*(Ye*
      Ye.adjoint()*Ye*Ye.adjoint()*TYe) - 8*(Ye*Ye.adjoint()*TYe*Ye.adjoint()*
      Ye) - 2*(Ye*Yv.adjoint()*Yv*Ye.adjoint()*TYe) - 4*(Ye*Yv.adjoint()*Yv*
      Yv.adjoint()*TYv) - 4*(Ye*Yv.adjoint()*TYv*Ye.adjoint()*Ye) - 4*(Ye*
      Yv.adjoint()*TYv*Yv.adjoint()*Yv) - 6*(TYe*Ye.adjoint()*Ye*Ye.adjoint()*
      Ye) - 4*(TYe*Yv.adjoint()*Yv*Ye.adjoint()*Ye) - 2*(TYe*Yv.adjoint()*Yv*
      Yv.adjoint()*Yv))).real();


   return beta_TYe;
}

/**
 * Calculates the three-loop beta function of TYe.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_TYe_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = ZEROMATRIX(3,3);


   return beta_TYe;
}

} // namespace flexiblesusy
