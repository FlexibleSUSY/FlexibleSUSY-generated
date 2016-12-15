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

// File generated at Thu 15 Dec 2016 12:59:55

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

   beta_TYe = (oneOver16PiSqr*(0.4*Ye*(9*MassB*Sqr(g1) + 5*(3*
      traceAdjYdTYd + traceAdjYeTYe + 3*MassWB*Sqr(g2))) + (3*traceYdAdjYd +
      traceYeAdjYe - 1.8*Sqr(g1) - 3*Sqr(g2))*TYe + 4*(Ye*Ye.adjoint()*TYe) + 2
      *(Ye*Yv.adjoint()*TYv) + 5*(TYe*Ye.adjoint()*Ye) + TYe*Yv.adjoint()*Yv))
      .real();


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

   beta_TYe = (twoLoop*(0.1*(-4*Ye*(135*Power(g1,4)*MassB + Sqr(g1)*(2*(
      traceAdjYdTYd - 3*traceAdjYeTYe - MassB*traceYdAdjYd + 3*MassB*
      traceYeAdjYe) + 9*(MassB + MassWB)*Sqr(g2)) + 5*(15*Power(g2,4)*MassWB +
      18*traceYdAdjYdTYdAdjYd + 3*traceYdAdjYuTYuAdjYd + 6*traceYeAdjYeTYeAdjYe
      + traceYeAdjYvTYvAdjYe + 3*traceYuAdjYdTYdAdjYu + traceYvAdjYeTYeAdjYv -
      16*(traceAdjYdTYd - MassG*traceYdAdjYd)*Sqr(g3))) + (135*Power(g1,4) + 2
      *Sqr(g1)*(-2*traceYdAdjYd + 6*traceYeAdjYe + 9*Sqr(g2)) + 5*(15*Power(g2,
      4) - 2*(9*traceYdAdjYdYdAdjYd + 3*traceYdAdjYuYuAdjYd + 3*
      traceYeAdjYeYeAdjYe + traceYeAdjYvYvAdjYe) + 32*traceYdAdjYd*Sqr(g3)))*
      TYe) - 6*(3*traceAdjYdTYd + traceAdjYeTYe + 2*MassWB*Sqr(g2))*(Ye*
      Ye.adjoint()*Ye) + (-4*(3*traceYdAdjYd + traceYeAdjYe) + 1.2*Sqr(g1) + 6*
      Sqr(g2))*(Ye*Ye.adjoint()*TYe) - 2*(3*traceAdjYuTYu + traceAdjYvTYv)*(Ye*
      Yv.adjoint()*Yv) - 2*(3*traceYuAdjYu + traceYvAdjYv)*(Ye*Yv.adjoint()*TYv
      ) + (-5*(3*traceYdAdjYd + traceYeAdjYe) - 1.2*Sqr(g1) + 12*Sqr(g2))*(TYe*
      Ye.adjoint()*Ye) + (-3*traceYuAdjYu - traceYvAdjYv)*(TYe*Yv.adjoint()*Yv)
      - 6*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*TYe) - 8*(Ye*Ye.adjoint()*TYe*
      Ye.adjoint()*Ye) - 2*(Ye*Yv.adjoint()*Yv*Ye.adjoint()*TYe) - 4*(Ye*
      Yv.adjoint()*Yv*Yv.adjoint()*TYv) - 4*(Ye*Yv.adjoint()*TYv*Ye.adjoint()*
      Ye) - 4*(Ye*Yv.adjoint()*TYv*Yv.adjoint()*Yv) - 6*(TYe*Ye.adjoint()*Ye*
      Ye.adjoint()*Ye) - 4*(TYe*Yv.adjoint()*Yv*Ye.adjoint()*Ye) - 2*(TYe*
      Yv.adjoint()*Yv*Yv.adjoint()*Yv))).real();


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
