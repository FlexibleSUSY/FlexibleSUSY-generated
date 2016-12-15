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

// File generated at Thu 15 Dec 2016 12:59:57

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
 * Calculates the one-loop beta function of TYu.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_TYu_one_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = (oneOver16PiSqr*(0.13333333333333333*Yu*(13*MassB*Sqr(g1) +
      5*(9*traceAdjYuTYu + 3*traceAdjYvTYv + 9*MassWB*Sqr(g2) + 16*MassG*Sqr(
      g3))) + (3*traceYuAdjYu + traceYvAdjYv - 0.8666666666666667*Sqr(g1) - 3*
      Sqr(g2) - 5.333333333333333*Sqr(g3))*TYu + 2*(Yu*Yd.adjoint()*TYd) + 4*(
      Yu*Yu.adjoint()*TYu) + TYu*Yd.adjoint()*Yd + 5*(TYu*Yu.adjoint()*Yu)))
      .real();


   return beta_TYu;
}

/**
 * Calculates the two-loop beta function of TYu.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_TYu_two_loop(const Soft_traces& soft_traces) const
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


   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = (twoLoop*(-0.008888888888888889*Yu*(2743*Power(g1,4)*MassB
      + 25*(-32*Power(g3,4)*MassG + 135*Power(g2,4)*MassWB + 9*(3*
      traceYdAdjYuTYuAdjYd + traceYeAdjYvTYvAdjYe + 3*traceYuAdjYdTYdAdjYu + 18
      *traceYuAdjYuTYuAdjYu + traceYvAdjYeTYeAdjYv + 6*traceYvAdjYvTYvAdjYv) -
      144*(traceAdjYuTYu - MassG*traceYuAdjYu)*Sqr(g3) + 72*(MassG + MassWB)*
      Sqr(g2)*Sqr(g3)) + 5*Sqr(g1)*(45*(MassB + MassWB)*Sqr(g2) + 4*(-9*
      traceAdjYuTYu + 9*MassB*traceYuAdjYu + 34*(MassB + MassG)*Sqr(g3)))) + (
      6.095555555555555*Power(g1,4) + 7.5*Power(g2,4) - 1.7777777777777777*
      Power(g3,4) - 3*traceYdAdjYuYuAdjYd - traceYeAdjYvYvAdjYe - 9*
      traceYuAdjYuYuAdjYu - 3*traceYvAdjYvYvAdjYv + 16*traceYuAdjYu*Sqr(g3) + 8
      *Sqr(g2)*Sqr(g3) + Sqr(g1)*(Sqr(g2) + 0.08888888888888889*(9*traceYuAdjYu
      + 34*Sqr(g3))))*TYu + (-2*(3*traceAdjYdTYd + traceAdjYeTYe) - 0.8*MassB*
      Sqr(g1))*(Yu*Yd.adjoint()*Yd) + (-2*(3*traceYdAdjYd + traceYeAdjYe) + 0.8
      *Sqr(g1))*(Yu*Yd.adjoint()*TYd) - 0.4*(2*MassB*Sqr(g1) + 15*(3*
      traceAdjYuTYu + traceAdjYvTYv + 2*MassWB*Sqr(g2)))*(Yu*Yu.adjoint()*Yu) +
      (-4*(3*traceYuAdjYu + traceYvAdjYv) + 1.2*Sqr(g1) + 6*Sqr(g2))*(Yu*
      Yu.adjoint()*TYu) + (-3*traceYdAdjYd - traceYeAdjYe + 0.4*Sqr(g1))*(TYu*
      Yd.adjoint()*Yd) + (-5*(3*traceYuAdjYu + traceYvAdjYv) + 12*Sqr(g2))*(TYu
      *Yu.adjoint()*Yu) - 4*(Yu*Yd.adjoint()*Yd*Yd.adjoint()*TYd) - 2*(Yu*
      Yd.adjoint()*Yd*Yu.adjoint()*TYu) - 4*(Yu*Yd.adjoint()*TYd*Yd.adjoint()*
      Yd) - 4*(Yu*Yd.adjoint()*TYd*Yu.adjoint()*Yu) - 6*(Yu*Yu.adjoint()*Yu*
      Yu.adjoint()*TYu) - 8*(Yu*Yu.adjoint()*TYu*Yu.adjoint()*Yu) - 2*(TYu*
      Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 4*(TYu*Yd.adjoint()*Yd*Yu.adjoint()*Yu
      ) - 6*(TYu*Yu.adjoint()*Yu*Yu.adjoint()*Yu))).real();


   return beta_TYu;
}

/**
 * Calculates the three-loop beta function of TYu.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_TYu_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = ZEROMATRIX(3,3);


   return beta_TYu;
}

} // namespace flexiblesusy
