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

// File generated at Fri 20 Oct 2017 08:58:50

#include "NUTNMSSM_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of TYu.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> NUTNMSSM_soft_parameters::calc_beta_TYu_1_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = (oneOver16PiSqr*((3*traceYuAdjYu + AbsSqr(Lambdax) -
      0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3))*TYu +
      0.13333333333333333*Yu*(45*traceAdjYuTYu + 13*MassB*Sqr(g1) + 45*MassWB*
      Sqr(g2) + 80*MassG*Sqr(g3) + 15*Conj(Lambdax)*TLambdax) + 2*(Yu*
      Yd.adjoint()*TYd) + 4*(Yu*Yu.adjoint()*TYu) + TYu*Yd.adjoint()*Yd + 5*(
      TYu*Yu.adjoint()*Yu))).real();


   return beta_TYu;
}

/**
 * Calculates the 2-loop beta function of TYu.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> NUTNMSSM_soft_parameters::calc_beta_TYu_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYuAdjYuTYuAdjYu = TRACE_STRUCT.traceYuAdjYuTYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = (twoLoop*(-0.008888888888888889*Yu*(2743*MassB*Quad(g1) +
      25*(27*(traceYdAdjYuTYuAdjYd + traceYuAdjYdTYdAdjYu + 6*
      traceYuAdjYuTYuAdjYu) + 135*MassWB*Quad(g2) - 32*MassG*Quad(g3) - 144*(
      traceAdjYuTYu - MassG*traceYuAdjYu)*Sqr(g3) + 72*(MassG + MassWB)*Sqr(g2)
      *Sqr(g3)) + 5*Sqr(g1)*(45*(MassB + MassWB)*Sqr(g2) + 4*(-9*traceAdjYuTYu
      + 9*MassB*traceYuAdjYu + 34*(MassB + MassG)*Sqr(g3)))) + (-3*
      traceYdAdjYuYuAdjYd - 9*traceYuAdjYuYuAdjYu + 6.095555555555555*Quad(g1)
      + 7.5*Quad(g2) - 1.7777777777777777*Quad(g3) + 16*traceYuAdjYu*Sqr(g3) +
      8*Sqr(g2)*Sqr(g3) + Sqr(g1)*(Sqr(g2) + 0.08888888888888889*(9*
      traceYuAdjYu + 34*Sqr(g3))))*TYu - 3*Lambdax*Sqr(Conj(Lambdax))*(Lambdax*
      TYu + 4*Yu*TLambdax) - Conj(Lambdax)*((3*traceYdAdjYd + traceYeAdjYe + 2*
      AbsSqr(Kappa))*Lambdax*TYu + 2*Yu*((3*traceAdjYdTYd + traceAdjYeTYe)*
      Lambdax + (3*traceYdAdjYd + traceYeAdjYe)*TLambdax + 2*Conj(Kappa)*(
      Lambdax*TKappa + Kappa*TLambdax))) - 0.4*(5*(3*traceAdjYdTYd +
      traceAdjYeTYe) + 2*MassB*Sqr(g1) + 5*Conj(Lambdax)*TLambdax)*(Yu*
      Yd.adjoint()*Yd) + (-2*(3*traceYdAdjYd + traceYeAdjYe) - 2*AbsSqr(Lambdax
      ) + 0.8*Sqr(g1))*(Yu*Yd.adjoint()*TYd) - 0.4*(45*traceAdjYuTYu + 2*MassB*
      Sqr(g1) + 30*MassWB*Sqr(g2) + 15*Conj(Lambdax)*TLambdax)*(Yu*Yu.adjoint()
      *Yu) + (-12*traceYuAdjYu - 4*AbsSqr(Lambdax) + 1.2*Sqr(g1) + 6*Sqr(g2))*(
      Yu*Yu.adjoint()*TYu) + (-3*traceYdAdjYd - traceYeAdjYe - AbsSqr(Lambdax)
      + 0.4*Sqr(g1))*(TYu*Yd.adjoint()*Yd) + (-15*traceYuAdjYu - 5*AbsSqr(
      Lambdax) + 12*Sqr(g2))*(TYu*Yu.adjoint()*Yu) - 4*(Yu*Yd.adjoint()*Yd*
      Yd.adjoint()*TYd) - 2*(Yu*Yd.adjoint()*Yd*Yu.adjoint()*TYu) - 4*(Yu*
      Yd.adjoint()*TYd*Yd.adjoint()*Yd) - 4*(Yu*Yd.adjoint()*TYd*Yu.adjoint()*
      Yu) - 6*(Yu*Yu.adjoint()*Yu*Yu.adjoint()*TYu) - 8*(Yu*Yu.adjoint()*TYu*
      Yu.adjoint()*Yu) - 2*(TYu*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 4*(TYu*
      Yd.adjoint()*Yd*Yu.adjoint()*Yu) - 6*(TYu*Yu.adjoint()*Yu*Yu.adjoint()*Yu
      ))).real();


   return beta_TYu;
}

/**
 * Calculates the 3-loop beta function of TYu.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> NUTNMSSM_soft_parameters::calc_beta_TYu_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = ZEROMATRIX(3,3);


   return beta_TYu;
}

} // namespace flexiblesusy
