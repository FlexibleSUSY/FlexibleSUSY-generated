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

// File generated at Tue 10 Oct 2017 21:34:52

#include "TMSSM_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of TYd.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_TYd_1_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = (oneOver16PiSqr*((3*traceYdAdjYd + traceYeAdjYe + 1.5*
      AbsSqr(Lambdax) - 0.4666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3))*TYd + 0.06666666666666667*Yd*(90*traceAdjYdTYd
      + 30*traceAdjYeTYe + 14*MassB*Sqr(g1) + 90*MassWB*Sqr(g2) + 160*MassG*
      Sqr(g3) + 45*Conj(Lambdax)*TLambdax) + 4*(Yd*Yd.adjoint()*TYd) + 2*(Yd*
      Yu.adjoint()*TYu) + 5*(TYd*Yd.adjoint()*Yd) + TYd*Yu.adjoint()*Yu)).real(
      );


   return beta_TYd;
}

/**
 * Calculates the 2-loop beta function of TYd.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_TYd_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
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


   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = (twoLoop*(-0.044444444444444446*Yd*(287*MassB*Quad(g1) + 5*
      (27*(6*traceYdAdjYdTYdAdjYd + traceYdAdjYuTYuAdjYd + 2*
      traceYeAdjYeTYeAdjYe + traceYuAdjYdTYdAdjYu) + 243*MassWB*Quad(g2) - 32*
      MassG*Quad(g3) - 144*(traceAdjYdTYd - MassG*traceYdAdjYd)*Sqr(g3) + 72*(
      MassG + MassWB)*Sqr(g2)*Sqr(g3)) + Sqr(g1)*(45*(MassB + MassWB)*Sqr(g2) +
      2*(9*(traceAdjYdTYd - 3*traceAdjYeTYe - MassB*traceYdAdjYd + 3*MassB*
      traceYeAdjYe) + 20*(MassB + MassG)*Sqr(g3)))) + (-9*traceYdAdjYdYdAdjYd -
      3*traceYdAdjYuYuAdjYd - 3*traceYeAdjYeYeAdjYe + 3.188888888888889*Quad(
      g1) + 13.5*Quad(g2) - 1.7777777777777777*Quad(g3) + Sqr(g1)*(-0.4*
      traceYdAdjYd + 1.2*traceYeAdjYe + Sqr(g2) + 0.8888888888888888*Sqr(g3)) +
      16*traceYdAdjYd*Sqr(g3) + 8*Sqr(g2)*Sqr(g3))*TYd - 3.75*Lambdax*Sqr(Conj
      (Lambdax))*(Lambdax*TYd + 4*Yd*TLambdax) - 1.5*Conj(Lambdax)*((3*
      traceYuAdjYu*Lambdax - 4*Lambdax*Sqr(g2))*TYd + 2*Yd*(3*traceAdjYuTYu*
      Lambdax + 4*MassWB*Lambdax*Sqr(g2) + 3*traceYuAdjYu*TLambdax - 4*Sqr(g2)*
      TLambdax)) + (-0.4*(4*MassB*Sqr(g1) + 15*(3*traceAdjYdTYd + traceAdjYeTYe
      + 2*MassWB*Sqr(g2))) - 9*Conj(Lambdax)*TLambdax)*(Yd*Yd.adjoint()*Yd) +
      (-12*traceYdAdjYd - 4*traceYeAdjYe - 6*AbsSqr(Lambdax) + 1.2*Sqr(g1) + 6*
      Sqr(g2))*(Yd*Yd.adjoint()*TYd) + (-6*traceAdjYuTYu - 1.6*MassB*Sqr(g1) -
      3*Conj(Lambdax)*TLambdax)*(Yd*Yu.adjoint()*Yu) + (-6*traceYuAdjYu - 3*
      AbsSqr(Lambdax) + 1.6*Sqr(g1))*(Yd*Yu.adjoint()*TYu) + (-15*traceYdAdjYd
      - 5*traceYeAdjYe - 7.5*AbsSqr(Lambdax) + 1.2*Sqr(g1) + 12*Sqr(g2))*(TYd*
      Yd.adjoint()*Yd) + (-3*traceYuAdjYu - 1.5*AbsSqr(Lambdax) + 0.8*Sqr(g1))*
      (TYd*Yu.adjoint()*Yu) - 6*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*TYd) - 8*(Yd*
      Yd.adjoint()*TYd*Yd.adjoint()*Yd) - 2*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*
      TYd) - 4*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*TYu) - 4*(Yd*Yu.adjoint()*TYu*
      Yd.adjoint()*Yd) - 4*(Yd*Yu.adjoint()*TYu*Yu.adjoint()*Yu) - 6*(TYd*
      Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 4*(TYd*Yu.adjoint()*Yu*Yd.adjoint()*Yd
      ) - 2*(TYd*Yu.adjoint()*Yu*Yu.adjoint()*Yu))).real();


   return beta_TYd;
}

/**
 * Calculates the 3-loop beta function of TYd.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_TYd_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = ZEROMATRIX(3,3);


   return beta_TYd;
}

} // namespace flexiblesusy
