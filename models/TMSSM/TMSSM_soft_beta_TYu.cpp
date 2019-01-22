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

// File generated at Tue 22 Jan 2019 16:50:15

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
 * Calculates the 1-loop beta function of TYu.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_TYu_1_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = (oneOver16PiSqr*(0.03333333333333333*(180*traceAdjYuTYu*Yu + 52*
      MassB*Yu*Sqr(g1) + 180*MassWB*Yu*Sqr(g2) + 320*MassG*Yu*Sqr(g3) + 90*
      traceYuAdjYu*TYu + 45*AbsSqr(Lambdax)*TYu - 26*Sqr(g1)*TYu - 90*Sqr(g2)*
      TYu - 160*Sqr(g3)*TYu + 90*Yu*Conj(Lambdax)*TLambdax) + 2*(Yu*Yd.adjoint(
      )*TYd) + 4*(Yu*Yu.adjoint()*TYu) + TYu*Yd.adjoint()*Yd + 5*(TYu*Yu.
      adjoint()*Yu))).real();


   return beta_TYu;
}

/**
 * Calculates the 2-loop beta function of TYu.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_TYu_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYuAdjYuTYuAdjYu = TRACE_STRUCT.traceYuAdjYuTYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = (twoLoop*(0.0011111111111111111*(-5400*traceYdAdjYuTYuAdjYd*Yu -
      5400*traceYuAdjYdTYdAdjYu*Yu - 32400*traceYuAdjYuTYuAdjYu*Yu - 8100*
      traceAdjYdTYd*Yu*AbsSqr(Lambdax) - 2700*traceAdjYeTYe*Yu*AbsSqr(Lambdax)
      - 21944*MassB*Yu*Quad(g1) - 48600*MassWB*Yu*Quad(g2) + 6400*MassG*Yu*Quad
      (g3) + 1440*traceAdjYuTYu*Yu*Sqr(g1) - 1440*MassB*traceYuAdjYu*Yu*Sqr(g1)
      - 10800*MassWB*Yu*AbsSqr(Lambdax)*Sqr(g2) - 1800*MassB*Yu*Sqr(g1)*Sqr(g2)
      - 1800*MassWB*Yu*Sqr(g1)*Sqr(g2) + 28800*traceAdjYuTYu*Yu*Sqr(g3) - 28800
      *MassG*traceYuAdjYu*Yu*Sqr(g3) - 5440*MassB*Yu*Sqr(g1)*Sqr(g3) - 5440*
      MassG*Yu*Sqr(g1)*Sqr(g3) - 14400*MassG*Yu*Sqr(g2)*Sqr(g3) - 14400*MassWB*
      Yu*Sqr(g2)*Sqr(g3) - 2700*traceYdAdjYuYuAdjYd*TYu - 8100*
      traceYuAdjYuYuAdjYu*TYu - 4050*traceYdAdjYd*AbsSqr(Lambdax)*TYu - 1350*
      traceYeAdjYe*AbsSqr(Lambdax)*TYu + 5486*Quad(g1)*TYu + 12150*Quad(g2)*TYu
       - 1600*Quad(g3)*TYu + 720*traceYuAdjYu*Sqr(g1)*TYu + 5400*AbsSqr(Lambdax
      )*Sqr(g2)*TYu + 900*Sqr(g1)*Sqr(g2)*TYu + 14400*traceYuAdjYu*Sqr(g3)*TYu
      + 2720*Sqr(g1)*Sqr(g3)*TYu + 7200*Sqr(g2)*Sqr(g3)*TYu - 3375*Sqr(Conj(
      Lambdax))*Sqr(Lambdax)*TYu - 8100*traceYdAdjYd*Yu*Conj(Lambdax)*TLambdax
      - 2700*traceYeAdjYe*Yu*Conj(Lambdax)*TLambdax + 10800*Yu*Conj(Lambdax)*
      Sqr(g2)*TLambdax - 13500*Yu*Lambdax*Sqr(Conj(Lambdax))*TLambdax) + 0.2*(-
      30*traceAdjYdTYd - 10*traceAdjYeTYe - 4*MassB*Sqr(g1) - 15*Conj(Lambdax)*
      TLambdax)*(Yu*Yd.adjoint()*Yd) + 0.2*(-30*traceYdAdjYd - 10*traceYeAdjYe
      - 15*AbsSqr(Lambdax) + 4*Sqr(g1))*(Yu*Yd.adjoint()*TYd) + 0.2*(-90*
      traceAdjYuTYu - 4*MassB*Sqr(g1) - 60*MassWB*Sqr(g2) - 45*Conj(Lambdax)*
      TLambdax)*(Yu*Yu.adjoint()*Yu) + 1.2*(-10*traceYuAdjYu - 5*AbsSqr(Lambdax
      ) + Sqr(g1) + 5*Sqr(g2))*(Yu*Yu.adjoint()*TYu) + 0.1*(-30*traceYdAdjYd -
      10*traceYeAdjYe - 15*AbsSqr(Lambdax) + 4*Sqr(g1))*(TYu*Yd.adjoint()*Yd) +
      1.5*(-10*traceYuAdjYu - 5*AbsSqr(Lambdax) + 8*Sqr(g2))*(TYu*Yu.adjoint()*
      Yu) - 4*(Yu*Yd.adjoint()*Yd*Yd.adjoint()*TYd) - 2*(Yu*Yd.adjoint()*Yd*Yu.
      adjoint()*TYu) - 4*(Yu*Yd.adjoint()*TYd*Yd.adjoint()*Yd) - 4*(Yu*Yd.
      adjoint()*TYd*Yu.adjoint()*Yu) - 6*(Yu*Yu.adjoint()*Yu*Yu.adjoint()*TYu)
      - 8*(Yu*Yu.adjoint()*TYu*Yu.adjoint()*Yu) - 2*(TYu*Yd.adjoint()*Yd*Yd.
      adjoint()*Yd) - 4*(TYu*Yd.adjoint()*Yd*Yu.adjoint()*Yu) - 6*(TYu*Yu.
      adjoint()*Yu*Yu.adjoint()*Yu))).real();


   return beta_TYu;
}

/**
 * Calculates the 3-loop beta function of TYu.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_TYu_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = ZEROMATRIX(3,3);


   return beta_TYu;
}

/**
 * Calculates the 4-loop beta function of TYu.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_TYu_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = ZEROMATRIX(3,3);


   return beta_TYu;
}

/**
 * Calculates the 5-loop beta function of TYu.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_TYu_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = ZEROMATRIX(3,3);


   return beta_TYu;
}

} // namespace flexiblesusy
