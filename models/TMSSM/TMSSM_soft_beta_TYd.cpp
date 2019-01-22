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

// File generated at Tue 22 Jan 2019 16:50:11

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

   beta_TYd = (oneOver16PiSqr*(0.03333333333333333*(180*traceAdjYdTYd*Yd + 60*
      traceAdjYeTYe*Yd + 28*MassB*Yd*Sqr(g1) + 180*MassWB*Yd*Sqr(g2) + 320*
      MassG*Yd*Sqr(g3) + 90*traceYdAdjYd*TYd + 30*traceYeAdjYe*TYd + 45*AbsSqr(
      Lambdax)*TYd - 14*Sqr(g1)*TYd - 90*Sqr(g2)*TYd - 160*Sqr(g3)*TYd + 90*Yd*
      Conj(Lambdax)*TLambdax) + 4*(Yd*Yd.adjoint()*TYd) + 2*(Yd*Yu.adjoint()*
      TYu) + 5*(TYd*Yd.adjoint()*Yd) + TYd*Yu.adjoint()*Yu)).real();


   return beta_TYd;
}

/**
 * Calculates the 2-loop beta function of TYd.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_TYd_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = (twoLoop*(0.005555555555555556*(-6480*traceYdAdjYdTYdAdjYd*Yd -
      1080*traceYdAdjYuTYuAdjYd*Yd - 2160*traceYeAdjYeTYeAdjYe*Yd - 1080*
      traceYuAdjYdTYdAdjYu*Yd - 1620*traceAdjYuTYu*Yd*AbsSqr(Lambdax) - 2296*
      MassB*Yd*Quad(g1) - 9720*MassWB*Yd*Quad(g2) + 1280*MassG*Yd*Quad(g3) -
      144*traceAdjYdTYd*Yd*Sqr(g1) + 432*traceAdjYeTYe*Yd*Sqr(g1) + 144*MassB*
      traceYdAdjYd*Yd*Sqr(g1) - 432*MassB*traceYeAdjYe*Yd*Sqr(g1) - 2160*MassWB
      *Yd*AbsSqr(Lambdax)*Sqr(g2) - 360*MassB*Yd*Sqr(g1)*Sqr(g2) - 360*MassWB*
      Yd*Sqr(g1)*Sqr(g2) + 5760*traceAdjYdTYd*Yd*Sqr(g3) - 5760*MassG*
      traceYdAdjYd*Yd*Sqr(g3) - 320*MassB*Yd*Sqr(g1)*Sqr(g3) - 320*MassG*Yd*Sqr
      (g1)*Sqr(g3) - 2880*MassG*Yd*Sqr(g2)*Sqr(g3) - 2880*MassWB*Yd*Sqr(g2)*Sqr
      (g3) - 1620*traceYdAdjYdYdAdjYd*TYd - 540*traceYdAdjYuYuAdjYd*TYd - 540*
      traceYeAdjYeYeAdjYe*TYd - 810*traceYuAdjYu*AbsSqr(Lambdax)*TYd + 574*Quad
      (g1)*TYd + 2430*Quad(g2)*TYd - 320*Quad(g3)*TYd - 72*traceYdAdjYd*Sqr(g1)
      *TYd + 216*traceYeAdjYe*Sqr(g1)*TYd + 1080*AbsSqr(Lambdax)*Sqr(g2)*TYd +
      180*Sqr(g1)*Sqr(g2)*TYd + 2880*traceYdAdjYd*Sqr(g3)*TYd + 160*Sqr(g1)*Sqr
      (g3)*TYd + 1440*Sqr(g2)*Sqr(g3)*TYd - 675*Sqr(Conj(Lambdax))*Sqr(Lambdax)
      *TYd - 1620*traceYuAdjYu*Yd*Conj(Lambdax)*TLambdax + 2160*Yd*Conj(Lambdax
      )*Sqr(g2)*TLambdax - 2700*Yd*Lambdax*Sqr(Conj(Lambdax))*TLambdax) + 0.2*(
      -90*traceAdjYdTYd - 30*traceAdjYeTYe - 8*MassB*Sqr(g1) - 60*MassWB*Sqr(g2
      ) - 45*Conj(Lambdax)*TLambdax)*(Yd*Yd.adjoint()*Yd) + 0.4*(-30*
      traceYdAdjYd - 10*traceYeAdjYe - 15*AbsSqr(Lambdax) + 3*Sqr(g1) + 15*Sqr(
      g2))*(Yd*Yd.adjoint()*TYd) + 0.2*(-30*traceAdjYuTYu - 8*MassB*Sqr(g1) -
      15*Conj(Lambdax)*TLambdax)*(Yd*Yu.adjoint()*Yu) + 0.2*(-30*traceYuAdjYu -
      15*AbsSqr(Lambdax) + 8*Sqr(g1))*(Yd*Yu.adjoint()*TYu) + 0.1*(-150*
      traceYdAdjYd - 50*traceYeAdjYe - 75*AbsSqr(Lambdax) + 12*Sqr(g1) + 120*
      Sqr(g2))*(TYd*Yd.adjoint()*Yd) + 0.1*(-30*traceYuAdjYu - 15*AbsSqr(
      Lambdax) + 8*Sqr(g1))*(TYd*Yu.adjoint()*Yu) - 6*(Yd*Yd.adjoint()*Yd*Yd.
      adjoint()*TYd) - 8*(Yd*Yd.adjoint()*TYd*Yd.adjoint()*Yd) - 2*(Yd*Yu.
      adjoint()*Yu*Yd.adjoint()*TYd) - 4*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*TYu)
      - 4*(Yd*Yu.adjoint()*TYu*Yd.adjoint()*Yd) - 4*(Yd*Yu.adjoint()*TYu*Yu.
      adjoint()*Yu) - 6*(TYd*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 4*(TYd*Yu.
      adjoint()*Yu*Yd.adjoint()*Yd) - 2*(TYd*Yu.adjoint()*Yu*Yu.adjoint()*Yu)))
      .real();


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

/**
 * Calculates the 4-loop beta function of TYd.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_TYd_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = ZEROMATRIX(3,3);


   return beta_TYd;
}

/**
 * Calculates the 5-loop beta function of TYd.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_TYd_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = ZEROMATRIX(3,3);


   return beta_TYd;
}

} // namespace flexiblesusy
