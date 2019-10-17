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

// File generated at Wed 16 Oct 2019 21:50:26

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
 * Calculates the 1-loop beta function of TYe.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_TYe_1_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = (oneOver16PiSqr*(0.1*(60*traceAdjYdTYd*Ye + 20*traceAdjYeTYe*Ye +
      36*MassB*Ye*Sqr(g1) + 60*MassWB*Ye*Sqr(g2) + 30*traceYdAdjYd*TYe + 10*
      traceYeAdjYe*TYe + 15*AbsSqr(Lambdax)*TYe - 18*Sqr(g1)*TYe - 30*Sqr(g2)*
      TYe + 30*Ye*Conj(Lambdax)*TLambdax) + 4*(Ye*Ye.adjoint()*TYe) + 5*(TYe*Ye
      .adjoint()*Ye))).real();


   return beta_TYe;
}

/**
 * Calculates the 2-loop beta function of TYe.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_TYe_2_loop(const Soft_traces& soft_traces) const
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


   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = (twoLoop*(0.05*(-720*traceYdAdjYdTYdAdjYd*Ye - 120*
      traceYdAdjYuTYuAdjYd*Ye - 240*traceYeAdjYeTYeAdjYe*Ye - 120*
      traceYuAdjYdTYdAdjYu*Ye - 180*traceAdjYuTYu*Ye*AbsSqr(Lambdax) - 1080*
      MassB*Ye*Quad(g1) - 1080*MassWB*Ye*Quad(g2) - 16*traceAdjYdTYd*Ye*Sqr(g1)
      + 48*traceAdjYeTYe*Ye*Sqr(g1) + 16*MassB*traceYdAdjYd*Ye*Sqr(g1) - 48*
      MassB*traceYeAdjYe*Ye*Sqr(g1) - 240*MassWB*Ye*AbsSqr(Lambdax)*Sqr(g2) -
      72*MassB*Ye*Sqr(g1)*Sqr(g2) - 72*MassWB*Ye*Sqr(g1)*Sqr(g2) + 640*
      traceAdjYdTYd*Ye*Sqr(g3) - 640*MassG*traceYdAdjYd*Ye*Sqr(g3) - 180*
      traceYdAdjYdYdAdjYd*TYe - 60*traceYdAdjYuYuAdjYd*TYe - 60*
      traceYeAdjYeYeAdjYe*TYe - 90*traceYuAdjYu*AbsSqr(Lambdax)*TYe + 270*Quad(
      g1)*TYe + 270*Quad(g2)*TYe - 8*traceYdAdjYd*Sqr(g1)*TYe + 24*traceYeAdjYe
      *Sqr(g1)*TYe + 120*AbsSqr(Lambdax)*Sqr(g2)*TYe + 36*Sqr(g1)*Sqr(g2)*TYe +
      320*traceYdAdjYd*Sqr(g3)*TYe - 75*Sqr(Conj(Lambdax))*Sqr(Lambdax)*TYe -
      180*traceYuAdjYu*Ye*Conj(Lambdax)*TLambdax + 240*Ye*Conj(Lambdax)*Sqr(g2)
      *TLambdax - 300*Ye*Lambdax*Sqr(Conj(Lambdax))*TLambdax) - 3*(6*
      traceAdjYdTYd + 2*traceAdjYeTYe + 4*MassWB*Sqr(g2) + 3*Conj(Lambdax)*
      TLambdax)*(Ye*Ye.adjoint()*Ye) + 0.4*(-30*traceYdAdjYd - 10*traceYeAdjYe
      - 15*AbsSqr(Lambdax) + 3*Sqr(g1) + 15*Sqr(g2))*(Ye*Ye.adjoint()*TYe) +
      0.1*(-150*traceYdAdjYd - 50*traceYeAdjYe - 75*AbsSqr(Lambdax) - 12*Sqr(g1
      ) + 120*Sqr(g2))*(TYe*Ye.adjoint()*Ye) - 6*(Ye*Ye.adjoint()*Ye*Ye.adjoint
      ()*TYe) - 8*(Ye*Ye.adjoint()*TYe*Ye.adjoint()*Ye) - 6*(TYe*Ye.adjoint()*
      Ye*Ye.adjoint()*Ye))).real();


   return beta_TYe;
}

/**
 * Calculates the 3-loop beta function of TYe.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_TYe_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = ZEROMATRIX(3,3);


   return beta_TYe;
}

/**
 * Calculates the 4-loop beta function of TYe.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_TYe_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = ZEROMATRIX(3,3);


   return beta_TYe;
}

/**
 * Calculates the 5-loop beta function of TYe.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_TYe_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = ZEROMATRIX(3,3);


   return beta_TYe;
}

} // namespace flexiblesusy
