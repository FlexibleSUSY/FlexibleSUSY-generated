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

// File generated at Wed 16 Oct 2019 21:50:41

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
 * Calculates the 1-loop beta function of mu2.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_mu2_1_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;


   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = (oneOver16PiSqr*(4*mHu2*(Yu*Yu.adjoint()) + 4*(TYu*(TYu).adjoint(
      )) + 2*(mu2*Yu*Yu.adjoint()) + 4*(Yu*mq2*Yu.adjoint()) + 2*(Yu*Yu.adjoint
      ()*mu2) - 0.26666666666666666*(3.872983346207417*g1*Tr11 + 8*AbsSqr(MassB
      )*Sqr(g1) + 40*AbsSqr(MassG)*Sqr(g3))*UNITMATRIX(3))).real();


   return beta_mu2;
}

/**
 * Calculates the 2-loop beta function of mu2.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_mu2_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr23 = TRACE_STRUCT.Tr23;


   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = (twoLoop*(-0.4*(30*traceconjTYuTpTYu + 30*tracemq2AdjYuYu + 30*
      tracemu2YuAdjYu + 60*mHu2*traceYuAdjYu + 15*mHd2*AbsSqr(Lambdax) + 30*
      mHu2*AbsSqr(Lambdax) + 15*mT2*AbsSqr(Lambdax) + 15*AbsSqr(TLambdax) + 2*
      mHu2*Sqr(g1) + 4*AbsSqr(MassB)*Sqr(g1) - 30*mHu2*Sqr(g2) - 60*AbsSqr(
      MassWB)*Sqr(g2))*(Yu*Yu.adjoint()) + 0.4*(-30*traceAdjYuTYu + 2*MassB*Sqr
      (g1) - 30*MassWB*Sqr(g2) - 15*Conj(Lambdax)*TLambdax)*(Yu*(TYu).adjoint()
      ) - 0.4*(30*traceconjTYuTpYu + 15*Conj(TLambdax)*Lambdax - 2*Conj(MassB)*
      Sqr(g1) + 30*Conj(MassWB)*Sqr(g2))*(TYu*Yu.adjoint()) - 0.4*(30*
      traceYuAdjYu + 15*AbsSqr(Lambdax) + 2*Sqr(g1) - 30*Sqr(g2))*(TYu*(TYu).
      adjoint()) + 0.2*(-30*traceYuAdjYu - 15*AbsSqr(Lambdax) - 2*Sqr(g1) + 30*
      Sqr(g2))*(mu2*Yu*Yu.adjoint()) - 0.4*(30*traceYuAdjYu + 15*AbsSqr(Lambdax
      ) + 2*Sqr(g1) - 30*Sqr(g2))*(Yu*mq2*Yu.adjoint()) + 0.2*(-30*traceYuAdjYu
       - 15*AbsSqr(Lambdax) - 2*Sqr(g1) + 30*Sqr(g2))*(Yu*Yu.adjoint()*mu2) - 4
      *(mHd2 + mHu2)*(Yu*Yd.adjoint()*Yd*Yu.adjoint()) - 4*(Yu*Yd.adjoint()*TYd
      *(TYu).adjoint()) - 8*mHu2*(Yu*Yu.adjoint()*Yu*Yu.adjoint()) - 4*(Yu*Yu.
      adjoint()*TYu*(TYu).adjoint()) - 4*(Yu*(TYd).adjoint()*TYd*Yu.adjoint())
      - 4*(Yu*(TYu).adjoint()*TYu*Yu.adjoint()) - 4*(TYu*Yd.adjoint()*Yd*(TYu).
      adjoint()) - 4*(TYu*Yu.adjoint()*Yu*(TYu).adjoint()) - 4*(TYu*(TYd).
      adjoint()*Yd*Yu.adjoint()) - 4*(TYu*(TYu).adjoint()*Yu*Yu.adjoint()) - 2*
      (mu2*Yu*Yd.adjoint()*Yd*Yu.adjoint()) - 2*(mu2*Yu*Yu.adjoint()*Yu*Yu.
      adjoint()) - 4*(Yu*mq2*Yd.adjoint()*Yd*Yu.adjoint()) - 4*(Yu*mq2*Yu.
      adjoint()*Yu*Yu.adjoint()) - 4*(Yu*Yd.adjoint()*md2*Yd*Yu.adjoint()) - 4*
      (Yu*Yd.adjoint()*Yd*mq2*Yu.adjoint()) - 2*(Yu*Yd.adjoint()*Yd*Yu.adjoint(
      )*mu2) - 4*(Yu*Yu.adjoint()*mu2*Yu*Yu.adjoint()) - 4*(Yu*Yu.adjoint()*Yu*
      mq2*Yu.adjoint()) - 2*(Yu*Yu.adjoint()*Yu*Yu.adjoint()*mu2) +
      0.07111111111111111*(-58.09475019311125*g1*Tr31 + 642*AbsSqr(MassB)*Quad(
      g1) + 150*Tr23*Quad(g3) - 600*AbsSqr(MassG)*Quad(g3) + 30*Tr2U111*Sqr(g1)
      + 160*AbsSqr(MassB)*Sqr(g1)*Sqr(g3) + 160*AbsSqr(MassG)*Sqr(g1)*Sqr(g3) +
      80*MassG*Conj(MassB)*Sqr(g1)*Sqr(g3) + 80*MassB*Conj(MassG)*Sqr(g1)*Sqr(
      g3))*UNITMATRIX(3))).real();


   return beta_mu2;
}

/**
 * Calculates the 3-loop beta function of mu2.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_mu2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = ZEROMATRIX(3,3);


   return beta_mu2;
}

/**
 * Calculates the 4-loop beta function of mu2.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_mu2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = ZEROMATRIX(3,3);


   return beta_mu2;
}

/**
 * Calculates the 5-loop beta function of mu2.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_mu2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = ZEROMATRIX(3,3);


   return beta_mu2;
}

} // namespace flexiblesusy
