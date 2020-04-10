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

// File generated at Fri 10 Apr 2020 20:18:33

#include "UMSSM_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of me2.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_me2_1_loop(const Soft_traces& soft_traces) const
{
   const auto Qe = INPUT(Qe);
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_me2;

   beta_me2 = (oneOver16PiSqr*(4*mHd2*(Ye*Ye.adjoint()) + 4*(TYe*(TYe).adjoint(
      )) + 2*(me2*Ye*Ye.adjoint()) + 4*(Ye*ml2*Ye.adjoint()) + 2*(Ye*Ye.adjoint
      ()*me2) + 0.4*(3.872983346207417*g1*Tr11 + 5*gp*Qe*Tr14 - 12*AbsSqr(MassB
      )*Sqr(g1) - 20*AbsSqr(MassU)*Sqr(gp)*Sqr(Qe))*UNITMATRIX(3))).real();


   return beta_me2;
}

/**
 * Calculates the 2-loop beta function of me2.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_me2_2_loop(const Soft_traces& soft_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qu = INPUT(Qu);
   const auto Qs = INPUT(Qs);
   const auto Qv = INPUT(Qv);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_me2;

   const Eigen::Matrix<double,3,3> beta_me2_1 = (UNITMATRIX(3)*(-0.8*twoLoop*(
      15*traceconjTYdTpTYd + 5*traceconjTYeTpTYe + 15*tracemd2YdAdjYd + 5*
      traceme2YeAdjYe + 5*traceml2AdjYeYe + 15*tracemq2AdjYdYd + 30*mHd2*
      traceYdAdjYd + 10*mHd2*traceYeAdjYe + 10*mHd2*AbsSqr(Lambdax) + 5*mHu2*
      AbsSqr(Lambdax) + 5*ms2*AbsSqr(Lambdax) + 5*AbsSqr(TLambdax) + 3*mHd2*Sqr
      (g1) + 6*AbsSqr(MassB)*Sqr(g1) - 15*mHd2*Sqr(g2) - 30*AbsSqr(MassWB)*Sqr(
      g2) + 10*mHd2*Sqr(gp)*Sqr(Qe) + 20*AbsSqr(MassU)*Sqr(gp)*Sqr(Qe) - 10*
      mHd2*Sqr(gp)*Sqr(QHd) - 20*AbsSqr(MassU)*Sqr(gp)*Sqr(QHd) - 10*mHd2*Sqr(
      gp)*Sqr(Ql) - 20*AbsSqr(MassU)*Sqr(gp)*Sqr(Ql))*(Ye*Ye.adjoint()) + 0.8*
      twoLoop*(-15*traceAdjYdTYd - 5*traceAdjYeTYe + 3*MassB*Sqr(g1) - 15*
      MassWB*Sqr(g2) + 10*MassU*Sqr(gp)*Sqr(Qe) - 10*MassU*Sqr(gp)*Sqr(QHd) -
      10*MassU*Sqr(gp)*Sqr(Ql) - 5*Conj(Lambdax)*TLambdax)*(Ye*(TYe).adjoint())
      - 0.8*twoLoop*(15*traceconjTYdTpYd + 5*traceconjTYeTpYe + 5*Conj(TLambdax
      )*Lambdax - 3*Conj(MassB)*Sqr(g1) + 15*Conj(MassWB)*Sqr(g2) - 10*Conj(
      MassU)*Sqr(gp)*Sqr(Qe) + 10*Conj(MassU)*Sqr(gp)*Sqr(QHd) + 10*Conj(MassU)
      *Sqr(gp)*Sqr(Ql))*(TYe*Ye.adjoint()) - 0.8*twoLoop*(15*traceYdAdjYd + 5*
      traceYeAdjYe + 5*AbsSqr(Lambdax) + 3*Sqr(g1) - 15*Sqr(g2) + 10*Sqr(gp)*
      Sqr(Qe) - 10*Sqr(gp)*Sqr(QHd) - 10*Sqr(gp)*Sqr(Ql))*(TYe*(TYe).adjoint())
      - 0.4*twoLoop*(15*traceYdAdjYd + 5*traceYeAdjYe + 5*AbsSqr(Lambdax) + 3*
      Sqr(g1) - 15*Sqr(g2) + 10*Sqr(gp)*Sqr(Qe) - 10*Sqr(gp)*Sqr(QHd) - 10*Sqr(
      gp)*Sqr(Ql))*(me2*Ye*Ye.adjoint()) - 0.8*twoLoop*(15*traceYdAdjYd + 5*
      traceYeAdjYe + 5*AbsSqr(Lambdax) + 3*Sqr(g1) - 15*Sqr(g2) + 10*Sqr(gp)*
      Sqr(Qe) - 10*Sqr(gp)*Sqr(QHd) - 10*Sqr(gp)*Sqr(Ql))*(Ye*ml2*Ye.adjoint())
      - 0.4*twoLoop*(15*traceYdAdjYd + 5*traceYeAdjYe + 5*AbsSqr(Lambdax) + 3*
      Sqr(g1) - 15*Sqr(g2) + 10*Sqr(gp)*Sqr(Qe) - 10*Sqr(gp)*Sqr(QHd) - 10*Sqr(
      gp)*Sqr(Ql))*(Ye*Ye.adjoint()*me2) - 8*mHd2*twoLoop*(Ye*Ye.adjoint()*Ye*
      Ye.adjoint()) - 4*twoLoop*(Ye*Ye.adjoint()*TYe*(TYe).adjoint()) - 4*
      twoLoop*(Ye*(TYe).adjoint()*TYe*Ye.adjoint()) - 4*(mHd2 + mHu2)*twoLoop*(
      Ye*Yv.conjugate()*Yv.transpose()*Ye.adjoint()) - 4*twoLoop*(Ye*Yv.
      conjugate()*(TYv).transpose()*(TYe).adjoint()) - 4*twoLoop*(Ye*TYv.
      conjugate()*(TYv).transpose()*Ye.adjoint()) - 4*twoLoop*(TYe*Ye.adjoint()
      *Ye*(TYe).adjoint()) - 4*twoLoop*(TYe*(TYe).adjoint()*Ye*Ye.adjoint()) -
      4*twoLoop*(TYe*Yv.conjugate()*Yv.transpose()*(TYe).adjoint()) - 4*twoLoop
      *(TYe*TYv.conjugate()*Yv.transpose()*Ye.adjoint()) - 2*twoLoop*(me2*Ye*Ye
      .adjoint()*Ye*Ye.adjoint()) - 2*twoLoop*(me2*Ye*Yv.conjugate()*Yv.
      transpose()*Ye.adjoint()) - 4*twoLoop*(Ye*ml2*Ye.adjoint()*Ye*Ye.adjoint(
      )) - 4*twoLoop*(Ye*ml2*Yv.conjugate()*Yv.transpose()*Ye.adjoint()) - 4*
      twoLoop*(Ye*Ye.adjoint()*me2*Ye*Ye.adjoint()) - 4*twoLoop*(Ye*Ye.adjoint(
      )*Ye*ml2*Ye.adjoint()) - 2*twoLoop*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*me2)
      - 4*twoLoop*(Ye*Yv.conjugate()*mvR2*Yv.transpose()*Ye.adjoint()) - 4*
      twoLoop*(Ye*Yv.conjugate()*Yv.transpose()*ml2*Ye.adjoint()) - 2*twoLoop*(
      Ye*Yv.conjugate()*Yv.transpose()*Ye.adjoint()*me2) + 0.32*twoLoop*(
      19.364916731037084*g1*gp*Qe*Tr2U114 + 19.364916731037084*g1*gp*Qe*Tr2U141
       + 19.364916731037084*g1*Tr31 + 25*gp*Qe*Tr34 + 351*AbsSqr(MassB)*Quad(g1
      ) + 15*Tr2U111*Sqr(g1) + 90*Qd*Qe*AbsSqr(MassB)*Sqr(g1)*Sqr(gp) + 45*
      MassU*Qd*Qe*Conj(MassB)*Sqr(g1)*Sqr(gp) + 25*Tr2U144*Sqr(gp)*Sqr(Qe))*
      UNITMATRIX(3))).real();
   const Eigen::Matrix<double,3,3> beta_me2_2 = (4.8*Qe*twoLoop*Sqr(gp)*(10*Qe*
      AbsSqr(MassB)*Sqr(g1) - 2*QHd*AbsSqr(MassB)*Sqr(g1) + 2*QHu*AbsSqr(MassB)
      *Sqr(g1) - 6*Ql*AbsSqr(MassB)*Sqr(g1) + 6*Qq*AbsSqr(MassB)*Sqr(g1) - 12*
      Qu*AbsSqr(MassB)*Sqr(g1) + 6*Qd*AbsSqr(MassU)*Sqr(g1) + 10*Qe*AbsSqr(
      MassU)*Sqr(g1) - 2*QHd*AbsSqr(MassU)*Sqr(g1) + 2*QHu*AbsSqr(MassU)*Sqr(g1
      ) - 6*Ql*AbsSqr(MassU)*Sqr(g1) + 6*Qq*AbsSqr(MassU)*Sqr(g1) - 12*Qu*
      AbsSqr(MassU)*Sqr(g1) + 5*MassU*Qe*Conj(MassB)*Sqr(g1) - MassU*QHd*Conj(
      MassB)*Sqr(g1) + MassU*QHu*Conj(MassB)*Sqr(g1) - 3*MassU*Ql*Conj(MassB)*
      Sqr(g1) + 3*MassU*Qq*Conj(MassB)*Sqr(g1) - 6*MassU*Qu*Conj(MassB)*Sqr(g1)
      + 3*MassB*Qd*Conj(MassU)*Sqr(g1) + 5*MassB*Qe*Conj(MassU)*Sqr(g1) - MassB
      *QHd*Conj(MassU)*Sqr(g1) + MassB*QHu*Conj(MassU)*Sqr(g1) - 3*MassB*Ql*
      Conj(MassU)*Sqr(g1) + 3*MassB*Qq*Conj(MassU)*Sqr(g1) - 6*MassB*Qu*Conj(
      MassU)*Sqr(g1) + 25*AbsSqr(MassU)*Cube(Qe)*Sqr(gp) + 45*Qe*AbsSqr(MassU)*
      Sqr(gp)*Sqr(Qd) + 10*Qe*AbsSqr(MassU)*Sqr(gp)*Sqr(QHd) + 10*Qe*AbsSqr(
      MassU)*Sqr(gp)*Sqr(QHu) + 30*Qe*AbsSqr(MassU)*Sqr(gp)*Sqr(Ql) + 90*Qe*
      AbsSqr(MassU)*Sqr(gp)*Sqr(Qq) + 5*Qe*AbsSqr(MassU)*Sqr(gp)*Sqr(Qs) + 45*
      Qe*AbsSqr(MassU)*Sqr(gp)*Sqr(Qu) + 15*Qe*AbsSqr(MassU)*Sqr(gp)*Sqr(Qv))*
      UNITMATRIX(3)).real();

   beta_me2 = beta_me2_1 + beta_me2_2;


   return beta_me2;
}

/**
 * Calculates the 3-loop beta function of me2.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_me2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_me2;

   beta_me2 = ZEROMATRIX(3,3);


   return beta_me2;
}

/**
 * Calculates the 4-loop beta function of me2.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_me2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_me2;

   beta_me2 = ZEROMATRIX(3,3);


   return beta_me2;
}

/**
 * Calculates the 5-loop beta function of me2.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_me2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_me2;

   beta_me2 = ZEROMATRIX(3,3);


   return beta_me2;
}

} // namespace flexiblesusy
