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

// File generated at Sun 4 Aug 2019 19:34:54

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
 * Calculates the 1-loop beta function of mvR2.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mvR2_1_loop(const Soft_traces& soft_traces) const
{
   const auto Qv = INPUT(Qv);
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_mvR2;

   beta_mvR2 = (oneOver16PiSqr*(4*mHu2*(Yv.transpose()*Yv.conjugate()) + 4*((
      TYv).transpose()*TYv.conjugate()) + 2*(mvR2*Yv.transpose()*Yv.conjugate()
      ) + 4*(Yv.transpose()*ml2*Yv.conjugate()) + 2*(Yv.transpose()*Yv.
      conjugate()*mvR2) - 2*gp*Qv*(-Tr14 + 4*gp*Qv*AbsSqr(MassU))*UNITMATRIX(3)
      )).real();


   return beta_mvR2;
}

/**
 * Calculates the 2-loop beta function of mvR2.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mvR2_2_loop(const Soft_traces& soft_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto Qv = INPUT(Qv);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qs = INPUT(Qs);
   const auto Qu = INPUT(Qu);
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double traceconjTYvTpTYv = TRACE_STRUCT.traceconjTYvTpTYv;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceYvAdjYvconjml2 = TRACE_STRUCT.traceYvAdjYvconjml2;
   const double traceYvconjmvR2AdjYv = TRACE_STRUCT.traceYvconjmvR2AdjYv;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYvTpYv = TRACE_STRUCT.traceconjTYvTpYv;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_mvR2;

   const Eigen::Matrix<double,3,3> beta_mvR2_1 = (UNITMATRIX(3)*(0.8*twoLoop*(-
      15*traceconjTYuTpTYu - 5*traceconjTYvTpTYv - 15*tracemq2AdjYuYu - 15*
      tracemu2YuAdjYu - 30*mHu2*traceYuAdjYu - 10*mHu2*traceYvAdjYv - 5*
      traceYvAdjYvconjml2 - 5*traceYvconjmvR2AdjYv - 5*mHd2*AbsSqr(Lambdax) -
      10*mHu2*AbsSqr(Lambdax) - 5*ms2*AbsSqr(Lambdax) - 5*AbsSqr(TLambdax) + 3*
      mHu2*Sqr(g1) + 6*AbsSqr(MassB)*Sqr(g1) + 15*mHu2*Sqr(g2) + 30*AbsSqr(
      MassWB)*Sqr(g2) + 10*mHu2*Sqr(gp)*Sqr(QHu) + 20*AbsSqr(MassU)*Sqr(gp)*Sqr
      (QHu) + 10*mHu2*Sqr(gp)*Sqr(Ql) + 20*AbsSqr(MassU)*Sqr(gp)*Sqr(Ql) - 10*
      mHu2*Sqr(gp)*Sqr(Qv) - 20*AbsSqr(MassU)*Sqr(gp)*Sqr(Qv))*(Yv.transpose()*
      Yv.conjugate()) - 0.8*twoLoop*(15*traceAdjYuTYu + 5*traceAdjYvTYv + 3*
      MassB*Sqr(g1) + 15*MassWB*Sqr(g2) + 10*MassU*Sqr(gp)*Sqr(QHu) + 10*MassU*
      Sqr(gp)*Sqr(Ql) - 10*MassU*Sqr(gp)*Sqr(Qv) + 5*Conj(Lambdax)*TLambdax)*(
      Yv.transpose()*TYv.conjugate()) - 0.8*twoLoop*(15*traceconjTYuTpYu + 5*
      traceconjTYvTpYv + 5*Conj(TLambdax)*Lambdax + 3*Conj(MassB)*Sqr(g1) + 15*
      Conj(MassWB)*Sqr(g2) + 10*Conj(MassU)*Sqr(gp)*Sqr(QHu) + 10*Conj(MassU)*
      Sqr(gp)*Sqr(Ql) - 10*Conj(MassU)*Sqr(gp)*Sqr(Qv))*((TYv).transpose()*Yv.
      conjugate()) + 0.8*twoLoop*(-15*traceYuAdjYu - 5*traceYvAdjYv - 5*AbsSqr(
      Lambdax) + 3*Sqr(g1) + 15*Sqr(g2) + 10*Sqr(gp)*Sqr(QHu) + 10*Sqr(gp)*Sqr(
      Ql) - 10*Sqr(gp)*Sqr(Qv))*((TYv).transpose()*TYv.conjugate()) + 0.4*
      twoLoop*(-15*traceYuAdjYu - 5*traceYvAdjYv - 5*AbsSqr(Lambdax) + 3*Sqr(g1
      ) + 15*Sqr(g2) + 10*Sqr(gp)*Sqr(QHu) + 10*Sqr(gp)*Sqr(Ql) - 10*Sqr(gp)*
      Sqr(Qv))*(mvR2*Yv.transpose()*Yv.conjugate()) + 0.8*twoLoop*(-15*
      traceYuAdjYu - 5*traceYvAdjYv - 5*AbsSqr(Lambdax) + 3*Sqr(g1) + 15*Sqr(g2
      ) + 10*Sqr(gp)*Sqr(QHu) + 10*Sqr(gp)*Sqr(Ql) - 10*Sqr(gp)*Sqr(Qv))*(Yv.
      transpose()*ml2*Yv.conjugate()) + 0.4*twoLoop*(-15*traceYuAdjYu - 5*
      traceYvAdjYv - 5*AbsSqr(Lambdax) + 3*Sqr(g1) + 15*Sqr(g2) + 10*Sqr(gp)*
      Sqr(QHu) + 10*Sqr(gp)*Sqr(Ql) - 10*Sqr(gp)*Sqr(Qv))*(Yv.transpose()*Yv.
      conjugate()*mvR2) - 4*(mHd2 + mHu2)*twoLoop*(Yv.transpose()*Ye.adjoint()*
      Ye*Yv.conjugate()) - 4*twoLoop*(Yv.transpose()*Ye.adjoint()*TYe*TYv.
      conjugate()) - 4*twoLoop*(Yv.transpose()*(TYe).adjoint()*TYe*Yv.conjugate
      ()) - 8*mHu2*twoLoop*(Yv.transpose()*Yv.conjugate()*Yv.transpose()*Yv.
      conjugate()) - 4*twoLoop*(Yv.transpose()*Yv.conjugate()*(TYv).transpose()
      *TYv.conjugate()) - 4*twoLoop*(Yv.transpose()*TYv.conjugate()*(TYv).
      transpose()*Yv.conjugate()) - 4*twoLoop*((TYv).transpose()*Ye.adjoint()*
      Ye*TYv.conjugate()) - 4*twoLoop*((TYv).transpose()*(TYe).adjoint()*Ye*Yv.
      conjugate()) - 4*twoLoop*((TYv).transpose()*Yv.conjugate()*Yv.transpose()
      *TYv.conjugate()) - 4*twoLoop*((TYv).transpose()*TYv.conjugate()*Yv.
      transpose()*Yv.conjugate()) - 2*twoLoop*(mvR2*Yv.transpose()*Ye.adjoint()
      *Ye*Yv.conjugate()) - 2*twoLoop*(mvR2*Yv.transpose()*Yv.conjugate()*Yv.
      transpose()*Yv.conjugate()) - 4*twoLoop*(Yv.transpose()*ml2*Ye.adjoint()*
      Ye*Yv.conjugate()) - 4*twoLoop*(Yv.transpose()*ml2*Yv.conjugate()*Yv.
      transpose()*Yv.conjugate()) - 4*twoLoop*(Yv.transpose()*Ye.adjoint()*me2*
      Ye*Yv.conjugate()) - 4*twoLoop*(Yv.transpose()*Ye.adjoint()*Ye*ml2*Yv.
      conjugate()) - 2*twoLoop*(Yv.transpose()*Ye.adjoint()*Ye*Yv.conjugate()*
      mvR2) - 4*twoLoop*(Yv.transpose()*Yv.conjugate()*mvR2*Yv.transpose()*Yv.
      conjugate()) - 4*twoLoop*(Yv.transpose()*Yv.conjugate()*Yv.transpose()*
      ml2*Yv.conjugate()) - 2*twoLoop*(Yv.transpose()*Yv.conjugate()*Yv.
      transpose()*Yv.conjugate()*mvR2) + 8*gp*Qv*twoLoop*(gp*Qv*Tr2U144 + Tr34
      + 27*Qv*AbsSqr(MassU)*Cube(gp)*Sqr(Qd) + 9*Qv*AbsSqr(MassU)*Cube(gp)*Sqr(
      Qe) + 6*Qv*AbsSqr(MassU)*Cube(gp)*Sqr(QHd) + 6*Qv*AbsSqr(MassU)*Cube(gp)*
      Sqr(QHu) + 18*Qv*AbsSqr(MassU)*Cube(gp)*Sqr(Ql) + 54*Qv*AbsSqr(MassU)*
      Cube(gp)*Sqr(Qq) + 3*Qv*AbsSqr(MassU)*Cube(gp)*Sqr(Qs))*UNITMATRIX(3))).
      real();
   const Eigen::Matrix<double,3,3> beta_mvR2_2 = (24*twoLoop*AbsSqr(MassU)*Quad
      (gp)*Sqr(Qv)*(9*Sqr(Qu) + 5*Sqr(Qv))*UNITMATRIX(3)).real();

   beta_mvR2 = beta_mvR2_1 + beta_mvR2_2;


   return beta_mvR2;
}

/**
 * Calculates the 3-loop beta function of mvR2.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mvR2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mvR2;

   beta_mvR2 = ZEROMATRIX(3,3);


   return beta_mvR2;
}

/**
 * Calculates the 4-loop beta function of mvR2.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mvR2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mvR2;

   beta_mvR2 = ZEROMATRIX(3,3);


   return beta_mvR2;
}

/**
 * Calculates the 5-loop beta function of mvR2.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mvR2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mvR2;

   beta_mvR2 = ZEROMATRIX(3,3);


   return beta_mvR2;
}

} // namespace flexiblesusy
