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

// File generated at Sun 28 Aug 2016 15:09:12

#include "UMSSM_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of mvR2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mvR2_one_loop(const Soft_traces& soft_traces) const
{
   const auto Qv = INPUT(Qv);
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_mvR2;

   beta_mvR2 = (oneOver16PiSqr*(4*mHu2*(Yv.transpose()*Yv.conjugate()) +
      4*((TYv).transpose()*TYv.conjugate()) + 2*(mvR2*Yv.transpose()*
      Yv.conjugate()) + 4*(Yv.transpose()*ml2*Yv.conjugate()) + 2*(Yv.transpose
      ()*Yv.conjugate()*mvR2) + 2*(gp*Qv*Tr14*UNITMATRIX(3) - 4*AbsSqr(MassU)*
      Sqr(gp)*Sqr(Qv)*UNITMATRIX(3)))).real();


   return beta_mvR2;
}

/**
 * Calculates the two-loop beta function of mvR2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mvR2_two_loop(const Soft_traces& soft_traces) const
{
   const auto Qv = INPUT(Qv);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
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

   const Eigen::Matrix<double,3,3> beta_mvR2_1 = (UNITMATRIX(3)*(-0.8*
      twoLoop*(15*traceconjTYuTpTYu + 5*traceconjTYvTpTYv + 15*tracemq2AdjYuYu
      + 15*tracemu2YuAdjYu + 30*mHu2*traceYuAdjYu + 10*mHu2*traceYvAdjYv + 5*
      traceYvAdjYvconjml2 + 5*traceYvconjmvR2AdjYv + 5*mHd2*AbsSqr(Lambdax) +
      10*mHu2*AbsSqr(Lambdax) + 5*ms2*AbsSqr(Lambdax) + 5*AbsSqr(TLambdax) - 3*
      mHu2*Sqr(g1) - 6*AbsSqr(MassB)*Sqr(g1) - 15*mHu2*Sqr(g2) - 30*AbsSqr(
      MassWB)*Sqr(g2) - 10*mHu2*Sqr(gp)*Sqr(QHu) - 10*mHu2*Sqr(gp)*Sqr(Ql) - 20
      *AbsSqr(MassU)*Sqr(gp)*(Sqr(QHu) + Sqr(Ql) - Sqr(Qv)) + 10*mHu2*Sqr(gp)*
      Sqr(Qv))*(Yv.transpose()*Yv.conjugate()) - 0.8*twoLoop*(3*MassB*Sqr(g1) +
      5*(3*traceAdjYuTYu + traceAdjYvTYv + 3*MassWB*Sqr(g2) + 2*MassU*Sqr(gp)*
      (Sqr(QHu) + Sqr(Ql) - Sqr(Qv))) + 5*Conj(Lambdax)*TLambdax)*(Yv.transpose
      ()*TYv.conjugate()) - 0.4*twoLoop*(30*traceconjTYuTpYu + 10*
      traceconjTYvTpYv + 10*Conj(TLambdax)*Lambdax + 6*Conj(MassB)*Sqr(g1) + 30
      *Conj(MassWB)*Sqr(g2) + 20*Conj(MassU)*Sqr(gp)*Sqr(QHu) + 20*Conj(MassU)*
      Sqr(gp)*Sqr(Ql) - 20*Conj(MassU)*Sqr(gp)*Sqr(Qv))*((TYv).transpose()*
      Yv.conjugate()) - 0.4*twoLoop*(30*traceYuAdjYu + 10*traceYvAdjYv + 10*
      AbsSqr(Lambdax) - 6*Sqr(g1) - 30*Sqr(g2) - 20*Sqr(gp)*Sqr(QHu) - 20*Sqr(
      gp)*Sqr(Ql) + 20*Sqr(gp)*Sqr(Qv))*((TYv).transpose()*TYv.conjugate()) -
      0.4*twoLoop*(15*traceYuAdjYu + 5*traceYvAdjYv + 5*AbsSqr(Lambdax) - 3*Sqr
      (g1) - 15*Sqr(g2) - 10*Sqr(gp)*Sqr(QHu) - 10*Sqr(gp)*Sqr(Ql) + 10*Sqr(gp)
      *Sqr(Qv))*(mvR2*Yv.transpose()*Yv.conjugate()) - 0.4*twoLoop*(30*
      traceYuAdjYu + 10*traceYvAdjYv + 10*AbsSqr(Lambdax) - 6*Sqr(g1) - 30*Sqr(
      g2) - 20*Sqr(gp)*Sqr(QHu) - 20*Sqr(gp)*Sqr(Ql) + 20*Sqr(gp)*Sqr(Qv))*(
      Yv.transpose()*ml2*Yv.conjugate()) - 0.4*twoLoop*(15*traceYuAdjYu + 5*
      traceYvAdjYv + 5*AbsSqr(Lambdax) - 3*Sqr(g1) - 15*Sqr(g2) - 10*Sqr(gp)*
      Sqr(QHu) - 10*Sqr(gp)*Sqr(Ql) + 10*Sqr(gp)*Sqr(Qv))*(Yv.transpose()*
      Yv.conjugate()*mvR2) - 0.4*(10*mHd2 + 10*mHu2)*twoLoop*(Yv.transpose()*
      Ye.adjoint()*Ye*Yv.conjugate()) - 4*twoLoop*(Yv.transpose()*Ye.adjoint()*
      TYe*TYv.conjugate()) - 4*twoLoop*(Yv.transpose()*(TYe).adjoint()*TYe*
      Yv.conjugate()) - 8*mHu2*twoLoop*(Yv.transpose()*Yv.conjugate()*
      Yv.transpose()*Yv.conjugate()) - 4*twoLoop*(Yv.transpose()*Yv.conjugate()
      *(TYv).transpose()*TYv.conjugate()) - 4*twoLoop*(Yv.transpose()*
      TYv.conjugate()*(TYv).transpose()*Yv.conjugate()) - 4*twoLoop*((TYv)
      .transpose()*Ye.adjoint()*Ye*TYv.conjugate()) - 4*twoLoop*((TYv)
      .transpose()*(TYe).adjoint()*Ye*Yv.conjugate()) - 4*twoLoop*((TYv)
      .transpose()*Yv.conjugate()*Yv.transpose()*TYv.conjugate()) - 4*twoLoop*(
      (TYv).transpose()*TYv.conjugate()*Yv.transpose()*Yv.conjugate()) - 2*
      twoLoop*(mvR2*Yv.transpose()*Ye.adjoint()*Ye*Yv.conjugate()) - 2*twoLoop*
      (mvR2*Yv.transpose()*Yv.conjugate()*Yv.transpose()*Yv.conjugate()) - 4*
      twoLoop*(Yv.transpose()*ml2*Ye.adjoint()*Ye*Yv.conjugate()) - 4*twoLoop*(
      Yv.transpose()*ml2*Yv.conjugate()*Yv.transpose()*Yv.conjugate()) - 4*
      twoLoop*(Yv.transpose()*Ye.adjoint()*me2*Ye*Yv.conjugate()) - 4*twoLoop*(
      Yv.transpose()*Ye.adjoint()*Ye*ml2*Yv.conjugate()) - 2*twoLoop*(
      Yv.transpose()*Ye.adjoint()*Ye*Yv.conjugate()*mvR2) - 4*twoLoop*(
      Yv.transpose()*Yv.conjugate()*mvR2*Yv.transpose()*Yv.conjugate()) - 4*
      twoLoop*(Yv.transpose()*Yv.conjugate()*Yv.transpose()*ml2*Yv.conjugate())
      - 2*twoLoop*(Yv.transpose()*Yv.conjugate()*Yv.transpose()*Yv.conjugate()
      *mvR2) - 0.4*twoLoop*(-20*gp*Qv*Tr34*UNITMATRIX(3) - 20*Tr2U144*Sqr(gp)*
      Sqr(Qv)*UNITMATRIX(3) - 540*Power(gp,4)*AbsSqr(MassU)*Sqr(Qd)*Sqr(Qv)*
      UNITMATRIX(3) - 180*Power(gp,4)*AbsSqr(MassU)*Sqr(Qe)*Sqr(Qv)*UNITMATRIX(
      3) - 120*Power(gp,4)*AbsSqr(MassU)*Sqr(QHd)*Sqr(Qv)*UNITMATRIX(3) - 120*
      Power(gp,4)*AbsSqr(MassU)*Sqr(QHu)*Sqr(Qv)*UNITMATRIX(3) - 360*Power(gp,4
      )*AbsSqr(MassU)*Sqr(Ql)*Sqr(Qv)*UNITMATRIX(3) - 1080*Power(gp,4)*AbsSqr(
      MassU)*Sqr(Qq)*Sqr(Qv)*UNITMATRIX(3) - 60*Power(gp,4)*AbsSqr(MassU)*Sqr(
      Qs)*Sqr(Qv)*UNITMATRIX(3)))).real();
   const Eigen::Matrix<double,3,3> beta_mvR2_2 = (24*Power(gp,4)*twoLoop*
      AbsSqr(MassU)*Sqr(Qv)*(9*Sqr(Qu) + 5*Sqr(Qv))*UNITMATRIX(3)).real();

   beta_mvR2 = beta_mvR2_1 + beta_mvR2_2;


   return beta_mvR2;
}

/**
 * Calculates the three-loop beta function of mvR2.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mvR2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mvR2;

   beta_mvR2 = ZEROMATRIX(3,3);


   return beta_mvR2;
}

} // namespace flexiblesusy
