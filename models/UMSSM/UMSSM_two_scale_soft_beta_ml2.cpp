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

// File generated at Mon 19 Sep 2016 09:59:19

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
 * Calculates the one-loop beta function of ml2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_ml2_one_loop(const Soft_traces& soft_traces) const
{
   const auto Ql = INPUT(Ql);
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_ml2;

   beta_ml2 = (oneOver16PiSqr*(2*mHd2*(Ye.adjoint()*Ye) + 2*((TYe)
      .adjoint()*TYe) + 2*mHu2*(Yv.conjugate()*Yv.transpose()) + 2*(
      TYv.conjugate()*(TYv).transpose()) + ml2*Ye.adjoint()*Ye + ml2*
      Yv.conjugate()*Yv.transpose() + 2*(Ye.adjoint()*me2*Ye) + Ye.adjoint()*Ye
      *ml2 + 2*(Yv.conjugate()*mvR2*Yv.transpose()) + Yv.conjugate()*
      Yv.transpose()*ml2 - 0.7745966692414834*g1*Tr11*UNITMATRIX(3) + 2*gp*Ql*
      Tr14*UNITMATRIX(3) - 1.2*AbsSqr(MassB)*Sqr(g1)*UNITMATRIX(3) - 6*AbsSqr(
      MassWB)*Sqr(g2)*UNITMATRIX(3) - 8*AbsSqr(MassU)*Sqr(gp)*Sqr(Ql)*
      UNITMATRIX(3))).real();


   return beta_ml2;
}

/**
 * Calculates the two-loop beta function of ml2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_ml2_two_loop(const Soft_traces& soft_traces) const
{
   const auto Ql = INPUT(Ql);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto Qd = INPUT(Qd);
   const auto QHu = INPUT(QHu);
   const auto Qq = INPUT(Qq);
   const auto Qu = INPUT(Qu);
   const auto Qv = INPUT(Qv);
   const auto Qs = INPUT(Qs);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double traceconjTYvTpTYv = TRACE_STRUCT.traceconjTYvTpTYv;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceYvAdjYvconjml2 = TRACE_STRUCT.traceYvAdjYvconjml2;
   const double traceYvconjmvR2AdjYv = TRACE_STRUCT.traceYvconjmvR2AdjYv;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYvTpYv = TRACE_STRUCT.traceconjTYvTpYv;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_ml2;

   const Eigen::Matrix<double,3,3> beta_ml2_1 = ((0.4*twoLoop*(-15*
      traceconjTYdTpTYd - 5*traceconjTYeTpTYe - 15*tracemd2YdAdjYd - 5*
      traceme2YeAdjYe - 5*traceml2AdjYeYe - 15*tracemq2AdjYdYd - 30*mHd2*
      traceYdAdjYd - 10*mHd2*traceYeAdjYe - 10*mHd2*AbsSqr(Lambdax) - 5*mHu2*
      AbsSqr(Lambdax) - 5*ms2*AbsSqr(Lambdax) + 6*mHd2*Sqr(g1) + 12*AbsSqr(
      MassB)*Sqr(g1) + 10*mHd2*Sqr(gp)*Sqr(Qe) + 10*mHd2*Sqr(gp)*Sqr(QHd) + 20*
      AbsSqr(MassU)*Sqr(gp)*(Sqr(Qe) + Sqr(QHd) - Sqr(Ql)) - 10*mHd2*Sqr(gp)*
      Sqr(Ql))*(Ye.adjoint()*Ye) - 0.4*twoLoop*(6*Conj(MassB)*Sqr(g1) + 5*(3*
      traceconjTYdTpYd + traceconjTYeTpYe + Conj(TLambdax)*Lambdax + 2*Conj(
      MassU)*Sqr(gp)*(Sqr(Qe) + Sqr(QHd) - Sqr(Ql))))*(Ye.adjoint()*TYe) + 0.2*
      twoLoop*(-30*traceAdjYdTYd - 10*traceAdjYeTYe - 12*MassB*Sqr(g1) - 20*
      MassU*Sqr(gp)*Sqr(Qe) - 20*MassU*Sqr(gp)*Sqr(QHd) + 20*MassU*Sqr(gp)*Sqr(
      Ql))*((TYe).adjoint()*Ye) + 0.2*twoLoop*(-30*traceYdAdjYd - 10*
      traceYeAdjYe - 10*AbsSqr(Lambdax) + 12*Sqr(g1) + 20*Sqr(gp)*Sqr(Qe) + 20*
      Sqr(gp)*Sqr(QHd) - 20*Sqr(gp)*Sqr(Ql))*((TYe).adjoint()*TYe) + 0.2*
      twoLoop*(-30*traceconjTYuTpTYu - 10*traceconjTYvTpTYv - 30*
      tracemq2AdjYuYu - 30*tracemu2YuAdjYu - 60*mHu2*traceYuAdjYu - 20*mHu2*
      traceYvAdjYv - 10*traceYvAdjYvconjml2 - 10*traceYvconjmvR2AdjYv - 10*mHd2
      *AbsSqr(Lambdax) - 20*mHu2*AbsSqr(Lambdax) - 10*ms2*AbsSqr(Lambdax) + 20*
      mHu2*Sqr(gp)*Sqr(QHu) + 40*AbsSqr(MassU)*Sqr(gp)*Sqr(QHu) - 20*mHu2*Sqr(
      gp)*Sqr(Ql) - 40*AbsSqr(MassU)*Sqr(gp)*Sqr(Ql) + 20*mHu2*Sqr(gp)*Sqr(Qv)
      + 40*AbsSqr(MassU)*Sqr(gp)*Sqr(Qv))*(Yv.conjugate()*Yv.transpose()) + 0.2
      *twoLoop*(-30*traceconjTYuTpYu - 10*traceconjTYvTpYv - 10*Conj(TLambdax)*
      Lambdax - 20*Conj(MassU)*Sqr(gp)*Sqr(QHu) + 20*Conj(MassU)*Sqr(gp)*Sqr(Ql
      ) - 20*Conj(MassU)*Sqr(gp)*Sqr(Qv))*(Yv.conjugate()*(TYv).transpose()) +
      0.2*twoLoop*(-30*traceAdjYuTYu - 10*traceAdjYvTYv - 20*MassU*Sqr(gp)*Sqr(
      QHu) + 20*MassU*Sqr(gp)*Sqr(Ql) - 20*MassU*Sqr(gp)*Sqr(Qv))*(
      TYv.conjugate()*Yv.transpose()) + 0.2*twoLoop*(-30*traceYuAdjYu - 10*
      traceYvAdjYv - 10*AbsSqr(Lambdax) + 20*Sqr(gp)*Sqr(QHu) - 20*Sqr(gp)*Sqr(
      Ql) + 20*Sqr(gp)*Sqr(Qv))*(TYv.conjugate()*(TYv).transpose()) + 0.2*
      twoLoop*(-15*traceYdAdjYd - 5*traceYeAdjYe - 5*AbsSqr(Lambdax) + 6*Sqr(g1
      ) + 10*Sqr(gp)*Sqr(Qe) + 10*Sqr(gp)*Sqr(QHd) - 10*Sqr(gp)*Sqr(Ql))*(ml2*
      Ye.adjoint()*Ye) + 0.2*twoLoop*(-15*traceYuAdjYu - 5*traceYvAdjYv - 5*
      AbsSqr(Lambdax) + 10*Sqr(gp)*Sqr(QHu) - 10*Sqr(gp)*Sqr(Ql) + 10*Sqr(gp)*
      Sqr(Qv))*(ml2*Yv.conjugate()*Yv.transpose()) + 0.2*twoLoop*(-30*
      traceYdAdjYd - 10*traceYeAdjYe - 10*AbsSqr(Lambdax) + 12*Sqr(g1) + 20*Sqr
      (gp)*Sqr(Qe) + 20*Sqr(gp)*Sqr(QHd) - 20*Sqr(gp)*Sqr(Ql))*(Ye.adjoint()*
      me2*Ye) + 0.2*twoLoop*(-15*traceYdAdjYd - 5*traceYeAdjYe - 5*AbsSqr(
      Lambdax) + 6*Sqr(g1) + 10*Sqr(gp)*Sqr(Qe) + 10*Sqr(gp)*Sqr(QHd) - 10*Sqr(
      gp)*Sqr(Ql))*(Ye.adjoint()*Ye*ml2))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_ml2_2 = (UNITMATRIX(3)*(-2*
      twoLoop*AbsSqr(TLambdax)*(Ye.adjoint()*Ye) - 2*twoLoop*Conj(Lambdax)*
      TLambdax*((TYe).adjoint()*Ye) - 2*twoLoop*AbsSqr(TLambdax)*(Yv.conjugate(
      )*Yv.transpose()) - 2*twoLoop*Conj(Lambdax)*TLambdax*(TYv.conjugate()*
      Yv.transpose()) + 2*twoLoop*(-3*traceYuAdjYu - traceYvAdjYv - AbsSqr(
      Lambdax) + 2*Sqr(gp)*(Sqr(QHu) - Sqr(Ql) + Sqr(Qv)))*(Yv.conjugate()*mvR2
      *Yv.transpose()) + twoLoop*(-3*traceYuAdjYu - traceYvAdjYv - AbsSqr(
      Lambdax) + 2*Sqr(gp)*(Sqr(QHu) - Sqr(Ql) + Sqr(Qv)))*(Yv.conjugate()*
      Yv.transpose()*ml2) - 8*mHd2*twoLoop*(Ye.adjoint()*Ye*Ye.adjoint()*Ye) -
      4*twoLoop*(Ye.adjoint()*Ye*(TYe).adjoint()*TYe) - 4*twoLoop*(Ye.adjoint()
      *TYe*(TYe).adjoint()*Ye) - 4*twoLoop*((TYe).adjoint()*Ye*Ye.adjoint()*TYe
      ) - 4*twoLoop*((TYe).adjoint()*TYe*Ye.adjoint()*Ye) - 8*mHu2*twoLoop*(
      Yv.conjugate()*Yv.transpose()*Yv.conjugate()*Yv.transpose()) - 4*twoLoop*
      (Yv.conjugate()*Yv.transpose()*TYv.conjugate()*(TYv).transpose()) - 4*
      twoLoop*(Yv.conjugate()*(TYv).transpose()*TYv.conjugate()*Yv.transpose())
      - 4*twoLoop*(TYv.conjugate()*Yv.transpose()*Yv.conjugate()*(TYv)
      .transpose()) - 4*twoLoop*(TYv.conjugate()*(TYv).transpose()*Yv.conjugate
      ()*Yv.transpose()) - 2*twoLoop*(ml2*Ye.adjoint()*Ye*Ye.adjoint()*Ye) - 2*
      twoLoop*(ml2*Yv.conjugate()*Yv.transpose()*Yv.conjugate()*Yv.transpose())
      - 4*twoLoop*(Ye.adjoint()*me2*Ye*Ye.adjoint()*Ye) - 4*twoLoop*(
      Ye.adjoint()*Ye*ml2*Ye.adjoint()*Ye) - 4*twoLoop*(Ye.adjoint()*Ye*
      Ye.adjoint()*me2*Ye) - 2*twoLoop*(Ye.adjoint()*Ye*Ye.adjoint()*Ye*ml2) -
      4*twoLoop*(Yv.conjugate()*mvR2*Yv.transpose()*Yv.conjugate()*Yv.transpose
      ()) - 4*twoLoop*(Yv.conjugate()*Yv.transpose()*ml2*Yv.conjugate()*
      Yv.transpose()) - 4*twoLoop*(Yv.conjugate()*Yv.transpose()*Yv.conjugate()
      *mvR2*Yv.transpose()) - 2*twoLoop*(Yv.conjugate()*Yv.transpose()*
      Yv.conjugate()*Yv.transpose()*ml2) + 0.04*twoLoop*(150*Power(g2,4)*Tr22*
      UNITMATRIX(3) - 77.45966692414834*g1*gp*Ql*Tr2U114*UNITMATRIX(3) -
      77.45966692414834*g1*gp*Ql*Tr2U141*UNITMATRIX(3) - 77.45966692414834*g1*
      Tr31*UNITMATRIX(3) + 200*gp*Ql*Tr34*UNITMATRIX(3) + 621*Power(g1,4)*
      AbsSqr(MassB)*UNITMATRIX(3) + 4800*Power(gp,4)*Power(Ql,4)*AbsSqr(MassU)*
      UNITMATRIX(3) + 825*Power(g2,4)*AbsSqr(MassWB)*UNITMATRIX(3) + 30*Tr2U111
      *Sqr(g1)*UNITMATRIX(3) + 90*AbsSqr(MassB)*Sqr(g1)*Sqr(g2)*UNITMATRIX(3) +
      90*AbsSqr(MassWB)*Sqr(g1)*Sqr(g2)*UNITMATRIX(3) + 45*MassWB*Conj(MassB)*
      Sqr(g1)*Sqr(g2)*UNITMATRIX(3) + 45*MassB*Conj(MassWB)*Sqr(g1)*Sqr(g2)*
      UNITMATRIX(3) - 360*Qd*Ql*AbsSqr(MassB)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) -
      360*Qe*Ql*AbsSqr(MassB)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) + 120*QHd*Ql*AbsSqr
      (MassB)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) - 120*QHu*Ql*AbsSqr(MassB)*Sqr(g1)*
      Sqr(gp)*UNITMATRIX(3) - 360*Ql*Qq*AbsSqr(MassB)*Sqr(g1)*Sqr(gp)*
      UNITMATRIX(3) + 720*Ql*Qu*AbsSqr(MassB)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) -
      360*Qd*Ql*AbsSqr(MassU)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) - 360*Qe*Ql*AbsSqr(
      MassU)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) + 120*QHd*Ql*AbsSqr(MassU)*Sqr(g1)*
      Sqr(gp)*UNITMATRIX(3) - 120*QHu*Ql*AbsSqr(MassU)*Sqr(g1)*Sqr(gp)*
      UNITMATRIX(3) - 360*Ql*Qq*AbsSqr(MassU)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) +
      720*Ql*Qu*AbsSqr(MassU)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) - 180*MassU*Qd*Ql*
      Conj(MassB)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) - 180*MassU*Qe*Ql*Conj(MassB)*
      Sqr(g1)*Sqr(gp)*UNITMATRIX(3) + 60*MassU*QHd*Ql*Conj(MassB)*Sqr(g1)*Sqr(
      gp)*UNITMATRIX(3) - 60*MassU*QHu*Ql*Conj(MassB)*Sqr(g1)*Sqr(gp)*
      UNITMATRIX(3) - 180*MassU*Ql*Qq*Conj(MassB)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3)
      + 360*MassU*Ql*Qu*Conj(MassB)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) - 180*MassB*
      Qd*Ql*Conj(MassU)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) - 180*MassB*Qe*Ql*Conj(
      MassU)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) + 60*MassB*QHd*Ql*Conj(MassU)*Sqr(g1
      )*Sqr(gp)*UNITMATRIX(3) - 60*MassB*QHu*Ql*Conj(MassU)*Sqr(g1)*Sqr(gp)*
      UNITMATRIX(3) - 180*MassB*Ql*Qq*Conj(MassU)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3)
      + 360*MassB*Ql*Qu*Conj(MassU)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) + 200*
      Tr2U144*Sqr(gp)*Sqr(Ql)*UNITMATRIX(3) + 480*AbsSqr(MassB)*Sqr(g1)*Sqr(gp)
      *Sqr(Ql)*UNITMATRIX(3) + 480*AbsSqr(MassU)*Sqr(g1)*Sqr(gp)*Sqr(Ql)*
      UNITMATRIX(3) + 240*MassU*Conj(MassB)*Sqr(g1)*Sqr(gp)*Sqr(Ql)*UNITMATRIX(
      3) + 240*MassB*Conj(MassU)*Sqr(g1)*Sqr(gp)*Sqr(Ql)*UNITMATRIX(3) + 600*
      AbsSqr(MassU)*Sqr(g2)*Sqr(gp)*Sqr(Ql)*UNITMATRIX(3) + 600*AbsSqr(MassWB)*
      Sqr(g2)*Sqr(gp)*Sqr(Ql)*UNITMATRIX(3) + 300*MassWB*Conj(MassU)*Sqr(g2)*
      Sqr(gp)*Sqr(Ql)*UNITMATRIX(3) + 300*MassU*Conj(MassWB)*Sqr(g2)*Sqr(gp)*
      Sqr(Ql)*UNITMATRIX(3) + 5400*Power(gp,4)*AbsSqr(MassU)*Sqr(Qd)*Sqr(Ql)*
      UNITMATRIX(3) + 1800*Power(gp,4)*AbsSqr(MassU)*Sqr(Qe)*Sqr(Ql)*UNITMATRIX
      (3) + 1200*Power(gp,4)*AbsSqr(MassU)*Sqr(QHd)*Sqr(Ql)*UNITMATRIX(3) +
      1200*Power(gp,4)*AbsSqr(MassU)*Sqr(QHu)*Sqr(Ql)*UNITMATRIX(3) + 10800*
      Power(gp,4)*AbsSqr(MassU)*Sqr(Ql)*Sqr(Qq)*UNITMATRIX(3) + 600*Power(gp,4)
      *AbsSqr(MassU)*Sqr(Ql)*Sqr(Qs)*UNITMATRIX(3) + 5400*Power(gp,4)*AbsSqr(
      MassU)*Sqr(Ql)*Sqr(Qu)*UNITMATRIX(3) + 1800*Power(gp,4)*AbsSqr(MassU)*Sqr
      (Ql)*Sqr(Qv)*UNITMATRIX(3)))).real();

   beta_ml2 = beta_ml2_1 + beta_ml2_2;


   return beta_ml2;
}

/**
 * Calculates the three-loop beta function of ml2.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_ml2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_ml2;

   beta_ml2 = ZEROMATRIX(3,3);


   return beta_ml2;
}

} // namespace flexiblesusy
