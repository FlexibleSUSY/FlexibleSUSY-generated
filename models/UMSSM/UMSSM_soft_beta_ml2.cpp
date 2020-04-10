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

// File generated at Fri 10 Apr 2020 20:18:24

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
 * Calculates the 1-loop beta function of ml2.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_ml2_1_loop(const Soft_traces& soft_traces) const
{
   const auto Ql = INPUT(Ql);
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_ml2;

   beta_ml2 = (oneOver16PiSqr*(2*mHd2*(Ye.adjoint()*Ye) + 2*((TYe).adjoint()*
      TYe) + 2*mHu2*(Yv.conjugate()*Yv.transpose()) + 2*(TYv.conjugate()*(TYv).
      transpose()) + ml2*Ye.adjoint()*Ye + ml2*Yv.conjugate()*Yv.transpose() +
      2*(Ye.adjoint()*me2*Ye) + Ye.adjoint()*Ye*ml2 + 2*(Yv.conjugate()*mvR2*Yv
      .transpose()) + Yv.conjugate()*Yv.transpose()*ml2 - 0.2*(
      3.872983346207417*g1*Tr11 - 10*gp*Ql*Tr14 + 6*AbsSqr(MassB)*Sqr(g1) + 30*
      AbsSqr(MassWB)*Sqr(g2) + 40*AbsSqr(MassU)*Sqr(gp)*Sqr(Ql))*UNITMATRIX(3))
      ).real();


   return beta_ml2;
}

/**
 * Calculates the 2-loop beta function of ml2.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_ml2_2_loop(const Soft_traces& soft_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto Ql = INPUT(Ql);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
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
      MassB)*Sqr(g1) + 10*mHd2*Sqr(gp)*Sqr(Qe) + 20*AbsSqr(MassU)*Sqr(gp)*Sqr(
      Qe) + 10*mHd2*Sqr(gp)*Sqr(QHd) + 20*AbsSqr(MassU)*Sqr(gp)*Sqr(QHd) - 10*
      mHd2*Sqr(gp)*Sqr(Ql) - 20*AbsSqr(MassU)*Sqr(gp)*Sqr(Ql))*(Ye.adjoint()*Ye
      ) - 0.4*twoLoop*(15*traceconjTYdTpYd + 5*traceconjTYeTpYe + 5*Conj(
      TLambdax)*Lambdax + 6*Conj(MassB)*Sqr(g1) + 10*Conj(MassU)*Sqr(gp)*Sqr(Qe
      ) + 10*Conj(MassU)*Sqr(gp)*Sqr(QHd) - 10*Conj(MassU)*Sqr(gp)*Sqr(Ql))*(Ye
      .adjoint()*TYe) - 0.4*twoLoop*(15*traceAdjYdTYd + 5*traceAdjYeTYe + 6*
      MassB*Sqr(g1) + 10*MassU*Sqr(gp)*Sqr(Qe) + 10*MassU*Sqr(gp)*Sqr(QHd) - 10
      *MassU*Sqr(gp)*Sqr(Ql))*((TYe).adjoint()*Ye) + 0.4*twoLoop*(-15*
      traceYdAdjYd - 5*traceYeAdjYe - 5*AbsSqr(Lambdax) + 6*Sqr(g1) + 10*Sqr(gp
      )*Sqr(Qe) + 10*Sqr(gp)*Sqr(QHd) - 10*Sqr(gp)*Sqr(Ql))*((TYe).adjoint()*
      TYe) + 2*twoLoop*(-3*traceconjTYuTpTYu - traceconjTYvTpTYv - 3*
      tracemq2AdjYuYu - 3*tracemu2YuAdjYu - 6*mHu2*traceYuAdjYu - 2*mHu2*
      traceYvAdjYv - traceYvAdjYvconjml2 - traceYvconjmvR2AdjYv - mHd2*AbsSqr(
      Lambdax) - 2*mHu2*AbsSqr(Lambdax) - ms2*AbsSqr(Lambdax) + 2*mHu2*Sqr(gp)*
      Sqr(QHu) + 4*AbsSqr(MassU)*Sqr(gp)*Sqr(QHu) - 2*mHu2*Sqr(gp)*Sqr(Ql) - 4*
      AbsSqr(MassU)*Sqr(gp)*Sqr(Ql) + 2*mHu2*Sqr(gp)*Sqr(Qv) + 4*AbsSqr(MassU)*
      Sqr(gp)*Sqr(Qv))*(Yv.conjugate()*Yv.transpose()) - 2*twoLoop*(3*
      traceconjTYuTpYu + traceconjTYvTpYv + Conj(TLambdax)*Lambdax + 2*Conj(
      MassU)*Sqr(gp)*Sqr(QHu) - 2*Conj(MassU)*Sqr(gp)*Sqr(Ql) + 2*Conj(MassU)*
      Sqr(gp)*Sqr(Qv))*(Yv.conjugate()*(TYv).transpose()) - 2*twoLoop*(3*
      traceAdjYuTYu + traceAdjYvTYv + 2*MassU*Sqr(gp)*Sqr(QHu) - 2*MassU*Sqr(gp
      )*Sqr(Ql) + 2*MassU*Sqr(gp)*Sqr(Qv))*(TYv.conjugate()*Yv.transpose()) + 2
      *twoLoop*(-3*traceYuAdjYu - traceYvAdjYv - AbsSqr(Lambdax) + 2*Sqr(gp)*
      Sqr(QHu) - 2*Sqr(gp)*Sqr(Ql) + 2*Sqr(gp)*Sqr(Qv))*(TYv.conjugate()*(TYv).
      transpose()) + 0.2*twoLoop*(-15*traceYdAdjYd - 5*traceYeAdjYe - 5*AbsSqr(
      Lambdax) + 6*Sqr(g1) + 10*Sqr(gp)*Sqr(Qe) + 10*Sqr(gp)*Sqr(QHd) - 10*Sqr(
      gp)*Sqr(Ql))*(ml2*Ye.adjoint()*Ye) + twoLoop*(-3*traceYuAdjYu -
      traceYvAdjYv - AbsSqr(Lambdax) + 2*Sqr(gp)*Sqr(QHu) - 2*Sqr(gp)*Sqr(Ql) +
      2*Sqr(gp)*Sqr(Qv))*(ml2*Yv.conjugate()*Yv.transpose()) + 0.4*twoLoop*(-15
      *traceYdAdjYd - 5*traceYeAdjYe - 5*AbsSqr(Lambdax) + 6*Sqr(g1) + 10*Sqr(
      gp)*Sqr(Qe) + 10*Sqr(gp)*Sqr(QHd) - 10*Sqr(gp)*Sqr(Ql))*(Ye.adjoint()*me2
      *Ye) + 0.2*twoLoop*(-15*traceYdAdjYd - 5*traceYeAdjYe - 5*AbsSqr(Lambdax)
      + 6*Sqr(g1) + 10*Sqr(gp)*Sqr(Qe) + 10*Sqr(gp)*Sqr(QHd) - 10*Sqr(gp)*Sqr(
      Ql))*(Ye.adjoint()*Ye*ml2))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_ml2_2 = (UNITMATRIX(3)*(-2*twoLoop*
      AbsSqr(TLambdax)*(Ye.adjoint()*Ye) - 2*twoLoop*Conj(Lambdax)*TLambdax*((
      TYe).adjoint()*Ye) - 2*twoLoop*AbsSqr(TLambdax)*(Yv.conjugate()*Yv.
      transpose()) - 2*twoLoop*Conj(Lambdax)*TLambdax*(TYv.conjugate()*Yv.
      transpose()) + 2*twoLoop*(-3*traceYuAdjYu - traceYvAdjYv - AbsSqr(Lambdax
      ) + 2*Sqr(gp)*Sqr(QHu) - 2*Sqr(gp)*Sqr(Ql) + 2*Sqr(gp)*Sqr(Qv))*(Yv.
      conjugate()*mvR2*Yv.transpose()) + twoLoop*(-3*traceYuAdjYu -
      traceYvAdjYv - AbsSqr(Lambdax) + 2*Sqr(gp)*Sqr(QHu) - 2*Sqr(gp)*Sqr(Ql) +
      2*Sqr(gp)*Sqr(Qv))*(Yv.conjugate()*Yv.transpose()*ml2) - 8*mHd2*twoLoop*(
      Ye.adjoint()*Ye*Ye.adjoint()*Ye) - 4*twoLoop*(Ye.adjoint()*Ye*(TYe).
      adjoint()*TYe) - 4*twoLoop*(Ye.adjoint()*TYe*(TYe).adjoint()*Ye) - 4*
      twoLoop*((TYe).adjoint()*Ye*Ye.adjoint()*TYe) - 4*twoLoop*((TYe).adjoint(
      )*TYe*Ye.adjoint()*Ye) - 8*mHu2*twoLoop*(Yv.conjugate()*Yv.transpose()*Yv
      .conjugate()*Yv.transpose()) - 4*twoLoop*(Yv.conjugate()*Yv.transpose()*
      TYv.conjugate()*(TYv).transpose()) - 4*twoLoop*(Yv.conjugate()*(TYv).
      transpose()*TYv.conjugate()*Yv.transpose()) - 4*twoLoop*(TYv.conjugate()*
      Yv.transpose()*Yv.conjugate()*(TYv).transpose()) - 4*twoLoop*(TYv.
      conjugate()*(TYv).transpose()*Yv.conjugate()*Yv.transpose()) - 2*twoLoop*
      (ml2*Ye.adjoint()*Ye*Ye.adjoint()*Ye) - 2*twoLoop*(ml2*Yv.conjugate()*Yv.
      transpose()*Yv.conjugate()*Yv.transpose()) - 4*twoLoop*(Ye.adjoint()*me2*
      Ye*Ye.adjoint()*Ye) - 4*twoLoop*(Ye.adjoint()*Ye*ml2*Ye.adjoint()*Ye) - 4
      *twoLoop*(Ye.adjoint()*Ye*Ye.adjoint()*me2*Ye) - 2*twoLoop*(Ye.adjoint()*
      Ye*Ye.adjoint()*Ye*ml2) - 4*twoLoop*(Yv.conjugate()*mvR2*Yv.transpose()*
      Yv.conjugate()*Yv.transpose()) - 4*twoLoop*(Yv.conjugate()*Yv.transpose()
      *ml2*Yv.conjugate()*Yv.transpose()) - 4*twoLoop*(Yv.conjugate()*Yv.
      transpose()*Yv.conjugate()*mvR2*Yv.transpose()) - 2*twoLoop*(Yv.conjugate
      ()*Yv.transpose()*Yv.conjugate()*Yv.transpose()*ml2) + 0.04*twoLoop*(-
      77.45966692414834*g1*gp*Ql*Tr2U114 - 77.45966692414834*g1*gp*Ql*Tr2U141 -
      77.45966692414834*g1*Tr31 + 200*gp*Ql*Tr34 + 621*AbsSqr(MassB)*Quad(g1) +
      150*Tr22*Quad(g2) + 825*AbsSqr(MassWB)*Quad(g2) + 4800*AbsSqr(MassU)*Quad
      (gp)*Quad(Ql) + 30*Tr2U111*Sqr(g1) + 90*AbsSqr(MassB)*Sqr(g1)*Sqr(g2) +
      90*AbsSqr(MassWB)*Sqr(g1)*Sqr(g2) + 45*MassWB*Conj(MassB)*Sqr(g1)*Sqr(g2)
      + 45*MassB*Conj(MassWB)*Sqr(g1)*Sqr(g2) - 360*Qd*Ql*AbsSqr(MassB)*Sqr(g1)
      *Sqr(gp) - 360*Qe*Ql*AbsSqr(MassB)*Sqr(g1)*Sqr(gp) + 120*QHd*Ql*AbsSqr(
      MassB)*Sqr(g1)*Sqr(gp) - 120*QHu*Ql*AbsSqr(MassB)*Sqr(g1)*Sqr(gp) - 360*
      Ql*Qq*AbsSqr(MassB)*Sqr(g1)*Sqr(gp) + 720*Ql*Qu*AbsSqr(MassB)*Sqr(g1)*Sqr
      (gp) - 360*Qd*Ql*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) - 360*Qe*Ql*AbsSqr(MassU)*
      Sqr(g1)*Sqr(gp) + 120*QHd*Ql*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) - 120*QHu*Ql*
      AbsSqr(MassU)*Sqr(g1)*Sqr(gp) - 360*Ql*Qq*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) +
      720*Ql*Qu*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) - 180*MassU*Qd*Ql*Conj(MassB)*Sqr
      (g1)*Sqr(gp) - 180*MassU*Qe*Ql*Conj(MassB)*Sqr(g1)*Sqr(gp) + 60*MassU*QHd
      *Ql*Conj(MassB)*Sqr(g1)*Sqr(gp) - 60*MassU*QHu*Ql*Conj(MassB)*Sqr(g1)*Sqr
      (gp) - 180*MassU*Ql*Qq*Conj(MassB)*Sqr(g1)*Sqr(gp) + 360*MassU*Ql*Qu*Conj
      (MassB)*Sqr(g1)*Sqr(gp) - 180*MassB*Qd*Ql*Conj(MassU)*Sqr(g1)*Sqr(gp) -
      180*MassB*Qe*Ql*Conj(MassU)*Sqr(g1)*Sqr(gp) + 60*MassB*QHd*Ql*Conj(MassU)
      *Sqr(g1)*Sqr(gp) - 60*MassB*QHu*Ql*Conj(MassU)*Sqr(g1)*Sqr(gp) - 180*
      MassB*Ql*Qq*Conj(MassU)*Sqr(g1)*Sqr(gp) + 360*MassB*Ql*Qu*Conj(MassU)*Sqr
      (g1)*Sqr(gp) + 200*Tr2U144*Sqr(gp)*Sqr(Ql) + 480*AbsSqr(MassB)*Sqr(g1)*
      Sqr(gp)*Sqr(Ql) + 480*AbsSqr(MassU)*Sqr(g1)*Sqr(gp)*Sqr(Ql) + 240*MassU*
      Conj(MassB)*Sqr(g1)*Sqr(gp)*Sqr(Ql) + 240*MassB*Conj(MassU)*Sqr(g1)*Sqr(
      gp)*Sqr(Ql) + 600*AbsSqr(MassU)*Sqr(g2)*Sqr(gp)*Sqr(Ql) + 600*AbsSqr(
      MassWB)*Sqr(g2)*Sqr(gp)*Sqr(Ql) + 300*MassWB*Conj(MassU)*Sqr(g2)*Sqr(gp)*
      Sqr(Ql) + 300*MassU*Conj(MassWB)*Sqr(g2)*Sqr(gp)*Sqr(Ql) + 5400*AbsSqr(
      MassU)*Quad(gp)*Sqr(Qd)*Sqr(Ql) + 1800*AbsSqr(MassU)*Quad(gp)*Sqr(Qe)*Sqr
      (Ql) + 1200*AbsSqr(MassU)*Quad(gp)*Sqr(QHd)*Sqr(Ql) + 1200*AbsSqr(MassU)*
      Quad(gp)*Sqr(QHu)*Sqr(Ql) + 10800*AbsSqr(MassU)*Quad(gp)*Sqr(Ql)*Sqr(Qq)
      + 600*AbsSqr(MassU)*Quad(gp)*Sqr(Ql)*Sqr(Qs) + 5400*AbsSqr(MassU)*Quad(gp
      )*Sqr(Ql)*Sqr(Qu) + 1800*AbsSqr(MassU)*Quad(gp)*Sqr(Ql)*Sqr(Qv))*
      UNITMATRIX(3))).real();

   beta_ml2 = beta_ml2_1 + beta_ml2_2;


   return beta_ml2;
}

/**
 * Calculates the 3-loop beta function of ml2.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_ml2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_ml2;

   beta_ml2 = ZEROMATRIX(3,3);


   return beta_ml2;
}

/**
 * Calculates the 4-loop beta function of ml2.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_ml2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_ml2;

   beta_ml2 = ZEROMATRIX(3,3);


   return beta_ml2;
}

/**
 * Calculates the 5-loop beta function of ml2.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_ml2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_ml2;

   beta_ml2 = ZEROMATRIX(3,3);


   return beta_ml2;
}

} // namespace flexiblesusy
