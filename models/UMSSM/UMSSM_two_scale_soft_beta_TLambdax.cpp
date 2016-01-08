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

// File generated at Fri 8 Jan 2016 12:29:36

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
 * Calculates the one-loop beta function of TLambdax.
 *
 * @return one-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_TLambdax_one_loop(const Soft_traces& soft_traces) const
{
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Qs = INPUT(Qs);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;


   double beta_TLambdax;

   beta_TLambdax = Re(oneOver16PiSqr*(0.4*Lambdax*(15*traceAdjYdTYd + 5*
      traceAdjYeTYe + 15*traceAdjYuTYu + 5*traceAdjYvTYv + 3*MassB*Sqr(g1) + 15
      *MassWB*Sqr(g2) + 10*MassU*Sqr(gp)*Sqr(QHd) + 10*MassU*Sqr(gp)*Sqr(QHu) +
      10*MassU*Sqr(gp)*Sqr(Qs)) + (3*traceYdAdjYd + traceYeAdjYe + 3*
      traceYuAdjYu + traceYvAdjYv + 12*AbsSqr(Lambdax) - 0.6*Sqr(g1) - 3*Sqr(g2
      ) - 2*Sqr(gp)*Sqr(QHd) - 2*Sqr(gp)*Sqr(QHu) - 2*Sqr(gp)*Sqr(Qs))*TLambdax
      ));


   return beta_TLambdax;
}

/**
 * Calculates the two-loop beta function of TLambdax.
 *
 * @return two-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_TLambdax_two_loop(const Soft_traces& soft_traces) const
{
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qs = INPUT(Qs);
   const auto Qu = INPUT(Qu);
   const auto Qv = INPUT(Qv);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYvAdjYvYvAdjYv = TRACE_STRUCT.traceYvAdjYvYvAdjYv;
   const double traceYvAdjYvTpYeconjYe =
      TRACE_STRUCT.traceYvAdjYvTpYeconjYe;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYuAdjYuTYuAdjYu = TRACE_STRUCT.traceYuAdjYuTYuAdjYu;
   const double traceYvAdjYvTYvAdjYv = TRACE_STRUCT.traceYvAdjYvTYvAdjYv;
   const double traceAdjYeTYeconjYvTpYv =
      TRACE_STRUCT.traceAdjYeTYeconjYvTpYv;
   const double traceAdjYvTpYeconjYeTYv =
      TRACE_STRUCT.traceAdjYvTpYeconjYeTYv;


   double beta_TLambdax;

   const double beta_TLambdax_1 = Re(-0.08*twoLoop*Lambdax*(261*Power(g1,
      4)*MassB + 5*Sqr(g1)*(9*(MassB + MassWB)*Sqr(g2) - 2*(-traceAdjYdTYd + 3*
      traceAdjYeTYe + 2*traceAdjYuTYu + 3*traceAdjYvTYv + MassB*traceYdAdjYd -
      3*MassB*traceYeAdjYe - 2*MassB*traceYuAdjYu + 3*(MassB + MassU)*Sqr(gp)*(
      3*Qd*(QHd - QHu) + 3*Qe*(QHd - QHu) + 2*QHd*QHu - 3*QHd*Ql + 3*QHu*Ql + 3
      *QHd*Qq - 3*QHu*Qq - 6*QHd*Qu + 6*QHu*Qu + 3*QHd*Qv - 3*QHu*Qv - 2*Sqr(
      QHd) - 2*Sqr(QHu)))) + 25*(15*Power(g2,4)*MassWB + 6*(MassU + MassWB)*Sqr
      (g2)*Sqr(gp)*(Sqr(QHd) + Sqr(QHu)) + 2*(traceAdjYeTYeconjYvTpYv +
      traceAdjYvTpYeconjYeTYv + 9*traceYdAdjYdTYdAdjYd + 3*traceYdAdjYuTYuAdjYd
      + 3*traceYeAdjYeTYeAdjYe + 3*traceYuAdjYdTYdAdjYu - 8*traceAdjYdTYd*Sqr(
      g3) - 8*traceAdjYuTYu*Sqr(g3) + 8*MassG*traceYdAdjYd*Sqr(g3) + 8*MassG*
      traceYuAdjYu*Sqr(g3) - Sqr(gp)*(3*(traceAdjYdTYd - MassU*traceYdAdjYd)*
      Sqr(Qd) + traceAdjYeTYe*Sqr(Qe) - MassU*traceYeAdjYe*Sqr(Qe) + (-3*
      traceAdjYdTYd - traceAdjYeTYe + 3*MassU*traceYdAdjYd + MassU*traceYeAdjYe
      )*Sqr(QHd) - 3*traceAdjYuTYu*Sqr(QHu) - traceAdjYvTYv*Sqr(QHu) + 3*MassU*
      traceYuAdjYu*Sqr(QHu) + traceAdjYeTYe*Sqr(Ql) + traceAdjYvTYv*Sqr(Ql) -
      MassU*traceYeAdjYe*Sqr(Ql) + 3*traceAdjYdTYd*Sqr(Qq) + 3*traceAdjYuTYu*
      Sqr(Qq) - 3*MassU*traceYdAdjYd*Sqr(Qq) - 3*MassU*traceYuAdjYu*Sqr(Qq) + 3
      *traceAdjYuTYu*Sqr(Qu) - 3*MassU*traceYuAdjYu*Sqr(Qu) + traceAdjYvTYv*Sqr
      (Qv)) + 2*Power(gp,4)*MassU*(4*Power(QHd,4) + 4*Power(QHu,4) + 3*Power(Qs
      ,4) + 4*Sqr(QHd)*Sqr(QHu) + 6*Sqr(QHd)*Sqr(Ql) + 6*Sqr(QHu)*Sqr(Ql) + 18*
      Sqr(QHd)*Sqr(Qq) + 18*Sqr(QHu)*Sqr(Qq) + 3*Sqr(QHd)*Sqr(Qs) + 3*Sqr(QHu)*
      Sqr(Qs) + 6*Sqr(Ql)*Sqr(Qs) + 18*Sqr(Qq)*Sqr(Qs) + 9*Sqr(Qd)*(Sqr(QHd) +
      Sqr(QHu) + Sqr(Qs)) + 3*Sqr(Qe)*(Sqr(QHd) + Sqr(QHu) + Sqr(Qs)) + 9*Sqr(
      QHd)*Sqr(Qu) + 9*Sqr(QHu)*Sqr(Qu) + 9*Sqr(Qs)*Sqr(Qu) + 3*Sqr(QHd)*Sqr(Qv
      ) + 3*Sqr(QHu)*Sqr(Qv) + 3*Sqr(Qs)*Sqr(Qv))))));
   const double beta_TLambdax_2 = Re(-0.02*twoLoop*(40*Lambdax*(45*
      traceYuAdjYuTYuAdjYu + 3*MassB*traceYvAdjYv*Sqr(g1) + 5*(3*
      traceYvAdjYvTYvAdjYv - MassU*traceYvAdjYv*Sqr(gp)*Sqr(QHu) + MassU*
      traceYvAdjYv*Sqr(gp)*Sqr(Ql) + MassU*traceYvAdjYv*Sqr(gp)*Sqr(Qv))) - (
      261*Power(g1,4) + 10*Sqr(g1)*(9*Sqr(g2) - 2*(traceYdAdjYd - 3*
      traceYeAdjYe - 2*traceYuAdjYu - 3*traceYvAdjYv + 3*Sqr(gp)*(3*Qd*QHd + 3*
      Qe*QHd - 3*Qd*QHu - 3*Qe*QHu + 2*QHd*QHu - 3*QHd*Ql + 3*QHu*Ql + 3*QHd*Qq
      - 3*QHu*Qq - 6*QHd*Qu + 6*QHu*Qu + 3*QHd*Qv - 3*QHu*Qv - 2*Sqr(QHd) - 2*
      Sqr(QHu)))) + 25*(15*Power(g2,4) + 12*Sqr(g2)*Sqr(gp)*(Sqr(QHd) + Sqr(QHu
      )) + 2*(-9*traceYdAdjYdYdAdjYd - 6*traceYdAdjYuYuAdjYd - 3*
      traceYeAdjYeYeAdjYe - 9*traceYuAdjYuYuAdjYu - 2*traceYvAdjYvTpYeconjYe -
      3*traceYvAdjYvYvAdjYv + 16*traceYdAdjYd*Sqr(g3) + 16*traceYuAdjYu*Sqr(g3)
      + 2*Sqr(gp)*(3*traceYdAdjYd*Sqr(Qd) + traceYeAdjYe*Sqr(Qe) - (3*
      traceYdAdjYd + traceYeAdjYe)*Sqr(QHd) - 3*traceYuAdjYu*Sqr(QHu) -
      traceYvAdjYv*Sqr(QHu) + traceYeAdjYe*Sqr(Ql) + traceYvAdjYv*Sqr(Ql) + 3*
      traceYdAdjYd*Sqr(Qq) + 3*traceYuAdjYu*Sqr(Qq) + 3*traceYuAdjYu*Sqr(Qu) +
      traceYvAdjYv*Sqr(Qv)) + 2*Power(gp,4)*(4*Power(QHd,4) + 4*Power(QHu,4) +
      3*Power(Qs,4) + 4*Sqr(QHd)*Sqr(QHu) + 6*Sqr(QHd)*Sqr(Ql) + 6*Sqr(QHu)*Sqr
      (Ql) + 18*Sqr(QHd)*Sqr(Qq) + 18*Sqr(QHu)*Sqr(Qq) + 3*Sqr(QHd)*Sqr(Qs) + 3
      *Sqr(QHu)*Sqr(Qs) + 6*Sqr(Ql)*Sqr(Qs) + 18*Sqr(Qq)*Sqr(Qs) + 9*Sqr(Qd)*(
      Sqr(QHd) + Sqr(QHu) + Sqr(Qs)) + 3*Sqr(Qe)*(Sqr(QHd) + Sqr(QHu) + Sqr(Qs)
      ) + 9*Sqr(QHd)*Sqr(Qu) + 9*Sqr(QHu)*Sqr(Qu) + 9*Sqr(Qs)*Sqr(Qu) + 3*Sqr(
      QHd)*Sqr(Qv) + 3*Sqr(QHu)*Sqr(Qv) + 3*Sqr(Qs)*Sqr(Qv)))))*TLambdax + 2500
      *Sqr(Conj(Lambdax))*Sqr(Lambdax)*TLambdax + 10*AbsSqr(Lambdax)*(2*Lambdax
      *(6*MassB*Sqr(g1) + 5*(3*(3*traceAdjYdTYd + traceAdjYeTYe + 3*
      traceAdjYuTYu + traceAdjYvTYv) + 6*MassWB*Sqr(g2) + 4*MassU*Sqr(gp)*(Sqr(
      QHd) + Sqr(QHu)))) - 3*(6*Sqr(g1) + 5*(-3*(3*traceYdAdjYd + traceYeAdjYe
      + 3*traceYuAdjYu + traceYvAdjYv) + 6*Sqr(g2) + 4*Sqr(gp)*(Sqr(QHd) + Sqr(
      QHu))))*TLambdax)));

   beta_TLambdax = beta_TLambdax_1 + beta_TLambdax_2;


   return beta_TLambdax;
}

/**
 * Calculates the three-loop beta function of TLambdax.
 *
 * @return three-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_TLambdax_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_TLambdax;

   beta_TLambdax = 0;


   return beta_TLambdax;
}

} // namespace flexiblesusy
