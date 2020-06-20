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
 * Calculates the 1-loop beta function of ms2.
 *
 * @return 1-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_ms2_1_loop(const Soft_traces& soft_traces) const
{
   const auto Qs = INPUT(Qs);
   const double Tr14 = TRACE_STRUCT.Tr14;


   double beta_ms2;

   beta_ms2 = Re(-2*(-(gp*Qs*Tr14) - 2*mHd2*AbsSqr(Lambdax) - 2*mHu2*AbsSqr(
      Lambdax) - 2*ms2*AbsSqr(Lambdax) - 2*AbsSqr(TLambdax) + 4*AbsSqr(MassU)*
      Sqr(gp)*Sqr(Qs)));


   return oneLoop * beta_ms2;
}

/**
 * Calculates the 2-loop beta function of ms2.
 *
 * @return 2-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_ms2_2_loop(const Soft_traces& soft_traces) const
{
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Qs = INPUT(Qs);
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
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
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double traceconjTYvTpYv = TRACE_STRUCT.traceconjTYvTpYv;
   const double traceconjTYvTpTYv = TRACE_STRUCT.traceconjTYvTpTYv;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceYvAdjYvconjml2 = TRACE_STRUCT.traceYvAdjYvconjml2;
   const double traceYvconjmvR2AdjYv = TRACE_STRUCT.traceYvconjmvR2AdjYv;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_ms2;

   beta_ms2 = Re(0.8*(10*gp*Qs*Tr34 - 15*traceconjTYdTpTYd*AbsSqr(Lambdax) - 5*
      traceconjTYeTpTYe*AbsSqr(Lambdax) - 15*traceconjTYuTpTYu*AbsSqr(Lambdax)
      - 5*traceconjTYvTpTYv*AbsSqr(Lambdax) - 15*tracemd2YdAdjYd*AbsSqr(Lambdax
      ) - 5*traceme2YeAdjYe*AbsSqr(Lambdax) - 5*traceml2AdjYeYe*AbsSqr(Lambdax)
      - 15*tracemq2AdjYdYd*AbsSqr(Lambdax) - 15*tracemq2AdjYuYu*AbsSqr(Lambdax)
      - 15*tracemu2YuAdjYu*AbsSqr(Lambdax) - 30*mHd2*traceYdAdjYd*AbsSqr(
      Lambdax) - 15*mHu2*traceYdAdjYd*AbsSqr(Lambdax) - 15*ms2*traceYdAdjYd*
      AbsSqr(Lambdax) - 10*mHd2*traceYeAdjYe*AbsSqr(Lambdax) - 5*mHu2*
      traceYeAdjYe*AbsSqr(Lambdax) - 5*ms2*traceYeAdjYe*AbsSqr(Lambdax) - 15*
      mHd2*traceYuAdjYu*AbsSqr(Lambdax) - 30*mHu2*traceYuAdjYu*AbsSqr(Lambdax)
      - 15*ms2*traceYuAdjYu*AbsSqr(Lambdax) - 5*mHd2*traceYvAdjYv*AbsSqr(
      Lambdax) - 10*mHu2*traceYvAdjYv*AbsSqr(Lambdax) - 5*ms2*traceYvAdjYv*
      AbsSqr(Lambdax) - 5*traceYvAdjYvconjml2*AbsSqr(Lambdax) - 5*
      traceYvconjmvR2AdjYv*AbsSqr(Lambdax) - 15*traceYdAdjYd*AbsSqr(TLambdax) -
      5*traceYeAdjYe*AbsSqr(TLambdax) - 15*traceYuAdjYu*AbsSqr(TLambdax) - 5*
      traceYvAdjYv*AbsSqr(TLambdax) - 40*AbsSqr(Lambdax)*AbsSqr(TLambdax) - 15*
      traceAdjYdTYd*Conj(TLambdax)*Lambdax - 5*traceAdjYeTYe*Conj(TLambdax)*
      Lambdax - 15*traceAdjYuTYu*Conj(TLambdax)*Lambdax - 5*traceAdjYvTYv*Conj(
      TLambdax)*Lambdax + 90*AbsSqr(MassU)*Quad(gp)*Quad(Qs) + 3*mHd2*AbsSqr(
      Lambdax)*Sqr(g1) + 3*mHu2*AbsSqr(Lambdax)*Sqr(g1) + 3*ms2*AbsSqr(Lambdax)
      *Sqr(g1) + 6*AbsSqr(MassB)*AbsSqr(Lambdax)*Sqr(g1) + 3*AbsSqr(TLambdax)*
      Sqr(g1) - 3*MassB*Conj(TLambdax)*Lambdax*Sqr(g1) + 15*mHd2*AbsSqr(Lambdax
      )*Sqr(g2) + 15*mHu2*AbsSqr(Lambdax)*Sqr(g2) + 15*ms2*AbsSqr(Lambdax)*Sqr(
      g2) + 30*AbsSqr(MassWB)*AbsSqr(Lambdax)*Sqr(g2) + 15*AbsSqr(TLambdax)*Sqr
      (g2) - 15*MassWB*Conj(TLambdax)*Lambdax*Sqr(g2) + 10*mHd2*AbsSqr(Lambdax)
      *Sqr(gp)*Sqr(QHd) + 10*mHu2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHd) + 10*ms2*
      AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHd) + 20*AbsSqr(MassU)*AbsSqr(Lambdax)*Sqr(
      gp)*Sqr(QHd) + 10*AbsSqr(TLambdax)*Sqr(gp)*Sqr(QHd) - 10*MassU*Conj(
      TLambdax)*Lambdax*Sqr(gp)*Sqr(QHd) + 10*mHd2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(
      QHu) + 10*mHu2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHu) + 10*ms2*AbsSqr(Lambdax)*
      Sqr(gp)*Sqr(QHu) + 20*AbsSqr(MassU)*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHu) + 10
      *AbsSqr(TLambdax)*Sqr(gp)*Sqr(QHu) - 10*MassU*Conj(TLambdax)*Lambdax*Sqr(
      gp)*Sqr(QHu) + 10*Tr2U144*Sqr(gp)*Sqr(Qs) - 10*mHd2*AbsSqr(Lambdax)*Sqr(
      gp)*Sqr(Qs) - 10*mHu2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs) - 10*ms2*AbsSqr(
      Lambdax)*Sqr(gp)*Sqr(Qs) - 20*AbsSqr(MassU)*AbsSqr(Lambdax)*Sqr(gp)*Sqr(
      Qs) - 10*AbsSqr(TLambdax)*Sqr(gp)*Sqr(Qs) + 10*MassU*Conj(TLambdax)*
      Lambdax*Sqr(gp)*Sqr(Qs) + 270*AbsSqr(MassU)*Quad(gp)*Sqr(Qd)*Sqr(Qs) + 90
      *AbsSqr(MassU)*Quad(gp)*Sqr(Qe)*Sqr(Qs) + 60*AbsSqr(MassU)*Quad(gp)*Sqr(
      QHd)*Sqr(Qs) + 60*AbsSqr(MassU)*Quad(gp)*Sqr(QHu)*Sqr(Qs) + 180*AbsSqr(
      MassU)*Quad(gp)*Sqr(Ql)*Sqr(Qs) + 540*AbsSqr(MassU)*Quad(gp)*Sqr(Qq)*Sqr(
      Qs) + 270*AbsSqr(MassU)*Quad(gp)*Sqr(Qs)*Sqr(Qu) + 90*AbsSqr(MassU)*Quad(
      gp)*Sqr(Qs)*Sqr(Qv) - 20*mHd2*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 20*mHu2*
      Sqr(Conj(Lambdax))*Sqr(Lambdax) - 20*ms2*Sqr(Conj(Lambdax))*Sqr(Lambdax)
      - 15*traceconjTYdTpYd*Conj(Lambdax)*TLambdax - 5*traceconjTYeTpYe*Conj(
      Lambdax)*TLambdax - 15*traceconjTYuTpYu*Conj(Lambdax)*TLambdax - 5*
      traceconjTYvTpYv*Conj(Lambdax)*TLambdax - 3*Conj(MassB)*Conj(Lambdax)*Sqr
      (g1)*TLambdax - 15*Conj(MassWB)*Conj(Lambdax)*Sqr(g2)*TLambdax - 10*Conj(
      MassU)*Conj(Lambdax)*Sqr(gp)*Sqr(QHd)*TLambdax - 10*Conj(MassU)*Conj(
      Lambdax)*Sqr(gp)*Sqr(QHu)*TLambdax + 10*Conj(MassU)*Conj(Lambdax)*Sqr(gp)
      *Sqr(Qs)*TLambdax));


   return twoLoop * beta_ms2;
}

/**
 * Calculates the 3-loop beta function of ms2.
 *
 * @return 3-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_ms2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_ms2;

   beta_ms2 = 0;


   return threeLoop * beta_ms2;
}

/**
 * Calculates the 4-loop beta function of ms2.
 *
 * @return 4-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_ms2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_ms2;

   beta_ms2 = 0;


   return fourLoop * beta_ms2;
}

/**
 * Calculates the 5-loop beta function of ms2.
 *
 * @return 5-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_ms2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_ms2;

   beta_ms2 = 0;


   return fiveLoop * beta_ms2;
}

} // namespace flexiblesusy
