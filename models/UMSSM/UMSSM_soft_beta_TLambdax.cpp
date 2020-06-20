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
 * Calculates the 1-loop beta function of TLambdax.
 *
 * @return 1-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_TLambdax_1_loop(const Soft_traces& soft_traces) const
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

   beta_TLambdax = Re(0.2*(30*traceAdjYdTYd*Lambdax + 10*traceAdjYeTYe*Lambdax
      + 30*traceAdjYuTYu*Lambdax + 10*traceAdjYvTYv*Lambdax + 6*MassB*Lambdax*
      Sqr(g1) + 30*MassWB*Lambdax*Sqr(g2) + 20*MassU*Lambdax*Sqr(gp)*Sqr(QHd) +
      20*MassU*Lambdax*Sqr(gp)*Sqr(QHu) + 20*MassU*Lambdax*Sqr(gp)*Sqr(Qs) + 15
      *traceYdAdjYd*TLambdax + 5*traceYeAdjYe*TLambdax + 15*traceYuAdjYu*
      TLambdax + 5*traceYvAdjYv*TLambdax + 60*AbsSqr(Lambdax)*TLambdax - 3*Sqr(
      g1)*TLambdax - 15*Sqr(g2)*TLambdax - 10*Sqr(gp)*Sqr(QHd)*TLambdax - 10*
      Sqr(gp)*Sqr(QHu)*TLambdax - 10*Sqr(gp)*Sqr(Qs)*TLambdax));


   return oneLoop * beta_TLambdax;
}

/**
 * Calculates the 2-loop beta function of TLambdax.
 *
 * @return 2-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_TLambdax_2_loop(const Soft_traces& soft_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto QHd = INPUT(QHd);
   const auto Qe = INPUT(Qe);
   const auto QHu = INPUT(QHu);
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
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYuAdjYuTYuAdjYu = TRACE_STRUCT.traceYuAdjYuTYuAdjYu;
   const double traceYvAdjYvYvAdjYv = TRACE_STRUCT.traceYvAdjYvYvAdjYv;
   const double traceYvAdjYvTYvAdjYv = TRACE_STRUCT.traceYvAdjYvTYvAdjYv;
   const double traceYvAdjYvTpYeconjYe = TRACE_STRUCT.traceYvAdjYvTpYeconjYe;
   const double traceAdjYeTYeconjYvTpYv = TRACE_STRUCT.traceAdjYeTYeconjYvTpYv;
   const double traceAdjYvTpYeconjYeTYv = TRACE_STRUCT.traceAdjYvTpYeconjYeTYv;


   double beta_TLambdax;

   const double beta_TLambdax_1 = Re(-0.08*Lambdax*(50*traceAdjYeTYeconjYvTpYv
      + 50*traceAdjYvTpYeconjYeTYv + 450*traceYdAdjYdTYdAdjYd + 150*
      traceYdAdjYuTYuAdjYd + 150*traceYeAdjYeTYeAdjYe + 150*
      traceYuAdjYdTYdAdjYu + 450*traceYuAdjYuTYuAdjYu + 150*
      traceYvAdjYvTYvAdjYv + 207*MassB*Quad(g1) + 375*MassWB*Quad(g2) + 400*
      MassU*Quad(gp)*Quad(QHd) + 400*MassU*Quad(gp)*Quad(QHu) + 300*MassU*Quad(
      gp)*Quad(Qs) + 10*traceAdjYdTYd*Sqr(g1) - 30*traceAdjYeTYe*Sqr(g1) - 20*
      traceAdjYuTYu*Sqr(g1) - 10*MassB*traceYdAdjYd*Sqr(g1) + 30*MassB*
      traceYeAdjYe*Sqr(g1) + 20*MassB*traceYuAdjYu*Sqr(g1) + 45*MassB*Sqr(g1)*
      Sqr(g2) + 45*MassWB*Sqr(g1)*Sqr(g2) - 400*traceAdjYdTYd*Sqr(g3) - 400*
      traceAdjYuTYu*Sqr(g3) + 400*MassG*traceYdAdjYd*Sqr(g3) + 400*MassG*
      traceYuAdjYu*Sqr(g3) - 90*MassB*Qd*QHd*Sqr(g1)*Sqr(gp) - 90*MassU*Qd*QHd*
      Sqr(g1)*Sqr(gp) - 90*MassB*Qe*QHd*Sqr(g1)*Sqr(gp) - 90*MassU*Qe*QHd*Sqr(
      g1)*Sqr(gp) + 90*MassB*Qd*QHu*Sqr(g1)*Sqr(gp) + 90*MassU*Qd*QHu*Sqr(g1)*
      Sqr(gp) + 90*MassB*Qe*QHu*Sqr(g1)*Sqr(gp) + 90*MassU*Qe*QHu*Sqr(g1)*Sqr(
      gp) - 60*MassB*QHd*QHu*Sqr(g1)*Sqr(gp) - 60*MassU*QHd*QHu*Sqr(g1)*Sqr(gp)
      + 90*MassB*QHd*Ql*Sqr(g1)*Sqr(gp) + 90*MassU*QHd*Ql*Sqr(g1)*Sqr(gp) - 90*
      MassB*QHu*Ql*Sqr(g1)*Sqr(gp) - 90*MassU*QHu*Ql*Sqr(g1)*Sqr(gp) - 90*MassB
      *QHd*Qq*Sqr(g1)*Sqr(gp) - 90*MassU*QHd*Qq*Sqr(g1)*Sqr(gp) + 90*MassB*QHu*
      Qq*Sqr(g1)*Sqr(gp) + 90*MassU*QHu*Qq*Sqr(g1)*Sqr(gp) + 180*MassB*QHd*Qu*
      Sqr(g1)*Sqr(gp) + 180*MassU*QHd*Qu*Sqr(g1)*Sqr(gp) - 180*MassB*QHu*Qu*Sqr
      (g1)*Sqr(gp) - 180*MassU*QHu*Qu*Sqr(g1)*Sqr(gp) - 150*traceAdjYdTYd*Sqr(
      gp)*Sqr(Qd) + 150*MassU*traceYdAdjYd*Sqr(gp)*Sqr(Qd) - 50*traceAdjYeTYe*
      Sqr(gp)*Sqr(Qe) + 50*MassU*traceYeAdjYe*Sqr(gp)*Sqr(Qe) + 150*
      traceAdjYdTYd*Sqr(gp)*Sqr(QHd) + 50*traceAdjYeTYe*Sqr(gp)*Sqr(QHd) - 150*
      MassU*traceYdAdjYd*Sqr(gp)*Sqr(QHd) - 50*MassU*traceYeAdjYe*Sqr(gp)*Sqr(
      QHd) + 60*MassB*Sqr(g1)*Sqr(gp)*Sqr(QHd) + 60*MassU*Sqr(g1)*Sqr(gp)*Sqr(
      QHd) + 150*MassU*Sqr(g2)*Sqr(gp)*Sqr(QHd) + 150*MassWB*Sqr(g2)*Sqr(gp)*
      Sqr(QHd) + 900*MassU*Quad(gp)*Sqr(Qd)*Sqr(QHd) + 300*MassU*Quad(gp)*Sqr(
      Qe)*Sqr(QHd) + 150*traceAdjYuTYu*Sqr(gp)*Sqr(QHu) + 50*traceAdjYvTYv*Sqr(
      gp)*Sqr(QHu) - 150*MassU*traceYuAdjYu*Sqr(gp)*Sqr(QHu) - 50*MassU*
      traceYvAdjYv*Sqr(gp)*Sqr(QHu) + 60*MassB*Sqr(g1)*Sqr(gp)*Sqr(QHu) + 60*
      MassU*Sqr(g1)*Sqr(gp)*Sqr(QHu) + 150*MassU*Sqr(g2)*Sqr(gp)*Sqr(QHu) + 150
      *MassWB*Sqr(g2)*Sqr(gp)*Sqr(QHu) + 900*MassU*Quad(gp)*Sqr(Qd)*Sqr(QHu) +
      300*MassU*Quad(gp)*Sqr(Qe)*Sqr(QHu) + 400*MassU*Quad(gp)*Sqr(QHd)*Sqr(QHu
      ) - 50*traceAdjYeTYe*Sqr(gp)*Sqr(Ql) - 50*traceAdjYvTYv*Sqr(gp)*Sqr(Ql) +
      50*MassU*traceYeAdjYe*Sqr(gp)*Sqr(Ql) + 50*MassU*traceYvAdjYv*Sqr(gp)*Sqr
      (Ql) + 600*MassU*Quad(gp)*Sqr(QHd)*Sqr(Ql) + 600*MassU*Quad(gp)*Sqr(QHu)*
      Sqr(Ql) - 150*traceAdjYdTYd*Sqr(gp)*Sqr(Qq) - 150*traceAdjYuTYu*Sqr(gp)*
      Sqr(Qq) + 150*MassU*traceYdAdjYd*Sqr(gp)*Sqr(Qq) + 150*MassU*traceYuAdjYu
      *Sqr(gp)*Sqr(Qq) + 1800*MassU*Quad(gp)*Sqr(QHd)*Sqr(Qq) + 1800*MassU*Quad
      (gp)*Sqr(QHu)*Sqr(Qq) + 900*MassU*Quad(gp)*Sqr(Qd)*Sqr(Qs) + 300*MassU*
      Quad(gp)*Sqr(Qe)*Sqr(Qs) + 300*MassU*Quad(gp)*Sqr(QHd)*Sqr(Qs) + 300*
      MassU*Quad(gp)*Sqr(QHu)*Sqr(Qs) + 600*MassU*Quad(gp)*Sqr(Ql)*Sqr(Qs) +
      1800*MassU*Quad(gp)*Sqr(Qq)*Sqr(Qs) - 150*traceAdjYuTYu*Sqr(gp)*Sqr(Qu) +
      150*MassU*traceYuAdjYu*Sqr(gp)*Sqr(Qu) + 900*MassU*Quad(gp)*Sqr(QHd)*Sqr(
      Qu) + 900*MassU*Quad(gp)*Sqr(QHu)*Sqr(Qu) + 900*MassU*Quad(gp)*Sqr(Qs)*
      Sqr(Qu) - 50*traceAdjYvTYv*Sqr(gp)*Sqr(Qv) + 50*MassU*traceYvAdjYv*Sqr(gp
      )*Sqr(Qv) + 300*MassU*Quad(gp)*Sqr(QHd)*Sqr(Qv) + 300*MassU*Quad(gp)*Sqr(
      QHu)*Sqr(Qv) + 300*MassU*Quad(gp)*Sqr(Qs)*Sqr(Qv)));
   const double beta_TLambdax_2 = Re(0.02*(-900*traceAdjYdTYd*Conj(Lambdax)*Sqr
      (Lambdax) - 300*traceAdjYeTYe*Conj(Lambdax)*Sqr(Lambdax) - 900*
      traceAdjYuTYu*Conj(Lambdax)*Sqr(Lambdax) - 300*traceAdjYvTYv*Conj(Lambdax
      )*Sqr(Lambdax) - 120*MassB*Conj(Lambdax)*Sqr(g1)*Sqr(Lambdax) - 600*
      MassWB*Conj(Lambdax)*Sqr(g2)*Sqr(Lambdax) - 400*MassU*Conj(Lambdax)*Sqr(
      gp)*Sqr(QHd)*Sqr(Lambdax) - 400*MassU*Conj(Lambdax)*Sqr(gp)*Sqr(QHu)*Sqr(
      Lambdax) - 450*traceYdAdjYdYdAdjYd*TLambdax - 300*traceYdAdjYuYuAdjYd*
      TLambdax - 150*traceYeAdjYeYeAdjYe*TLambdax - 450*traceYuAdjYuYuAdjYu*
      TLambdax - 100*traceYvAdjYvTpYeconjYe*TLambdax - 150*traceYvAdjYvYvAdjYv*
      TLambdax - 1350*traceYdAdjYd*AbsSqr(Lambdax)*TLambdax - 450*traceYeAdjYe*
      AbsSqr(Lambdax)*TLambdax - 1350*traceYuAdjYu*AbsSqr(Lambdax)*TLambdax -
      450*traceYvAdjYv*AbsSqr(Lambdax)*TLambdax + 207*Quad(g1)*TLambdax + 375*
      Quad(g2)*TLambdax + 400*Quad(gp)*Quad(QHd)*TLambdax + 400*Quad(gp)*Quad(
      QHu)*TLambdax + 300*Quad(gp)*Quad(Qs)*TLambdax - 20*traceYdAdjYd*Sqr(g1)*
      TLambdax + 60*traceYeAdjYe*Sqr(g1)*TLambdax + 40*traceYuAdjYu*Sqr(g1)*
      TLambdax + 180*AbsSqr(Lambdax)*Sqr(g1)*TLambdax + 900*AbsSqr(Lambdax)*Sqr
      (g2)*TLambdax + 90*Sqr(g1)*Sqr(g2)*TLambdax + 800*traceYdAdjYd*Sqr(g3)*
      TLambdax + 800*traceYuAdjYu*Sqr(g3)*TLambdax - 180*Qd*QHd*Sqr(g1)*Sqr(gp)
      *TLambdax - 180*Qe*QHd*Sqr(g1)*Sqr(gp)*TLambdax + 180*Qd*QHu*Sqr(g1)*Sqr(
      gp)*TLambdax + 180*Qe*QHu*Sqr(g1)*Sqr(gp)*TLambdax - 120*QHd*QHu*Sqr(g1)*
      Sqr(gp)*TLambdax + 180*QHd*Ql*Sqr(g1)*Sqr(gp)*TLambdax - 180*QHu*Ql*Sqr(
      g1)*Sqr(gp)*TLambdax - 180*QHd*Qq*Sqr(g1)*Sqr(gp)*TLambdax + 180*QHu*Qq*
      Sqr(g1)*Sqr(gp)*TLambdax + 360*QHd*Qu*Sqr(g1)*Sqr(gp)*TLambdax - 360*QHu*
      Qu*Sqr(g1)*Sqr(gp)*TLambdax + 300*traceYdAdjYd*Sqr(gp)*Sqr(Qd)*TLambdax +
      100*traceYeAdjYe*Sqr(gp)*Sqr(Qe)*TLambdax - 300*traceYdAdjYd*Sqr(gp)*Sqr(
      QHd)*TLambdax - 100*traceYeAdjYe*Sqr(gp)*Sqr(QHd)*TLambdax + 600*AbsSqr(
      Lambdax)*Sqr(gp)*Sqr(QHd)*TLambdax + 120*Sqr(g1)*Sqr(gp)*Sqr(QHd)*
      TLambdax + 300*Sqr(g2)*Sqr(gp)*Sqr(QHd)*TLambdax + 900*Quad(gp)*Sqr(Qd)*
      Sqr(QHd)*TLambdax + 300*Quad(gp)*Sqr(Qe)*Sqr(QHd)*TLambdax - 300*
      traceYuAdjYu*Sqr(gp)*Sqr(QHu)*TLambdax - 100*traceYvAdjYv*Sqr(gp)*Sqr(QHu
      )*TLambdax + 600*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHu)*TLambdax + 120*Sqr(g1)*
      Sqr(gp)*Sqr(QHu)*TLambdax + 300*Sqr(g2)*Sqr(gp)*Sqr(QHu)*TLambdax + 900*
      Quad(gp)*Sqr(Qd)*Sqr(QHu)*TLambdax + 300*Quad(gp)*Sqr(Qe)*Sqr(QHu)*
      TLambdax + 400*Quad(gp)*Sqr(QHd)*Sqr(QHu)*TLambdax + 100*traceYeAdjYe*Sqr
      (gp)*Sqr(Ql)*TLambdax + 100*traceYvAdjYv*Sqr(gp)*Sqr(Ql)*TLambdax + 600*
      Quad(gp)*Sqr(QHd)*Sqr(Ql)*TLambdax + 600*Quad(gp)*Sqr(QHu)*Sqr(Ql)*
      TLambdax + 300*traceYdAdjYd*Sqr(gp)*Sqr(Qq)*TLambdax + 300*traceYuAdjYu*
      Sqr(gp)*Sqr(Qq)*TLambdax + 1800*Quad(gp)*Sqr(QHd)*Sqr(Qq)*TLambdax + 1800
      *Quad(gp)*Sqr(QHu)*Sqr(Qq)*TLambdax + 900*Quad(gp)*Sqr(Qd)*Sqr(Qs)*
      TLambdax + 300*Quad(gp)*Sqr(Qe)*Sqr(Qs)*TLambdax + 300*Quad(gp)*Sqr(QHd)*
      Sqr(Qs)*TLambdax + 300*Quad(gp)*Sqr(QHu)*Sqr(Qs)*TLambdax + 600*Quad(gp)*
      Sqr(Ql)*Sqr(Qs)*TLambdax + 1800*Quad(gp)*Sqr(Qq)*Sqr(Qs)*TLambdax + 300*
      traceYuAdjYu*Sqr(gp)*Sqr(Qu)*TLambdax + 900*Quad(gp)*Sqr(QHd)*Sqr(Qu)*
      TLambdax + 900*Quad(gp)*Sqr(QHu)*Sqr(Qu)*TLambdax + 900*Quad(gp)*Sqr(Qs)*
      Sqr(Qu)*TLambdax + 100*traceYvAdjYv*Sqr(gp)*Sqr(Qv)*TLambdax + 300*Quad(
      gp)*Sqr(QHd)*Sqr(Qv)*TLambdax + 300*Quad(gp)*Sqr(QHu)*Sqr(Qv)*TLambdax +
      300*Quad(gp)*Sqr(Qs)*Sqr(Qv)*TLambdax - 2500*Sqr(Conj(Lambdax))*Sqr(
      Lambdax)*TLambdax));

   beta_TLambdax = beta_TLambdax_1 + beta_TLambdax_2;


   return twoLoop * beta_TLambdax;
}

/**
 * Calculates the 3-loop beta function of TLambdax.
 *
 * @return 3-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_TLambdax_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_TLambdax;

   beta_TLambdax = 0;


   return threeLoop * beta_TLambdax;
}

/**
 * Calculates the 4-loop beta function of TLambdax.
 *
 * @return 4-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_TLambdax_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_TLambdax;

   beta_TLambdax = 0;


   return fourLoop * beta_TLambdax;
}

/**
 * Calculates the 5-loop beta function of TLambdax.
 *
 * @return 5-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_TLambdax_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_TLambdax;

   beta_TLambdax = 0;


   return fiveLoop * beta_TLambdax;
}

} // namespace flexiblesusy
