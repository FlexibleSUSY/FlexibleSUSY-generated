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

// File generated at Tue 22 Jan 2019 17:29:31

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
 * Calculates the 1-loop beta function of MassB.
 *
 * @return 1-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_MassB_1_loop(const Soft_traces& soft_traces) const
{


   double beta_MassB;

   beta_MassB = Re(13.2*MassB*oneOver16PiSqr*Sqr(g1));


   return beta_MassB;
}

/**
 * Calculates the 2-loop beta function of MassB.
 *
 * @return 2-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_MassB_2_loop(const Soft_traces& soft_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qu = INPUT(Qu);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;


   double beta_MassB;

   beta_MassB = Re(0.08*twoLoop*Sqr(g1)*(70*traceAdjYdTYd + 90*traceAdjYeTYe +
      130*traceAdjYuTYu + 30*traceAdjYvTYv - 70*MassB*traceYdAdjYd - 90*MassB*
      traceYeAdjYe - 130*MassB*traceYuAdjYu - 30*MassB*traceYvAdjYv - 30*MassB*
      AbsSqr(Lambdax) + 398*MassB*Sqr(g1) + 135*MassB*Sqr(g2) + 135*MassWB*Sqr(
      g2) + 440*MassB*Sqr(g3) + 440*MassG*Sqr(g3) + 60*MassB*Sqr(gp)*Sqr(Qd) +
      60*MassU*Sqr(gp)*Sqr(Qd) + 180*MassB*Sqr(gp)*Sqr(Qe) + 180*MassU*Sqr(gp)*
      Sqr(Qe) + 30*MassB*Sqr(gp)*Sqr(QHd) + 30*MassU*Sqr(gp)*Sqr(QHd) + 30*
      MassB*Sqr(gp)*Sqr(QHu) + 30*MassU*Sqr(gp)*Sqr(QHu) + 90*MassB*Sqr(gp)*Sqr
      (Ql) + 90*MassU*Sqr(gp)*Sqr(Ql) + 30*MassB*Sqr(gp)*Sqr(Qq) + 30*MassU*Sqr
      (gp)*Sqr(Qq) + 240*MassB*Sqr(gp)*Sqr(Qu) + 240*MassU*Sqr(gp)*Sqr(Qu) + 30
      *Conj(Lambdax)*TLambdax));


   return beta_MassB;
}

/**
 * Calculates the 3-loop beta function of MassB.
 *
 * @return 3-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_MassB_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MassB;

   beta_MassB = 0;


   return beta_MassB;
}

/**
 * Calculates the 4-loop beta function of MassB.
 *
 * @return 4-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_MassB_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MassB;

   beta_MassB = 0;


   return beta_MassB;
}

/**
 * Calculates the 5-loop beta function of MassB.
 *
 * @return 5-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_MassB_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MassB;

   beta_MassB = 0;


   return beta_MassB;
}

} // namespace flexiblesusy
