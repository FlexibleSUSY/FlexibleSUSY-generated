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


#include "MRSSM2_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of MDBS.
 *
 * @return 1-loop beta function
 */
double MRSSM2_soft_parameters::calc_beta_MDBS_1_loop(const Soft_traces& soft_traces) const
{


   double beta_MDBS;

   beta_MDBS = Re(0.4*MDBS*(5*AbsSqr(LamSD) + 5*AbsSqr(LamSU) + 18*Sqr(g1)));


   return oneLoop * beta_MDBS;
}

/**
 * Calculates the 2-loop beta function of MDBS.
 *
 * @return 2-loop beta function
 */
double MRSSM2_soft_parameters::calc_beta_MDBS_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_MDBS;

   beta_MDBS = Re(0.04*MDBS*(-150*traceYdAdjYd*AbsSqr(LamSD) - 50*traceYeAdjYe*
      AbsSqr(LamSD) - 150*traceYuAdjYu*AbsSqr(LamSU) - 150*AbsSqr(LamSD)*AbsSqr
      (LamTD) - 150*AbsSqr(LamSU)*AbsSqr(LamTU) + 208*Quad(g1) - 70*
      traceYdAdjYd*Sqr(g1) - 90*traceYeAdjYe*Sqr(g1) - 130*traceYuAdjYu*Sqr(g1)
      - 45*AbsSqr(LamTD)*Sqr(g1) - 45*AbsSqr(LamTU)*Sqr(g1) + 150*AbsSqr(LamSD)
      *Sqr(g2) + 150*AbsSqr(LamSU)*Sqr(g2) + 180*Sqr(g1)*Sqr(g2) + 440*Sqr(g1)*
      Sqr(g3) - 100*Sqr(LamSD)*Sqr(Conj(LamSD)) - 100*Sqr(LamSU)*Sqr(Conj(LamSU
      ))));


   return twoLoop * beta_MDBS;
}

/**
 * Calculates the 3-loop beta function of MDBS.
 *
 * @return 3-loop beta function
 */
double MRSSM2_soft_parameters::calc_beta_MDBS_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MDBS;

   beta_MDBS = 0;


   return threeLoop * beta_MDBS;
}

/**
 * Calculates the 4-loop beta function of MDBS.
 *
 * @return 4-loop beta function
 */
double MRSSM2_soft_parameters::calc_beta_MDBS_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MDBS;

   beta_MDBS = 0;


   return fourLoop * beta_MDBS;
}

/**
 * Calculates the 5-loop beta function of MDBS.
 *
 * @return 5-loop beta function
 */
double MRSSM2_soft_parameters::calc_beta_MDBS_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MDBS;

   beta_MDBS = 0;


   return fiveLoop * beta_MDBS;
}

} // namespace flexiblesusy
