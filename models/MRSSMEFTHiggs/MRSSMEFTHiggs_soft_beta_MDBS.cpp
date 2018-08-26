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

// File generated at Sun 26 Aug 2018 14:25:29

#include "MRSSMEFTHiggs_soft_parameters.hpp"
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
double MRSSMEFTHiggs_soft_parameters::calc_beta_MDBS_1_loop(const Soft_traces& soft_traces) const
{


   double beta_MDBS;

   beta_MDBS = Re(0.4*MDBS*oneOver16PiSqr*(5*AbsSqr(LamSD) + 5*AbsSqr(LamSU) +
      18*Sqr(g1)));


   return beta_MDBS;
}

/**
 * Calculates the 2-loop beta function of MDBS.
 *
 * @return 2-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_MDBS_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_MDBS;

   beta_MDBS = Re(-0.04*MDBS*twoLoop*(-150*AbsSqr(LamSU)*(-traceYuAdjYu -
      AbsSqr(LamTU) + Sqr(g2)) - 50*AbsSqr(LamSD)*(-3*traceYdAdjYd -
      traceYeAdjYe - 3*AbsSqr(LamTD) + 3*Sqr(g2)) + Sqr(g1)*(70*traceYdAdjYd +
      90*traceYeAdjYe + 130*traceYuAdjYu + 45*AbsSqr(LamTD) + 45*AbsSqr(LamTU)
      - 208*Sqr(g1) - 180*Sqr(g2) - 440*Sqr(g3)) + 100*Sqr(LamSD)*Sqr(Conj(
      LamSD)) + 100*Sqr(LamSU)*Sqr(Conj(LamSU))));


   return beta_MDBS;
}

/**
 * Calculates the 3-loop beta function of MDBS.
 *
 * @return 3-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_MDBS_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MDBS;

   beta_MDBS = 0;


   return beta_MDBS;
}

/**
 * Calculates the 4-loop beta function of MDBS.
 *
 * @return 4-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_MDBS_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MDBS;

   beta_MDBS = 0;


   return beta_MDBS;
}

} // namespace flexiblesusy
