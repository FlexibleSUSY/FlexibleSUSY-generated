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

// File generated at Sun 4 Aug 2019 17:35:59

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
 * Calculates the 1-loop beta function of BMuD.
 *
 * @return 1-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_BMuD_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_BMuD;

   beta_BMuD = Re(0.2*oneOver16PiSqr*(15*traceYdAdjYd*BMuD + 5*traceYeAdjYe*
      BMuD + 30*AbsSqr(LamSD)*BMuD + 15*AbsSqr(LamTD)*BMuD + 20*LamSD*BMuU*Conj
      (LamSU) - 3*BMuD*Sqr(g1) - 15*BMuD*Sqr(g2)));


   return beta_BMuD;
}

/**
 * Calculates the 2-loop beta function of BMuD.
 *
 * @return 2-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_BMuD_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_BMuD;

   beta_BMuD = Re(0.1*twoLoop*(-90*traceYdAdjYdYdAdjYd*BMuD - 30*
      traceYdAdjYuYuAdjYd*BMuD - 30*traceYeAdjYeYeAdjYe*BMuD - 150*traceYdAdjYd
      *AbsSqr(LamSD)*BMuD - 50*traceYeAdjYe*AbsSqr(LamSD)*BMuD - 40*AbsSqr(
      LamSD)*AbsSqr(LamSU)*BMuD - 45*traceYdAdjYd*AbsSqr(LamTD)*BMuD - 15*
      traceYeAdjYe*AbsSqr(LamTD)*BMuD - 180*AbsSqr(LamSD)*AbsSqr(LamTD)*BMuD -
      30*AbsSqr(LamTD)*AbsSqr(LamTU)*BMuD - 120*LamSD*traceYuAdjYu*BMuU*Conj(
      LamSU) - 120*LamSD*AbsSqr(LamTU)*BMuU*Conj(LamSU) + 45*BMuD*Quad(g1) +
      165*BMuD*Quad(g2) - 4*traceYdAdjYd*BMuD*Sqr(g1) + 12*traceYeAdjYe*BMuD*
      Sqr(g1) + 72*AbsSqr(LamSD)*BMuD*Sqr(g1) + 72*LamSD*BMuU*Conj(LamSU)*Sqr(
      g1) + 360*AbsSqr(LamSD)*BMuD*Sqr(g2) + 120*AbsSqr(LamTD)*BMuD*Sqr(g2) +
      360*LamSD*BMuU*Conj(LamSU)*Sqr(g2) + 18*BMuD*Sqr(g1)*Sqr(g2) + 160*
      traceYdAdjYd*BMuD*Sqr(g3) - 140*BMuD*Sqr(LamSD)*Sqr(Conj(LamSD)) - 80*
      LamSD*LamSU*BMuU*Sqr(Conj(LamSU)) - 75*BMuD*Sqr(LamTD)*Sqr(Conj(LamTD))))
      ;


   return beta_BMuD;
}

/**
 * Calculates the 3-loop beta function of BMuD.
 *
 * @return 3-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_BMuD_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_BMuD;

   beta_BMuD = 0;


   return beta_BMuD;
}

/**
 * Calculates the 4-loop beta function of BMuD.
 *
 * @return 4-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_BMuD_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_BMuD;

   beta_BMuD = 0;


   return beta_BMuD;
}

/**
 * Calculates the 5-loop beta function of BMuD.
 *
 * @return 5-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_BMuD_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_BMuD;

   beta_BMuD = 0;


   return beta_BMuD;
}

} // namespace flexiblesusy
