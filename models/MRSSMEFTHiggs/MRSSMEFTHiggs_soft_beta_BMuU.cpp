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
 * Calculates the 1-loop beta function of BMuU.
 *
 * @return 1-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_BMuU_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_BMuU;

   beta_BMuU = Re(0.2*(15*traceYuAdjYu*BMuU + 30*AbsSqr(LamSU)*BMuU + 15*AbsSqr
      (LamTU)*BMuU + 20*LamSU*BMuD*Conj(LamSD) - 3*BMuU*Sqr(g1) - 15*BMuU*Sqr(
      g2)));


   return oneLoop * beta_BMuU;
}

/**
 * Calculates the 2-loop beta function of BMuU.
 *
 * @return 2-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_BMuU_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_BMuU;

   beta_BMuU = Re(0.1*(-30*traceYdAdjYuYuAdjYd*BMuU - 90*traceYuAdjYuYuAdjYu*
      BMuU - 150*traceYuAdjYu*AbsSqr(LamSU)*BMuU - 40*AbsSqr(LamSD)*AbsSqr(
      LamSU)*BMuU - 45*traceYuAdjYu*AbsSqr(LamTU)*BMuU - 180*AbsSqr(LamSU)*
      AbsSqr(LamTU)*BMuU - 30*AbsSqr(LamTD)*AbsSqr(LamTU)*BMuU - 120*LamSU*
      traceYdAdjYd*BMuD*Conj(LamSD) - 40*LamSU*traceYeAdjYe*BMuD*Conj(LamSD) -
      120*LamSU*AbsSqr(LamTD)*BMuD*Conj(LamSD) + 45*BMuU*Quad(g1) + 165*BMuU*
      Quad(g2) + 8*traceYuAdjYu*BMuU*Sqr(g1) + 24*AbsSqr(LamSU)*BMuU*Sqr(g1) +
      24*LamSU*BMuD*Conj(LamSD)*Sqr(g1) + 120*AbsSqr(LamSU)*BMuU*Sqr(g2) + 120*
      AbsSqr(LamTU)*BMuU*Sqr(g2) + 120*LamSU*BMuD*Conj(LamSD)*Sqr(g2) + 18*BMuU
      *Sqr(g1)*Sqr(g2) + 160*traceYuAdjYu*BMuU*Sqr(g3) - 80*LamSD*LamSU*BMuD*
      Sqr(Conj(LamSD)) - 140*BMuU*Sqr(LamSU)*Sqr(Conj(LamSU)) - 75*BMuU*Sqr(
      LamTU)*Sqr(Conj(LamTU))));


   return twoLoop * beta_BMuU;
}

/**
 * Calculates the 3-loop beta function of BMuU.
 *
 * @return 3-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_BMuU_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_BMuU;

   beta_BMuU = 0;


   return threeLoop * beta_BMuU;
}

/**
 * Calculates the 4-loop beta function of BMuU.
 *
 * @return 4-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_BMuU_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_BMuU;

   beta_BMuU = 0;


   return fourLoop * beta_BMuU;
}

/**
 * Calculates the 5-loop beta function of BMuU.
 *
 * @return 5-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_BMuU_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_BMuU;

   beta_BMuU = 0;


   return fiveLoop * beta_BMuU;
}

} // namespace flexiblesusy
