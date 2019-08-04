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

// File generated at Sun 4 Aug 2019 17:12:06

#include "CE6SSM_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of mHp2.
 *
 * @return 1-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_mHp2_1_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   double beta_mHp2;

   beta_mHp2 = Re(0.2*oneOver16PiSqr*(-3.872983346207417*g1*Tr11 +
      3.1622776601683795*gN*Tr14 - 6*AbsSqr(MassB)*Sqr(g1) - 30*AbsSqr(MassWB)*
      Sqr(g2) - 4*AbsSqr(MassBp)*Sqr(gN)));


   return beta_mHp2;
}

/**
 * Calculates the 2-loop beta function of mHp2.
 *
 * @return 2-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_mHp2_2_loop(const Soft_traces& soft_traces) const
{
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_mHp2;

   beta_mHp2 = Re(0.04*twoLoop*(-24.49489742783178*g1*gN*Tr2U114 -
      24.49489742783178*g1*gN*Tr2U141 - 77.45966692414834*g1*Tr31 +
      63.24555320336759*gN*Tr34 + 891*AbsSqr(MassB)*Quad(g1) + 150*Tr22*Quad(g2
      ) + 2175*AbsSqr(MassWB)*Quad(g2) + 576*AbsSqr(MassBp)*Quad(gN) + 30*
      Tr2U111*Sqr(g1) + 90*AbsSqr(MassB)*Sqr(g1)*Sqr(g2) + 90*AbsSqr(MassWB)*
      Sqr(g1)*Sqr(g2) + 45*MassWB*Conj(MassB)*Sqr(g1)*Sqr(g2) + 45*MassB*Conj(
      MassWB)*Sqr(g1)*Sqr(g2) + 20*Tr2U144*Sqr(gN) + 36*AbsSqr(MassB)*Sqr(g1)*
      Sqr(gN) + 36*AbsSqr(MassBp)*Sqr(g1)*Sqr(gN) + 18*MassBp*Conj(MassB)*Sqr(
      g1)*Sqr(gN) + 18*MassB*Conj(MassBp)*Sqr(g1)*Sqr(gN) + 60*AbsSqr(MassBp)*
      Sqr(g2)*Sqr(gN) + 60*AbsSqr(MassWB)*Sqr(g2)*Sqr(gN) + 30*MassWB*Conj(
      MassBp)*Sqr(g2)*Sqr(gN) + 30*MassBp*Conj(MassWB)*Sqr(g2)*Sqr(gN)));


   return beta_mHp2;
}

/**
 * Calculates the 3-loop beta function of mHp2.
 *
 * @return 3-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_mHp2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHp2;

   beta_mHp2 = 0;


   return beta_mHp2;
}

/**
 * Calculates the 4-loop beta function of mHp2.
 *
 * @return 4-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_mHp2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHp2;

   beta_mHp2 = 0;


   return beta_mHp2;
}

/**
 * Calculates the 5-loop beta function of mHp2.
 *
 * @return 5-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_mHp2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHp2;

   beta_mHp2 = 0;


   return beta_mHp2;
}

} // namespace flexiblesusy
