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


#include "SplitMSSM_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of MassWB.
 *
 * @return 1-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_MassWB_1_loop(const Soft_traces& soft_traces) const
{


   double beta_MassWB;

   beta_MassWB = Re(4*g2d*g2u*Mu - 12*MassWB*Sqr(g2) + MassWB*Sqr(g2d) + MassWB
      *Sqr(g2u));


   return oneLoop * beta_MassWB;
}

/**
 * Calculates the 2-loop beta function of MassWB.
 *
 * @return 2-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_MassWB_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_MassWB;

   beta_MassWB = Re(0.008333333333333333*(-240*g2d*g2u*gYd*gYu*MassB - 1440*g2d
      *g2u*traceYdAdjYd*Mu - 480*g2d*g2u*traceYeAdjYe*Mu - 1440*g2d*g2u*
      traceYuAdjYu*Mu - 1080*g2u*Cube(g2d)*Mu - 1080*g2d*Cube(g2u)*Mu - 9320*
      MassWB*Quad(g2) - 435*MassWB*Quad(g2d) - 435*MassWB*Quad(g2u) + 576*g2d*
      g2u*Mu*Sqr(g1) + 5760*g2d*g2u*Mu*Sqr(g2) - 540*MassWB*traceYdAdjYd*Sqr(
      g2d) - 180*MassWB*traceYeAdjYe*Sqr(g2d) - 540*MassWB*traceYuAdjYu*Sqr(g2d
      ) + 153*MassWB*Sqr(g1)*Sqr(g2d) + 1365*MassWB*Sqr(g2)*Sqr(g2d) - 540*
      MassWB*traceYdAdjYd*Sqr(g2u) - 180*MassWB*traceYeAdjYe*Sqr(g2u) - 540*
      MassWB*traceYuAdjYu*Sqr(g2u) + 153*MassWB*Sqr(g1)*Sqr(g2u) + 1365*MassWB*
      Sqr(g2)*Sqr(g2u) - 1260*MassWB*Sqr(g2d)*Sqr(g2u) - 360*g2d*g2u*Mu*Sqr(gYd
      ) + 120*MassB*Sqr(g2d)*Sqr(gYd) - 105*MassWB*Sqr(g2d)*Sqr(gYd) - 90*
      MassWB*Sqr(g2u)*Sqr(gYd) - 360*g2d*g2u*Mu*Sqr(gYu) - 90*MassWB*Sqr(g2d)*
      Sqr(gYu) + 120*MassB*Sqr(g2u)*Sqr(gYu) - 105*MassWB*Sqr(g2u)*Sqr(gYu)));


   return twoLoop * beta_MassWB;
}

/**
 * Calculates the 3-loop beta function of MassWB.
 *
 * @return 3-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_MassWB_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MassWB;

   beta_MassWB = 0;


   return threeLoop * beta_MassWB;
}

/**
 * Calculates the 4-loop beta function of MassWB.
 *
 * @return 4-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_MassWB_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MassWB;

   beta_MassWB = 0;


   return fourLoop * beta_MassWB;
}

/**
 * Calculates the 5-loop beta function of MassWB.
 *
 * @return 5-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_MassWB_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MassWB;

   beta_MassWB = 0;


   return fiveLoop * beta_MassWB;
}

} // namespace flexiblesusy
