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
 * Calculates the 1-loop beta function of Mu.
 *
 * @return 1-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_Mu_1_loop(const Soft_traces& soft_traces) const
{


   double beta_Mu;

   beta_Mu = Re(0.05*(20*gYd*gYu*MassB + 60*g2d*g2u*MassWB - 18*Mu*Sqr(g1) - 90
      *Mu*Sqr(g2) + 15*Mu*Sqr(g2d) + 15*Mu*Sqr(g2u) + 5*Mu*Sqr(gYd) + 5*Mu*Sqr(
      gYu)));


   return oneLoop * beta_Mu;
}

/**
 * Calculates the 2-loop beta function of Mu.
 *
 * @return 2-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_Mu_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_Mu;

   beta_Mu = Re(0.00125*(-2400*gYd*gYu*MassB*traceYdAdjYd - 7200*g2d*g2u*MassWB
      *traceYdAdjYd - 800*gYd*gYu*MassB*traceYeAdjYe - 2400*g2d*g2u*MassWB*
      traceYeAdjYe - 2400*gYd*gYu*MassB*traceYuAdjYu - 7200*g2d*g2u*MassWB*
      traceYuAdjYu - 6000*g2u*MassWB*Cube(g2d) - 6000*g2d*MassWB*Cube(g2u) -
      1200*gYu*MassB*Cube(gYd) - 1200*gYd*MassB*Cube(gYu) + 2400*g2d*g2u*gYd*
      gYu*Mu + 2718*Mu*Quad(g1) - 21050*Mu*Quad(g2) - 1500*Mu*Quad(g2d) - 1500*
      Mu*Quad(g2u) - 200*Mu*Quad(gYd) - 200*Mu*Quad(gYu) + 720*gYd*gYu*MassB*
      Sqr(g1) + 2160*g2d*g2u*MassWB*Sqr(g1) + 3600*gYd*gYu*MassB*Sqr(g2) +
      34800*g2d*g2u*MassWB*Sqr(g2) - 540*Mu*Sqr(g1)*Sqr(g2) - 1200*gYd*gYu*
      MassB*Sqr(g2d) - 2700*traceYdAdjYd*Mu*Sqr(g2d) - 900*traceYeAdjYe*Mu*Sqr(
      g2d) - 2700*traceYuAdjYu*Mu*Sqr(g2d) + 495*Mu*Sqr(g1)*Sqr(g2d) + 9075*Mu*
      Sqr(g2)*Sqr(g2d) - 1200*gYd*gYu*MassB*Sqr(g2u) - 2700*traceYdAdjYd*Mu*Sqr
      (g2u) - 900*traceYeAdjYe*Mu*Sqr(g2u) - 2700*traceYuAdjYu*Mu*Sqr(g2u) +
      495*Mu*Sqr(g1)*Sqr(g2u) + 9075*Mu*Sqr(g2)*Sqr(g2u) - 9000*Mu*Sqr(g2d)*Sqr
      (g2u) - 1200*g2d*g2u*MassWB*Sqr(gYd) - 900*traceYdAdjYd*Mu*Sqr(gYd) - 300
      *traceYeAdjYe*Mu*Sqr(gYd) - 900*traceYuAdjYu*Mu*Sqr(gYd) + 165*Mu*Sqr(g1)
      *Sqr(gYd) + 825*Mu*Sqr(g2)*Sqr(gYd) - 900*Mu*Sqr(g2d)*Sqr(gYd) - 900*Mu*
      Sqr(g2u)*Sqr(gYd) - 1200*g2d*g2u*MassWB*Sqr(gYu) - 900*traceYdAdjYd*Mu*
      Sqr(gYu) - 300*traceYeAdjYe*Mu*Sqr(gYu) - 900*traceYuAdjYu*Mu*Sqr(gYu) +
      165*Mu*Sqr(g1)*Sqr(gYu) + 825*Mu*Sqr(g2)*Sqr(gYu) - 900*Mu*Sqr(g2d)*Sqr(
      gYu) - 900*Mu*Sqr(g2u)*Sqr(gYu) - 1600*Mu*Sqr(gYd)*Sqr(gYu)));


   return twoLoop * beta_Mu;
}

/**
 * Calculates the 3-loop beta function of Mu.
 *
 * @return 3-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_Mu_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Mu;

   beta_Mu = 0;


   return threeLoop * beta_Mu;
}

/**
 * Calculates the 4-loop beta function of Mu.
 *
 * @return 4-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_Mu_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Mu;

   beta_Mu = 0;


   return fourLoop * beta_Mu;
}

/**
 * Calculates the 5-loop beta function of Mu.
 *
 * @return 5-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_Mu_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Mu;

   beta_Mu = 0;


   return fiveLoop * beta_Mu;
}

} // namespace flexiblesusy
