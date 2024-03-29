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
 * Calculates the 1-loop beta function of MassB.
 *
 * @return 1-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_MassB_1_loop(const Soft_traces& soft_traces) const
{


   double beta_MassB;

   beta_MassB = Re(4*gYd*gYu*Mu + MassB*Sqr(gYd) + MassB*Sqr(gYu));


   return oneLoop * beta_MassB;
}

/**
 * Calculates the 2-loop beta function of MassB.
 *
 * @return 2-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_MassB_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_MassB;

   beta_MassB = Re(0.025*(-240*g2d*g2u*gYd*gYu*MassWB - 480*gYd*gYu*
      traceYdAdjYd*Mu - 160*gYd*gYu*traceYeAdjYe*Mu - 480*gYd*gYu*traceYuAdjYu*
      Mu - 120*gYu*Cube(gYd)*Mu - 120*gYd*Cube(gYu)*Mu + 5*MassB*Quad(gYd) + 5*
      MassB*Quad(gYu) + 192*gYd*gYu*Mu*Sqr(g1) + 960*gYd*gYu*Mu*Sqr(g2) - 360*
      gYd*gYu*Mu*Sqr(g2d) - 360*gYd*gYu*Mu*Sqr(g2u) - 180*MassB*traceYdAdjYd*
      Sqr(gYd) - 60*MassB*traceYeAdjYe*Sqr(gYd) - 180*MassB*traceYuAdjYu*Sqr(
      gYd) + 51*MassB*Sqr(g1)*Sqr(gYd) + 255*MassB*Sqr(g2)*Sqr(gYd) - 105*MassB
      *Sqr(g2d)*Sqr(gYd) + 120*MassWB*Sqr(g2d)*Sqr(gYd) - 90*MassB*Sqr(g2u)*Sqr
      (gYd) - 180*MassB*traceYdAdjYd*Sqr(gYu) - 60*MassB*traceYeAdjYe*Sqr(gYu)
      - 180*MassB*traceYuAdjYu*Sqr(gYu) + 51*MassB*Sqr(g1)*Sqr(gYu) + 255*MassB
      *Sqr(g2)*Sqr(gYu) - 90*MassB*Sqr(g2d)*Sqr(gYu) - 105*MassB*Sqr(g2u)*Sqr(
      gYu) + 120*MassWB*Sqr(g2u)*Sqr(gYu) - 140*MassB*Sqr(gYd)*Sqr(gYu)));


   return twoLoop * beta_MassB;
}

/**
 * Calculates the 3-loop beta function of MassB.
 *
 * @return 3-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_MassB_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MassB;

   beta_MassB = 0;


   return threeLoop * beta_MassB;
}

/**
 * Calculates the 4-loop beta function of MassB.
 *
 * @return 4-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_MassB_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MassB;

   beta_MassB = 0;


   return fourLoop * beta_MassB;
}

/**
 * Calculates the 5-loop beta function of MassB.
 *
 * @return 5-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_MassB_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MassB;

   beta_MassB = 0;


   return fiveLoop * beta_MassB;
}

} // namespace flexiblesusy
