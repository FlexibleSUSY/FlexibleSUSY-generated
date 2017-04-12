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

// File generated at Wed 12 Apr 2017 11:14:25

#include "SplitMSSM_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of MassB.
 *
 * @return one-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_MassB_one_loop(const Soft_traces& soft_traces) const
{


   double beta_MassB;

   beta_MassB = Re(oneOver16PiSqr*(4*gYd*gYu*Mu + MassB*Sqr(gYd) + MassB*
      Sqr(gYu)));


   return beta_MassB;
}

/**
 * Calculates the two-loop beta function of MassB.
 *
 * @return two-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_MassB_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_MassB;

   beta_MassB = Re(0.025*twoLoop*(3*Sqr(g1)*(64*gYd*gYu*Mu + 17*MassB*Sqr
      (gYd) + 17*MassB*Sqr(gYu)) - 5*(-(Power(gYd,4)*MassB) - Power(gYu,4)*
      MassB + 48*g2d*g2u*gYd*gYu*MassWB + 8*Power(gYd,3)*gYu*Mu + 8*gYd*Power(
      gYu,3)*Mu + 24*gYd*gYu*Mu*Sqr(g2u) + 36*MassB*traceYdAdjYd*Sqr(gYd) + 12*
      MassB*traceYeAdjYe*Sqr(gYd) + 36*MassB*traceYuAdjYu*Sqr(gYd) + 18*MassB*
      Sqr(g2u)*Sqr(gYd) + 36*MassB*traceYdAdjYd*Sqr(gYu) + 12*MassB*
      traceYeAdjYe*Sqr(gYu) + 36*MassB*traceYuAdjYu*Sqr(gYu) + 21*MassB*Sqr(g2u
      )*Sqr(gYu) - 24*MassWB*Sqr(g2u)*Sqr(gYu) + 28*MassB*Sqr(gYd)*Sqr(gYu) + 3
      *Sqr(g2d)*(8*gYd*gYu*Mu + (7*MassB - 8*MassWB)*Sqr(gYd) + 6*MassB*Sqr(gYu
      )) - 3*Sqr(g2)*(64*gYd*gYu*Mu + 17*MassB*Sqr(gYd) + 17*MassB*Sqr(gYu)))))
      ;


   return beta_MassB;
}

/**
 * Calculates the three-loop beta function of MassB.
 *
 * @return three-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_MassB_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MassB;

   beta_MassB = 0;


   return beta_MassB;
}

} // namespace flexiblesusy
