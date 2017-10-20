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

// File generated at Fri 20 Oct 2017 08:35:32

#include "HGTHDMIIMSSMBC_soft_parameters.hpp"
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
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_MassB_1_loop(const Soft_traces& soft_traces) const
{


   double beta_MassB;

   beta_MassB = Re(MassB*oneOver16PiSqr*(Sqr(g1dp) + Sqr(g2up)));


   return beta_MassB;
}

/**
 * Calculates the 2-loop beta function of MassB.
 *
 * @return 2-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_MassB_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_MassB;

   beta_MassB = Re(0.025*twoLoop*(51*MassB*Sqr(g1)*(Sqr(g1dp) + Sqr(g2up)
      ) - 5*(-(MassB*Quad(g1dp)) + 3*(7*MassB - 8*MassWB)*Sqr(g1d)*Sqr(g1dp) +
      3*MassB*Sqr(g1dp)*(4*(3*traceYdAdjYd + traceYeAdjYe) - 17*Sqr(g2)) - Sqr(
      g2up)*(51*MassB*Sqr(g2) - 3*(7*MassB - 8*MassWB)*Sqr(g2u) + MassB*(-36*
      traceYuAdjYu + Sqr(g2up))))));


   return beta_MassB;
}

/**
 * Calculates the 3-loop beta function of MassB.
 *
 * @return 3-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_MassB_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MassB;

   beta_MassB = 0;


   return beta_MassB;
}

} // namespace flexiblesusy
