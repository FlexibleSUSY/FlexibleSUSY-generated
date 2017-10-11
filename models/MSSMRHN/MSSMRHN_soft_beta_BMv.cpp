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

// File generated at Tue 10 Oct 2017 22:48:57

#include "MSSMRHN_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of BMv.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_BMv_1_loop(const Soft_traces& soft_traces) const
{


   Eigen::Matrix<double,3,3> beta_BMv;

   beta_BMv = (oneOver16PiSqr*(4*(Mv*Yv.conjugate()*(TYv).transpose()) +
      2*(Yv*Yv.adjoint()*BMv) + 2*(BMv*Yv.conjugate()*Yv.transpose()) + 4*(TYv*
      Yv.adjoint()*Mv))).real();


   return beta_BMv;
}

/**
 * Calculates the 2-loop beta function of BMv.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_BMv_2_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   Eigen::Matrix<double,3,3> beta_BMv;

   beta_BMv = (twoLoop*(-0.8*(3*MassB*Sqr(g1) + 5*(3*traceAdjYuTYu +
      traceAdjYvTYv + 3*MassWB*Sqr(g2)))*(Mv*Yv.conjugate()*Yv.transpose()) + (
      -12*traceYuAdjYu - 4*traceYvAdjYv + 2.4*Sqr(g1) + 12*Sqr(g2))*(Mv*
      Yv.conjugate()*(TYv).transpose()) - 0.8*(3*MassB*Sqr(g1) + 5*(3*
      traceAdjYuTYu + traceAdjYvTYv + 3*MassWB*Sqr(g2)))*(Yv*Yv.adjoint()*Mv) +
      (-2*(3*traceYuAdjYu + traceYvAdjYv) + 1.2*Sqr(g1) + 6*Sqr(g2))*(Yv*
      Yv.adjoint()*BMv) + (-2*(3*traceYuAdjYu + traceYvAdjYv) + 1.2*Sqr(g1) + 6
      *Sqr(g2))*(BMv*Yv.conjugate()*Yv.transpose()) + (-12*traceYuAdjYu - 4*
      traceYvAdjYv + 2.4*Sqr(g1) + 12*Sqr(g2))*(TYv*Yv.adjoint()*Mv) - 4*(Mv*
      Yv.conjugate()*Ye.transpose()*Ye.conjugate()*(TYv).transpose()) - 4*(Mv*
      Yv.conjugate()*Yv.transpose()*Yv.conjugate()*(TYv).transpose()) - 4*(Mv*
      Yv.conjugate()*(TYe).transpose()*Ye.conjugate()*Yv.transpose()) - 4*(Mv*
      Yv.conjugate()*(TYv).transpose()*Yv.conjugate()*Yv.transpose()) - 2*(Yv*
      Ye.adjoint()*Ye*Yv.adjoint()*BMv) - 4*(Yv*Ye.adjoint()*TYe*Yv.adjoint()*
      Mv) - 2*(Yv*Yv.adjoint()*Yv*Yv.adjoint()*BMv) - 4*(Yv*Yv.adjoint()*TYv*
      Yv.adjoint()*Mv) - 2*(BMv*Yv.conjugate()*Ye.transpose()*Ye.conjugate()*
      Yv.transpose()) - 2*(BMv*Yv.conjugate()*Yv.transpose()*Yv.conjugate()*
      Yv.transpose()) - 4*(TYv*Ye.adjoint()*Ye*Yv.adjoint()*Mv) - 4*(TYv*
      Yv.adjoint()*Yv*Yv.adjoint()*Mv))).real();


   return beta_BMv;
}

/**
 * Calculates the 3-loop beta function of BMv.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_BMv_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_BMv;

   beta_BMv = ZEROMATRIX(3,3);


   return beta_BMv;
}

} // namespace flexiblesusy
