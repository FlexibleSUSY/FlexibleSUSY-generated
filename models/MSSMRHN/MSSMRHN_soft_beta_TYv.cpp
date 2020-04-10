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

// File generated at Fri 10 Apr 2020 20:35:12

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
 * Calculates the 1-loop beta function of TYv.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_TYv_1_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   Eigen::Matrix<double,3,3> beta_TYv;

   beta_TYv = (oneOver16PiSqr*(0.2*(30*traceAdjYuTYu*Yv + 10*traceAdjYvTYv*Yv +
      6*MassB*Yv*Sqr(g1) + 30*MassWB*Yv*Sqr(g2) + 15*traceYuAdjYu*TYv + 5*
      traceYvAdjYv*TYv - 3*Sqr(g1)*TYv - 15*Sqr(g2)*TYv) + 2*(Yv*Ye.adjoint()*
      TYe) + 4*(Yv*Yv.adjoint()*TYv) + TYv*Ye.adjoint()*Ye + 5*(TYv*Yv.adjoint(
      )*Yv))).real();


   return beta_TYv;
}

/**
 * Calculates the 2-loop beta function of TYv.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_TYv_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYvTYvAdjYe = TRACE_STRUCT.traceYeAdjYvTYvAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYuAdjYuTYuAdjYu = TRACE_STRUCT.traceYuAdjYuTYuAdjYu;
   const double traceYvAdjYeTYeAdjYv = TRACE_STRUCT.traceYvAdjYeTYeAdjYv;
   const double traceYvAdjYvTYvAdjYv = TRACE_STRUCT.traceYvAdjYvTYvAdjYv;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYvYvAdjYe = TRACE_STRUCT.traceYeAdjYvYvAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYvAdjYvYvAdjYv = TRACE_STRUCT.traceYvAdjYvYvAdjYv;


   Eigen::Matrix<double,3,3> beta_TYv;

   beta_TYv = (twoLoop*(0.02*(-300*traceYdAdjYuTYuAdjYd*Yv - 100*
      traceYeAdjYvTYvAdjYe*Yv - 300*traceYuAdjYdTYdAdjYu*Yv - 1800*
      traceYuAdjYuTYuAdjYu*Yv - 100*traceYvAdjYeTYeAdjYv*Yv - 600*
      traceYvAdjYvTYvAdjYv*Yv - 828*MassB*Yv*Quad(g1) - 1500*MassWB*Yv*Quad(g2)
      + 80*traceAdjYuTYu*Yv*Sqr(g1) - 80*MassB*traceYuAdjYu*Yv*Sqr(g1) - 180*
      MassB*Yv*Sqr(g1)*Sqr(g2) - 180*MassWB*Yv*Sqr(g1)*Sqr(g2) + 1600*
      traceAdjYuTYu*Yv*Sqr(g3) - 1600*MassG*traceYuAdjYu*Yv*Sqr(g3) - 150*
      traceYdAdjYuYuAdjYd*TYv - 50*traceYeAdjYvYvAdjYe*TYv - 450*
      traceYuAdjYuYuAdjYu*TYv - 150*traceYvAdjYvYvAdjYv*TYv + 207*Quad(g1)*TYv
      + 375*Quad(g2)*TYv + 40*traceYuAdjYu*Sqr(g1)*TYv + 90*Sqr(g1)*Sqr(g2)*TYv
       + 800*traceYuAdjYu*Sqr(g3)*TYv) - 0.4*(15*traceAdjYdTYd + 5*
      traceAdjYeTYe + 6*MassB*Sqr(g1))*(Yv*Ye.adjoint()*Ye) + 0.4*(-15*
      traceYdAdjYd - 5*traceYeAdjYe + 6*Sqr(g1))*(Yv*Ye.adjoint()*TYe) - 1.2*(
      15*traceAdjYuTYu + 5*traceAdjYvTYv + 2*MassB*Sqr(g1) + 10*MassWB*Sqr(g2))
      *(Yv*Yv.adjoint()*Yv) + 0.4*(-30*traceYuAdjYu - 10*traceYvAdjYv + 3*Sqr(
      g1) + 15*Sqr(g2))*(Yv*Yv.adjoint()*TYv) + 0.2*(-15*traceYdAdjYd - 5*
      traceYeAdjYe + 6*Sqr(g1))*(TYv*Ye.adjoint()*Ye) + 0.2*(-75*traceYuAdjYu -
      25*traceYvAdjYv + 12*Sqr(g1) + 60*Sqr(g2))*(TYv*Yv.adjoint()*Yv) - 4*(Yv*
      Ye.adjoint()*Ye*Ye.adjoint()*TYe) - 2*(Yv*Ye.adjoint()*Ye*Yv.adjoint()*
      TYv) - 4*(Yv*Ye.adjoint()*TYe*Ye.adjoint()*Ye) - 4*(Yv*Ye.adjoint()*TYe*
      Yv.adjoint()*Yv) - 6*(Yv*Yv.adjoint()*Yv*Yv.adjoint()*TYv) - 8*(Yv*Yv.
      adjoint()*TYv*Yv.adjoint()*Yv) - 2*(TYv*Ye.adjoint()*Ye*Ye.adjoint()*Ye)
      - 4*(TYv*Ye.adjoint()*Ye*Yv.adjoint()*Yv) - 6*(TYv*Yv.adjoint()*Yv*Yv.
      adjoint()*Yv))).real();


   return beta_TYv;
}

/**
 * Calculates the 3-loop beta function of TYv.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_TYv_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYv;

   beta_TYv = ZEROMATRIX(3,3);


   return beta_TYv;
}

/**
 * Calculates the 4-loop beta function of TYv.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_TYv_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYv;

   beta_TYv = ZEROMATRIX(3,3);


   return beta_TYv;
}

/**
 * Calculates the 5-loop beta function of TYv.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_TYv_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYv;

   beta_TYv = ZEROMATRIX(3,3);


   return beta_TYv;
}

} // namespace flexiblesusy
