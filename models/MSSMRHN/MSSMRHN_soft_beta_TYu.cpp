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
 * Calculates the 1-loop beta function of TYu.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_TYu_1_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = (0.06666666666666667*(90*traceAdjYuTYu*Yu + 30*traceAdjYvTYv*Yu +
      26*MassB*Yu*Sqr(g1) + 90*MassWB*Yu*Sqr(g2) + 160*MassG*Yu*Sqr(g3) + 45*
      traceYuAdjYu*TYu + 15*traceYvAdjYv*TYu - 13*Sqr(g1)*TYu - 45*Sqr(g2)*TYu
      - 80*Sqr(g3)*TYu) + 2*(Yu*Yd.adjoint()*TYd) + 4*(Yu*Yu.adjoint()*TYu) +
      TYu*Yd.adjoint()*Yd + 5*(TYu*Yu.adjoint()*Yu)).real();


   return oneLoop * beta_TYu;
}

/**
 * Calculates the 2-loop beta function of TYu.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_TYu_2_loop(const Soft_traces& soft_traces) const
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


   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = (0.0022222222222222222*(-2700*traceYdAdjYuTYuAdjYd*Yu - 900*
      traceYeAdjYvTYvAdjYe*Yu - 2700*traceYuAdjYdTYdAdjYu*Yu - 16200*
      traceYuAdjYuTYuAdjYu*Yu - 900*traceYvAdjYeTYeAdjYv*Yu - 5400*
      traceYvAdjYvTYvAdjYv*Yu - 10972*MassB*Yu*Quad(g1) - 13500*MassWB*Yu*Quad(
      g2) + 3200*MassG*Yu*Quad(g3) + 720*traceAdjYuTYu*Yu*Sqr(g1) - 720*MassB*
      traceYuAdjYu*Yu*Sqr(g1) - 900*MassB*Yu*Sqr(g1)*Sqr(g2) - 900*MassWB*Yu*
      Sqr(g1)*Sqr(g2) + 14400*traceAdjYuTYu*Yu*Sqr(g3) - 14400*MassG*
      traceYuAdjYu*Yu*Sqr(g3) - 2720*MassB*Yu*Sqr(g1)*Sqr(g3) - 2720*MassG*Yu*
      Sqr(g1)*Sqr(g3) - 7200*MassG*Yu*Sqr(g2)*Sqr(g3) - 7200*MassWB*Yu*Sqr(g2)*
      Sqr(g3) - 1350*traceYdAdjYuYuAdjYd*TYu - 450*traceYeAdjYvYvAdjYe*TYu -
      4050*traceYuAdjYuYuAdjYu*TYu - 1350*traceYvAdjYvYvAdjYv*TYu + 2743*Quad(
      g1)*TYu + 3375*Quad(g2)*TYu - 800*Quad(g3)*TYu + 360*traceYuAdjYu*Sqr(g1)
      *TYu + 450*Sqr(g1)*Sqr(g2)*TYu + 7200*traceYuAdjYu*Sqr(g3)*TYu + 1360*Sqr
      (g1)*Sqr(g3)*TYu + 3600*Sqr(g2)*Sqr(g3)*TYu) - 0.4*(15*traceAdjYdTYd + 5*
      traceAdjYeTYe + 2*MassB*Sqr(g1))*(Yu*Yd.adjoint()*Yd) + 0.4*(-15*
      traceYdAdjYd - 5*traceYeAdjYe + 2*Sqr(g1))*(Yu*Yd.adjoint()*TYd) - 0.4*(
      45*traceAdjYuTYu + 15*traceAdjYvTYv + 2*MassB*Sqr(g1) + 30*MassWB*Sqr(g2)
      )*(Yu*Yu.adjoint()*Yu) + 0.4*(-30*traceYuAdjYu - 10*traceYvAdjYv + 3*Sqr(
      g1) + 15*Sqr(g2))*(Yu*Yu.adjoint()*TYu) + 0.2*(-15*traceYdAdjYd - 5*
      traceYeAdjYe + 2*Sqr(g1))*(TYu*Yd.adjoint()*Yd) + (-15*traceYuAdjYu - 5*
      traceYvAdjYv + 12*Sqr(g2))*(TYu*Yu.adjoint()*Yu) - 4*(Yu*Yd.adjoint()*Yd*
      Yd.adjoint()*TYd) - 2*(Yu*Yd.adjoint()*Yd*Yu.adjoint()*TYu) - 4*(Yu*Yd.
      adjoint()*TYd*Yd.adjoint()*Yd) - 4*(Yu*Yd.adjoint()*TYd*Yu.adjoint()*Yu)
      - 6*(Yu*Yu.adjoint()*Yu*Yu.adjoint()*TYu) - 8*(Yu*Yu.adjoint()*TYu*Yu.
      adjoint()*Yu) - 2*(TYu*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 4*(TYu*Yd.
      adjoint()*Yd*Yu.adjoint()*Yu) - 6*(TYu*Yu.adjoint()*Yu*Yu.adjoint()*Yu)).
      real();


   return twoLoop * beta_TYu;
}

/**
 * Calculates the 3-loop beta function of TYu.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_TYu_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = ZEROMATRIX(3,3);


   return threeLoop * beta_TYu;
}

/**
 * Calculates the 4-loop beta function of TYu.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_TYu_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = ZEROMATRIX(3,3);


   return fourLoop * beta_TYu;
}

/**
 * Calculates the 5-loop beta function of TYu.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_TYu_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = ZEROMATRIX(3,3);


   return fiveLoop * beta_TYu;
}

} // namespace flexiblesusy
