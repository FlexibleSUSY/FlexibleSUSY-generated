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
 * Calculates the 1-loop beta function of TYd.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_TYd_1_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = (0.06666666666666667*(90*traceAdjYdTYd*Yd + 30*traceAdjYeTYe*Yd +
      14*MassB*Yd*Sqr(g1) + 90*MassWB*Yd*Sqr(g2) + 160*MassG*Yd*Sqr(g3) + 45*
      traceYdAdjYd*TYd + 15*traceYeAdjYe*TYd - 7*Sqr(g1)*TYd - 45*Sqr(g2)*TYd -
      80*Sqr(g3)*TYd) + 4*(Yd*Yd.adjoint()*TYd) + 2*(Yd*Yu.adjoint()*TYu) + 5*(
      TYd*Yd.adjoint()*Yd) + TYd*Yu.adjoint()*Yu).real();


   return oneLoop * beta_TYd;
}

/**
 * Calculates the 2-loop beta function of TYd.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_TYd_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYeAdjYvTYvAdjYe = TRACE_STRUCT.traceYeAdjYvTYvAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYvAdjYeTYeAdjYv = TRACE_STRUCT.traceYvAdjYeTYeAdjYv;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYeAdjYvYvAdjYe = TRACE_STRUCT.traceYeAdjYvYvAdjYe;


   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = (0.011111111111111112*(-3240*traceYdAdjYdTYdAdjYd*Yd - 540*
      traceYdAdjYuTYuAdjYd*Yd - 1080*traceYeAdjYeTYeAdjYe*Yd - 180*
      traceYeAdjYvTYvAdjYe*Yd - 540*traceYuAdjYdTYdAdjYu*Yd - 180*
      traceYvAdjYeTYeAdjYv*Yd - 1148*MassB*Yd*Quad(g1) - 2700*MassWB*Yd*Quad(g2
      ) + 640*MassG*Yd*Quad(g3) - 72*traceAdjYdTYd*Yd*Sqr(g1) + 216*
      traceAdjYeTYe*Yd*Sqr(g1) + 72*MassB*traceYdAdjYd*Yd*Sqr(g1) - 216*MassB*
      traceYeAdjYe*Yd*Sqr(g1) - 180*MassB*Yd*Sqr(g1)*Sqr(g2) - 180*MassWB*Yd*
      Sqr(g1)*Sqr(g2) + 2880*traceAdjYdTYd*Yd*Sqr(g3) - 2880*MassG*traceYdAdjYd
      *Yd*Sqr(g3) - 160*MassB*Yd*Sqr(g1)*Sqr(g3) - 160*MassG*Yd*Sqr(g1)*Sqr(g3)
      - 1440*MassG*Yd*Sqr(g2)*Sqr(g3) - 1440*MassWB*Yd*Sqr(g2)*Sqr(g3) - 810*
      traceYdAdjYdYdAdjYd*TYd - 270*traceYdAdjYuYuAdjYd*TYd - 270*
      traceYeAdjYeYeAdjYe*TYd - 90*traceYeAdjYvYvAdjYe*TYd + 287*Quad(g1)*TYd +
      675*Quad(g2)*TYd - 160*Quad(g3)*TYd - 36*traceYdAdjYd*Sqr(g1)*TYd + 108*
      traceYeAdjYe*Sqr(g1)*TYd + 90*Sqr(g1)*Sqr(g2)*TYd + 1440*traceYdAdjYd*Sqr
      (g3)*TYd + 80*Sqr(g1)*Sqr(g3)*TYd + 720*Sqr(g2)*Sqr(g3)*TYd) - 0.4*(45*
      traceAdjYdTYd + 15*traceAdjYeTYe + 4*MassB*Sqr(g1) + 30*MassWB*Sqr(g2))*(
      Yd*Yd.adjoint()*Yd) + 0.4*(-30*traceYdAdjYd - 10*traceYeAdjYe + 3*Sqr(g1)
      + 15*Sqr(g2))*(Yd*Yd.adjoint()*TYd) - 0.4*(15*traceAdjYuTYu + 5*
      traceAdjYvTYv + 4*MassB*Sqr(g1))*(Yd*Yu.adjoint()*Yu) + 0.4*(-15*
      traceYuAdjYu - 5*traceYvAdjYv + 4*Sqr(g1))*(Yd*Yu.adjoint()*TYu) + 0.2*(-
      75*traceYdAdjYd - 25*traceYeAdjYe + 6*Sqr(g1) + 60*Sqr(g2))*(TYd*Yd.
      adjoint()*Yd) + 0.2*(-15*traceYuAdjYu - 5*traceYvAdjYv + 4*Sqr(g1))*(TYd*
      Yu.adjoint()*Yu) - 6*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*TYd) - 8*(Yd*Yd.
      adjoint()*TYd*Yd.adjoint()*Yd) - 2*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*TYd)
      - 4*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*TYu) - 4*(Yd*Yu.adjoint()*TYu*Yd.
      adjoint()*Yd) - 4*(Yd*Yu.adjoint()*TYu*Yu.adjoint()*Yu) - 6*(TYd*Yd.
      adjoint()*Yd*Yd.adjoint()*Yd) - 4*(TYd*Yu.adjoint()*Yu*Yd.adjoint()*Yd) -
      2*(TYd*Yu.adjoint()*Yu*Yu.adjoint()*Yu)).real();


   return twoLoop * beta_TYd;
}

/**
 * Calculates the 3-loop beta function of TYd.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_TYd_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = ZEROMATRIX(3,3);


   return threeLoop * beta_TYd;
}

/**
 * Calculates the 4-loop beta function of TYd.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_TYd_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = ZEROMATRIX(3,3);


   return fourLoop * beta_TYd;
}

/**
 * Calculates the 5-loop beta function of TYd.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_soft_parameters::calc_beta_TYd_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = ZEROMATRIX(3,3);


   return fiveLoop * beta_TYd;
}

} // namespace flexiblesusy
