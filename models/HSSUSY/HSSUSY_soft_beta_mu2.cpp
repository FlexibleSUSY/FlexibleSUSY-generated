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

// File generated at Wed 16 Oct 2019 21:34:00

#include "HSSUSY_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of mu2.
 *
 * @return 1-loop beta function
 */
double HSSUSY_soft_parameters::calc_beta_mu2_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_mu2;

   beta_mu2 = Re(-0.1*mu2*oneOver16PiSqr*(-60*traceYdAdjYd - 20*traceYeAdjYe -
      60*traceYuAdjYu - 60*Lambdax + 9*Sqr(g1) + 45*Sqr(g2)));


   return beta_mu2;
}

/**
 * Calculates the 2-loop beta function of mu2.
 *
 * @return 2-loop beta function
 */
double HSSUSY_soft_parameters::calc_beta_mu2_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_mu2;

   beta_mu2 = Re(0.0025*mu2*twoLoop*(-5400*traceYdAdjYdYdAdjYd - 8400*
      traceYdAdjYuYuAdjYd - 1800*traceYeAdjYeYeAdjYe - 5400*traceYuAdjYuYuAdjYu
       - 14400*traceYdAdjYd*Lambdax - 4800*traceYeAdjYe*Lambdax - 14400*
      traceYuAdjYu*Lambdax + 1671*Quad(g1) - 3625*Quad(g2) + 500*traceYdAdjYd*
      Sqr(g1) + 1500*traceYeAdjYe*Sqr(g1) + 1700*traceYuAdjYu*Sqr(g1) + 2880*
      Lambdax*Sqr(g1) + 4500*traceYdAdjYd*Sqr(g2) + 1500*traceYeAdjYe*Sqr(g2) +
      4500*traceYuAdjYu*Sqr(g2) + 14400*Lambdax*Sqr(g2) + 450*Sqr(g1)*Sqr(g2) +
      16000*traceYdAdjYd*Sqr(g3) + 16000*traceYuAdjYu*Sqr(g3) - 6000*Sqr(
      Lambdax)));


   return beta_mu2;
}

/**
 * Calculates the 3-loop beta function of mu2.
 *
 * @return 3-loop beta function
 */
double HSSUSY_soft_parameters::calc_beta_mu2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mu2;

   beta_mu2 = Re(16.756629209125986*mu2*threeLoop*(15.30737457986509*Cube(
      Lambdax) + 1.*Power6(g1) + 36.01244566362596*Power6(g2) +
      18.429131386172795*Power6(Yd(2,2)) + 0.4139679292425927*Power6(Ye(2,2)) +
      18.429131386172795*Power6(Yu(2,2)) - 3.927137010147313*Lambdax*Quad(g1) -
      3.850086285654542*Lambdax*Quad(g2) + 20.73175259456566*Lambdax*Quad(Yd(2,
      2)) + 9.05898764799418*Lambdax*Quad(Ye(2,2)) + 20.73175259456566*Lambdax*
      Quad(Yu(2,2)) + 1.1852954889988687*Quad(g2)*Sqr(g1) + 1.4222858389686859*
      Quad(Yd(2,2))*Sqr(g1) - 1.7965342456882272*Quad(Ye(2,2))*Sqr(g1) -
      0.8960859852603394*Quad(Yu(2,2))*Sqr(g1) + 1.167061165496437*Quad(g1)*Sqr
      (g2) - 0.4570467771716553*Quad(Yd(2,2))*Sqr(g2) + 0.3847519367275216*Quad
      (Ye(2,2))*Sqr(g2) - 0.4570467771716553*Quad(Yu(2,2))*Sqr(g2) -
      2.2572004440579008*Lambdax*Sqr(g1)*Sqr(g2) - 0.5001699678078098*Quad(g1)*
      Sqr(g3) - 3.410249780507793*Quad(g2)*Sqr(g3) - 24.974054450461615*Quad(Yd
      (2,2))*Sqr(g3) - 24.974054450461615*Quad(Yu(2,2))*Sqr(g3) -
      1.1507215629992067*Sqr(g1)*Sqr(Lambdax) - 5.753607814996034*Sqr(g2)*Sqr(
      Lambdax) - 1.4129200768537742*Quad(g1)*Sqr(Yd(2,2)) - 12.249064902749408*
      Quad(g2)*Sqr(Yd(2,2)) + 21.30308732502087*Quad(g3)*Sqr(Yd(2,2)) +
      8.593613799222506*Quad(Ye(2,2))*Sqr(Yd(2,2)) + 35.35455070551097*Quad(Yu(
      2,2))*Sqr(Yd(2,2)) - 3.9225017102152724*Lambdax*Sqr(g1)*Sqr(Yd(2,2)) -
      19.012859486691834*Lambdax*Sqr(g2)*Sqr(Yd(2,2)) - 0.9679005125426392*Sqr(
      g1)*Sqr(g2)*Sqr(Yd(2,2)) + 4.7971925150761985*Lambdax*Sqr(g3)*Sqr(Yd(2,2)
      ) - 0.24969029303485452*Sqr(g1)*Sqr(g3)*Sqr(Yd(2,2)) + 0.903779089067895*
      Sqr(g2)*Sqr(g3)*Sqr(Yd(2,2)) + 4.431082115224105*Sqr(Lambdax)*Sqr(Yd(2,2)
      ) - 2.4518611783899384*Quad(g1)*Sqr(Ye(2,2)) - 4.083021634249803*Quad(g2)
      *Sqr(Ye(2,2)) + 8.593613799222506*Quad(Yd(2,2))*Sqr(Ye(2,2)) +
      8.593613799222506*Quad(Yu(2,2))*Sqr(Ye(2,2)) - 0.9077345271487414*Lambdax
      *Sqr(g1)*Sqr(Ye(2,2)) - 6.337619828897279*Lambdax*Sqr(g2)*Sqr(Ye(2,2)) +
      1.2654587331364668*Sqr(g1)*Sqr(g2)*Sqr(Ye(2,2)) + 1.4770273717413684*Sqr(
      Lambdax)*Sqr(Ye(2,2)) - 6.445210349416881*Lambdax*Sqr(Yd(2,2))*Sqr(Ye(2,2
      )) - 0.32226051747084405*Sqr(g1)*Sqr(Yd(2,2))*Sqr(Ye(2,2)) -
      1.6113025873542206*Sqr(g2)*Sqr(Yd(2,2))*Sqr(Ye(2,2)) - 3.3087112184463003
      *Quad(g1)*Sqr(Yu(2,2)) - 12.249064902749408*Quad(g2)*Sqr(Yu(2,2)) +
      21.30308732502087*Quad(g3)*Sqr(Yu(2,2)) + 35.35455070551097*Quad(Yd(2,2))
      *Sqr(Yu(2,2)) + 8.593613799222506*Quad(Ye(2,2))*Sqr(Yu(2,2)) -
      3.5627122715845574*Lambdax*Sqr(g1)*Sqr(Yu(2,2)) - 19.012859486691834*
      Lambdax*Sqr(g2)*Sqr(Yu(2,2)) + 1.369048888979986*Sqr(g1)*Sqr(g2)*Sqr(Yu(2
      ,2)) + 4.7971925150761985*Lambdax*Sqr(g3)*Sqr(Yu(2,2)) +
      1.0416480395104468*Sqr(g1)*Sqr(g3)*Sqr(Yu(2,2)) + 0.903779089067895*Sqr(
      g2)*Sqr(g3)*Sqr(Yu(2,2)) + 4.431082115224105*Sqr(Lambdax)*Sqr(Yu(2,2)) -
      24.89428427856406*Lambdax*Sqr(Yd(2,2))*Sqr(Yu(2,2)) - 2.858120668905674*
      Sqr(g1)*Sqr(Yd(2,2))*Sqr(Yu(2,2)) - 7.499185617963107*Sqr(g2)*Sqr(Yd(2,2)
      )*Sqr(Yu(2,2)) - 1.9930895579602697*Sqr(g3)*Sqr(Yd(2,2))*Sqr(Yu(2,2)) -
      6.445210349416881*Lambdax*Sqr(Ye(2,2))*Sqr(Yu(2,2)) - 0.32226051747084405
      *Sqr(g1)*Sqr(Ye(2,2))*Sqr(Yu(2,2)) - 1.6113025873542206*Sqr(g2)*Sqr(Ye(2,
      2))*Sqr(Yu(2,2)) + 1.253235345719949*Sqr(Yd(2,2))*Sqr(Ye(2,2))*Sqr(Yu(2,2
      ))));


   return beta_mu2;
}

/**
 * Calculates the 4-loop beta function of mu2.
 *
 * @return 4-loop beta function
 */
double HSSUSY_soft_parameters::calc_beta_mu2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mu2;

   beta_mu2 = 0;


   return beta_mu2;
}

/**
 * Calculates the 5-loop beta function of mu2.
 *
 * @return 5-loop beta function
 */
double HSSUSY_soft_parameters::calc_beta_mu2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mu2;

   beta_mu2 = 0;


   return beta_mu2;
}

} // namespace flexiblesusy
