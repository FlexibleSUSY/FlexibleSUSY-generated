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

// File generated at Mon 9 May 2016 12:03:51

#include "SM_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of mu2.
 *
 * @return one-loop beta function
 */
double SM_soft_parameters::calc_beta_mu2_one_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_mu2;

   beta_mu2 = Re(oneOver16PiSqr*(6*mu2*traceYdAdjYd + 2*mu2*traceYeAdjYe
      + 6*mu2*traceYuAdjYu + 6*mu2*Lambdax - 0.9*mu2*Sqr(g1) - 4.5*mu2*Sqr(g2))
      );


   return beta_mu2;
}

/**
 * Calculates the two-loop beta function of mu2.
 *
 * @return two-loop beta function
 */
double SM_soft_parameters::calc_beta_mu2_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_mu2;

   beta_mu2 = Re(0.0025*mu2*twoLoop*(1671*Power(g1,4) - 3625*Power(g2,4)
      - 5400*traceYdAdjYdYdAdjYd - 8400*traceYdAdjYuYuAdjYd - 1800*
      traceYeAdjYeYeAdjYe - 5400*traceYuAdjYuYuAdjYu - 14400*traceYuAdjYu*
      Lambdax + 1700*traceYuAdjYu*Sqr(g1) + 2880*Lambdax*Sqr(g1) + 4500*
      traceYuAdjYu*Sqr(g2) + 14400*Lambdax*Sqr(g2) + 450*Sqr(g1)*Sqr(g2) + 300*
      traceYeAdjYe*(-16*Lambdax + 5*Sqr(g1) + 5*Sqr(g2)) + 16000*traceYuAdjYu*
      Sqr(g3) + 100*traceYdAdjYd*(-144*Lambdax + 5*Sqr(g1) + 45*Sqr(g2) + 160*
      Sqr(g3)) - 6000*Sqr(Lambdax)));


   return beta_mu2;
}

/**
 * Calculates the three-loop beta function of mu2.
 *
 * @return three-loop beta function
 */
double SM_soft_parameters::calc_beta_mu2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mu2;

   beta_mu2 = Re(2*mu2*threeLoop*(8.378314604562993*Power(g1,6) +
      301.7235994495886*Power(g2,6) + 128.25*Power(Lambdax,3) +
      178.48396765750311*Power(g3,4)*Sqr(Yd(2,2)) + 40.19238810996319*Lambdax*
      Sqr(g3)*Sqr(Yd(2,2)) + 37.125*Sqr(Lambdax)*Sqr(Yd(2,2)) + 12.375*Sqr(
      Lambdax)*Sqr(Ye(2,2)) - 54.*Lambdax*Sqr(Yd(2,2))*Sqr(Ye(2,2)) + Power(g2,
      4)*(-32.25723415592714*Lambdax - 28.572145541236182*Sqr(g3) -
      102.62651936694535*Sqr(Yd(2,2)) - 34.208839788981784*Sqr(Ye(2,2)) -
      102.62651936694535*Sqr(Yu(2,2))) + Power(g1,4)*(-32.90278936623708*
      Lambdax + 9.778005607297105*Sqr(g2) - 4.190581346047974*Sqr(g3) -
      11.837888914984243*Sqr(Yd(2,2)) - 20.542464319265452*Sqr(Ye(2,2)) -
      27.721423523790055*Sqr(Yu(2,2))) + 178.48396765750311*Power(g3,4)*Sqr(Yu(
      2,2)) + 40.19238810996319*Lambdax*Sqr(g3)*Sqr(Yu(2,2)) + 37.125*Sqr(
      Lambdax)*Sqr(Yu(2,2)) - 208.57214554123618*Lambdax*Sqr(Yd(2,2))*Sqr(Yu(2,
      2)) - 16.698731351660527*Sqr(g3)*Sqr(Yd(2,2))*Sqr(Yu(2,2)) - 54.*Lambdax*
      Sqr(Ye(2,2))*Sqr(Yu(2,2)) + 10.5*Sqr(Yd(2,2))*Sqr(Ye(2,2))*Sqr(Yu(2,2)) +
      173.69714554123618*Lambdax*Power(Yd(2,2),4) - 209.24048513745396*Sqr(g3)
      *Power(Yd(2,2),4) + 72.*Sqr(Ye(2,2))*Power(Yd(2,2),4) + 296.2115485137454
      *Sqr(Yu(2,2))*Power(Yd(2,2),4) + 154.40506064218175*Power(Yd(2,2),6) +
      75.89904851374538*Lambdax*Power(Ye(2,2),4) + 72.*Sqr(Yd(2,2))*Power(Ye(2,
      2),4) + 72.*Sqr(Yu(2,2))*Power(Ye(2,2),4) + 3.4683535473939138*Power(Ye(2
      ,2),6) + 173.69714554123618*Lambdax*Power(Yu(2,2),4) - 209.24048513745396
      *Sqr(g3)*Power(Yu(2,2),4) + 296.2115485137454*Sqr(Yd(2,2))*Power(Yu(2,2),
      4) + 72.*Sqr(Ye(2,2))*Power(Yu(2,2),4) + 154.40506064218175*Power(Yu(2,2)
      ,6) + Sqr(g1)*(9.930778506201856*Power(g2,4) - 9.641107277061808*Sqr(
      Lambdax) - 32.86395336511993*Lambdax*Sqr(Yd(2,2)) - 7.605285445876382*
      Lambdax*Sqr(Ye(2,2)) - 2.7*Sqr(Yd(2,2))*Sqr(Ye(2,2)) - 29.849524256872694
      *Lambdax*Sqr(Yu(2,2)) - 23.946234141895758*Sqr(Yd(2,2))*Sqr(Yu(2,2)) -
      2.7*Sqr(Ye(2,2))*Sqr(Yu(2,2)) + Sqr(g3)*(-2.091983828751535*Sqr(Yd(2,2))
      + 8.727254982244787*Sqr(Yu(2,2))) + Sqr(g2)*(-18.911535445876382*Lambdax
      - 8.109375*Sqr(Yd(2,2)) + 10.602411385309043*Sqr(Ye(2,2)) +
      11.470322300901756*Sqr(Yu(2,2))) + 11.916358216494471*Power(Yd(2,2),4) -
      15.051929108247235*Power(Ye(2,2),4) - 7.507690297250921*Power(Yu(2,2),4))
      + Sqr(g2)*(-48.205536385309045*Sqr(Lambdax) - 53.09857277061809*Lambdax*
      Sqr(Ye(2,2)) + Sqr(Yd(2,2))*(-159.29571831185424*Lambdax -
      13.500000000000004*Sqr(Ye(2,2)) - 62.83053638530905*Sqr(Yu(2,2))) -
      159.29571831185424*Lambdax*Sqr(Yu(2,2)) - 13.500000000000004*Sqr(Ye(2,2))
      *Sqr(Yu(2,2)) + Sqr(g3)*(7.572145541236182*Sqr(Yd(2,2)) +
      7.572145541236182*Sqr(Yu(2,2))) - 3.829281688145728*Power(Yd(2,2),4) +
      3.223572770618091*Power(Ye(2,2),4) - 3.829281688145728*Power(Yu(2,2),4)))
      );


   return beta_mu2;
}

} // namespace flexiblesusy
