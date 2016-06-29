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

// File generated at Wed 29 Jun 2016 12:03:44

#include "UMSSM_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of md2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_md2_one_loop(const Soft_traces& soft_traces) const
{
   const auto Qd = INPUT(Qd);
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = (oneOver16PiSqr*(4*mHd2*(Yd*Yd.adjoint()) + 4*(TYd*(TYd)
      .adjoint()) + 2*(md2*Yd*Yd.adjoint()) + 4*(Yd*mq2*Yd.adjoint()) + 2*(Yd*
      Yd.adjoint()*md2) + 0.5163977794943222*g1*Tr11*UNITMATRIX(3) + 2*gp*Qd*
      Tr14*UNITMATRIX(3) - 0.5333333333333333*AbsSqr(MassB)*Sqr(g1)*UNITMATRIX(
      3) - 10.666666666666666*AbsSqr(MassG)*Sqr(g3)*UNITMATRIX(3) - 8*AbsSqr(
      MassU)*Sqr(gp)*Sqr(Qd)*UNITMATRIX(3))).real();


   return beta_md2;
}

/**
 * Calculates the two-loop beta function of md2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_md2_two_loop(const Soft_traces& soft_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qs = INPUT(Qs);
   const auto Qu = INPUT(Qu);
   const auto Qv = INPUT(Qv);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr23 = TRACE_STRUCT.Tr23;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_md2;

   const Eigen::Matrix<double,3,3> beta_md2_1 = (-0.008888888888888889*
      twoLoop*UNITMATRIX(3)*(-90*(-15*traceconjTYdTpTYd - 5*traceconjTYeTpTYe -
      15*tracemd2YdAdjYd - 5*traceme2YeAdjYe - 5*traceml2AdjYeYe - 15*
      tracemq2AdjYdYd - 30*mHd2*traceYdAdjYd - 10*mHd2*traceYeAdjYe - 10*mHd2*
      AbsSqr(Lambdax) - 5*mHu2*AbsSqr(Lambdax) - 5*ms2*AbsSqr(Lambdax) - 5*
      AbsSqr(TLambdax) + mHd2*Sqr(g1) + 2*AbsSqr(MassB)*Sqr(g1) + 15*mHd2*Sqr(
      g2) + 30*AbsSqr(MassWB)*Sqr(g2) - 10*mHd2*Sqr(gp)*Sqr(Qd) + 10*mHd2*Sqr(
      gp)*Sqr(QHd) - 20*AbsSqr(MassU)*Sqr(gp)*(Sqr(Qd) - Sqr(QHd) - Sqr(Qq)) +
      10*mHd2*Sqr(gp)*Sqr(Qq))*(Yd*Yd.adjoint()) + 90*(MassB*Sqr(g1) + 5*(3*
      traceAdjYdTYd + traceAdjYeTYe + 3*MassWB*Sqr(g2) + 2*MassU*Sqr(gp)*(-Sqr(
      Qd) + Sqr(QHd) + Sqr(Qq))) + 5*Conj(Lambdax)*TLambdax)*(Yd*(TYd).adjoint(
      )) + 1350*traceconjTYdTpYd*(TYd*Yd.adjoint()) + 450*traceconjTYeTpYe*(TYd
      *Yd.adjoint()) + 450*Conj(TLambdax)*Lambdax*(TYd*Yd.adjoint()) + 90*Conj(
      MassB)*Sqr(g1)*(TYd*Yd.adjoint()) + 1350*Conj(MassWB)*Sqr(g2)*(TYd*
      Yd.adjoint()) - 900*Conj(MassU)*Sqr(gp)*Sqr(Qd)*(TYd*Yd.adjoint()) + 900*
      Conj(MassU)*Sqr(gp)*Sqr(QHd)*(TYd*Yd.adjoint()) + 900*Conj(MassU)*Sqr(gp)
      *Sqr(Qq)*(TYd*Yd.adjoint()) + 1350*traceYdAdjYd*(TYd*(TYd).adjoint()) +
      450*traceYeAdjYe*(TYd*(TYd).adjoint()) + 450*AbsSqr(Lambdax)*(TYd*(TYd)
      .adjoint()) - 90*Sqr(g1)*(TYd*(TYd).adjoint()) - 1350*Sqr(g2)*(TYd*(TYd)
      .adjoint()) + 900*Sqr(gp)*Sqr(Qd)*(TYd*(TYd).adjoint()) - 900*Sqr(gp)*Sqr
      (QHd)*(TYd*(TYd).adjoint()) - 900*Sqr(gp)*Sqr(Qq)*(TYd*(TYd).adjoint()) +
      675*traceYdAdjYd*(md2*Yd*Yd.adjoint()) + 225*traceYeAdjYe*(md2*Yd*
      Yd.adjoint()) + 225*AbsSqr(Lambdax)*(md2*Yd*Yd.adjoint()) - 45*Sqr(g1)*(
      md2*Yd*Yd.adjoint()) - 675*Sqr(g2)*(md2*Yd*Yd.adjoint()) + 450*Sqr(gp)*
      Sqr(Qd)*(md2*Yd*Yd.adjoint()) - 450*Sqr(gp)*Sqr(QHd)*(md2*Yd*Yd.adjoint()
      ) - 450*Sqr(gp)*Sqr(Qq)*(md2*Yd*Yd.adjoint()) + 1350*traceYdAdjYd*(Yd*mq2
      *Yd.adjoint()) + 450*traceYeAdjYe*(Yd*mq2*Yd.adjoint()) + 450*AbsSqr(
      Lambdax)*(Yd*mq2*Yd.adjoint()) - 90*Sqr(g1)*(Yd*mq2*Yd.adjoint()) - 1350*
      Sqr(g2)*(Yd*mq2*Yd.adjoint()) + 900*Sqr(gp)*Sqr(Qd)*(Yd*mq2*Yd.adjoint())
      - 900*Sqr(gp)*Sqr(QHd)*(Yd*mq2*Yd.adjoint()) - 900*Sqr(gp)*Sqr(Qq)*(Yd*
      mq2*Yd.adjoint()) + 675*traceYdAdjYd*(Yd*Yd.adjoint()*md2) + 225*
      traceYeAdjYe*(Yd*Yd.adjoint()*md2) + 225*AbsSqr(Lambdax)*(Yd*Yd.adjoint()
      *md2) - 45*Sqr(g1)*(Yd*Yd.adjoint()*md2) - 675*Sqr(g2)*(Yd*Yd.adjoint()*
      md2) + 450*Sqr(gp)*Sqr(Qd)*(Yd*Yd.adjoint()*md2) - 450*Sqr(gp)*Sqr(QHd)*(
      Yd*Yd.adjoint()*md2) - 450*Sqr(gp)*Sqr(Qq)*(Yd*Yd.adjoint()*md2) + 900*
      mHd2*(Yd*Yd.adjoint()*Yd*Yd.adjoint()) + 450*(Yd*Yd.adjoint()*TYd*(TYd)
      .adjoint()) + 450*mHd2*(Yd*Yu.adjoint()*Yu*Yd.adjoint()) + 450*mHu2*(Yd*
      Yu.adjoint()*Yu*Yd.adjoint()) + 450*(Yd*Yu.adjoint()*TYu*(TYd).adjoint())
      + 450*(Yd*(TYd).adjoint()*TYd*Yd.adjoint()) + 450*(Yd*(TYu).adjoint()*
      TYu*Yd.adjoint()) + 450*(TYd*Yd.adjoint()*Yd*(TYd).adjoint()) + 450*(TYd*
      Yu.adjoint()*Yu*(TYd).adjoint()) + 450*(TYd*(TYd).adjoint()*Yd*Yd.adjoint
      ()) + 450*(TYd*(TYu).adjoint()*Yu*Yd.adjoint()) + 225*(md2*Yd*Yd.adjoint(
      )*Yd*Yd.adjoint()) + 225*(md2*Yd*Yu.adjoint()*Yu*Yd.adjoint()) + 450*(Yd*
      mq2*Yd.adjoint()*Yd*Yd.adjoint()) + 450*(Yd*mq2*Yu.adjoint()*Yu*
      Yd.adjoint()) + 450*(Yd*Yd.adjoint()*md2*Yd*Yd.adjoint()) + 450*(Yd*
      Yd.adjoint()*Yd*mq2*Yd.adjoint()) + 225*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*
      md2) + 450*(Yd*Yu.adjoint()*mu2*Yu*Yd.adjoint()) + 450*(Yd*Yu.adjoint()*
      Yu*mq2*Yd.adjoint()) + 225*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*md2) - 1200*
      Power(g3,4)*Tr23*UNITMATRIX(3) - 232.379000772445*g1*gp*Qd*Tr2U114*
      UNITMATRIX(3) - 232.379000772445*g1*gp*Qd*Tr2U141*UNITMATRIX(3) -
      232.379000772445*g1*Tr31*UNITMATRIX(3) - 900*gp*Qd*Tr34*UNITMATRIX(3) -
      1212*Power(g1,4)*AbsSqr(MassB)*UNITMATRIX(3) - 60*Tr2U111*Sqr(g1)*
      UNITMATRIX(3) - 320*AbsSqr(MassB)*Sqr(g1)*Sqr(g3)*UNITMATRIX(3) - 900*
      Tr2U144*Sqr(gp)*Sqr(Qd)*UNITMATRIX(3))).real();
   const Eigen::Matrix<double,3,3> beta_md2_2 = (0.17777777777777778*
      twoLoop*(Conj(MassB)*Sqr(g1)*(8*MassG*Sqr(g3) + 3*(2*MassB + MassU)*Qd*(
      11*Qd + 3*(3*Qe - QHd + QHu - 3*Ql + 3*Qq - 6*Qu))*Sqr(gp)) + 8*Conj(
      MassG)*Sqr(g3)*((MassB + 2*MassG)*Sqr(g1) + 15*(-2*MassG*Sqr(g3) + (2*
      MassG + MassU)*Sqr(gp)*Sqr(Qd))) + 3*Qd*Conj(MassU)*Sqr(gp)*((MassB + 2*
      MassU)*(11*Qd + 3*(3*Qe - QHd + QHu - 3*Ql + 3*Qq - 6*Qu))*Sqr(g1) + 5*Qd
      *(8*(MassG + 2*MassU)*Sqr(g3) + 9*MassU*Sqr(gp)*(11*Sqr(Qd) + 3*Sqr(Qe) +
      2*Sqr(QHd) + 2*Sqr(QHu) + 6*Sqr(Ql) + 18*Sqr(Qq) + Sqr(Qs) + 9*Sqr(Qu) +
      3*Sqr(Qv)))))*UNITMATRIX(3)).real();

   beta_md2 = beta_md2_1 + beta_md2_2;


   return beta_md2;
}

/**
 * Calculates the three-loop beta function of md2.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_md2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = ZEROMATRIX(3,3);


   return beta_md2;
}

} // namespace flexiblesusy
