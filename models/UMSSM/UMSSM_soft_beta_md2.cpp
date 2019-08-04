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

// File generated at Sun 4 Aug 2019 19:34:47

#include "UMSSM_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of md2.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_md2_1_loop(const Soft_traces& soft_traces) const
{
   const auto Qd = INPUT(Qd);
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = (oneOver16PiSqr*(4*mHd2*(Yd*Yd.adjoint()) + 4*(TYd*(TYd).adjoint(
      )) + 2*(md2*Yd*Yd.adjoint()) + 4*(Yd*mq2*Yd.adjoint()) + 2*(Yd*Yd.adjoint
      ()*md2) + 0.13333333333333333*(3.872983346207417*g1*Tr11 + 15*gp*Qd*Tr14
      - 4*AbsSqr(MassB)*Sqr(g1) - 80*AbsSqr(MassG)*Sqr(g3) - 60*AbsSqr(MassU)*
      Sqr(gp)*Sqr(Qd))*UNITMATRIX(3))).real();


   return beta_md2;
}

/**
 * Calculates the 2-loop beta function of md2.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_md2_2_loop(const Soft_traces& soft_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qu = INPUT(Qu);
   const auto Qs = INPUT(Qs);
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

   const Eigen::Matrix<double,3,3> beta_md2_1 = (UNITMATRIX(3)*(0.8*twoLoop*(-
      15*traceconjTYdTpTYd - 5*traceconjTYeTpTYe - 15*tracemd2YdAdjYd - 5*
      traceme2YeAdjYe - 5*traceml2AdjYeYe - 15*tracemq2AdjYdYd - 30*mHd2*
      traceYdAdjYd - 10*mHd2*traceYeAdjYe - 10*mHd2*AbsSqr(Lambdax) - 5*mHu2*
      AbsSqr(Lambdax) - 5*ms2*AbsSqr(Lambdax) - 5*AbsSqr(TLambdax) + mHd2*Sqr(
      g1) + 2*AbsSqr(MassB)*Sqr(g1) + 15*mHd2*Sqr(g2) + 30*AbsSqr(MassWB)*Sqr(
      g2) - 10*mHd2*Sqr(gp)*Sqr(Qd) - 20*AbsSqr(MassU)*Sqr(gp)*Sqr(Qd) + 10*
      mHd2*Sqr(gp)*Sqr(QHd) + 20*AbsSqr(MassU)*Sqr(gp)*Sqr(QHd) + 10*mHd2*Sqr(
      gp)*Sqr(Qq) + 20*AbsSqr(MassU)*Sqr(gp)*Sqr(Qq))*(Yd*Yd.adjoint()) - 0.8*
      twoLoop*(15*traceAdjYdTYd + 5*traceAdjYeTYe + MassB*Sqr(g1) + 15*MassWB*
      Sqr(g2) - 10*MassU*Sqr(gp)*Sqr(Qd) + 10*MassU*Sqr(gp)*Sqr(QHd) + 10*MassU
      *Sqr(gp)*Sqr(Qq) + 5*Conj(Lambdax)*TLambdax)*(Yd*(TYd).adjoint()) - 0.8*
      twoLoop*(15*traceconjTYdTpYd + 5*traceconjTYeTpYe + 5*Conj(TLambdax)*
      Lambdax + Conj(MassB)*Sqr(g1) + 15*Conj(MassWB)*Sqr(g2) - 10*Conj(MassU)*
      Sqr(gp)*Sqr(Qd) + 10*Conj(MassU)*Sqr(gp)*Sqr(QHd) + 10*Conj(MassU)*Sqr(gp
      )*Sqr(Qq))*(TYd*Yd.adjoint()) + 0.8*twoLoop*(-15*traceYdAdjYd - 5*
      traceYeAdjYe - 5*AbsSqr(Lambdax) + Sqr(g1) + 15*Sqr(g2) - 10*Sqr(gp)*Sqr(
      Qd) + 10*Sqr(gp)*Sqr(QHd) + 10*Sqr(gp)*Sqr(Qq))*(TYd*(TYd).adjoint()) +
      0.4*twoLoop*(-15*traceYdAdjYd - 5*traceYeAdjYe - 5*AbsSqr(Lambdax) + Sqr(
      g1) + 15*Sqr(g2) - 10*Sqr(gp)*Sqr(Qd) + 10*Sqr(gp)*Sqr(QHd) + 10*Sqr(gp)*
      Sqr(Qq))*(md2*Yd*Yd.adjoint()) + 0.8*twoLoop*(-15*traceYdAdjYd - 5*
      traceYeAdjYe - 5*AbsSqr(Lambdax) + Sqr(g1) + 15*Sqr(g2) - 10*Sqr(gp)*Sqr(
      Qd) + 10*Sqr(gp)*Sqr(QHd) + 10*Sqr(gp)*Sqr(Qq))*(Yd*mq2*Yd.adjoint()) +
      0.4*twoLoop*(-15*traceYdAdjYd - 5*traceYeAdjYe - 5*AbsSqr(Lambdax) + Sqr(
      g1) + 15*Sqr(g2) - 10*Sqr(gp)*Sqr(Qd) + 10*Sqr(gp)*Sqr(QHd) + 10*Sqr(gp)*
      Sqr(Qq))*(Yd*Yd.adjoint()*md2) - 8*mHd2*twoLoop*(Yd*Yd.adjoint()*Yd*Yd.
      adjoint()) - 4*twoLoop*(Yd*Yd.adjoint()*TYd*(TYd).adjoint()) - 4*(mHd2 +
      mHu2)*twoLoop*(Yd*Yu.adjoint()*Yu*Yd.adjoint()) - 4*twoLoop*(Yd*Yu.
      adjoint()*TYu*(TYd).adjoint()) - 4*twoLoop*(Yd*(TYd).adjoint()*TYd*Yd.
      adjoint()) - 4*twoLoop*(Yd*(TYu).adjoint()*TYu*Yd.adjoint()) - 4*twoLoop*
      (TYd*Yd.adjoint()*Yd*(TYd).adjoint()) - 4*twoLoop*(TYd*Yu.adjoint()*Yu*(
      TYd).adjoint()) - 4*twoLoop*(TYd*(TYd).adjoint()*Yd*Yd.adjoint()) - 4*
      twoLoop*(TYd*(TYu).adjoint()*Yu*Yd.adjoint()) - 2*twoLoop*(md2*Yd*Yd.
      adjoint()*Yd*Yd.adjoint()) - 2*twoLoop*(md2*Yd*Yu.adjoint()*Yu*Yd.adjoint
      ()) - 4*twoLoop*(Yd*mq2*Yd.adjoint()*Yd*Yd.adjoint()) - 4*twoLoop*(Yd*mq2
      *Yu.adjoint()*Yu*Yd.adjoint()) - 4*twoLoop*(Yd*Yd.adjoint()*md2*Yd*Yd.
      adjoint()) - 4*twoLoop*(Yd*Yd.adjoint()*Yd*mq2*Yd.adjoint()) - 2*twoLoop*
      (Yd*Yd.adjoint()*Yd*Yd.adjoint()*md2) - 4*twoLoop*(Yd*Yu.adjoint()*mu2*Yu
      *Yd.adjoint()) - 4*twoLoop*(Yd*Yu.adjoint()*Yu*mq2*Yd.adjoint()) - 2*
      twoLoop*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*md2) + 0.035555555555555556*
      twoLoop*(58.09475019311125*g1*gp*Qd*Tr2U114 + 58.09475019311125*g1*gp*Qd*
      Tr2U141 + 58.09475019311125*g1*Tr31 + 225*gp*Qd*Tr34 + 303*AbsSqr(MassB)*
      Quad(g1) + 300*Tr23*Quad(g3) + 15*Tr2U111*Sqr(g1) + 80*AbsSqr(MassB)*Sqr(
      g1)*Sqr(g3) + 225*Tr2U144*Sqr(gp)*Sqr(Qd))*UNITMATRIX(3))).real();
   const Eigen::Matrix<double,3,3> beta_md2_2 = (0.17777777777777778*twoLoop*(-
      240*AbsSqr(MassG)*Quad(g3) + 1485*AbsSqr(MassU)*Quad(gp)*Quad(Qd) + 16*
      AbsSqr(MassG)*Sqr(g1)*Sqr(g3) + 8*MassG*Conj(MassB)*Sqr(g1)*Sqr(g3) + 8*
      MassB*Conj(MassG)*Sqr(g1)*Sqr(g3) + 54*Qd*Qe*AbsSqr(MassB)*Sqr(g1)*Sqr(gp
      ) - 18*Qd*QHd*AbsSqr(MassB)*Sqr(g1)*Sqr(gp) + 18*Qd*QHu*AbsSqr(MassB)*Sqr
      (g1)*Sqr(gp) - 54*Qd*Ql*AbsSqr(MassB)*Sqr(g1)*Sqr(gp) + 54*Qd*Qq*AbsSqr(
      MassB)*Sqr(g1)*Sqr(gp) - 108*Qd*Qu*AbsSqr(MassB)*Sqr(g1)*Sqr(gp) + 54*Qd*
      Qe*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) - 18*Qd*QHd*AbsSqr(MassU)*Sqr(g1)*Sqr(gp
      ) + 18*Qd*QHu*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) - 54*Qd*Ql*AbsSqr(MassU)*Sqr(
      g1)*Sqr(gp) + 54*Qd*Qq*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) - 108*Qd*Qu*AbsSqr(
      MassU)*Sqr(g1)*Sqr(gp) + 27*MassU*Qd*Qe*Conj(MassB)*Sqr(g1)*Sqr(gp) - 9*
      MassU*Qd*QHd*Conj(MassB)*Sqr(g1)*Sqr(gp) + 9*MassU*Qd*QHu*Conj(MassB)*Sqr
      (g1)*Sqr(gp) - 27*MassU*Qd*Ql*Conj(MassB)*Sqr(g1)*Sqr(gp) + 27*MassU*Qd*
      Qq*Conj(MassB)*Sqr(g1)*Sqr(gp) - 54*MassU*Qd*Qu*Conj(MassB)*Sqr(g1)*Sqr(
      gp) + 27*MassB*Qd*Qe*Conj(MassU)*Sqr(g1)*Sqr(gp) - 9*MassB*Qd*QHd*Conj(
      MassU)*Sqr(g1)*Sqr(gp) + 9*MassB*Qd*QHu*Conj(MassU)*Sqr(g1)*Sqr(gp) - 27*
      MassB*Qd*Ql*Conj(MassU)*Sqr(g1)*Sqr(gp) + 27*MassB*Qd*Qq*Conj(MassU)*Sqr(
      g1)*Sqr(gp) - 54*MassB*Qd*Qu*Conj(MassU)*Sqr(g1)*Sqr(gp) + 66*AbsSqr(
      MassB)*Sqr(g1)*Sqr(gp)*Sqr(Qd) + 66*AbsSqr(MassU)*Sqr(g1)*Sqr(gp)*Sqr(Qd)
      + 33*MassU*Conj(MassB)*Sqr(g1)*Sqr(gp)*Sqr(Qd) + 33*MassB*Conj(MassU)*Sqr
      (g1)*Sqr(gp)*Sqr(Qd) + 240*AbsSqr(MassG)*Sqr(g3)*Sqr(gp)*Sqr(Qd) + 240*
      AbsSqr(MassU)*Sqr(g3)*Sqr(gp)*Sqr(Qd) + 120*MassU*Conj(MassG)*Sqr(g3)*Sqr
      (gp)*Sqr(Qd) + 120*MassG*Conj(MassU)*Sqr(g3)*Sqr(gp)*Sqr(Qd) + 405*AbsSqr
      (MassU)*Quad(gp)*Sqr(Qd)*Sqr(Qe) + 270*AbsSqr(MassU)*Quad(gp)*Sqr(Qd)*Sqr
      (QHd) + 270*AbsSqr(MassU)*Quad(gp)*Sqr(Qd)*Sqr(QHu) + 810*AbsSqr(MassU)*
      Quad(gp)*Sqr(Qd)*Sqr(Ql) + 2430*AbsSqr(MassU)*Quad(gp)*Sqr(Qd)*Sqr(Qq) +
      135*AbsSqr(MassU)*Quad(gp)*Sqr(Qd)*Sqr(Qs) + 1215*AbsSqr(MassU)*Quad(gp)*
      Sqr(Qd)*Sqr(Qu) + 405*AbsSqr(MassU)*Quad(gp)*Sqr(Qd)*Sqr(Qv))*UNITMATRIX(
      3)).real();

   beta_md2 = beta_md2_1 + beta_md2_2;


   return beta_md2;
}

/**
 * Calculates the 3-loop beta function of md2.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_md2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = ZEROMATRIX(3,3);


   return beta_md2;
}

/**
 * Calculates the 4-loop beta function of md2.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_md2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = ZEROMATRIX(3,3);


   return beta_md2;
}

/**
 * Calculates the 5-loop beta function of md2.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_md2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = ZEROMATRIX(3,3);


   return beta_md2;
}

} // namespace flexiblesusy
