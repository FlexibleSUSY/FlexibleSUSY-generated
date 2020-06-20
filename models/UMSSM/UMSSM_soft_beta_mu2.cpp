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
 * Calculates the 1-loop beta function of mu2.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mu2_1_loop(const Soft_traces& soft_traces) const
{
   const auto Qu = INPUT(Qu);
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = (4*mHu2*(Yu*Yu.adjoint()) + 4*(TYu*(TYu).adjoint()) + 2*(mu2*Yu*
      Yu.adjoint()) + 4*(Yu*mq2*Yu.adjoint()) + 2*(Yu*Yu.adjoint()*mu2) -
      0.13333333333333333*(7.745966692414834*g1*Tr11 - 15*gp*Qu*Tr14 + 16*
      AbsSqr(MassB)*Sqr(g1) + 80*AbsSqr(MassG)*Sqr(g3) + 60*AbsSqr(MassU)*Sqr(
      gp)*Sqr(Qu))*UNITMATRIX(3)).real();


   return oneLoop * beta_mu2;
}

/**
 * Calculates the 2-loop beta function of mu2.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mu2_2_loop(const Soft_traces& soft_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto Qu = INPUT(Qu);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qs = INPUT(Qs);
   const auto Qv = INPUT(Qv);
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double traceconjTYvTpTYv = TRACE_STRUCT.traceconjTYvTpTYv;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceYvAdjYvconjml2 = TRACE_STRUCT.traceYvAdjYvconjml2;
   const double traceYvconjmvR2AdjYv = TRACE_STRUCT.traceYvconjmvR2AdjYv;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYvTpYv = TRACE_STRUCT.traceconjTYvTpYv;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr23 = TRACE_STRUCT.Tr23;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_mu2;

   const Eigen::Matrix<double,3,3> beta_mu2_1 = (UNITMATRIX(3)*(-0.8*(15*
      traceconjTYuTpTYu + 5*traceconjTYvTpTYv + 15*tracemq2AdjYuYu + 15*
      tracemu2YuAdjYu + 30*mHu2*traceYuAdjYu + 10*mHu2*traceYvAdjYv + 5*
      traceYvAdjYvconjml2 + 5*traceYvconjmvR2AdjYv + 5*mHd2*AbsSqr(Lambdax) +
      10*mHu2*AbsSqr(Lambdax) + 5*ms2*AbsSqr(Lambdax) + 5*AbsSqr(TLambdax) +
      mHu2*Sqr(g1) + 2*AbsSqr(MassB)*Sqr(g1) - 15*mHu2*Sqr(g2) - 30*AbsSqr(
      MassWB)*Sqr(g2) - 10*mHu2*Sqr(gp)*Sqr(QHu) - 20*AbsSqr(MassU)*Sqr(gp)*Sqr
      (QHu) - 10*mHu2*Sqr(gp)*Sqr(Qq) - 20*AbsSqr(MassU)*Sqr(gp)*Sqr(Qq) + 10*
      mHu2*Sqr(gp)*Sqr(Qu) + 20*AbsSqr(MassU)*Sqr(gp)*Sqr(Qu))*(Yu*Yu.adjoint()
      ) + 0.8*(-15*traceAdjYuTYu - 5*traceAdjYvTYv + MassB*Sqr(g1) - 15*MassWB*
      Sqr(g2) - 10*MassU*Sqr(gp)*Sqr(QHu) - 10*MassU*Sqr(gp)*Sqr(Qq) + 10*MassU
      *Sqr(gp)*Sqr(Qu) - 5*Conj(Lambdax)*TLambdax)*(Yu*(TYu).adjoint()) - 0.8*(
      15*traceconjTYuTpYu + 5*traceconjTYvTpYv + 5*Conj(TLambdax)*Lambdax -
      Conj(MassB)*Sqr(g1) + 15*Conj(MassWB)*Sqr(g2) + 10*Conj(MassU)*Sqr(gp)*
      Sqr(QHu) + 10*Conj(MassU)*Sqr(gp)*Sqr(Qq) - 10*Conj(MassU)*Sqr(gp)*Sqr(Qu
      ))*(TYu*Yu.adjoint()) - 0.8*(15*traceYuAdjYu + 5*traceYvAdjYv + 5*AbsSqr(
      Lambdax) + Sqr(g1) - 15*Sqr(g2) - 10*Sqr(gp)*Sqr(QHu) - 10*Sqr(gp)*Sqr(Qq
      ) + 10*Sqr(gp)*Sqr(Qu))*(TYu*(TYu).adjoint()) - 0.4*(15*traceYuAdjYu + 5*
      traceYvAdjYv + 5*AbsSqr(Lambdax) + Sqr(g1) - 15*Sqr(g2) - 10*Sqr(gp)*Sqr(
      QHu) - 10*Sqr(gp)*Sqr(Qq) + 10*Sqr(gp)*Sqr(Qu))*(mu2*Yu*Yu.adjoint()) -
      0.8*(15*traceYuAdjYu + 5*traceYvAdjYv + 5*AbsSqr(Lambdax) + Sqr(g1) - 15*
      Sqr(g2) - 10*Sqr(gp)*Sqr(QHu) - 10*Sqr(gp)*Sqr(Qq) + 10*Sqr(gp)*Sqr(Qu))*
      (Yu*mq2*Yu.adjoint()) - 0.4*(15*traceYuAdjYu + 5*traceYvAdjYv + 5*AbsSqr(
      Lambdax) + Sqr(g1) - 15*Sqr(g2) - 10*Sqr(gp)*Sqr(QHu) - 10*Sqr(gp)*Sqr(Qq
      ) + 10*Sqr(gp)*Sqr(Qu))*(Yu*Yu.adjoint()*mu2) - 4*(mHd2 + mHu2)*(Yu*Yd.
      adjoint()*Yd*Yu.adjoint()) - 4*(Yu*Yd.adjoint()*TYd*(TYu).adjoint()) - 8*
      mHu2*(Yu*Yu.adjoint()*Yu*Yu.adjoint()) - 4*(Yu*Yu.adjoint()*TYu*(TYu).
      adjoint()) - 4*(Yu*(TYd).adjoint()*TYd*Yu.adjoint()) - 4*(Yu*(TYu).
      adjoint()*TYu*Yu.adjoint()) - 4*(TYu*Yd.adjoint()*Yd*(TYu).adjoint()) - 4
      *(TYu*Yu.adjoint()*Yu*(TYu).adjoint()) - 4*(TYu*(TYd).adjoint()*Yd*Yu.
      adjoint()) - 4*(TYu*(TYu).adjoint()*Yu*Yu.adjoint()) - 2*(mu2*Yu*Yd.
      adjoint()*Yd*Yu.adjoint()) - 2*(mu2*Yu*Yu.adjoint()*Yu*Yu.adjoint()) - 4*
      (Yu*mq2*Yd.adjoint()*Yd*Yu.adjoint()) - 4*(Yu*mq2*Yu.adjoint()*Yu*Yu.
      adjoint()) - 4*(Yu*Yd.adjoint()*md2*Yd*Yu.adjoint()) - 4*(Yu*Yd.adjoint()
      *Yd*mq2*Yu.adjoint()) - 2*(Yu*Yd.adjoint()*Yd*Yu.adjoint()*mu2) - 4*(Yu*
      Yu.adjoint()*mu2*Yu*Yu.adjoint()) - 4*(Yu*Yu.adjoint()*Yu*mq2*Yu.adjoint(
      )) - 2*(Yu*Yu.adjoint()*Yu*Yu.adjoint()*mu2) + 0.035555555555555556*(-
      116.1895003862225*g1*gp*Qu*Tr2U114 - 116.1895003862225*g1*gp*Qu*Tr2U141 -
      116.1895003862225*g1*Tr31 + 225*gp*Qu*Tr34 + 1284*AbsSqr(MassB)*Quad(g1)
      + 300*Tr23*Quad(g3) + 60*Tr2U111*Sqr(g1) + 320*AbsSqr(MassB)*Sqr(g1)*Sqr(
      g3) + 225*Tr2U144*Sqr(gp)*Sqr(Qu))*UNITMATRIX(3))).real();
   const Eigen::Matrix<double,3,3> beta_mu2_2 = (0.17777777777777778*(-240*
      AbsSqr(MassG)*Quad(g3) + 1485*AbsSqr(MassU)*Quad(gp)*Quad(Qu) + 64*AbsSqr
      (MassG)*Sqr(g1)*Sqr(g3) + 32*MassG*Conj(MassB)*Sqr(g1)*Sqr(g3) + 32*MassB
      *Conj(MassG)*Sqr(g1)*Sqr(g3) - 108*Qd*Qu*AbsSqr(MassB)*Sqr(g1)*Sqr(gp) -
      108*Qe*Qu*AbsSqr(MassB)*Sqr(g1)*Sqr(gp) + 36*QHd*Qu*AbsSqr(MassB)*Sqr(g1)
      *Sqr(gp) - 36*QHu*Qu*AbsSqr(MassB)*Sqr(g1)*Sqr(gp) + 108*Ql*Qu*AbsSqr(
      MassB)*Sqr(g1)*Sqr(gp) - 108*Qq*Qu*AbsSqr(MassB)*Sqr(g1)*Sqr(gp) - 108*Qd
      *Qu*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) - 108*Qe*Qu*AbsSqr(MassU)*Sqr(g1)*Sqr(
      gp) + 36*QHd*Qu*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) - 36*QHu*Qu*AbsSqr(MassU)*
      Sqr(g1)*Sqr(gp) + 108*Ql*Qu*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) - 108*Qq*Qu*
      AbsSqr(MassU)*Sqr(g1)*Sqr(gp) - 54*MassU*Qd*Qu*Conj(MassB)*Sqr(g1)*Sqr(gp
      ) - 54*MassU*Qe*Qu*Conj(MassB)*Sqr(g1)*Sqr(gp) + 18*MassU*QHd*Qu*Conj(
      MassB)*Sqr(g1)*Sqr(gp) - 18*MassU*QHu*Qu*Conj(MassB)*Sqr(g1)*Sqr(gp) + 54
      *MassU*Ql*Qu*Conj(MassB)*Sqr(g1)*Sqr(gp) - 54*MassU*Qq*Qu*Conj(MassB)*Sqr
      (g1)*Sqr(gp) - 54*MassB*Qd*Qu*Conj(MassU)*Sqr(g1)*Sqr(gp) - 54*MassB*Qe*
      Qu*Conj(MassU)*Sqr(g1)*Sqr(gp) + 18*MassB*QHd*Qu*Conj(MassU)*Sqr(g1)*Sqr(
      gp) - 18*MassB*QHu*Qu*Conj(MassU)*Sqr(g1)*Sqr(gp) + 54*MassB*Ql*Qu*Conj(
      MassU)*Sqr(g1)*Sqr(gp) - 54*MassB*Qq*Qu*Conj(MassU)*Sqr(g1)*Sqr(gp) + 264
      *AbsSqr(MassB)*Sqr(g1)*Sqr(gp)*Sqr(Qu) + 264*AbsSqr(MassU)*Sqr(g1)*Sqr(gp
      )*Sqr(Qu) + 132*MassU*Conj(MassB)*Sqr(g1)*Sqr(gp)*Sqr(Qu) + 132*MassB*
      Conj(MassU)*Sqr(g1)*Sqr(gp)*Sqr(Qu) + 240*AbsSqr(MassG)*Sqr(g3)*Sqr(gp)*
      Sqr(Qu) + 240*AbsSqr(MassU)*Sqr(g3)*Sqr(gp)*Sqr(Qu) + 120*MassU*Conj(
      MassG)*Sqr(g3)*Sqr(gp)*Sqr(Qu) + 120*MassG*Conj(MassU)*Sqr(g3)*Sqr(gp)*
      Sqr(Qu) + 1215*AbsSqr(MassU)*Quad(gp)*Sqr(Qd)*Sqr(Qu) + 405*AbsSqr(MassU)
      *Quad(gp)*Sqr(Qe)*Sqr(Qu) + 270*AbsSqr(MassU)*Quad(gp)*Sqr(QHd)*Sqr(Qu) +
      270*AbsSqr(MassU)*Quad(gp)*Sqr(QHu)*Sqr(Qu) + 810*AbsSqr(MassU)*Quad(gp)*
      Sqr(Ql)*Sqr(Qu) + 2430*AbsSqr(MassU)*Quad(gp)*Sqr(Qq)*Sqr(Qu) + 135*
      AbsSqr(MassU)*Quad(gp)*Sqr(Qs)*Sqr(Qu) + 405*AbsSqr(MassU)*Quad(gp)*Sqr(
      Qu)*Sqr(Qv))*UNITMATRIX(3)).real();

   beta_mu2 = beta_mu2_1 + beta_mu2_2;


   return twoLoop * beta_mu2;
}

/**
 * Calculates the 3-loop beta function of mu2.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mu2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = ZEROMATRIX(3,3);


   return threeLoop * beta_mu2;
}

/**
 * Calculates the 4-loop beta function of mu2.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mu2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = ZEROMATRIX(3,3);


   return fourLoop * beta_mu2;
}

/**
 * Calculates the 5-loop beta function of mu2.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mu2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = ZEROMATRIX(3,3);


   return fiveLoop * beta_mu2;
}

} // namespace flexiblesusy
