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

// File generated at Fri 8 Jan 2016 15:12:48

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
 * Calculates the one-loop beta function of mvR2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mvR2_one_loop(const Soft_traces& soft_traces) const
{
   const auto Qv = INPUT(Qv);
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_mvR2;

   beta_mvR2 = (oneOver16PiSqr*(4*mHu2*(Yv.transpose()*Yv.conjugate()) +
      4*((TYv).transpose()*TYv.conjugate()) + 2*(mvR2*Yv.transpose()*
      Yv.conjugate()) + 4*(Yv.transpose()*ml2*Yv.conjugate()) + 2*(Yv.transpose
      ()*Yv.conjugate()*mvR2) + 1.5491933384829668*g1*Tr11*UNITMATRIX(3) + 2*gp
      *Qv*Tr14*UNITMATRIX(3) - 4.8*AbsSqr(MassB)*Sqr(g1)*UNITMATRIX(3) - 8*
      AbsSqr(MassU)*Sqr(gp)*Sqr(Qv)*UNITMATRIX(3))).real();


   return beta_mvR2;
}

/**
 * Calculates the two-loop beta function of mvR2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mvR2_two_loop(const Soft_traces& soft_traces) const
{
   const auto Qv = INPUT(Qv);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto Qq = INPUT(Qq);
   const auto Qu = INPUT(Qu);
   const auto Qs = INPUT(Qs);
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
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_mvR2;

   const Eigen::Matrix<double,3,3> beta_mvR2_1 = (-0.08*twoLoop*
      UNITMATRIX(3)*(10*(15*traceconjTYuTpTYu + 5*traceconjTYvTpTYv + 15*
      tracemq2AdjYuYu + 15*tracemu2YuAdjYu + 30*mHu2*traceYuAdjYu + 10*mHu2*
      traceYvAdjYv + 5*traceYvAdjYvconjml2 + 5*traceYvconjmvR2AdjYv + 5*mHd2*
      AbsSqr(Lambdax) + 10*mHu2*AbsSqr(Lambdax) + 5*ms2*AbsSqr(Lambdax) + 5*
      AbsSqr(TLambdax) + 3*mHu2*Sqr(g1) + 6*AbsSqr(MassB)*Sqr(g1) - 15*mHu2*Sqr
      (g2) - 30*AbsSqr(MassWB)*Sqr(g2) - 10*mHu2*Sqr(gp)*Sqr(QHu) - 10*mHu2*Sqr
      (gp)*Sqr(Ql) - 20*AbsSqr(MassU)*Sqr(gp)*(Sqr(QHu) + Sqr(Ql) - Sqr(Qv)) +
      10*mHu2*Sqr(gp)*Sqr(Qv))*(Yv.transpose()*Yv.conjugate()) - 10*(3*MassB*
      Sqr(g1) - 5*(3*traceAdjYuTYu + traceAdjYvTYv + 3*MassWB*Sqr(g2) + 2*MassU
      *Sqr(gp)*(Sqr(QHu) + Sqr(Ql) - Sqr(Qv))) - 5*Conj(Lambdax)*TLambdax)*(
      Yv.transpose()*TYv.conjugate()) + 150*traceconjTYuTpYu*((TYv).transpose()
      *Yv.conjugate()) + 50*traceconjTYvTpYv*((TYv).transpose()*Yv.conjugate())
      + 50*Conj(TLambdax)*Lambdax*((TYv).transpose()*Yv.conjugate()) - 30*Conj
      (MassB)*Sqr(g1)*((TYv).transpose()*Yv.conjugate()) + 150*Conj(MassWB)*Sqr
      (g2)*((TYv).transpose()*Yv.conjugate()) + 100*Conj(MassU)*Sqr(gp)*Sqr(QHu
      )*((TYv).transpose()*Yv.conjugate()) + 100*Conj(MassU)*Sqr(gp)*Sqr(Ql)*((
      TYv).transpose()*Yv.conjugate()) - 100*Conj(MassU)*Sqr(gp)*Sqr(Qv)*((TYv)
      .transpose()*Yv.conjugate()) + 150*traceYuAdjYu*((TYv).transpose()*
      TYv.conjugate()) + 50*traceYvAdjYv*((TYv).transpose()*TYv.conjugate()) +
      50*AbsSqr(Lambdax)*((TYv).transpose()*TYv.conjugate()) + 30*Sqr(g1)*((TYv
      ).transpose()*TYv.conjugate()) - 150*Sqr(g2)*((TYv).transpose()*
      TYv.conjugate()) - 100*Sqr(gp)*Sqr(QHu)*((TYv).transpose()*TYv.conjugate(
      )) - 100*Sqr(gp)*Sqr(Ql)*((TYv).transpose()*TYv.conjugate()) + 100*Sqr(gp
      )*Sqr(Qv)*((TYv).transpose()*TYv.conjugate()) + 75*traceYuAdjYu*(mvR2*
      Yv.transpose()*Yv.conjugate()) + 25*traceYvAdjYv*(mvR2*Yv.transpose()*
      Yv.conjugate()) + 25*AbsSqr(Lambdax)*(mvR2*Yv.transpose()*Yv.conjugate())
      + 15*Sqr(g1)*(mvR2*Yv.transpose()*Yv.conjugate()) - 75*Sqr(g2)*(mvR2*
      Yv.transpose()*Yv.conjugate()) - 50*Sqr(gp)*Sqr(QHu)*(mvR2*Yv.transpose()
      *Yv.conjugate()) - 50*Sqr(gp)*Sqr(Ql)*(mvR2*Yv.transpose()*Yv.conjugate()
      ) + 50*Sqr(gp)*Sqr(Qv)*(mvR2*Yv.transpose()*Yv.conjugate()) + 150*
      traceYuAdjYu*(Yv.transpose()*ml2*Yv.conjugate()) + 50*traceYvAdjYv*(
      Yv.transpose()*ml2*Yv.conjugate()) + 50*AbsSqr(Lambdax)*(Yv.transpose()*
      ml2*Yv.conjugate()) + 30*Sqr(g1)*(Yv.transpose()*ml2*Yv.conjugate()) -
      150*Sqr(g2)*(Yv.transpose()*ml2*Yv.conjugate()) - 100*Sqr(gp)*Sqr(QHu)*(
      Yv.transpose()*ml2*Yv.conjugate()) - 100*Sqr(gp)*Sqr(Ql)*(Yv.transpose()*
      ml2*Yv.conjugate()) + 100*Sqr(gp)*Sqr(Qv)*(Yv.transpose()*ml2*
      Yv.conjugate()) + 75*traceYuAdjYu*(Yv.transpose()*Yv.conjugate()*mvR2) +
      25*traceYvAdjYv*(Yv.transpose()*Yv.conjugate()*mvR2) + 25*AbsSqr(Lambdax)
      *(Yv.transpose()*Yv.conjugate()*mvR2) + 15*Sqr(g1)*(Yv.transpose()*
      Yv.conjugate()*mvR2) - 75*Sqr(g2)*(Yv.transpose()*Yv.conjugate()*mvR2) -
      50*Sqr(gp)*Sqr(QHu)*(Yv.transpose()*Yv.conjugate()*mvR2) - 50*Sqr(gp)*Sqr
      (Ql)*(Yv.transpose()*Yv.conjugate()*mvR2) + 50*Sqr(gp)*Sqr(Qv)*(
      Yv.transpose()*Yv.conjugate()*mvR2) + 50*mHd2*(Yv.transpose()*Ye.adjoint(
      )*Ye*Yv.conjugate()) + 50*mHu2*(Yv.transpose()*Ye.adjoint()*Ye*
      Yv.conjugate()) + 50*(Yv.transpose()*Ye.adjoint()*TYe*TYv.conjugate()) +
      50*(Yv.transpose()*(TYe).adjoint()*TYe*Yv.conjugate()) + 100*mHu2*(
      Yv.transpose()*Yv.conjugate()*Yv.transpose()*Yv.conjugate()) + 50*(
      Yv.transpose()*Yv.conjugate()*(TYv).transpose()*TYv.conjugate()) + 50*(
      Yv.transpose()*TYv.conjugate()*(TYv).transpose()*Yv.conjugate()) + 50*((
      TYv).transpose()*Ye.adjoint()*Ye*TYv.conjugate()) + 50*((TYv).transpose()
      *(TYe).adjoint()*Ye*Yv.conjugate()) + 50*((TYv).transpose()*Yv.conjugate(
      )*Yv.transpose()*TYv.conjugate()) + 50*((TYv).transpose()*TYv.conjugate()
      *Yv.transpose()*Yv.conjugate()) + 25*(mvR2*Yv.transpose()*Ye.adjoint()*Ye
      *Yv.conjugate()) + 25*(mvR2*Yv.transpose()*Yv.conjugate()*Yv.transpose()*
      Yv.conjugate()) + 50*(Yv.transpose()*ml2*Ye.adjoint()*Ye*Yv.conjugate())
      + 50*(Yv.transpose()*ml2*Yv.conjugate()*Yv.transpose()*Yv.conjugate()) +
      50*(Yv.transpose()*Ye.adjoint()*me2*Ye*Yv.conjugate()) + 50*(Yv.transpose
      ()*Ye.adjoint()*Ye*ml2*Yv.conjugate()) + 25*(Yv.transpose()*Ye.adjoint()*
      Ye*Yv.conjugate()*mvR2) + 50*(Yv.transpose()*Yv.conjugate()*mvR2*
      Yv.transpose()*Yv.conjugate()) + 50*(Yv.transpose()*Yv.conjugate()*
      Yv.transpose()*ml2*Yv.conjugate()) + 25*(Yv.transpose()*Yv.conjugate()*
      Yv.transpose()*Yv.conjugate()*mvR2) - 77.45966692414834*g1*gp*Qv*Tr2U114*
      UNITMATRIX(3) - 77.45966692414834*g1*gp*Qv*Tr2U141*UNITMATRIX(3) -
      77.45966692414834*g1*Tr31*UNITMATRIX(3) - 100*gp*Qv*Tr34*UNITMATRIX(3) -
      1728*Power(g1,4)*AbsSqr(MassB)*UNITMATRIX(3) - 60*Tr2U111*Sqr(g1)*
      UNITMATRIX(3) - 360*Qd*Qv*AbsSqr(MassB)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) -
      180*MassU*Qd*Qv*Conj(MassB)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) - 100*Tr2U144*
      Sqr(gp)*Sqr(Qv)*UNITMATRIX(3))).real();
   const Eigen::Matrix<double,3,3> beta_mvR2_2 = (4.8*Qv*twoLoop*Sqr(gp)*
      ((2*MassB + MassU)*(3*Qe - QHd + QHu - 3*Ql + 3*Qq - 6*Qu + 5*Qv)*Conj(
      MassB)*Sqr(g1) + Conj(MassU)*((MassB + 2*MassU)*(3*Qd + 3*Qe - QHd + QHu
      - 3*Ql + 3*Qq - 6*Qu + 5*Qv)*Sqr(g1) + 5*MassU*Qv*Sqr(gp)*(9*Sqr(Qd) + 3*
      Sqr(Qe) + 2*Sqr(QHd) + 2*Sqr(QHu) + 6*Sqr(Ql) + 18*Sqr(Qq) + Sqr(Qs) + 9*
      Sqr(Qu) + 5*Sqr(Qv))))*UNITMATRIX(3)).real();

   beta_mvR2 = beta_mvR2_1 + beta_mvR2_2;


   return beta_mvR2;
}

/**
 * Calculates the three-loop beta function of mvR2.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mvR2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mvR2;

   beta_mvR2 = ZEROMATRIX(3,3);


   return beta_mvR2;
}

} // namespace flexiblesusy
