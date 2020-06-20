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


#include "UMSSM_susy_parameters.hpp"
#include "wrappers.hpp"
#include "functors.hpp"

#include <iostream>

namespace flexiblesusy {

#define CLASSNAME UMSSM_susy_parameters
#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces
#define TRACE_STRUCT_TYPE Susy_traces
#define CALCULATE_TRACES(l) calc_susy_traces(l);

const int UMSSM_susy_parameters::numberOfParameters;

UMSSM_susy_parameters::UMSSM_susy_parameters(const UMSSM_input_parameters& input_)
   : input(input_)
{
   set_number_of_parameters(numberOfParameters);
}

UMSSM_susy_parameters::UMSSM_susy_parameters(
   double scale_, int loops_, int thresholds_,
   const UMSSM_input_parameters& input_
   , const Eigen::Matrix<double,3,3>& Yd_, const Eigen::Matrix<double,3,3>& Ye_,
   double Lambdax_, const Eigen::Matrix<double,3,3>& Yv_, const Eigen::Matrix<
   double,3,3>& Yu_, double g1_, double g2_, double g3_, double gp_, double vd_
   , double vu_, double vS_
)
   : Beta_function()
   , Yd(Yd_), Ye(Ye_), Lambdax(Lambdax_), Yv(Yv_), Yu(Yu_), g1(g1_), g2(g2_), g3(
   g3_), gp(gp_), vd(vd_), vu(vu_), vS(vS_)
   , input(input_)
{
   set_number_of_parameters(numberOfParameters);
   set_scale(scale_);
   set_loops(loops_);
   set_thresholds(thresholds_);
}

Eigen::ArrayXd UMSSM_susy_parameters::beta() const
{
   return calc_beta().get().unaryExpr(Chop<double>(get_zero_threshold()));
}

UMSSM_susy_parameters UMSSM_susy_parameters::calc_beta(int loops) const
{
   Eigen::Matrix<double,3,3> beta_Yd = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_Ye = Eigen::Matrix<double,3,3>::Zero();
   double beta_Lambdax = 0.;
   Eigen::Matrix<double,3,3> beta_Yv = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_Yu = Eigen::Matrix<double,3,3>::Zero();
   double beta_g1 = 0.;
   double beta_g2 = 0.;
   double beta_g3 = 0.;
   double beta_gp = 0.;
   double beta_vd = 0.;
   double beta_vu = 0.;
   double beta_vS = 0.;

   if (loops > 0) {
      const auto TRACE_STRUCT = CALCULATE_TRACES(loops);

      beta_Yd += calc_beta_Yd_1_loop(TRACE_STRUCT);
      beta_Ye += calc_beta_Ye_1_loop(TRACE_STRUCT);
      beta_Lambdax += calc_beta_Lambdax_1_loop(TRACE_STRUCT);
      beta_Yv += calc_beta_Yv_1_loop(TRACE_STRUCT);
      beta_Yu += calc_beta_Yu_1_loop(TRACE_STRUCT);
      beta_g1 += calc_beta_g1_1_loop(TRACE_STRUCT);
      beta_g2 += calc_beta_g2_1_loop(TRACE_STRUCT);
      beta_g3 += calc_beta_g3_1_loop(TRACE_STRUCT);
      beta_gp += calc_beta_gp_1_loop(TRACE_STRUCT);
      beta_vd += calc_beta_vd_1_loop(TRACE_STRUCT);
      beta_vu += calc_beta_vu_1_loop(TRACE_STRUCT);
      beta_vS += calc_beta_vS_1_loop(TRACE_STRUCT);

      if (loops > 1) {
         beta_Yd += calc_beta_Yd_2_loop(TRACE_STRUCT);
         beta_Ye += calc_beta_Ye_2_loop(TRACE_STRUCT);
         beta_Lambdax += calc_beta_Lambdax_2_loop(TRACE_STRUCT);
         beta_Yv += calc_beta_Yv_2_loop(TRACE_STRUCT);
         beta_Yu += calc_beta_Yu_2_loop(TRACE_STRUCT);
         beta_g1 += calc_beta_g1_2_loop(TRACE_STRUCT);
         beta_g2 += calc_beta_g2_2_loop(TRACE_STRUCT);
         beta_g3 += calc_beta_g3_2_loop(TRACE_STRUCT);
         beta_gp += calc_beta_gp_2_loop(TRACE_STRUCT);
         beta_vd += calc_beta_vd_2_loop(TRACE_STRUCT);
         beta_vu += calc_beta_vu_2_loop(TRACE_STRUCT);
         beta_vS += calc_beta_vS_2_loop(TRACE_STRUCT);

         if (loops > 2) {

            if (loops > 3) {

               if (loops > 4) {

               }
            }
         }
      }
   }


   return UMSSM_susy_parameters(get_scale(), loops, get_thresholds(), input,
                    beta_Yd, beta_Ye, beta_Lambdax, beta_Yv, beta_Yu, beta_g1, beta_g2, beta_g3, beta_gp, beta_vd, beta_vu, beta_vS);
}

UMSSM_susy_parameters UMSSM_susy_parameters::calc_beta() const
{
   return calc_beta(get_loops());
}

void UMSSM_susy_parameters::clear()
{
   reset();
   Yd = Eigen::Matrix<double,3,3>::Zero();
   Ye = Eigen::Matrix<double,3,3>::Zero();
   Lambdax = 0.;
   Yv = Eigen::Matrix<double,3,3>::Zero();
   Yu = Eigen::Matrix<double,3,3>::Zero();
   g1 = 0.;
   g2 = 0.;
   g3 = 0.;
   gp = 0.;
   vd = 0.;
   vu = 0.;
   vS = 0.;

}

Eigen::Matrix<double,3,3> CLASSNAME::get_SqSq() const
{
   Eigen::Matrix<double,3,3> anomDim;
   const auto Qq = INPUT(Qq);
   const auto QHu = INPUT(QHu);
   const auto Qu = INPUT(Qu);
   const auto Qd = INPUT(Qd);
   const auto QHd = INPUT(QHd);
   const auto Qe = INPUT(Qe);
   const auto Ql = INPUT(Ql);
   const auto Qs = INPUT(Qs);
   const auto Qv = INPUT(Qv);

   anomDim = (oneOver16PiSqr*(Yd.adjoint()*Yd + Yu.adjoint()*Yu -
      0.03333333333333333*(Sqr(g1) + 45*Sqr(g2) + 80*Sqr(g3) + 60*Sqr(gp)*Sqr(
      Qq))*UNITMATRIX(3))).real();

   if (get_loops() > 1) {
      anomDim += (twoLoop*(-(AbsSqr(Lambdax)*(Yu.adjoint()*Yu)) + 0.8*Sqr(g1)*(
         Yu.adjoint()*Yu) + 2*Sqr(gp)*Sqr(QHu)*(Yu.adjoint()*Yu) - 2*Sqr(gp)*
         Sqr(Qq)*(Yu.adjoint()*Yu) + 2*Sqr(gp)*Sqr(Qu)*(Yu.adjoint()*Yu) - 2*(
         Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 2*(Yu.adjoint()*Yu*Yu.adjoint()*Yu)
         + Yd.adjoint()*Yd*(-AbsSqr(Lambdax) + 0.4*Sqr(g1) + 2*Sqr(gp)*Sqr(Qd)
         + 2*Sqr(gp)*Sqr(QHd) - 2*Sqr(gp)*Sqr(Qq) - 3*(Yd*Yd.adjoint()).trace()
         - (Ye*Ye.adjoint()).trace()) - 3*(Yu.adjoint()*Yu)*(Yu*Yu.adjoint()).
         trace() - Yu.adjoint()*Yu*(Yv*Yv.adjoint()).trace() +
         0.0011111111111111111*(199*Quad(g1) + 10*Sqr(g1)*(9*Sqr(g2) + 16*Sqr(
         g3) + 12*Qq*(9*Qd + 9*Qe - 3*QHd + 3*QHu - 9*Ql + 10*Qq - 18*Qu)*Sqr(
         gp)) + 25*(135*Quad(g2) + 72*Sqr(g2)*(4*Sqr(g3) + 3*Sqr(gp)*Sqr(Qq)) +
         8*(-4*Quad(g3) + 48*Sqr(g3)*Sqr(gp)*Sqr(Qq) + 9*Quad(gp)*Sqr(Qq)*(9*
         Sqr(Qd) + 3*Sqr(Qe) + 2*Sqr(QHd) + 2*Sqr(QHu) + 6*Sqr(Ql) + 20*Sqr(Qq)
         + Sqr(Qs) + 9*Sqr(Qu) + 3*Sqr(Qv)))))*UNITMATRIX(3))).real();
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SlSl() const
{
   Eigen::Matrix<double,3,3> anomDim;
   const auto Ql = INPUT(Ql);
   const auto QHu = INPUT(QHu);
   const auto Qv = INPUT(Qv);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto Qd = INPUT(Qd);
   const auto Qq = INPUT(Qq);
   const auto Qu = INPUT(Qu);
   const auto Qs = INPUT(Qs);

   anomDim = (oneOver16PiSqr*(Ye.adjoint()*Ye + Yv.conjugate()*Yv.transpose() -
      0.1*(3*Sqr(g1) + 15*Sqr(g2) + 20*Sqr(gp)*Sqr(Ql))*UNITMATRIX(3))).real();

   if (get_loops() > 1) {
      anomDim += (twoLoop*(-(AbsSqr(Lambdax)*(Yv.conjugate()*Yv.transpose())) +
         2*Sqr(gp)*Sqr(QHu)*(Yv.conjugate()*Yv.transpose()) - 2*Sqr(gp)*Sqr(Ql)
         *(Yv.conjugate()*Yv.transpose()) + 2*Sqr(gp)*Sqr(Qv)*(Yv.conjugate()*
         Yv.transpose()) - 2*(Ye.adjoint()*Ye*Ye.adjoint()*Ye) - 2*(Yv.
         conjugate()*Yv.transpose()*Yv.conjugate()*Yv.transpose()) + Ye.adjoint
         ()*Ye*(-AbsSqr(Lambdax) + 1.2*Sqr(g1) + 2*Sqr(gp)*Sqr(Qe) + 2*Sqr(gp)*
         Sqr(QHd) - 2*Sqr(gp)*Sqr(Ql) - 3*(Yd*Yd.adjoint()).trace() - (Ye*Ye.
         adjoint()).trace()) - 3*(Yv.conjugate()*Yv.transpose())*(Yu*Yu.adjoint
         ()).trace() - Yv.conjugate()*Yv.transpose()*(Yv*Yv.adjoint()).trace()
         + 0.01*(207*Quad(g1) + 30*Sqr(g1)*(3*Sqr(g2) + 4*Ql*(-3*Qd - 3*Qe +
         QHd - QHu + 4*Ql - 3*Qq + 6*Qu)*Sqr(gp)) + 25*(15*Quad(g2) + 24*Sqr(g2
         )*Sqr(gp)*Sqr(Ql) + 8*Quad(gp)*Sqr(Ql)*(9*Sqr(Qd) + 3*Sqr(Qe) + 2*Sqr(
         QHd) + 2*Sqr(QHu) + 8*Sqr(Ql) + 18*Sqr(Qq) + Sqr(Qs) + 9*Sqr(Qu) + 3*
         Sqr(Qv))))*UNITMATRIX(3))).real();
   }

   return anomDim;
}

double CLASSNAME::get_SHdSHd() const
{
   double anomDim = 0;
   const auto QHd = INPUT(QHd);
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qs = INPUT(Qs);
   const auto Qu = INPUT(Qu);
   const auto Qv = INPUT(Qv);

   anomDim = Re(oneOver16PiSqr*(AbsSqr(Lambdax) - 0.3*Sqr(g1) - 1.5*Sqr(g2) - 2
      *Sqr(gp)*Sqr(QHd) + 3*(Yd*Yd.adjoint()).trace() + (Ye*Ye.adjoint()).trace
      ()));

   if (get_loops() > 1) {
      anomDim += Re(twoLoop*(2.07*Quad(g1) + 3.75*Quad(g2) + 8*Quad(gp)*Quad(
         QHd) + 0.9*Sqr(g1)*Sqr(g2) - 3.6*Qd*QHd*Sqr(g1)*Sqr(gp) - 3.6*Qe*QHd*
         Sqr(g1)*Sqr(gp) - 1.2*QHd*QHu*Sqr(g1)*Sqr(gp) + 3.6*QHd*Ql*Sqr(g1)*Sqr
         (gp) - 3.6*QHd*Qq*Sqr(g1)*Sqr(gp) + 7.2*QHd*Qu*Sqr(g1)*Sqr(gp) + 2.4*
         Sqr(g1)*Sqr(gp)*Sqr(QHd) + 6*Sqr(g2)*Sqr(gp)*Sqr(QHd) + 18*Quad(gp)*
         Sqr(Qd)*Sqr(QHd) + 6*Quad(gp)*Sqr(Qe)*Sqr(QHd) + 4*Quad(gp)*Sqr(QHd)*
         Sqr(QHu) + 12*Quad(gp)*Sqr(QHd)*Sqr(Ql) + 36*Quad(gp)*Sqr(QHd)*Sqr(Qq)
         + 2*Quad(gp)*Sqr(QHd)*Sqr(Qs) + 18*Quad(gp)*Sqr(QHd)*Sqr(Qu) + 6*Quad(
         gp)*Sqr(QHd)*Sqr(Qv) - 3*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 0.4*(Sqr(g1
         ) - 5*(8*Sqr(g3) + 3*Sqr(gp)*(Sqr(Qd) - Sqr(QHd) + Sqr(Qq))))*(Yd*Yd.
         adjoint()).trace() + 1.2*Sqr(g1)*(Ye*Ye.adjoint()).trace() + 2*Sqr(gp)
         *Sqr(Qe)*(Ye*Ye.adjoint()).trace() - 2*Sqr(gp)*Sqr(QHd)*(Ye*Ye.adjoint
         ()).trace() + 2*Sqr(gp)*Sqr(Ql)*(Ye*Ye.adjoint()).trace() - AbsSqr(
         Lambdax)*(2*Sqr(gp)*(Sqr(QHd) - Sqr(QHu) - Sqr(Qs)) + 3*(Yu*Yu.adjoint
         ()).trace() + (Yv*Yv.adjoint()).trace()) - 9*(Yd*Yd.adjoint()*Yd*Yd.
         adjoint()).trace() - 3*(Yd*Yu.adjoint()*Yu*Yd.adjoint()).trace() - 3*(
         Ye*Ye.adjoint()*Ye*Ye.adjoint()).trace() - (Yv*Yv.adjoint()*Ye.
         transpose()*Ye.conjugate()).trace()));
   }

   return anomDim;
}

double CLASSNAME::get_SHuSHu() const
{
   double anomDim = 0;
   const auto QHu = INPUT(QHu);
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qs = INPUT(Qs);
   const auto Qu = INPUT(Qu);
   const auto Qv = INPUT(Qv);

   anomDim = Re(oneOver16PiSqr*(AbsSqr(Lambdax) - 0.3*Sqr(g1) - 1.5*Sqr(g2) - 2
      *Sqr(gp)*Sqr(QHu) + 3*(Yu*Yu.adjoint()).trace() + (Yv*Yv.adjoint()).trace
      ()));

   if (get_loops() > 1) {
      anomDim += Re(twoLoop*(2.07*Quad(g1) + 3.75*Quad(g2) + 8*Quad(gp)*Quad(
         QHu) + 0.9*Sqr(g1)*Sqr(g2) + 3.6*Qd*QHu*Sqr(g1)*Sqr(gp) + 3.6*Qe*QHu*
         Sqr(g1)*Sqr(gp) - 1.2*QHd*QHu*Sqr(g1)*Sqr(gp) - 3.6*QHu*Ql*Sqr(g1)*Sqr
         (gp) + 3.6*QHu*Qq*Sqr(g1)*Sqr(gp) - 7.2*QHu*Qu*Sqr(g1)*Sqr(gp) + 2.4*
         Sqr(g1)*Sqr(gp)*Sqr(QHu) + 6*Sqr(g2)*Sqr(gp)*Sqr(QHu) + 18*Quad(gp)*
         Sqr(Qd)*Sqr(QHu) + 6*Quad(gp)*Sqr(Qe)*Sqr(QHu) + 4*Quad(gp)*Sqr(QHd)*
         Sqr(QHu) + 12*Quad(gp)*Sqr(QHu)*Sqr(Ql) + 36*Quad(gp)*Sqr(QHu)*Sqr(Qq)
         + 2*Quad(gp)*Sqr(QHu)*Sqr(Qs) + 18*Quad(gp)*Sqr(QHu)*Sqr(Qu) + 6*Quad(
         gp)*Sqr(QHu)*Sqr(Qv) - 3*Sqr(Conj(Lambdax))*Sqr(Lambdax) + AbsSqr(
         Lambdax)*(2*Sqr(gp)*(Sqr(QHd) - Sqr(QHu) + Sqr(Qs)) - 3*(Yd*Yd.adjoint
         ()).trace() - (Ye*Ye.adjoint()).trace()) + 0.4*(2*Sqr(g1) + 40*Sqr(g3)
         + 15*Sqr(gp)*(-Sqr(QHu) + Sqr(Qq) + Sqr(Qu)))*(Yu*Yu.adjoint()).trace(
         ) - 2*Sqr(gp)*Sqr(QHu)*(Yv*Yv.adjoint()).trace() + 2*Sqr(gp)*Sqr(Ql)*(
         Yv*Yv.adjoint()).trace() + 2*Sqr(gp)*Sqr(Qv)*(Yv*Yv.adjoint()).trace()
         - 3*(Yd*Yu.adjoint()*Yu*Yd.adjoint()).trace() - 9*(Yu*Yu.adjoint()*Yu*
         Yu.adjoint()).trace() - 3*(Yv*Yv.adjoint()*Yv*Yv.adjoint()).trace() -
         (Yv*Yv.adjoint()*Ye.transpose()*Ye.conjugate()).trace()));
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SdRSdR() const
{
   Eigen::Matrix<double,3,3> anomDim;
   const auto Qd = INPUT(Qd);
   const auto QHd = INPUT(QHd);
   const auto Qq = INPUT(Qq);
   const auto Qe = INPUT(Qe);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qu = INPUT(Qu);
   const auto Qs = INPUT(Qs);
   const auto Qv = INPUT(Qv);

   anomDim = (oneOver16PiSqr*(2*(Yd.conjugate()*Yd.transpose()) -
      0.13333333333333333*(Sqr(g1) + 20*Sqr(g3) + 15*Sqr(gp)*Sqr(Qd))*
      UNITMATRIX(3))).real();

   if (get_loops() > 1) {
      anomDim += (twoLoop*(-2*(Yd.conjugate()*Yd.transpose()*Yd.conjugate()*Yd.
         transpose() + Yd.conjugate()*Yu.transpose()*Yu.conjugate()*Yd.
         transpose()) + Yd.conjugate()*Yd.transpose()*(-2*AbsSqr(Lambdax) + 0.4
         *Sqr(g1) + 6*Sqr(g2) - 4*Sqr(gp)*Sqr(Qd) + 4*Sqr(gp)*Sqr(QHd) + 4*Sqr(
         gp)*Sqr(Qq) - 6*(Yd*Yd.adjoint()).trace() - 2*(Ye*Ye.adjoint()).trace(
         )) + 0.008888888888888889*(101*Quad(g1) + 10*Sqr(g1)*(8*Sqr(g3) + 3*Qd
         *(11*Qd + 9*Qe - 3*QHd + 3*QHu - 9*Ql + 9*Qq - 18*Qu)*Sqr(gp)) - 25*(4
         *Quad(g3) - 48*Sqr(g3)*Sqr(gp)*Sqr(Qd) - 9*Quad(gp)*Sqr(Qd)*(11*Sqr(Qd
         ) + 3*Sqr(Qe) + 2*Sqr(QHd) + 2*Sqr(QHu) + 6*Sqr(Ql) + 18*Sqr(Qq) + Sqr
         (Qs) + 9*Sqr(Qu) + 3*Sqr(Qv))))*UNITMATRIX(3))).real();
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SuRSuR() const
{
   Eigen::Matrix<double,3,3> anomDim;
   const auto Qu = INPUT(Qu);
   const auto QHu = INPUT(QHu);
   const auto Qq = INPUT(Qq);
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto Ql = INPUT(Ql);
   const auto Qs = INPUT(Qs);
   const auto Qv = INPUT(Qv);

   anomDim = (oneOver16PiSqr*(2*(Yu.conjugate()*Yu.transpose()) -
      0.13333333333333333*(4*Sqr(g1) + 20*Sqr(g3) + 15*Sqr(gp)*Sqr(Qu))*
      UNITMATRIX(3))).real();

   if (get_loops() > 1) {
      anomDim += (twoLoop*(-2*(Yu.conjugate()*Yd.transpose()*Yd.conjugate()*Yu.
         transpose() + Yu.conjugate()*Yu.transpose()*Yu.conjugate()*Yu.
         transpose()) + Yu.conjugate()*Yu.transpose()*(-2*AbsSqr(Lambdax) - 0.4
         *Sqr(g1) + 6*Sqr(g2) + 4*Sqr(gp)*Sqr(QHu) + 4*Sqr(gp)*Sqr(Qq) - 4*Sqr(
         gp)*Sqr(Qu) - 6*(Yu*Yu.adjoint()).trace() - 2*(Yv*Yv.adjoint()).trace(
         )) + 0.008888888888888889*(428*Quad(g1) + 20*Sqr(g1)*(16*Sqr(g3) - 3*(
         9*Qd + 9*Qe - 3*QHd + 3*QHu - 9*Ql + 9*Qq - 22*Qu)*Qu*Sqr(gp)) - 25*(4
         *Quad(g3) - 48*Sqr(g3)*Sqr(gp)*Sqr(Qu) - 9*Quad(gp)*Sqr(Qu)*(9*Sqr(Qd)
         + 3*Sqr(Qe) + 2*Sqr(QHd) + 2*Sqr(QHu) + 6*Sqr(Ql) + 18*Sqr(Qq) + Sqr(
         Qs) + 11*Sqr(Qu) + 3*Sqr(Qv))))*UNITMATRIX(3))).real();
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SeRSeR() const
{
   Eigen::Matrix<double,3,3> anomDim;
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto Ql = INPUT(Ql);
   const auto Qd = INPUT(Qd);
   const auto QHu = INPUT(QHu);
   const auto Qq = INPUT(Qq);
   const auto Qu = INPUT(Qu);
   const auto Qs = INPUT(Qs);
   const auto Qv = INPUT(Qv);

   anomDim = (oneOver16PiSqr*(2*(Ye.conjugate()*Ye.transpose()) - 0.4*(3*Sqr(g1
      ) + 5*Sqr(gp)*Sqr(Qe))*UNITMATRIX(3))).real();

   if (get_loops() > 1) {
      anomDim += (twoLoop*(-2*(Ye.conjugate()*Yv*Yv.adjoint()*Ye.transpose() +
         Ye.conjugate()*Ye.transpose()*Ye.conjugate()*Ye.transpose()) + Ye.
         conjugate()*Ye.transpose()*(-2*AbsSqr(Lambdax) - 1.2*Sqr(g1) + 6*Sqr(
         g2) - 4*Sqr(gp)*Sqr(Qe) + 4*Sqr(gp)*Sqr(QHd) + 4*Sqr(gp)*Sqr(Ql) - 6*(
         Yd*Yd.adjoint()).trace() - 2*(Ye*Ye.adjoint()).trace()) + 0.08*(117*
         Quad(g1) + 30*Qe*(3*Qd + 5*Qe - QHd + QHu - 3*Ql + 3*Qq - 6*Qu)*Sqr(g1
         )*Sqr(gp) + 25*Quad(gp)*Sqr(Qe)*(9*Sqr(Qd) + 5*Sqr(Qe) + 2*Sqr(QHd) +
         2*Sqr(QHu) + 6*Sqr(Ql) + 18*Sqr(Qq) + Sqr(Qs) + 9*Sqr(Qu) + 3*Sqr(Qv))
         )*UNITMATRIX(3))).real();
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SvRSvR() const
{
   Eigen::Matrix<double,3,3> anomDim;
   const auto Qv = INPUT(Qv);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto Qq = INPUT(Qq);
   const auto Qs = INPUT(Qs);
   const auto Qu = INPUT(Qu);

   anomDim = (oneOver16PiSqr*(2*(Yv.adjoint()*Yv) - 2*Sqr(gp)*Sqr(Qv)*
      UNITMATRIX(3))).real();

   if (get_loops() > 1) {
      anomDim += (twoLoop*(-2*(Yv.adjoint()*Yv*Yv.adjoint()*Yv + Yv.adjoint()*
         Ye.transpose()*Ye.conjugate()*Yv) + Yv.adjoint()*Yv*(-2*AbsSqr(Lambdax
         ) + 1.2*Sqr(g1) + 6*Sqr(g2) + 4*Sqr(gp)*Sqr(QHu) + 4*Sqr(gp)*Sqr(Ql) -
         4*Sqr(gp)*Sqr(Qv) - 6*(Yu*Yu.adjoint()).trace() - 2*(Yv*Yv.adjoint()).
         trace()) + 2*Quad(gp)*Sqr(Qv)*(9*Sqr(Qd) + 3*Sqr(Qe) + 2*Sqr(QHd) + 2*
         Sqr(QHu) + 6*Sqr(Ql) + 18*Sqr(Qq) + Sqr(Qs) + 9*Sqr(Qu) + 5*Sqr(Qv))*
         UNITMATRIX(3))).real();
   }

   return anomDim;
}

double CLASSNAME::get_SsRSsR() const
{
   double anomDim = 0;
   const auto Qs = INPUT(Qs);
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qu = INPUT(Qu);
   const auto Qv = INPUT(Qv);

   anomDim = Re(oneOver16PiSqr*(2*AbsSqr(Lambdax) - 2*Sqr(gp)*Sqr(Qs)));

   if (get_loops() > 1) {
      anomDim += Re(twoLoop*(2*Quad(gp)*Sqr(Qs)*(9*Sqr(Qd) + 3*Sqr(Qe) + 2*Sqr(
         QHd) + 2*Sqr(QHu) + 6*Sqr(Ql) + 18*Sqr(Qq) + 3*Sqr(Qs) + 9*Sqr(Qu) + 3
         *Sqr(Qv)) - 4*Sqr(Conj(Lambdax))*Sqr(Lambdax) + 0.4*AbsSqr(Lambdax)*(3
         *Sqr(g1) + 15*Sqr(g2) + 10*Sqr(gp)*Sqr(QHd) + 10*Sqr(gp)*Sqr(QHu) - 10
         *Sqr(gp)*Sqr(Qs) - 15*(Yd*Yd.adjoint()).trace() - 5*(Ye*Ye.adjoint()).
         trace() - 15*(Yu*Yu.adjoint()).trace() - 5*(Yv*Yv.adjoint()).trace()))
         );
   }

   return anomDim;
}



Eigen::ArrayXd UMSSM_susy_parameters::get() const
{
   Eigen::ArrayXd pars(numberOfParameters);

   pars(0) = Yd(0,0);
   pars(1) = Yd(0,1);
   pars(2) = Yd(0,2);
   pars(3) = Yd(1,0);
   pars(4) = Yd(1,1);
   pars(5) = Yd(1,2);
   pars(6) = Yd(2,0);
   pars(7) = Yd(2,1);
   pars(8) = Yd(2,2);
   pars(9) = Ye(0,0);
   pars(10) = Ye(0,1);
   pars(11) = Ye(0,2);
   pars(12) = Ye(1,0);
   pars(13) = Ye(1,1);
   pars(14) = Ye(1,2);
   pars(15) = Ye(2,0);
   pars(16) = Ye(2,1);
   pars(17) = Ye(2,2);
   pars(18) = Lambdax;
   pars(19) = Yv(0,0);
   pars(20) = Yv(0,1);
   pars(21) = Yv(0,2);
   pars(22) = Yv(1,0);
   pars(23) = Yv(1,1);
   pars(24) = Yv(1,2);
   pars(25) = Yv(2,0);
   pars(26) = Yv(2,1);
   pars(27) = Yv(2,2);
   pars(28) = Yu(0,0);
   pars(29) = Yu(0,1);
   pars(30) = Yu(0,2);
   pars(31) = Yu(1,0);
   pars(32) = Yu(1,1);
   pars(33) = Yu(1,2);
   pars(34) = Yu(2,0);
   pars(35) = Yu(2,1);
   pars(36) = Yu(2,2);
   pars(37) = g1;
   pars(38) = g2;
   pars(39) = g3;
   pars(40) = gp;
   pars(41) = vd;
   pars(42) = vu;
   pars(43) = vS;


   return pars;
}

void UMSSM_susy_parameters::print() const
{
   this->print(std::cerr);
}

void UMSSM_susy_parameters::print(std::ostream& ostr) const
{
   ostr << "susy parameters at Q = " << get_scale() << ":\n";
   ostr << "Yd = " << Yd << '\n';
   ostr << "Ye = " << Ye << '\n';
   ostr << "Lambdax = " << Lambdax << '\n';
   ostr << "Yv = " << Yv << '\n';
   ostr << "Yu = " << Yu << '\n';
   ostr << "g1 = " << g1 << '\n';
   ostr << "g2 = " << g2 << '\n';
   ostr << "g3 = " << g3 << '\n';
   ostr << "gp = " << gp << '\n';
   ostr << "vd = " << vd << '\n';
   ostr << "vu = " << vu << '\n';
   ostr << "vS = " << vS << '\n';

}

void UMSSM_susy_parameters::set(const Eigen::ArrayXd& pars)
{
   Yd(0,0) = pars(0);
   Yd(0,1) = pars(1);
   Yd(0,2) = pars(2);
   Yd(1,0) = pars(3);
   Yd(1,1) = pars(4);
   Yd(1,2) = pars(5);
   Yd(2,0) = pars(6);
   Yd(2,1) = pars(7);
   Yd(2,2) = pars(8);
   Ye(0,0) = pars(9);
   Ye(0,1) = pars(10);
   Ye(0,2) = pars(11);
   Ye(1,0) = pars(12);
   Ye(1,1) = pars(13);
   Ye(1,2) = pars(14);
   Ye(2,0) = pars(15);
   Ye(2,1) = pars(16);
   Ye(2,2) = pars(17);
   Lambdax = pars(18);
   Yv(0,0) = pars(19);
   Yv(0,1) = pars(20);
   Yv(0,2) = pars(21);
   Yv(1,0) = pars(22);
   Yv(1,1) = pars(23);
   Yv(1,2) = pars(24);
   Yv(2,0) = pars(25);
   Yv(2,1) = pars(26);
   Yv(2,2) = pars(27);
   Yu(0,0) = pars(28);
   Yu(0,1) = pars(29);
   Yu(0,2) = pars(30);
   Yu(1,0) = pars(31);
   Yu(1,1) = pars(32);
   Yu(1,2) = pars(33);
   Yu(2,0) = pars(34);
   Yu(2,1) = pars(35);
   Yu(2,2) = pars(36);
   g1 = pars(37);
   g2 = pars(38);
   g3 = pars(39);
   gp = pars(40);
   vd = pars(41);
   vu = pars(42);
   vS = pars(43);

}

const UMSSM_input_parameters& UMSSM_susy_parameters::get_input() const
{
   return input;
}

UMSSM_input_parameters& UMSSM_susy_parameters::get_input()
{
   return input;
}

void UMSSM_susy_parameters::set_input_parameters(const UMSSM_input_parameters& input_)
{
   input = input_;
}

UMSSM_susy_parameters::Susy_traces UMSSM_susy_parameters::calc_susy_traces(int loops) const
{
   Susy_traces susy_traces;

   if (loops > 0) {
      

      TRACE_STRUCT.traceYdAdjYd = Re((Yd*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYeAdjYe = Re((Ye*Ye.adjoint()).trace());
      TRACE_STRUCT.traceYuAdjYu = Re((Yu*Yu.adjoint()).trace());
      TRACE_STRUCT.traceYvAdjYv = Re((Yv*Yv.adjoint()).trace());

   }

   if (loops > 1) {
      TRACE_STRUCT.traceYdAdjYdYdAdjYd = Re((Yd*Yd.adjoint()*Yd*Yd.adjoint()).trace()
         );
      TRACE_STRUCT.traceYdAdjYuYuAdjYd = Re((Yd*Yu.adjoint()*Yu*Yd.adjoint()).trace()
         );
      TRACE_STRUCT.traceYeAdjYeYeAdjYe = Re((Ye*Ye.adjoint()*Ye*Ye.adjoint()).trace()
         );
      TRACE_STRUCT.traceYuAdjYuYuAdjYu = Re((Yu*Yu.adjoint()*Yu*Yu.adjoint()).trace()
         );
      TRACE_STRUCT.traceYvAdjYvYvAdjYv = Re((Yv*Yv.adjoint()*Yv*Yv.adjoint()).trace()
         );
      TRACE_STRUCT.traceYvAdjYvTpYeconjYe = Re((Yv*Yv.adjoint()*Ye.transpose()*Ye.
         conjugate()).trace());

   }

   if (loops > 2) {

   }

   return susy_traces;
}

std::ostream& operator<<(std::ostream& ostr, const UMSSM_susy_parameters& susy_pars)
{
   susy_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
