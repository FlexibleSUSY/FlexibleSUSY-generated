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

// File generated at Mon 27 Feb 2017 13:21:54

#include "MRSSMtower_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"
#include "functors.hpp"

#include <iostream>

namespace flexiblesusy {

#define CLASSNAME MRSSMtower_susy_parameters
#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces
#define TRACE_STRUCT_TYPE Susy_traces
#define CALCULATE_TRACES() calc_susy_traces(TRACE_STRUCT);

const int MRSSMtower_susy_parameters::numberOfParameters;

MRSSMtower_susy_parameters::MRSSMtower_susy_parameters(const MRSSMtower_input_parameters& input_)
   : Beta_function()
   , Yd(Eigen::Matrix<double,3,3>::Zero()), Ye(Eigen::Matrix<double,3,3>::Zero(
   )), LamTD(0), LamTU(0), LamSD(0), LamSU(0), Yu(Eigen::Matrix<double,3,3>
   ::Zero()), Mu(0), MuD(0), MuU(0), g1(0), g2(0), g3(0), vd(0), vu(0), vT(0),
   vS(0)

   , input(input_)
{
   set_number_of_parameters(numberOfParameters);
}

MRSSMtower_susy_parameters::MRSSMtower_susy_parameters(
   double scale_, unsigned loops_, unsigned thresholds_,
   const MRSSMtower_input_parameters& input_
   , const Eigen::Matrix<double,3,3>& Yd_, const Eigen::Matrix<double,3,3>& Ye_
   , double LamTD_, double LamTU_, double LamSD_, double LamSU_, const
   Eigen::Matrix<double,3,3>& Yu_, double Mu_, double MuD_, double MuU_, double
   g1_, double g2_, double g3_, double vd_, double vu_, double vT_, double vS_

)
   : Beta_function()
   , Yd(Yd_), Ye(Ye_), LamTD(LamTD_), LamTU(LamTU_), LamSD(LamSD_), LamSU(
   LamSU_), Yu(Yu_), Mu(Mu_), MuD(MuD_), MuU(MuU_), g1(g1_), g2(g2_), g3(g3_),
   vd(vd_), vu(vu_), vT(vT_), vS(vS_)

   , input(input_)
{
   set_number_of_parameters(numberOfParameters);
   set_scale(scale_);
   set_loops(loops_);
   set_thresholds(thresholds_);
}

Eigen::ArrayXd MRSSMtower_susy_parameters::beta() const
{
   return calc_beta().get().unaryExpr(Chop<double>(get_zero_threshold()));
}

MRSSMtower_susy_parameters MRSSMtower_susy_parameters::calc_beta() const
{
   Eigen::Matrix<double,3,3> beta_Yd = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_Ye = Eigen::Matrix<double,3,3>::Zero();
   double beta_LamTD = 0.;
   double beta_LamTU = 0.;
   double beta_LamSD = 0.;
   double beta_LamSU = 0.;
   Eigen::Matrix<double,3,3> beta_Yu = Eigen::Matrix<double,3,3>::Zero();
   double beta_Mu = 0.;
   double beta_MuD = 0.;
   double beta_MuU = 0.;
   double beta_g1 = 0.;
   double beta_g2 = 0.;
   double beta_g3 = 0.;
   double beta_vd = 0.;
   double beta_vu = 0.;
   double beta_vT = 0.;
   double beta_vS = 0.;

   if (get_loops() > 0) {
      TRACE_STRUCT_TYPE TRACE_STRUCT;
      CALCULATE_TRACES();

      beta_Yd += calc_beta_Yd_one_loop(TRACE_STRUCT);
      beta_Ye += calc_beta_Ye_one_loop(TRACE_STRUCT);
      beta_LamTD += calc_beta_LamTD_one_loop(TRACE_STRUCT);
      beta_LamTU += calc_beta_LamTU_one_loop(TRACE_STRUCT);
      beta_LamSD += calc_beta_LamSD_one_loop(TRACE_STRUCT);
      beta_LamSU += calc_beta_LamSU_one_loop(TRACE_STRUCT);
      beta_Yu += calc_beta_Yu_one_loop(TRACE_STRUCT);
      beta_Mu += calc_beta_Mu_one_loop(TRACE_STRUCT);
      beta_MuD += calc_beta_MuD_one_loop(TRACE_STRUCT);
      beta_MuU += calc_beta_MuU_one_loop(TRACE_STRUCT);
      beta_g1 += calc_beta_g1_one_loop(TRACE_STRUCT);
      beta_g2 += calc_beta_g2_one_loop(TRACE_STRUCT);
      beta_g3 += calc_beta_g3_one_loop(TRACE_STRUCT);
      beta_vd += calc_beta_vd_one_loop(TRACE_STRUCT);
      beta_vu += calc_beta_vu_one_loop(TRACE_STRUCT);
      beta_vT += calc_beta_vT_one_loop(TRACE_STRUCT);
      beta_vS += calc_beta_vS_one_loop(TRACE_STRUCT);

      if (get_loops() > 1) {
         beta_Yd += calc_beta_Yd_two_loop(TRACE_STRUCT);
         beta_Ye += calc_beta_Ye_two_loop(TRACE_STRUCT);
         beta_LamTD += calc_beta_LamTD_two_loop(TRACE_STRUCT);
         beta_LamTU += calc_beta_LamTU_two_loop(TRACE_STRUCT);
         beta_LamSD += calc_beta_LamSD_two_loop(TRACE_STRUCT);
         beta_LamSU += calc_beta_LamSU_two_loop(TRACE_STRUCT);
         beta_Yu += calc_beta_Yu_two_loop(TRACE_STRUCT);
         beta_Mu += calc_beta_Mu_two_loop(TRACE_STRUCT);
         beta_MuD += calc_beta_MuD_two_loop(TRACE_STRUCT);
         beta_MuU += calc_beta_MuU_two_loop(TRACE_STRUCT);
         beta_g1 += calc_beta_g1_two_loop(TRACE_STRUCT);
         beta_g2 += calc_beta_g2_two_loop(TRACE_STRUCT);
         beta_g3 += calc_beta_g3_two_loop(TRACE_STRUCT);
         beta_vd += calc_beta_vd_two_loop(TRACE_STRUCT);
         beta_vu += calc_beta_vu_two_loop(TRACE_STRUCT);
         beta_vT += calc_beta_vT_two_loop(TRACE_STRUCT);
         beta_vS += calc_beta_vS_two_loop(TRACE_STRUCT);

         if (get_loops() > 2) {

         }
      }
   }


   return MRSSMtower_susy_parameters(get_scale(), get_loops(), get_thresholds(), input,
                    beta_Yd, beta_Ye, beta_LamTD, beta_LamTU, beta_LamSD, beta_LamSU, beta_Yu, beta_Mu, beta_MuD, beta_MuU, beta_g1, beta_g2, beta_g3, beta_vd, beta_vu, beta_vT, beta_vS);
}

MRSSMtower_susy_parameters MRSSMtower_susy_parameters::calc_beta(unsigned loops) const
{
   MRSSMtower_susy_parameters p(*this);
   p.set_loops(loops);

   return p.calc_beta();
}

void MRSSMtower_susy_parameters::clear()
{
   reset();
   Yd = Eigen::Matrix<double,3,3>::Zero();
   Ye = Eigen::Matrix<double,3,3>::Zero();
   LamTD = 0.;
   LamTU = 0.;
   LamSD = 0.;
   LamSU = 0.;
   Yu = Eigen::Matrix<double,3,3>::Zero();
   Mu = 0.;
   MuD = 0.;
   MuU = 0.;
   g1 = 0.;
   g2 = 0.;
   g3 = 0.;
   vd = 0.;
   vu = 0.;
   vT = 0.;
   vS = 0.;

}

Eigen::Matrix<double,3,3> CLASSNAME::get_SqSq() const
{
   Eigen::Matrix<double,3,3> anomDim;

   anomDim = (oneOver16PiSqr*(Yd.adjoint()*Yd + Yu.adjoint()*Yu -
      0.03333333333333333*(Sqr(g1) + 45*Sqr(g2) + 80*Sqr(g3))*UNITMATRIX(3)))
      .real();

   if (get_loops() > 1) {
      anomDim += (twoLoop*(-(AbsSqr(LamSU)*(Yu.adjoint()*Yu)) - 1.5*
         AbsSqr(LamTU)*(Yu.adjoint()*Yu) + 0.8*Sqr(g1)*(Yu.adjoint()*Yu) - 2*(
         Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 2*(Yu.adjoint()*Yu*Yu.adjoint()*Yu)
         + Yd.adjoint()*Yd*(-AbsSqr(LamSD) - 1.5*AbsSqr(LamTD) + 0.4*Sqr(g1) -
         3*(Yd*Yd.adjoint()).trace() - (Ye*Ye.adjoint()).trace()) - 3*(
         Yu.adjoint()*Yu)*(Yu*Yu.adjoint()).trace() + (0.2411111111111111*Power
         (g1,4) + 8.25*Power(g2,4) + 7.111111111111111*Power(g3,4) + 8*Sqr(g2)*
         Sqr(g3) + 0.011111111111111112*Sqr(g1)*(9*Sqr(g2) + 16*Sqr(g3)))*
         UNITMATRIX(3))).real();
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SlSl() const
{
   Eigen::Matrix<double,3,3> anomDim;

   anomDim = (oneOver16PiSqr*(Ye.adjoint()*Ye - 0.3*(Sqr(g1) + 5*Sqr(g2))
      *UNITMATRIX(3))).real();

   if (get_loops() > 1) {
      anomDim += (twoLoop*(-2*(Ye.adjoint()*Ye*Ye.adjoint()*Ye) +
         Ye.adjoint()*Ye*(-AbsSqr(LamSD) - 1.5*AbsSqr(LamTD) + 1.2*Sqr(g1) - 3*
         (Yd*Yd.adjoint()).trace() - (Ye*Ye.adjoint()).trace()) + 0.15*(15*
         Power(g1,4) + 55*Power(g2,4) + 6*Sqr(g1)*Sqr(g2))*UNITMATRIX(3))).real
         ();
   }

   return anomDim;
}

double CLASSNAME::get_SHdSHd() const
{
   double anomDim = 0;

   anomDim = Re(oneOver16PiSqr*(AbsSqr(LamSD) + 1.5*AbsSqr(LamTD) - 0.3*
      Sqr(g1) - 1.5*Sqr(g2) + 3*(Yd*Yd.adjoint()).trace() + (Ye*Ye.adjoint())
      .trace()));

   if (get_loops() > 1) {
      anomDim += Re(twoLoop*(2.25*Power(g1,4) + 8.25*Power(g2,4) -
         AbsSqr(LamSD)*(2*AbsSqr(LamSU) + 3*AbsSqr(LamTD)) + 0.9*Sqr(g1)*Sqr(g2
         ) + Conj(LamTD)*(-1.5*LamTD*AbsSqr(LamTU) + 6*LamTD*Sqr(g2)) - 3*Sqr(
         LamSD)*Sqr(Conj(LamSD)) - 3.75*Sqr(LamTD)*Sqr(Conj(LamTD)) - 0.4*Sqr(
         g1)*(Yd*Yd.adjoint()).trace() + 16*Sqr(g3)*(Yd*Yd.adjoint()).trace() +
         1.2*Sqr(g1)*(Ye*Ye.adjoint()).trace() - 9*(Yd*Yd.adjoint()*Yd*
         Yd.adjoint()).trace() - 3*(Yd*Yu.adjoint()*Yu*Yd.adjoint()).trace() -
         3*(Ye*Ye.adjoint()*Ye*Ye.adjoint()).trace()));
   }

   return anomDim;
}

double CLASSNAME::get_SHuSHu() const
{
   double anomDim = 0;

   anomDim = Re(oneOver16PiSqr*(AbsSqr(LamSU) - 0.3*(-5*AbsSqr(LamTU) +
      Sqr(g1) + 5*Sqr(g2) - 10*(Yu*Yu.adjoint()).trace())));

   if (get_loops() > 1) {
      anomDim += Re(twoLoop*(2.25*Power(g1,4) + 8.25*Power(g2,4) - 2*
         AbsSqr(LamSD)*AbsSqr(LamSU) - 3*AbsSqr(LamSU)*AbsSqr(LamTU) - 1.5*
         AbsSqr(LamTD)*AbsSqr(LamTU) + 6*AbsSqr(LamTU)*Sqr(g2) + 0.9*Sqr(g1)*
         Sqr(g2) - 3*Sqr(LamSU)*Sqr(Conj(LamSU)) - 3.75*Sqr(LamTU)*Sqr(Conj(
         LamTU)) + 0.8*Sqr(g1)*(Yu*Yu.adjoint()).trace() + 16*Sqr(g3)*(Yu*
         Yu.adjoint()).trace() - 3*(Yd*Yu.adjoint()*Yu*Yd.adjoint()).trace() -
         9*(Yu*Yu.adjoint()*Yu*Yu.adjoint()).trace()));
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SdRSdR() const
{
   Eigen::Matrix<double,3,3> anomDim;

   anomDim = (oneOver16PiSqr*(2*(Yd.conjugate()*Yd.transpose()) -
      0.13333333333333333*(Sqr(g1) + 20*Sqr(g3))*UNITMATRIX(3))).real();

   if (get_loops() > 1) {
      anomDim += (twoLoop*(-2*(Yd.conjugate()*Yd.transpose()*
         Yd.conjugate()*Yd.transpose() + Yd.conjugate()*Yu.transpose()*
         Yu.conjugate()*Yd.transpose()) + Yd.conjugate()*Yd.transpose()*(-2*
         AbsSqr(LamSD) - 3*AbsSqr(LamTD) + 0.4*Sqr(g1) + 6*Sqr(g2) - 6*(Yd*
         Yd.adjoint()).trace() - 2*(Ye*Ye.adjoint()).trace()) +
         0.08888888888888889*(11*Power(g1,4) + 80*Power(g3,4) + 8*Sqr(g1)*Sqr(
         g3))*UNITMATRIX(3))).real();
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SuRSuR() const
{
   Eigen::Matrix<double,3,3> anomDim;

   anomDim = (oneOver16PiSqr*(2*(Yu.conjugate()*Yu.transpose()) -
      0.5333333333333333*(Sqr(g1) + 5*Sqr(g3))*UNITMATRIX(3))).real();

   if (get_loops() > 1) {
      anomDim += (twoLoop*(-2*(Yu.conjugate()*Yd.transpose()*
         Yd.conjugate()*Yu.transpose() + Yu.conjugate()*Yu.transpose()*
         Yu.conjugate()*Yu.transpose()) + Yu.conjugate()*Yu.transpose()*(-2*
         AbsSqr(LamSU) - 3*AbsSqr(LamTU) - 0.4*Sqr(g1) + 6*Sqr(g2) - 6*(Yu*
         Yu.adjoint()).trace()) + 0.14222222222222222*(29*Power(g1,4) + 50*
         Power(g3,4) + 20*Sqr(g1)*Sqr(g3))*UNITMATRIX(3))).real();
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SeRSeR() const
{
   Eigen::Matrix<double,3,3> anomDim;

   anomDim = (oneOver16PiSqr*(2*(Ye.conjugate()*Ye.transpose()) - 1.2*Sqr
      (g1)*UNITMATRIX(3))).real();

   if (get_loops() > 1) {
      anomDim += (twoLoop*(-2*(Ye.conjugate()*Ye.transpose()*
         Ye.conjugate()*Ye.transpose()) + Ye.conjugate()*Ye.transpose()*(-2*
         AbsSqr(LamSD) - 3*AbsSqr(LamTD) - 1.2*Sqr(g1) + 6*Sqr(g2) - 6*(Yd*
         Yd.adjoint()).trace() - 2*(Ye*Ye.adjoint()).trace()) + 10.08*Power(g1,
         4)*UNITMATRIX(3))).real();
   }

   return anomDim;
}

double CLASSNAME::get_SsSs() const
{
   double anomDim = 0;

   anomDim = Re(2*oneOver16PiSqr*(AbsSqr(LamSD) + AbsSqr(LamSU)));

   if (get_loops() > 1) {
      anomDim += Re(twoLoop*(-4*Sqr(LamSD)*Sqr(Conj(LamSD)) + Conj(
         LamSD)*(-6*LamSD*AbsSqr(LamTD) + 1.2*LamSD*Sqr(g1) + 6*LamSD*Sqr(g2) -
         6*LamSD*(Yd*Yd.adjoint()).trace() - 2*LamSD*(Ye*Ye.adjoint()).trace()
         ) - 0.4*AbsSqr(LamSU)*(10*AbsSqr(LamSU) - 3*(-5*AbsSqr(LamTU) + Sqr(g1
         ) + 5*Sqr(g2) - 5*(Yu*Yu.adjoint()).trace()))));
   }

   return anomDim;
}

double CLASSNAME::get_STST() const
{
   double anomDim = 0;

   anomDim = Re(oneOver16PiSqr*(AbsSqr(LamTD) + AbsSqr(LamTU) - 4*Sqr(g2)
      ));

   if (get_loops() > 1) {
      anomDim += Re(twoLoop*(32*Power(g2,4) - 3*Sqr(LamTD)*Sqr(Conj(
         LamTD)) - 3*Sqr(LamTU)*Sqr(Conj(LamTU)) + 0.2*AbsSqr(LamTD)*(-10*
         AbsSqr(LamSD) + 3*Sqr(g1) - 5*Sqr(g2) - 15*(Yd*Yd.adjoint()).trace() -
         5*(Ye*Ye.adjoint()).trace()) + Conj(LamTU)*(-2*LamTU*AbsSqr(LamSU) +
         0.6*LamTU*Sqr(g1) - LamTU*Sqr(g2) - 3*LamTU*(Yu*Yu.adjoint()).trace())
         ));
   }

   return anomDim;
}

double CLASSNAME::get_SOcSOc() const
{
   double anomDim = 0;

   anomDim = Re(-6*oneOver16PiSqr*Sqr(g3));

   if (get_loops() > 1) {
      anomDim += Re(36*Power(g3,4)*twoLoop);
   }

   return anomDim;
}

double CLASSNAME::get_SRdSRd() const
{
   double anomDim = 0;

   anomDim = Re(oneOver16PiSqr*(AbsSqr(LamSD) - 0.3*(-5*AbsSqr(LamTD) +
      Sqr(g1) + 5*Sqr(g2))));

   if (get_loops() > 1) {
      anomDim += Re(twoLoop*(-3*Sqr(LamSD)*Sqr(Conj(LamSD)) + 0.15*(15
         *Power(g1,4) + 55*Power(g2,4) + 6*Sqr(g1)*Sqr(g2) - 25*Sqr(LamTD)*Sqr(
         Conj(LamTD)) + 10*AbsSqr(LamTD)*(-AbsSqr(LamTU) + 4*Sqr(g2) - 3*(Yd*
         Yd.adjoint()).trace() - (Ye*Ye.adjoint()).trace())) - AbsSqr(LamSD)*(2
         *AbsSqr(LamSU) + 3*AbsSqr(LamTD) + 3*(Yd*Yd.adjoint()).trace() + (Ye*
         Ye.adjoint()).trace())));
   }

   return anomDim;
}

double CLASSNAME::get_SRuSRu() const
{
   double anomDim = 0;

   anomDim = Re(oneOver16PiSqr*(AbsSqr(LamSU) - 0.3*(-5*AbsSqr(LamTU) +
      Sqr(g1) + 5*Sqr(g2))));

   if (get_loops() > 1) {
      anomDim += Re(twoLoop*(-2*AbsSqr(LamSD)*AbsSqr(LamSU) + 0.15*(15
         *Power(g1,4) + 55*Power(g2,4) + 6*Sqr(g1)*Sqr(g2) - 20*Sqr(LamSU)*Sqr(
         Conj(LamSU)) - 25*Sqr(LamTU)*Sqr(Conj(LamTU)) + 10*AbsSqr(LamTU)*(
         -AbsSqr(LamTD) + 4*Sqr(g2) - 3*(Yu*Yu.adjoint()).trace()) - 20*AbsSqr(
         LamSU)*(AbsSqr(LamTU) + (Yu*Yu.adjoint()).trace()))));
   }

   return anomDim;
}


Eigen::ArrayXd MRSSMtower_susy_parameters::get() const
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
   pars(18) = LamTD;
   pars(19) = LamTU;
   pars(20) = LamSD;
   pars(21) = LamSU;
   pars(22) = Yu(0,0);
   pars(23) = Yu(0,1);
   pars(24) = Yu(0,2);
   pars(25) = Yu(1,0);
   pars(26) = Yu(1,1);
   pars(27) = Yu(1,2);
   pars(28) = Yu(2,0);
   pars(29) = Yu(2,1);
   pars(30) = Yu(2,2);
   pars(31) = Mu;
   pars(32) = MuD;
   pars(33) = MuU;
   pars(34) = g1;
   pars(35) = g2;
   pars(36) = g3;
   pars(37) = vd;
   pars(38) = vu;
   pars(39) = vT;
   pars(40) = vS;


   return pars;
}

void MRSSMtower_susy_parameters::print(std::ostream& ostr) const
{
   ostr << "susy parameters at Q = " << get_scale() << ":\n";
   ostr << "Yd = " << Yd << '\n';
   ostr << "Ye = " << Ye << '\n';
   ostr << "LamTD = " << LamTD << '\n';
   ostr << "LamTU = " << LamTU << '\n';
   ostr << "LamSD = " << LamSD << '\n';
   ostr << "LamSU = " << LamSU << '\n';
   ostr << "Yu = " << Yu << '\n';
   ostr << "Mu = " << Mu << '\n';
   ostr << "MuD = " << MuD << '\n';
   ostr << "MuU = " << MuU << '\n';
   ostr << "g1 = " << g1 << '\n';
   ostr << "g2 = " << g2 << '\n';
   ostr << "g3 = " << g3 << '\n';
   ostr << "vd = " << vd << '\n';
   ostr << "vu = " << vu << '\n';
   ostr << "vT = " << vT << '\n';
   ostr << "vS = " << vS << '\n';

}

void MRSSMtower_susy_parameters::set(const Eigen::ArrayXd& pars)
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
   LamTD = pars(18);
   LamTU = pars(19);
   LamSD = pars(20);
   LamSU = pars(21);
   Yu(0,0) = pars(22);
   Yu(0,1) = pars(23);
   Yu(0,2) = pars(24);
   Yu(1,0) = pars(25);
   Yu(1,1) = pars(26);
   Yu(1,2) = pars(27);
   Yu(2,0) = pars(28);
   Yu(2,1) = pars(29);
   Yu(2,2) = pars(30);
   Mu = pars(31);
   MuD = pars(32);
   MuU = pars(33);
   g1 = pars(34);
   g2 = pars(35);
   g3 = pars(36);
   vd = pars(37);
   vu = pars(38);
   vT = pars(39);
   vS = pars(40);

}

const MRSSMtower_input_parameters& MRSSMtower_susy_parameters::get_input() const
{
   return input;
}

MRSSMtower_input_parameters& MRSSMtower_susy_parameters::get_input()
{
   return input;
}

void MRSSMtower_susy_parameters::set_input_parameters(const MRSSMtower_input_parameters& input_)
{
   input = input_;
}

void MRSSMtower_susy_parameters::calc_susy_traces(Susy_traces& susy_traces) const
{
   if (get_loops() > 0) {
      TRACE_STRUCT.traceYdAdjYd = Re((Yd*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYeAdjYe = Re((Ye*Ye.adjoint()).trace());
      TRACE_STRUCT.traceYuAdjYu = Re((Yu*Yu.adjoint()).trace());

   }

   if (get_loops() > 1) {
      TRACE_STRUCT.traceYdAdjYdYdAdjYd = Re((Yd*Yd.adjoint()*Yd*Yd.adjoint())
         .trace());
      TRACE_STRUCT.traceYdAdjYuYuAdjYd = Re((Yd*Yu.adjoint()*Yu*Yd.adjoint())
         .trace());
      TRACE_STRUCT.traceYeAdjYeYeAdjYe = Re((Ye*Ye.adjoint()*Ye*Ye.adjoint())
         .trace());
      TRACE_STRUCT.traceYuAdjYuYuAdjYu = Re((Yu*Yu.adjoint()*Yu*Yu.adjoint())
         .trace());

   }

   if (get_loops() > 2) {

   }
}

std::ostream& operator<<(std::ostream& ostr, const MRSSMtower_susy_parameters& susy_pars)
{
   susy_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
