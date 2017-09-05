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

// File generated at Tue 5 Sep 2017 10:32:14

#include "SplitMSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"
#include "functors.hpp"

#include <iostream>

namespace flexiblesusy {

#define CLASSNAME SplitMSSM_susy_parameters
#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces
#define TRACE_STRUCT_TYPE Susy_traces
#define CALCULATE_TRACES() calc_susy_traces(TRACE_STRUCT);

const int SplitMSSM_susy_parameters::numberOfParameters;

SplitMSSM_susy_parameters::SplitMSSM_susy_parameters(const SplitMSSM_input_parameters& input_)
   : Beta_function()
   , g1(0), g2(0), g3(0), Lambdax(0), Yu(Eigen::Matrix<double,3,3>::Zero()), Yd
   (Eigen::Matrix<double,3,3>::Zero()), Ye(Eigen::Matrix<double,3,3>::Zero()),
   gYd(0), g2d(0), gYu(0), g2u(0)

   , input(input_)
{
   set_number_of_parameters(numberOfParameters);
}

SplitMSSM_susy_parameters::SplitMSSM_susy_parameters(
   double scale_, unsigned loops_, unsigned thresholds_,
   const SplitMSSM_input_parameters& input_
   , double g1_, double g2_, double g3_, double Lambdax_, const Eigen::Matrix<
   double,3,3>& Yu_, const Eigen::Matrix<double,3,3>& Yd_, const Eigen::Matrix<
   double,3,3>& Ye_, double gYd_, double g2d_, double gYu_, double g2u_

)
   : Beta_function()
   , g1(g1_), g2(g2_), g3(g3_), Lambdax(Lambdax_), Yu(Yu_), Yd(Yd_), Ye(Ye_),
   gYd(gYd_), g2d(g2d_), gYu(gYu_), g2u(g2u_)

   , input(input_)
{
   set_number_of_parameters(numberOfParameters);
   set_scale(scale_);
   set_loops(loops_);
   set_thresholds(thresholds_);
}

Eigen::ArrayXd SplitMSSM_susy_parameters::beta() const
{
   return calc_beta().get().unaryExpr(Chop<double>(get_zero_threshold()));
}

SplitMSSM_susy_parameters SplitMSSM_susy_parameters::calc_beta() const
{
   double beta_g1 = 0.;
   double beta_g2 = 0.;
   double beta_g3 = 0.;
   double beta_Lambdax = 0.;
   Eigen::Matrix<double,3,3> beta_Yu = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_Yd = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_Ye = Eigen::Matrix<double,3,3>::Zero();
   double beta_gYd = 0.;
   double beta_g2d = 0.;
   double beta_gYu = 0.;
   double beta_g2u = 0.;

   if (get_loops() > 0) {
      TRACE_STRUCT_TYPE TRACE_STRUCT;
      CALCULATE_TRACES();

      beta_g1 += calc_beta_g1_one_loop(TRACE_STRUCT);
      beta_g2 += calc_beta_g2_one_loop(TRACE_STRUCT);
      beta_g3 += calc_beta_g3_one_loop(TRACE_STRUCT);
      beta_Lambdax += calc_beta_Lambdax_one_loop(TRACE_STRUCT);
      beta_Yu += calc_beta_Yu_one_loop(TRACE_STRUCT);
      beta_Yd += calc_beta_Yd_one_loop(TRACE_STRUCT);
      beta_Ye += calc_beta_Ye_one_loop(TRACE_STRUCT);
      beta_gYd += calc_beta_gYd_one_loop(TRACE_STRUCT);
      beta_g2d += calc_beta_g2d_one_loop(TRACE_STRUCT);
      beta_gYu += calc_beta_gYu_one_loop(TRACE_STRUCT);
      beta_g2u += calc_beta_g2u_one_loop(TRACE_STRUCT);

      if (get_loops() > 1) {
         beta_g1 += calc_beta_g1_two_loop(TRACE_STRUCT);
         beta_g2 += calc_beta_g2_two_loop(TRACE_STRUCT);
         beta_g3 += calc_beta_g3_two_loop(TRACE_STRUCT);
         beta_Lambdax += calc_beta_Lambdax_two_loop(TRACE_STRUCT);
         beta_Yu += calc_beta_Yu_two_loop(TRACE_STRUCT);
         beta_Yd += calc_beta_Yd_two_loop(TRACE_STRUCT);
         beta_Ye += calc_beta_Ye_two_loop(TRACE_STRUCT);
         beta_gYd += calc_beta_gYd_two_loop(TRACE_STRUCT);
         beta_g2d += calc_beta_g2d_two_loop(TRACE_STRUCT);
         beta_gYu += calc_beta_gYu_two_loop(TRACE_STRUCT);
         beta_g2u += calc_beta_g2u_two_loop(TRACE_STRUCT);

         if (get_loops() > 2) {

         }
      }
   }


   return SplitMSSM_susy_parameters(get_scale(), get_loops(), get_thresholds(), input,
                    beta_g1, beta_g2, beta_g3, beta_Lambdax, beta_Yu, beta_Yd, beta_Ye, beta_gYd, beta_g2d, beta_gYu, beta_g2u);
}

SplitMSSM_susy_parameters SplitMSSM_susy_parameters::calc_beta(unsigned loops) const
{
   SplitMSSM_susy_parameters p(*this);
   p.set_loops(loops);

   return p.calc_beta();
}

void SplitMSSM_susy_parameters::clear()
{
   reset();
   g1 = 0.;
   g2 = 0.;
   g3 = 0.;
   Lambdax = 0.;
   Yu = Eigen::Matrix<double,3,3>::Zero();
   Yd = Eigen::Matrix<double,3,3>::Zero();
   Ye = Eigen::Matrix<double,3,3>::Zero();
   gYd = 0.;
   g2d = 0.;
   gYu = 0.;
   g2u = 0.;

}



Eigen::ArrayXd SplitMSSM_susy_parameters::get() const
{
   Eigen::ArrayXd pars(numberOfParameters);

   pars(0) = g1;
   pars(1) = g2;
   pars(2) = g3;
   pars(3) = Lambdax;
   pars(4) = Yu(0,0);
   pars(5) = Yu(0,1);
   pars(6) = Yu(0,2);
   pars(7) = Yu(1,0);
   pars(8) = Yu(1,1);
   pars(9) = Yu(1,2);
   pars(10) = Yu(2,0);
   pars(11) = Yu(2,1);
   pars(12) = Yu(2,2);
   pars(13) = Yd(0,0);
   pars(14) = Yd(0,1);
   pars(15) = Yd(0,2);
   pars(16) = Yd(1,0);
   pars(17) = Yd(1,1);
   pars(18) = Yd(1,2);
   pars(19) = Yd(2,0);
   pars(20) = Yd(2,1);
   pars(21) = Yd(2,2);
   pars(22) = Ye(0,0);
   pars(23) = Ye(0,1);
   pars(24) = Ye(0,2);
   pars(25) = Ye(1,0);
   pars(26) = Ye(1,1);
   pars(27) = Ye(1,2);
   pars(28) = Ye(2,0);
   pars(29) = Ye(2,1);
   pars(30) = Ye(2,2);
   pars(31) = gYd;
   pars(32) = g2d;
   pars(33) = gYu;
   pars(34) = g2u;


   return pars;
}

void SplitMSSM_susy_parameters::print(std::ostream& ostr) const
{
   ostr << "susy parameters at Q = " << get_scale() << ":\n";
   ostr << "g1 = " << g1 << '\n';
   ostr << "g2 = " << g2 << '\n';
   ostr << "g3 = " << g3 << '\n';
   ostr << "Lambdax = " << Lambdax << '\n';
   ostr << "Yu = " << Yu << '\n';
   ostr << "Yd = " << Yd << '\n';
   ostr << "Ye = " << Ye << '\n';
   ostr << "gYd = " << gYd << '\n';
   ostr << "g2d = " << g2d << '\n';
   ostr << "gYu = " << gYu << '\n';
   ostr << "g2u = " << g2u << '\n';

}

void SplitMSSM_susy_parameters::set(const Eigen::ArrayXd& pars)
{
   g1 = pars(0);
   g2 = pars(1);
   g3 = pars(2);
   Lambdax = pars(3);
   Yu(0,0) = pars(4);
   Yu(0,1) = pars(5);
   Yu(0,2) = pars(6);
   Yu(1,0) = pars(7);
   Yu(1,1) = pars(8);
   Yu(1,2) = pars(9);
   Yu(2,0) = pars(10);
   Yu(2,1) = pars(11);
   Yu(2,2) = pars(12);
   Yd(0,0) = pars(13);
   Yd(0,1) = pars(14);
   Yd(0,2) = pars(15);
   Yd(1,0) = pars(16);
   Yd(1,1) = pars(17);
   Yd(1,2) = pars(18);
   Yd(2,0) = pars(19);
   Yd(2,1) = pars(20);
   Yd(2,2) = pars(21);
   Ye(0,0) = pars(22);
   Ye(0,1) = pars(23);
   Ye(0,2) = pars(24);
   Ye(1,0) = pars(25);
   Ye(1,1) = pars(26);
   Ye(1,2) = pars(27);
   Ye(2,0) = pars(28);
   Ye(2,1) = pars(29);
   Ye(2,2) = pars(30);
   gYd = pars(31);
   g2d = pars(32);
   gYu = pars(33);
   g2u = pars(34);

}

const SplitMSSM_input_parameters& SplitMSSM_susy_parameters::get_input() const
{
   return input;
}

SplitMSSM_input_parameters& SplitMSSM_susy_parameters::get_input()
{
   return input;
}

void SplitMSSM_susy_parameters::set_input_parameters(const SplitMSSM_input_parameters& input_)
{
   input = input_;
}

void SplitMSSM_susy_parameters::calc_susy_traces(Susy_traces& susy_traces) const
{
   if (get_loops() > 0) {
      TRACE_STRUCT.traceYdAdjYd = Re((Yd*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYeAdjYe = Re((Ye*Ye.adjoint()).trace());
      TRACE_STRUCT.traceYuAdjYu = Re((Yu*Yu.adjoint()).trace());
      TRACE_STRUCT.traceYdAdjYdYdAdjYd = Re((Yd*Yd.adjoint()*Yd*Yd.adjoint())
         .trace());
      TRACE_STRUCT.traceYeAdjYeYeAdjYe = Re((Ye*Ye.adjoint()*Ye*Ye.adjoint())
         .trace());
      TRACE_STRUCT.traceYuAdjYuYuAdjYu = Re((Yu*Yu.adjoint()*Yu*Yu.adjoint())
         .trace());

   }

   if (get_loops() > 1) {
      TRACE_STRUCT.traceYdAdjYuYuAdjYd = Re((Yd*Yu.adjoint()*Yu*Yd.adjoint())
         .trace());
      TRACE_STRUCT.traceYdAdjYdYdAdjYdYdAdjYd = Re((Yd*Yd.adjoint()*Yd*Yd.adjoint(
         )*Yd*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYdAdjYdYdAdjYuYuAdjYd = Re((Yd*Yd.adjoint()*Yd*Yu.adjoint(
         )*Yu*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYdAdjYuYuAdjYdYdAdjYd = Re((Yd*Yu.adjoint()*Yu*Yd.adjoint(
         )*Yd*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYdAdjYuYuAdjYuYuAdjYd = Re((Yd*Yu.adjoint()*Yu*Yu.adjoint(
         )*Yu*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYeAdjYeYeAdjYeYeAdjYe = Re((Ye*Ye.adjoint()*Ye*Ye.adjoint(
         )*Ye*Ye.adjoint()).trace());
      TRACE_STRUCT.traceYuAdjYuYuAdjYuYuAdjYu = Re((Yu*Yu.adjoint()*Yu*Yu.adjoint(
         )*Yu*Yu.adjoint()).trace());

   }

   if (get_loops() > 2) {

   }
}

std::ostream& operator<<(std::ostream& ostr, const SplitMSSM_susy_parameters& susy_pars)
{
   susy_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
