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

// File generated at Wed 12 Apr 2017 10:59:41

#include "HGTHDMIIMSSMBC_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"
#include "functors.hpp"

#include <iostream>

namespace flexiblesusy {

#define CLASSNAME HGTHDMIIMSSMBC_susy_parameters
#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces
#define TRACE_STRUCT_TYPE Susy_traces
#define CALCULATE_TRACES() calc_susy_traces(TRACE_STRUCT);

const int HGTHDMIIMSSMBC_susy_parameters::numberOfParameters;

HGTHDMIIMSSMBC_susy_parameters::HGTHDMIIMSSMBC_susy_parameters(const HGTHDMIIMSSMBC_input_parameters& input_)
   : Beta_function()
   , g1(0), g2(0), g3(0), Lambda6(0), Lambda5(0), Lambda7(0), Lambda1(0),
   Lambda4(0), Lambda3(0), Lambda2(0), Yu(Eigen::Matrix<double,3,3>::Zero()),
   Yd(Eigen::Matrix<double,3,3>::Zero()), Ye(Eigen::Matrix<double,3,3>::Zero())
   , g1dp(0), g1d(0), g2up(0), g2u(0)

   , input(input_)
{
   set_number_of_parameters(numberOfParameters);
}

HGTHDMIIMSSMBC_susy_parameters::HGTHDMIIMSSMBC_susy_parameters(
   double scale_, unsigned loops_, unsigned thresholds_,
   const HGTHDMIIMSSMBC_input_parameters& input_
   , double g1_, double g2_, double g3_, double Lambda6_, double Lambda5_,
   double Lambda7_, double Lambda1_, double Lambda4_, double Lambda3_, double
   Lambda2_, const Eigen::Matrix<double,3,3>& Yu_, const Eigen::Matrix<double,3
   ,3>& Yd_, const Eigen::Matrix<double,3,3>& Ye_, double g1dp_, double g1d_,
   double g2up_, double g2u_

)
   : Beta_function()
   , g1(g1_), g2(g2_), g3(g3_), Lambda6(Lambda6_), Lambda5(Lambda5_), Lambda7(
   Lambda7_), Lambda1(Lambda1_), Lambda4(Lambda4_), Lambda3(Lambda3_), Lambda2(
   Lambda2_), Yu(Yu_), Yd(Yd_), Ye(Ye_), g1dp(g1dp_), g1d(g1d_), g2up(g2up_),
   g2u(g2u_)

   , input(input_)
{
   set_number_of_parameters(numberOfParameters);
   set_scale(scale_);
   set_loops(loops_);
   set_thresholds(thresholds_);
}

Eigen::ArrayXd HGTHDMIIMSSMBC_susy_parameters::beta() const
{
   return calc_beta().get().unaryExpr(Chop<double>(get_zero_threshold()));
}

HGTHDMIIMSSMBC_susy_parameters HGTHDMIIMSSMBC_susy_parameters::calc_beta() const
{
   double beta_g1 = 0.;
   double beta_g2 = 0.;
   double beta_g3 = 0.;
   double beta_Lambda6 = 0.;
   double beta_Lambda5 = 0.;
   double beta_Lambda7 = 0.;
   double beta_Lambda1 = 0.;
   double beta_Lambda4 = 0.;
   double beta_Lambda3 = 0.;
   double beta_Lambda2 = 0.;
   Eigen::Matrix<double,3,3> beta_Yu = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_Yd = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_Ye = Eigen::Matrix<double,3,3>::Zero();
   double beta_g1dp = 0.;
   double beta_g1d = 0.;
   double beta_g2up = 0.;
   double beta_g2u = 0.;

   if (get_loops() > 0) {
      TRACE_STRUCT_TYPE TRACE_STRUCT;
      CALCULATE_TRACES();

      beta_g1 += calc_beta_g1_one_loop(TRACE_STRUCT);
      beta_g2 += calc_beta_g2_one_loop(TRACE_STRUCT);
      beta_g3 += calc_beta_g3_one_loop(TRACE_STRUCT);
      beta_Lambda6 += calc_beta_Lambda6_one_loop(TRACE_STRUCT);
      beta_Lambda5 += calc_beta_Lambda5_one_loop(TRACE_STRUCT);
      beta_Lambda7 += calc_beta_Lambda7_one_loop(TRACE_STRUCT);
      beta_Lambda1 += calc_beta_Lambda1_one_loop(TRACE_STRUCT);
      beta_Lambda4 += calc_beta_Lambda4_one_loop(TRACE_STRUCT);
      beta_Lambda3 += calc_beta_Lambda3_one_loop(TRACE_STRUCT);
      beta_Lambda2 += calc_beta_Lambda2_one_loop(TRACE_STRUCT);
      beta_Yu += calc_beta_Yu_one_loop(TRACE_STRUCT);
      beta_Yd += calc_beta_Yd_one_loop(TRACE_STRUCT);
      beta_Ye += calc_beta_Ye_one_loop(TRACE_STRUCT);
      beta_g1dp += calc_beta_g1dp_one_loop(TRACE_STRUCT);
      beta_g1d += calc_beta_g1d_one_loop(TRACE_STRUCT);
      beta_g2up += calc_beta_g2up_one_loop(TRACE_STRUCT);
      beta_g2u += calc_beta_g2u_one_loop(TRACE_STRUCT);

      if (get_loops() > 1) {
         beta_g1 += calc_beta_g1_two_loop(TRACE_STRUCT);
         beta_g2 += calc_beta_g2_two_loop(TRACE_STRUCT);
         beta_g3 += calc_beta_g3_two_loop(TRACE_STRUCT);
         beta_Lambda6 += calc_beta_Lambda6_two_loop(TRACE_STRUCT);
         beta_Lambda5 += calc_beta_Lambda5_two_loop(TRACE_STRUCT);
         beta_Lambda7 += calc_beta_Lambda7_two_loop(TRACE_STRUCT);
         beta_Lambda1 += calc_beta_Lambda1_two_loop(TRACE_STRUCT);
         beta_Lambda4 += calc_beta_Lambda4_two_loop(TRACE_STRUCT);
         beta_Lambda3 += calc_beta_Lambda3_two_loop(TRACE_STRUCT);
         beta_Lambda2 += calc_beta_Lambda2_two_loop(TRACE_STRUCT);
         beta_Yu += calc_beta_Yu_two_loop(TRACE_STRUCT);
         beta_Yd += calc_beta_Yd_two_loop(TRACE_STRUCT);
         beta_Ye += calc_beta_Ye_two_loop(TRACE_STRUCT);
         beta_g1dp += calc_beta_g1dp_two_loop(TRACE_STRUCT);
         beta_g1d += calc_beta_g1d_two_loop(TRACE_STRUCT);
         beta_g2up += calc_beta_g2up_two_loop(TRACE_STRUCT);
         beta_g2u += calc_beta_g2u_two_loop(TRACE_STRUCT);

         if (get_loops() > 2) {

         }
      }
   }


   return HGTHDMIIMSSMBC_susy_parameters(get_scale(), get_loops(), get_thresholds(), input,
                    beta_g1, beta_g2, beta_g3, beta_Lambda6, beta_Lambda5, beta_Lambda7, beta_Lambda1, beta_Lambda4, beta_Lambda3, beta_Lambda2, beta_Yu, beta_Yd, beta_Ye, beta_g1dp, beta_g1d, beta_g2up, beta_g2u);
}

HGTHDMIIMSSMBC_susy_parameters HGTHDMIIMSSMBC_susy_parameters::calc_beta(unsigned loops) const
{
   HGTHDMIIMSSMBC_susy_parameters p(*this);
   p.set_loops(loops);

   return p.calc_beta();
}

void HGTHDMIIMSSMBC_susy_parameters::clear()
{
   reset();
   g1 = 0.;
   g2 = 0.;
   g3 = 0.;
   Lambda6 = 0.;
   Lambda5 = 0.;
   Lambda7 = 0.;
   Lambda1 = 0.;
   Lambda4 = 0.;
   Lambda3 = 0.;
   Lambda2 = 0.;
   Yu = Eigen::Matrix<double,3,3>::Zero();
   Yd = Eigen::Matrix<double,3,3>::Zero();
   Ye = Eigen::Matrix<double,3,3>::Zero();
   g1dp = 0.;
   g1d = 0.;
   g2up = 0.;
   g2u = 0.;

}



Eigen::ArrayXd HGTHDMIIMSSMBC_susy_parameters::get() const
{
   Eigen::ArrayXd pars(numberOfParameters);

   pars(0) = g1;
   pars(1) = g2;
   pars(2) = g3;
   pars(3) = Lambda6;
   pars(4) = Lambda5;
   pars(5) = Lambda7;
   pars(6) = Lambda1;
   pars(7) = Lambda4;
   pars(8) = Lambda3;
   pars(9) = Lambda2;
   pars(10) = Yu(0,0);
   pars(11) = Yu(0,1);
   pars(12) = Yu(0,2);
   pars(13) = Yu(1,0);
   pars(14) = Yu(1,1);
   pars(15) = Yu(1,2);
   pars(16) = Yu(2,0);
   pars(17) = Yu(2,1);
   pars(18) = Yu(2,2);
   pars(19) = Yd(0,0);
   pars(20) = Yd(0,1);
   pars(21) = Yd(0,2);
   pars(22) = Yd(1,0);
   pars(23) = Yd(1,1);
   pars(24) = Yd(1,2);
   pars(25) = Yd(2,0);
   pars(26) = Yd(2,1);
   pars(27) = Yd(2,2);
   pars(28) = Ye(0,0);
   pars(29) = Ye(0,1);
   pars(30) = Ye(0,2);
   pars(31) = Ye(1,0);
   pars(32) = Ye(1,1);
   pars(33) = Ye(1,2);
   pars(34) = Ye(2,0);
   pars(35) = Ye(2,1);
   pars(36) = Ye(2,2);
   pars(37) = g1dp;
   pars(38) = g1d;
   pars(39) = g2up;
   pars(40) = g2u;


   return pars;
}

void HGTHDMIIMSSMBC_susy_parameters::print(std::ostream& ostr) const
{
   ostr << "susy parameters at Q = " << get_scale() << ":\n";
   ostr << "g1 = " << g1 << '\n';
   ostr << "g2 = " << g2 << '\n';
   ostr << "g3 = " << g3 << '\n';
   ostr << "Lambda6 = " << Lambda6 << '\n';
   ostr << "Lambda5 = " << Lambda5 << '\n';
   ostr << "Lambda7 = " << Lambda7 << '\n';
   ostr << "Lambda1 = " << Lambda1 << '\n';
   ostr << "Lambda4 = " << Lambda4 << '\n';
   ostr << "Lambda3 = " << Lambda3 << '\n';
   ostr << "Lambda2 = " << Lambda2 << '\n';
   ostr << "Yu = " << Yu << '\n';
   ostr << "Yd = " << Yd << '\n';
   ostr << "Ye = " << Ye << '\n';
   ostr << "g1dp = " << g1dp << '\n';
   ostr << "g1d = " << g1d << '\n';
   ostr << "g2up = " << g2up << '\n';
   ostr << "g2u = " << g2u << '\n';

}

void HGTHDMIIMSSMBC_susy_parameters::set(const Eigen::ArrayXd& pars)
{
   g1 = pars(0);
   g2 = pars(1);
   g3 = pars(2);
   Lambda6 = pars(3);
   Lambda5 = pars(4);
   Lambda7 = pars(5);
   Lambda1 = pars(6);
   Lambda4 = pars(7);
   Lambda3 = pars(8);
   Lambda2 = pars(9);
   Yu(0,0) = pars(10);
   Yu(0,1) = pars(11);
   Yu(0,2) = pars(12);
   Yu(1,0) = pars(13);
   Yu(1,1) = pars(14);
   Yu(1,2) = pars(15);
   Yu(2,0) = pars(16);
   Yu(2,1) = pars(17);
   Yu(2,2) = pars(18);
   Yd(0,0) = pars(19);
   Yd(0,1) = pars(20);
   Yd(0,2) = pars(21);
   Yd(1,0) = pars(22);
   Yd(1,1) = pars(23);
   Yd(1,2) = pars(24);
   Yd(2,0) = pars(25);
   Yd(2,1) = pars(26);
   Yd(2,2) = pars(27);
   Ye(0,0) = pars(28);
   Ye(0,1) = pars(29);
   Ye(0,2) = pars(30);
   Ye(1,0) = pars(31);
   Ye(1,1) = pars(32);
   Ye(1,2) = pars(33);
   Ye(2,0) = pars(34);
   Ye(2,1) = pars(35);
   Ye(2,2) = pars(36);
   g1dp = pars(37);
   g1d = pars(38);
   g2up = pars(39);
   g2u = pars(40);

}

const HGTHDMIIMSSMBC_input_parameters& HGTHDMIIMSSMBC_susy_parameters::get_input() const
{
   return input;
}

HGTHDMIIMSSMBC_input_parameters& HGTHDMIIMSSMBC_susy_parameters::get_input()
{
   return input;
}

void HGTHDMIIMSSMBC_susy_parameters::set_input_parameters(const HGTHDMIIMSSMBC_input_parameters& input_)
{
   input = input_;
}

void HGTHDMIIMSSMBC_susy_parameters::calc_susy_traces(Susy_traces& susy_traces) const
{
   if (get_loops() > 0) {
      TRACE_STRUCT.traceYdAdjYd = Re((Yd*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYeAdjYe = Re((Ye*Ye.adjoint()).trace());
      TRACE_STRUCT.traceYuAdjYu = Re((Yu*Yu.adjoint()).trace());
      TRACE_STRUCT.traceYdAdjYdYdAdjYd = Re((Yd*Yd.adjoint()*Yd*Yd.adjoint())
         .trace());
      TRACE_STRUCT.traceYdAdjYuYuAdjYd = Re((Yd*Yu.adjoint()*Yu*Yd.adjoint())
         .trace());
      TRACE_STRUCT.traceYeAdjYeYeAdjYe = Re((Ye*Ye.adjoint()*Ye*Ye.adjoint())
         .trace());
      TRACE_STRUCT.traceYuAdjYuYuAdjYu = Re((Yu*Yu.adjoint()*Yu*Yu.adjoint())
         .trace());

   }

   if (get_loops() > 1) {
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

std::ostream& operator<<(std::ostream& ostr, const HGTHDMIIMSSMBC_susy_parameters& susy_pars)
{
   susy_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
