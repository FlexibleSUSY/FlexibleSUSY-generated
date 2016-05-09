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

// File generated at Mon 9 May 2016 11:57:13

#include "THDMIIMSSMBC_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"
#include "functors.hpp"

#include <iostream>

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces
#define TRACE_STRUCT_TYPE Soft_traces
#define CALCULATE_TRACES() calc_soft_traces(TRACE_STRUCT);

const int THDMIIMSSMBC_soft_parameters::numberOfParameters;

THDMIIMSSMBC_soft_parameters::THDMIIMSSMBC_soft_parameters(const THDMIIMSSMBC_input_parameters& input_)
   : THDMIIMSSMBC_susy_parameters(input_)
   , M122(0), M112(0), M222(0), v1(0), v2(0)

{
   set_number_of_parameters(numberOfParameters);
}

THDMIIMSSMBC_soft_parameters::THDMIIMSSMBC_soft_parameters(
   const THDMIIMSSMBC_susy_parameters& susy_model
   , double M122_, double M112_, double M222_, double v1_, double v2_

)
   : THDMIIMSSMBC_susy_parameters(susy_model)
   , M122(M122_), M112(M112_), M222(M222_), v1(v1_), v2(v2_)

{
   set_number_of_parameters(numberOfParameters);
}

Eigen::ArrayXd THDMIIMSSMBC_soft_parameters::beta() const
{
   return calc_beta().get().unaryExpr(Chop<double>(get_zero_threshold()));
}

THDMIIMSSMBC_soft_parameters THDMIIMSSMBC_soft_parameters::calc_beta() const
{
   double beta_M122 = 0.;
   double beta_M112 = 0.;
   double beta_M222 = 0.;
   double beta_v1 = 0.;
   double beta_v2 = 0.;

   if (get_loops() > 0) {
      TRACE_STRUCT_TYPE TRACE_STRUCT;
      CALCULATE_TRACES();

      beta_M122 += calc_beta_M122_one_loop(TRACE_STRUCT);
      beta_M112 += calc_beta_M112_one_loop(TRACE_STRUCT);
      beta_M222 += calc_beta_M222_one_loop(TRACE_STRUCT);
      beta_v1 += calc_beta_v1_one_loop(TRACE_STRUCT);
      beta_v2 += calc_beta_v2_one_loop(TRACE_STRUCT);

      if (get_loops() > 1) {
         beta_M122 += calc_beta_M122_two_loop(TRACE_STRUCT);
         beta_M112 += calc_beta_M112_two_loop(TRACE_STRUCT);
         beta_M222 += calc_beta_M222_two_loop(TRACE_STRUCT);
         beta_v1 += calc_beta_v1_two_loop(TRACE_STRUCT);
         beta_v2 += calc_beta_v2_two_loop(TRACE_STRUCT);

         if (get_loops() > 2) {

         }
      }
   }


   const THDMIIMSSMBC_susy_parameters susy_betas(THDMIIMSSMBC_susy_parameters::calc_beta());

   return THDMIIMSSMBC_soft_parameters(susy_betas, beta_M122, beta_M112, beta_M222, beta_v1, beta_v2);
}

void THDMIIMSSMBC_soft_parameters::clear()
{
   THDMIIMSSMBC_susy_parameters::clear();

   M122 = 0.;
   M112 = 0.;
   M222 = 0.;
   v1 = 0.;
   v2 = 0.;

}

Eigen::ArrayXd THDMIIMSSMBC_soft_parameters::get() const
{
   Eigen::ArrayXd pars(THDMIIMSSMBC_susy_parameters::get());
   pars.conservativeResize(numberOfParameters);

   pars(37) = M122;
   pars(38) = M112;
   pars(39) = M222;
   pars(40) = v1;
   pars(41) = v2;


   return pars;
}

void THDMIIMSSMBC_soft_parameters::print(std::ostream& ostr) const
{
   THDMIIMSSMBC_susy_parameters::print(ostr);
   ostr << "soft parameters at Q = " << get_scale() << ":\n";
   ostr << "M122 = " << M122 << '\n';
   ostr << "M112 = " << M112 << '\n';
   ostr << "M222 = " << M222 << '\n';
   ostr << "v1 = " << v1 << '\n';
   ostr << "v2 = " << v2 << '\n';

}

void THDMIIMSSMBC_soft_parameters::set(const Eigen::ArrayXd& pars)
{
   THDMIIMSSMBC_susy_parameters::set(pars);

   M122 = pars(37);
   M112 = pars(38);
   M222 = pars(39);
   v1 = pars(40);
   v2 = pars(41);

}

void THDMIIMSSMBC_soft_parameters::calc_soft_traces(Soft_traces& soft_traces) const
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

std::ostream& operator<<(std::ostream& ostr, const THDMIIMSSMBC_soft_parameters& soft_pars)
{
   soft_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
