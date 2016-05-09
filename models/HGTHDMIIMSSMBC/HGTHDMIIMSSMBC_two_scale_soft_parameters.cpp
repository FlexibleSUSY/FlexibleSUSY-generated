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

// File generated at Mon 9 May 2016 12:04:21

#include "HGTHDMIIMSSMBC_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"
#include "functors.hpp"

#include <iostream>

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces
#define TRACE_STRUCT_TYPE Soft_traces
#define CALCULATE_TRACES() calc_soft_traces(TRACE_STRUCT);

const int HGTHDMIIMSSMBC_soft_parameters::numberOfParameters;

HGTHDMIIMSSMBC_soft_parameters::HGTHDMIIMSSMBC_soft_parameters(const HGTHDMIIMSSMBC_input_parameters& input_)
   : HGTHDMIIMSSMBC_susy_parameters(input_)
   , MassB(0), MassG(0), MassWB(0), Mu(0), M122(0), M112(0), M222(0), v1(0), v2
   (0)

{
   set_number_of_parameters(numberOfParameters);
}

HGTHDMIIMSSMBC_soft_parameters::HGTHDMIIMSSMBC_soft_parameters(
   const HGTHDMIIMSSMBC_susy_parameters& susy_model
   , double MassB_, double MassG_, double MassWB_, double Mu_, double M122_,
   double M112_, double M222_, double v1_, double v2_

)
   : HGTHDMIIMSSMBC_susy_parameters(susy_model)
   , MassB(MassB_), MassG(MassG_), MassWB(MassWB_), Mu(Mu_), M122(M122_), M112(
   M112_), M222(M222_), v1(v1_), v2(v2_)

{
   set_number_of_parameters(numberOfParameters);
}

Eigen::ArrayXd HGTHDMIIMSSMBC_soft_parameters::beta() const
{
   return calc_beta().get().unaryExpr(Chop<double>(get_zero_threshold()));
}

HGTHDMIIMSSMBC_soft_parameters HGTHDMIIMSSMBC_soft_parameters::calc_beta() const
{
   double beta_MassB = 0.;
   double beta_MassG = 0.;
   double beta_MassWB = 0.;
   double beta_Mu = 0.;
   double beta_M122 = 0.;
   double beta_M112 = 0.;
   double beta_M222 = 0.;
   double beta_v1 = 0.;
   double beta_v2 = 0.;

   if (get_loops() > 0) {
      TRACE_STRUCT_TYPE TRACE_STRUCT;
      CALCULATE_TRACES();

      beta_MassB += calc_beta_MassB_one_loop(TRACE_STRUCT);
      beta_MassG += calc_beta_MassG_one_loop(TRACE_STRUCT);
      beta_MassWB += calc_beta_MassWB_one_loop(TRACE_STRUCT);
      beta_Mu += calc_beta_Mu_one_loop(TRACE_STRUCT);
      beta_M122 += calc_beta_M122_one_loop(TRACE_STRUCT);
      beta_M112 += calc_beta_M112_one_loop(TRACE_STRUCT);
      beta_M222 += calc_beta_M222_one_loop(TRACE_STRUCT);
      beta_v1 += calc_beta_v1_one_loop(TRACE_STRUCT);
      beta_v2 += calc_beta_v2_one_loop(TRACE_STRUCT);

      if (get_loops() > 1) {
         beta_MassB += calc_beta_MassB_two_loop(TRACE_STRUCT);
         beta_MassG += calc_beta_MassG_two_loop(TRACE_STRUCT);
         beta_MassWB += calc_beta_MassWB_two_loop(TRACE_STRUCT);
         beta_Mu += calc_beta_Mu_two_loop(TRACE_STRUCT);
         beta_M122 += calc_beta_M122_two_loop(TRACE_STRUCT);
         beta_M112 += calc_beta_M112_two_loop(TRACE_STRUCT);
         beta_M222 += calc_beta_M222_two_loop(TRACE_STRUCT);
         beta_v1 += calc_beta_v1_two_loop(TRACE_STRUCT);
         beta_v2 += calc_beta_v2_two_loop(TRACE_STRUCT);

         if (get_loops() > 2) {

         }
      }
   }


   const HGTHDMIIMSSMBC_susy_parameters susy_betas(HGTHDMIIMSSMBC_susy_parameters::calc_beta());

   return HGTHDMIIMSSMBC_soft_parameters(susy_betas, beta_MassB, beta_MassG, beta_MassWB, beta_Mu, beta_M122, beta_M112, beta_M222, beta_v1, beta_v2);
}

void HGTHDMIIMSSMBC_soft_parameters::clear()
{
   HGTHDMIIMSSMBC_susy_parameters::clear();

   MassB = 0.;
   MassG = 0.;
   MassWB = 0.;
   Mu = 0.;
   M122 = 0.;
   M112 = 0.;
   M222 = 0.;
   v1 = 0.;
   v2 = 0.;

}

Eigen::ArrayXd HGTHDMIIMSSMBC_soft_parameters::get() const
{
   Eigen::ArrayXd pars(HGTHDMIIMSSMBC_susy_parameters::get());
   pars.conservativeResize(numberOfParameters);

   pars(41) = MassB;
   pars(42) = MassG;
   pars(43) = MassWB;
   pars(44) = Mu;
   pars(45) = M122;
   pars(46) = M112;
   pars(47) = M222;
   pars(48) = v1;
   pars(49) = v2;


   return pars;
}

void HGTHDMIIMSSMBC_soft_parameters::print(std::ostream& ostr) const
{
   HGTHDMIIMSSMBC_susy_parameters::print(ostr);
   ostr << "soft parameters at Q = " << get_scale() << ":\n";
   ostr << "MassB = " << MassB << '\n';
   ostr << "MassG = " << MassG << '\n';
   ostr << "MassWB = " << MassWB << '\n';
   ostr << "Mu = " << Mu << '\n';
   ostr << "M122 = " << M122 << '\n';
   ostr << "M112 = " << M112 << '\n';
   ostr << "M222 = " << M222 << '\n';
   ostr << "v1 = " << v1 << '\n';
   ostr << "v2 = " << v2 << '\n';

}

void HGTHDMIIMSSMBC_soft_parameters::set(const Eigen::ArrayXd& pars)
{
   HGTHDMIIMSSMBC_susy_parameters::set(pars);

   MassB = pars(41);
   MassG = pars(42);
   MassWB = pars(43);
   Mu = pars(44);
   M122 = pars(45);
   M112 = pars(46);
   M222 = pars(47);
   v1 = pars(48);
   v2 = pars(49);

}

void HGTHDMIIMSSMBC_soft_parameters::calc_soft_traces(Soft_traces& soft_traces) const
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

std::ostream& operator<<(std::ostream& ostr, const HGTHDMIIMSSMBC_soft_parameters& soft_pars)
{
   soft_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy