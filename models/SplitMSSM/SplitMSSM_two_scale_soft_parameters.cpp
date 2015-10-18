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

// File generated at Sun 18 Oct 2015 11:39:46

#include "SplitMSSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"
#include "functors.hpp"

#include <iostream>

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

SplitMSSM_soft_parameters::SplitMSSM_soft_parameters(const SplitMSSM_input_parameters& input_)
   : SplitMSSM_susy_parameters(input_)
   , MassB(0), MassG(0), MassWB(0), Mu(0), mu2(0), v(0)

{
   set_number_of_parameters(numberOfParameters);
}

SplitMSSM_soft_parameters::SplitMSSM_soft_parameters(
   const SplitMSSM_susy_parameters& susy_model
   , double MassB_, double MassG_, double MassWB_, double Mu_, double mu2_,
   double v_

)
   : SplitMSSM_susy_parameters(susy_model)
   , MassB(MassB_), MassG(MassG_), MassWB(MassWB_), Mu(Mu_), mu2(mu2_), v(v_)

{
   set_number_of_parameters(numberOfParameters);
}

Eigen::ArrayXd SplitMSSM_soft_parameters::beta() const
{
   return calc_beta().get().unaryExpr(Chop<double>(get_zero_threshold()));
}

SplitMSSM_soft_parameters SplitMSSM_soft_parameters::calc_beta() const
{
   Soft_traces soft_traces;
   calc_soft_traces(soft_traces);

   double beta_MassB(calc_beta_MassB_one_loop(TRACE_STRUCT));
   double beta_MassG(calc_beta_MassG_one_loop(TRACE_STRUCT));
   double beta_MassWB(calc_beta_MassWB_one_loop(TRACE_STRUCT));
   double beta_Mu(calc_beta_Mu_one_loop(TRACE_STRUCT));
   double beta_mu2(calc_beta_mu2_one_loop(TRACE_STRUCT));
   double beta_v(calc_beta_v_one_loop(TRACE_STRUCT));

   if (get_loops() > 1) {
      beta_MassB += calc_beta_MassB_two_loop(TRACE_STRUCT);
      beta_MassG += calc_beta_MassG_two_loop(TRACE_STRUCT);
      beta_MassWB += calc_beta_MassWB_two_loop(TRACE_STRUCT);
      beta_Mu += calc_beta_Mu_two_loop(TRACE_STRUCT);
      beta_mu2 += calc_beta_mu2_two_loop(TRACE_STRUCT);
      beta_v += calc_beta_v_two_loop(TRACE_STRUCT);

      if (get_loops() > 2) {

      }
   }


   const SplitMSSM_susy_parameters susy_betas(SplitMSSM_susy_parameters::calc_beta());

   return SplitMSSM_soft_parameters(susy_betas, beta_MassB, beta_MassG, beta_MassWB, beta_Mu, beta_mu2, beta_v);
}

void SplitMSSM_soft_parameters::clear()
{
   SplitMSSM_susy_parameters::clear();

   MassB = 0.;
   MassG = 0.;
   MassWB = 0.;
   Mu = 0.;
   mu2 = 0.;
   v = 0.;

}

Eigen::ArrayXd SplitMSSM_soft_parameters::get() const
{
   Eigen::ArrayXd pars(SplitMSSM_susy_parameters::get());
   pars.conservativeResize(numberOfParameters);

   pars(35) = MassB;
   pars(36) = MassG;
   pars(37) = MassWB;
   pars(38) = Mu;
   pars(39) = mu2;
   pars(40) = v;


   return pars;
}

void SplitMSSM_soft_parameters::print(std::ostream& ostr) const
{
   SplitMSSM_susy_parameters::print(ostr);
   ostr << "soft parameters:\n";
   ostr << "MassB = " << MassB << '\n';
   ostr << "MassG = " << MassG << '\n';
   ostr << "MassWB = " << MassWB << '\n';
   ostr << "Mu = " << Mu << '\n';
   ostr << "mu2 = " << mu2 << '\n';
   ostr << "v = " << v << '\n';

}

void SplitMSSM_soft_parameters::set(const Eigen::ArrayXd& pars)
{
   SplitMSSM_susy_parameters::set(pars);

   MassB = pars(35);
   MassG = pars(36);
   MassWB = pars(37);
   Mu = pars(38);
   mu2 = pars(39);
   v = pars(40);

}

void SplitMSSM_soft_parameters::calc_soft_traces(Soft_traces& soft_traces) const
{
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

std::ostream& operator<<(std::ostream& ostr, const SplitMSSM_soft_parameters& soft_pars)
{
   soft_pars.print(std::cout);
   return ostr;
}

} // namespace flexiblesusy
