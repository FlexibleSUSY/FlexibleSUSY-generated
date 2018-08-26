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

// File generated at Sun 26 Aug 2018 14:09:50

#include "SplitMSSM_soft_parameters.hpp"
#include "config.h"
#include "global_thread_pool.hpp"
#include "wrappers.hpp"
#include "functors.hpp"

#include <iostream>

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces
#define TRACE_STRUCT_TYPE Soft_traces
#define CALCULATE_TRACES(l) calc_soft_traces(l);

const int SplitMSSM_soft_parameters::numberOfParameters;

SplitMSSM_soft_parameters::SplitMSSM_soft_parameters(const SplitMSSM_input_parameters& input_)
   : SplitMSSM_susy_parameters(input_)
{
   set_number_of_parameters(numberOfParameters);
}

SplitMSSM_soft_parameters::SplitMSSM_soft_parameters(
   const SplitMSSM_susy_parameters& susy_model
   , double MassB_, double MassG_, double MassWB_, double Mu_, double mu2_, double
    v_
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

SplitMSSM_soft_parameters SplitMSSM_soft_parameters::calc_beta(int loops) const
{
   double beta_MassB = 0.;
   double beta_MassG = 0.;
   double beta_MassWB = 0.;
   double beta_Mu = 0.;
   double beta_mu2 = 0.;
   double beta_v = 0.;

   if (loops > 0) {
      const auto TRACE_STRUCT = CALCULATE_TRACES(loops);

      beta_MassB += calc_beta_MassB_1_loop(TRACE_STRUCT);
      beta_MassG += calc_beta_MassG_1_loop(TRACE_STRUCT);
      beta_MassWB += calc_beta_MassWB_1_loop(TRACE_STRUCT);
      beta_Mu += calc_beta_Mu_1_loop(TRACE_STRUCT);
      beta_mu2 += calc_beta_mu2_1_loop(TRACE_STRUCT);
      beta_v += calc_beta_v_1_loop(TRACE_STRUCT);

      if (loops > 1) {
         beta_MassB += calc_beta_MassB_2_loop(TRACE_STRUCT);
         beta_MassG += calc_beta_MassG_2_loop(TRACE_STRUCT);
         beta_MassWB += calc_beta_MassWB_2_loop(TRACE_STRUCT);
         beta_Mu += calc_beta_Mu_2_loop(TRACE_STRUCT);
         beta_mu2 += calc_beta_mu2_2_loop(TRACE_STRUCT);
         beta_v += calc_beta_v_2_loop(TRACE_STRUCT);

         if (loops > 2) {
         #ifdef ENABLE_THREADS
            {


            }
         #else
         #endif

            if (loops > 3) {

            }
         }
      }
   }


   const SplitMSSM_susy_parameters susy_betas(SplitMSSM_susy_parameters::calc_beta(loops));

   return SplitMSSM_soft_parameters(susy_betas, beta_MassB, beta_MassG, beta_MassWB, beta_Mu, beta_mu2, beta_v);
}

SplitMSSM_soft_parameters SplitMSSM_soft_parameters::calc_beta() const
{
   return calc_beta(get_loops());
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
   ostr << "soft parameters at Q = " << get_scale() << ":\n";
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

SplitMSSM_soft_parameters::Soft_traces SplitMSSM_soft_parameters::calc_soft_traces(int loops) const
{
   Soft_traces soft_traces;

   if (loops > 0) {
      

      TRACE_STRUCT.traceYdAdjYd = Re((Yd*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYeAdjYe = Re((Ye*Ye.adjoint()).trace());
      TRACE_STRUCT.traceYuAdjYu = Re((Yu*Yu.adjoint()).trace());

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

   }

   if (loops > 2) {

   }

   return soft_traces;
}

std::ostream& operator<<(std::ostream& ostr, const SplitMSSM_soft_parameters& soft_pars)
{
   soft_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
