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

#ifndef SplitMSSM_soft_parameters_H
#define SplitMSSM_soft_parameters_H

#include "SplitMSSM_susy_parameters.hpp"

#include <iosfwd>

namespace flexiblesusy {

#ifdef TRACE_STRUCT_TYPE
   #undef TRACE_STRUCT_TYPE
#endif
#define TRACE_STRUCT_TYPE Soft_traces

class SplitMSSM_soft_parameters : public SplitMSSM_susy_parameters {
public:
   explicit SplitMSSM_soft_parameters(const SplitMSSM_input_parameters& input_ = SplitMSSM_input_parameters());
   SplitMSSM_soft_parameters(const SplitMSSM_susy_parameters& , double MassB_, double MassG_, double MassWB_, double Mu_, double mu2_, double
    v_);
   SplitMSSM_soft_parameters(const SplitMSSM_soft_parameters&) = default;
   SplitMSSM_soft_parameters(SplitMSSM_soft_parameters&&) = default;
   virtual ~SplitMSSM_soft_parameters() = default;
   SplitMSSM_soft_parameters& operator=(const SplitMSSM_soft_parameters&) = default;
   SplitMSSM_soft_parameters& operator=(SplitMSSM_soft_parameters&&) = default;

   virtual Eigen::ArrayXd beta() const override;
   virtual Eigen::ArrayXd get() const override;
   virtual void print(std::ostream&) const override;
   virtual void set(const Eigen::ArrayXd&) override;

   SplitMSSM_soft_parameters calc_beta() const;
   SplitMSSM_soft_parameters calc_beta(int) const;
   virtual void clear() override;

   void set_MassB(double MassB_) { MassB = MassB_; }
   void set_MassG(double MassG_) { MassG = MassG_; }
   void set_MassWB(double MassWB_) { MassWB = MassWB_; }
   void set_Mu(double Mu_) { Mu = Mu_; }
   void set_mu2(double mu2_) { mu2 = mu2_; }
   void set_v(double v_) { v = v_; }

   double get_MassB() const { return MassB; }
   double get_MassG() const { return MassG; }
   double get_MassWB() const { return MassWB; }
   double get_Mu() const { return Mu; }
   double get_mu2() const { return mu2; }
   double get_v() const { return v; }


protected:
   double MassB{};
   double MassG{};
   double MassWB{};
   double Mu{};
   double mu2{};
   double v{};


private:
   static const int numberOfParameters = 41;

   struct Soft_traces {
      double traceYdAdjYd{};
      double traceYeAdjYe{};
      double traceYuAdjYu{};
      double traceYdAdjYdYdAdjYd{};
      double traceYdAdjYuYuAdjYd{};
      double traceYeAdjYeYeAdjYe{};
      double traceYuAdjYuYuAdjYu{};

   };
   Soft_traces calc_soft_traces(int) const;

   double calc_beta_MassB_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassB_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassB_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassB_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassG_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassG_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassG_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassG_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassWB_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassWB_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassWB_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassWB_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Mu_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Mu_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Mu_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Mu_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mu2_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mu2_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mu2_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mu2_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_4_loop(const TRACE_STRUCT_TYPE&) const;

};

std::ostream& operator<<(std::ostream&, const SplitMSSM_soft_parameters&);

} // namespace flexiblesusy

#undef TRACE_STRUCT_TYPE

#endif
