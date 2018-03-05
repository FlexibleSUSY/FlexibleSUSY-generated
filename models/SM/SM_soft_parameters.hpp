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

// File generated at Mon 5 Mar 2018 17:44:56

#ifndef SM_soft_parameters_H
#define SM_soft_parameters_H

#include "SM_susy_parameters.hpp"

#include <iosfwd>

namespace flexiblesusy {

#ifdef TRACE_STRUCT_TYPE
   #undef TRACE_STRUCT_TYPE
#endif
#define TRACE_STRUCT_TYPE Soft_traces

class SM_soft_parameters : public SM_susy_parameters {
public:
   explicit SM_soft_parameters(const SM_input_parameters& input_ = SM_input_parameters());
   SM_soft_parameters(const SM_susy_parameters& , double mu2_, double v_
);
   SM_soft_parameters(const SM_soft_parameters&) = default;
   SM_soft_parameters(SM_soft_parameters&&) = default;
   virtual ~SM_soft_parameters() = default;
   SM_soft_parameters& operator=(const SM_soft_parameters&) = default;
   SM_soft_parameters& operator=(SM_soft_parameters&&) = default;

   virtual Eigen::ArrayXd beta() const override;
   virtual Eigen::ArrayXd get() const override;
   virtual void print(std::ostream&) const override;
   virtual void set(const Eigen::ArrayXd&) override;

   SM_soft_parameters calc_beta() const;
   SM_soft_parameters calc_beta(int) const;
   virtual void clear() override;

   void set_mu2(double mu2_) { mu2 = mu2_; }
   void set_v(double v_) { v = v_; }

   double get_mu2() const { return mu2; }
   double get_v() const { return v; }


protected:
   double mu2{};
   double v{};


private:
   static const int numberOfParameters = 33;

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

   double calc_beta_mu2_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mu2_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mu2_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mu2_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_4_loop(const TRACE_STRUCT_TYPE&) const;

};

std::ostream& operator<<(std::ostream&, const SM_soft_parameters&);

} // namespace flexiblesusy

#undef TRACE_STRUCT_TYPE

#endif
