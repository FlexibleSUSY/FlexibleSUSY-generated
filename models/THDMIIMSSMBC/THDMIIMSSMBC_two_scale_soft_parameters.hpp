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

// File generated at Tue 12 Jul 2016 10:31:03

#ifndef THDMIIMSSMBC_TWO_SCALE_soft_parameters_H
#define THDMIIMSSMBC_TWO_SCALE_soft_parameters_H

#include "THDMIIMSSMBC_two_scale_susy_parameters.hpp"

#include <iosfwd>

namespace flexiblesusy {

#ifdef TRACE_STRUCT_TYPE
   #undef TRACE_STRUCT_TYPE
#endif
#define TRACE_STRUCT_TYPE Soft_traces

class THDMIIMSSMBC_soft_parameters : public THDMIIMSSMBC_susy_parameters {
public:
   explicit THDMIIMSSMBC_soft_parameters(const THDMIIMSSMBC_input_parameters& input_ = THDMIIMSSMBC_input_parameters());
   THDMIIMSSMBC_soft_parameters(const THDMIIMSSMBC_susy_parameters& , double M122_, double M112_, double M222_, double v1_, double v2_
);
   virtual ~THDMIIMSSMBC_soft_parameters() {}
   virtual Eigen::ArrayXd beta() const;
   virtual Eigen::ArrayXd get() const;
   virtual void print(std::ostream&) const;
   virtual void set(const Eigen::ArrayXd&);

   THDMIIMSSMBC_soft_parameters calc_beta() const;
   virtual void clear();

   void set_M122(double M122_) { M122 = M122_; }
   void set_M112(double M112_) { M112 = M112_; }
   void set_M222(double M222_) { M222 = M222_; }
   void set_v1(double v1_) { v1 = v1_; }
   void set_v2(double v2_) { v2 = v2_; }

   double get_M122() const { return M122; }
   double get_M112() const { return M112; }
   double get_M222() const { return M222; }
   double get_v1() const { return v1; }
   double get_v2() const { return v2; }


protected:
   double M122;
   double M112;
   double M222;
   double v1;
   double v2;


private:
   static const int numberOfParameters = 42;

   struct Soft_traces {
      double traceYdAdjYd;
      double traceYeAdjYe;
      double traceYuAdjYu;
      double traceYdAdjYdYdAdjYd;
      double traceYdAdjYuYuAdjYd;
      double traceYeAdjYeYeAdjYe;
      double traceYuAdjYuYuAdjYu;

   };
   void calc_soft_traces(Soft_traces&) const;

   double calc_beta_M122_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M122_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M122_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M112_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M112_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M112_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M222_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M222_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_M222_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v1_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v1_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v1_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v2_three_loop(const TRACE_STRUCT_TYPE&) const;

};

std::ostream& operator<<(std::ostream&, const THDMIIMSSMBC_soft_parameters&);

} // namespace flexiblesusy

#undef TRACE_STRUCT_TYPE

#endif
