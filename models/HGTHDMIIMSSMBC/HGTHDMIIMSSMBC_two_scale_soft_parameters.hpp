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

// File generated at Sun 28 Aug 2016 15:02:42

#ifndef HGTHDMIIMSSMBC_TWO_SCALE_soft_parameters_H
#define HGTHDMIIMSSMBC_TWO_SCALE_soft_parameters_H

#include "HGTHDMIIMSSMBC_two_scale_susy_parameters.hpp"

#include <iosfwd>

namespace flexiblesusy {

#ifdef TRACE_STRUCT_TYPE
   #undef TRACE_STRUCT_TYPE
#endif
#define TRACE_STRUCT_TYPE Soft_traces

class HGTHDMIIMSSMBC_soft_parameters : public HGTHDMIIMSSMBC_susy_parameters {
public:
   explicit HGTHDMIIMSSMBC_soft_parameters(const HGTHDMIIMSSMBC_input_parameters& input_ = HGTHDMIIMSSMBC_input_parameters());
   HGTHDMIIMSSMBC_soft_parameters(const HGTHDMIIMSSMBC_susy_parameters& , double MassB_, double MassG_, double MassWB_, double Mu_, double M122_,
   double M112_, double M222_, double v1_, double v2_
);
   virtual ~HGTHDMIIMSSMBC_soft_parameters() {}
   virtual Eigen::ArrayXd beta() const;
   virtual Eigen::ArrayXd get() const;
   virtual void print(std::ostream&) const;
   virtual void set(const Eigen::ArrayXd&);

   HGTHDMIIMSSMBC_soft_parameters calc_beta() const;
   HGTHDMIIMSSMBC_soft_parameters calc_beta(unsigned) const;
   virtual void clear();

   void set_MassB(double MassB_) { MassB = MassB_; }
   void set_MassG(double MassG_) { MassG = MassG_; }
   void set_MassWB(double MassWB_) { MassWB = MassWB_; }
   void set_Mu(double Mu_) { Mu = Mu_; }
   void set_M122(double M122_) { M122 = M122_; }
   void set_M112(double M112_) { M112 = M112_; }
   void set_M222(double M222_) { M222 = M222_; }
   void set_v1(double v1_) { v1 = v1_; }
   void set_v2(double v2_) { v2 = v2_; }

   double get_MassB() const { return MassB; }
   double get_MassG() const { return MassG; }
   double get_MassWB() const { return MassWB; }
   double get_Mu() const { return Mu; }
   double get_M122() const { return M122; }
   double get_M112() const { return M112; }
   double get_M222() const { return M222; }
   double get_v1() const { return v1; }
   double get_v2() const { return v2; }


protected:
   double MassB;
   double MassG;
   double MassWB;
   double Mu;
   double M122;
   double M112;
   double M222;
   double v1;
   double v2;


private:
   static const int numberOfParameters = 50;

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

   double calc_beta_MassB_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassB_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassB_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassG_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassG_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassG_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassWB_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassWB_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassWB_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Mu_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Mu_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Mu_three_loop(const TRACE_STRUCT_TYPE&) const;
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

std::ostream& operator<<(std::ostream&, const HGTHDMIIMSSMBC_soft_parameters&);

} // namespace flexiblesusy

#undef TRACE_STRUCT_TYPE

#endif
