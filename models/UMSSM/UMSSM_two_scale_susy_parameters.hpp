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

// File generated at Tue 24 Feb 2015 17:35:12

#ifndef UMSSM_TWO_SCALE_susy_parameters_H
#define UMSSM_TWO_SCALE_susy_parameters_H

#include "betafunction.hpp"
#include "UMSSM_input_parameters.hpp"

#include <iosfwd>
#include <string>
#include <vector>
#include <Eigen/Core>

namespace flexiblesusy {

#ifdef TRACE_STRUCT_TYPE
   #undef TRACE_STRUCT_TYPE
#endif
#define TRACE_STRUCT_TYPE Susy_traces

class UMSSM_susy_parameters : public Beta_function {
public:
   explicit UMSSM_susy_parameters(const UMSSM_input_parameters& input_ = UMSSM_input_parameters());
   UMSSM_susy_parameters(double scale_, double loops_, double thresholds_, const UMSSM_input_parameters& input_, const Eigen::Matrix<double,3,3>& Yd_, const Eigen::Matrix<double,3,3>& Ye_
   , double Lambdax_, const Eigen::Matrix<double,3,3>& Yu_, double g1_, double
   g2_, double g3_, double gp_, double vd_, double vu_, double vS_
);
   virtual ~UMSSM_susy_parameters() {}
   virtual Eigen::ArrayXd beta() const;
   virtual const Eigen::ArrayXd get() const;
   virtual void print(std::ostream&) const;
   virtual void set(const Eigen::ArrayXd&);
   const UMSSM_input_parameters& get_input() const;
   void set_input_parameters(const UMSSM_input_parameters&);

   UMSSM_susy_parameters calc_beta() const;
   virtual void clear();

   void set_Yd(const Eigen::Matrix<double,3,3>& Yd_) { Yd = Yd_; }
   void set_Yd(int i, int k, double value) { Yd(i,k) = value; }
   void set_Ye(const Eigen::Matrix<double,3,3>& Ye_) { Ye = Ye_; }
   void set_Ye(int i, int k, double value) { Ye(i,k) = value; }
   void set_Lambdax(double Lambdax_) { Lambdax = Lambdax_; }
   void set_Yu(const Eigen::Matrix<double,3,3>& Yu_) { Yu = Yu_; }
   void set_Yu(int i, int k, double value) { Yu(i,k) = value; }
   void set_g1(double g1_) { g1 = g1_; }
   void set_g2(double g2_) { g2 = g2_; }
   void set_g3(double g3_) { g3 = g3_; }
   void set_gp(double gp_) { gp = gp_; }
   void set_vd(double vd_) { vd = vd_; }
   void set_vu(double vu_) { vu = vu_; }
   void set_vS(double vS_) { vS = vS_; }

   const Eigen::Matrix<double,3,3>& get_Yd() const { return Yd; }
   double get_Yd(int i, int k) const { return Yd(i,k); }
   const Eigen::Matrix<double,3,3>& get_Ye() const { return Ye; }
   double get_Ye(int i, int k) const { return Ye(i,k); }
   double get_Lambdax() const { return Lambdax; }
   const Eigen::Matrix<double,3,3>& get_Yu() const { return Yu; }
   double get_Yu(int i, int k) const { return Yu(i,k); }
   double get_g1() const { return g1; }
   double get_g2() const { return g2; }
   double get_g3() const { return g3; }
   double get_gp() const { return gp; }
   double get_vd() const { return vd; }
   double get_vu() const { return vu; }
   double get_vS() const { return vS; }

   Eigen::Matrix<double,3,3> get_SqSq() const;
   Eigen::Matrix<double,3,3> get_SlSl() const;
   double get_SHdSHd() const;
   double get_SHuSHu() const;
   Eigen::Matrix<double,3,3> get_SdRSdR() const;
   Eigen::Matrix<double,3,3> get_SuRSuR() const;
   Eigen::Matrix<double,3,3> get_SeRSeR() const;
   double get_SsRSsR() const;


protected:
   Eigen::Matrix<double,3,3> Yd;
   Eigen::Matrix<double,3,3> Ye;
   double Lambdax;
   Eigen::Matrix<double,3,3> Yu;
   double g1;
   double g2;
   double g3;
   double gp;
   double vd;
   double vu;
   double vS;

   UMSSM_input_parameters input;

private:
   static const int numberOfParameters = 35;

   struct Susy_traces {
      double traceYdAdjYd;
      double traceYeAdjYe;
      double traceYuAdjYu;
      double traceYdAdjYdYdAdjYd;
      double traceYdAdjYuYuAdjYd;
      double traceYeAdjYeYeAdjYe;
      double traceYuAdjYuYuAdjYu;

   };
   void calc_susy_traces(Susy_traces&) const;

   Eigen::Matrix<double,3,3> calc_beta_Yd_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yd_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Ye_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Ye_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambdax_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambdax_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yu_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yu_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g1_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g1_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g3_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g3_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_gp_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_gp_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_vd_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_vd_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_vu_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_vu_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_vS_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_vS_two_loop(const TRACE_STRUCT_TYPE&) const;

};

std::ostream& operator<<(std::ostream&, const UMSSM_susy_parameters&);

} // namespace flexiblesusy

#endif
