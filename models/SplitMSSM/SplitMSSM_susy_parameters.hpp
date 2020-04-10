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

// File generated at Fri 10 Apr 2020 19:53:06

#ifndef SplitMSSM_susy_parameters_H
#define SplitMSSM_susy_parameters_H

#include "betafunction.hpp"
#include "SplitMSSM_input_parameters.hpp"

#include <iosfwd>
#include <string>
#include <vector>
#include <Eigen/Core>

namespace flexiblesusy {

#ifdef TRACE_STRUCT_TYPE
   #undef TRACE_STRUCT_TYPE
#endif
#define TRACE_STRUCT_TYPE Susy_traces

class SplitMSSM_susy_parameters : public Beta_function {
public:
   explicit SplitMSSM_susy_parameters(const SplitMSSM_input_parameters& input_ = SplitMSSM_input_parameters());
   SplitMSSM_susy_parameters(double scale_, int loops_, int thresholds_, const SplitMSSM_input_parameters& input_, double g1_, double g2_, double g3_, double Lambdax_, const Eigen::Matrix<
   double,3,3>& Yu_, const Eigen::Matrix<double,3,3>& Yd_, const Eigen::Matrix<
   double,3,3>& Ye_, double gYd_, double g2d_, double gYu_, double g2u_);
   SplitMSSM_susy_parameters(const SplitMSSM_susy_parameters&) = default;
   SplitMSSM_susy_parameters(SplitMSSM_susy_parameters&&) = default;
   virtual ~SplitMSSM_susy_parameters() = default;
   SplitMSSM_susy_parameters& operator=(const SplitMSSM_susy_parameters&) = default;
   SplitMSSM_susy_parameters& operator=(SplitMSSM_susy_parameters&&) = default;

   virtual Eigen::ArrayXd beta() const override;
   virtual Eigen::ArrayXd get() const override;
   virtual void print(std::ostream&) const;
   virtual void set(const Eigen::ArrayXd&) override;
   const SplitMSSM_input_parameters& get_input() const;
   SplitMSSM_input_parameters& get_input();
   void set_input_parameters(const SplitMSSM_input_parameters&);

   SplitMSSM_susy_parameters calc_beta() const;
   SplitMSSM_susy_parameters calc_beta(int) const;
   virtual void clear();

   void set_g1(double g1_) { g1 = g1_; }
   void set_g2(double g2_) { g2 = g2_; }
   void set_g3(double g3_) { g3 = g3_; }
   void set_Lambdax(double Lambdax_) { Lambdax = Lambdax_; }
   void set_Yu(const Eigen::Matrix<double,3,3>& Yu_) { Yu = Yu_; }
   void set_Yu(int i, int k, const double& value) { Yu(i,k) = value; }
   void set_Yd(const Eigen::Matrix<double,3,3>& Yd_) { Yd = Yd_; }
   void set_Yd(int i, int k, const double& value) { Yd(i,k) = value; }
   void set_Ye(const Eigen::Matrix<double,3,3>& Ye_) { Ye = Ye_; }
   void set_Ye(int i, int k, const double& value) { Ye(i,k) = value; }
   void set_gYd(double gYd_) { gYd = gYd_; }
   void set_g2d(double g2d_) { g2d = g2d_; }
   void set_gYu(double gYu_) { gYu = gYu_; }
   void set_g2u(double g2u_) { g2u = g2u_; }

   double get_g1() const { return g1; }
   double get_g2() const { return g2; }
   double get_g3() const { return g3; }
   double get_Lambdax() const { return Lambdax; }
   const Eigen::Matrix<double,3,3>& get_Yu() const { return Yu; }
   double get_Yu(int i, int k) const { return Yu(i,k); }
   const Eigen::Matrix<double,3,3>& get_Yd() const { return Yd; }
   double get_Yd(int i, int k) const { return Yd(i,k); }
   const Eigen::Matrix<double,3,3>& get_Ye() const { return Ye; }
   double get_Ye(int i, int k) const { return Ye(i,k); }
   double get_gYd() const { return gYd; }
   double get_g2d() const { return g2d; }
   double get_gYu() const { return gYu; }
   double get_g2u() const { return g2u; }



protected:
   double g1{};
   double g2{};
   double g3{};
   double Lambdax{};
   Eigen::Matrix<double,3,3> Yu{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> Yd{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> Ye{Eigen::Matrix<double,3,3>::Zero()};
   double gYd{};
   double g2d{};
   double gYu{};
   double g2u{};

   SplitMSSM_input_parameters input{};

private:
   static const int numberOfParameters = 35;

   struct Susy_traces {
      double traceYdAdjYd{};
      double traceYeAdjYe{};
      double traceYuAdjYu{};
      double traceYdAdjYdYdAdjYd{};
      double traceYeAdjYeYeAdjYe{};
      double traceYuAdjYuYuAdjYu{};
      double traceYdAdjYuYuAdjYd{};
      double traceYdAdjYdYdAdjYdYdAdjYd{};
      double traceYdAdjYdYdAdjYuYuAdjYd{};
      double traceYdAdjYuYuAdjYdYdAdjYd{};
      double traceYdAdjYuYuAdjYuYuAdjYd{};
      double traceYeAdjYeYeAdjYeYeAdjYe{};
      double traceYuAdjYuYuAdjYuYuAdjYu{};

   };
   Susy_traces calc_susy_traces(int) const;

   double calc_beta_g1_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g1_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g1_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g1_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g1_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g2_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g2_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g2_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g2_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g2_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g3_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g3_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g3_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g3_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g3_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambdax_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambdax_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambdax_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambdax_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambdax_5_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yu_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yu_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yu_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yu_4_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yu_5_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yd_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yd_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yd_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yd_4_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yd_5_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Ye_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Ye_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Ye_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Ye_4_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_Ye_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_gYd_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_gYd_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_gYd_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_gYd_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_gYd_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g2d_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g2d_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g2d_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g2d_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g2d_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_gYu_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_gYu_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_gYu_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_gYu_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_gYu_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g2u_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g2u_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g2u_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g2u_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g2u_5_loop(const TRACE_STRUCT_TYPE&) const;

};

std::ostream& operator<<(std::ostream&, const SplitMSSM_susy_parameters&);

#undef TRACE_STRUCT_TYPE

} // namespace flexiblesusy

#endif
