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

// File generated at Sun 4 Aug 2019 19:57:44

#ifndef MSSMRHN_soft_parameters_H
#define MSSMRHN_soft_parameters_H

#include "MSSMRHN_susy_parameters.hpp"

#include <iosfwd>

namespace flexiblesusy {

#ifdef TRACE_STRUCT_TYPE
   #undef TRACE_STRUCT_TYPE
#endif
#define TRACE_STRUCT_TYPE Soft_traces

class MSSMRHN_soft_parameters : public MSSMRHN_susy_parameters {
public:
   explicit MSSMRHN_soft_parameters(const MSSMRHN_input_parameters& input_ = MSSMRHN_input_parameters());
   MSSMRHN_soft_parameters(const MSSMRHN_susy_parameters& , const Eigen::Matrix<double,3,3>& TYd_, const Eigen::Matrix<double,3,3>& TYe_,
   const Eigen::Matrix<double,3,3>& TYu_, const Eigen::Matrix<double,3,3>& TYv_
   , double BMu_, const Eigen::Matrix<double,3,3>& BMv_, const Eigen::Matrix<
   double,3,3>& mq2_, const Eigen::Matrix<double,3,3>& ml2_, double mHd2_,
   double mHu2_, const Eigen::Matrix<double,3,3>& md2_, const Eigen::Matrix<
   double,3,3>& mu2_, const Eigen::Matrix<double,3,3>& me2_, const Eigen::
   Matrix<double,3,3>& mv2_, double MassB_, double MassWB_, double MassG_);
   MSSMRHN_soft_parameters(const MSSMRHN_soft_parameters&) = default;
   MSSMRHN_soft_parameters(MSSMRHN_soft_parameters&&) = default;
   virtual ~MSSMRHN_soft_parameters() = default;
   MSSMRHN_soft_parameters& operator=(const MSSMRHN_soft_parameters&) = default;
   MSSMRHN_soft_parameters& operator=(MSSMRHN_soft_parameters&&) = default;

   virtual Eigen::ArrayXd beta() const override;
   virtual Eigen::ArrayXd get() const override;
   virtual void print(std::ostream&) const override;
   virtual void set(const Eigen::ArrayXd&) override;

   MSSMRHN_soft_parameters calc_beta() const;
   MSSMRHN_soft_parameters calc_beta(int) const;
   virtual void clear() override;

   void set_TYd(const Eigen::Matrix<double,3,3>& TYd_) { TYd = TYd_; }
   void set_TYd(int i, int k, const double& value) { TYd(i,k) = value; }
   void set_TYe(const Eigen::Matrix<double,3,3>& TYe_) { TYe = TYe_; }
   void set_TYe(int i, int k, const double& value) { TYe(i,k) = value; }
   void set_TYu(const Eigen::Matrix<double,3,3>& TYu_) { TYu = TYu_; }
   void set_TYu(int i, int k, const double& value) { TYu(i,k) = value; }
   void set_TYv(const Eigen::Matrix<double,3,3>& TYv_) { TYv = TYv_; }
   void set_TYv(int i, int k, const double& value) { TYv(i,k) = value; }
   void set_BMu(double BMu_) { BMu = BMu_; }
   void set_BMv(const Eigen::Matrix<double,3,3>& BMv_) { BMv = BMv_; }
   void set_BMv(int i, int k, const double& value) { BMv(i,k) = value; }
   void set_mq2(const Eigen::Matrix<double,3,3>& mq2_) { mq2 = mq2_; }
   void set_mq2(int i, int k, const double& value) { mq2(i,k) = value; }
   void set_ml2(const Eigen::Matrix<double,3,3>& ml2_) { ml2 = ml2_; }
   void set_ml2(int i, int k, const double& value) { ml2(i,k) = value; }
   void set_mHd2(double mHd2_) { mHd2 = mHd2_; }
   void set_mHu2(double mHu2_) { mHu2 = mHu2_; }
   void set_md2(const Eigen::Matrix<double,3,3>& md2_) { md2 = md2_; }
   void set_md2(int i, int k, const double& value) { md2(i,k) = value; }
   void set_mu2(const Eigen::Matrix<double,3,3>& mu2_) { mu2 = mu2_; }
   void set_mu2(int i, int k, const double& value) { mu2(i,k) = value; }
   void set_me2(const Eigen::Matrix<double,3,3>& me2_) { me2 = me2_; }
   void set_me2(int i, int k, const double& value) { me2(i,k) = value; }
   void set_mv2(const Eigen::Matrix<double,3,3>& mv2_) { mv2 = mv2_; }
   void set_mv2(int i, int k, const double& value) { mv2(i,k) = value; }
   void set_MassB(double MassB_) { MassB = MassB_; }
   void set_MassWB(double MassWB_) { MassWB = MassWB_; }
   void set_MassG(double MassG_) { MassG = MassG_; }

   const Eigen::Matrix<double,3,3>& get_TYd() const { return TYd; }
   double get_TYd(int i, int k) const { return TYd(i,k); }
   const Eigen::Matrix<double,3,3>& get_TYe() const { return TYe; }
   double get_TYe(int i, int k) const { return TYe(i,k); }
   const Eigen::Matrix<double,3,3>& get_TYu() const { return TYu; }
   double get_TYu(int i, int k) const { return TYu(i,k); }
   const Eigen::Matrix<double,3,3>& get_TYv() const { return TYv; }
   double get_TYv(int i, int k) const { return TYv(i,k); }
   double get_BMu() const { return BMu; }
   const Eigen::Matrix<double,3,3>& get_BMv() const { return BMv; }
   double get_BMv(int i, int k) const { return BMv(i,k); }
   const Eigen::Matrix<double,3,3>& get_mq2() const { return mq2; }
   double get_mq2(int i, int k) const { return mq2(i,k); }
   const Eigen::Matrix<double,3,3>& get_ml2() const { return ml2; }
   double get_ml2(int i, int k) const { return ml2(i,k); }
   double get_mHd2() const { return mHd2; }
   double get_mHu2() const { return mHu2; }
   const Eigen::Matrix<double,3,3>& get_md2() const { return md2; }
   double get_md2(int i, int k) const { return md2(i,k); }
   const Eigen::Matrix<double,3,3>& get_mu2() const { return mu2; }
   double get_mu2(int i, int k) const { return mu2(i,k); }
   const Eigen::Matrix<double,3,3>& get_me2() const { return me2; }
   double get_me2(int i, int k) const { return me2(i,k); }
   const Eigen::Matrix<double,3,3>& get_mv2() const { return mv2; }
   double get_mv2(int i, int k) const { return mv2(i,k); }
   double get_MassB() const { return MassB; }
   double get_MassWB() const { return MassWB; }
   double get_MassG() const { return MassG; }


protected:
   Eigen::Matrix<double,3,3> TYd{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> TYe{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> TYu{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> TYv{Eigen::Matrix<double,3,3>::Zero()};
   double BMu{};
   Eigen::Matrix<double,3,3> BMv{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> mq2{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> ml2{Eigen::Matrix<double,3,3>::Zero()};
   double mHd2{};
   double mHu2{};
   Eigen::Matrix<double,3,3> md2{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> mu2{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> me2{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> mv2{Eigen::Matrix<double,3,3>::Zero()};
   double MassB{};
   double MassWB{};
   double MassG{};


private:
   static const int numberOfParameters = 156;

   struct Soft_traces {
      double traceAdjYdTYd{};
      double traceAdjYeTYe{};
      double traceYdAdjYd{};
      double traceYeAdjYe{};
      double traceYdAdjYdTYdAdjYd{};
      double traceYdAdjYuTYuAdjYd{};
      double traceYeAdjYeTYeAdjYe{};
      double traceYeAdjYvTYvAdjYe{};
      double traceYuAdjYdTYdAdjYu{};
      double traceYvAdjYeTYeAdjYv{};
      double traceAdjYuTYu{};
      double traceAdjYvTYv{};
      double traceYuAdjYu{};
      double traceYvAdjYv{};
      double traceYdAdjYdYdAdjYd{};
      double traceYdAdjYuYuAdjYd{};
      double traceYeAdjYeYeAdjYe{};
      double traceYeAdjYvYvAdjYe{};
      double traceYuAdjYuTYuAdjYu{};
      double traceYvAdjYvTYvAdjYv{};
      double traceYuAdjYuYuAdjYu{};
      double traceYvAdjYvYvAdjYv{};
      double traceconjTYdTpTYd{};
      double traceconjTYeTpTYe{};
      double tracemd2YdAdjYd{};
      double traceme2YeAdjYe{};
      double traceml2AdjYeYe{};
      double tracemq2AdjYdYd{};
      double traceconjTYdTpYd{};
      double traceconjTYeTpYe{};
      double traceconjTYuTpTYu{};
      double traceconjTYvTpTYv{};
      double traceml2AdjYvYv{};
      double tracemq2AdjYuYu{};
      double tracemu2YuAdjYu{};
      double tracemv2YvAdjYv{};
      double traceconjTYuTpYu{};
      double traceconjTYvTpYv{};
      double traceYdAdjYdTYdAdjTYd{};
      double traceYdAdjYuTYuAdjTYd{};
      double traceYdAdjTYdTYdAdjYd{};
      double traceYdAdjTYuTYuAdjYd{};
      double traceYeAdjYeTYeAdjTYe{};
      double traceYeAdjYvTYvAdjTYe{};
      double traceYeAdjTYeTYeAdjYe{};
      double traceYeAdjTYvTYvAdjYe{};
      double traceYuAdjYdTYdAdjTYu{};
      double traceYuAdjTYdTYdAdjYu{};
      double traceYvAdjYeTYeAdjTYv{};
      double traceYvAdjTYeTYeAdjYv{};
      double tracemd2YdAdjYdYdAdjYd{};
      double tracemd2YdAdjYuYuAdjYd{};
      double traceme2YeAdjYeYeAdjYe{};
      double traceme2YeAdjYvYvAdjYe{};
      double traceml2AdjYeYeAdjYeYe{};
      double traceml2AdjYeYeAdjYvYv{};
      double traceml2AdjYvYvAdjYeYe{};
      double tracemq2AdjYdYdAdjYdYd{};
      double tracemq2AdjYdYdAdjYuYu{};
      double tracemq2AdjYuYuAdjYdYd{};
      double tracemu2YuAdjYdYdAdjYu{};
      double tracemv2YvAdjYeYeAdjYv{};
      double traceYuAdjYuTYuAdjTYu{};
      double traceYuAdjTYuTYuAdjYu{};
      double traceYvAdjYvTYvAdjTYv{};
      double traceYvAdjTYvTYvAdjYv{};
      double traceml2AdjYvYvAdjYvYv{};
      double tracemq2AdjYuYuAdjYuYu{};
      double tracemu2YuAdjYuYuAdjYu{};
      double tracemv2YvAdjYvYvAdjYv{};
      double Tr11{};
      double Tr2U111{};
      double Tr31{};
      double Tr22{};
      double Tr23{};

   };
   Soft_traces calc_soft_traces(int) const;

   Eigen::Matrix<double,3,3> calc_beta_TYd_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYd_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYd_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYd_4_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYd_5_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYe_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYe_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYe_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYe_4_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYe_5_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYu_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYu_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYu_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYu_4_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYu_5_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYv_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYv_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYv_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYv_4_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYv_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMu_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMu_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMu_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMu_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMu_5_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_BMv_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_BMv_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_BMv_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_BMv_4_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_BMv_5_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mq2_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mq2_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mq2_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mq2_4_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mq2_5_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_ml2_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_ml2_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_ml2_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_ml2_4_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_ml2_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHd2_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHd2_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHd2_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHd2_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHd2_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHu2_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHu2_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHu2_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHu2_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHu2_5_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_md2_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_md2_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_md2_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_md2_4_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_md2_5_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mu2_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mu2_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mu2_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mu2_4_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mu2_5_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_me2_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_me2_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_me2_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_me2_4_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_me2_5_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mv2_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mv2_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mv2_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mv2_4_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mv2_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassB_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassB_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassB_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassB_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassB_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassWB_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassWB_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassWB_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassWB_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassWB_5_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassG_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassG_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassG_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassG_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassG_5_loop(const TRACE_STRUCT_TYPE&) const;

};

std::ostream& operator<<(std::ostream&, const MSSMRHN_soft_parameters&);

} // namespace flexiblesusy

#undef TRACE_STRUCT_TYPE

#endif
