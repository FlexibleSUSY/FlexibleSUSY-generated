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

// File generated at Fri 20 Oct 2017 08:34:13

#ifndef MRSSMEFTHiggs_soft_parameters_H
#define MRSSMEFTHiggs_soft_parameters_H

#include "MRSSMEFTHiggs_susy_parameters.hpp"

#include <iosfwd>

namespace flexiblesusy {

#ifdef TRACE_STRUCT_TYPE
   #undef TRACE_STRUCT_TYPE
#endif
#define TRACE_STRUCT_TYPE Soft_traces

class MRSSMEFTHiggs_soft_parameters : public MRSSMEFTHiggs_susy_parameters {
public:
   explicit MRSSMEFTHiggs_soft_parameters(const MRSSMEFTHiggs_input_parameters& input_ = MRSSMEFTHiggs_input_parameters());
   MRSSMEFTHiggs_soft_parameters(const MRSSMEFTHiggs_susy_parameters& , double BMu_, double BMuD_, double BMuU_, const Eigen::Matrix<double,3,3>&
   mq2_, const Eigen::Matrix<double,3,3>& ml2_, double mHd2_, double mHu2_,
   const Eigen::Matrix<double,3,3>& md2_, const Eigen::Matrix<double,3,3>& mu2_
   , const Eigen::Matrix<double,3,3>& me2_, double mS2_, double mT2_, double
   moc2_, double mRd2_, double mRu2_, double MDBS_, double MDWBT_, double
   MDGoc_
);
   MRSSMEFTHiggs_soft_parameters(const MRSSMEFTHiggs_soft_parameters&) = default;
   MRSSMEFTHiggs_soft_parameters(MRSSMEFTHiggs_soft_parameters&&) = default;
   virtual ~MRSSMEFTHiggs_soft_parameters() = default;
   MRSSMEFTHiggs_soft_parameters& operator=(const MRSSMEFTHiggs_soft_parameters&) = default;
   MRSSMEFTHiggs_soft_parameters& operator=(MRSSMEFTHiggs_soft_parameters&&) = default;

   virtual Eigen::ArrayXd beta() const override;
   virtual Eigen::ArrayXd get() const override;
   virtual void print(std::ostream&) const override;
   virtual void set(const Eigen::ArrayXd&) override;

   MRSSMEFTHiggs_soft_parameters calc_beta() const;
   MRSSMEFTHiggs_soft_parameters calc_beta(int) const;
   virtual void clear() override;

   void set_BMu(double BMu_) { BMu = BMu_; }
   void set_BMuD(double BMuD_) { BMuD = BMuD_; }
   void set_BMuU(double BMuU_) { BMuU = BMuU_; }
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
   void set_mS2(double mS2_) { mS2 = mS2_; }
   void set_mT2(double mT2_) { mT2 = mT2_; }
   void set_moc2(double moc2_) { moc2 = moc2_; }
   void set_mRd2(double mRd2_) { mRd2 = mRd2_; }
   void set_mRu2(double mRu2_) { mRu2 = mRu2_; }
   void set_MDBS(double MDBS_) { MDBS = MDBS_; }
   void set_MDWBT(double MDWBT_) { MDWBT = MDWBT_; }
   void set_MDGoc(double MDGoc_) { MDGoc = MDGoc_; }

   double get_BMu() const { return BMu; }
   double get_BMuD() const { return BMuD; }
   double get_BMuU() const { return BMuU; }
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
   double get_mS2() const { return mS2; }
   double get_mT2() const { return mT2; }
   double get_moc2() const { return moc2; }
   double get_mRd2() const { return mRd2; }
   double get_mRu2() const { return mRu2; }
   double get_MDBS() const { return MDBS; }
   double get_MDWBT() const { return MDWBT; }
   double get_MDGoc() const { return MDGoc; }


protected:
   double BMu{};
   double BMuD{};
   double BMuU{};
   Eigen::Matrix<double,3,3> mq2{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> ml2{Eigen::Matrix<double,3,3>::Zero()};
   double mHd2{};
   double mHu2{};
   Eigen::Matrix<double,3,3> md2{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> mu2{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> me2{Eigen::Matrix<double,3,3>::Zero()};
   double mS2{};
   double mT2{};
   double moc2{};
   double mRd2{};
   double mRu2{};
   double MDBS{};
   double MDWBT{};
   double MDGoc{};


private:
   static const int numberOfParameters = 99;

   struct Soft_traces {
      double traceYdAdjYd{};
      double traceYeAdjYe{};
      double traceYuAdjYu{};
      double traceYdAdjYdYdAdjYd{};
      double traceYdAdjYuYuAdjYd{};
      double traceYeAdjYeYeAdjYe{};
      double traceYuAdjYuYuAdjYu{};
      double tracemd2YdAdjYd{};
      double traceme2YeAdjYe{};
      double traceml2AdjYeYe{};
      double tracemq2AdjYdYd{};
      double tracemq2AdjYuYu{};
      double tracemu2YuAdjYu{};
      double tracemd2YdAdjYdYdAdjYd{};
      double tracemd2YdAdjYuYuAdjYd{};
      double traceme2YeAdjYeYeAdjYe{};
      double traceml2AdjYeYeAdjYeYe{};
      double tracemq2AdjYdYdAdjYdYd{};
      double tracemq2AdjYdYdAdjYuYu{};
      double tracemq2AdjYuYuAdjYdYd{};
      double tracemu2YuAdjYdYdAdjYu{};
      double tracemq2AdjYuYuAdjYuYu{};
      double tracemu2YuAdjYuYuAdjYu{};
      double Tr11{};
      double Tr2U111{};
      double Tr31{};
      double Tr22{};
      double Tr23{};

   };
   Soft_traces calc_soft_traces(int) const;

   double calc_beta_BMu_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMu_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMu_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMuD_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMuD_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMuD_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMuU_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMuU_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMuU_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mq2_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mq2_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mq2_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_ml2_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_ml2_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_ml2_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHd2_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHd2_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHd2_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHu2_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHu2_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHu2_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_md2_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_md2_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_md2_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mu2_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mu2_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mu2_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_me2_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_me2_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_me2_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mS2_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mS2_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mS2_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mT2_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mT2_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mT2_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_moc2_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_moc2_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_moc2_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mRd2_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mRd2_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mRd2_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mRu2_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mRu2_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mRu2_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MDBS_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MDBS_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MDBS_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MDWBT_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MDWBT_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MDWBT_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MDGoc_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MDGoc_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MDGoc_3_loop(const TRACE_STRUCT_TYPE&) const;

};

std::ostream& operator<<(std::ostream&, const MRSSMEFTHiggs_soft_parameters&);

} // namespace flexiblesusy

#undef TRACE_STRUCT_TYPE

#endif
