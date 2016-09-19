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

// File generated at Mon 19 Sep 2016 09:53:53

#ifndef MRSSM_TWO_SCALE_soft_parameters_H
#define MRSSM_TWO_SCALE_soft_parameters_H

#include "MRSSM_two_scale_susy_parameters.hpp"

#include <iosfwd>

namespace flexiblesusy {

#ifdef TRACE_STRUCT_TYPE
   #undef TRACE_STRUCT_TYPE
#endif
#define TRACE_STRUCT_TYPE Soft_traces

class MRSSM_soft_parameters : public MRSSM_susy_parameters {
public:
   explicit MRSSM_soft_parameters(const MRSSM_input_parameters& input_ = MRSSM_input_parameters());
   MRSSM_soft_parameters(const MRSSM_susy_parameters& , double BMu_, double BMuD_, double BMuU_, const Eigen::Matrix<double,3,3>&
   mq2_, const Eigen::Matrix<double,3,3>& ml2_, double mHd2_, double mHu2_,
   const Eigen::Matrix<double,3,3>& md2_, const Eigen::Matrix<double,3,3>& mu2_
   , const Eigen::Matrix<double,3,3>& me2_, double mS2_, double mT2_, double
   moc2_, double mRd2_, double mRu2_, double MDBS_, double MDWBT_, double
   MDGoc_
);
   virtual ~MRSSM_soft_parameters() {}
   virtual Eigen::ArrayXd beta() const;
   virtual Eigen::ArrayXd get() const;
   virtual void print(std::ostream&) const;
   virtual void set(const Eigen::ArrayXd&);

   MRSSM_soft_parameters calc_beta() const;
   MRSSM_soft_parameters calc_beta(unsigned) const;
   virtual void clear();

   void set_BMu(double BMu_) { BMu = BMu_; }
   void set_BMuD(double BMuD_) { BMuD = BMuD_; }
   void set_BMuU(double BMuU_) { BMuU = BMuU_; }
   void set_mq2(const Eigen::Matrix<double,3,3>& mq2_) { mq2 = mq2_; }
   void set_mq2(int i, int k, double value) { mq2(i,k) = value; }
   void set_ml2(const Eigen::Matrix<double,3,3>& ml2_) { ml2 = ml2_; }
   void set_ml2(int i, int k, double value) { ml2(i,k) = value; }
   void set_mHd2(double mHd2_) { mHd2 = mHd2_; }
   void set_mHu2(double mHu2_) { mHu2 = mHu2_; }
   void set_md2(const Eigen::Matrix<double,3,3>& md2_) { md2 = md2_; }
   void set_md2(int i, int k, double value) { md2(i,k) = value; }
   void set_mu2(const Eigen::Matrix<double,3,3>& mu2_) { mu2 = mu2_; }
   void set_mu2(int i, int k, double value) { mu2(i,k) = value; }
   void set_me2(const Eigen::Matrix<double,3,3>& me2_) { me2 = me2_; }
   void set_me2(int i, int k, double value) { me2(i,k) = value; }
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
   double BMu;
   double BMuD;
   double BMuU;
   Eigen::Matrix<double,3,3> mq2;
   Eigen::Matrix<double,3,3> ml2;
   double mHd2;
   double mHu2;
   Eigen::Matrix<double,3,3> md2;
   Eigen::Matrix<double,3,3> mu2;
   Eigen::Matrix<double,3,3> me2;
   double mS2;
   double mT2;
   double moc2;
   double mRd2;
   double mRu2;
   double MDBS;
   double MDWBT;
   double MDGoc;


private:
   static const int numberOfParameters = 99;

   struct Soft_traces {
      double traceYdAdjYd;
      double traceYeAdjYe;
      double traceYuAdjYu;
      double traceYdAdjYdYdAdjYd;
      double traceYdAdjYuYuAdjYd;
      double traceYeAdjYeYeAdjYe;
      double traceYuAdjYuYuAdjYu;
      double tracemd2YdAdjYd;
      double traceme2YeAdjYe;
      double traceml2AdjYeYe;
      double tracemq2AdjYdYd;
      double tracemq2AdjYuYu;
      double tracemu2YuAdjYu;
      double tracemd2YdAdjYdYdAdjYd;
      double tracemd2YdAdjYuYuAdjYd;
      double traceme2YeAdjYeYeAdjYe;
      double traceml2AdjYeYeAdjYeYe;
      double tracemq2AdjYdYdAdjYdYd;
      double tracemq2AdjYdYdAdjYuYu;
      double tracemq2AdjYuYuAdjYdYd;
      double tracemu2YuAdjYdYdAdjYu;
      double tracemq2AdjYuYuAdjYuYu;
      double tracemu2YuAdjYuYuAdjYu;
      double Tr11;
      double Tr2U111;
      double Tr31;
      double Tr22;
      double Tr23;

   };
   void calc_soft_traces(Soft_traces&) const;

   double calc_beta_BMu_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMu_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMu_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMuD_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMuD_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMuD_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMuU_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMuU_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMuU_three_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mq2_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mq2_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mq2_three_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_ml2_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_ml2_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_ml2_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHd2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHd2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHd2_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHu2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHu2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHu2_three_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_md2_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_md2_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_md2_three_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mu2_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mu2_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mu2_three_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_me2_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_me2_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_me2_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mS2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mS2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mS2_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mT2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mT2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mT2_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_moc2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_moc2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_moc2_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mRd2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mRd2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mRd2_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mRu2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mRu2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mRu2_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MDBS_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MDBS_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MDBS_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MDWBT_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MDWBT_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MDWBT_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MDGoc_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MDGoc_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MDGoc_three_loop(const TRACE_STRUCT_TYPE&) const;

};

std::ostream& operator<<(std::ostream&, const MRSSM_soft_parameters&);

} // namespace flexiblesusy

#undef TRACE_STRUCT_TYPE

#endif
