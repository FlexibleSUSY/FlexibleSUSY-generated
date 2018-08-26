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

// File generated at Sun 26 Aug 2018 14:01:50

/**
 * @file CE6SSM_semi_analytic_solutions.hpp
 * @brief contains class for computing the semi-analytic solutions
 *
 * This file was generated at Sun 26 Aug 2018 14:01:50 with FlexibleSUSY
 * 2.2.0 (git commit: 8489097de2d6938a6da0149378457b5ad13d9425) and SARAH 4.13.0 .
 */

#ifndef CE6SSM_SEMI_ANALYTIC_SOLUTIONS_H
#define CE6SSM_SEMI_ANALYTIC_SOLUTIONS_H

#include "CE6SSM_mass_eigenstates.hpp"

#include <Eigen/Core>
#include <Eigen/SVD>

#include <map>
#include <vector>

namespace flexiblesusy {

/**
 * @class CE6SSM_semi_analytic_solutions
 * @brief class with routines for computing the semi-analytic solutions
 */
class CE6SSM_semi_analytic_solutions {
public:
   CE6SSM_semi_analytic_solutions();
   CE6SSM_semi_analytic_solutions(
      const CE6SSM_semi_analytic_solutions&) = default;
   CE6SSM_semi_analytic_solutions(
      CE6SSM_semi_analytic_solutions&&) = default;
   ~CE6SSM_semi_analytic_solutions() = default;
   CE6SSM_semi_analytic_solutions& operator=(
      const CE6SSM_semi_analytic_solutions&) = default;
   CE6SSM_semi_analytic_solutions& operator=(
      CE6SSM_semi_analytic_solutions&&) = default;

   /**
    * @brief returns the scale where the boundary conditions are imposed
    *
    * @return the boundary condition scale
    */
   double get_input_scale() const { return input_scale; }

   /**
    * @brief returns the scale at which the coefficients are calculated
    *
    * @return the scale where the coefficients apply
    */
   double get_output_scale() const { return output_scale; }

   /**
    * @brief sets the scale where the boundary conditions are imposed
    *
    * @param[in] s the value of the boundary condition scale
    */
   void set_input_scale(double s) { input_scale = s; }

   /**
    * @brief sets the scale at which the coefficients are calculated
    * @param[in] s the scale to calculate the coefficients at
    */
   void set_output_scale(double s) { output_scale = s; }

   double get_MassBCoeff1() const { return MassBCoeff1; }
   double get_MassBCoeff2() const { return MassBCoeff2; }
   double get_MassBpCoeff1() const { return MassBpCoeff1; }
   double get_MassBpCoeff2() const { return MassBpCoeff2; }
   double get_MassGCoeff1() const { return MassGCoeff1; }
   double get_MassGCoeff2() const { return MassGCoeff2; }
   double get_MassWBCoeff1() const { return MassWBCoeff1; }
   double get_MassWBCoeff2() const { return MassWBCoeff2; }
   const Eigen::Matrix<double,3,3>& get_TYdCoeff1() const { return TYdCoeff1; }
   double get_TYdCoeff1(int i, int k) const { return TYdCoeff1(i,k); }
   const Eigen::Matrix<double,3,3>& get_TYdCoeff2() const { return TYdCoeff2; }
   double get_TYdCoeff2(int i, int k) const { return TYdCoeff2(i,k); }
   const Eigen::Matrix<double,3,3>& get_TYeCoeff1() const { return TYeCoeff1; }
   double get_TYeCoeff1(int i, int k) const { return TYeCoeff1(i,k); }
   const Eigen::Matrix<double,3,3>& get_TYeCoeff2() const { return TYeCoeff2; }
   double get_TYeCoeff2(int i, int k) const { return TYeCoeff2(i,k); }
   const Eigen::Matrix<double,3,3>& get_TYuCoeff1() const { return TYuCoeff1; }
   double get_TYuCoeff1(int i, int k) const { return TYuCoeff1(i,k); }
   const Eigen::Matrix<double,3,3>& get_TYuCoeff2() const { return TYuCoeff2; }
   double get_TYuCoeff2(int i, int k) const { return TYuCoeff2(i,k); }
   const Eigen::Matrix<double,3,3>& get_TKappaCoeff1() const { return TKappaCoeff1
      ; }
   double get_TKappaCoeff1(int i, int k) const { return TKappaCoeff1(i,k); }
   const Eigen::Matrix<double,3,3>& get_TKappaCoeff2() const { return TKappaCoeff2
      ; }
   double get_TKappaCoeff2(int i, int k) const { return TKappaCoeff2(i,k); }
   double get_TLambdaxCoeff1() const { return TLambdaxCoeff1; }
   double get_TLambdaxCoeff2() const { return TLambdaxCoeff2; }
   const Eigen::Matrix<double,2,2>& get_TLambda12Coeff1() const { return
      TLambda12Coeff1; }
   double get_TLambda12Coeff1(int i, int k) const { return TLambda12Coeff1(i,k); }
   const Eigen::Matrix<double,2,2>& get_TLambda12Coeff2() const { return
      TLambda12Coeff2; }
   double get_TLambda12Coeff2(int i, int k) const { return TLambda12Coeff2(i,k); }
   const Eigen::Matrix<double,3,3>& get_md2Coeff1() const { return md2Coeff1; }
   double get_md2Coeff1(int i, int k) const { return md2Coeff1(i,k); }
   const Eigen::Matrix<double,3,3>& get_md2Coeff2() const { return md2Coeff2; }
   double get_md2Coeff2(int i, int k) const { return md2Coeff2(i,k); }
   const Eigen::Matrix<double,3,3>& get_md2Coeff3() const { return md2Coeff3; }
   double get_md2Coeff3(int i, int k) const { return md2Coeff3(i,k); }
   const Eigen::Matrix<double,3,3>& get_md2Coeff4() const { return md2Coeff4; }
   double get_md2Coeff4(int i, int k) const { return md2Coeff4(i,k); }
   const Eigen::Matrix<double,3,3>& get_mDx2Coeff1() const { return mDx2Coeff1; }
   double get_mDx2Coeff1(int i, int k) const { return mDx2Coeff1(i,k); }
   const Eigen::Matrix<double,3,3>& get_mDx2Coeff2() const { return mDx2Coeff2; }
   double get_mDx2Coeff2(int i, int k) const { return mDx2Coeff2(i,k); }
   const Eigen::Matrix<double,3,3>& get_mDx2Coeff3() const { return mDx2Coeff3; }
   double get_mDx2Coeff3(int i, int k) const { return mDx2Coeff3(i,k); }
   const Eigen::Matrix<double,3,3>& get_mDx2Coeff4() const { return mDx2Coeff4; }
   double get_mDx2Coeff4(int i, int k) const { return mDx2Coeff4(i,k); }
   const Eigen::Matrix<double,3,3>& get_mDxbar2Coeff1() const { return
      mDxbar2Coeff1; }
   double get_mDxbar2Coeff1(int i, int k) const { return mDxbar2Coeff1(i,k); }
   const Eigen::Matrix<double,3,3>& get_mDxbar2Coeff2() const { return
      mDxbar2Coeff2; }
   double get_mDxbar2Coeff2(int i, int k) const { return mDxbar2Coeff2(i,k); }
   const Eigen::Matrix<double,3,3>& get_mDxbar2Coeff3() const { return
      mDxbar2Coeff3; }
   double get_mDxbar2Coeff3(int i, int k) const { return mDxbar2Coeff3(i,k); }
   const Eigen::Matrix<double,3,3>& get_mDxbar2Coeff4() const { return
      mDxbar2Coeff4; }
   double get_mDxbar2Coeff4(int i, int k) const { return mDxbar2Coeff4(i,k); }
   const Eigen::Matrix<double,3,3>& get_me2Coeff1() const { return me2Coeff1; }
   double get_me2Coeff1(int i, int k) const { return me2Coeff1(i,k); }
   const Eigen::Matrix<double,3,3>& get_me2Coeff2() const { return me2Coeff2; }
   double get_me2Coeff2(int i, int k) const { return me2Coeff2(i,k); }
   const Eigen::Matrix<double,3,3>& get_me2Coeff3() const { return me2Coeff3; }
   double get_me2Coeff3(int i, int k) const { return me2Coeff3(i,k); }
   const Eigen::Matrix<double,3,3>& get_me2Coeff4() const { return me2Coeff4; }
   double get_me2Coeff4(int i, int k) const { return me2Coeff4(i,k); }
   const Eigen::Matrix<double,2,2>& get_mH1I2Coeff1() const { return mH1I2Coeff1;
      }
   double get_mH1I2Coeff1(int i, int k) const { return mH1I2Coeff1(i,k); }
   const Eigen::Matrix<double,2,2>& get_mH1I2Coeff2() const { return mH1I2Coeff2;
      }
   double get_mH1I2Coeff2(int i, int k) const { return mH1I2Coeff2(i,k); }
   const Eigen::Matrix<double,2,2>& get_mH1I2Coeff3() const { return mH1I2Coeff3;
      }
   double get_mH1I2Coeff3(int i, int k) const { return mH1I2Coeff3(i,k); }
   const Eigen::Matrix<double,2,2>& get_mH1I2Coeff4() const { return mH1I2Coeff4;
      }
   double get_mH1I2Coeff4(int i, int k) const { return mH1I2Coeff4(i,k); }
   const Eigen::Matrix<double,2,2>& get_mH2I2Coeff1() const { return mH2I2Coeff1;
      }
   double get_mH2I2Coeff1(int i, int k) const { return mH2I2Coeff1(i,k); }
   const Eigen::Matrix<double,2,2>& get_mH2I2Coeff2() const { return mH2I2Coeff2;
      }
   double get_mH2I2Coeff2(int i, int k) const { return mH2I2Coeff2(i,k); }
   const Eigen::Matrix<double,2,2>& get_mH2I2Coeff3() const { return mH2I2Coeff3;
      }
   double get_mH2I2Coeff3(int i, int k) const { return mH2I2Coeff3(i,k); }
   const Eigen::Matrix<double,2,2>& get_mH2I2Coeff4() const { return mH2I2Coeff4;
      }
   double get_mH2I2Coeff4(int i, int k) const { return mH2I2Coeff4(i,k); }
   double get_mHd2Coeff1() const { return mHd2Coeff1; }
   double get_mHd2Coeff2() const { return mHd2Coeff2; }
   double get_mHd2Coeff3() const { return mHd2Coeff3; }
   double get_mHd2Coeff4() const { return mHd2Coeff4; }
   double get_mHp2Coeff1() const { return mHp2Coeff1; }
   double get_mHp2Coeff2() const { return mHp2Coeff2; }
   double get_mHp2Coeff3() const { return mHp2Coeff3; }
   double get_mHp2Coeff4() const { return mHp2Coeff4; }
   double get_mHpbar2Coeff1() const { return mHpbar2Coeff1; }
   double get_mHpbar2Coeff2() const { return mHpbar2Coeff2; }
   double get_mHpbar2Coeff3() const { return mHpbar2Coeff3; }
   double get_mHpbar2Coeff4() const { return mHpbar2Coeff4; }
   double get_mHu2Coeff1() const { return mHu2Coeff1; }
   double get_mHu2Coeff2() const { return mHu2Coeff2; }
   double get_mHu2Coeff3() const { return mHu2Coeff3; }
   double get_mHu2Coeff4() const { return mHu2Coeff4; }
   const Eigen::Matrix<double,3,3>& get_ml2Coeff1() const { return ml2Coeff1; }
   double get_ml2Coeff1(int i, int k) const { return ml2Coeff1(i,k); }
   const Eigen::Matrix<double,3,3>& get_ml2Coeff2() const { return ml2Coeff2; }
   double get_ml2Coeff2(int i, int k) const { return ml2Coeff2(i,k); }
   const Eigen::Matrix<double,3,3>& get_ml2Coeff3() const { return ml2Coeff3; }
   double get_ml2Coeff3(int i, int k) const { return ml2Coeff3(i,k); }
   const Eigen::Matrix<double,3,3>& get_ml2Coeff4() const { return ml2Coeff4; }
   double get_ml2Coeff4(int i, int k) const { return ml2Coeff4(i,k); }
   const Eigen::Matrix<double,3,3>& get_mq2Coeff1() const { return mq2Coeff1; }
   double get_mq2Coeff1(int i, int k) const { return mq2Coeff1(i,k); }
   const Eigen::Matrix<double,3,3>& get_mq2Coeff2() const { return mq2Coeff2; }
   double get_mq2Coeff2(int i, int k) const { return mq2Coeff2(i,k); }
   const Eigen::Matrix<double,3,3>& get_mq2Coeff3() const { return mq2Coeff3; }
   double get_mq2Coeff3(int i, int k) const { return mq2Coeff3(i,k); }
   const Eigen::Matrix<double,3,3>& get_mq2Coeff4() const { return mq2Coeff4; }
   double get_mq2Coeff4(int i, int k) const { return mq2Coeff4(i,k); }
   double get_ms2Coeff1() const { return ms2Coeff1; }
   double get_ms2Coeff2() const { return ms2Coeff2; }
   double get_ms2Coeff3() const { return ms2Coeff3; }
   double get_ms2Coeff4() const { return ms2Coeff4; }
   const Eigen::Matrix<double,2,2>& get_msI2Coeff1() const { return msI2Coeff1; }
   double get_msI2Coeff1(int i, int k) const { return msI2Coeff1(i,k); }
   const Eigen::Matrix<double,2,2>& get_msI2Coeff2() const { return msI2Coeff2; }
   double get_msI2Coeff2(int i, int k) const { return msI2Coeff2(i,k); }
   const Eigen::Matrix<double,2,2>& get_msI2Coeff3() const { return msI2Coeff3; }
   double get_msI2Coeff3(int i, int k) const { return msI2Coeff3(i,k); }
   const Eigen::Matrix<double,2,2>& get_msI2Coeff4() const { return msI2Coeff4; }
   double get_msI2Coeff4(int i, int k) const { return msI2Coeff4(i,k); }
   const Eigen::Matrix<double,3,3>& get_mu2Coeff1() const { return mu2Coeff1; }
   double get_mu2Coeff1(int i, int k) const { return mu2Coeff1(i,k); }
   const Eigen::Matrix<double,3,3>& get_mu2Coeff2() const { return mu2Coeff2; }
   double get_mu2Coeff2(int i, int k) const { return mu2Coeff2(i,k); }
   const Eigen::Matrix<double,3,3>& get_mu2Coeff3() const { return mu2Coeff3; }
   double get_mu2Coeff3(int i, int k) const { return mu2Coeff3(i,k); }
   const Eigen::Matrix<double,3,3>& get_mu2Coeff4() const { return mu2Coeff4; }
   double get_mu2Coeff4(int i, int k) const { return mu2Coeff4(i,k); }
   double get_BMuPrCoeff1() const { return BMuPrCoeff1; }
   double get_BMuPrCoeff2() const { return BMuPrCoeff2; }
   double get_BMuPrCoeff3() const { return BMuPrCoeff3; }

   /**
    * @brief calculates semi-analytic coefficients for a model
    *
    * @param[in] model model to calculate semi-analytic coefficients for
    */
   void calculate_coefficients(const CE6SSM_mass_eigenstates&);

   /**
    * @brief applies the semi-analytic solutions to the model
    *
    * @param[out] model model to apply semi-analytic solutions to
    */
   void evaluate_solutions(CE6SSM_mass_eigenstates&) const;

private:
   struct Boundary_values {
      double Azero{};
      double m12{};
      double m0Sq{};
      double BMuPrimeInput{};
      double MuPr{};

      
   };

   struct Model_data {
      Boundary_values boundary_values{};
      CE6SSM_mass_eigenstates model{};
      std::vector<int> basis_sets{};
   };

   using Data_vector_t = std::vector<const Model_data*>;

   double input_scale{0.};  ///< scale at which boundary conditions hold
   double output_scale{0.}; ///< scale at which coefficients are calculated
   std::array<Model_data,7> trial_data{};

   // semi-analytic solution coefficients
   double MassBCoeff1{};
   double MassBCoeff2{};
   double MassBpCoeff1{};
   double MassBpCoeff2{};
   double MassGCoeff1{};
   double MassGCoeff2{};
   double MassWBCoeff1{};
   double MassWBCoeff2{};
   Eigen::Matrix<double,3,3> TYdCoeff1{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> TYdCoeff2{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> TYeCoeff1{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> TYeCoeff2{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> TYuCoeff1{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> TYuCoeff2{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> TKappaCoeff1{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> TKappaCoeff2{Eigen::Matrix<double,3,3>::Zero()};
   double TLambdaxCoeff1{};
   double TLambdaxCoeff2{};
   Eigen::Matrix<double,2,2> TLambda12Coeff1{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> TLambda12Coeff2{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,3,3> md2Coeff1{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> md2Coeff2{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> md2Coeff3{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> md2Coeff4{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> mDx2Coeff1{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> mDx2Coeff2{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> mDx2Coeff3{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> mDx2Coeff4{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> mDxbar2Coeff1{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> mDxbar2Coeff2{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> mDxbar2Coeff3{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> mDxbar2Coeff4{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> me2Coeff1{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> me2Coeff2{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> me2Coeff3{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> me2Coeff4{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,2,2> mH1I2Coeff1{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> mH1I2Coeff2{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> mH1I2Coeff3{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> mH1I2Coeff4{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> mH2I2Coeff1{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> mH2I2Coeff2{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> mH2I2Coeff3{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> mH2I2Coeff4{Eigen::Matrix<double,2,2>::Zero()};
   double mHd2Coeff1{};
   double mHd2Coeff2{};
   double mHd2Coeff3{};
   double mHd2Coeff4{};
   double mHp2Coeff1{};
   double mHp2Coeff2{};
   double mHp2Coeff3{};
   double mHp2Coeff4{};
   double mHpbar2Coeff1{};
   double mHpbar2Coeff2{};
   double mHpbar2Coeff3{};
   double mHpbar2Coeff4{};
   double mHu2Coeff1{};
   double mHu2Coeff2{};
   double mHu2Coeff3{};
   double mHu2Coeff4{};
   Eigen::Matrix<double,3,3> ml2Coeff1{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> ml2Coeff2{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> ml2Coeff3{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> ml2Coeff4{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> mq2Coeff1{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> mq2Coeff2{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> mq2Coeff3{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> mq2Coeff4{Eigen::Matrix<double,3,3>::Zero()};
   double ms2Coeff1{};
   double ms2Coeff2{};
   double ms2Coeff3{};
   double ms2Coeff4{};
   Eigen::Matrix<double,2,2> msI2Coeff1{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> msI2Coeff2{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> msI2Coeff3{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> msI2Coeff4{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,3,3> mu2Coeff1{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> mu2Coeff2{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> mu2Coeff3{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> mu2Coeff4{Eigen::Matrix<double,3,3>::Zero()};
   double BMuPrCoeff1{};
   double BMuPrCoeff2{};
   double BMuPrCoeff3{};


   void initialize_trial_values();
   void calculate_trial_data(const CE6SSM_mass_eigenstates&);
   void set_to_boundary_values(CE6SSM_mass_eigenstates&,
                               const Boundary_values&) const;
   CE6SSM_mass_eigenstates run_to_output_scale(
      const CE6SSM_mass_eigenstates&, const Boundary_values&) const;
   std::map<int,Data_vector_t> create_datasets() const;

   template <class MatrixType, int BasisSize, class BasisEvaluator>
   Eigen::JacobiSVD<MatrixType> create_solver(
      const Data_vector_t&, const BasisEvaluator&) const;

   void calculate_MassB_coefficients(const Eigen::JacobiSVD<Eigen::MatrixXd>&, const Data_vector_t&);
   void calculate_MassBp_coefficients(const Eigen::JacobiSVD<Eigen::MatrixXd>&, const Data_vector_t&);
   void calculate_MassG_coefficients(const Eigen::JacobiSVD<Eigen::MatrixXd>&, const Data_vector_t&);
   void calculate_MassWB_coefficients(const Eigen::JacobiSVD<Eigen::MatrixXd>&, const Data_vector_t&);
   void calculate_TYd_coefficients(const Eigen::JacobiSVD<Eigen::MatrixXd>&, const Data_vector_t&);
   void calculate_TYe_coefficients(const Eigen::JacobiSVD<Eigen::MatrixXd>&, const Data_vector_t&);
   void calculate_TYu_coefficients(const Eigen::JacobiSVD<Eigen::MatrixXd>&, const Data_vector_t&);
   void calculate_TKappa_coefficients(const Eigen::JacobiSVD<Eigen::MatrixXd>&, const Data_vector_t&);
   void calculate_TLambdax_coefficients(const Eigen::JacobiSVD<Eigen::MatrixXd>&, const Data_vector_t&);
   void calculate_TLambda12_coefficients(const Eigen::JacobiSVD<Eigen::MatrixXd>&, const Data_vector_t&);
   void calculate_md2_coefficients(const Eigen::JacobiSVD<Eigen::MatrixXd>&, const Data_vector_t&);
   void calculate_mDx2_coefficients(const Eigen::JacobiSVD<Eigen::MatrixXd>&, const Data_vector_t&);
   void calculate_mDxbar2_coefficients(const Eigen::JacobiSVD<Eigen::MatrixXd>&, const Data_vector_t&);
   void calculate_me2_coefficients(const Eigen::JacobiSVD<Eigen::MatrixXd>&, const Data_vector_t&);
   void calculate_mH1I2_coefficients(const Eigen::JacobiSVD<Eigen::MatrixXd>&, const Data_vector_t&);
   void calculate_mH2I2_coefficients(const Eigen::JacobiSVD<Eigen::MatrixXd>&, const Data_vector_t&);
   void calculate_mHd2_coefficients(const Eigen::JacobiSVD<Eigen::MatrixXd>&, const Data_vector_t&);
   void calculate_mHp2_coefficients(const Eigen::JacobiSVD<Eigen::MatrixXd>&, const Data_vector_t&);
   void calculate_mHpbar2_coefficients(const Eigen::JacobiSVD<Eigen::MatrixXd>&, const Data_vector_t&);
   void calculate_mHu2_coefficients(const Eigen::JacobiSVD<Eigen::MatrixXd>&, const Data_vector_t&);
   void calculate_ml2_coefficients(const Eigen::JacobiSVD<Eigen::MatrixXd>&, const Data_vector_t&);
   void calculate_mq2_coefficients(const Eigen::JacobiSVD<Eigen::MatrixXd>&, const Data_vector_t&);
   void calculate_ms2_coefficients(const Eigen::JacobiSVD<Eigen::MatrixXd>&, const Data_vector_t&);
   void calculate_msI2_coefficients(const Eigen::JacobiSVD<Eigen::MatrixXd>&, const Data_vector_t&);
   void calculate_mu2_coefficients(const Eigen::JacobiSVD<Eigen::MatrixXd>&, const Data_vector_t&);
   void calculate_BMuPr_coefficients(const Eigen::JacobiSVD<Eigen::MatrixXd>&, const Data_vector_t&);

};

template <class MatrixType, int BasisSize, class BasisEvaluator>
Eigen::JacobiSVD<MatrixType> CE6SSM_semi_analytic_solutions::create_solver(
   const CE6SSM_semi_analytic_solutions::Data_vector_t& data,
   const BasisEvaluator& ev) const
{
   const std::size_t n = data.size();
   MatrixType lhs(n, BasisSize);
   for (std::size_t i = 0; i < n; ++i) {
      lhs.row(i) = ev(data[i]->boundary_values);
   }
   return Eigen::JacobiSVD<MatrixType>(
      lhs, Eigen::ComputeThinU | Eigen::ComputeThinV);
}

} // namespace flexiblesusy

#endif
