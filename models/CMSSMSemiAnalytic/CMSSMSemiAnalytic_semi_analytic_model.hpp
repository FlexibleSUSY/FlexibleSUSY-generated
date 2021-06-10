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


/**
 * @file CMSSMSemiAnalytic_semi_analytic_model.hpp
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the semi_analytic solver by solving EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated with FlexibleSUSY 2.6.0 and SARAH 4.14.5 .
 */

#ifndef CMSSMSemiAnalytic_SEMI_ANALYTIC_MODEL_H
#define CMSSMSemiAnalytic_SEMI_ANALYTIC_MODEL_H

#include "CMSSMSemiAnalytic_model.hpp"
#include "CMSSMSemiAnalytic_model_slha.hpp"
#include "CMSSMSemiAnalytic_semi_analytic_solutions.hpp"

#include "model.hpp"

namespace flexiblesusy {

class Semi_analytic;
/**
 * @class CMSSMSemiAnalytic<Semi_analytic>
 * @brief model class with routines for determining masses and mixings and EWSB
 */
template<>
class CMSSMSemiAnalytic<Semi_analytic> : public Model, public CMSSMSemiAnalytic_slha {
public:
   explicit CMSSMSemiAnalytic(const CMSSMSemiAnalytic_input_parameters& input_ = CMSSMSemiAnalytic_input_parameters(), bool do_convert_masses_to_slha = true);
   explicit CMSSMSemiAnalytic(const CMSSMSemiAnalytic_slha&, bool do_convert_masses_to_slha = true);
   CMSSMSemiAnalytic(const CMSSMSemiAnalytic&) = default;
   CMSSMSemiAnalytic(CMSSMSemiAnalytic&&) = default;
   virtual ~CMSSMSemiAnalytic() = default;
   CMSSMSemiAnalytic& operator=(const CMSSMSemiAnalytic&) = default;
   CMSSMSemiAnalytic& operator=(CMSSMSemiAnalytic&&) = default;

   // interface functions
   virtual void calculate_spectrum();
   virtual void clear_problems();
   virtual std::string name() const;
   virtual void run_to(double scale, double eps = -1.0);
   virtual void print(std::ostream&) const;
   virtual void set_precision(double);

   /**
    * @brief returns the current values of the semi-analytic coefficients
    *
    * @return the current set of semi-analytic solutions
    */
   const CMSSMSemiAnalytic_semi_analytic_solutions& get_semi_analytic_solutions() const;

   /**
    * @brief returns the current values of the semi-analytic coefficients
    *
    * @return the current set of semi-analytic solutions
    */
   CMSSMSemiAnalytic_semi_analytic_solutions& get_semi_analytic_solutions();

   /**
    * @brief calculates the semi-analytic solutions for the soft parameters
    *
    * @param[in] input_scale the scale where the boundary conditions hold
    */
   void calculate_semi_analytic_solutions(double);

   double get_MassBCoeff1() const { return solutions.get_MassBCoeff1(); }
   double get_MassBCoeff2() const { return solutions.get_MassBCoeff2(); }

   double get_MassGCoeff1() const { return solutions.get_MassGCoeff1(); }
   double get_MassGCoeff2() const { return solutions.get_MassGCoeff2(); }

   double get_MassWBCoeff1() const { return solutions.get_MassWBCoeff1(); }
   double get_MassWBCoeff2() const { return solutions.get_MassWBCoeff2(); }

   const Eigen::Matrix<double,3,3>& get_TYdCoeff1() const { return solutions.
      get_TYdCoeff1(); }
   const Eigen::Matrix<double,3,3>& get_TYdCoeff2() const { return solutions.
      get_TYdCoeff2(); }

   const Eigen::Matrix<double,3,3>& get_TYeCoeff1() const { return solutions.
      get_TYeCoeff1(); }
   const Eigen::Matrix<double,3,3>& get_TYeCoeff2() const { return solutions.
      get_TYeCoeff2(); }

   const Eigen::Matrix<double,3,3>& get_TYuCoeff1() const { return solutions.
      get_TYuCoeff1(); }
   const Eigen::Matrix<double,3,3>& get_TYuCoeff2() const { return solutions.
      get_TYuCoeff2(); }

   const Eigen::Matrix<double,3,3>& get_md2Coeff1() const { return solutions.
      get_md2Coeff1(); }
   const Eigen::Matrix<double,3,3>& get_md2Coeff2() const { return solutions.
      get_md2Coeff2(); }
   const Eigen::Matrix<double,3,3>& get_md2Coeff3() const { return solutions.
      get_md2Coeff3(); }
   const Eigen::Matrix<double,3,3>& get_md2Coeff4() const { return solutions.
      get_md2Coeff4(); }

   const Eigen::Matrix<double,3,3>& get_me2Coeff1() const { return solutions.
      get_me2Coeff1(); }
   const Eigen::Matrix<double,3,3>& get_me2Coeff2() const { return solutions.
      get_me2Coeff2(); }
   const Eigen::Matrix<double,3,3>& get_me2Coeff3() const { return solutions.
      get_me2Coeff3(); }
   const Eigen::Matrix<double,3,3>& get_me2Coeff4() const { return solutions.
      get_me2Coeff4(); }

   double get_mHd2Coeff1() const { return solutions.get_mHd2Coeff1(); }
   double get_mHd2Coeff2() const { return solutions.get_mHd2Coeff2(); }
   double get_mHd2Coeff3() const { return solutions.get_mHd2Coeff3(); }
   double get_mHd2Coeff4() const { return solutions.get_mHd2Coeff4(); }

   double get_mHu2Coeff1() const { return solutions.get_mHu2Coeff1(); }
   double get_mHu2Coeff2() const { return solutions.get_mHu2Coeff2(); }
   double get_mHu2Coeff3() const { return solutions.get_mHu2Coeff3(); }
   double get_mHu2Coeff4() const { return solutions.get_mHu2Coeff4(); }

   const Eigen::Matrix<double,3,3>& get_ml2Coeff1() const { return solutions.
      get_ml2Coeff1(); }
   const Eigen::Matrix<double,3,3>& get_ml2Coeff2() const { return solutions.
      get_ml2Coeff2(); }
   const Eigen::Matrix<double,3,3>& get_ml2Coeff3() const { return solutions.
      get_ml2Coeff3(); }
   const Eigen::Matrix<double,3,3>& get_ml2Coeff4() const { return solutions.
      get_ml2Coeff4(); }

   const Eigen::Matrix<double,3,3>& get_mq2Coeff1() const { return solutions.
      get_mq2Coeff1(); }
   const Eigen::Matrix<double,3,3>& get_mq2Coeff2() const { return solutions.
      get_mq2Coeff2(); }
   const Eigen::Matrix<double,3,3>& get_mq2Coeff3() const { return solutions.
      get_mq2Coeff3(); }
   const Eigen::Matrix<double,3,3>& get_mq2Coeff4() const { return solutions.
      get_mq2Coeff4(); }

   const Eigen::Matrix<double,3,3>& get_mu2Coeff1() const { return solutions.
      get_mu2Coeff1(); }
   const Eigen::Matrix<double,3,3>& get_mu2Coeff2() const { return solutions.
      get_mu2Coeff2(); }
   const Eigen::Matrix<double,3,3>& get_mu2Coeff3() const { return solutions.
      get_mu2Coeff3(); }
   const Eigen::Matrix<double,3,3>& get_mu2Coeff4() const { return solutions.
      get_mu2Coeff4(); }

   double get_BMuCoeff1() const { return solutions.get_BMuCoeff1(); }
   double get_BMuCoeff2() const { return solutions.get_BMuCoeff2(); }
   double get_BMuCoeff3() const { return solutions.get_BMuCoeff3(); }

private:
   /// semi-analytic solutions for the model soft parameters
   CMSSMSemiAnalytic_semi_analytic_solutions solutions{};
};

std::ostream& operator<<(std::ostream&, const CMSSMSemiAnalytic<Semi_analytic>&);

} // namespace flexiblesusy

#endif
