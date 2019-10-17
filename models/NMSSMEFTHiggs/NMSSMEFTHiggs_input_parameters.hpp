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

// File generated at Wed 16 Oct 2019 19:36:33

#ifndef NMSSMEFTHiggs_INPUT_PARAMETERS_H
#define NMSSMEFTHiggs_INPUT_PARAMETERS_H

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct NMSSMEFTHiggs_input_parameters {
   double MSUSY{};
   double M1Input{};
   double M2Input{};
   double M3Input{};
   double MuInput{};
   double TanBeta{};
   double LambdaInput{};
   double KappaInput{};
   double ALambdaInput{};
   double AKappaInput{};
   Eigen::Matrix<double,3,3> mq2Input{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> mu2Input{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> md2Input{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> ml2Input{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> me2Input{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> AuInput{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> AdInput{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> AeInput{Eigen::Matrix<double,3,3>::Zero()};


   Eigen::ArrayXd get() const;
   void set(const Eigen::ArrayXd&);
};

std::ostream& operator<<(std::ostream&, const NMSSMEFTHiggs_input_parameters&);

} // namespace flexiblesusy

#endif
