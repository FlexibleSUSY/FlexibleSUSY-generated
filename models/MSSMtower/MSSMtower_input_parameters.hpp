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

// File generated at Wed 12 Apr 2017 10:30:06

#ifndef MSSMtower_INPUT_PARAMETERS_H
#define MSSMtower_INPUT_PARAMETERS_H

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct MSSMtower_input_parameters {
   double MSUSY;
   double M1Input;
   double M2Input;
   double M3Input;
   double MuInput;
   double mAInput;
   double TanBeta;
   Eigen::Matrix<double,3,3> mq2Input;
   Eigen::Matrix<double,3,3> mu2Input;
   Eigen::Matrix<double,3,3> md2Input;
   Eigen::Matrix<double,3,3> ml2Input;
   Eigen::Matrix<double,3,3> me2Input;
   Eigen::Matrix<double,3,3> AuInput;
   Eigen::Matrix<double,3,3> AdInput;
   Eigen::Matrix<double,3,3> AeInput;

   MSSMtower_input_parameters()
      : MSUSY(0), M1Input(0), M2Input(0), M3Input(0), MuInput(0), mAInput(0),
   TanBeta(0), mq2Input(Eigen::Matrix<double,3,3>::Zero()), mu2Input(
   Eigen::Matrix<double,3,3>::Zero()), md2Input(Eigen::Matrix<double,3,3>::Zero
   ()), ml2Input(Eigen::Matrix<double,3,3>::Zero()), me2Input(Eigen::Matrix<
   double,3,3>::Zero()), AuInput(Eigen::Matrix<double,3,3>::Zero()), AdInput(
   Eigen::Matrix<double,3,3>::Zero()), AeInput(Eigen::Matrix<double,3,3>::Zero(
   ))

   {}

   Eigen::ArrayXd get() const;
   void set(const Eigen::ArrayXd&);
};

std::ostream& operator<<(std::ostream&, const MSSMtower_input_parameters&);

} // namespace flexiblesusy

#endif
