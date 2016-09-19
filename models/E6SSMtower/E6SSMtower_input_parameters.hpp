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

// File generated at Mon 19 Sep 2016 09:41:17

#ifndef E6SSMtower_INPUT_PARAMETERS_H
#define E6SSMtower_INPUT_PARAMETERS_H

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct E6SSMtower_input_parameters {
   double MSUSY;
   double M1Input;
   double M2Input;
   double M3Input;
   double MuInput;
   double mAInput;
   double TanBeta;
   double LambdaInput;
   double gNInput;
   double M4Input;
   double mHp2Input;
   double mHpbar2Input;
   double MuPrInput;
   double BMuPrInput;
   Eigen::Matrix<double,3,3> mq2Input;
   Eigen::Matrix<double,3,3> mu2Input;
   Eigen::Matrix<double,3,3> md2Input;
   Eigen::Matrix<double,3,3> ml2Input;
   Eigen::Matrix<double,3,3> me2Input;
   Eigen::Matrix<double,3,3> AuInput;
   Eigen::Matrix<double,3,3> AdInput;
   Eigen::Matrix<double,3,3> AeInput;
   Eigen::Matrix<double,2,2> Lambda12Input;
   Eigen::Matrix<double,2,2> ALambda12Input;
   Eigen::Matrix<double,3,3> KappaInput;
   Eigen::Matrix<double,3,3> AKappaInput;
   Eigen::Matrix<double,3,3> mDx2Input;
   Eigen::Matrix<double,3,3> mDxbar2Input;
   Eigen::Matrix<double,2,2> mH1I2Input;
   Eigen::Matrix<double,2,2> mH2I2Input;
   Eigen::Matrix<double,2,2> msI2Input;

   E6SSMtower_input_parameters()
      : MSUSY(0), M1Input(0), M2Input(0), M3Input(0), MuInput(0), mAInput(0),
   TanBeta(0), LambdaInput(0), gNInput(0), M4Input(0), mHp2Input(0),
   mHpbar2Input(0), MuPrInput(0), BMuPrInput(0), mq2Input(Eigen::Matrix<double,
   3,3>::Zero()), mu2Input(Eigen::Matrix<double,3,3>::Zero()), md2Input(
   Eigen::Matrix<double,3,3>::Zero()), ml2Input(Eigen::Matrix<double,3,3>::Zero
   ()), me2Input(Eigen::Matrix<double,3,3>::Zero()), AuInput(Eigen::Matrix<
   double,3,3>::Zero()), AdInput(Eigen::Matrix<double,3,3>::Zero()), AeInput(
   Eigen::Matrix<double,3,3>::Zero()), Lambda12Input(Eigen::Matrix<double,2,2>
   ::Zero()), ALambda12Input(Eigen::Matrix<double,2,2>::Zero()), KappaInput(
   Eigen::Matrix<double,3,3>::Zero()), AKappaInput(Eigen::Matrix<double,3,3>
   ::Zero()), mDx2Input(Eigen::Matrix<double,3,3>::Zero()), mDxbar2Input(
   Eigen::Matrix<double,3,3>::Zero()), mH1I2Input(Eigen::Matrix<double,2,2>
   ::Zero()), mH2I2Input(Eigen::Matrix<double,2,2>::Zero()), msI2Input(
   Eigen::Matrix<double,2,2>::Zero())

   {}

   Eigen::ArrayXd get() const;
   void set(const Eigen::ArrayXd&);
};

std::ostream& operator<<(std::ostream&, const E6SSMtower_input_parameters&);

} // namespace flexiblesusy

#endif
