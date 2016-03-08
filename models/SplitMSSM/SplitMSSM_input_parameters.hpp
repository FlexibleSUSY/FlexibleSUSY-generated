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

// File generated at Tue 8 Mar 2016 16:04:18

#ifndef SplitMSSM_INPUT_PARAMETERS_H
#define SplitMSSM_INPUT_PARAMETERS_H

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct SplitMSSM_input_parameters {
   double MSUSY;
   double M1Input;
   double M2Input;
   double M3Input;
   double MuInput;
   double mAInput;
   double MEWSB;
   double AtInput;
   double TanBeta;
   Eigen::Matrix<double,3,3> msq2;
   Eigen::Matrix<double,3,3> msu2;
   Eigen::Matrix<double,3,3> msd2;
   Eigen::Matrix<double,3,3> msl2;
   Eigen::Matrix<double,3,3> mse2;

   SplitMSSM_input_parameters()
      : MSUSY(0), M1Input(0), M2Input(0), M3Input(0), MuInput(0), mAInput(0),
   MEWSB(0), AtInput(0), TanBeta(0), msq2(Eigen::Matrix<double,3,3>::Zero()),
   msu2(Eigen::Matrix<double,3,3>::Zero()), msd2(Eigen::Matrix<double,3,3>
   ::Zero()), msl2(Eigen::Matrix<double,3,3>::Zero()), mse2(Eigen::Matrix<
   double,3,3>::Zero())

   {}

   Eigen::ArrayXd get() const;
   void set(const Eigen::ArrayXd&);
};

std::ostream& operator<<(std::ostream&, const SplitMSSM_input_parameters&);

} // namespace flexiblesusy

#endif
