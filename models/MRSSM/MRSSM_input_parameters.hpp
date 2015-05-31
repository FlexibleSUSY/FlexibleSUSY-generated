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

// File generated at Sun 31 May 2015 12:24:58

#ifndef MRSSM_INPUT_PARAMETERS_H
#define MRSSM_INPUT_PARAMETERS_H

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct MRSSM_input_parameters {
   double TanBeta;
   double LamTDInput;
   double LamTUInput;
   double LamSDInput;
   double LamSUInput;
   double MuInput;
   double MuDInput;
   double MuUInput;
   double vTInput;
   double vSInput;
   double BMuInput;
   double BMuDInput;
   double BMuUInput;
   Eigen::Matrix<double,3,3> mq2Input;
   Eigen::Matrix<double,3,3> ml2Input;
   Eigen::Matrix<double,3,3> md2Input;
   Eigen::Matrix<double,3,3> mu2Input;
   Eigen::Matrix<double,3,3> me2Input;
   double moc2Input;
   double mRd2Input;
   double mRu2Input;
   double MDBSInput;
   double MDWBTInput;
   double MDGocInput;

   MRSSM_input_parameters()
      : TanBeta(0), LamTDInput(0), LamTUInput(0), LamSDInput(0), LamSUInput(0),
   MuInput(0), MuDInput(0), MuUInput(0), vTInput(0), vSInput(0), BMuInput(0),
   BMuDInput(0), BMuUInput(0), mq2Input(Eigen::Matrix<double,3,3>::Zero()),
   ml2Input(Eigen::Matrix<double,3,3>::Zero()), md2Input(Eigen::Matrix<double,3
   ,3>::Zero()), mu2Input(Eigen::Matrix<double,3,3>::Zero()), me2Input(
   Eigen::Matrix<double,3,3>::Zero()), moc2Input(0), mRd2Input(0), mRu2Input(0)
   , MDBSInput(0), MDWBTInput(0), MDGocInput(0)

   {}
};

std::ostream& operator<<(std::ostream&, const MRSSM_input_parameters&);

} // namespace flexiblesusy

#endif
