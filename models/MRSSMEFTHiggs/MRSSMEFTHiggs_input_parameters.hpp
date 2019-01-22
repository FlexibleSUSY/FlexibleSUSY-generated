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

// File generated at Tue 22 Jan 2019 13:53:44

#ifndef MRSSMEFTHiggs_INPUT_PARAMETERS_H
#define MRSSMEFTHiggs_INPUT_PARAMETERS_H

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct MRSSMEFTHiggs_input_parameters {
   double TanBeta{};
   double MS{};
   double LamTDInput{};
   double LamTUInput{};
   double LamSDInput{};
   double LamSUInput{};
   double MuDInput{};
   double MuUInput{};
   double BMuInput{};
   Eigen::Matrix<double,3,3> mq2Input{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> ml2Input{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> md2Input{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> mu2Input{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> me2Input{Eigen::Matrix<double,3,3>::Zero()};
   double mS2Input{};
   double mT2Input{};
   double moc2Input{};
   double mRd2Input{};
   double mRu2Input{};
   double MDBSInput{};
   double MDWBTInput{};
   double MDGocInput{};


   Eigen::ArrayXd get() const;
   void set(const Eigen::ArrayXd&);
};

std::ostream& operator<<(std::ostream&, const MRSSMEFTHiggs_input_parameters&);

} // namespace flexiblesusy

#endif
