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


#ifndef MSSMRHN_INPUT_PARAMETERS_H
#define MSSMRHN_INPUT_PARAMETERS_H

#include <complex>
#include <iosfwd>
#include <Eigen/Core>

namespace flexiblesusy {

struct MSSMRHN_input_parameters {
   double m0{};
   double m12{};
   double TanBeta{};
   int SignMu{1};
   double Azero{};
   Eigen::Matrix<double,3,3> BMvInput{Eigen::Matrix<double,3,3>::Zero()};


   Eigen::ArrayXd get() const;
   void set(const Eigen::ArrayXd&);
};

std::ostream& operator<<(std::ostream&, const MSSMRHN_input_parameters&);

} // namespace flexiblesusy

#endif
