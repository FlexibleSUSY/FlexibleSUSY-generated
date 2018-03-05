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

// File generated at Mon 5 Mar 2018 16:09:02

#ifndef CE6SSM_INPUT_PARAMETERS_H
#define CE6SSM_INPUT_PARAMETERS_H

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct CE6SSM_input_parameters {
   double TanBeta{};
   double m0SqGuess{};
   double m12Guess{};
   double AzeroGuess{};
   double LambdaInput{};
   double KappaInput{};
   double MuPrimeInput{};
   double BMuPrimeInput{};
   double vsInput{};
   double Lambda12Input{};


   Eigen::ArrayXd get() const;
   void set(const Eigen::ArrayXd&);
};

std::ostream& operator<<(std::ostream&, const CE6SSM_input_parameters&);

} // namespace flexiblesusy

#endif
