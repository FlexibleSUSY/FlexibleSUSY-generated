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

// File generated at Sun 4 Aug 2019 19:44:56

#ifndef lowNMSSMTanBetaAtMZ_INPUT_PARAMETERS_H
#define lowNMSSMTanBetaAtMZ_INPUT_PARAMETERS_H

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct lowNMSSMTanBetaAtMZ_input_parameters {
   double TanBeta{};
   double Qin{};
   double M1Input{};
   double M2Input{};
   double M3Input{};
   double AtInput{};
   double AbInput{};
   double ATauInput{};
   double ml1Input{};
   double ml2Input{};
   double ml3Input{};
   double me1Input{};
   double me2Input{};
   double me3Input{};
   double mq1Input{};
   double mq2Input{};
   double mq3Input{};
   double md1Input{};
   double md2Input{};
   double md3Input{};
   double mu1Input{};
   double mu2Input{};
   double mu3Input{};
   double LambdaInput{};
   double KappaInput{};
   double ALambdaInput{};
   double AKappaInput{};
   double MuEffInput{};


   Eigen::ArrayXd get() const;
   void set(const Eigen::ArrayXd&);
};

std::ostream& operator<<(std::ostream&, const lowNMSSMTanBetaAtMZ_input_parameters&);

} // namespace flexiblesusy

#endif
