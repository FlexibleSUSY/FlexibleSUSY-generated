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
 * @file MRSSM2_unitarity.cpp
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#include "MRSSM2_unitarity.hpp"
#include "MRSSM2_mass_eigenstates.hpp"
#include "cxx_qft/MRSSM2_fields.hpp"

#include "sum.hpp"
#include "wrappers.hpp"

#include <Eigen/Eigenvalues>

namespace flexiblesusy {
namespace MRSSM2_unitarity {

static constexpr int size = 0;

namespace {

inline
double Sqrt2(int i, int j) {
   // 1/Sqrt[2] or 1
   return i==j ? oneOverSqrt2 : 1;
}

inline
bool is_unitarity_fulfilled(double eigenval) {
   // best eigenvalue <= 1/2 then the point is allowed
   return eigenval <= 0.5;
}

inline
double best_eigenvalue(Eigen::MatrixXcd const& m) {
   // max(|re(eigenvalues(a0))|)
   return m.eigenvalues().unaryExpr(
         [](std::complex<double> const& el) {return std::real(el);}
      ).cwiseAbs().maxCoeff();
}

}

// s -> infinity limit
UnitarityInfiniteS
max_scattering_eigenvalue_infinite_s(MRSSM2_mass_eigenstates const& model) {
   Eigen::MatrixXcd matrix = Eigen::MatrixXcd::Zero(size, size);

   using namespace MRSSM2_cxx_diagrams::fields;


   // reinstate factor 1/(16*Pi) that was removed from temp values
   matrix = matrix.unaryExpr([](std::complex<double> const& el) { return oneOver16Pi*el; });

   const double best_eigenval = best_eigenvalue(matrix);
   return {is_unitarity_fulfilled(best_eigenval), model.get_scale(), best_eigenval, matrix};
}

} // namespace MRSSM2_unitarity
} // namespace flexiblesusy
