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
 * @file cxx_qft/MSSMNoFVatMGUT_fields.hpp
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#ifndef MSSMNoFVatMGUT_CXXQFT_FIELDS_H
#define MSSMNoFVatMGUT_CXXQFT_FIELDS_H

#include "cxx_qft/fields.hpp"

namespace flexiblesusy {
namespace MSSMNoFVatMGUT_cxx_diagrams {
namespace fields {

using cxx_diagrams::fields::ParticleType;
using cxx_diagrams::fields::ParticleColorRep;
using cxx_diagrams::fields::conj;
using cxx_diagrams::fields::bar;

struct VG {
   static constexpr auto particleType = ParticleType::vector;
   static constexpr auto colorRep = ParticleColorRep::octet;
   static constexpr auto massless = true;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 0;
   using lorentz_conjugate = VG;
};

struct gG {
   static constexpr auto particleType = ParticleType::ghost;
   static constexpr auto colorRep = ParticleColorRep::octet;
   static constexpr auto massless = true;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 0;
   using lorentz_conjugate = typename bar<gG>::type;
};

struct Glu {
   static constexpr auto particleType = ParticleType::fermion;
   static constexpr auto colorRep = ParticleColorRep::octet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, false>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 0;
   using lorentz_conjugate = Glu;
};

struct VP {
   static constexpr auto particleType = ParticleType::vector;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = true;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 0;
   using lorentz_conjugate = VP;
};

struct VZ {
   static constexpr auto particleType = ParticleType::vector;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 0;
   using lorentz_conjugate = VZ;
};

struct gP {
   static constexpr auto particleType = ParticleType::ghost;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = true;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 0;
   using lorentz_conjugate = typename bar<gP>::type;
};

struct gZ {
   static constexpr auto particleType = ParticleType::ghost;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 0;
   using lorentz_conjugate = typename bar<gZ>::type;
};

struct gWm {
   static constexpr auto particleType = ParticleType::ghost;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = -1;
   using lorentz_conjugate = typename bar<gWm>::type;
};

struct gWmC {
   static constexpr auto particleType = ParticleType::ghost;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 1;
   using lorentz_conjugate = typename bar<gWmC>::type;
};

struct Fd {
   static constexpr auto particleType = ParticleType::fermion;
   static constexpr auto colorRep = ParticleColorRep::triplet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = -0.3333333333333333;
   using lorentz_conjugate = typename bar<Fd>::type;
};

struct Fs {
   static constexpr auto particleType = ParticleType::fermion;
   static constexpr auto colorRep = ParticleColorRep::triplet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = -0.3333333333333333;
   using lorentz_conjugate = typename bar<Fs>::type;
};

struct Fb {
   static constexpr auto particleType = ParticleType::fermion;
   static constexpr auto colorRep = ParticleColorRep::triplet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = -0.3333333333333333;
   using lorentz_conjugate = typename bar<Fb>::type;
};

struct Fu {
   static constexpr auto particleType = ParticleType::fermion;
   static constexpr auto colorRep = ParticleColorRep::triplet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 0.6666666666666666;
   using lorentz_conjugate = typename bar<Fu>::type;
};

struct Fc {
   static constexpr auto particleType = ParticleType::fermion;
   static constexpr auto colorRep = ParticleColorRep::triplet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 0.6666666666666666;
   using lorentz_conjugate = typename bar<Fc>::type;
};

struct Ft {
   static constexpr auto particleType = ParticleType::fermion;
   static constexpr auto colorRep = ParticleColorRep::triplet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 0.6666666666666666;
   using lorentz_conjugate = typename bar<Ft>::type;
};

struct Fve {
   static constexpr auto particleType = ParticleType::fermion;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 0;
   using lorentz_conjugate = typename bar<Fve>::type;
};

struct Fvm {
   static constexpr auto particleType = ParticleType::fermion;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 0;
   using lorentz_conjugate = typename bar<Fvm>::type;
};

struct Fvt {
   static constexpr auto particleType = ParticleType::fermion;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 0;
   using lorentz_conjugate = typename bar<Fvt>::type;
};

struct Fe {
   static constexpr auto particleType = ParticleType::fermion;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = -1;
   using lorentz_conjugate = typename bar<Fe>::type;
};

struct Fm {
   static constexpr auto particleType = ParticleType::fermion;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = -1;
   using lorentz_conjugate = typename bar<Fm>::type;
};

struct Ftau {
   static constexpr auto particleType = ParticleType::fermion;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = -1;
   using lorentz_conjugate = typename bar<Ftau>::type;
};

struct SveL {
   static constexpr auto particleType = ParticleType::scalar;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, false>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 0;
   using lorentz_conjugate = typename conj<SveL>::type;
};

struct SvmL {
   static constexpr auto particleType = ParticleType::scalar;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, false>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 0;
   using lorentz_conjugate = typename conj<SvmL>::type;
};

struct SvtL {
   static constexpr auto particleType = ParticleType::scalar;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, false>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 0;
   using lorentz_conjugate = typename conj<SvtL>::type;
};

struct Sd {
   static constexpr auto particleType = ParticleType::scalar;
   static constexpr auto colorRep = ParticleColorRep::triplet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 2>
   >;
   static constexpr int numberOfGenerations = 2;
   using sm_flags = boost::mpl::vector_c<bool, false, false>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electricCharge = -0.3333333333333333;
   using lorentz_conjugate = typename conj<Sd>::type;
};

struct Su {
   static constexpr auto particleType = ParticleType::scalar;
   static constexpr auto colorRep = ParticleColorRep::triplet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 2>
   >;
   static constexpr int numberOfGenerations = 2;
   using sm_flags = boost::mpl::vector_c<bool, false, false>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electricCharge = 0.6666666666666666;
   using lorentz_conjugate = typename conj<Su>::type;
};

struct Se {
   static constexpr auto particleType = ParticleType::scalar;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 2>
   >;
   static constexpr int numberOfGenerations = 2;
   using sm_flags = boost::mpl::vector_c<bool, false, false>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electricCharge = -1;
   using lorentz_conjugate = typename conj<Se>::type;
};

struct Sm {
   static constexpr auto particleType = ParticleType::scalar;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 2>
   >;
   static constexpr int numberOfGenerations = 2;
   using sm_flags = boost::mpl::vector_c<bool, false, false>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electricCharge = -1;
   using lorentz_conjugate = typename conj<Sm>::type;
};

struct Stau {
   static constexpr auto particleType = ParticleType::scalar;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 2>
   >;
   static constexpr int numberOfGenerations = 2;
   using sm_flags = boost::mpl::vector_c<bool, false, false>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electricCharge = -1;
   using lorentz_conjugate = typename conj<Stau>::type;
};

struct Ss {
   static constexpr auto particleType = ParticleType::scalar;
   static constexpr auto colorRep = ParticleColorRep::triplet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 2>
   >;
   static constexpr int numberOfGenerations = 2;
   using sm_flags = boost::mpl::vector_c<bool, false, false>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electricCharge = -0.3333333333333333;
   using lorentz_conjugate = typename conj<Ss>::type;
};

struct Sc {
   static constexpr auto particleType = ParticleType::scalar;
   static constexpr auto colorRep = ParticleColorRep::triplet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 2>
   >;
   static constexpr int numberOfGenerations = 2;
   using sm_flags = boost::mpl::vector_c<bool, false, false>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electricCharge = 0.6666666666666666;
   using lorentz_conjugate = typename conj<Sc>::type;
};

struct Sb {
   static constexpr auto particleType = ParticleType::scalar;
   static constexpr auto colorRep = ParticleColorRep::triplet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 2>
   >;
   static constexpr int numberOfGenerations = 2;
   using sm_flags = boost::mpl::vector_c<bool, false, false>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electricCharge = -0.3333333333333333;
   using lorentz_conjugate = typename conj<Sb>::type;
};

struct St {
   static constexpr auto particleType = ParticleType::scalar;
   static constexpr auto colorRep = ParticleColorRep::triplet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 2>
   >;
   static constexpr int numberOfGenerations = 2;
   using sm_flags = boost::mpl::vector_c<bool, false, false>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electricCharge = 0.6666666666666666;
   using lorentz_conjugate = typename conj<St>::type;
};

struct hh {
   static constexpr auto particleType = ParticleType::scalar;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 2>
   >;
   static constexpr int numberOfGenerations = 2;
   using sm_flags = boost::mpl::vector_c<bool, false, false>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electricCharge = 0;
   using lorentz_conjugate = hh;
};

struct Ah {
   static constexpr auto particleType = ParticleType::scalar;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 2>
   >;
   static constexpr int numberOfGenerations = 2;
   using sm_flags = boost::mpl::vector_c<bool, true, false>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electricCharge = 0;
   using lorentz_conjugate = Ah;
};

struct Hpm {
   static constexpr auto particleType = ParticleType::scalar;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 2>
   >;
   static constexpr int numberOfGenerations = 2;
   using sm_flags = boost::mpl::vector_c<bool, true, false>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electricCharge = -1;
   using lorentz_conjugate = typename conj<Hpm>::type;
};

struct Chi {
   static constexpr auto particleType = ParticleType::fermion;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 4>
   >;
   static constexpr int numberOfGenerations = 4;
   using sm_flags = boost::mpl::vector_c<bool, false, false, false, false>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electricCharge = 0;
   using lorentz_conjugate = Chi;
};

struct Cha {
   static constexpr auto particleType = ParticleType::fermion;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 2>
   >;
   static constexpr int numberOfGenerations = 2;
   using sm_flags = boost::mpl::vector_c<bool, false, false>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electricCharge = -1;
   using lorentz_conjugate = typename bar<Cha>::type;
};

struct VWm {
   static constexpr auto particleType = ParticleType::vector;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = -1;
   using lorentz_conjugate = typename conj<VWm>::type;
};

// Named fields
using Electron = Fe;

using scalars = boost::mpl::vector<SveL, SvmL, SvtL, Sd, Su, Se, Sm, Stau, Ss, Sc, Sb, St, hh, Ah, Hpm>;
using fermions = boost::mpl::vector<Glu, Fd, Fs, Fb, Fu, Fc, Ft, Fve, Fvm, Fvt, Fe, Fm, Ftau, Chi, Cha>;
using vectors = boost::mpl::vector<VG, VP, VZ, VWm>;
using ghosts = boost::mpl::vector<gG, gP, gZ, gWm, gWmC>;

} // namespace fields

} // namespace MSSMNoFVatMGUT_cxx_diagrams

namespace cxx_diagrams {
namespace fields{

// Fields that are their own Lorentz conjugates.
template<> struct conj<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::VG> { using type = flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::VG; };
template<> struct bar<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Glu> { using type = flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Glu; };
template<> struct conj<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::VP> { using type = flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::VP; };
template<> struct conj<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::VZ> { using type = flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::VZ; };
template<> struct conj<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::hh> { using type = flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::hh; };
template<> struct conj<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Ah> { using type = flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Ah; };
template<> struct bar<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Chi> { using type = flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Chi; };



template<>
struct is_vector<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::VG > : public std::true_type {};
template<>
struct is_ghost<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::gG > : public std::true_type {};
template<>
struct is_fermion<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Glu > : public std::true_type {};
template<>
struct is_vector<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::VP > : public std::true_type {};
template<>
struct is_vector<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::VZ > : public std::true_type {};
template<>
struct is_ghost<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::gP > : public std::true_type {};
template<>
struct is_ghost<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::gZ > : public std::true_type {};
template<>
struct is_ghost<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::gWm > : public std::true_type {};
template<>
struct is_ghost<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::gWmC > : public std::true_type {};
template<>
struct is_fermion<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Fd > : public std::true_type {};
template<>
struct is_fermion<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Fs > : public std::true_type {};
template<>
struct is_fermion<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Fb > : public std::true_type {};
template<>
struct is_fermion<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Fu > : public std::true_type {};
template<>
struct is_fermion<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Fc > : public std::true_type {};
template<>
struct is_fermion<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Ft > : public std::true_type {};
template<>
struct is_fermion<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Fve > : public std::true_type {};
template<>
struct is_fermion<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Fvm > : public std::true_type {};
template<>
struct is_fermion<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Fvt > : public std::true_type {};
template<>
struct is_fermion<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Fe > : public std::true_type {};
template<>
struct is_fermion<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Fm > : public std::true_type {};
template<>
struct is_fermion<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Ftau > : public std::true_type {};
template<>
struct is_scalar<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::SveL > : public std::true_type {};
template<>
struct is_scalar<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::SvmL > : public std::true_type {};
template<>
struct is_scalar<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::SvtL > : public std::true_type {};
template<>
struct is_scalar<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Sd > : public std::true_type {};
template<>
struct is_scalar<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Su > : public std::true_type {};
template<>
struct is_scalar<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Se > : public std::true_type {};
template<>
struct is_scalar<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Sm > : public std::true_type {};
template<>
struct is_scalar<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Stau > : public std::true_type {};
template<>
struct is_scalar<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Ss > : public std::true_type {};
template<>
struct is_scalar<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Sc > : public std::true_type {};
template<>
struct is_scalar<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Sb > : public std::true_type {};
template<>
struct is_scalar<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::St > : public std::true_type {};
template<>
struct is_scalar<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::hh > : public std::true_type {};
template<>
struct is_scalar<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Ah > : public std::true_type {};
template<>
struct is_scalar<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Hpm > : public std::true_type {};
template<>
struct is_fermion<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Chi > : public std::true_type {};
template<>
struct is_fermion<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::Cha > : public std::true_type {};
template<>
struct is_vector<flexiblesusy::MSSMNoFVatMGUT_cxx_diagrams::fields::VWm > : public std::true_type {};


} // namespace fields
} // namespace cxx_diagrams

} // namespace flexiblesusy

#endif
