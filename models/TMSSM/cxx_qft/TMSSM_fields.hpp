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
 * @file cxx_qft/TMSSM_fields.hpp
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#ifndef TMSSM_CXXQFT_FIELDS_H
#define TMSSM_CXXQFT_FIELDS_H

#include "cxx_qft/fields.hpp"

namespace flexiblesusy {
namespace TMSSM_cxx_diagrams {
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

struct Fv {
   static constexpr auto particleType = ParticleType::fermion;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = true;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 3>
   >;
   static constexpr int numberOfGenerations = 3;
   using sm_flags = boost::mpl::vector_c<bool, true, true, true>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electricCharge = 0;
   using lorentz_conjugate = typename bar<Fv>::type;
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

struct Sd {
   static constexpr auto particleType = ParticleType::scalar;
   static constexpr auto colorRep = ParticleColorRep::triplet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 6>
   >;
   static constexpr int numberOfGenerations = 6;
   using sm_flags = boost::mpl::vector_c<bool, false, false, false, false, false, false>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electricCharge = -0.3333333333333333;
   using lorentz_conjugate = typename conj<Sd>::type;
};

struct Sv {
   static constexpr auto particleType = ParticleType::scalar;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 3>
   >;
   static constexpr int numberOfGenerations = 3;
   using sm_flags = boost::mpl::vector_c<bool, false, false, false>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electricCharge = 0;
   using lorentz_conjugate = typename conj<Sv>::type;
};

struct Su {
   static constexpr auto particleType = ParticleType::scalar;
   static constexpr auto colorRep = ParticleColorRep::triplet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 6>
   >;
   static constexpr int numberOfGenerations = 6;
   using sm_flags = boost::mpl::vector_c<bool, false, false, false, false, false, false>;
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
     boost::mpl::vector_c<int, 6>
   >;
   static constexpr int numberOfGenerations = 6;
   using sm_flags = boost::mpl::vector_c<bool, false, false, false, false, false, false>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electricCharge = -1;
   using lorentz_conjugate = typename conj<Se>::type;
};

struct hh {
   static constexpr auto particleType = ParticleType::scalar;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 3>
   >;
   static constexpr int numberOfGenerations = 3;
   using sm_flags = boost::mpl::vector_c<bool, false, false, false>;
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
     boost::mpl::vector_c<int, 3>
   >;
   static constexpr int numberOfGenerations = 3;
   using sm_flags = boost::mpl::vector_c<bool, true, false, false>;
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
     boost::mpl::vector_c<int, 4>
   >;
   static constexpr int numberOfGenerations = 4;
   using sm_flags = boost::mpl::vector_c<bool, true, false, false, false>;
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
     boost::mpl::vector_c<int, 5>
   >;
   static constexpr int numberOfGenerations = 5;
   using sm_flags = boost::mpl::vector_c<bool, false, false, false, false, false>;
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
     boost::mpl::vector_c<int, 3>
   >;
   static constexpr int numberOfGenerations = 3;
   using sm_flags = boost::mpl::vector_c<bool, false, false, false>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electricCharge = -1;
   using lorentz_conjugate = typename bar<Cha>::type;
};

struct Fe {
   static constexpr auto particleType = ParticleType::fermion;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 3>
   >;
   static constexpr int numberOfGenerations = 3;
   using sm_flags = boost::mpl::vector_c<bool, true, true, true>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electricCharge = -1;
   using lorentz_conjugate = typename bar<Fe>::type;
};

struct Fd {
   static constexpr auto particleType = ParticleType::fermion;
   static constexpr auto colorRep = ParticleColorRep::triplet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 3>
   >;
   static constexpr int numberOfGenerations = 3;
   using sm_flags = boost::mpl::vector_c<bool, true, true, true>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electricCharge = -0.3333333333333333;
   using lorentz_conjugate = typename bar<Fd>::type;
};

struct Fu {
   static constexpr auto particleType = ParticleType::fermion;
   static constexpr auto colorRep = ParticleColorRep::triplet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 3>
   >;
   static constexpr int numberOfGenerations = 3;
   using sm_flags = boost::mpl::vector_c<bool, true, true, true>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electricCharge = 0.6666666666666666;
   using lorentz_conjugate = typename bar<Fu>::type;
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

using scalars = boost::mpl::vector<Sd, Sv, Su, Se, hh, Ah, Hpm>;
using fermions = boost::mpl::vector<Glu, Fv, Chi, Cha, Fe, Fd, Fu>;
using vectors = boost::mpl::vector<VG, VP, VZ, VWm>;
using ghosts = boost::mpl::vector<gG, gP, gZ, gWm, gWmC>;

} // namespace fields

} // namespace TMSSM_cxx_diagrams

namespace cxx_diagrams {
namespace fields{

// Fields that are their own Lorentz conjugates.
template<> struct conj<flexiblesusy::TMSSM_cxx_diagrams::fields::VG> { using type = flexiblesusy::TMSSM_cxx_diagrams::fields::VG; };
template<> struct bar<flexiblesusy::TMSSM_cxx_diagrams::fields::Glu> { using type = flexiblesusy::TMSSM_cxx_diagrams::fields::Glu; };
template<> struct conj<flexiblesusy::TMSSM_cxx_diagrams::fields::VP> { using type = flexiblesusy::TMSSM_cxx_diagrams::fields::VP; };
template<> struct conj<flexiblesusy::TMSSM_cxx_diagrams::fields::VZ> { using type = flexiblesusy::TMSSM_cxx_diagrams::fields::VZ; };
template<> struct conj<flexiblesusy::TMSSM_cxx_diagrams::fields::hh> { using type = flexiblesusy::TMSSM_cxx_diagrams::fields::hh; };
template<> struct conj<flexiblesusy::TMSSM_cxx_diagrams::fields::Ah> { using type = flexiblesusy::TMSSM_cxx_diagrams::fields::Ah; };
template<> struct bar<flexiblesusy::TMSSM_cxx_diagrams::fields::Chi> { using type = flexiblesusy::TMSSM_cxx_diagrams::fields::Chi; };



template<>
struct is_vector<flexiblesusy::TMSSM_cxx_diagrams::fields::VG > : public std::true_type {};
template<>
struct is_ghost<flexiblesusy::TMSSM_cxx_diagrams::fields::gG > : public std::true_type {};
template<>
struct is_fermion<flexiblesusy::TMSSM_cxx_diagrams::fields::Glu > : public std::true_type {};
template<>
struct is_fermion<flexiblesusy::TMSSM_cxx_diagrams::fields::Fv > : public std::true_type {};
template<>
struct is_vector<flexiblesusy::TMSSM_cxx_diagrams::fields::VP > : public std::true_type {};
template<>
struct is_vector<flexiblesusy::TMSSM_cxx_diagrams::fields::VZ > : public std::true_type {};
template<>
struct is_ghost<flexiblesusy::TMSSM_cxx_diagrams::fields::gP > : public std::true_type {};
template<>
struct is_ghost<flexiblesusy::TMSSM_cxx_diagrams::fields::gZ > : public std::true_type {};
template<>
struct is_ghost<flexiblesusy::TMSSM_cxx_diagrams::fields::gWm > : public std::true_type {};
template<>
struct is_ghost<flexiblesusy::TMSSM_cxx_diagrams::fields::gWmC > : public std::true_type {};
template<>
struct is_scalar<flexiblesusy::TMSSM_cxx_diagrams::fields::Sd > : public std::true_type {};
template<>
struct is_scalar<flexiblesusy::TMSSM_cxx_diagrams::fields::Sv > : public std::true_type {};
template<>
struct is_scalar<flexiblesusy::TMSSM_cxx_diagrams::fields::Su > : public std::true_type {};
template<>
struct is_scalar<flexiblesusy::TMSSM_cxx_diagrams::fields::Se > : public std::true_type {};
template<>
struct is_scalar<flexiblesusy::TMSSM_cxx_diagrams::fields::hh > : public std::true_type {};
template<>
struct is_scalar<flexiblesusy::TMSSM_cxx_diagrams::fields::Ah > : public std::true_type {};
template<>
struct is_scalar<flexiblesusy::TMSSM_cxx_diagrams::fields::Hpm > : public std::true_type {};
template<>
struct is_fermion<flexiblesusy::TMSSM_cxx_diagrams::fields::Chi > : public std::true_type {};
template<>
struct is_fermion<flexiblesusy::TMSSM_cxx_diagrams::fields::Cha > : public std::true_type {};
template<>
struct is_fermion<flexiblesusy::TMSSM_cxx_diagrams::fields::Fe > : public std::true_type {};
template<>
struct is_fermion<flexiblesusy::TMSSM_cxx_diagrams::fields::Fd > : public std::true_type {};
template<>
struct is_fermion<flexiblesusy::TMSSM_cxx_diagrams::fields::Fu > : public std::true_type {};
template<>
struct is_vector<flexiblesusy::TMSSM_cxx_diagrams::fields::VWm > : public std::true_type {};


} // namespace fields
} // namespace cxx_diagrams

} // namespace flexiblesusy

#endif
