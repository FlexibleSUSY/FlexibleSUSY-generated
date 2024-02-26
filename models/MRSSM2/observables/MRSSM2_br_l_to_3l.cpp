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
 * @file MRSSM2_br_l_to_3l.cpp
 *
 * This file was generated at Sun 25 Feb 2024 20:35:27 with FlexibleSUSY
 * 2.8.0 and SARAH 4.15.1
 */

#include <valarray>
#include <complex>

#include "MRSSM2_mass_eigenstates.hpp"
#include "cxx_qft/MRSSM2_qft.hpp"

#include "MRSSM2_br_l_to_3l.hpp"
#include "MRSSM2_FFV_form_factors.hpp"

#define NPF(name) fix_tensors_sign<name>(model, nI, nO, nA);
#define F4F(lepton, g) embed_photon<lepton, Photon>(model, form_factors, g);

#include "loop_libraries/loop_library.hpp"
#include "cxx_qft/MRSSM2_npointfunctions_wilsoncoeffs.hpp"
#include "concatenate.hpp"
#include <limits>
#include <type_traits>
#include <boost/fusion/include/at_key.hpp>
#include "wrappers.hpp"

namespace flexiblesusy {

namespace MRSSM2_cxx_diagrams {
namespace npointfunctions {


class nPointFeFeFeFe_3917882094413845 : public correlation_function_context<4,0>
{
   using generic_sum_base = correlation_function_context<4,0>;

   template<class GenericFieldMap>
   struct subexpression_base :
   generic_sum_base, index_map_interface<GenericFieldMap>
   {
      subexpression_base( const subexpression_base & ) = default;

      subexpression_base( const generic_sum_base &gsb,
         const typename field_index_map<GenericFieldMap>::type &fim ) :
      generic_sum_base( gsb ), index_map_interface<GenericFieldMap>( fim )
      {}
   }; // End of subexpression_base<GenericFieldMap>

   struct GenericF7Key {};
   struct GenericF8Key {};
   struct GenericS5Key {};
   struct GenericS6Key {};
   struct GenericF6Key {};
   struct GenericS7Key {};
   struct GenericS8Key {};

   template<class GenericFieldMap>
   struct genericSum1_impl : generic_sum_base {
      genericSum1_impl( const generic_sum_base &base ) :
      generic_sum_base( base ) {
      } // End of constructor genericSum1_impl
   
      std::array<std::complex<double>,10> operator()( void ) {
         using boost::mpl::at;
         using boost::fusion::at_key;
         using gF7 = typename at<GenericFieldMap,GenericF7Key>::type;
         using gF8 = typename at<GenericFieldMap,GenericF8Key>::type;
         using gS5 = typename at<GenericFieldMap,GenericS5Key>::type;
         using gS6 = typename at<GenericFieldMap,GenericS6Key>::type;
   
         std::array<int, 1> i1 {this->external_indices(0)};
         std::array<int, 1> i2 {this->external_indices(1)};
         std::array<int, 1> i3 {this->external_indices(2)};
         std::array<int, 1> i4 {this->external_indices(3)};
         
         typename field_index_map<GenericFieldMap>::type index_map;
         const context_with_vertices &context = *this;
         std::complex<double> S_LL = 0.0;
         std::complex<double> S_LR = 0.0;
         std::complex<double> S_RL = 0.0;
         std::complex<double> S_RR = 0.0;
         std::complex<double> V_LL = 0.0;
         std::complex<double> V_LR = 0.0;
         std::complex<double> V_RL = 0.0;
         std::complex<double> V_RR = 0.0;
         std::complex<double> T_LL = 0.0;
         std::complex<double> T_RR = 0.0;
   
         // Shorter aliases for large types
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _gF7 = typename cxx_diagrams::bar<gF7>::type;
         using _gF8 = typename cxx_diagrams::bar<gF8>::type;
         using _gS5 = typename cxx_diagrams::conj<gS5>::type;
         using _gS6 = typename cxx_diagrams::conj<gS6>::type;
         using Fe = fields::Fe;
         using Fe = fields::Fe;
         
         // Aliases for loop function identifiers
         using looplibrary::bb0;
         using looplibrary::cc0;
         
         double mF7, mF8, mS5, mS6;
         std::complex<double> g1, g2, g3, g4, g5, g6, g7, g8;
         std::array<int, 8> z{};
         looplibrary::Bcoeff_t l1 {};
         looplibrary::Ccoeff_t l2 {};
         // Start of summation over generic fields.
         for( const auto &iF7 : index_range<gF7>() ) {
         at_key<GenericF7Key>( index_map ) = iF7;
         for( const auto &iF8 : index_range<gF8>() ) {
         at_key<GenericF8Key>( index_map ) = iF8;
         for( const auto &iS5 : index_range<gS5>() ) {
         at_key<GenericS5Key>( index_map ) = iS5;
         for( const auto &iS6 : index_range<gS6>() ) {
         at_key<GenericS6Key>( index_map ) = iS6;
         
            mF7 = context.mass<gF7>(iF7);
            mF8 = context.mass<gF8>(iF8);
            mS5 = context.mass<gS5>(iS5);
            mS6 = context.mass<gS6>(iS6);
            
            g1 = NPF_L(_Fe, _gF8, gS6) NPF_I(i3, iF8, iS6);
            g2 = NPF_R(_Fe, _gF8, gS6) NPF_I(i3, iF8, iS6);
            g3 = NPF_L(_Fe, Fe, _gS5) NPF_I(i4, i2, iS5);
            g4 = NPF_R(_Fe, Fe, _gS5) NPF_I(i4, i2, iS5);
            g5 = NPF_L(_gF7, Fe, _gS6) NPF_I(iF7, i1, iS6);
            g6 = NPF_R(_gF7, Fe, _gS6) NPF_I(iF7, i1, iS6);
            g7 = NPF_L(gF8, gF7, gS5) NPF_I(iF8, iF7, iS5);
            g8 = NPF_R(gF8, gF7, gS5) NPF_I(iF8, iF7, iS5);
            
            z[0] = (std::abs(g1) < std::numeric_limits<double>::epsilon())?0:1;
            z[1] = (std::abs(g2) < std::numeric_limits<double>::epsilon())?0:1;
            z[2] = (std::abs(g3) < std::numeric_limits<double>::epsilon())?0:1;
            z[3] = (std::abs(g4) < std::numeric_limits<double>::epsilon())?0:1;
            z[4] = (std::abs(g5) < std::numeric_limits<double>::epsilon())?0:1;
            z[5] = (std::abs(g6) < std::numeric_limits<double>::epsilon())?0:1;
            z[6] = (std::abs(g7) < std::numeric_limits<double>::epsilon())?0:1;
            z[7] = (std::abs(g8) < std::numeric_limits<double>::epsilon())?0:1;
            
            if( (z[2] + z[3])*(z[0]*z[4] + z[1]*z[5])*(z[6] + z[7]) == 0 ) continue;
            Loop_library::get().B(l1,0,Sqr(mF7),Sqr(mF8),Sqr(context.scale()));
            Loop_library::get().C(l2,0,0,0,Sqr(mF8),Sqr(mF7),Sqr(mS6),Sqr(context.scale()));
         
            S_LL += (0.006332573977646111*g1*g3*g5*(g7*l2[cc0]*mF7*mF8 + g8*(l1[bb0] + l2[cc0]*Sqr(mS6))))/Sqr(mS5);
            S_LR += (0.006332573977646111*g1*g4*g5*(g7*l2[cc0]*mF7*mF8 + g8*(l1[bb0] + l2[cc0]*Sqr(mS6))))/Sqr(mS5);
            S_RL += (0.006332573977646111*g2*g3*g6*(g8*l2[cc0]*mF7*mF8 + g7*(l1[bb0] + l2[cc0]*Sqr(mS6))))/Sqr(mS5);
            S_RR += (0.006332573977646111*g2*g4*g6*(g8*l2[cc0]*mF7*mF8 + g7*(l1[bb0] + l2[cc0]*Sqr(mS6))))/Sqr(mS5);
            V_LL += 0;
            V_LR += 0;
            V_RL += 0;
            V_RR += 0;
            T_LL += 0;
            T_RR += 0;
         
         }}}} // End of summation over generic fields
   
         return {S_LL, S_LR, S_RL, S_RR, V_LL, V_LR, V_RL, V_RR, T_LL, T_RR};
      } // End of operator()( void )
   }; // End of struct genericSum1_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum1( void ) {
      using GenericKeys = boost::mpl::vector< GenericF7Key, GenericF8Key, GenericS5Key, GenericS6Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Cha1>::type, fields::Cha1, fields::hh, fields::Sv>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Chi>::type, fields::Chi, fields::hh, fields::Se>,
         boost::mpl::vector<fields::Chi, typename cxx_diagrams::bar<fields::Chi>::type, fields::hh, fields::Se>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Cha1>::type, fields::Cha1, fields::Ah, fields::Sv>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Chi>::type, fields::Chi, fields::Ah, fields::Se>,
         boost::mpl::vector<fields::Chi, typename cxx_diagrams::bar<fields::Chi>::type, fields::Ah, fields::Se>
         >;
      using combinatorial_factors = boost::mpl::vector<
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>
         >;
      using color_factors = boost::mpl::vector<
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>
         >;
      return accumulate_generic<
         GenericKeys,
         GenericInsertions,
         combinatorial_factors,
         color_factors,
         boost::mpl::int_<10>,
         genericSum1_impl
         >( *this );
   } // End of function genericSum1()
   
   template<class GenericFieldMap>
   struct genericSum2_impl : generic_sum_base {
      genericSum2_impl( const generic_sum_base &base ) :
      generic_sum_base( base ) {
      } // End of constructor genericSum2_impl
   
      std::array<std::complex<double>,10> operator()( void ) {
         using boost::mpl::at;
         using boost::fusion::at_key;
         using gF6 = typename at<GenericFieldMap,GenericF6Key>::type;
         using gS5 = typename at<GenericFieldMap,GenericS5Key>::type;
         using gS7 = typename at<GenericFieldMap,GenericS7Key>::type;
         using gS8 = typename at<GenericFieldMap,GenericS8Key>::type;
   
         std::array<int, 1> i1 {this->external_indices(0)};
         std::array<int, 1> i2 {this->external_indices(1)};
         std::array<int, 1> i3 {this->external_indices(2)};
         std::array<int, 1> i4 {this->external_indices(3)};
         
         typename field_index_map<GenericFieldMap>::type index_map;
         const context_with_vertices &context = *this;
         std::complex<double> S_LL = 0.0;
         std::complex<double> S_LR = 0.0;
         std::complex<double> S_RL = 0.0;
         std::complex<double> S_RR = 0.0;
         std::complex<double> V_LL = 0.0;
         std::complex<double> V_LR = 0.0;
         std::complex<double> V_RL = 0.0;
         std::complex<double> V_RR = 0.0;
         std::complex<double> T_LL = 0.0;
         std::complex<double> T_RR = 0.0;
   
         // Shorter aliases for large types
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _gF6 = typename cxx_diagrams::bar<gF6>::type;
         using _gS5 = typename cxx_diagrams::conj<gS5>::type;
         using _gS7 = typename cxx_diagrams::conj<gS7>::type;
         using _gS8 = typename cxx_diagrams::conj<gS8>::type;
         using Fe = fields::Fe;
         using Fe = fields::Fe;
         
         // Aliases for loop function identifiers
         using looplibrary::cc0;
         
         double mF6, mS5, mS7, mS8;
         std::complex<double> g1, g2, g3, g4, g5, g6, g7;
         std::array<int, 7> z{};
         looplibrary::Ccoeff_t l1 {};
         // Start of summation over generic fields.
         for( const auto &iF6 : index_range<gF6>() ) {
         at_key<GenericF6Key>( index_map ) = iF6;
         for( const auto &iS5 : index_range<gS5>() ) {
         at_key<GenericS5Key>( index_map ) = iS5;
         for( const auto &iS7 : index_range<gS7>() ) {
         at_key<GenericS7Key>( index_map ) = iS7;
         for( const auto &iS8 : index_range<gS8>() ) {
         at_key<GenericS8Key>( index_map ) = iS8;
         
            mF6 = context.mass<gF6>(iF6);
            mS5 = context.mass<gS5>(iS5);
            mS7 = context.mass<gS7>(iS7);
            mS8 = context.mass<gS8>(iS8);
            
            g1 = NPF_L(_Fe, gF6, _gS8) NPF_I(i3, iF6, iS8);
            g2 = NPF_R(_Fe, gF6, _gS8) NPF_I(i3, iF6, iS8);
            g3 = NPF_L(_Fe, Fe, _gS5) NPF_I(i4, i2, iS5);
            g4 = NPF_R(_Fe, Fe, _gS5) NPF_I(i4, i2, iS5);
            g5 = NPF_L(_gF6, Fe, _gS7) NPF_I(iF6, i1, iS7);
            g6 = NPF_R(_gF6, Fe, _gS7) NPF_I(iF6, i1, iS7);
            g7 = NPF_S(gS5, gS7, gS8) NPF_I(iS5, iS7, iS8);
            
            z[0] = (std::abs(g1) < std::numeric_limits<double>::epsilon())?0:1;
            z[1] = (std::abs(g2) < std::numeric_limits<double>::epsilon())?0:1;
            z[2] = (std::abs(g3) < std::numeric_limits<double>::epsilon())?0:1;
            z[3] = (std::abs(g4) < std::numeric_limits<double>::epsilon())?0:1;
            z[4] = (std::abs(g5) < std::numeric_limits<double>::epsilon())?0:1;
            z[5] = (std::abs(g6) < std::numeric_limits<double>::epsilon())?0:1;
            z[6] = (std::abs(g7) < std::numeric_limits<double>::epsilon())?0:1;
            
            if( (z[2] + z[3])*(z[0]*z[4] + z[1]*z[5])*z[6] == 0 ) continue;
            Loop_library::get().C(l1,0,0,0,Sqr(mF6),Sqr(mS8),Sqr(mS7),Sqr(context.scale()));
         
            S_LL += (0.006332573977646111*g1*g3*g5*g7*l1[cc0]*mF6)/Sqr(mS5);
            S_LR += (0.006332573977646111*g1*g4*g5*g7*l1[cc0]*mF6)/Sqr(mS5);
            S_RL += (0.006332573977646111*g2*g3*g6*g7*l1[cc0]*mF6)/Sqr(mS5);
            S_RR += (0.006332573977646111*g2*g4*g6*g7*l1[cc0]*mF6)/Sqr(mS5);
            V_LL += 0;
            V_LR += 0;
            V_RL += 0;
            V_RR += 0;
            T_LL += 0;
            T_RR += 0;
         
         }}}} // End of summation over generic fields
   
         return {S_LL, S_LR, S_RL, S_RR, V_LL, V_LR, V_RL, V_RR, T_LL, T_RR};
      } // End of operator()( void )
   }; // End of struct genericSum2_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum2( void ) {
      using GenericKeys = boost::mpl::vector< GenericF6Key, GenericS5Key, GenericS7Key, GenericS8Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<fields::Fv, fields::hh, fields::Hpm, typename cxx_diagrams::conj<fields::Hpm>::type>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Cha1>::type, fields::hh, fields::Sv, typename cxx_diagrams::conj<fields::Sv>::type>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Chi>::type, fields::hh, fields::Se, typename cxx_diagrams::conj<fields::Se>::type>,
         boost::mpl::vector<fields::Chi, fields::hh, fields::Se, typename cxx_diagrams::conj<fields::Se>::type>,
         boost::mpl::vector<fields::Fv, fields::Ah, fields::Hpm, typename cxx_diagrams::conj<fields::Hpm>::type>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Cha1>::type, fields::Ah, fields::Sv, typename cxx_diagrams::conj<fields::Sv>::type>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Chi>::type, fields::Ah, fields::Se, typename cxx_diagrams::conj<fields::Se>::type>,
         boost::mpl::vector<fields::Chi, fields::Ah, fields::Se, typename cxx_diagrams::conj<fields::Se>::type>
         >;
      using combinatorial_factors = boost::mpl::vector<
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>
         >;
      using color_factors = boost::mpl::vector<
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>
         >;
      return accumulate_generic<
         GenericKeys,
         GenericInsertions,
         combinatorial_factors,
         color_factors,
         boost::mpl::int_<10>,
         genericSum2_impl
         >( *this );
   } // End of function genericSum2()
   
   template<class GenericFieldMap>
   struct genericSum3_impl : generic_sum_base {
      genericSum3_impl( const generic_sum_base &base ) :
      generic_sum_base( base ) {
      } // End of constructor genericSum3_impl
   
      std::array<std::complex<double>,10> operator()( void ) {
         using boost::mpl::at;
         using boost::fusion::at_key;
         using gF6 = typename at<GenericFieldMap,GenericF6Key>::type;
         using gF7 = typename at<GenericFieldMap,GenericF7Key>::type;
         using gS5 = typename at<GenericFieldMap,GenericS5Key>::type;
         using gS8 = typename at<GenericFieldMap,GenericS8Key>::type;
   
         std::array<int, 1> i1 {this->external_indices(0)};
         std::array<int, 1> i2 {this->external_indices(1)};
         std::array<int, 1> i3 {this->external_indices(2)};
         std::array<int, 1> i4 {this->external_indices(3)};
         
         typename field_index_map<GenericFieldMap>::type index_map;
         const context_with_vertices &context = *this;
         std::complex<double> S_LL = 0.0;
         std::complex<double> S_LR = 0.0;
         std::complex<double> S_RL = 0.0;
         std::complex<double> S_RR = 0.0;
         std::complex<double> V_LL = 0.0;
         std::complex<double> V_LR = 0.0;
         std::complex<double> V_RL = 0.0;
         std::complex<double> V_RR = 0.0;
         std::complex<double> T_LL = 0.0;
         std::complex<double> T_RR = 0.0;
   
         // Shorter aliases for large types
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _gF6 = typename cxx_diagrams::bar<gF6>::type;
         using _gF7 = typename cxx_diagrams::bar<gF7>::type;
         using _gS5 = typename cxx_diagrams::conj<gS5>::type;
         using _gS8 = typename cxx_diagrams::conj<gS8>::type;
         using Fe = fields::Fe;
         using Fe = fields::Fe;
         
         // Aliases for loop function identifiers
         using looplibrary::bb0;
         
         double mF6, mF7, mS5, mS8;
         std::complex<double> g1, g2, g3, g4, g5, g6, g7, g8;
         std::array<int, 8> z{};
         looplibrary::Bcoeff_t l1 {};
         // Start of summation over generic fields.
         for( const auto &iF6 : index_range<gF6>() ) {
         at_key<GenericF6Key>( index_map ) = iF6;
         if( (std::is_same_v<gF6,fields::Fe> || std::is_same_v<gF6,typename cxx_diagrams::bar<fields::Fe>::type>) && iF6 == i3 ) continue;
         for( const auto &iF7 : index_range<gF7>() ) {
         at_key<GenericF7Key>( index_map ) = iF7;
         for( const auto &iS5 : index_range<gS5>() ) {
         at_key<GenericS5Key>( index_map ) = iS5;
         for( const auto &iS8 : index_range<gS8>() ) {
         at_key<GenericS8Key>( index_map ) = iS8;
         
            mF6 = context.mass<gF6>(iF6);
            mF7 = context.mass<gF7>(iF7);
            mS5 = context.mass<gS5>(iS5);
            mS8 = context.mass<gS8>(iS8);
            
            g1 = NPF_L(_Fe, _gF7, _gS8) NPF_I(i3, iF7, iS8);
            g2 = NPF_R(_Fe, _gF7, _gS8) NPF_I(i3, iF7, iS8);
            g3 = NPF_L(_Fe, Fe, gS5) NPF_I(i4, i2, iS5);
            g4 = NPF_R(_Fe, Fe, gS5) NPF_I(i4, i2, iS5);
            g5 = NPF_L(_gF6, Fe, _gS5) NPF_I(iF6, i1, iS5);
            g6 = NPF_R(_gF6, Fe, _gS5) NPF_I(iF6, i1, iS5);
            g7 = NPF_L(gF7, gF6, gS8) NPF_I(iF7, iF6, iS8);
            g8 = NPF_R(gF7, gF6, gS8) NPF_I(iF7, iF6, iS8);
            
            z[0] = (std::abs(g1) < std::numeric_limits<double>::epsilon())?0:1;
            z[1] = (std::abs(g2) < std::numeric_limits<double>::epsilon())?0:1;
            z[2] = (std::abs(g3) < std::numeric_limits<double>::epsilon())?0:1;
            z[3] = (std::abs(g4) < std::numeric_limits<double>::epsilon())?0:1;
            z[4] = (std::abs(g5) < std::numeric_limits<double>::epsilon())?0:1;
            z[5] = (std::abs(g6) < std::numeric_limits<double>::epsilon())?0:1;
            z[6] = (std::abs(g7) < std::numeric_limits<double>::epsilon())?0:1;
            z[7] = (std::abs(g8) < std::numeric_limits<double>::epsilon())?0:1;
            
            if( (z[2] + z[3])*(z[0]*z[4]*z[6] + z[1]*z[5]*z[7]) == 0 ) continue;
            Loop_library::get().B(l1,0,Sqr(mF7),Sqr(mS8),Sqr(context.scale()));
         
            S_LL += (-0.006332573977646111*g1*g3*g5*g7*l1[bb0]*mF7)/(mF6*Sqr(mS5));
            S_LR += (-0.006332573977646111*g1*g4*g5*g7*l1[bb0]*mF7)/(mF6*Sqr(mS5));
            S_RL += (-0.006332573977646111*g2*g3*g6*g8*l1[bb0]*mF7)/(mF6*Sqr(mS5));
            S_RR += (-0.006332573977646111*g2*g4*g6*g8*l1[bb0]*mF7)/(mF6*Sqr(mS5));
            V_LL += 0;
            V_LR += 0;
            V_RL += 0;
            V_RR += 0;
            T_LL += 0;
            T_RR += 0;
         
         }}}} // End of summation over generic fields
   
         return {S_LL, S_LR, S_RL, S_RR, V_LL, V_LR, V_RL, V_RR, T_LL, T_RR};
      } // End of operator()( void )
   }; // End of struct genericSum3_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum3( void ) {
      using GenericKeys = boost::mpl::vector< GenericF6Key, GenericF7Key, GenericS5Key, GenericS8Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<fields::Fe, typename cxx_diagrams::bar<fields::Fv>::type, fields::hh, typename cxx_diagrams::conj<fields::Hpm>::type>,
         boost::mpl::vector<fields::Fe, fields::Cha1, fields::hh, typename cxx_diagrams::conj<fields::Sv>::type>,
         boost::mpl::vector<fields::Fe, typename cxx_diagrams::bar<fields::Chi>::type, fields::hh, typename cxx_diagrams::conj<fields::Se>::type>,
         boost::mpl::vector<fields::Fe, fields::Chi, fields::hh, typename cxx_diagrams::conj<fields::Se>::type>,
         boost::mpl::vector<fields::Fe, typename cxx_diagrams::bar<fields::Fv>::type, fields::Ah, typename cxx_diagrams::conj<fields::Hpm>::type>,
         boost::mpl::vector<fields::Fe, fields::Cha1, fields::Ah, typename cxx_diagrams::conj<fields::Sv>::type>,
         boost::mpl::vector<fields::Fe, typename cxx_diagrams::bar<fields::Chi>::type, fields::Ah, typename cxx_diagrams::conj<fields::Se>::type>,
         boost::mpl::vector<fields::Fe, fields::Chi, fields::Ah, typename cxx_diagrams::conj<fields::Se>::type>
         >;
      using combinatorial_factors = boost::mpl::vector<
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>
         >;
      using color_factors = boost::mpl::vector<
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>
         >;
      return accumulate_generic<
         GenericKeys,
         GenericInsertions,
         combinatorial_factors,
         color_factors,
         boost::mpl::int_<10>,
         genericSum3_impl
         >( *this );
   } // End of function genericSum3()
   
   template<class GenericFieldMap>
   struct genericSum4_impl : generic_sum_base {
      genericSum4_impl( const generic_sum_base &base ) :
      generic_sum_base( base ) {
      } // End of constructor genericSum4_impl
   
      std::array<std::complex<double>,10> operator()( void ) {
         using boost::mpl::at;
         using boost::fusion::at_key;
         using gF6 = typename at<GenericFieldMap,GenericF6Key>::type;
         using gF7 = typename at<GenericFieldMap,GenericF7Key>::type;
         using gS5 = typename at<GenericFieldMap,GenericS5Key>::type;
         using gS8 = typename at<GenericFieldMap,GenericS8Key>::type;
   
         std::array<int, 1> i1 {this->external_indices(0)};
         std::array<int, 1> i2 {this->external_indices(1)};
         std::array<int, 1> i3 {this->external_indices(2)};
         std::array<int, 1> i4 {this->external_indices(3)};
         
         typename field_index_map<GenericFieldMap>::type index_map;
         const context_with_vertices &context = *this;
         std::complex<double> S_LL = 0.0;
         std::complex<double> S_LR = 0.0;
         std::complex<double> S_RL = 0.0;
         std::complex<double> S_RR = 0.0;
         std::complex<double> V_LL = 0.0;
         std::complex<double> V_LR = 0.0;
         std::complex<double> V_RL = 0.0;
         std::complex<double> V_RR = 0.0;
         std::complex<double> T_LL = 0.0;
         std::complex<double> T_RR = 0.0;
   
         // Shorter aliases for large types
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _gF6 = typename cxx_diagrams::bar<gF6>::type;
         using _gF7 = typename cxx_diagrams::bar<gF7>::type;
         using _gS5 = typename cxx_diagrams::conj<gS5>::type;
         using _gS8 = typename cxx_diagrams::conj<gS8>::type;
         using Fe = fields::Fe;
         using Fe = fields::Fe;
         
         // Aliases for loop function identifiers
         using looplibrary::bb0;
         using looplibrary::bb1;
         
         double mFe1, mF6, mF7, mS5, mS8;
         std::complex<double> g1, g2, g3, g4, g5, g6, g7, g8;
         std::array<int, 8> z{};
         looplibrary::Bcoeff_t l1 {};
         // Start of summation over generic fields.
         for( const auto &iF6 : index_range<gF6>() ) {
         at_key<GenericF6Key>( index_map ) = iF6;
         if( (std::is_same_v<gF6,fields::Fe> || std::is_same_v<gF6,typename cxx_diagrams::bar<fields::Fe>::type>) && iF6 == i1 ) continue;
         for( const auto &iF7 : index_range<gF7>() ) {
         at_key<GenericF7Key>( index_map ) = iF7;
         for( const auto &iS5 : index_range<gS5>() ) {
         at_key<GenericS5Key>( index_map ) = iS5;
         for( const auto &iS8 : index_range<gS8>() ) {
         at_key<GenericS8Key>( index_map ) = iS8;
         
            mFe1 = context.mass<Fe>(i1);
            mF6 = context.mass<gF6>(iF6);
            mF7 = context.mass<gF7>(iF7);
            mS5 = context.mass<gS5>(iS5);
            mS8 = context.mass<gS8>(iS8);
            
            g1 = NPF_L(_Fe, _gF6, gS5) NPF_I(i3, iF6, iS5);
            g2 = NPF_R(_Fe, _gF6, gS5) NPF_I(i3, iF6, iS5);
            g3 = NPF_L(_Fe, Fe, _gS5) NPF_I(i4, i2, iS5);
            g4 = NPF_R(_Fe, Fe, _gS5) NPF_I(i4, i2, iS5);
            g5 = NPF_L(_gF7, Fe, _gS8) NPF_I(iF7, i1, iS8);
            g6 = NPF_R(_gF7, Fe, _gS8) NPF_I(iF7, i1, iS8);
            g7 = NPF_L(gF6, gF7, gS8) NPF_I(iF6, iF7, iS8);
            g8 = NPF_R(gF6, gF7, gS8) NPF_I(iF6, iF7, iS8);
            
            z[0] = (std::abs(g1) < std::numeric_limits<double>::epsilon())?0:1;
            z[1] = (std::abs(g2) < std::numeric_limits<double>::epsilon())?0:1;
            z[2] = (std::abs(g3) < std::numeric_limits<double>::epsilon())?0:1;
            z[3] = (std::abs(g4) < std::numeric_limits<double>::epsilon())?0:1;
            z[4] = (std::abs(g5) < std::numeric_limits<double>::epsilon())?0:1;
            z[5] = (std::abs(g6) < std::numeric_limits<double>::epsilon())?0:1;
            z[6] = (std::abs(g7) < std::numeric_limits<double>::epsilon())?0:1;
            z[7] = (std::abs(g8) < std::numeric_limits<double>::epsilon())?0:1;
            
            if( (z[0] + z[1])*(z[2] + z[3])*(z[4] + z[5])*(z[6] + z[7]) == 0 ) continue;
            Loop_library::get().B(l1,0,Sqr(mF7),Sqr(mS8),Sqr(context.scale()));
         
            S_LL += (0.006332573977646111*g1*g3*(g5*g7*l1[bb0]*mF6*mF7 - g6*g7*l1[bb1]*mF6*mFe1 + g6*g8*l1[bb0]*mF7*mFe1 - g5*g8*l1[bb1]*Sqr(mFe1)))/((Sqr(mF6) - Sqr(mFe1))*(Sqr(mFe1) - Sqr(mS5)));
            S_LR += (0.006332573977646111*g1*g4*(g5*g7*l1[bb0]*mF6*mF7 - g6*g7*l1[bb1]*mF6*mFe1 + g6*g8*l1[bb0]*mF7*mFe1 - g5*g8*l1[bb1]*Sqr(mFe1)))/((Sqr(mF6) - Sqr(mFe1))*(Sqr(mFe1) - Sqr(mS5)));
            S_RL += (0.006332573977646111*g2*g3*(g6*g8*l1[bb0]*mF6*mF7 - g5*g8*l1[bb1]*mF6*mFe1 + g5*g7*l1[bb0]*mF7*mFe1 - g6*g7*l1[bb1]*Sqr(mFe1)))/((Sqr(mF6) - Sqr(mFe1))*(Sqr(mFe1) - Sqr(mS5)));
            S_RR += (0.006332573977646111*g2*g4*(g6*g8*l1[bb0]*mF6*mF7 - g5*g8*l1[bb1]*mF6*mFe1 + g5*g7*l1[bb0]*mF7*mFe1 - g6*g7*l1[bb1]*Sqr(mFe1)))/((Sqr(mF6) - Sqr(mFe1))*(Sqr(mFe1) - Sqr(mS5)));
            V_LL += 0;
            V_LR += 0;
            V_RL += 0;
            V_RR += 0;
            T_LL += 0;
            T_RR += 0;
         
         }}}} // End of summation over generic fields
   
         return {S_LL, S_LR, S_RL, S_RR, V_LL, V_LR, V_RL, V_RR, T_LL, T_RR};
      } // End of operator()( void )
   }; // End of struct genericSum4_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum4( void ) {
      using GenericKeys = boost::mpl::vector< GenericF6Key, GenericF7Key, GenericS5Key, GenericS8Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, fields::Fv, fields::hh, fields::Hpm>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, typename cxx_diagrams::bar<fields::Cha1>::type, fields::hh, fields::Sv>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, typename cxx_diagrams::bar<fields::Chi>::type, fields::hh, fields::Se>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, fields::Chi, fields::hh, fields::Se>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, fields::Fv, fields::Ah, fields::Hpm>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, typename cxx_diagrams::bar<fields::Cha1>::type, fields::Ah, fields::Sv>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, typename cxx_diagrams::bar<fields::Chi>::type, fields::Ah, fields::Se>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, fields::Chi, fields::Ah, fields::Se>
         >;
      using combinatorial_factors = boost::mpl::vector<
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>
         >;
      using color_factors = boost::mpl::vector<
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>
         >;
      return accumulate_generic<
         GenericKeys,
         GenericInsertions,
         combinatorial_factors,
         color_factors,
         boost::mpl::int_<10>,
         genericSum4_impl
         >( *this );
   } // End of function genericSum4()

   public:
   nPointFeFeFeFe_3917882094413845( const MRSSM2_mass_eigenstates &model,
   const std::array<int, 4> &indices,
   const std::array<Eigen::Vector4d, 0> &momenta
    ) :
   correlation_function_context<4,0> { model, indices, momenta }
   {}

   std::array<std::complex<double>,10> calculate( void ) {
      std::array<std::complex<double>,10> genericSummation;
      constexpr int coeffsLength = genericSummation.size();
      const auto genericsum1 = genericSum1();
      const auto genericsum2 = genericSum2();
      const auto genericsum3 = genericSum3();
      const auto genericsum4 = genericSum4();
   
      for ( std::size_t i=0; i<coeffsLength; i++ ) {
         genericSummation[i] += genericsum1[i]+genericsum2[i]+genericsum3[i]+genericsum4[i];
      }
      return genericSummation;
   } // End of calculate()
}; // End of nPointFeFeFeFe_3917882094413845

std::array<std::complex<double>, 10> scalars1loop(const MRSSM2_mass_eigenstates &model,
const std::array<int, 4> &indices,
const std::array<Eigen::Vector4d, 0> &momenta
) {
   nPointFeFeFeFe_3917882094413845 helper{model, indices, momenta};
   return helper.calculate();
}
class nPointFeFeFeFe_3917882099255229 : public correlation_function_context<4,0>
{
   using generic_sum_base = correlation_function_context<4,0>;

   template<class GenericFieldMap>
   struct subexpression_base :
   generic_sum_base, index_map_interface<GenericFieldMap>
   {
      subexpression_base( const subexpression_base & ) = default;

      subexpression_base( const generic_sum_base &gsb,
         const typename field_index_map<GenericFieldMap>::type &fim ) :
      generic_sum_base( gsb ), index_map_interface<GenericFieldMap>( fim )
      {}
   }; // End of subexpression_base<GenericFieldMap>

   struct GenericF7Key {};
   struct GenericF8Key {};
   struct GenericS6Key {};
   struct GenericV5Key {};
   struct GenericF6Key {};
   struct GenericS7Key {};
   struct GenericS8Key {};

   template<class GenericFieldMap>
   struct genericSum1_impl : generic_sum_base {
      genericSum1_impl( const generic_sum_base &base ) :
      generic_sum_base( base ) {
      } // End of constructor genericSum1_impl
   
      std::array<std::complex<double>,10> operator()( void ) {
         using boost::mpl::at;
         using boost::fusion::at_key;
         using gF7 = typename at<GenericFieldMap,GenericF7Key>::type;
         using gF8 = typename at<GenericFieldMap,GenericF8Key>::type;
         using gS6 = typename at<GenericFieldMap,GenericS6Key>::type;
         using gV5 = typename at<GenericFieldMap,GenericV5Key>::type;
   
         std::array<int, 1> i1 {this->external_indices(0)};
         std::array<int, 1> i2 {this->external_indices(1)};
         std::array<int, 1> i3 {this->external_indices(2)};
         std::array<int, 1> i4 {this->external_indices(3)};
         
         typename field_index_map<GenericFieldMap>::type index_map;
         const context_with_vertices &context = *this;
         std::complex<double> S_LL = 0.0;
         std::complex<double> S_LR = 0.0;
         std::complex<double> S_RL = 0.0;
         std::complex<double> S_RR = 0.0;
         std::complex<double> V_LL = 0.0;
         std::complex<double> V_LR = 0.0;
         std::complex<double> V_RL = 0.0;
         std::complex<double> V_RR = 0.0;
         std::complex<double> T_LL = 0.0;
         std::complex<double> T_RR = 0.0;
   
         // Shorter aliases for large types
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _gF7 = typename cxx_diagrams::bar<gF7>::type;
         using _gF8 = typename cxx_diagrams::bar<gF8>::type;
         using _gS6 = typename cxx_diagrams::conj<gS6>::type;
         using _gV5 = typename cxx_diagrams::conj<gV5>::type;
         using Fe = fields::Fe;
         using Fe = fields::Fe;
         
         // Aliases for loop function identifiers
         using looplibrary::bb0;
         using looplibrary::cc0;
         using looplibrary::cc00;
         
         double mF7, mF8, mS6, mV5;
         std::complex<double> g1, g2, g3, g4, g5, g6, g7, g8;
         std::array<int, 8> z{};
         looplibrary::Bcoeff_t l1 {};
         looplibrary::Ccoeff_t l2 {};
         // Start of summation over generic fields.
         for( const auto &iF7 : index_range<gF7>() ) {
         at_key<GenericF7Key>( index_map ) = iF7;
         for( const auto &iF8 : index_range<gF8>() ) {
         at_key<GenericF8Key>( index_map ) = iF8;
         for( const auto &iS6 : index_range<gS6>() ) {
         at_key<GenericS6Key>( index_map ) = iS6;
         for( const auto &iV5 : index_range<gV5>() ) {
         at_key<GenericV5Key>( index_map ) = iV5;
         
            mF7 = context.mass<gF7>(iF7);
            mF8 = context.mass<gF8>(iF8);
            mS6 = context.mass<gS6>(iS6);
            mV5 = context.mass<gV5>(iV5);
            
            g1 = NPF_L(_Fe, _gF8, gS6) NPF_I(i3, iF8, iS6);
            g2 = NPF_R(_Fe, _gF8, gS6) NPF_I(i3, iF8, iS6);
            g3 = NPF_L(_Fe, Fe, _gV5) NPF_I(i4, i2, iV5);
            g4 = NPF_R(_Fe, Fe, _gV5) NPF_I(i4, i2, iV5);
            g5 = NPF_L(_gF7, Fe, _gS6) NPF_I(iF7, i1, iS6);
            g6 = NPF_R(_gF7, Fe, _gS6) NPF_I(iF7, i1, iS6);
            g7 = NPF_L(gF8, gF7, gV5) NPF_I(iF8, iF7, iV5);
            g8 = NPF_R(gF8, gF7, gV5) NPF_I(iF8, iF7, iV5);
            
            z[0] = (std::abs(g1) < std::numeric_limits<double>::epsilon())?0:1;
            z[1] = (std::abs(g2) < std::numeric_limits<double>::epsilon())?0:1;
            z[2] = (std::abs(g3) < std::numeric_limits<double>::epsilon())?0:1;
            z[3] = (std::abs(g4) < std::numeric_limits<double>::epsilon())?0:1;
            z[4] = (std::abs(g5) < std::numeric_limits<double>::epsilon())?0:1;
            z[5] = (std::abs(g6) < std::numeric_limits<double>::epsilon())?0:1;
            z[6] = (std::abs(g7) < std::numeric_limits<double>::epsilon())?0:1;
            z[7] = (std::abs(g8) < std::numeric_limits<double>::epsilon())?0:1;
            
            if( (z[2] + z[3])*(z[1]*z[4] + z[0]*z[5])*(z[6] + z[7]) == 0 ) continue;
            Loop_library::get().B(l1,0,Sqr(mF7),Sqr(mF8),Sqr(context.scale()));
            Loop_library::get().C(l2,0,0,0,Sqr(mF8),Sqr(mF7),Sqr(mS6),Sqr(context.scale()));
         
            S_LL += 0;
            S_LR += 0;
            S_RL += 0;
            S_RR += 0;
            V_LL += (0.006332573977646111*g2*g3*g5*(-(g7*l2[cc0]*mF7*mF8) + g8*(l1[bb0] - 2*l2[cc00] + l2[cc0]*Sqr(mS6))))/Sqr(mV5);
            V_LR += (0.006332573977646111*g2*g4*g5*(-(g7*l2[cc0]*mF7*mF8) + g8*(l1[bb0] - 2*l2[cc00] + l2[cc0]*Sqr(mS6))))/Sqr(mV5);
            V_RL += (0.006332573977646111*g1*g3*g6*(-(g8*l2[cc0]*mF7*mF8) + g7*(l1[bb0] - 2*l2[cc00] + l2[cc0]*Sqr(mS6))))/Sqr(mV5);
            V_RR += (0.006332573977646111*g1*g4*g6*(-(g8*l2[cc0]*mF7*mF8) + g7*(l1[bb0] - 2*l2[cc00] + l2[cc0]*Sqr(mS6))))/Sqr(mV5);
            T_LL += 0;
            T_RR += 0;
         
         }}}} // End of summation over generic fields
   
         return {S_LL, S_LR, S_RL, S_RR, V_LL, V_LR, V_RL, V_RR, T_LL, T_RR};
      } // End of operator()( void )
   }; // End of struct genericSum1_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum1( void ) {
      using GenericKeys = boost::mpl::vector< GenericF7Key, GenericF8Key, GenericS6Key, GenericV5Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<fields::Fv, typename cxx_diagrams::bar<fields::Fv>::type, fields::Hpm, fields::VZ>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Cha1>::type, fields::Cha1, fields::Sv, fields::VZ>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Chi>::type, fields::Chi, fields::Se, fields::VZ>,
         boost::mpl::vector<fields::Chi, typename cxx_diagrams::bar<fields::Chi>::type, fields::Se, fields::VZ>
         >;
      using combinatorial_factors = boost::mpl::vector<
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>
         >;
      using color_factors = boost::mpl::vector<
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>
         >;
      return accumulate_generic<
         GenericKeys,
         GenericInsertions,
         combinatorial_factors,
         color_factors,
         boost::mpl::int_<10>,
         genericSum1_impl
         >( *this );
   } // End of function genericSum1()
   
   template<class GenericFieldMap>
   struct genericSum2_impl : generic_sum_base {
      genericSum2_impl( const generic_sum_base &base ) :
      generic_sum_base( base ) {
      } // End of constructor genericSum2_impl
   
      std::array<std::complex<double>,10> operator()( void ) {
         using boost::mpl::at;
         using boost::fusion::at_key;
         using gF6 = typename at<GenericFieldMap,GenericF6Key>::type;
         using gS7 = typename at<GenericFieldMap,GenericS7Key>::type;
         using gS8 = typename at<GenericFieldMap,GenericS8Key>::type;
         using gV5 = typename at<GenericFieldMap,GenericV5Key>::type;
   
         std::array<int, 1> i1 {this->external_indices(0)};
         std::array<int, 1> i2 {this->external_indices(1)};
         std::array<int, 1> i3 {this->external_indices(2)};
         std::array<int, 1> i4 {this->external_indices(3)};
         
         typename field_index_map<GenericFieldMap>::type index_map;
         const context_with_vertices &context = *this;
         std::complex<double> S_LL = 0.0;
         std::complex<double> S_LR = 0.0;
         std::complex<double> S_RL = 0.0;
         std::complex<double> S_RR = 0.0;
         std::complex<double> V_LL = 0.0;
         std::complex<double> V_LR = 0.0;
         std::complex<double> V_RL = 0.0;
         std::complex<double> V_RR = 0.0;
         std::complex<double> T_LL = 0.0;
         std::complex<double> T_RR = 0.0;
   
         // Shorter aliases for large types
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _gF6 = typename cxx_diagrams::bar<gF6>::type;
         using _gS7 = typename cxx_diagrams::conj<gS7>::type;
         using _gS8 = typename cxx_diagrams::conj<gS8>::type;
         using _gV5 = typename cxx_diagrams::conj<gV5>::type;
         using Fe = fields::Fe;
         using Fe = fields::Fe;
         
         // Aliases for loop function identifiers
         using looplibrary::cc00;
         
         double mF6, mS7, mS8, mV5;
         std::complex<double> g1, g2, g3, g4, g5, g6, g7;
         std::array<int, 7> z{};
         looplibrary::Ccoeff_t l1 {};
         // Start of summation over generic fields.
         for( const auto &iF6 : index_range<gF6>() ) {
         at_key<GenericF6Key>( index_map ) = iF6;
         for( const auto &iS7 : index_range<gS7>() ) {
         at_key<GenericS7Key>( index_map ) = iS7;
         for( const auto &iS8 : index_range<gS8>() ) {
         at_key<GenericS8Key>( index_map ) = iS8;
         for( const auto &iV5 : index_range<gV5>() ) {
         at_key<GenericV5Key>( index_map ) = iV5;
         
            mF6 = context.mass<gF6>(iF6);
            mS7 = context.mass<gS7>(iS7);
            mS8 = context.mass<gS8>(iS8);
            mV5 = context.mass<gV5>(iV5);
            
            g1 = NPF_L(_Fe, gF6, _gS8) NPF_I(i3, iF6, iS8);
            g2 = NPF_R(_Fe, gF6, _gS8) NPF_I(i3, iF6, iS8);
            g3 = NPF_L(_Fe, Fe, _gV5) NPF_I(i4, i2, iV5);
            g4 = NPF_R(_Fe, Fe, _gV5) NPF_I(i4, i2, iV5);
            g5 = NPF_L(_gF6, Fe, _gS7) NPF_I(iF6, i1, iS7);
            g6 = NPF_R(_gF6, Fe, _gS7) NPF_I(iF6, i1, iS7);
            g7 = NPF_SSV(gS7, gS8, gV5) NPF_D(0, 1) NPF_I(iS7, iS8, iV5);
            
            z[0] = (std::abs(g1) < std::numeric_limits<double>::epsilon())?0:1;
            z[1] = (std::abs(g2) < std::numeric_limits<double>::epsilon())?0:1;
            z[2] = (std::abs(g3) < std::numeric_limits<double>::epsilon())?0:1;
            z[3] = (std::abs(g4) < std::numeric_limits<double>::epsilon())?0:1;
            z[4] = (std::abs(g5) < std::numeric_limits<double>::epsilon())?0:1;
            z[5] = (std::abs(g6) < std::numeric_limits<double>::epsilon())?0:1;
            z[6] = (std::abs(g7) < std::numeric_limits<double>::epsilon())?0:1;
            
            if( (z[2] + z[3])*(z[1]*z[4] + z[0]*z[5])*z[6] == 0 ) continue;
            Loop_library::get().C(l1,0,0,0,Sqr(mF6),Sqr(mS8),Sqr(mS7),Sqr(context.scale()));
         
            S_LL += 0;
            S_LR += 0;
            S_RL += 0;
            S_RR += 0;
            V_LL += (0.012665147955292222*g2*g3*g5*g7*l1[cc00])/Sqr(mV5);
            V_LR += (0.012665147955292222*g2*g4*g5*g7*l1[cc00])/Sqr(mV5);
            V_RL += (0.012665147955292222*g1*g3*g6*g7*l1[cc00])/Sqr(mV5);
            V_RR += (0.012665147955292222*g1*g4*g6*g7*l1[cc00])/Sqr(mV5);
            T_LL += 0;
            T_RR += 0;
         
         }}}} // End of summation over generic fields
   
         return {S_LL, S_LR, S_RL, S_RR, V_LL, V_LR, V_RL, V_RR, T_LL, T_RR};
      } // End of operator()( void )
   }; // End of struct genericSum2_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum2( void ) {
      using GenericKeys = boost::mpl::vector< GenericF6Key, GenericS7Key, GenericS8Key, GenericV5Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<fields::Fv, fields::Hpm, typename cxx_diagrams::conj<fields::Hpm>::type, fields::VZ>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Cha1>::type, fields::Sv, typename cxx_diagrams::conj<fields::Sv>::type, fields::VZ>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Chi>::type, fields::Se, typename cxx_diagrams::conj<fields::Se>::type, fields::VZ>,
         boost::mpl::vector<fields::Chi, fields::Se, typename cxx_diagrams::conj<fields::Se>::type, fields::VZ>
         >;
      using combinatorial_factors = boost::mpl::vector<
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>
         >;
      using color_factors = boost::mpl::vector<
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>
         >;
      return accumulate_generic<
         GenericKeys,
         GenericInsertions,
         combinatorial_factors,
         color_factors,
         boost::mpl::int_<10>,
         genericSum2_impl
         >( *this );
   } // End of function genericSum2()
   
   template<class GenericFieldMap>
   struct genericSum3_impl : generic_sum_base {
      genericSum3_impl( const generic_sum_base &base ) :
      generic_sum_base( base ) {
      } // End of constructor genericSum3_impl
   
      std::array<std::complex<double>,10> operator()( void ) {
         using boost::mpl::at;
         using boost::fusion::at_key;
         using gF6 = typename at<GenericFieldMap,GenericF6Key>::type;
         using gF7 = typename at<GenericFieldMap,GenericF7Key>::type;
         using gS8 = typename at<GenericFieldMap,GenericS8Key>::type;
         using gV5 = typename at<GenericFieldMap,GenericV5Key>::type;
   
         std::array<int, 1> i1 {this->external_indices(0)};
         std::array<int, 1> i2 {this->external_indices(1)};
         std::array<int, 1> i3 {this->external_indices(2)};
         std::array<int, 1> i4 {this->external_indices(3)};
         
         typename field_index_map<GenericFieldMap>::type index_map;
         const context_with_vertices &context = *this;
         std::complex<double> S_LL = 0.0;
         std::complex<double> S_LR = 0.0;
         std::complex<double> S_RL = 0.0;
         std::complex<double> S_RR = 0.0;
         std::complex<double> V_LL = 0.0;
         std::complex<double> V_LR = 0.0;
         std::complex<double> V_RL = 0.0;
         std::complex<double> V_RR = 0.0;
         std::complex<double> T_LL = 0.0;
         std::complex<double> T_RR = 0.0;
   
         // Shorter aliases for large types
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _gF6 = typename cxx_diagrams::bar<gF6>::type;
         using _gF7 = typename cxx_diagrams::bar<gF7>::type;
         using _gS8 = typename cxx_diagrams::conj<gS8>::type;
         using _gV5 = typename cxx_diagrams::conj<gV5>::type;
         using Fe = fields::Fe;
         using Fe = fields::Fe;
         
         // Aliases for loop function identifiers
         using looplibrary::bb0;
         
         double mF6, mF7, mS8, mV5;
         std::complex<double> g1, g2, g3, g4, g5, g6, g7, g8;
         std::array<int, 8> z{};
         looplibrary::Bcoeff_t l1 {};
         // Start of summation over generic fields.
         for( const auto &iF6 : index_range<gF6>() ) {
         at_key<GenericF6Key>( index_map ) = iF6;
         if( (std::is_same_v<gF6,fields::Fe> || std::is_same_v<gF6,typename cxx_diagrams::bar<fields::Fe>::type>) && iF6 == i3 ) continue;
         for( const auto &iF7 : index_range<gF7>() ) {
         at_key<GenericF7Key>( index_map ) = iF7;
         for( const auto &iS8 : index_range<gS8>() ) {
         at_key<GenericS8Key>( index_map ) = iS8;
         for( const auto &iV5 : index_range<gV5>() ) {
         at_key<GenericV5Key>( index_map ) = iV5;
         
            mF6 = context.mass<gF6>(iF6);
            mF7 = context.mass<gF7>(iF7);
            mS8 = context.mass<gS8>(iS8);
            mV5 = context.mass<gV5>(iV5);
            
            g1 = NPF_L(_Fe, _gF7, _gS8) NPF_I(i3, iF7, iS8);
            g2 = NPF_R(_Fe, _gF7, _gS8) NPF_I(i3, iF7, iS8);
            g3 = NPF_L(_Fe, Fe, gV5) NPF_I(i4, i2, iV5);
            g4 = NPF_R(_Fe, Fe, gV5) NPF_I(i4, i2, iV5);
            g5 = NPF_L(_gF6, Fe, _gV5) NPF_I(iF6, i1, iV5);
            g6 = NPF_R(_gF6, Fe, _gV5) NPF_I(iF6, i1, iV5);
            g7 = NPF_L(gF7, gF6, gS8) NPF_I(iF7, iF6, iS8);
            g8 = NPF_R(gF7, gF6, gS8) NPF_I(iF7, iF6, iS8);
            
            z[0] = (std::abs(g1) < std::numeric_limits<double>::epsilon())?0:1;
            z[1] = (std::abs(g2) < std::numeric_limits<double>::epsilon())?0:1;
            z[2] = (std::abs(g3) < std::numeric_limits<double>::epsilon())?0:1;
            z[3] = (std::abs(g4) < std::numeric_limits<double>::epsilon())?0:1;
            z[4] = (std::abs(g5) < std::numeric_limits<double>::epsilon())?0:1;
            z[5] = (std::abs(g6) < std::numeric_limits<double>::epsilon())?0:1;
            z[6] = (std::abs(g7) < std::numeric_limits<double>::epsilon())?0:1;
            z[7] = (std::abs(g8) < std::numeric_limits<double>::epsilon())?0:1;
            
            if( (z[2] + z[3])*(z[0]*z[5]*z[6] + z[1]*z[4]*z[7]) == 0 ) continue;
            Loop_library::get().B(l1,0,Sqr(mF7),Sqr(mS8),Sqr(context.scale()));
         
            S_LL += 0;
            S_LR += 0;
            S_RL += 0;
            S_RR += 0;
            V_LL += (0.006332573977646111*g2*g3*g5*g8*l1[bb0]*mF7)/(mF6*Sqr(mV5));
            V_LR += (0.006332573977646111*g2*g4*g5*g8*l1[bb0]*mF7)/(mF6*Sqr(mV5));
            V_RL += (0.006332573977646111*g1*g3*g6*g7*l1[bb0]*mF7)/(mF6*Sqr(mV5));
            V_RR += (0.006332573977646111*g1*g4*g6*g7*l1[bb0]*mF7)/(mF6*Sqr(mV5));
            T_LL += 0;
            T_RR += 0;
         
         }}}} // End of summation over generic fields
   
         return {S_LL, S_LR, S_RL, S_RR, V_LL, V_LR, V_RL, V_RR, T_LL, T_RR};
      } // End of operator()( void )
   }; // End of struct genericSum3_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum3( void ) {
      using GenericKeys = boost::mpl::vector< GenericF6Key, GenericF7Key, GenericS8Key, GenericV5Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<fields::Fe, typename cxx_diagrams::bar<fields::Fv>::type, typename cxx_diagrams::conj<fields::Hpm>::type, fields::VZ>,
         boost::mpl::vector<fields::Fe, fields::Cha1, typename cxx_diagrams::conj<fields::Sv>::type, fields::VZ>,
         boost::mpl::vector<fields::Fe, typename cxx_diagrams::bar<fields::Chi>::type, typename cxx_diagrams::conj<fields::Se>::type, fields::VZ>,
         boost::mpl::vector<fields::Fe, fields::Chi, typename cxx_diagrams::conj<fields::Se>::type, fields::VZ>
         >;
      using combinatorial_factors = boost::mpl::vector<
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>
         >;
      using color_factors = boost::mpl::vector<
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>
         >;
      return accumulate_generic<
         GenericKeys,
         GenericInsertions,
         combinatorial_factors,
         color_factors,
         boost::mpl::int_<10>,
         genericSum3_impl
         >( *this );
   } // End of function genericSum3()
   
   template<class GenericFieldMap>
   struct genericSum4_impl : generic_sum_base {
      genericSum4_impl( const generic_sum_base &base ) :
      generic_sum_base( base ) {
      } // End of constructor genericSum4_impl
   
      std::array<std::complex<double>,10> operator()( void ) {
         using boost::mpl::at;
         using boost::fusion::at_key;
         using gF6 = typename at<GenericFieldMap,GenericF6Key>::type;
         using gF7 = typename at<GenericFieldMap,GenericF7Key>::type;
         using gS8 = typename at<GenericFieldMap,GenericS8Key>::type;
         using gV5 = typename at<GenericFieldMap,GenericV5Key>::type;
   
         std::array<int, 1> i1 {this->external_indices(0)};
         std::array<int, 1> i2 {this->external_indices(1)};
         std::array<int, 1> i3 {this->external_indices(2)};
         std::array<int, 1> i4 {this->external_indices(3)};
         
         typename field_index_map<GenericFieldMap>::type index_map;
         const context_with_vertices &context = *this;
         std::complex<double> S_LL = 0.0;
         std::complex<double> S_LR = 0.0;
         std::complex<double> S_RL = 0.0;
         std::complex<double> S_RR = 0.0;
         std::complex<double> V_LL = 0.0;
         std::complex<double> V_LR = 0.0;
         std::complex<double> V_RL = 0.0;
         std::complex<double> V_RR = 0.0;
         std::complex<double> T_LL = 0.0;
         std::complex<double> T_RR = 0.0;
   
         // Shorter aliases for large types
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _gF6 = typename cxx_diagrams::bar<gF6>::type;
         using _gF7 = typename cxx_diagrams::bar<gF7>::type;
         using _gS8 = typename cxx_diagrams::conj<gS8>::type;
         using _gV5 = typename cxx_diagrams::conj<gV5>::type;
         using Fe = fields::Fe;
         using Fe = fields::Fe;
         
         // Aliases for loop function identifiers
         using looplibrary::bb0;
         using looplibrary::bb1;
         
         double mFe1, mF6, mF7, mS8, mV5;
         std::complex<double> g1, g2, g3, g4, g5, g6, g7, g8;
         std::array<int, 8> z{};
         looplibrary::Bcoeff_t l1 {};
         // Start of summation over generic fields.
         for( const auto &iF6 : index_range<gF6>() ) {
         at_key<GenericF6Key>( index_map ) = iF6;
         if( (std::is_same_v<gF6,fields::Fe> || std::is_same_v<gF6,typename cxx_diagrams::bar<fields::Fe>::type>) && iF6 == i1 ) continue;
         for( const auto &iF7 : index_range<gF7>() ) {
         at_key<GenericF7Key>( index_map ) = iF7;
         for( const auto &iS8 : index_range<gS8>() ) {
         at_key<GenericS8Key>( index_map ) = iS8;
         for( const auto &iV5 : index_range<gV5>() ) {
         at_key<GenericV5Key>( index_map ) = iV5;
         
            mFe1 = context.mass<Fe>(i1);
            mF6 = context.mass<gF6>(iF6);
            mF7 = context.mass<gF7>(iF7);
            mS8 = context.mass<gS8>(iS8);
            mV5 = context.mass<gV5>(iV5);
            
            g1 = NPF_L(_Fe, _gF6, gV5) NPF_I(i3, iF6, iV5);
            g2 = NPF_R(_Fe, _gF6, gV5) NPF_I(i3, iF6, iV5);
            g3 = NPF_L(_Fe, Fe, _gV5) NPF_I(i4, i2, iV5);
            g4 = NPF_R(_Fe, Fe, _gV5) NPF_I(i4, i2, iV5);
            g5 = NPF_L(_gF7, Fe, _gS8) NPF_I(iF7, i1, iS8);
            g6 = NPF_R(_gF7, Fe, _gS8) NPF_I(iF7, i1, iS8);
            g7 = NPF_L(gF6, gF7, gS8) NPF_I(iF6, iF7, iS8);
            g8 = NPF_R(gF6, gF7, gS8) NPF_I(iF6, iF7, iS8);
            
            z[0] = (std::abs(g1) < std::numeric_limits<double>::epsilon())?0:1;
            z[1] = (std::abs(g2) < std::numeric_limits<double>::epsilon())?0:1;
            z[2] = (std::abs(g3) < std::numeric_limits<double>::epsilon())?0:1;
            z[3] = (std::abs(g4) < std::numeric_limits<double>::epsilon())?0:1;
            z[4] = (std::abs(g5) < std::numeric_limits<double>::epsilon())?0:1;
            z[5] = (std::abs(g6) < std::numeric_limits<double>::epsilon())?0:1;
            z[6] = (std::abs(g7) < std::numeric_limits<double>::epsilon())?0:1;
            z[7] = (std::abs(g8) < std::numeric_limits<double>::epsilon())?0:1;
            
            if( (z[0] + z[1])*(z[2] + z[3])*(z[4] + z[5])*(z[6] + z[7]) == 0 ) continue;
            Loop_library::get().B(l1,0,Sqr(mF7),Sqr(mS8),Sqr(context.scale()));
         
            S_LL += 0;
            S_LR += 0;
            S_RL += 0;
            S_RR += 0;
            V_LL += (0.006332573977646111*g1*g3*(-(g5*g7*l1[bb0]*mF6*mF7) + g6*g7*l1[bb1]*mF6*mFe1 - g6*g8*l1[bb0]*mF7*mFe1 + g5*g8*l1[bb1]*Sqr(mFe1)))/((Sqr(mF6) - Sqr(mFe1))*(Sqr(mFe1) - Sqr(mV5)));
            V_LR += (0.006332573977646111*g1*g4*(-(g5*g7*l1[bb0]*mF6*mF7) + g6*g7*l1[bb1]*mF6*mFe1 - g6*g8*l1[bb0]*mF7*mFe1 + g5*g8*l1[bb1]*Sqr(mFe1)))/((Sqr(mF6) - Sqr(mFe1))*(Sqr(mFe1) - Sqr(mV5)));
            V_RL += (0.006332573977646111*g2*g3*(-(g6*g8*l1[bb0]*mF6*mF7) + g5*g8*l1[bb1]*mF6*mFe1 - g5*g7*l1[bb0]*mF7*mFe1 + g6*g7*l1[bb1]*Sqr(mFe1)))/((Sqr(mF6) - Sqr(mFe1))*(Sqr(mFe1) - Sqr(mV5)));
            V_RR += (0.006332573977646111*g2*g4*(-(g6*g8*l1[bb0]*mF6*mF7) + g5*g8*l1[bb1]*mF6*mFe1 - g5*g7*l1[bb0]*mF7*mFe1 + g6*g7*l1[bb1]*Sqr(mFe1)))/((Sqr(mF6) - Sqr(mFe1))*(Sqr(mFe1) - Sqr(mV5)));
            T_LL += 0;
            T_RR += 0;
         
         }}}} // End of summation over generic fields
   
         return {S_LL, S_LR, S_RL, S_RR, V_LL, V_LR, V_RL, V_RR, T_LL, T_RR};
      } // End of operator()( void )
   }; // End of struct genericSum4_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum4( void ) {
      using GenericKeys = boost::mpl::vector< GenericF6Key, GenericF7Key, GenericS8Key, GenericV5Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, fields::Fv, fields::Hpm, fields::VZ>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, typename cxx_diagrams::bar<fields::Cha1>::type, fields::Sv, fields::VZ>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, typename cxx_diagrams::bar<fields::Chi>::type, fields::Se, fields::VZ>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, fields::Chi, fields::Se, fields::VZ>
         >;
      using combinatorial_factors = boost::mpl::vector<
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>
         >;
      using color_factors = boost::mpl::vector<
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>
         >;
      return accumulate_generic<
         GenericKeys,
         GenericInsertions,
         combinatorial_factors,
         color_factors,
         boost::mpl::int_<10>,
         genericSum4_impl
         >( *this );
   } // End of function genericSum4()

   public:
   nPointFeFeFeFe_3917882099255229( const MRSSM2_mass_eigenstates &model,
   const std::array<int, 4> &indices,
   const std::array<Eigen::Vector4d, 0> &momenta
    ) :
   correlation_function_context<4,0> { model, indices, momenta }
   {}

   std::array<std::complex<double>,10> calculate( void ) {
      std::array<std::complex<double>,10> genericSummation;
      constexpr int coeffsLength = genericSummation.size();
      const auto genericsum1 = genericSum1();
      const auto genericsum2 = genericSum2();
      const auto genericsum3 = genericSum3();
      const auto genericsum4 = genericSum4();
   
      for ( std::size_t i=0; i<coeffsLength; i++ ) {
         genericSummation[i] += genericsum1[i]+genericsum2[i]+genericsum3[i]+genericsum4[i];
      }
      return genericSummation;
   } // End of calculate()
}; // End of nPointFeFeFeFe_3917882099255229

std::array<std::complex<double>, 10> vectors1loop(const MRSSM2_mass_eigenstates &model,
const std::array<int, 4> &indices,
const std::array<Eigen::Vector4d, 0> &momenta
) {
   nPointFeFeFeFe_3917882099255229 helper{model, indices, momenta};
   return helper.calculate();
}
class nPointFeFeFeFe_3917882127499752 : public correlation_function_context<4,0>
{
   using generic_sum_base = correlation_function_context<4,0>;

   template<class GenericFieldMap>
   struct subexpression_base :
   generic_sum_base, index_map_interface<GenericFieldMap>
   {
      subexpression_base( const subexpression_base & ) = default;

      subexpression_base( const generic_sum_base &gsb,
         const typename field_index_map<GenericFieldMap>::type &fim ) :
      generic_sum_base( gsb ), index_map_interface<GenericFieldMap>( fim )
      {}
   }; // End of subexpression_base<GenericFieldMap>

   struct GenericF5Key {};
   struct GenericF8Key {};
   struct GenericS6Key {};
   struct GenericS7Key {};
   struct GenericF6Key {};
   struct GenericF7Key {};
   struct GenericS5Key {};
   struct GenericS8Key {};

   template<class GenericFieldMap>
   struct genericSum1_impl : generic_sum_base {
      genericSum1_impl( const generic_sum_base &base ) :
      generic_sum_base( base ) {
      } // End of constructor genericSum1_impl
   
      std::array<std::complex<double>,10> operator()( void ) {
         using boost::mpl::at;
         using boost::fusion::at_key;
         using gF5 = typename at<GenericFieldMap,GenericF5Key>::type;
         using gF8 = typename at<GenericFieldMap,GenericF8Key>::type;
         using gS6 = typename at<GenericFieldMap,GenericS6Key>::type;
         using gS7 = typename at<GenericFieldMap,GenericS7Key>::type;
   
         std::array<int, 1> i1 {this->external_indices(0)};
         std::array<int, 1> i2 {this->external_indices(1)};
         std::array<int, 1> i3 {this->external_indices(2)};
         std::array<int, 1> i4 {this->external_indices(3)};
         
         typename field_index_map<GenericFieldMap>::type index_map;
         const context_with_vertices &context = *this;
         std::complex<double> S_LL = 0.0;
         std::complex<double> S_LR = 0.0;
         std::complex<double> S_RL = 0.0;
         std::complex<double> S_RR = 0.0;
         std::complex<double> V_LL = 0.0;
         std::complex<double> V_LR = 0.0;
         std::complex<double> V_RL = 0.0;
         std::complex<double> V_RR = 0.0;
         std::complex<double> T_LL = 0.0;
         std::complex<double> T_RR = 0.0;
   
         // Shorter aliases for large types
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _gF5 = typename cxx_diagrams::bar<gF5>::type;
         using _gF8 = typename cxx_diagrams::bar<gF8>::type;
         using _gS6 = typename cxx_diagrams::conj<gS6>::type;
         using _gS7 = typename cxx_diagrams::conj<gS7>::type;
         using Fe = fields::Fe;
         using Fe = fields::Fe;
         
         // Aliases for loop function identifiers
         using looplibrary::dd0;
         using looplibrary::dd00;
         
         double mF5, mF8, mS6, mS7;
         std::complex<double> g1, g2, g3, g4, g5, g6, g7, g8;
         std::array<int, 8> z{};
         looplibrary::Dcoeff_t l1 {};
         // Start of summation over generic fields.
         for( const auto &iF5 : index_range<gF5>() ) {
         at_key<GenericF5Key>( index_map ) = iF5;
         for( const auto &iF8 : index_range<gF8>() ) {
         at_key<GenericF8Key>( index_map ) = iF8;
         for( const auto &iS6 : index_range<gS6>() ) {
         at_key<GenericS6Key>( index_map ) = iS6;
         for( const auto &iS7 : index_range<gS7>() ) {
         at_key<GenericS7Key>( index_map ) = iS7;
         
            mF5 = context.mass<gF5>(iF5);
            mF8 = context.mass<gF8>(iF8);
            mS6 = context.mass<gS6>(iS6);
            mS7 = context.mass<gS7>(iS7);
            
            g1 = NPF_L(_Fe, gF8, gS7) NPF_I(i4, iF8, iS7);
            g2 = NPF_R(_Fe, gF8, gS7) NPF_I(i4, iF8, iS7);
            g3 = NPF_L(_gF5, Fe, _gS6) NPF_I(iF5, i1, iS6);
            g4 = NPF_R(_gF5, Fe, _gS6) NPF_I(iF5, i1, iS6);
            g5 = NPF_L(_gF8, _Fe, gS6) NPF_I(iF8, i3, iS6);
            g6 = NPF_R(_gF8, _Fe, gS6) NPF_I(iF8, i3, iS6);
            g7 = NPF_L(Fe, gF5, _gS7) NPF_I(i2, iF5, iS7);
            g8 = NPF_R(Fe, gF5, _gS7) NPF_I(i2, iF5, iS7);
            
            z[0] = (std::abs(g1) < std::numeric_limits<double>::epsilon())?0:1;
            z[1] = (std::abs(g2) < std::numeric_limits<double>::epsilon())?0:1;
            z[2] = (std::abs(g3) < std::numeric_limits<double>::epsilon())?0:1;
            z[3] = (std::abs(g4) < std::numeric_limits<double>::epsilon())?0:1;
            z[4] = (std::abs(g5) < std::numeric_limits<double>::epsilon())?0:1;
            z[5] = (std::abs(g6) < std::numeric_limits<double>::epsilon())?0:1;
            z[6] = (std::abs(g7) < std::numeric_limits<double>::epsilon())?0:1;
            z[7] = (std::abs(g8) < std::numeric_limits<double>::epsilon())?0:1;
            
            if( z[0]*(z[2]*z[4]*z[6] + z[3]*z[5]*z[6] + z[3]*z[4]*z[7] + z[2]*z[5]*z[7]) + z[1]*(z[3]*z[4]*z[6] + z[2]*z[5]*z[6] + z[2]*z[4]*z[7] + z[3]*z[5]*z[7]) == 0 ) continue;
            Loop_library::get().D(l1,0,0,0,0,0,0,Sqr(mF8),Sqr(mF5),Sqr(mS7),Sqr(mS6),Sqr(context.scale()));
         
            S_LL += 0.0031662869888230555*g1*g3*g5*g7*l1[dd0]*mF5*mF8;
            S_LR += 0.012665147955292222*g2*g3*g5*g8*l1[dd00];
            S_RL += 0.012665147955292222*g1*g4*g6*g7*l1[dd00];
            S_RR += 0.0031662869888230555*g2*g4*g6*g8*l1[dd0]*mF5*mF8;
            V_LL += -0.0031662869888230555*g2*g3*g6*g7*l1[dd0]*mF5*mF8;
            V_LR += 0.006332573977646111*g1*g3*g6*g8*l1[dd00];
            V_RL += 0.006332573977646111*g2*g4*g5*g7*l1[dd00];
            V_RR += -0.0031662869888230555*g1*g4*g5*g8*l1[dd0]*mF5*mF8;
            T_LL += 0.0007915717472057639*g1*g3*g5*g7*l1[dd0]*mF5*mF8;
            T_RR += 0.0007915717472057639*g2*g4*g6*g8*l1[dd0]*mF5*mF8;
         
         }}}} // End of summation over generic fields
   
         return {S_LL, S_LR, S_RL, S_RR, V_LL, V_LR, V_RL, V_RR, T_LL, T_RR};
      } // End of operator()( void )
   }; // End of struct genericSum1_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum1( void ) {
      using GenericKeys = boost::mpl::vector< GenericF5Key, GenericF8Key, GenericS6Key, GenericS7Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Chi>::type, typename cxx_diagrams::bar<fields::Chi>::type, fields::Se, fields::Se>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Chi>::type, fields::Chi, fields::Se, fields::Se>,
         boost::mpl::vector<fields::Chi, typename cxx_diagrams::bar<fields::Chi>::type, fields::Se, fields::Se>,
         boost::mpl::vector<fields::Chi, fields::Chi, fields::Se, fields::Se>
         >;
      using combinatorial_factors = boost::mpl::vector<
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>
         >;
      using color_factors = boost::mpl::vector<
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>
         >;
      return accumulate_generic<
         GenericKeys,
         GenericInsertions,
         combinatorial_factors,
         color_factors,
         boost::mpl::int_<10>,
         genericSum1_impl
         >( *this );
   } // End of function genericSum1()
   
   template<class GenericFieldMap>
   struct genericSum2_impl : generic_sum_base {
      genericSum2_impl( const generic_sum_base &base ) :
      generic_sum_base( base ) {
      } // End of constructor genericSum2_impl
   
      std::array<std::complex<double>,10> operator()( void ) {
         using boost::mpl::at;
         using boost::fusion::at_key;
         using gF5 = typename at<GenericFieldMap,GenericF5Key>::type;
         using gF8 = typename at<GenericFieldMap,GenericF8Key>::type;
         using gS6 = typename at<GenericFieldMap,GenericS6Key>::type;
         using gS7 = typename at<GenericFieldMap,GenericS7Key>::type;
   
         std::array<int, 1> i1 {this->external_indices(0)};
         std::array<int, 1> i2 {this->external_indices(1)};
         std::array<int, 1> i3 {this->external_indices(2)};
         std::array<int, 1> i4 {this->external_indices(3)};
         
         typename field_index_map<GenericFieldMap>::type index_map;
         const context_with_vertices &context = *this;
         std::complex<double> S_LL = 0.0;
         std::complex<double> S_LR = 0.0;
         std::complex<double> S_RL = 0.0;
         std::complex<double> S_RR = 0.0;
         std::complex<double> V_LL = 0.0;
         std::complex<double> V_LR = 0.0;
         std::complex<double> V_RL = 0.0;
         std::complex<double> V_RR = 0.0;
         std::complex<double> T_LL = 0.0;
         std::complex<double> T_RR = 0.0;
   
         // Shorter aliases for large types
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _gF5 = typename cxx_diagrams::bar<gF5>::type;
         using _gF8 = typename cxx_diagrams::bar<gF8>::type;
         using _gS6 = typename cxx_diagrams::conj<gS6>::type;
         using _gS7 = typename cxx_diagrams::conj<gS7>::type;
         using Fe = fields::Fe;
         using Fe = fields::Fe;
         
         // Aliases for loop function identifiers
         using looplibrary::dd0;
         using looplibrary::dd00;
         
         double mF5, mF8, mS6, mS7;
         std::complex<double> g1, g2, g3, g4, g5, g6, g7, g8;
         std::array<int, 8> z{};
         looplibrary::Dcoeff_t l1 {};
         // Start of summation over generic fields.
         for( const auto &iF5 : index_range<gF5>() ) {
         at_key<GenericF5Key>( index_map ) = iF5;
         for( const auto &iF8 : index_range<gF8>() ) {
         at_key<GenericF8Key>( index_map ) = iF8;
         for( const auto &iS6 : index_range<gS6>() ) {
         at_key<GenericS6Key>( index_map ) = iS6;
         for( const auto &iS7 : index_range<gS7>() ) {
         at_key<GenericS7Key>( index_map ) = iS7;
         
            mF5 = context.mass<gF5>(iF5);
            mF8 = context.mass<gF8>(iF8);
            mS6 = context.mass<gS6>(iS6);
            mS7 = context.mass<gS7>(iS7);
            
            g1 = NPF_L(_Fe, gF8, gS6) NPF_I(i4, iF8, iS6);
            g2 = NPF_R(_Fe, gF8, gS6) NPF_I(i4, iF8, iS6);
            g3 = NPF_L(_gF5, Fe, _gS6) NPF_I(iF5, i1, iS6);
            g4 = NPF_R(_gF5, Fe, _gS6) NPF_I(iF5, i1, iS6);
            g5 = NPF_L(_gF8, _Fe, gS7) NPF_I(iF8, i3, iS7);
            g6 = NPF_R(_gF8, _Fe, gS7) NPF_I(iF8, i3, iS7);
            g7 = NPF_L(Fe, gF5, _gS7) NPF_I(i2, iF5, iS7);
            g8 = NPF_R(Fe, gF5, _gS7) NPF_I(i2, iF5, iS7);
            
            z[0] = (std::abs(g1) < std::numeric_limits<double>::epsilon())?0:1;
            z[1] = (std::abs(g2) < std::numeric_limits<double>::epsilon())?0:1;
            z[2] = (std::abs(g3) < std::numeric_limits<double>::epsilon())?0:1;
            z[3] = (std::abs(g4) < std::numeric_limits<double>::epsilon())?0:1;
            z[4] = (std::abs(g5) < std::numeric_limits<double>::epsilon())?0:1;
            z[5] = (std::abs(g6) < std::numeric_limits<double>::epsilon())?0:1;
            z[6] = (std::abs(g7) < std::numeric_limits<double>::epsilon())?0:1;
            z[7] = (std::abs(g8) < std::numeric_limits<double>::epsilon())?0:1;
            
            if( z[0]*(z[2]*z[4]*z[6] + z[3]*z[5]*z[6] + z[3]*z[4]*z[7] + z[2]*z[5]*z[7]) + z[1]*(z[3]*z[4]*z[6] + z[2]*z[5]*z[6] + z[2]*z[4]*z[7] + z[3]*z[5]*z[7]) == 0 ) continue;
            Loop_library::get().D(l1,0,0,0,0,0,0,Sqr(mF8),Sqr(mF5),Sqr(mS7),Sqr(mS6),Sqr(context.scale()));
         
            S_LL += 0.0031662869888230555*g1*g3*g5*g7*l1[dd0]*mF5*mF8;
            S_LR += -0.012665147955292222*g2*g3*g5*g8*l1[dd00];
            S_RL += -0.012665147955292222*g1*g4*g6*g7*l1[dd00];
            S_RR += 0.0031662869888230555*g2*g4*g6*g8*l1[dd0]*mF5*mF8;
            V_LL += -0.0031662869888230555*g2*g3*g6*g7*l1[dd0]*mF5*mF8;
            V_LR += -0.006332573977646111*g1*g3*g6*g8*l1[dd00];
            V_RL += -0.006332573977646111*g2*g4*g5*g7*l1[dd00];
            V_RR += -0.0031662869888230555*g1*g4*g5*g8*l1[dd0]*mF5*mF8;
            T_LL += 0.0007915717472057639*g1*g3*g5*g7*l1[dd0]*mF5*mF8;
            T_RR += 0.0007915717472057639*g2*g4*g6*g8*l1[dd0]*mF5*mF8;
         
         }}}} // End of summation over generic fields
   
         return {S_LL, S_LR, S_RL, S_RR, V_LL, V_LR, V_RL, V_RR, T_LL, T_RR};
      } // End of operator()( void )
   }; // End of struct genericSum2_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum2( void ) {
      using GenericKeys = boost::mpl::vector< GenericF5Key, GenericF8Key, GenericS6Key, GenericS7Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Chi>::type, typename cxx_diagrams::bar<fields::Chi>::type, fields::Se, fields::Se>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Chi>::type, fields::Chi, fields::Se, fields::Se>,
         boost::mpl::vector<fields::Chi, typename cxx_diagrams::bar<fields::Chi>::type, fields::Se, fields::Se>,
         boost::mpl::vector<fields::Chi, fields::Chi, fields::Se, fields::Se>
         >;
      using combinatorial_factors = boost::mpl::vector<
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>
         >;
      using color_factors = boost::mpl::vector<
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>
         >;
      return accumulate_generic<
         GenericKeys,
         GenericInsertions,
         combinatorial_factors,
         color_factors,
         boost::mpl::int_<10>,
         genericSum2_impl
         >( *this );
   } // End of function genericSum2()
   
   template<class GenericFieldMap>
   struct genericSum3_impl : generic_sum_base {
      genericSum3_impl( const generic_sum_base &base ) :
      generic_sum_base( base ) {
      } // End of constructor genericSum3_impl
   
      std::array<std::complex<double>,10> operator()( void ) {
         using boost::mpl::at;
         using boost::fusion::at_key;
         using gF5 = typename at<GenericFieldMap,GenericF5Key>::type;
         using gF8 = typename at<GenericFieldMap,GenericF8Key>::type;
         using gS6 = typename at<GenericFieldMap,GenericS6Key>::type;
         using gS7 = typename at<GenericFieldMap,GenericS7Key>::type;
   
         std::array<int, 1> i1 {this->external_indices(0)};
         std::array<int, 1> i2 {this->external_indices(1)};
         std::array<int, 1> i3 {this->external_indices(2)};
         std::array<int, 1> i4 {this->external_indices(3)};
         
         typename field_index_map<GenericFieldMap>::type index_map;
         const context_with_vertices &context = *this;
         std::complex<double> S_LL = 0.0;
         std::complex<double> S_LR = 0.0;
         std::complex<double> S_RL = 0.0;
         std::complex<double> S_RR = 0.0;
         std::complex<double> V_LL = 0.0;
         std::complex<double> V_LR = 0.0;
         std::complex<double> V_RL = 0.0;
         std::complex<double> V_RR = 0.0;
         std::complex<double> T_LL = 0.0;
         std::complex<double> T_RR = 0.0;
   
         // Shorter aliases for large types
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _gF5 = typename cxx_diagrams::bar<gF5>::type;
         using _gF8 = typename cxx_diagrams::bar<gF8>::type;
         using _gS6 = typename cxx_diagrams::conj<gS6>::type;
         using _gS7 = typename cxx_diagrams::conj<gS7>::type;
         using Fe = fields::Fe;
         using Fe = fields::Fe;
         
         // Aliases for loop function identifiers
         using looplibrary::dd0;
         using looplibrary::dd00;
         
         double mF5, mF8, mS6, mS7;
         std::complex<double> g1, g2, g3, g4, g5, g6, g7, g8;
         std::array<int, 8> z{};
         looplibrary::Dcoeff_t l1 {};
         // Start of summation over generic fields.
         for( const auto &iF5 : index_range<gF5>() ) {
         at_key<GenericF5Key>( index_map ) = iF5;
         for( const auto &iF8 : index_range<gF8>() ) {
         at_key<GenericF8Key>( index_map ) = iF8;
         for( const auto &iS6 : index_range<gS6>() ) {
         at_key<GenericS6Key>( index_map ) = iS6;
         for( const auto &iS7 : index_range<gS7>() ) {
         at_key<GenericS7Key>( index_map ) = iS7;
         
            mF5 = context.mass<gF5>(iF5);
            mF8 = context.mass<gF8>(iF8);
            mS6 = context.mass<gS6>(iS6);
            mS7 = context.mass<gS7>(iS7);
            
            g1 = NPF_L(_Fe, gF5, gS7) NPF_I(i3, iF5, iS7);
            g2 = NPF_R(_Fe, gF5, gS7) NPF_I(i3, iF5, iS7);
            g3 = NPF_L(_Fe, gF8, gS6) NPF_I(i4, iF8, iS6);
            g4 = NPF_R(_Fe, gF8, gS6) NPF_I(i4, iF8, iS6);
            g5 = NPF_L(_gF5, Fe, _gS6) NPF_I(iF5, i1, iS6);
            g6 = NPF_R(_gF5, Fe, _gS6) NPF_I(iF5, i1, iS6);
            g7 = NPF_L(_gF8, Fe, _gS7) NPF_I(iF8, i2, iS7);
            g8 = NPF_R(_gF8, Fe, _gS7) NPF_I(iF8, i2, iS7);
            
            z[0] = (std::abs(g1) < std::numeric_limits<double>::epsilon())?0:1;
            z[1] = (std::abs(g2) < std::numeric_limits<double>::epsilon())?0:1;
            z[2] = (std::abs(g3) < std::numeric_limits<double>::epsilon())?0:1;
            z[3] = (std::abs(g4) < std::numeric_limits<double>::epsilon())?0:1;
            z[4] = (std::abs(g5) < std::numeric_limits<double>::epsilon())?0:1;
            z[5] = (std::abs(g6) < std::numeric_limits<double>::epsilon())?0:1;
            z[6] = (std::abs(g7) < std::numeric_limits<double>::epsilon())?0:1;
            z[7] = (std::abs(g8) < std::numeric_limits<double>::epsilon())?0:1;
            
            if( z[0]*(z[2]*z[4]*z[6] + z[3]*z[5]*z[6] + z[3]*z[4]*z[7] + z[2]*z[5]*z[7]) + z[1]*(z[3]*z[4]*z[6] + z[2]*z[5]*z[6] + z[2]*z[4]*z[7] + z[3]*z[5]*z[7]) == 0 ) continue;
            Loop_library::get().D(l1,0,0,0,0,0,0,Sqr(mF8),Sqr(mF5),Sqr(mS7),Sqr(mS6),Sqr(context.scale()));
         
            S_LL += -0.006332573977646111*g1*g3*g5*g7*l1[dd0]*mF5*mF8;
            S_LR += -0.006332573977646111*g1*g4*g5*g8*l1[dd0]*mF5*mF8;
            S_RL += -0.006332573977646111*g2*g3*g6*g7*l1[dd0]*mF5*mF8;
            S_RR += -0.006332573977646111*g2*g4*g6*g8*l1[dd0]*mF5*mF8;
            V_LL += -0.006332573977646111*g2*g4*g5*g7*l1[dd00];
            V_LR += -0.006332573977646111*g2*g3*g5*g8*l1[dd00];
            V_RL += -0.006332573977646111*g1*g4*g6*g7*l1[dd00];
            V_RR += -0.006332573977646111*g1*g3*g6*g8*l1[dd00];
            T_LL += 0;
            T_RR += 0;
         
         }}}} // End of summation over generic fields
   
         return {S_LL, S_LR, S_RL, S_RR, V_LL, V_LR, V_RL, V_RR, T_LL, T_RR};
      } // End of operator()( void )
   }; // End of struct genericSum3_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum3( void ) {
      using GenericKeys = boost::mpl::vector< GenericF5Key, GenericF8Key, GenericS6Key, GenericS7Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<fields::Fv, fields::Fv, fields::Hpm, fields::Hpm>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Cha1>::type, typename cxx_diagrams::bar<fields::Cha1>::type, fields::Sv, fields::Sv>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Chi>::type, typename cxx_diagrams::bar<fields::Chi>::type, fields::Se, fields::Se>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Chi>::type, fields::Chi, fields::Se, fields::Se>,
         boost::mpl::vector<fields::Chi, typename cxx_diagrams::bar<fields::Chi>::type, fields::Se, fields::Se>,
         boost::mpl::vector<fields::Chi, fields::Chi, fields::Se, fields::Se>
         >;
      using combinatorial_factors = boost::mpl::vector<
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>
         >;
      using color_factors = boost::mpl::vector<
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>
         >;
      return accumulate_generic<
         GenericKeys,
         GenericInsertions,
         combinatorial_factors,
         color_factors,
         boost::mpl::int_<10>,
         genericSum3_impl
         >( *this );
   } // End of function genericSum3()
   
   template<class GenericFieldMap>
   struct genericSum4_impl : generic_sum_base {
      genericSum4_impl( const generic_sum_base &base ) :
      generic_sum_base( base ) {
      } // End of constructor genericSum4_impl
   
      std::array<std::complex<double>,10> operator()( void ) {
         using boost::mpl::at;
         using boost::fusion::at_key;
         using gF6 = typename at<GenericFieldMap,GenericF6Key>::type;
         using gF7 = typename at<GenericFieldMap,GenericF7Key>::type;
         using gS5 = typename at<GenericFieldMap,GenericS5Key>::type;
         using gS8 = typename at<GenericFieldMap,GenericS8Key>::type;
   
         std::array<int, 1> i1 {this->external_indices(0)};
         std::array<int, 1> i2 {this->external_indices(1)};
         std::array<int, 1> i3 {this->external_indices(2)};
         std::array<int, 1> i4 {this->external_indices(3)};
         
         typename field_index_map<GenericFieldMap>::type index_map;
         const context_with_vertices &context = *this;
         std::complex<double> S_LL = 0.0;
         std::complex<double> S_LR = 0.0;
         std::complex<double> S_RL = 0.0;
         std::complex<double> S_RR = 0.0;
         std::complex<double> V_LL = 0.0;
         std::complex<double> V_LR = 0.0;
         std::complex<double> V_RL = 0.0;
         std::complex<double> V_RR = 0.0;
         std::complex<double> T_LL = 0.0;
         std::complex<double> T_RR = 0.0;
   
         // Shorter aliases for large types
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _gF6 = typename cxx_diagrams::bar<gF6>::type;
         using _gF7 = typename cxx_diagrams::bar<gF7>::type;
         using _gS5 = typename cxx_diagrams::conj<gS5>::type;
         using _gS8 = typename cxx_diagrams::conj<gS8>::type;
         using Fe = fields::Fe;
         using Fe = fields::Fe;
         
         // Aliases for loop function identifiers
         using looplibrary::dd0;
         using looplibrary::dd00;
         
         double mF6, mF7, mS5, mS8;
         std::complex<double> g1, g2, g3, g4, g5, g6, g7, g8;
         std::array<int, 8> z{};
         looplibrary::Dcoeff_t l1 {};
         // Start of summation over generic fields.
         for( const auto &iF6 : index_range<gF6>() ) {
         at_key<GenericF6Key>( index_map ) = iF6;
         for( const auto &iF7 : index_range<gF7>() ) {
         at_key<GenericF7Key>( index_map ) = iF7;
         for( const auto &iS5 : index_range<gS5>() ) {
         at_key<GenericS5Key>( index_map ) = iS5;
         for( const auto &iS8 : index_range<gS8>() ) {
         at_key<GenericS8Key>( index_map ) = iS8;
         
            mF6 = context.mass<gF6>(iF6);
            mF7 = context.mass<gF7>(iF7);
            mS5 = context.mass<gS5>(iS5);
            mS8 = context.mass<gS8>(iS8);
            
            g1 = NPF_L(_Fe, gF7, gS5) NPF_I(i3, iF7, iS5);
            g2 = NPF_R(_Fe, gF7, gS5) NPF_I(i3, iF7, iS5);
            g3 = NPF_L(_Fe, gF6, gS8) NPF_I(i4, iF6, iS8);
            g4 = NPF_R(_Fe, gF6, gS8) NPF_I(i4, iF6, iS8);
            g5 = NPF_L(_gF6, Fe, _gS5) NPF_I(iF6, i1, iS5);
            g6 = NPF_R(_gF6, Fe, _gS5) NPF_I(iF6, i1, iS5);
            g7 = NPF_L(_gF7, Fe, _gS8) NPF_I(iF7, i2, iS8);
            g8 = NPF_R(_gF7, Fe, _gS8) NPF_I(iF7, i2, iS8);
            
            z[0] = (std::abs(g1) < std::numeric_limits<double>::epsilon())?0:1;
            z[1] = (std::abs(g2) < std::numeric_limits<double>::epsilon())?0:1;
            z[2] = (std::abs(g3) < std::numeric_limits<double>::epsilon())?0:1;
            z[3] = (std::abs(g4) < std::numeric_limits<double>::epsilon())?0:1;
            z[4] = (std::abs(g5) < std::numeric_limits<double>::epsilon())?0:1;
            z[5] = (std::abs(g6) < std::numeric_limits<double>::epsilon())?0:1;
            z[6] = (std::abs(g7) < std::numeric_limits<double>::epsilon())?0:1;
            z[7] = (std::abs(g8) < std::numeric_limits<double>::epsilon())?0:1;
            
            if( z[0]*(z[2]*z[4]*z[6] + z[3]*z[5]*z[6] + z[3]*z[4]*z[7] + z[2]*z[5]*z[7]) + z[1]*(z[3]*z[4]*z[6] + z[2]*z[5]*z[6] + z[2]*z[4]*z[7] + z[3]*z[5]*z[7]) == 0 ) continue;
            Loop_library::get().D(l1,0,0,0,0,0,0,Sqr(mF7),Sqr(mF6),Sqr(mS8),Sqr(mS5),Sqr(context.scale()));
         
            S_LL += 0.0031662869888230555*g1*g3*g5*g7*l1[dd0]*mF6*mF7;
            S_LR += 0.012665147955292222*g1*g4*g5*g8*l1[dd00];
            S_RL += 0.012665147955292222*g2*g3*g6*g7*l1[dd00];
            S_RR += 0.0031662869888230555*g2*g4*g6*g8*l1[dd0]*mF6*mF7;
            V_LL += -0.006332573977646111*g2*g4*g5*g7*l1[dd00];
            V_LR += 0.0031662869888230555*g2*g3*g5*g8*l1[dd0]*mF6*mF7;
            V_RL += 0.0031662869888230555*g1*g4*g6*g7*l1[dd0]*mF6*mF7;
            V_RR += -0.006332573977646111*g1*g3*g6*g8*l1[dd00];
            T_LL += -0.0007915717472057639*g1*g3*g5*g7*l1[dd0]*mF6*mF7;
            T_RR += -0.0007915717472057639*g2*g4*g6*g8*l1[dd0]*mF6*mF7;
         
         }}}} // End of summation over generic fields
   
         return {S_LL, S_LR, S_RL, S_RR, V_LL, V_LR, V_RL, V_RR, T_LL, T_RR};
      } // End of operator()( void )
   }; // End of struct genericSum4_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum4( void ) {
      using GenericKeys = boost::mpl::vector< GenericF6Key, GenericF7Key, GenericS5Key, GenericS8Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<fields::Fv, fields::Fv, fields::Hpm, fields::Hpm>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Cha1>::type, typename cxx_diagrams::bar<fields::Cha1>::type, fields::Sv, fields::Sv>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Chi>::type, typename cxx_diagrams::bar<fields::Chi>::type, fields::Se, fields::Se>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Chi>::type, fields::Chi, fields::Se, fields::Se>,
         boost::mpl::vector<fields::Chi, typename cxx_diagrams::bar<fields::Chi>::type, fields::Se, fields::Se>,
         boost::mpl::vector<fields::Chi, fields::Chi, fields::Se, fields::Se>
         >;
      using combinatorial_factors = boost::mpl::vector<
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>
         >;
      using color_factors = boost::mpl::vector<
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>
         >;
      return accumulate_generic<
         GenericKeys,
         GenericInsertions,
         combinatorial_factors,
         color_factors,
         boost::mpl::int_<10>,
         genericSum4_impl
         >( *this );
   } // End of function genericSum4()

   public:
   nPointFeFeFeFe_3917882127499752( const MRSSM2_mass_eigenstates &model,
   const std::array<int, 4> &indices,
   const std::array<Eigen::Vector4d, 0> &momenta
    ) :
   correlation_function_context<4,0> { model, indices, momenta }
   {}

   std::array<std::complex<double>,10> calculate( void ) {
      std::array<std::complex<double>,10> genericSummation;
      constexpr int coeffsLength = genericSummation.size();
      const auto genericsum1 = genericSum1();
      const auto genericsum2 = genericSum2();
      const auto genericsum3 = genericSum3();
      const auto genericsum4 = genericSum4();
   
      for ( std::size_t i=0; i<coeffsLength; i++ ) {
         genericSummation[i] += genericsum1[i]+genericsum2[i]+genericsum3[i]+genericsum4[i];
      }
      return genericSummation;
   } // End of calculate()
}; // End of nPointFeFeFeFe_3917882127499752

std::array<std::complex<double>, 10> boxes1loop(const MRSSM2_mass_eigenstates &model,
const std::array<int, 4> &indices,
const std::array<Eigen::Vector4d, 0> &momenta
) {
   nPointFeFeFeFe_3917882127499752 helper{model, indices, momenta};
   return helper.calculate();
}

} // namespace npointfunctions
} // namespace MRSSM2_cxx_diagrams

namespace {

std::valarray<std::complex<double>> zero(
   int generationIndex1,
   int generationIndex2,
   const MRSSM2_mass_eigenstates&,
   bool){
   std::valarray<std::complex<double>> res {0.0, 0.0, 0.0, 0.0};
   return res;
}

std::array<std::complex<double>,10> zero(
   const MRSSM2_mass_eigenstates&,
   const std::array<int,4>&,
   const std::array<Eigen::Vector4d,0>&){
   std::array<std::complex<double>,10> res{};
   return res;
}

double get_MSUSY(const MRSSM2_mass_eigenstates_interface& model)
{
   return Min(model.get_MSRdp(), model.get_MSRum(), model.get_MSd().tail<6>().
      minCoeff(), model.get_MSu().tail<6>().minCoeff(), model.get_MSe().tail<6>().
      minCoeff(), model.get_MHpm().tail<3>().minCoeff(), model.get_MCha1().tail<2>
      ().minCoeff(), model.get_MCha2().tail<2>().minCoeff());
}

void run_to_MSUSY(MRSSM2_mass_eigenstates& model)
{
   const double precision_goal = model.get_precision();
   const int Nmax = static_cast<int>(-std::log10(precision_goal) * 10);
   int it = 0;
   double precision = 0.;
   double MSUSY_old = 0., MSUSY_new = 0.;

   VERBOSE_MSG("MRSSM2_br_l_to_3l: iterate to run to MSUSY ...");

   do {
      MSUSY_new = get_MSUSY(model);

      if (MSUSY_new <= 0.)
         return;

      VERBOSE_MSG("MRSSM2_br_l_to_3l:    running to new MSUSY = " << MSUSY_new);

      model.run_to(MSUSY_new);
      model.calculate_DRbar_masses();

      precision = MaxRelDiff(MSUSY_old, MSUSY_new);
      MSUSY_old = MSUSY_new;
   } while (precision > precision_goal && ++it < Nmax);

   VERBOSE_MSG("MRSSM2_br_l_to_3l: iteration finished. MSUSY = " << MSUSY_new);

   if (precision > precision_goal) {
      WARNING("MRSSM2_br_l_to_3l: run_to_MSUSY(): "
              "Iteration did not converge after " << Nmax
              << " iterations (precision goal: " << precision_goal << ")");
   }
}

} // anonymous namespace

using namespace MRSSM2_cxx_diagrams;
using namespace MRSSM2_FFV_form_factors;

namespace MRSSM2_br_l_to_3l {

typedef std::valarray<std::complex<double>> (*ffv)(
   int, int, const MRSSM2_mass_eigenstates&, bool);

typedef std::array<std::complex<double>,10> (*npf)(
   const MRSSM2_mass_eigenstates&,
   const std::array<int,4>&,
   const std::array<Eigen::Vector4d,0>&);

/**
 * @param[in] nI     Generation index of incoming lepton.
 * @return Total decay width in GeV, according to PDG data.
 */
double get_total_width(int nI) {
   // https://pdg.lbl.gov/2020/reviews/rpp2020-rev-phys-constants.pdf
   constexpr double hbar = 6.582119569e-25, // [GeV*s]
   // https://pdg.lbl.gov/2020/listings/rpp2020-list-muon.pdf
                    muon = 2.1969811e-6,    // [s]
   // https://pdg.lbl.gov/2020/listings/rpp2020-list-tau.pdf
                    tau  = 290.3e-15;       // [s]
   switch (nI) {
      case 1: return hbar/muon;
      case 2: return hbar/tau;
      default: throw std::invalid_argument("Unrecognized lepton");
   }
}

/**
 * @tparam    Lepton Type of a lepton field.
 * @tparam    Photon Type of a photon field.
 * @param[in] model  Mass eigenstates.
 * @param[in] nA     Generation index.
 * @return Left part of ffv coupling, multiplied by -i.
 */
template <class Lepton, class Photon>
std::complex<double> left(const MRSSM2_mass_eigenstates& model, int nA) {
    context_base context {model};
    using vertex = Vertex<typename Lepton::lorentz_conjugate, Lepton, Photon>;
    std::array<int, 2> indices {nA, nA};
    const auto value =  vertex::evaluate(indices, context);
    return value.left();
}

/**
 * @tparam    Lepton Type of a lepton field.
 * @tparam    Photon Type of a photon field.
 * @param[in] model  Mass eigenstates.
 * @param[in] nA     Generation index.
 * @return Right part of ffv coupling, multiplied by -i.
 */
template <class Lepton, class Photon>
std::complex<double> right(const MRSSM2_mass_eigenstates& model, int nA) {
    context_base context {model};
    using vertex = Vertex<typename Lepton::lorentz_conjugate, Lepton, Photon>;
    std::array<int, 2> indices {nA, nA};
    const auto value =  vertex::evaluate(indices, context);
    return value.right();
}

/**
 * @tparam    L     Type of a lepton field.
 * @tparam    A     Type of a photon field.
 * @tparam    T     Type of a form factors.
 * @param[in] model Mass eigenstates.
 * @param[in] ff    Lepton-photon form factors.
 * @param[in] g     Generation index for leptons.
 * @return Set of four-fermion coefficients from photon penguin amplitudes
 *         (without overall i; with appropriate embedding).
 */
template <class L, class A, class T>
Eigen::Array<std::complex<double>,10,1> embed_photon(
   const MRSSM2_mass_eigenstates& model, const T& ff, int g) {
   // Get quark-photon couplings (without i, as everywhere):
   const auto lL =  left<L, A>(model, g);
   const auto lR = right<L, A>(model, g);
   Eigen::Array<std::complex<double>,10,1> res{};
   // Minus comes from the form_factors embedding into four-fermion amplitude:
   res[4] = - ff[0] * lL;
   res[5] = - ff[0] * lR;
   res[6] = - ff[1] * lL;
   res[7] = - ff[1] * lR;
   return res;
};

/**
 * @tparam    Name  Function name for npf function.
 * @param[in] model Mass eigenstates.
 * @param[in] nI    Generation index for incoming lepton.
 * @param[in] nO    Generation index for outgoing lepton.
 * @param[in] nA    Generation index for lepton pair.
 * @return Set of four-fermion coefficients for non-photonic amplitudes
 *         (without overall i; with fixed signs for tensor operators).
 */
template <npf Name>
Eigen::Array<std::complex<double>,10,1> fix_tensors_sign(
   const MRSSM2_mass_eigenstates& model, int nI, int nO, int nA) {
   const auto npf = Name(model,
         std::array<int,4>{nI, nA, nO, nA},
         std::array<Eigen::Vector4d, 0>{});
   Eigen::Array<std::complex<double>,10,1> res(npf.data());
   res[8] = - res[8];
   res[9] = - res[9];
   return res;
};

/**
 * @tparam    A      Type of input sets of coefficients.
 * @param[in] photon T-coefficients of all penguins.
 * @param[in] res    T-coefficients of all penguins.
 * @param[in] boxes  Coefficients of all boxes.
 * @return Wilson coefficients of four fermion operators.
 */
template<class A>
Eigen::Array<std::complex<double>,10,1> fierz(A photon, A rest, A box) {
   A t_channel = photon + rest;
   A tu_channels;
   tu_channels << t_channel[0]/2 - 6*t_channel[8],
                  t_channel[1] - 2*t_channel[5],
                  t_channel[2] - 2*t_channel[6],
                  t_channel[3]/2 - 6*t_channel[9],
                  2*t_channel[4],
                  t_channel[5] - t_channel[1]/2,
                  t_channel[6] - t_channel[2]/2,
                  2*t_channel[7],
                  1.5*t_channel[8] - t_channel[0]/8,
                  1.5*t_channel[9] - t_channel[3]/8;
   tu_channels = tu_channels + box;
   A res{};
   res << 2*tu_channels[0],
          0,
          0,
          2*tu_channels[3],
          tu_channels[4]/2,
          tu_channels[5],
          tu_channels[6],
          tu_channels[7]/2,
          0,
          0;
   return res;
}

struct LeptonOperators {
   std::complex<double> SLL, SLR, SRL, SRR,
                        VLL, VLR, VRL, VRR,
                        TLL, TRR;
   template<class T>
   LeptonOperators(T coeffs) {
      SLL = coeffs[0];
      SLR = coeffs[1];
      SRL = coeffs[2];
      SRR = coeffs[3];
      VLL = coeffs[4];
      VLR = coeffs[5];
      VRL = coeffs[6];
      VRR = coeffs[7];
      TLL = coeffs[8];
      TRR = coeffs[9];
   }
};

/**
 * @param[in] g Generation index for leptons.
 * @param[in] qedqcd Reference to low-energy data.
 * @return A mass of a lepton for a given generation number.
 */
double get_pole_mass(int g, const softsusy::QedQcd& qedqcd) {
   switch (g) {
      case 0: return qedqcd.displayMass(softsusy::mElectron);
      case 1: return qedqcd.displayMass(softsusy::mMuon);
      case 2: return qedqcd.displayMass(softsusy::mTau);
      default: throw std::invalid_argument("Unrecognized lepton");
   }
}

/**
 * @tparam    T1     Type of a dipole coefficient.
 * @tparam    T2     Type of a four-lepton coefficients.
 * @param[in] nI     Generation of decaying lepton.
 * @param[in] nA     Generation of any lepton in pair.
 * @param[in] model  Mass eigenstates.
 * @param[in] qedqcd Reference to low-energy data.
 * @return Width of decay with leptons of the same generation in the end.
 */
template<class T1, class T2>
double width_same(int nI, int nA, T1 DL, T1 DR, T2 coeffs,
   const MRSSM2_mass_eigenstates& model, const softsusy::QedQcd& qedqcd) {

   context_base context {model};
   const double e = unit_charge(context);
   const auto mass_heavy = get_pole_mass(nI, qedqcd);
   const auto mass_light = get_pole_mass(nA, qedqcd);
   const auto r = mass_heavy / mass_light;
   const LeptonOperators C(coeffs);

   using std::norm;
   using std::real;
   using std::conj;
   using std::log;

   // 1802.06803: wrong eq. (3.12) for X term;
   // 1702.03020: wrong X term for their covariant derivative close to eq. (A.1).
   // hep-ph/0004025, hep-ph/9909265: same;
   double res = e*e * (norm(DL) + norm(DR)) * (-11. + 8.*log(r)) / 192.;
   res += e * (2.*real(C.VLL*conj(DR)) + real(C.VLR*conj(DR))) / 192.;
   res += e * (2.*real(C.VRR*conj(DL)) + real(C.VRL*conj(DL))) / 192.;
   res += (norm(C.SLL) + 16.*norm(C.VLL) + 8.*norm(C.VLR)) / 12288.;
   res += (norm(C.SRR) + 16.*norm(C.VRR) + 8.*norm(C.VRL)) / 12288.;
   return res * pow(mass_heavy, 5) / pow(Pi, 3);
}

/**
 * @tparam    T1     Type of a dipole coefficient.
 * @tparam    T2     Type of a four-lepton coefficients.
 * @param[in] nI     Generation of decaying lepton.
 * @param[in] nA     Generation of any lepton in pair.
 * @param[in] model  Mass eigenstates.
 * @param[in] qedqcd Reference to low-energy data.
 * @return Width of decay with leptons of different generations in the end.
 */
template<class T1, class T2>
double width_diff(int nI, int nA, T1 DL, T1 DR, T2 coeffs,
   const MRSSM2_mass_eigenstates& model, const softsusy::QedQcd& qedqcd) {

   context_base context {model};
   const double e = unit_charge(context);
   const auto mass_heavy = get_pole_mass(nI, qedqcd);
   const auto mass_light = get_pole_mass(nA, qedqcd);
   const auto r = mass_heavy / mass_light;
   const LeptonOperators C(coeffs);

   using std::norm;
   using std::real;
   using std::conj;
   using std::log;

   // 1802.06803: wrong eq. (3.12) - they take heaviest mass of final states.
   //             correct eq. (3.11);
   // 1212.5939: wrong r in eq. (3.10) see eq. (4.9) of hep-ph/9403398, where
   //            text explanation is given;
   double res = e*e * (norm(DL) + norm(DR)) * (-3. + 2.*log(r)) / 48.;
   res += e * (real(C.VLL*conj(DR)) + real(C.VLR*conj(DR))) / 192.;
   res += e * (real(C.VRR*conj(DL)) + real(C.VRL*conj(DL))) / 192.;
   res += (norm(C.SLL) + norm(C.SLR) + norm(C.SRL) + norm(C.SRR)) / 6144.;
   res += (norm(C.VLL) + norm(C.VLR) + norm(C.VRL) + norm(C.VRR)) / 1536.;
   res += (norm(C.TLL) + norm(C.TRR)) / 128.;
   return res * pow(mass_heavy, 5) / pow(Pi, 3);
}

/**
 * @tparam    Lepton  Type of a lepton field.
 * @tparam    Photon  Type of a photon field.
 * @tparam    Factor  Function name for photon t-penguins.
 * @tparam    Scalars Function name for scalar t-penguins.
 * @tparam    Vectors Function name for vector t-penguins.
 * @tparam    Boxes   Function name for all box diagrams.
 * @param[in] nI      Generation of decaying lepton.
 * @param[in] nO      Generation of lepton, which can be separated from pair.
 * @param[in] nA      Generation of lepton and antilepton.
 * @param[in] model   Mass eigenstates.
 * @param[in] qedqcd  Reference to low-energy data.
 * @return Observable value and Wilson coefficients used to derive it.
 */
template<class Lepton, class Photon, ffv Factor, npf Scalars, npf Vectors, npf Boxes>
Eigen::Array<std::complex<double>,13,1> forge(int nI, int nO, int nA,
   const MRSSM2_mass_eigenstates& model_, const softsusy::QedQcd& qedqcd) {

   MRSSM2_mass_eigenstates model(model_);
   try {
      run_to_MSUSY(model);
      model.get_physical().clear();
      model.solve_ewsb();
   } catch (const Error& e) {
      ERROR("MRSSM2_br_l_to_3l:" << e.what_detailed());
      return std::numeric_limits<Eigen::Array<std::complex<double>,13,1>>::quiet_NaN();
   }

   context_base context {model};

   const auto form_factors = Factor(nI, nO, model, false);
   const auto photon_amp = F4F(Lepton, nA);

   // Full amplitude is calculated with the following convention (for the case
   // 4!=3, otherwise Fierz transformations should be used hep-ph/0412245, note
   // the half for every sigma-sigma summation in eqs. 40 and 42):
   //    <out:4,3| T exp(i L_full dx) |in:2,1> = i * npf,
   // which means, that it has to be matched with -1 * C_XY, because
   //    <out:4,3| T exp(i L_low dx) |in:2,1> = -i * C_XY [3 X 1]*[4 Y 2].
   auto t_amp = NPF(Scalars);
   t_amp = t_amp + NPF(Vectors);
   const auto box_amp = NPF(Boxes);

   // If final leptons are the same, then we need to add u-penguins as well. If
   // not, then we can just sum them (we always calculate all boxes in npf).
   auto CXYl = (nO == nA) ? fierz(photon_amp, t_amp, box_amp)
                          : photon_amp + t_amp + box_amp;
   // Matching
   const auto DL = - 0.5 * form_factors[2];
   const auto DR = - 0.5 * form_factors[3];
   CXYl = - CXYl;

   // @todo Running of coefficients

   double partial_width;
   if (nO == nA) partial_width = width_same(nI, nA, DL, DR, CXYl, model, qedqcd);
   else          partial_width = width_diff(nI, nA, DL, DR, CXYl, model, qedqcd);

   const double total_width = get_total_width(nI);
   Eigen::Array<std::complex<double>,13,1> res;

   res << partial_width/total_width,
          DL,
          DR,
          CXYl[0],
          CXYl[1],
          CXYl[2],
          CXYl[3],
          CXYl[4],
          CXYl[5],
          CXYl[6],
          CXYl[7],
          CXYl[8],
          CXYl[9];

   return res;
}

Eigen::Array<std::complex<double>,13,1> calculate_Fe_to_FeFebarFe_for_All_1loop(int nI, int nO, int nA, const MRSSM2_mass_eigenstates& model, const softsusy::QedQcd& qedqcd) {
      return forge<
         fields::Fe, fields::VP,
         calculate_Fe_Fe_VP_form_factors,
         npointfunctions::boxes1loop,
         npointfunctions::scalars1loop,
         npointfunctions::vectors1loop
      >(nI, nO, nA, model, qedqcd);
}

} // namespace MRSSM2_br_l_to_3l
} // namespace flexiblesusy
