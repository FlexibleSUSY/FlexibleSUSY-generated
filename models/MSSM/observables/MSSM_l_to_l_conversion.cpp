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
 * @file MSSM_l_to_l_conversion.cpp
 *
 * This file was generated at Sun 25 Feb 2024 21:19:21 with FlexibleSUSY
 * 2.8.0 and SARAH 4.15.1
 */

#include <valarray>
#include <complex>

#include "MSSM_mass_eigenstates.hpp"
#include "cxx_qft/MSSM_qft.hpp"

#include "MSSM_l_to_l_conversion.hpp"
#include "observables/l_to_l_conversion/settings.hpp"
#include "MSSM_FFV_form_factors.hpp"

#define F4F(quark, g) embed_photon<quark, Photon>(model, form_factors, g)
#define NPF(name, g) fix_tensors_sign<name>(model, in, out, g)
#define GS(N,Q) parameters.get(LToLConversion_settings::scalar_##N##Q)*m.N/m.Q
#define GV(N,Q) parameters.get(LToLConversion_settings::vector_##N##Q)
#define GT(N,Q) parameters.get(LToLConversion_settings::tensor_##N##Q)

#include "loop_libraries/loop_library.hpp"
#include "cxx_qft/MSSM_npointfunctions_wilsoncoeffs.hpp"
#include "concatenate.hpp"
#include <limits>
#include <type_traits>
#include <boost/fusion/include/at_key.hpp>
#include "wrappers.hpp"

namespace flexiblesusy {

namespace MSSM_cxx_diagrams {
namespace npointfunctions {

class nPointFeFuFeFu_3917884760616542 : public correlation_function_context<4,0>
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
   struct GenericV5Key {};
   struct GenericF5Key {};

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
         using _Fu = typename cxx_diagrams::bar<fields::Fu>::type;
         using _gF7 = typename cxx_diagrams::bar<gF7>::type;
         using _gF8 = typename cxx_diagrams::bar<gF8>::type;
         using _gS5 = typename cxx_diagrams::conj<gS5>::type;
         using _gS6 = typename cxx_diagrams::conj<gS6>::type;
         using Fe = fields::Fe;
         using Fu = fields::Fu;
         
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
            g3 = NPF_L(_Fu, Fu, _gS5) NPF_I(i4, i2, iS5);
            g4 = NPF_R(_Fu, Fu, _gS5) NPF_I(i4, i2, iS5);
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
         boost::mpl::vector<fields::Cha, typename cxx_diagrams::bar<fields::Cha>::type, fields::hh, fields::Sv>,
         boost::mpl::vector<fields::Chi, fields::Chi, fields::hh, fields::Se>,
         boost::mpl::vector<fields::Cha, typename cxx_diagrams::bar<fields::Cha>::type, fields::Ah, fields::Sv>,
         boost::mpl::vector<fields::Chi, fields::Chi, fields::Ah, fields::Se>
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
         using _Fu = typename cxx_diagrams::bar<fields::Fu>::type;
         using _gF6 = typename cxx_diagrams::bar<gF6>::type;
         using _gS5 = typename cxx_diagrams::conj<gS5>::type;
         using _gS7 = typename cxx_diagrams::conj<gS7>::type;
         using _gS8 = typename cxx_diagrams::conj<gS8>::type;
         using Fe = fields::Fe;
         using Fu = fields::Fu;
         
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
            g3 = NPF_L(_Fu, Fu, _gS5) NPF_I(i4, i2, iS5);
            g4 = NPF_R(_Fu, Fu, _gS5) NPF_I(i4, i2, iS5);
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
         boost::mpl::vector<fields::Cha, fields::hh, fields::Sv, typename cxx_diagrams::conj<fields::Sv>::type>,
         boost::mpl::vector<fields::Chi, fields::hh, fields::Se, typename cxx_diagrams::conj<fields::Se>::type>,
         boost::mpl::vector<fields::Fv, fields::Ah, fields::Hpm, typename cxx_diagrams::conj<fields::Hpm>::type>,
         boost::mpl::vector<fields::Chi, fields::Ah, fields::Se, typename cxx_diagrams::conj<fields::Se>::type>
         >;
      using combinatorial_factors = boost::mpl::vector<
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
         using _Fu = typename cxx_diagrams::bar<fields::Fu>::type;
         using _gF7 = typename cxx_diagrams::bar<gF7>::type;
         using _gF8 = typename cxx_diagrams::bar<gF8>::type;
         using _gS6 = typename cxx_diagrams::conj<gS6>::type;
         using _gV5 = typename cxx_diagrams::conj<gV5>::type;
         using Fe = fields::Fe;
         using Fu = fields::Fu;
         
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
            g3 = NPF_L(_Fu, Fu, _gV5) NPF_I(i4, i2, iV5);
            g4 = NPF_R(_Fu, Fu, _gV5) NPF_I(i4, i2, iV5);
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
   }; // End of struct genericSum3_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum3( void ) {
      using GenericKeys = boost::mpl::vector< GenericF7Key, GenericF8Key, GenericS6Key, GenericV5Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<fields::Fv, typename cxx_diagrams::bar<fields::Fv>::type, fields::Hpm, fields::VZ>,
         boost::mpl::vector<fields::Cha, typename cxx_diagrams::bar<fields::Cha>::type, fields::Sv, fields::VZ>,
         boost::mpl::vector<fields::Chi, fields::Chi, fields::Se, fields::VZ>
         >;
      using combinatorial_factors = boost::mpl::vector<
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>
         >;
      using color_factors = boost::mpl::vector<
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
         using _Fu = typename cxx_diagrams::bar<fields::Fu>::type;
         using _gF6 = typename cxx_diagrams::bar<gF6>::type;
         using _gS7 = typename cxx_diagrams::conj<gS7>::type;
         using _gS8 = typename cxx_diagrams::conj<gS8>::type;
         using _gV5 = typename cxx_diagrams::conj<gV5>::type;
         using Fe = fields::Fe;
         using Fu = fields::Fu;
         
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
            g3 = NPF_L(_Fu, Fu, _gV5) NPF_I(i4, i2, iV5);
            g4 = NPF_R(_Fu, Fu, _gV5) NPF_I(i4, i2, iV5);
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
   }; // End of struct genericSum4_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum4( void ) {
      using GenericKeys = boost::mpl::vector< GenericF6Key, GenericS7Key, GenericS8Key, GenericV5Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<fields::Fv, fields::Hpm, typename cxx_diagrams::conj<fields::Hpm>::type, fields::VZ>,
         boost::mpl::vector<fields::Cha, fields::Sv, typename cxx_diagrams::conj<fields::Sv>::type, fields::VZ>,
         boost::mpl::vector<fields::Chi, fields::Se, typename cxx_diagrams::conj<fields::Se>::type, fields::VZ>
         >;
      using combinatorial_factors = boost::mpl::vector<
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>
         >;
      using color_factors = boost::mpl::vector<
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
   
   template<class GenericFieldMap>
   struct genericSum5_impl : generic_sum_base {
      genericSum5_impl( const generic_sum_base &base ) :
      generic_sum_base( base ) {
      } // End of constructor genericSum5_impl
   
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
         using _Fu = typename cxx_diagrams::bar<fields::Fu>::type;
         using _gF5 = typename cxx_diagrams::bar<gF5>::type;
         using _gF8 = typename cxx_diagrams::bar<gF8>::type;
         using _gS6 = typename cxx_diagrams::conj<gS6>::type;
         using _gS7 = typename cxx_diagrams::conj<gS7>::type;
         using Fe = fields::Fe;
         using Fu = fields::Fu;
         
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
            
            g1 = NPF_L(_Fu, gF8, gS7) NPF_I(i4, iF8, iS7);
            g2 = NPF_R(_Fu, gF8, gS7) NPF_I(i4, iF8, iS7);
            g3 = NPF_L(_gF5, Fe, _gS6) NPF_I(iF5, i1, iS6);
            g4 = NPF_R(_gF5, Fe, _gS6) NPF_I(iF5, i1, iS6);
            g5 = NPF_L(_gF8, _Fe, gS6) NPF_I(iF8, i3, iS6);
            g6 = NPF_R(_gF8, _Fe, gS6) NPF_I(iF8, i3, iS6);
            g7 = NPF_L(Fu, gF5, _gS7) NPF_I(i2, iF5, iS7);
            g8 = NPF_R(Fu, gF5, _gS7) NPF_I(i2, iF5, iS7);
            
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
   }; // End of struct genericSum5_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum5( void ) {
      using GenericKeys = boost::mpl::vector< GenericF5Key, GenericF8Key, GenericS6Key, GenericS7Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<fields::Chi, fields::Chi, fields::Se, fields::Su>,
         boost::mpl::vector<fields::Cha, typename cxx_diagrams::bar<fields::Cha>::type, fields::Sv, fields::Sd>
         >;
      using combinatorial_factors = boost::mpl::vector<
         boost::mpl::int_<1>,
         boost::mpl::int_<1>
         >;
      using color_factors = boost::mpl::vector<
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>
         >;
      return accumulate_generic<
         GenericKeys,
         GenericInsertions,
         combinatorial_factors,
         color_factors,
         boost::mpl::int_<10>,
         genericSum5_impl
         >( *this );
   } // End of function genericSum5()
   
   template<class GenericFieldMap>
   struct genericSum6_impl : generic_sum_base {
      genericSum6_impl( const generic_sum_base &base ) :
      generic_sum_base( base ) {
      } // End of constructor genericSum6_impl
   
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
         using _Fu = typename cxx_diagrams::bar<fields::Fu>::type;
         using _gF6 = typename cxx_diagrams::bar<gF6>::type;
         using _gF7 = typename cxx_diagrams::bar<gF7>::type;
         using _gS5 = typename cxx_diagrams::conj<gS5>::type;
         using _gS8 = typename cxx_diagrams::conj<gS8>::type;
         using Fe = fields::Fe;
         using Fu = fields::Fu;
         
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
            
            g1 = NPF_L(_Fe, gF6, _gS8) NPF_I(i3, iF6, iS8);
            g2 = NPF_R(_Fe, gF6, _gS8) NPF_I(i3, iF6, iS8);
            g3 = NPF_L(_Fu, gF7, gS8) NPF_I(i4, iF7, iS8);
            g4 = NPF_R(_Fu, gF7, gS8) NPF_I(i4, iF7, iS8);
            g5 = NPF_L(_gF6, Fe, _gS5) NPF_I(iF6, i1, iS5);
            g6 = NPF_R(_gF6, Fe, _gS5) NPF_I(iF6, i1, iS5);
            g7 = NPF_L(_gF7, Fu, gS5) NPF_I(iF7, i2, iS5);
            g8 = NPF_R(_gF7, Fu, gS5) NPF_I(iF7, i2, iS5);
            
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
         
            S_LL += -0.006332573977646111*g1*g3*g5*g7*l1[dd0]*mF6*mF7;
            S_LR += -0.006332573977646111*g1*g4*g5*g8*l1[dd0]*mF6*mF7;
            S_RL += -0.006332573977646111*g2*g3*g6*g7*l1[dd0]*mF6*mF7;
            S_RR += -0.006332573977646111*g2*g4*g6*g8*l1[dd0]*mF6*mF7;
            V_LL += 0.006332573977646111*g2*g4*g5*g7*l1[dd00];
            V_LR += 0.006332573977646111*g2*g3*g5*g8*l1[dd00];
            V_RL += 0.006332573977646111*g1*g4*g6*g7*l1[dd00];
            V_RR += 0.006332573977646111*g1*g3*g6*g8*l1[dd00];
            T_LL += 0;
            T_RR += 0;
         
         }}}} // End of summation over generic fields
   
         return {S_LL, S_LR, S_RL, S_RR, V_LL, V_LR, V_RL, V_RR, T_LL, T_RR};
      } // End of operator()( void )
   }; // End of struct genericSum6_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum6( void ) {
      using GenericKeys = boost::mpl::vector< GenericF6Key, GenericF7Key, GenericS5Key, GenericS8Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<fields::Fv, fields::Fd, fields::Hpm, typename cxx_diagrams::conj<fields::Hpm>::type>
         >;
      using combinatorial_factors = boost::mpl::vector<
         boost::mpl::int_<1>
         >;
      using color_factors = boost::mpl::vector<
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>
         >;
      return accumulate_generic<
         GenericKeys,
         GenericInsertions,
         combinatorial_factors,
         color_factors,
         boost::mpl::int_<10>,
         genericSum6_impl
         >( *this );
   } // End of function genericSum6()
   
   template<class GenericFieldMap>
   struct genericSum7_impl : generic_sum_base {
      genericSum7_impl( const generic_sum_base &base ) :
      generic_sum_base( base ) {
      } // End of constructor genericSum7_impl
   
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
         using _Fu = typename cxx_diagrams::bar<fields::Fu>::type;
         using _gF6 = typename cxx_diagrams::bar<gF6>::type;
         using _gF7 = typename cxx_diagrams::bar<gF7>::type;
         using _gS5 = typename cxx_diagrams::conj<gS5>::type;
         using _gS8 = typename cxx_diagrams::conj<gS8>::type;
         using Fe = fields::Fe;
         using Fu = fields::Fu;
         
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
            g3 = NPF_L(_Fu, gF6, gS8) NPF_I(i4, iF6, iS8);
            g4 = NPF_R(_Fu, gF6, gS8) NPF_I(i4, iF6, iS8);
            g5 = NPF_L(_gF6, Fe, _gS5) NPF_I(iF6, i1, iS5);
            g6 = NPF_R(_gF6, Fe, _gS5) NPF_I(iF6, i1, iS5);
            g7 = NPF_L(_gF7, Fu, _gS8) NPF_I(iF7, i2, iS8);
            g8 = NPF_R(_gF7, Fu, _gS8) NPF_I(iF7, i2, iS8);
            
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
   }; // End of struct genericSum7_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum7( void ) {
      using GenericKeys = boost::mpl::vector< GenericF6Key, GenericF7Key, GenericS5Key, GenericS8Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<fields::Chi, fields::Chi, fields::Se, fields::Su>
         >;
      using combinatorial_factors = boost::mpl::vector<
         boost::mpl::int_<1>
         >;
      using color_factors = boost::mpl::vector<
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>
         >;
      return accumulate_generic<
         GenericKeys,
         GenericInsertions,
         combinatorial_factors,
         color_factors,
         boost::mpl::int_<10>,
         genericSum7_impl
         >( *this );
   } // End of function genericSum7()
   
   template<class GenericFieldMap>
   struct genericSum8_impl : generic_sum_base {
      genericSum8_impl( const generic_sum_base &base ) :
      generic_sum_base( base ) {
      } // End of constructor genericSum8_impl
   
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
         using _Fu = typename cxx_diagrams::bar<fields::Fu>::type;
         using _gF6 = typename cxx_diagrams::bar<gF6>::type;
         using _gF7 = typename cxx_diagrams::bar<gF7>::type;
         using _gS5 = typename cxx_diagrams::conj<gS5>::type;
         using _gS8 = typename cxx_diagrams::conj<gS8>::type;
         using Fe = fields::Fe;
         using Fu = fields::Fu;
         
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
            g3 = NPF_L(_Fu, Fu, gS5) NPF_I(i4, i2, iS5);
            g4 = NPF_R(_Fu, Fu, gS5) NPF_I(i4, i2, iS5);
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
   }; // End of struct genericSum8_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum8( void ) {
      using GenericKeys = boost::mpl::vector< GenericF6Key, GenericF7Key, GenericS5Key, GenericS8Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<fields::Fe, typename cxx_diagrams::bar<fields::Fv>::type, fields::hh, typename cxx_diagrams::conj<fields::Hpm>::type>,
         boost::mpl::vector<fields::Fe, typename cxx_diagrams::bar<fields::Cha>::type, fields::hh, typename cxx_diagrams::conj<fields::Sv>::type>,
         boost::mpl::vector<fields::Fe, fields::Chi, fields::hh, typename cxx_diagrams::conj<fields::Se>::type>,
         boost::mpl::vector<fields::Fe, typename cxx_diagrams::bar<fields::Fv>::type, fields::Ah, typename cxx_diagrams::conj<fields::Hpm>::type>,
         boost::mpl::vector<fields::Fe, typename cxx_diagrams::bar<fields::Cha>::type, fields::Ah, typename cxx_diagrams::conj<fields::Sv>::type>,
         boost::mpl::vector<fields::Fe, fields::Chi, fields::Ah, typename cxx_diagrams::conj<fields::Se>::type>
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
         genericSum8_impl
         >( *this );
   } // End of function genericSum8()
   
   template<class GenericFieldMap>
   struct genericSum9_impl : generic_sum_base {
      genericSum9_impl( const generic_sum_base &base ) :
      generic_sum_base( base ) {
      } // End of constructor genericSum9_impl
   
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
         using _Fu = typename cxx_diagrams::bar<fields::Fu>::type;
         using _gF6 = typename cxx_diagrams::bar<gF6>::type;
         using _gF7 = typename cxx_diagrams::bar<gF7>::type;
         using _gS8 = typename cxx_diagrams::conj<gS8>::type;
         using _gV5 = typename cxx_diagrams::conj<gV5>::type;
         using Fe = fields::Fe;
         using Fu = fields::Fu;
         
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
            g3 = NPF_L(_Fu, Fu, gV5) NPF_I(i4, i2, iV5);
            g4 = NPF_R(_Fu, Fu, gV5) NPF_I(i4, i2, iV5);
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
   }; // End of struct genericSum9_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum9( void ) {
      using GenericKeys = boost::mpl::vector< GenericF6Key, GenericF7Key, GenericS8Key, GenericV5Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<fields::Fe, typename cxx_diagrams::bar<fields::Fv>::type, typename cxx_diagrams::conj<fields::Hpm>::type, fields::VZ>,
         boost::mpl::vector<fields::Fe, typename cxx_diagrams::bar<fields::Cha>::type, typename cxx_diagrams::conj<fields::Sv>::type, fields::VZ>,
         boost::mpl::vector<fields::Fe, fields::Chi, typename cxx_diagrams::conj<fields::Se>::type, fields::VZ>
         >;
      using combinatorial_factors = boost::mpl::vector<
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>
         >;
      using color_factors = boost::mpl::vector<
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
         genericSum9_impl
         >( *this );
   } // End of function genericSum9()
   
   template<class GenericFieldMap>
   struct genericSum10_impl : generic_sum_base {
      genericSum10_impl( const generic_sum_base &base ) :
      generic_sum_base( base ) {
      } // End of constructor genericSum10_impl
   
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
         using _Fu = typename cxx_diagrams::bar<fields::Fu>::type;
         using _gF6 = typename cxx_diagrams::bar<gF6>::type;
         using _gF7 = typename cxx_diagrams::bar<gF7>::type;
         using _gS5 = typename cxx_diagrams::conj<gS5>::type;
         using _gS8 = typename cxx_diagrams::conj<gS8>::type;
         using Fe = fields::Fe;
         using Fu = fields::Fu;
         
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
            g3 = NPF_L(_Fu, Fu, _gS5) NPF_I(i4, i2, iS5);
            g4 = NPF_R(_Fu, Fu, _gS5) NPF_I(i4, i2, iS5);
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
   }; // End of struct genericSum10_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum10( void ) {
      using GenericKeys = boost::mpl::vector< GenericF6Key, GenericF7Key, GenericS5Key, GenericS8Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, fields::Fv, fields::hh, fields::Hpm>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, fields::Cha, fields::hh, fields::Sv>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, fields::Chi, fields::hh, fields::Se>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, fields::Fv, fields::Ah, fields::Hpm>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, fields::Cha, fields::Ah, fields::Sv>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, fields::Chi, fields::Ah, fields::Se>
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
         genericSum10_impl
         >( *this );
   } // End of function genericSum10()
   
   template<class GenericFieldMap>
   struct genericSum11_impl : generic_sum_base {
      genericSum11_impl( const generic_sum_base &base ) :
      generic_sum_base( base ) {
      } // End of constructor genericSum11_impl
   
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
         using _Fu = typename cxx_diagrams::bar<fields::Fu>::type;
         using _gF6 = typename cxx_diagrams::bar<gF6>::type;
         using _gF7 = typename cxx_diagrams::bar<gF7>::type;
         using _gS8 = typename cxx_diagrams::conj<gS8>::type;
         using _gV5 = typename cxx_diagrams::conj<gV5>::type;
         using Fe = fields::Fe;
         using Fu = fields::Fu;
         
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
            g3 = NPF_L(_Fu, Fu, _gV5) NPF_I(i4, i2, iV5);
            g4 = NPF_R(_Fu, Fu, _gV5) NPF_I(i4, i2, iV5);
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
   }; // End of struct genericSum11_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum11( void ) {
      using GenericKeys = boost::mpl::vector< GenericF6Key, GenericF7Key, GenericS8Key, GenericV5Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, fields::Fv, fields::Hpm, fields::VZ>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, fields::Cha, fields::Sv, fields::VZ>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, fields::Chi, fields::Se, fields::VZ>
         >;
      using combinatorial_factors = boost::mpl::vector<
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>
         >;
      using color_factors = boost::mpl::vector<
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
         genericSum11_impl
         >( *this );
   } // End of function genericSum11()

   public:
   nPointFeFuFeFu_3917884760616542( const MSSM_mass_eigenstates &model,
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
      const auto genericsum5 = genericSum5();
      const auto genericsum6 = genericSum6();
      const auto genericsum7 = genericSum7();
      const auto genericsum8 = genericSum8();
      const auto genericsum9 = genericSum9();
      const auto genericsum10 = genericSum10();
      const auto genericsum11 = genericSum11();
   
      for ( std::size_t i=0; i<coeffsLength; i++ ) {
         genericSummation[i] += genericsum1[i]+genericsum2[i]+genericsum3[i]+genericsum4[i]+genericsum5[i]+genericsum6[i]+genericsum7[i]+genericsum8[i]+genericsum9[i]+genericsum10[i]+genericsum11[i];
      }
      return genericSummation;
   } // End of calculate()
}; // End of nPointFeFuFeFu_3917884760616542

std::array<std::complex<double>, 10> conversion_FeFu_to_FeFu_All1loop(const MSSM_mass_eigenstates &model,
const std::array<int, 4> &indices,
const std::array<Eigen::Vector4d, 0> &momenta
) {
   nPointFeFuFeFu_3917884760616542 helper{model, indices, momenta};
   return helper.calculate();
}

class nPointFeFdFeFd_3917884761161352 : public correlation_function_context<4,0>
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
   struct GenericV5Key {};
   struct GenericF5Key {};

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
         using _Fd = typename cxx_diagrams::bar<fields::Fd>::type;
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _gF7 = typename cxx_diagrams::bar<gF7>::type;
         using _gF8 = typename cxx_diagrams::bar<gF8>::type;
         using _gS5 = typename cxx_diagrams::conj<gS5>::type;
         using _gS6 = typename cxx_diagrams::conj<gS6>::type;
         using Fd = fields::Fd;
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
            
            g1 = NPF_L(_Fd, Fd, _gS5) NPF_I(i4, i2, iS5);
            g2 = NPF_R(_Fd, Fd, _gS5) NPF_I(i4, i2, iS5);
            g3 = NPF_L(_Fe, _gF8, gS6) NPF_I(i3, iF8, iS6);
            g4 = NPF_R(_Fe, _gF8, gS6) NPF_I(i3, iF8, iS6);
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
            
            if( (z[0] + z[1])*(z[2]*z[4] + z[3]*z[5])*(z[6] + z[7]) == 0 ) continue;
            Loop_library::get().B(l1,0,Sqr(mF7),Sqr(mF8),Sqr(context.scale()));
            Loop_library::get().C(l2,0,0,0,Sqr(mF8),Sqr(mF7),Sqr(mS6),Sqr(context.scale()));
         
            S_LL += (0.006332573977646111*g1*g3*g5*(g7*l2[cc0]*mF7*mF8 + g8*(l1[bb0] + l2[cc0]*Sqr(mS6))))/Sqr(mS5);
            S_LR += (0.006332573977646111*g2*g3*g5*(g7*l2[cc0]*mF7*mF8 + g8*(l1[bb0] + l2[cc0]*Sqr(mS6))))/Sqr(mS5);
            S_RL += (0.006332573977646111*g1*g4*g6*(g8*l2[cc0]*mF7*mF8 + g7*(l1[bb0] + l2[cc0]*Sqr(mS6))))/Sqr(mS5);
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
         boost::mpl::vector<fields::Cha, typename cxx_diagrams::bar<fields::Cha>::type, fields::hh, fields::Sv>,
         boost::mpl::vector<fields::Chi, fields::Chi, fields::hh, fields::Se>,
         boost::mpl::vector<fields::Cha, typename cxx_diagrams::bar<fields::Cha>::type, fields::Ah, fields::Sv>,
         boost::mpl::vector<fields::Chi, fields::Chi, fields::Ah, fields::Se>
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
         using _Fd = typename cxx_diagrams::bar<fields::Fd>::type;
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _gF6 = typename cxx_diagrams::bar<gF6>::type;
         using _gS5 = typename cxx_diagrams::conj<gS5>::type;
         using _gS7 = typename cxx_diagrams::conj<gS7>::type;
         using _gS8 = typename cxx_diagrams::conj<gS8>::type;
         using Fd = fields::Fd;
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
            
            g1 = NPF_L(_Fd, Fd, _gS5) NPF_I(i4, i2, iS5);
            g2 = NPF_R(_Fd, Fd, _gS5) NPF_I(i4, i2, iS5);
            g3 = NPF_L(_Fe, gF6, _gS8) NPF_I(i3, iF6, iS8);
            g4 = NPF_R(_Fe, gF6, _gS8) NPF_I(i3, iF6, iS8);
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
            
            if( (z[0] + z[1])*(z[2]*z[4] + z[3]*z[5])*z[6] == 0 ) continue;
            Loop_library::get().C(l1,0,0,0,Sqr(mF6),Sqr(mS8),Sqr(mS7),Sqr(context.scale()));
         
            S_LL += (0.006332573977646111*g1*g3*g5*g7*l1[cc0]*mF6)/Sqr(mS5);
            S_LR += (0.006332573977646111*g2*g3*g5*g7*l1[cc0]*mF6)/Sqr(mS5);
            S_RL += (0.006332573977646111*g1*g4*g6*g7*l1[cc0]*mF6)/Sqr(mS5);
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
         boost::mpl::vector<fields::Cha, fields::hh, fields::Sv, typename cxx_diagrams::conj<fields::Sv>::type>,
         boost::mpl::vector<fields::Chi, fields::hh, fields::Se, typename cxx_diagrams::conj<fields::Se>::type>,
         boost::mpl::vector<fields::Fv, fields::Ah, fields::Hpm, typename cxx_diagrams::conj<fields::Hpm>::type>,
         boost::mpl::vector<fields::Chi, fields::Ah, fields::Se, typename cxx_diagrams::conj<fields::Se>::type>
         >;
      using combinatorial_factors = boost::mpl::vector<
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
         using _Fd = typename cxx_diagrams::bar<fields::Fd>::type;
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _gF7 = typename cxx_diagrams::bar<gF7>::type;
         using _gF8 = typename cxx_diagrams::bar<gF8>::type;
         using _gS6 = typename cxx_diagrams::conj<gS6>::type;
         using _gV5 = typename cxx_diagrams::conj<gV5>::type;
         using Fd = fields::Fd;
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
            
            g1 = NPF_L(_Fd, Fd, _gV5) NPF_I(i4, i2, iV5);
            g2 = NPF_R(_Fd, Fd, _gV5) NPF_I(i4, i2, iV5);
            g3 = NPF_L(_Fe, _gF8, gS6) NPF_I(i3, iF8, iS6);
            g4 = NPF_R(_Fe, _gF8, gS6) NPF_I(i3, iF8, iS6);
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
            
            if( (z[0] + z[1])*(z[3]*z[4] + z[2]*z[5])*(z[6] + z[7]) == 0 ) continue;
            Loop_library::get().B(l1,0,Sqr(mF7),Sqr(mF8),Sqr(context.scale()));
            Loop_library::get().C(l2,0,0,0,Sqr(mF8),Sqr(mF7),Sqr(mS6),Sqr(context.scale()));
         
            S_LL += 0;
            S_LR += 0;
            S_RL += 0;
            S_RR += 0;
            V_LL += (0.006332573977646111*g1*g4*g5*(-(g7*l2[cc0]*mF7*mF8) + g8*(l1[bb0] - 2*l2[cc00] + l2[cc0]*Sqr(mS6))))/Sqr(mV5);
            V_LR += (0.006332573977646111*g2*g4*g5*(-(g7*l2[cc0]*mF7*mF8) + g8*(l1[bb0] - 2*l2[cc00] + l2[cc0]*Sqr(mS6))))/Sqr(mV5);
            V_RL += (0.006332573977646111*g1*g3*g6*(-(g8*l2[cc0]*mF7*mF8) + g7*(l1[bb0] - 2*l2[cc00] + l2[cc0]*Sqr(mS6))))/Sqr(mV5);
            V_RR += (0.006332573977646111*g2*g3*g6*(-(g8*l2[cc0]*mF7*mF8) + g7*(l1[bb0] - 2*l2[cc00] + l2[cc0]*Sqr(mS6))))/Sqr(mV5);
            T_LL += 0;
            T_RR += 0;
         
         }}}} // End of summation over generic fields
   
         return {S_LL, S_LR, S_RL, S_RR, V_LL, V_LR, V_RL, V_RR, T_LL, T_RR};
      } // End of operator()( void )
   }; // End of struct genericSum3_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum3( void ) {
      using GenericKeys = boost::mpl::vector< GenericF7Key, GenericF8Key, GenericS6Key, GenericV5Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<fields::Fv, typename cxx_diagrams::bar<fields::Fv>::type, fields::Hpm, fields::VZ>,
         boost::mpl::vector<fields::Cha, typename cxx_diagrams::bar<fields::Cha>::type, fields::Sv, fields::VZ>,
         boost::mpl::vector<fields::Chi, fields::Chi, fields::Se, fields::VZ>
         >;
      using combinatorial_factors = boost::mpl::vector<
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>
         >;
      using color_factors = boost::mpl::vector<
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
         using _Fd = typename cxx_diagrams::bar<fields::Fd>::type;
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _gF6 = typename cxx_diagrams::bar<gF6>::type;
         using _gS7 = typename cxx_diagrams::conj<gS7>::type;
         using _gS8 = typename cxx_diagrams::conj<gS8>::type;
         using _gV5 = typename cxx_diagrams::conj<gV5>::type;
         using Fd = fields::Fd;
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
            
            g1 = NPF_L(_Fd, Fd, _gV5) NPF_I(i4, i2, iV5);
            g2 = NPF_R(_Fd, Fd, _gV5) NPF_I(i4, i2, iV5);
            g3 = NPF_L(_Fe, gF6, _gS8) NPF_I(i3, iF6, iS8);
            g4 = NPF_R(_Fe, gF6, _gS8) NPF_I(i3, iF6, iS8);
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
            
            if( (z[0] + z[1])*(z[3]*z[4] + z[2]*z[5])*z[6] == 0 ) continue;
            Loop_library::get().C(l1,0,0,0,Sqr(mF6),Sqr(mS8),Sqr(mS7),Sqr(context.scale()));
         
            S_LL += 0;
            S_LR += 0;
            S_RL += 0;
            S_RR += 0;
            V_LL += (0.012665147955292222*g1*g4*g5*g7*l1[cc00])/Sqr(mV5);
            V_LR += (0.012665147955292222*g2*g4*g5*g7*l1[cc00])/Sqr(mV5);
            V_RL += (0.012665147955292222*g1*g3*g6*g7*l1[cc00])/Sqr(mV5);
            V_RR += (0.012665147955292222*g2*g3*g6*g7*l1[cc00])/Sqr(mV5);
            T_LL += 0;
            T_RR += 0;
         
         }}}} // End of summation over generic fields
   
         return {S_LL, S_LR, S_RL, S_RR, V_LL, V_LR, V_RL, V_RR, T_LL, T_RR};
      } // End of operator()( void )
   }; // End of struct genericSum4_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum4( void ) {
      using GenericKeys = boost::mpl::vector< GenericF6Key, GenericS7Key, GenericS8Key, GenericV5Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<fields::Fv, fields::Hpm, typename cxx_diagrams::conj<fields::Hpm>::type, fields::VZ>,
         boost::mpl::vector<fields::Cha, fields::Sv, typename cxx_diagrams::conj<fields::Sv>::type, fields::VZ>,
         boost::mpl::vector<fields::Chi, fields::Se, typename cxx_diagrams::conj<fields::Se>::type, fields::VZ>
         >;
      using combinatorial_factors = boost::mpl::vector<
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>
         >;
      using color_factors = boost::mpl::vector<
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
   
   template<class GenericFieldMap>
   struct genericSum5_impl : generic_sum_base {
      genericSum5_impl( const generic_sum_base &base ) :
      generic_sum_base( base ) {
      } // End of constructor genericSum5_impl
   
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
         using _Fd = typename cxx_diagrams::bar<fields::Fd>::type;
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _gF5 = typename cxx_diagrams::bar<gF5>::type;
         using _gF8 = typename cxx_diagrams::bar<gF8>::type;
         using _gS6 = typename cxx_diagrams::conj<gS6>::type;
         using _gS7 = typename cxx_diagrams::conj<gS7>::type;
         using Fd = fields::Fd;
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
            
            g1 = NPF_L(_Fd, gF8, gS7) NPF_I(i4, iF8, iS7);
            g2 = NPF_R(_Fd, gF8, gS7) NPF_I(i4, iF8, iS7);
            g3 = NPF_L(_gF5, Fe, _gS6) NPF_I(iF5, i1, iS6);
            g4 = NPF_R(_gF5, Fe, _gS6) NPF_I(iF5, i1, iS6);
            g5 = NPF_L(_gF8, _Fe, gS6) NPF_I(iF8, i3, iS6);
            g6 = NPF_R(_gF8, _Fe, gS6) NPF_I(iF8, i3, iS6);
            g7 = NPF_L(Fd, gF5, _gS7) NPF_I(i2, iF5, iS7);
            g8 = NPF_R(Fd, gF5, _gS7) NPF_I(i2, iF5, iS7);
            
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
   }; // End of struct genericSum5_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum5( void ) {
      using GenericKeys = boost::mpl::vector< GenericF5Key, GenericF8Key, GenericS6Key, GenericS7Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<fields::Chi, fields::Chi, fields::Se, fields::Sd>
         >;
      using combinatorial_factors = boost::mpl::vector<
         boost::mpl::int_<1>
         >;
      using color_factors = boost::mpl::vector<
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>
         >;
      return accumulate_generic<
         GenericKeys,
         GenericInsertions,
         combinatorial_factors,
         color_factors,
         boost::mpl::int_<10>,
         genericSum5_impl
         >( *this );
   } // End of function genericSum5()
   
   template<class GenericFieldMap>
   struct genericSum6_impl : generic_sum_base {
      genericSum6_impl( const generic_sum_base &base ) :
      generic_sum_base( base ) {
      } // End of constructor genericSum6_impl
   
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
         using _Fd = typename cxx_diagrams::bar<fields::Fd>::type;
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _gF5 = typename cxx_diagrams::bar<gF5>::type;
         using _gF8 = typename cxx_diagrams::bar<gF8>::type;
         using _gS6 = typename cxx_diagrams::conj<gS6>::type;
         using _gS7 = typename cxx_diagrams::conj<gS7>::type;
         using Fd = fields::Fd;
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
            
            g1 = NPF_L(_Fd, gF8, gS6) NPF_I(i4, iF8, iS6);
            g2 = NPF_R(_Fd, gF8, gS6) NPF_I(i4, iF8, iS6);
            g3 = NPF_L(_Fe, gF5, gS7) NPF_I(i3, iF5, iS7);
            g4 = NPF_R(_Fe, gF5, gS7) NPF_I(i3, iF5, iS7);
            g5 = NPF_L(_gF5, Fe, _gS6) NPF_I(iF5, i1, iS6);
            g6 = NPF_R(_gF5, Fe, _gS6) NPF_I(iF5, i1, iS6);
            g7 = NPF_L(_gF8, Fd, _gS7) NPF_I(iF8, i2, iS7);
            g8 = NPF_R(_gF8, Fd, _gS7) NPF_I(iF8, i2, iS7);
            
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
            S_LR += -0.006332573977646111*g2*g3*g5*g8*l1[dd0]*mF5*mF8;
            S_RL += -0.006332573977646111*g1*g4*g6*g7*l1[dd0]*mF5*mF8;
            S_RR += -0.006332573977646111*g2*g4*g6*g8*l1[dd0]*mF5*mF8;
            V_LL += -0.006332573977646111*g2*g4*g5*g7*l1[dd00];
            V_LR += -0.006332573977646111*g1*g4*g5*g8*l1[dd00];
            V_RL += -0.006332573977646111*g2*g3*g6*g7*l1[dd00];
            V_RR += -0.006332573977646111*g1*g3*g6*g8*l1[dd00];
            T_LL += 0;
            T_RR += 0;
         
         }}}} // End of summation over generic fields
   
         return {S_LL, S_LR, S_RL, S_RR, V_LL, V_LR, V_RL, V_RR, T_LL, T_RR};
      } // End of operator()( void )
   }; // End of struct genericSum6_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum6( void ) {
      using GenericKeys = boost::mpl::vector< GenericF5Key, GenericF8Key, GenericS6Key, GenericS7Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<fields::Fv, fields::Fu, fields::Hpm, fields::Hpm>
         >;
      using combinatorial_factors = boost::mpl::vector<
         boost::mpl::int_<1>
         >;
      using color_factors = boost::mpl::vector<
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>
         >;
      return accumulate_generic<
         GenericKeys,
         GenericInsertions,
         combinatorial_factors,
         color_factors,
         boost::mpl::int_<10>,
         genericSum6_impl
         >( *this );
   } // End of function genericSum6()
   
   template<class GenericFieldMap>
   struct genericSum7_impl : generic_sum_base {
      genericSum7_impl( const generic_sum_base &base ) :
      generic_sum_base( base ) {
      } // End of constructor genericSum7_impl
   
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
         using _Fd = typename cxx_diagrams::bar<fields::Fd>::type;
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _gF6 = typename cxx_diagrams::bar<gF6>::type;
         using _gF7 = typename cxx_diagrams::bar<gF7>::type;
         using _gS5 = typename cxx_diagrams::conj<gS5>::type;
         using _gS8 = typename cxx_diagrams::conj<gS8>::type;
         using Fd = fields::Fd;
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
            
            g1 = NPF_L(_Fd, gF6, gS8) NPF_I(i4, iF6, iS8);
            g2 = NPF_R(_Fd, gF6, gS8) NPF_I(i4, iF6, iS8);
            g3 = NPF_L(_Fe, gF7, gS5) NPF_I(i3, iF7, iS5);
            g4 = NPF_R(_Fe, gF7, gS5) NPF_I(i3, iF7, iS5);
            g5 = NPF_L(_gF6, Fe, _gS5) NPF_I(iF6, i1, iS5);
            g6 = NPF_R(_gF6, Fe, _gS5) NPF_I(iF6, i1, iS5);
            g7 = NPF_L(_gF7, Fd, _gS8) NPF_I(iF7, i2, iS8);
            g8 = NPF_R(_gF7, Fd, _gS8) NPF_I(iF7, i2, iS8);
            
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
            S_LR += 0.012665147955292222*g2*g3*g5*g8*l1[dd00];
            S_RL += 0.012665147955292222*g1*g4*g6*g7*l1[dd00];
            S_RR += 0.0031662869888230555*g2*g4*g6*g8*l1[dd0]*mF6*mF7;
            V_LL += -0.006332573977646111*g2*g4*g5*g7*l1[dd00];
            V_LR += 0.0031662869888230555*g1*g4*g5*g8*l1[dd0]*mF6*mF7;
            V_RL += 0.0031662869888230555*g2*g3*g6*g7*l1[dd0]*mF6*mF7;
            V_RR += -0.006332573977646111*g1*g3*g6*g8*l1[dd00];
            T_LL += -0.0007915717472057639*g1*g3*g5*g7*l1[dd0]*mF6*mF7;
            T_RR += -0.0007915717472057639*g2*g4*g6*g8*l1[dd0]*mF6*mF7;
         
         }}}} // End of summation over generic fields
   
         return {S_LL, S_LR, S_RL, S_RR, V_LL, V_LR, V_RL, V_RR, T_LL, T_RR};
      } // End of operator()( void )
   }; // End of struct genericSum7_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum7( void ) {
      using GenericKeys = boost::mpl::vector< GenericF6Key, GenericF7Key, GenericS5Key, GenericS8Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<fields::Cha, fields::Cha, fields::Sv, fields::Su>,
         boost::mpl::vector<fields::Chi, fields::Chi, fields::Se, fields::Sd>
         >;
      using combinatorial_factors = boost::mpl::vector<
         boost::mpl::int_<1>,
         boost::mpl::int_<1>
         >;
      using color_factors = boost::mpl::vector<
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>,
         detail::complex_helper<detail::ratio_helper<1, 1>,detail::ratio_helper<0, 1>>
         >;
      return accumulate_generic<
         GenericKeys,
         GenericInsertions,
         combinatorial_factors,
         color_factors,
         boost::mpl::int_<10>,
         genericSum7_impl
         >( *this );
   } // End of function genericSum7()
   
   template<class GenericFieldMap>
   struct genericSum8_impl : generic_sum_base {
      genericSum8_impl( const generic_sum_base &base ) :
      generic_sum_base( base ) {
      } // End of constructor genericSum8_impl
   
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
         using _Fd = typename cxx_diagrams::bar<fields::Fd>::type;
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _gF6 = typename cxx_diagrams::bar<gF6>::type;
         using _gF7 = typename cxx_diagrams::bar<gF7>::type;
         using _gS5 = typename cxx_diagrams::conj<gS5>::type;
         using _gS8 = typename cxx_diagrams::conj<gS8>::type;
         using Fd = fields::Fd;
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
            
            g1 = NPF_L(_Fd, Fd, gS5) NPF_I(i4, i2, iS5);
            g2 = NPF_R(_Fd, Fd, gS5) NPF_I(i4, i2, iS5);
            g3 = NPF_L(_Fe, _gF7, _gS8) NPF_I(i3, iF7, iS8);
            g4 = NPF_R(_Fe, _gF7, _gS8) NPF_I(i3, iF7, iS8);
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
            
            if( (z[0] + z[1])*(z[2]*z[4]*z[6] + z[3]*z[5]*z[7]) == 0 ) continue;
            Loop_library::get().B(l1,0,Sqr(mF7),Sqr(mS8),Sqr(context.scale()));
         
            S_LL += (-0.006332573977646111*g1*g3*g5*g7*l1[bb0]*mF7)/(mF6*Sqr(mS5));
            S_LR += (-0.006332573977646111*g2*g3*g5*g7*l1[bb0]*mF7)/(mF6*Sqr(mS5));
            S_RL += (-0.006332573977646111*g1*g4*g6*g8*l1[bb0]*mF7)/(mF6*Sqr(mS5));
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
   }; // End of struct genericSum8_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum8( void ) {
      using GenericKeys = boost::mpl::vector< GenericF6Key, GenericF7Key, GenericS5Key, GenericS8Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<fields::Fe, typename cxx_diagrams::bar<fields::Fv>::type, fields::hh, typename cxx_diagrams::conj<fields::Hpm>::type>,
         boost::mpl::vector<fields::Fe, typename cxx_diagrams::bar<fields::Cha>::type, fields::hh, typename cxx_diagrams::conj<fields::Sv>::type>,
         boost::mpl::vector<fields::Fe, fields::Chi, fields::hh, typename cxx_diagrams::conj<fields::Se>::type>,
         boost::mpl::vector<fields::Fe, typename cxx_diagrams::bar<fields::Fv>::type, fields::Ah, typename cxx_diagrams::conj<fields::Hpm>::type>,
         boost::mpl::vector<fields::Fe, typename cxx_diagrams::bar<fields::Cha>::type, fields::Ah, typename cxx_diagrams::conj<fields::Sv>::type>,
         boost::mpl::vector<fields::Fe, fields::Chi, fields::Ah, typename cxx_diagrams::conj<fields::Se>::type>
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
         genericSum8_impl
         >( *this );
   } // End of function genericSum8()
   
   template<class GenericFieldMap>
   struct genericSum9_impl : generic_sum_base {
      genericSum9_impl( const generic_sum_base &base ) :
      generic_sum_base( base ) {
      } // End of constructor genericSum9_impl
   
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
         using _Fd = typename cxx_diagrams::bar<fields::Fd>::type;
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _gF6 = typename cxx_diagrams::bar<gF6>::type;
         using _gF7 = typename cxx_diagrams::bar<gF7>::type;
         using _gS8 = typename cxx_diagrams::conj<gS8>::type;
         using _gV5 = typename cxx_diagrams::conj<gV5>::type;
         using Fd = fields::Fd;
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
            
            g1 = NPF_L(_Fd, Fd, gV5) NPF_I(i4, i2, iV5);
            g2 = NPF_R(_Fd, Fd, gV5) NPF_I(i4, i2, iV5);
            g3 = NPF_L(_Fe, _gF7, _gS8) NPF_I(i3, iF7, iS8);
            g4 = NPF_R(_Fe, _gF7, _gS8) NPF_I(i3, iF7, iS8);
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
            
            if( (z[0] + z[1])*(z[2]*z[5]*z[6] + z[3]*z[4]*z[7]) == 0 ) continue;
            Loop_library::get().B(l1,0,Sqr(mF7),Sqr(mS8),Sqr(context.scale()));
         
            S_LL += 0;
            S_LR += 0;
            S_RL += 0;
            S_RR += 0;
            V_LL += (0.006332573977646111*g1*g4*g5*g8*l1[bb0]*mF7)/(mF6*Sqr(mV5));
            V_LR += (0.006332573977646111*g2*g4*g5*g8*l1[bb0]*mF7)/(mF6*Sqr(mV5));
            V_RL += (0.006332573977646111*g1*g3*g6*g7*l1[bb0]*mF7)/(mF6*Sqr(mV5));
            V_RR += (0.006332573977646111*g2*g3*g6*g7*l1[bb0]*mF7)/(mF6*Sqr(mV5));
            T_LL += 0;
            T_RR += 0;
         
         }}}} // End of summation over generic fields
   
         return {S_LL, S_LR, S_RL, S_RR, V_LL, V_LR, V_RL, V_RR, T_LL, T_RR};
      } // End of operator()( void )
   }; // End of struct genericSum9_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum9( void ) {
      using GenericKeys = boost::mpl::vector< GenericF6Key, GenericF7Key, GenericS8Key, GenericV5Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<fields::Fe, typename cxx_diagrams::bar<fields::Fv>::type, typename cxx_diagrams::conj<fields::Hpm>::type, fields::VZ>,
         boost::mpl::vector<fields::Fe, typename cxx_diagrams::bar<fields::Cha>::type, typename cxx_diagrams::conj<fields::Sv>::type, fields::VZ>,
         boost::mpl::vector<fields::Fe, fields::Chi, typename cxx_diagrams::conj<fields::Se>::type, fields::VZ>
         >;
      using combinatorial_factors = boost::mpl::vector<
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>
         >;
      using color_factors = boost::mpl::vector<
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
         genericSum9_impl
         >( *this );
   } // End of function genericSum9()
   
   template<class GenericFieldMap>
   struct genericSum10_impl : generic_sum_base {
      genericSum10_impl( const generic_sum_base &base ) :
      generic_sum_base( base ) {
      } // End of constructor genericSum10_impl
   
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
         using _Fd = typename cxx_diagrams::bar<fields::Fd>::type;
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _gF6 = typename cxx_diagrams::bar<gF6>::type;
         using _gF7 = typename cxx_diagrams::bar<gF7>::type;
         using _gS5 = typename cxx_diagrams::conj<gS5>::type;
         using _gS8 = typename cxx_diagrams::conj<gS8>::type;
         using Fd = fields::Fd;
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
            
            g1 = NPF_L(_Fd, Fd, _gS5) NPF_I(i4, i2, iS5);
            g2 = NPF_R(_Fd, Fd, _gS5) NPF_I(i4, i2, iS5);
            g3 = NPF_L(_Fe, _gF6, gS5) NPF_I(i3, iF6, iS5);
            g4 = NPF_R(_Fe, _gF6, gS5) NPF_I(i3, iF6, iS5);
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
            S_LR += (0.006332573977646111*g2*g3*(g5*g7*l1[bb0]*mF6*mF7 - g6*g7*l1[bb1]*mF6*mFe1 + g6*g8*l1[bb0]*mF7*mFe1 - g5*g8*l1[bb1]*Sqr(mFe1)))/((Sqr(mF6) - Sqr(mFe1))*(Sqr(mFe1) - Sqr(mS5)));
            S_RL += (0.006332573977646111*g1*g4*(g6*g8*l1[bb0]*mF6*mF7 - g5*g8*l1[bb1]*mF6*mFe1 + g5*g7*l1[bb0]*mF7*mFe1 - g6*g7*l1[bb1]*Sqr(mFe1)))/((Sqr(mF6) - Sqr(mFe1))*(Sqr(mFe1) - Sqr(mS5)));
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
   }; // End of struct genericSum10_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum10( void ) {
      using GenericKeys = boost::mpl::vector< GenericF6Key, GenericF7Key, GenericS5Key, GenericS8Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, fields::Fv, fields::hh, fields::Hpm>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, fields::Cha, fields::hh, fields::Sv>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, fields::Chi, fields::hh, fields::Se>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, fields::Fv, fields::Ah, fields::Hpm>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, fields::Cha, fields::Ah, fields::Sv>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, fields::Chi, fields::Ah, fields::Se>
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
         genericSum10_impl
         >( *this );
   } // End of function genericSum10()
   
   template<class GenericFieldMap>
   struct genericSum11_impl : generic_sum_base {
      genericSum11_impl( const generic_sum_base &base ) :
      generic_sum_base( base ) {
      } // End of constructor genericSum11_impl
   
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
         using _Fd = typename cxx_diagrams::bar<fields::Fd>::type;
         using _Fe = typename cxx_diagrams::bar<fields::Fe>::type;
         using _gF6 = typename cxx_diagrams::bar<gF6>::type;
         using _gF7 = typename cxx_diagrams::bar<gF7>::type;
         using _gS8 = typename cxx_diagrams::conj<gS8>::type;
         using _gV5 = typename cxx_diagrams::conj<gV5>::type;
         using Fd = fields::Fd;
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
            
            g1 = NPF_L(_Fd, Fd, _gV5) NPF_I(i4, i2, iV5);
            g2 = NPF_R(_Fd, Fd, _gV5) NPF_I(i4, i2, iV5);
            g3 = NPF_L(_Fe, _gF6, gV5) NPF_I(i3, iF6, iV5);
            g4 = NPF_R(_Fe, _gF6, gV5) NPF_I(i3, iF6, iV5);
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
            V_LR += (0.006332573977646111*g2*g3*(-(g5*g7*l1[bb0]*mF6*mF7) + g6*g7*l1[bb1]*mF6*mFe1 - g6*g8*l1[bb0]*mF7*mFe1 + g5*g8*l1[bb1]*Sqr(mFe1)))/((Sqr(mF6) - Sqr(mFe1))*(Sqr(mFe1) - Sqr(mV5)));
            V_RL += (0.006332573977646111*g1*g4*(-(g6*g8*l1[bb0]*mF6*mF7) + g5*g8*l1[bb1]*mF6*mFe1 - g5*g7*l1[bb0]*mF7*mFe1 + g6*g7*l1[bb1]*Sqr(mFe1)))/((Sqr(mF6) - Sqr(mFe1))*(Sqr(mFe1) - Sqr(mV5)));
            V_RR += (0.006332573977646111*g2*g4*(-(g6*g8*l1[bb0]*mF6*mF7) + g5*g8*l1[bb1]*mF6*mFe1 - g5*g7*l1[bb0]*mF7*mFe1 + g6*g7*l1[bb1]*Sqr(mFe1)))/((Sqr(mF6) - Sqr(mFe1))*(Sqr(mFe1) - Sqr(mV5)));
            T_LL += 0;
            T_RR += 0;
         
         }}}} // End of summation over generic fields
   
         return {S_LL, S_LR, S_RL, S_RR, V_LL, V_LR, V_RL, V_RR, T_LL, T_RR};
      } // End of operator()( void )
   }; // End of struct genericSum11_impl<GenericFieldMap>
   
   std::array<std::complex<double>,10> genericSum11( void ) {
      using GenericKeys = boost::mpl::vector< GenericF6Key, GenericF7Key, GenericS8Key, GenericV5Key >;
      using GenericInsertions = boost::mpl::vector<
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, fields::Fv, fields::Hpm, fields::VZ>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, fields::Cha, fields::Sv, fields::VZ>,
         boost::mpl::vector<typename cxx_diagrams::bar<fields::Fe>::type, fields::Chi, fields::Se, fields::VZ>
         >;
      using combinatorial_factors = boost::mpl::vector<
         boost::mpl::int_<1>,
         boost::mpl::int_<1>,
         boost::mpl::int_<1>
         >;
      using color_factors = boost::mpl::vector<
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
         genericSum11_impl
         >( *this );
   } // End of function genericSum11()

   public:
   nPointFeFdFeFd_3917884761161352( const MSSM_mass_eigenstates &model,
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
      const auto genericsum5 = genericSum5();
      const auto genericsum6 = genericSum6();
      const auto genericsum7 = genericSum7();
      const auto genericsum8 = genericSum8();
      const auto genericsum9 = genericSum9();
      const auto genericsum10 = genericSum10();
      const auto genericsum11 = genericSum11();
   
      for ( std::size_t i=0; i<coeffsLength; i++ ) {
         genericSummation[i] += genericsum1[i]+genericsum2[i]+genericsum3[i]+genericsum4[i]+genericsum5[i]+genericsum6[i]+genericsum7[i]+genericsum8[i]+genericsum9[i]+genericsum10[i]+genericsum11[i];
      }
      return genericSummation;
   } // End of calculate()
}; // End of nPointFeFdFeFd_3917884761161352

std::array<std::complex<double>, 10> conversion_FeFd_to_FeFd_All1loop(const MSSM_mass_eigenstates &model,
const std::array<int, 4> &indices,
const std::array<Eigen::Vector4d, 0> &momenta
) {
   nPointFeFdFeFd_3917884761161352 helper{model, indices, momenta};
   return helper.calculate();
}

} // namespace npointfunctions
} // namespace MSSM_cxx_diagrams

namespace {

std::valarray<std::complex<double>> zero(
   int generationIndex1,
   int generationIndex2,
   const MSSM_mass_eigenstates&,
   bool){
   std::valarray<std::complex<double>> res {0.0, 0.0, 0.0, 0.0};
   return res;
}

double get_MSUSY(const MSSM_mass_eigenstates_interface& model)
{
   return Min(model.get_MSd().tail<6>().minCoeff(), model.get_MSu().tail<6>().
      minCoeff(), model.get_MSe().tail<6>().minCoeff(), model.get_MHpm().tail<1>()
      .minCoeff(), model.get_MCha().tail<2>().minCoeff());
}

void run_to_MSUSY(MSSM_mass_eigenstates& model)
{
   const double precision_goal = model.get_precision();
   const int Nmax = static_cast<int>(-std::log10(precision_goal) * 10);
   int it = 0;
   double precision = 0.;
   double MSUSY_old = 0., MSUSY_new = 0.;

   VERBOSE_MSG("MSSM_l_to_l_conversion: iterate to run to MSUSY ...");

   do {
      MSUSY_new = get_MSUSY(model);

      if (MSUSY_new <= 0.)
         return;

      VERBOSE_MSG("MSSM_l_to_l_conversion:    running to new MSUSY = " << MSUSY_new);

      model.run_to(MSUSY_new);
      model.calculate_DRbar_masses();

      precision = MaxRelDiff(MSUSY_old, MSUSY_new);
      MSUSY_old = MSUSY_new;
   } while (precision > precision_goal && ++it < Nmax);

   VERBOSE_MSG("MSSM_l_to_l_conversion: iteration finished. MSUSY = " << MSUSY_new);

   if (precision > precision_goal) {
      WARNING("MSSM_l_to_l_conversion: run_to_MSUSY(): "
              "Iteration did not converge after " << Nmax
              << " iterations (precision goal: " << precision_goal << ")");
   }
}

} // anonymous namespace

using namespace MSSM_cxx_diagrams;
using namespace MSSM_FFV_form_factors;

namespace MSSM_l_to_l_conversion {

struct overlap_integrals {
   double D;
   double Vp;
   double Vn;
   double Sp;
   double Sn;
};

template <class Lepton, class Up, class Down>
struct Mass {
   // PDG 2018 data
   double p = 0.938272081, n = 0.939565413;
   // Other masses
   double l, u, c, d, s, b;
   Mass(const context_base& context, int in) {
      l = context.mass<Lepton>({in});
      u = context.mass<Up>({0});
      c = context.mass<Up>({1});
      d = context.mass<Down>({0});
      s = context.mass<Down>({1});
      b = context.mass<Down>({2});
   }
};

overlap_integrals get_overlap_integrals(const Nucleus N, const softsusy::QedQcd& qedqcd) {
   overlap_integrals res;

   // get muon pole mass from input slha file
   const auto muon_pole_mass_5o2 = pow(qedqcd.displayMass(softsusy::mMuon), 5./2.);

   // Tab. 2 of hep-ph/0203110
   switch (N) {
      case Nucleus::Au:
         res.D  = 0.1670 * muon_pole_mass_5o2;
         res.Vp = 0.0859 * muon_pole_mass_5o2;
         res.Vn = 0.1080 * muon_pole_mass_5o2;
         res.Sp = 0.0523 * muon_pole_mass_5o2;
         res.Sn = 0.0610 * muon_pole_mass_5o2;
         return res;
      case Nucleus::Al:
         res.D  = 0.0357 * muon_pole_mass_5o2;
         res.Vp = 0.0159 * muon_pole_mass_5o2;
         res.Vn = 0.0169 * muon_pole_mass_5o2;
         res.Sp = 0.0153 * muon_pole_mass_5o2;
         res.Sn = 0.0163 * muon_pole_mass_5o2;
         return res;
      default:
         throw std::invalid_argument("Unknown nucleus");
   }
}

double get_capture_rate (const Nucleus N) {
   switch (N) {
      case Nucleus::Au:
         return 8.84868e-18;
      case Nucleus::Al:
         return 4.64079e-19;
      default:
         throw std::invalid_argument("Unknown nucleus");
   }
}

typedef std::valarray<std::complex<double>>  (*ffv_function)
   (int, int, const MSSM_mass_eigenstates&, bool);

typedef std::array<std::complex<double>, 10> (*npf_function)
   (const MSSM_mass_eigenstates&, const std::array<int,4>&, const std::array<Eigen::Vector4d,0>&);

/**
 * @tparam    Q     Type of a quark field.
 * @tparam    A     Type of a photon field.
 * @param[in] model Mass eigenstates.
 * @param[in] g     Generation index for quarks.
 * @return Left part of qqv coupling, multiplied by -i.
 */
template <class Q, class A>
std::complex<double> left(const MSSM_mass_eigenstates& model, int g) {
    context_base context {model};
    using vertex = Vertex<typename Q::lorentz_conjugate, Q, A>;
    std::array<int, 2> indices {g, g};
    const auto value =  vertex::evaluate(indices, context);
    return value.left();
}

/**
 * @tparam    Q     Type of a quark field.
 * @tparam    A     Type of a photon field.
 * @param[in] model Mass eigenstates.
 * @param[in] g     Generation index for quarks.
 * @return Right part of qqv coupling, multiplied by -i.
 */
template <class Q, class A>
std::complex<double> right(const MSSM_mass_eigenstates& model, int g) {
    context_base context {model};
    using vertex = Vertex<typename Q::lorentz_conjugate, Q, A>;
    std::array<int, 2> indices {g, g};
    const auto value =  vertex::evaluate(indices, context);
    return value.right();
}

/**
 * @tparam    Q     Type of a quark field.
 * @tparam    A     Type of a photon field.
 * @tparam    T     Type of a form factors.
 * @param[in] model Mass eigenstates.
 * @param[in] ff    Lepton-photon form factors.
 * @param[in] g     Generation index for quarks.
 * @return Set of four-fermion coefficients from photon penguin amplitudes
 *         (without overall i; with appropriate embedding).
 */
template <class Q, class A, class T>
Eigen::Array<std::complex<double>,10,1> embed_photon(
   const MSSM_mass_eigenstates& model, const T& ff, int g) {
   // Get quark-photon couplings (without i, as everywhere):
   const auto qL =  left<Q, A>(model, g);
   const auto qR = right<Q, A>(model, g);
   Eigen::Array<std::complex<double>,10,1> res{};
   // Term from eq. (t.1) will contribute to four-fermion vector coefficients.
   // Minus comes from the form_factors embedding into four-fermion amplitude:
   res[4] = - ff[0] * qL;
   res[5] = - ff[0] * qR;
   res[6] = - ff[1] * qL;
   res[7] = - ff[1] * qR;
   return res;
};

/**
 * @tparam    Name  Function name for npf function.
 * @param[in] model Mass eigenstates.
 * @param[in] in    Generation index for incoming lepton.
 * @param[in] out   Generation index for outgoing lepton.
 * @param[in] g     Generation index for quarks.
 * @return Set of four-fermion coefficients for non-photonic amplitudes
 *         (without overall i; with fixed signs for tensor operators).
 */
template <npf_function Name>
Eigen::Array<std::complex<double>,10,1> fix_tensors_sign(
   const MSSM_mass_eigenstates& model, int in, int out, int g) {
   const auto npf = Name(model,
         std::array<int,4>{in, g, out, g},
         std::array<Eigen::Vector4d, 0>{});
   Eigen::Array<std::complex<double>,10,1> res(npf.data());
   res[8] = - res[8];
   res[9] = - res[9];
   return res;
};

/**
 * Form factors are defined via the following formula (q = pj - pi > 0; see
 * eq. (3.4) of 1902.06650 (up to electric charge);
 * pi - momenta[going from blob] of outgoing lepton,
 * pj - momenta[going into blob] of incoming lepton):
 *    <pi, q| T exp(i L_MSSM dx)|pj> =
 *    i * ubari
 *          q^2 gamma_mu (A1_X * P_X)                                  (t.1)
 *          + (zero after embedding term)                              (t.2)
 *          + i mj sigma_munu q^nu (A2_X * P_X)                        (t.3)
 *        uj e^*mu
 * @note Form factors below are ordered as A1_L, A1_R, A2_L, A2_R.
 * @note NPF function returns such result, that G4 = i * npf.
 * @note Match to a low energy effective model with the same covariant derivative
 *       for photon, as in MSSM, using eq. (t.3).
 *       Minus comes from the form_factors embedding into four-fermion
 *       amplitude, because we use descending ordering for external fermions
 *          <out:q4,l3| T exp(i L_MSSM dx) |in:q2,l1> =: G4
 * @note For four fermion coefficient matching minus comes from the G4
 *       definition, because
 *          <out:q4,l3| T exp(i L_low dx) |in:q2,l1> = -i * C_XY [3X1]*[4Y2].
 * @tparam    Lepton   Type of a lepton fields.
 * @tparam    Up       Type of a up-quark fields.
 * @tparam    Down     Type of a down-quark fields.
 * @tparam    Photon   Type of a photon field.
 * @tparam    photon   Function name for photon form factors.
 * @tparam    npf_up   Function name for npf function for up-quark contribution.
 * @tparam    npf_down Function name for npf function for down-quark
 *                     contribution.
 * @param[in] in     Generation index for incoming lepton.
 * @param[in] out    Generation index for outgoing lepton.
 * @param[in] model  Mass eigenstates.
 * @param[in] parameters Parameters for the observable calculation.
 * @param[in] qedqcd Reference to low-energy data.
 * @return Set of (l to l conversion) process and Wilson coefficients.
 */
template<class Lepton, class Up, class Down, class Photon,
   ffv_function photon, npf_function npf_up, npf_function npf_down>
Eigen::Array<std::complex<double>, 13, 1> forge_conversion(int in, int out,
   const MSSM_l_to_l_conversion::Nucleus nucleus,
   const MSSM_mass_eigenstates& model_,
   const LToLConversion_settings& parameters,
   const softsusy::QedQcd& qedqcd) {

   MSSM_mass_eigenstates model(model_);
   try {
      run_to_MSUSY(model);
      model.get_physical().clear();
      model.solve_ewsb();
   } catch (const Error& e) {
      ERROR("MSSM_l_to_l_conversion:" << e.what_detailed());
      return std::numeric_limits<Eigen::Array<std::complex<double>,13,1>>::quiet_NaN();
   }

   context_base context {model};

   const auto form_factors = photon(in, out, model, false);
   const auto photon_u = F4F(Up,   0);
   const auto photon_d = F4F(Down, 0);
   const auto photon_s = F4F(Down, 1);

   const auto npf_u = NPF(npf_up,   0);
   const auto npf_d = NPF(npf_down, 0);
   const auto npf_s = NPF(npf_down, 1);

   // Matching
   const auto DL = - 0.5 * form_factors[2];
   const auto DR = - 0.5 * form_factors[3];
   const auto CXYu = - photon_u - npf_u;
   const auto CXYd = - photon_d - npf_d;
   const auto CXYs = - photon_s - npf_d;

   auto CSLu = ( CXYu[0] + CXYu[1] )/2.;
   auto CSRu = ( CXYu[2] + CXYu[3] )/2.;
   auto CSLd = ( CXYd[0] + CXYd[1] )/2.;
   auto CSRd = ( CXYd[2] + CXYd[3] )/2.;
   auto CSLs = ( CXYs[0] + CXYs[1] )/2.;
   auto CSRs = ( CXYs[2] + CXYs[3] )/2.;

   auto CVLu = ( CXYu[4] + CXYu[5] )/2.;
   auto CVRu = ( CXYu[6] + CXYu[7] )/2.;
   auto CVLd = ( CXYd[4] + CXYd[5] )/2.;
   auto CVRd = ( CXYd[6] + CXYd[7] )/2.;

   auto CTLu = CXYu[8], CTRu = CXYu[9];
   auto CTLd = CXYd[8], CTRd = CXYd[9];
   auto CTLs = CXYs[8], CTRs = CXYs[9];

   // Vector contribution.
   auto gpLV = ( GV(p,u)*CVLu + GV(p,d)*CVLd );
   auto gpRV = ( GV(p,u)*CVRu + GV(p,d)*CVRd );
   auto gnLV = ( GV(n,u)*CVLu + GV(n,d)*CVLd );
   auto gnRV = ( GV(n,u)*CVRu + GV(n,d)*CVRd );

   // Scalar contribution from scalar coefficients.
   Mass<Lepton, Up, Down> m(context, in);
   auto gpLS = ( GS(p,u)*CSLu + GS(p,d)*CSLd + GS(p,s)*CSLs );
   auto gpRS = ( GS(p,u)*CSRu + GS(p,d)*CSRd + GS(p,s)*CSRs );
   auto gnLS = ( GS(n,u)*CSLu + GS(n,d)*CSLd + GS(n,s)*CSLs );
   auto gnRS = ( GS(n,u)*CSRu + GS(n,d)*CSRd + GS(n,s)*CSRs );

   if (parameters.get(LToLConversion_settings::include_tensor_contribution)) {
      gpLS += (m.l/m.p)*(GT(p,u)*CTLu + GT(p,d)*CTLd + GT(p,s)*CTLs);
      gpRS += (m.l/m.p)*(GT(p,u)*CTRu + GT(p,d)*CTRd + GT(p,s)*CTRs);
      gnLS += (m.l/m.n)*(GT(n,u)*CTLu + GT(n,d)*CTLd + GT(n,s)*CTLs);
      gnRS += (m.l/m.n)*(GT(n,u)*CTRu + GT(n,d)*CTRd + GT(n,s)*CTRs);
   }

   if (parameters.get(LToLConversion_settings::include_gluonic_contribution)) {
      // Only scalar contributions is needed - no photon; minus from matching.
      const auto CXYc = - NPF(npf_up,   1);
      const auto CXYb = - NPF(npf_down, 2);

      auto CSLc = ( CXYc[0] + CXYc[1] )/2.;
      auto CSRc = ( CXYc[2] + CXYc[3] )/2.;
      auto CSLb = ( CXYb[0] + CXYb[1] )/2.;
      auto CSRb = ( CXYb[2] + CXYb[3] )/2.;

      const double Ggp = - 8.*Pi/9 * (1 - m.u/m.p * GS(p,u)
                                        - m.d/m.p * GS(p,d)
                                        - m.s/m.p * GS(p,s));
      const double Ggn = - 8.*Pi/9 * (1 - m.u/m.n * GS(n,u)
                                        - m.d/m.n * GS(n,d)
                                        - m.s/m.n * GS(n,s));

      // Scalar contribution from gluonic coefficients.
      gpLS += - 1./(12*Pi) * m.p*Ggp*(CSLc/m.c + CSLb/m.b);
      gpRS += - 1./(12*Pi) * m.p*Ggp*(CSRc/m.c + CSRb/m.b);
      gnLS += - 1./(12*Pi) * m.n*Ggn*(CSLc/m.c + CSLb/m.b);
      gnRS += - 1./(12*Pi) * m.n*Ggn*(CSRc/m.c + CSRb/m.b);
   }

   const auto ff = get_overlap_integrals(nucleus, qedqcd);
   const auto left  = DR*ff.D/4. - gpLV*ff.Vp - gnLV*ff.Vn - gpRS*ff.Sp - gnRS*ff.Sn;
   const auto right = DL*ff.D/4. - gpRV*ff.Vp - gnRV*ff.Vn - gpLS*ff.Sp - gnLS*ff.Sn;

   // eq. (3.31) of 1902.06650 (with some coefficient redefinitions)
   const double conversion_rate = 4*(std::norm(left) + std::norm(right));
   const double capture_rate = get_capture_rate(nucleus);
   Eigen::Array<std::complex<double>,13,1> res;
   res << conversion_rate/capture_rate,
          DL,      DR,
          CXYu[0], CXYu[1], CXYu[2], CXYu[3],
          CXYu[4], CXYu[5], CXYu[6], CXYu[7],
          CXYu[8], CXYu[9];
   return res;
}

Eigen::Array<std::complex<double>,13,1> calculate_FeFe_forAll_1loop(int in, int out, const MSSM_l_to_l_conversion::Nucleus n, const MSSM_mass_eigenstates& model, const LToLConversion_settings& parameters, const softsusy::QedQcd& qedqcd) {
         return forge_conversion<
            fields::Fe, fields::Fu, fields::Fd, fields::VP,
            calculate_Fe_Fe_VP_form_factors,
            npointfunctions::conversion_FeFu_to_FeFu_All1loop,
            npointfunctions::conversion_FeFd_to_FeFd_All1loop
         >(in, out, n, model, parameters, qedqcd);
}

} // namespace MSSM_l_to_l_conversion
} // namespace flexiblesusy
