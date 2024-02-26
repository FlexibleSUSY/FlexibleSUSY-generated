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
 * @file cxx_qft/CE6SSM_vertices.hpp
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#ifndef CE6SSM_CXXQFT_VERTICES_H
#define CE6SSM_CXXQFT_VERTICES_H

#include "multiindex.hpp"
#include "numerics2.hpp"

#include "CE6SSM_fields.hpp"
#include "cxx_qft/vertices.hpp"

#include <array>
#include <algorithm>

#include <boost/mpl/erase.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/joint_view.hpp>
#include <boost/mpl/vector.hpp>

namespace flexiblesusy {
namespace CE6SSM_cxx_diagrams {

namespace detail {
template<class... Fields> struct VertexImpl;
} // namespace detail

template <class... Fields>
struct Vertex {
   using index_bounds = typename boost::mpl::fold<
      boost::mpl::vector<Fields...>,
      boost::mpl::pair<boost::mpl::vector<>, boost::mpl::vector<>>,
      boost::mpl::pair<
         boost::mpl::joint_view<
            boost::mpl::first<boost::mpl::_1>,
            boost::mpl::first<meta::index_bounds<boost::mpl::_2>>
         >,
         boost::mpl::joint_view<
            boost::mpl::second<boost::mpl::_1>,
            boost::mpl::second<meta::index_bounds<boost::mpl::_2>>
         >
      >
   >::type;
   using indices_type = std::array<int,
      cxx_diagrams::detail::total_number_of_field_indices<
         boost::mpl::vector<Fields...>
      >::value
   >;
   using vertex_type = decltype(
      detail::VertexImpl<Fields...>::evaluate(
         std::declval<indices_type>(),
         std::declval<context_base>()
      )
   );

   template <int FieldIndex>
   static typename field_indices<typename boost::mpl::at_c<
      boost::mpl::vector<Fields...>, FieldIndex>::type
   >::type indices_of_field(const indices_type& indices)
   {
      using namespace boost::mpl;
      using fields = vector<Fields...>;

      using result_type = typename field_indices<
         typename boost::mpl::at_c<fields, FieldIndex>::type
      >::type;

      using preceeding_fields = typename erase<fields,
         typename advance<
            typename begin<fields>::type,
            int_<FieldIndex>
         >::type,
         typename end<fields>::type
      >::type;

      constexpr int offset =
         cxx_diagrams::detail::total_number_of_field_indices<preceeding_fields>::value;
      constexpr int length = std::tuple_size<result_type>::value;

      result_type result_indices;
      std::copy(indices.begin() + offset,
         indices.begin() + offset + length,
         result_indices.begin()
      );

      return result_indices;
   }

   static vertex_type
   evaluate(const indices_type& indices, const context_base& context)
   {
      return detail::VertexImpl<Fields...>::evaluate(indices, context);
   }
};

struct context_base;

namespace detail {
template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::Ah, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::ChaI>::type, CE6SSM_cxx_diagrams::fields::ChaI>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::Ah, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Cha>::type, CE6SSM_cxx_diagrams::fields::Cha>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::Ah, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fd>::type, CE6SSM_cxx_diagrams::fields::Fd>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::Ah, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::FDX>::type, CE6SSM_cxx_diagrams::fields::FDX>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::Ah, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fe>::type, CE6SSM_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::Ah, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fu>::type, CE6SSM_cxx_diagrams::fields::Fu>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::Ah, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Hpm>::type, CE6SSM_cxx_diagrams::fields::Hpm>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::Ah, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Sd>::type, CE6SSM_cxx_diagrams::fields::Sd>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::Ah, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SDX>::type, CE6SSM_cxx_diagrams::fields::SDX>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::Ah, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Se>::type, CE6SSM_cxx_diagrams::fields::Se>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::Ah, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHIp>::type, CE6SSM_cxx_diagrams::fields::SHIp>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::Ah, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Su>::type, CE6SSM_cxx_diagrams::fields::Su>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::Chi, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Se>::type, CE6SSM_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::Fe, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fe>::type, CE6SSM_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::Fe, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fe>::type, CE6SSM_cxx_diagrams::fields::VZ>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::hh, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::ChaI>::type, CE6SSM_cxx_diagrams::fields::ChaI>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::hh, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Cha>::type, CE6SSM_cxx_diagrams::fields::Cha>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::hh, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fd>::type, CE6SSM_cxx_diagrams::fields::Fd>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::hh, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::FDX>::type, CE6SSM_cxx_diagrams::fields::FDX>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::hh, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fe>::type, CE6SSM_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::hh, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fu>::type, CE6SSM_cxx_diagrams::fields::Fu>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::hh, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Hpm>::type, CE6SSM_cxx_diagrams::fields::Hpm>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::hh, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Sd>::type, CE6SSM_cxx_diagrams::fields::Sd>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::hh, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SDX>::type, CE6SSM_cxx_diagrams::fields::SDX>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::hh, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Se>::type, CE6SSM_cxx_diagrams::fields::Se>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::hh, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHIp>::type, CE6SSM_cxx_diagrams::fields::SHIp>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::hh, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHpp>::type, CE6SSM_cxx_diagrams::fields::SHpp>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::hh, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Su>::type, CE6SSM_cxx_diagrams::fields::Su>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::hh, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::VWm>::type, CE6SSM_cxx_diagrams::fields::VWm>
{
   static cxx_diagrams::InverseMetricVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::VP, CE6SSM_cxx_diagrams::fields::ChaI, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::ChaI>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::VP, CE6SSM_cxx_diagrams::fields::Cha, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Cha>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::VP, CE6SSM_cxx_diagrams::fields::Fd, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fd>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::VP, CE6SSM_cxx_diagrams::fields::FDX, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::FDX>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::VP, CE6SSM_cxx_diagrams::fields::Fe, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fe>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::VP, CE6SSM_cxx_diagrams::fields::Fu, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fu>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::VP, CE6SSM_cxx_diagrams::fields::Hpm, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Hpm>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::VP, CE6SSM_cxx_diagrams::fields::Sd, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Sd>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::VP, CE6SSM_cxx_diagrams::fields::SDX, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SDX>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::VP, CE6SSM_cxx_diagrams::fields::Se, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Se>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::VP, CE6SSM_cxx_diagrams::fields::SHIp, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHIp>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::VP, CE6SSM_cxx_diagrams::fields::SHpp, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHpp>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::VP, CE6SSM_cxx_diagrams::fields::Su, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Su>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::VP, CE6SSM_cxx_diagrams::fields::VWm, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::VWm>::type>
{
   static cxx_diagrams::TripleVectorVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::VP, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fe>::type, CE6SSM_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::VZ, CE6SSM_cxx_diagrams::fields::ChaI, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::ChaI>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::VZ, CE6SSM_cxx_diagrams::fields::Cha, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Cha>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::VZ, CE6SSM_cxx_diagrams::fields::Fd, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fd>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::VZ, CE6SSM_cxx_diagrams::fields::FDX, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::FDX>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::VZ, CE6SSM_cxx_diagrams::fields::Fe, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fe>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::VZ, CE6SSM_cxx_diagrams::fields::Fu, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fu>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::VZp, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fe>::type, CE6SSM_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CE6SSM_cxx_diagrams::fields::VZ, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fe>::type, CE6SSM_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Cha>::type, CE6SSM_cxx_diagrams::fields::Cha, CE6SSM_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fe>::type, CE6SSM_cxx_diagrams::fields::Cha, CE6SSM_cxx_diagrams::fields::Sv>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fe>::type, CE6SSM_cxx_diagrams::fields::Fe, CE6SSM_cxx_diagrams::fields::Ah>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fe>::type, CE6SSM_cxx_diagrams::fields::Fe, CE6SSM_cxx_diagrams::fields::hh>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fe>::type, CE6SSM_cxx_diagrams::fields::Fe, CE6SSM_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fe>::type, CE6SSM_cxx_diagrams::fields::Fe, CE6SSM_cxx_diagrams::fields::VZp>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fe>::type, CE6SSM_cxx_diagrams::fields::Fe, CE6SSM_cxx_diagrams::fields::VZ>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fe>::type, CE6SSM_cxx_diagrams::fields::Hpm, CE6SSM_cxx_diagrams::fields::Fv>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fe>::type, CE6SSM_cxx_diagrams::fields::Se, CE6SSM_cxx_diagrams::fields::Chi>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fe>::type, CE6SSM_cxx_diagrams::fields::VWm, CE6SSM_cxx_diagrams::fields::Fv>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fv>::type, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Hpm>::type, CE6SSM_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fv>::type, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::VWm>::type, CE6SSM_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Hpm>::type, CE6SSM_cxx_diagrams::fields::Hpm, CE6SSM_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Hpm>::type, CE6SSM_cxx_diagrams::fields::VWm, CE6SSM_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::InverseMetricVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Se>::type, CE6SSM_cxx_diagrams::fields::Se, CE6SSM_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Sv>::type, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Cha>::type, CE6SSM_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::VWm>::type, CE6SSM_cxx_diagrams::fields::Hpm, CE6SSM_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::InverseMetricVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::VWm>::type, CE6SSM_cxx_diagrams::fields::VWm, CE6SSM_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::TripleVectorVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};



cxx_diagrams::ChiralVertex unit_charge(const context_base& context);
} // namespace detail

inline double unit_charge(const context_base& context)
{
   return -(detail::unit_charge(context).left().real() /
            fields::Electron::electricCharge);
}

} // namespace CE6SSM_cxx_diagrams
} // namespace flexiblesusy

#endif
