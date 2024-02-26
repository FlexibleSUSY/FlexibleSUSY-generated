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
 * @file cxx_qft/MRSSM2_vertices.hpp
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#ifndef MRSSM2_CXXQFT_VERTICES_H
#define MRSSM2_CXXQFT_VERTICES_H

#include "multiindex.hpp"
#include "numerics2.hpp"

#include "MRSSM2_fields.hpp"
#include "cxx_qft/vertices.hpp"

#include <array>
#include <algorithm>

#include <boost/mpl/erase.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/joint_view.hpp>
#include <boost/mpl/vector.hpp>

namespace flexiblesusy {
namespace MRSSM2_cxx_diagrams {

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
template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::Hpm, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Hpm>::type>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::Se, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Se>::type>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::Sv, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Sv>::type>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha1>::type, MRSSM2_cxx_diagrams::fields::Cha1>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha2>::type, MRSSM2_cxx_diagrams::fields::Cha2>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fd>::type, MRSSM2_cxx_diagrams::fields::Fd>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type, MRSSM2_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fu>::type, MRSSM2_cxx_diagrams::fields::Fu>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Hpm>::type, MRSSM2_cxx_diagrams::fields::Hpm>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Sd>::type, MRSSM2_cxx_diagrams::fields::Sd>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Se>::type, MRSSM2_cxx_diagrams::fields::Se>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type, MRSSM2_cxx_diagrams::fields::SRdp>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRum>::type, MRSSM2_cxx_diagrams::fields::SRum>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Su>::type, MRSSM2_cxx_diagrams::fields::Su>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Cha1, MRSSM2_cxx_diagrams::fields::Fd, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Su>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Cha1, MRSSM2_cxx_diagrams::fields::Fe, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Sv>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Cha1, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha1>::type, MRSSM2_cxx_diagrams::fields::Ah>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Cha1, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha1>::type, MRSSM2_cxx_diagrams::fields::hh>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Cha1, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha1>::type, MRSSM2_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Cha1, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha1>::type, MRSSM2_cxx_diagrams::fields::VZ>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Chi, MRSSM2_cxx_diagrams::fields::Fd, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Sd>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Chi, MRSSM2_cxx_diagrams::fields::Fe, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Se>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Chi, MRSSM2_cxx_diagrams::fields::Fu, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Su>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Chi, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type, MRSSM2_cxx_diagrams::fields::Ah>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Chi, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type, MRSSM2_cxx_diagrams::fields::hh>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Chi, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type, MRSSM2_cxx_diagrams::fields::VZ>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Chi, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type, MRSSM2_cxx_diagrams::fields::Se>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Chi, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Se>::type, MRSSM2_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Fd, MRSSM2_cxx_diagrams::fields::Chi, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Sd>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Fd, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Sd>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Fe, MRSSM2_cxx_diagrams::fields::Chi, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Se>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Fe, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Se>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Fe, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type, MRSSM2_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Fe, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type, MRSSM2_cxx_diagrams::fields::VZ>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Fu, MRSSM2_cxx_diagrams::fields::Chi, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Su>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Fu, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha1>::type, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Sd>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Fu, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Su>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::Hpm, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Hpm>::type>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::Se, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Se>::type>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::Sv, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Sv>::type>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha1>::type, MRSSM2_cxx_diagrams::fields::Cha1>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha2>::type, MRSSM2_cxx_diagrams::fields::Cha2>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fd>::type, MRSSM2_cxx_diagrams::fields::Fd>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type, MRSSM2_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fu>::type, MRSSM2_cxx_diagrams::fields::Fu>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Hpm>::type, MRSSM2_cxx_diagrams::fields::Hpm>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Sd>::type, MRSSM2_cxx_diagrams::fields::Sd>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Se>::type, MRSSM2_cxx_diagrams::fields::Se>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type, MRSSM2_cxx_diagrams::fields::SRdp>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRum>::type, MRSSM2_cxx_diagrams::fields::SRum>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Su>::type, MRSSM2_cxx_diagrams::fields::Su>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::VWm>::type, MRSSM2_cxx_diagrams::fields::VWm>
{
   static cxx_diagrams::InverseMetricVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Hpm, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Hpm>::type, MRSSM2_cxx_diagrams::fields::VZ>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Se, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Se>::type, MRSSM2_cxx_diagrams::fields::VZ>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::Sv, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Sv>::type, MRSSM2_cxx_diagrams::fields::VZ>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::VP, MRSSM2_cxx_diagrams::fields::Cha1, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha1>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::VP, MRSSM2_cxx_diagrams::fields::Cha2, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha2>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::VP, MRSSM2_cxx_diagrams::fields::Fd, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fd>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::VP, MRSSM2_cxx_diagrams::fields::Fe, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::VP, MRSSM2_cxx_diagrams::fields::Fu, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fu>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::VP, MRSSM2_cxx_diagrams::fields::Hpm, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Hpm>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::VP, MRSSM2_cxx_diagrams::fields::Sd, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Sd>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::VP, MRSSM2_cxx_diagrams::fields::Se, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Se>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::VP, MRSSM2_cxx_diagrams::fields::SRdp, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::VP, MRSSM2_cxx_diagrams::fields::SRum, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRum>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::VP, MRSSM2_cxx_diagrams::fields::Su, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Su>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::VP, MRSSM2_cxx_diagrams::fields::VWm, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::VWm>::type>
{
   static cxx_diagrams::TripleVectorVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::VP, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type, MRSSM2_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::VZ, MRSSM2_cxx_diagrams::fields::Cha1, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha1>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::VZ, MRSSM2_cxx_diagrams::fields::Cha2, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha2>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::VZ, MRSSM2_cxx_diagrams::fields::Fd, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fd>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::VZ, MRSSM2_cxx_diagrams::fields::Fe, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::VZ, MRSSM2_cxx_diagrams::fields::Fu, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fu>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MRSSM2_cxx_diagrams::fields::VZ, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type, MRSSM2_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha1>::type, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type, MRSSM2_cxx_diagrams::fields::Sv>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type, MRSSM2_cxx_diagrams::fields::Chi, MRSSM2_cxx_diagrams::fields::Ah>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type, MRSSM2_cxx_diagrams::fields::Chi, MRSSM2_cxx_diagrams::fields::hh>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type, MRSSM2_cxx_diagrams::fields::Chi, MRSSM2_cxx_diagrams::fields::VZ>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type, MRSSM2_cxx_diagrams::fields::Fd, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Sd>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type, MRSSM2_cxx_diagrams::fields::Fe, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Se>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type, MRSSM2_cxx_diagrams::fields::Fu, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Su>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type, MRSSM2_cxx_diagrams::fields::Se>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Se>::type, MRSSM2_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fd>::type, MRSSM2_cxx_diagrams::fields::Chi, MRSSM2_cxx_diagrams::fields::Sd>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fd>::type, MRSSM2_cxx_diagrams::fields::Fd, MRSSM2_cxx_diagrams::fields::Ah>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fd>::type, MRSSM2_cxx_diagrams::fields::Fd, MRSSM2_cxx_diagrams::fields::hh>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fd>::type, MRSSM2_cxx_diagrams::fields::Fd, MRSSM2_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fd>::type, MRSSM2_cxx_diagrams::fields::Fd, MRSSM2_cxx_diagrams::fields::VZ>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fd>::type, MRSSM2_cxx_diagrams::fields::Fu, MRSSM2_cxx_diagrams::fields::Hpm>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fd>::type, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha1>::type, MRSSM2_cxx_diagrams::fields::Su>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fd>::type, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type, MRSSM2_cxx_diagrams::fields::Sd>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type, MRSSM2_cxx_diagrams::fields::Chi, MRSSM2_cxx_diagrams::fields::Se>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type, MRSSM2_cxx_diagrams::fields::Fe, MRSSM2_cxx_diagrams::fields::Ah>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type, MRSSM2_cxx_diagrams::fields::Fe, MRSSM2_cxx_diagrams::fields::hh>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type, MRSSM2_cxx_diagrams::fields::Fe, MRSSM2_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type, MRSSM2_cxx_diagrams::fields::Fe, MRSSM2_cxx_diagrams::fields::VZ>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type, MRSSM2_cxx_diagrams::fields::Fv, MRSSM2_cxx_diagrams::fields::Hpm>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type, MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::Fv>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type, MRSSM2_cxx_diagrams::fields::Se, MRSSM2_cxx_diagrams::fields::Chi>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type, MRSSM2_cxx_diagrams::fields::Se, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type, MRSSM2_cxx_diagrams::fields::VWm, MRSSM2_cxx_diagrams::fields::Fv>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha1>::type, MRSSM2_cxx_diagrams::fields::Sv>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type, MRSSM2_cxx_diagrams::fields::Se>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fu>::type, MRSSM2_cxx_diagrams::fields::Cha1, MRSSM2_cxx_diagrams::fields::Sd>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fu>::type, MRSSM2_cxx_diagrams::fields::Chi, MRSSM2_cxx_diagrams::fields::Su>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fu>::type, MRSSM2_cxx_diagrams::fields::Fd, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Hpm>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fu>::type, MRSSM2_cxx_diagrams::fields::Fu, MRSSM2_cxx_diagrams::fields::Ah>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fu>::type, MRSSM2_cxx_diagrams::fields::Fu, MRSSM2_cxx_diagrams::fields::hh>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fu>::type, MRSSM2_cxx_diagrams::fields::Fu, MRSSM2_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fu>::type, MRSSM2_cxx_diagrams::fields::Fu, MRSSM2_cxx_diagrams::fields::VZ>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fu>::type, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type, MRSSM2_cxx_diagrams::fields::Su>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fv>::type, MRSSM2_cxx_diagrams::fields::Fe, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Hpm>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fv>::type, MRSSM2_cxx_diagrams::fields::Fv, MRSSM2_cxx_diagrams::fields::VZ>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fv>::type, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Hpm>::type, MRSSM2_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fv>::type, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::VWm>::type, MRSSM2_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Hpm>::type, MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Hpm>::type, MRSSM2_cxx_diagrams::fields::VWm, MRSSM2_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::InverseMetricVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Se>::type, MRSSM2_cxx_diagrams::fields::Se, MRSSM2_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Sv>::type, MRSSM2_cxx_diagrams::fields::Cha1, MRSSM2_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::VWm>::type, MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::InverseMetricVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::VWm>::type, MRSSM2_cxx_diagrams::fields::VWm, MRSSM2_cxx_diagrams::fields::VP>
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

} // namespace MRSSM2_cxx_diagrams
} // namespace flexiblesusy

#endif
