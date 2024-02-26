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
 * @file cxx_qft/MSSMNoFVHimalaya_vertices.hpp
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#ifndef MSSMNoFVHimalaya_CXXQFT_VERTICES_H
#define MSSMNoFVHimalaya_CXXQFT_VERTICES_H

#include "multiindex.hpp"
#include "numerics2.hpp"

#include "MSSMNoFVHimalaya_fields.hpp"
#include "cxx_qft/vertices.hpp"

#include <array>
#include <algorithm>

#include <boost/mpl/erase.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/joint_view.hpp>
#include <boost/mpl/vector.hpp>

namespace flexiblesusy {
namespace MSSMNoFVHimalaya_cxx_diagrams {

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
template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Cha>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Cha>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fb>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Fb>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fc>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Fc>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fd>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Fd>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fe>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Fm>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fs>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Fs>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Ftau>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Ftau>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Ft>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Ft>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fu>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Fu>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Hpm>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Hpm>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Sb>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Sb>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Sc>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Sc>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Sd>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Sd>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Se>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Se>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Sm>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Sm>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Ss>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Ss>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Stau>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Stau>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::St>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::St>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::Ah, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Su>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Su>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::Chi, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Sm>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Fm>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::Fm, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::Fm, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::VZ>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Cha>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Cha>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fb>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Fb>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fc>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Fc>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fd>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Fd>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fe>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Fm>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fs>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Fs>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Ftau>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Ftau>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Ft>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Ft>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fu>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Fu>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Hpm>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Hpm>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Sb>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Sb>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Sc>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Sc>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Sd>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Sd>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Se>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Se>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Sm>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Sm>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Ss>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Ss>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Stau>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Stau>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::St>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::St>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Su>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Su>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::hh, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::VWm>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::VWm>
{
   static cxx_diagrams::InverseMetricVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VP, MSSMNoFVHimalaya_cxx_diagrams::fields::Cha, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Cha>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VP, MSSMNoFVHimalaya_cxx_diagrams::fields::Fb, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fb>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VP, MSSMNoFVHimalaya_cxx_diagrams::fields::Fc, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fc>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VP, MSSMNoFVHimalaya_cxx_diagrams::fields::Fd, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fd>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VP, MSSMNoFVHimalaya_cxx_diagrams::fields::Fe, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fe>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VP, MSSMNoFVHimalaya_cxx_diagrams::fields::Fm, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VP, MSSMNoFVHimalaya_cxx_diagrams::fields::Fs, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fs>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VP, MSSMNoFVHimalaya_cxx_diagrams::fields::Ftau, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Ftau>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VP, MSSMNoFVHimalaya_cxx_diagrams::fields::Ft, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Ft>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VP, MSSMNoFVHimalaya_cxx_diagrams::fields::Fu, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fu>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VP, MSSMNoFVHimalaya_cxx_diagrams::fields::Hpm, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Hpm>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VP, MSSMNoFVHimalaya_cxx_diagrams::fields::Sb, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Sb>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VP, MSSMNoFVHimalaya_cxx_diagrams::fields::Sc, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Sc>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VP, MSSMNoFVHimalaya_cxx_diagrams::fields::Sd, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Sd>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VP, MSSMNoFVHimalaya_cxx_diagrams::fields::Se, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Se>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VP, MSSMNoFVHimalaya_cxx_diagrams::fields::Sm, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Sm>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VP, MSSMNoFVHimalaya_cxx_diagrams::fields::Ss, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Ss>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VP, MSSMNoFVHimalaya_cxx_diagrams::fields::Stau, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Stau>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VP, MSSMNoFVHimalaya_cxx_diagrams::fields::St, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::St>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VP, MSSMNoFVHimalaya_cxx_diagrams::fields::Su, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Su>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VP, MSSMNoFVHimalaya_cxx_diagrams::fields::VWm, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::VWm>::type>
{
   static cxx_diagrams::TripleVectorVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VP, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Fm>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VZ, MSSMNoFVHimalaya_cxx_diagrams::fields::Cha, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Cha>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VZ, MSSMNoFVHimalaya_cxx_diagrams::fields::Fb, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fb>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VZ, MSSMNoFVHimalaya_cxx_diagrams::fields::Fc, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fc>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VZ, MSSMNoFVHimalaya_cxx_diagrams::fields::Fd, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fd>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VZ, MSSMNoFVHimalaya_cxx_diagrams::fields::Fe, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fe>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VZ, MSSMNoFVHimalaya_cxx_diagrams::fields::Fm, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VZ, MSSMNoFVHimalaya_cxx_diagrams::fields::Fs, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fs>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VZ, MSSMNoFVHimalaya_cxx_diagrams::fields::Ftau, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Ftau>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VZ, MSSMNoFVHimalaya_cxx_diagrams::fields::Ft, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Ft>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VZ, MSSMNoFVHimalaya_cxx_diagrams::fields::Fu, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fu>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<MSSMNoFVHimalaya_cxx_diagrams::fields::VZ, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Fm>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Cha>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Cha, MSSMNoFVHimalaya_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Cha, MSSMNoFVHimalaya_cxx_diagrams::fields::SvmL>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Fm, MSSMNoFVHimalaya_cxx_diagrams::fields::Ah>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Fm, MSSMNoFVHimalaya_cxx_diagrams::fields::hh>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Fm, MSSMNoFVHimalaya_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Fm, MSSMNoFVHimalaya_cxx_diagrams::fields::VZ>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Hpm, MSSMNoFVHimalaya_cxx_diagrams::fields::Fvm>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Sm, MSSMNoFVHimalaya_cxx_diagrams::fields::Chi>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fm>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::VWm, MSSMNoFVHimalaya_cxx_diagrams::fields::Fvm>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fvm>::type, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Hpm>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Fm>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Fvm>::type, typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::VWm>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Fm>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Hpm>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Hpm, MSSMNoFVHimalaya_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Hpm>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::VWm, MSSMNoFVHimalaya_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::InverseMetricVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::Sm>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Sm, MSSMNoFVHimalaya_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::SvmL>::type, typename MSSMNoFVHimalaya_cxx_diagrams::fields::bar<MSSMNoFVHimalaya_cxx_diagrams::fields::Cha>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Fm>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::VWm>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::Hpm, MSSMNoFVHimalaya_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::InverseMetricVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<typename MSSMNoFVHimalaya_cxx_diagrams::fields::conj<MSSMNoFVHimalaya_cxx_diagrams::fields::VWm>::type, MSSMNoFVHimalaya_cxx_diagrams::fields::VWm, MSSMNoFVHimalaya_cxx_diagrams::fields::VP>
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

} // namespace MSSMNoFVHimalaya_cxx_diagrams
} // namespace flexiblesusy

#endif
