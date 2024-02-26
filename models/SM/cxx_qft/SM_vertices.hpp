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
 * @file cxx_qft/SM_vertices.hpp
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#ifndef SM_CXXQFT_VERTICES_H
#define SM_CXXQFT_VERTICES_H

#include "multiindex.hpp"
#include "numerics2.hpp"

#include "SM_fields.hpp"
#include "cxx_qft/vertices.hpp"

#include <array>
#include <algorithm>

#include <boost/mpl/erase.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/joint_view.hpp>
#include <boost/mpl/vector.hpp>

namespace flexiblesusy {
namespace SM_cxx_diagrams {

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
template<> struct VertexImpl<SM_cxx_diagrams::fields::Ah, typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fd>::type, SM_cxx_diagrams::fields::Fd>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<SM_cxx_diagrams::fields::Ah, typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fe>::type, SM_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<SM_cxx_diagrams::fields::Ah, typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fu>::type, SM_cxx_diagrams::fields::Fu>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<SM_cxx_diagrams::fields::Fe, typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fe>::type, SM_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<SM_cxx_diagrams::fields::Fe, typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fe>::type, SM_cxx_diagrams::fields::VZ>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<SM_cxx_diagrams::fields::hh, typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fd>::type, SM_cxx_diagrams::fields::Fd>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<SM_cxx_diagrams::fields::hh, typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fe>::type, SM_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<SM_cxx_diagrams::fields::hh, typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fu>::type, SM_cxx_diagrams::fields::Fu>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<SM_cxx_diagrams::fields::hh, typename SM_cxx_diagrams::fields::conj<SM_cxx_diagrams::fields::Hp>::type, SM_cxx_diagrams::fields::Hp>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<SM_cxx_diagrams::fields::hh, typename SM_cxx_diagrams::fields::conj<SM_cxx_diagrams::fields::VWp>::type, SM_cxx_diagrams::fields::VWp>
{
   static cxx_diagrams::InverseMetricVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<SM_cxx_diagrams::fields::Hp, typename SM_cxx_diagrams::fields::conj<SM_cxx_diagrams::fields::Hp>::type, SM_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<SM_cxx_diagrams::fields::Hp, typename SM_cxx_diagrams::fields::conj<SM_cxx_diagrams::fields::VWp>::type, SM_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::InverseMetricVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<SM_cxx_diagrams::fields::VP, SM_cxx_diagrams::fields::Fd, typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fd>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<SM_cxx_diagrams::fields::VP, SM_cxx_diagrams::fields::Fe, typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fe>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<SM_cxx_diagrams::fields::VP, SM_cxx_diagrams::fields::Fu, typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fu>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<SM_cxx_diagrams::fields::VP, SM_cxx_diagrams::fields::Hp, typename SM_cxx_diagrams::fields::conj<SM_cxx_diagrams::fields::Hp>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<SM_cxx_diagrams::fields::VP, SM_cxx_diagrams::fields::VWp, typename SM_cxx_diagrams::fields::conj<SM_cxx_diagrams::fields::VWp>::type>
{
   static cxx_diagrams::TripleVectorVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<SM_cxx_diagrams::fields::VP, typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fe>::type, SM_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<SM_cxx_diagrams::fields::VWp, typename SM_cxx_diagrams::fields::conj<SM_cxx_diagrams::fields::Hp>::type, SM_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::InverseMetricVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<SM_cxx_diagrams::fields::VWp, typename SM_cxx_diagrams::fields::conj<SM_cxx_diagrams::fields::VWp>::type, SM_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::TripleVectorVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<SM_cxx_diagrams::fields::VZ, SM_cxx_diagrams::fields::Fd, typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fd>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<SM_cxx_diagrams::fields::VZ, SM_cxx_diagrams::fields::Fe, typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fe>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<SM_cxx_diagrams::fields::VZ, SM_cxx_diagrams::fields::Fu, typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fu>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<SM_cxx_diagrams::fields::VZ, typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fe>::type, SM_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fd>::type, SM_cxx_diagrams::fields::Fd, SM_cxx_diagrams::fields::Ah>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fd>::type, SM_cxx_diagrams::fields::Fd, SM_cxx_diagrams::fields::hh>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fd>::type, SM_cxx_diagrams::fields::Fd, SM_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fd>::type, SM_cxx_diagrams::fields::Fd, SM_cxx_diagrams::fields::VZ>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fd>::type, SM_cxx_diagrams::fields::Fu, typename SM_cxx_diagrams::fields::conj<SM_cxx_diagrams::fields::Hp>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fe>::type, SM_cxx_diagrams::fields::Fe, SM_cxx_diagrams::fields::Ah>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fe>::type, SM_cxx_diagrams::fields::Fe, SM_cxx_diagrams::fields::hh>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fe>::type, SM_cxx_diagrams::fields::Fe, SM_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fe>::type, SM_cxx_diagrams::fields::Fe, SM_cxx_diagrams::fields::VZ>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fe>::type, SM_cxx_diagrams::fields::Fv, typename SM_cxx_diagrams::fields::conj<SM_cxx_diagrams::fields::Hp>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fe>::type, typename SM_cxx_diagrams::fields::conj<SM_cxx_diagrams::fields::Hp>::type, SM_cxx_diagrams::fields::Fv>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fe>::type, typename SM_cxx_diagrams::fields::conj<SM_cxx_diagrams::fields::VWp>::type, SM_cxx_diagrams::fields::Fv>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fu>::type, SM_cxx_diagrams::fields::Fd, SM_cxx_diagrams::fields::Hp>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fu>::type, SM_cxx_diagrams::fields::Fu, SM_cxx_diagrams::fields::Ah>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fu>::type, SM_cxx_diagrams::fields::Fu, SM_cxx_diagrams::fields::hh>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fu>::type, SM_cxx_diagrams::fields::Fu, SM_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fu>::type, SM_cxx_diagrams::fields::Fu, SM_cxx_diagrams::fields::VZ>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fv>::type, SM_cxx_diagrams::fields::Fe, SM_cxx_diagrams::fields::Hp>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fv>::type, SM_cxx_diagrams::fields::Fv, SM_cxx_diagrams::fields::VZ>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fv>::type, SM_cxx_diagrams::fields::Hp, SM_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fv>::type, SM_cxx_diagrams::fields::VWp, SM_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename SM_cxx_diagrams::fields::conj<SM_cxx_diagrams::fields::Hp>::type, SM_cxx_diagrams::fields::Hp, SM_cxx_diagrams::fields::VZ>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};



cxx_diagrams::ChiralVertex unit_charge(const context_base& context);
} // namespace detail

inline double unit_charge(const context_base& context)
{
   return -(detail::unit_charge(context).left().real() /
            fields::Electron::electricCharge);
}

} // namespace SM_cxx_diagrams
} // namespace flexiblesusy

#endif
