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
 * @file cxx_qft/CMSSMSemiAnalytic_vertices.hpp
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#ifndef CMSSMSemiAnalytic_CXXQFT_VERTICES_H
#define CMSSMSemiAnalytic_CXXQFT_VERTICES_H

#include "multiindex.hpp"
#include "numerics2.hpp"

#include "CMSSMSemiAnalytic_fields.hpp"
#include "cxx_qft/vertices.hpp"

#include <array>
#include <algorithm>

#include <boost/mpl/erase.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/joint_view.hpp>
#include <boost/mpl/vector.hpp>

namespace flexiblesusy {
namespace CMSSMSemiAnalytic_cxx_diagrams {

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
template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::Ah, typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Cha>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Cha>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::Ah, typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Fd>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Fd>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::Ah, typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Fe>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::Ah, typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Fu>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Fu>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::Ah, typename CMSSMSemiAnalytic_cxx_diagrams::fields::conj<CMSSMSemiAnalytic_cxx_diagrams::fields::Hpm>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Hpm>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::Ah, typename CMSSMSemiAnalytic_cxx_diagrams::fields::conj<CMSSMSemiAnalytic_cxx_diagrams::fields::Sd>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Sd>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::Ah, typename CMSSMSemiAnalytic_cxx_diagrams::fields::conj<CMSSMSemiAnalytic_cxx_diagrams::fields::Se>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Se>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::Ah, typename CMSSMSemiAnalytic_cxx_diagrams::fields::conj<CMSSMSemiAnalytic_cxx_diagrams::fields::Su>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Su>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::Chi, typename CMSSMSemiAnalytic_cxx_diagrams::fields::conj<CMSSMSemiAnalytic_cxx_diagrams::fields::Se>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::Fe, typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Fe>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::Fe, typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Fe>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::VZ>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::hh, typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Cha>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Cha>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::hh, typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Fd>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Fd>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::hh, typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Fe>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::hh, typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Fu>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Fu>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::hh, typename CMSSMSemiAnalytic_cxx_diagrams::fields::conj<CMSSMSemiAnalytic_cxx_diagrams::fields::Hpm>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Hpm>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::hh, typename CMSSMSemiAnalytic_cxx_diagrams::fields::conj<CMSSMSemiAnalytic_cxx_diagrams::fields::Sd>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Sd>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::hh, typename CMSSMSemiAnalytic_cxx_diagrams::fields::conj<CMSSMSemiAnalytic_cxx_diagrams::fields::Se>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Se>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::hh, typename CMSSMSemiAnalytic_cxx_diagrams::fields::conj<CMSSMSemiAnalytic_cxx_diagrams::fields::Su>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Su>
{
   static cxx_diagrams::ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::hh, typename CMSSMSemiAnalytic_cxx_diagrams::fields::conj<CMSSMSemiAnalytic_cxx_diagrams::fields::VWm>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::VWm>
{
   static cxx_diagrams::InverseMetricVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::VP, CMSSMSemiAnalytic_cxx_diagrams::fields::Cha, typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Cha>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::VP, CMSSMSemiAnalytic_cxx_diagrams::fields::Fd, typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Fd>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::VP, CMSSMSemiAnalytic_cxx_diagrams::fields::Fe, typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Fe>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::VP, CMSSMSemiAnalytic_cxx_diagrams::fields::Fu, typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Fu>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::VP, CMSSMSemiAnalytic_cxx_diagrams::fields::Hpm, typename CMSSMSemiAnalytic_cxx_diagrams::fields::conj<CMSSMSemiAnalytic_cxx_diagrams::fields::Hpm>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::VP, CMSSMSemiAnalytic_cxx_diagrams::fields::Sd, typename CMSSMSemiAnalytic_cxx_diagrams::fields::conj<CMSSMSemiAnalytic_cxx_diagrams::fields::Sd>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::VP, CMSSMSemiAnalytic_cxx_diagrams::fields::Se, typename CMSSMSemiAnalytic_cxx_diagrams::fields::conj<CMSSMSemiAnalytic_cxx_diagrams::fields::Se>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::VP, CMSSMSemiAnalytic_cxx_diagrams::fields::Su, typename CMSSMSemiAnalytic_cxx_diagrams::fields::conj<CMSSMSemiAnalytic_cxx_diagrams::fields::Su>::type>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::VP, CMSSMSemiAnalytic_cxx_diagrams::fields::VWm, typename CMSSMSemiAnalytic_cxx_diagrams::fields::conj<CMSSMSemiAnalytic_cxx_diagrams::fields::VWm>::type>
{
   static cxx_diagrams::TripleVectorVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::VP, typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Fe>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::VZ, CMSSMSemiAnalytic_cxx_diagrams::fields::Cha, typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Cha>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::VZ, CMSSMSemiAnalytic_cxx_diagrams::fields::Fd, typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Fd>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::VZ, CMSSMSemiAnalytic_cxx_diagrams::fields::Fe, typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Fe>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::VZ, CMSSMSemiAnalytic_cxx_diagrams::fields::Fu, typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Fu>::type>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<CMSSMSemiAnalytic_cxx_diagrams::fields::VZ, typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Fe>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Cha>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Cha, CMSSMSemiAnalytic_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Fe>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Cha, CMSSMSemiAnalytic_cxx_diagrams::fields::Sv>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Fe>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Fe, CMSSMSemiAnalytic_cxx_diagrams::fields::Ah>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Fe>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Fe, CMSSMSemiAnalytic_cxx_diagrams::fields::hh>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Fe>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Fe, CMSSMSemiAnalytic_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Fe>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Fe, CMSSMSemiAnalytic_cxx_diagrams::fields::VZ>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Fe>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Hpm, CMSSMSemiAnalytic_cxx_diagrams::fields::Fv>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Fe>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Se, CMSSMSemiAnalytic_cxx_diagrams::fields::Chi>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Fe>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::VWm, CMSSMSemiAnalytic_cxx_diagrams::fields::Fv>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Fv>::type, typename CMSSMSemiAnalytic_cxx_diagrams::fields::conj<CMSSMSemiAnalytic_cxx_diagrams::fields::Hpm>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Fv>::type, typename CMSSMSemiAnalytic_cxx_diagrams::fields::conj<CMSSMSemiAnalytic_cxx_diagrams::fields::VWm>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CMSSMSemiAnalytic_cxx_diagrams::fields::conj<CMSSMSemiAnalytic_cxx_diagrams::fields::Hpm>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Hpm, CMSSMSemiAnalytic_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CMSSMSemiAnalytic_cxx_diagrams::fields::conj<CMSSMSemiAnalytic_cxx_diagrams::fields::Hpm>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::VWm, CMSSMSemiAnalytic_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::InverseMetricVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CMSSMSemiAnalytic_cxx_diagrams::fields::conj<CMSSMSemiAnalytic_cxx_diagrams::fields::Se>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Se, CMSSMSemiAnalytic_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CMSSMSemiAnalytic_cxx_diagrams::fields::conj<CMSSMSemiAnalytic_cxx_diagrams::fields::Sv>::type, typename CMSSMSemiAnalytic_cxx_diagrams::fields::bar<CMSSMSemiAnalytic_cxx_diagrams::fields::Cha>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Fe>
{
   static cxx_diagrams::ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CMSSMSemiAnalytic_cxx_diagrams::fields::conj<CMSSMSemiAnalytic_cxx_diagrams::fields::VWm>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::Hpm, CMSSMSemiAnalytic_cxx_diagrams::fields::VP>
{
   static cxx_diagrams::InverseMetricVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<typename CMSSMSemiAnalytic_cxx_diagrams::fields::conj<CMSSMSemiAnalytic_cxx_diagrams::fields::VWm>::type, CMSSMSemiAnalytic_cxx_diagrams::fields::VWm, CMSSMSemiAnalytic_cxx_diagrams::fields::VP>
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

} // namespace CMSSMSemiAnalytic_cxx_diagrams
} // namespace flexiblesusy

#endif
