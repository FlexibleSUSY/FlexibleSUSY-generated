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
 * @file cxx_qft/MSSMEFTHiggs_vertices.hpp
 *
 * This file was generated with FlexibleSUSY 2.6.0 and SARAH 4.14.5 .
 */

#ifndef MSSMEFTHiggs_CXXQFT_VERTICES_H
#define MSSMEFTHiggs_CXXQFT_VERTICES_H

#include "multiindex.hpp"
#include "numerics2.hpp"

#include "MSSMEFTHiggs_fields.hpp"

#include <array>
#include <algorithm>

#include <boost/mpl/erase.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/joint_view.hpp>
#include <boost/mpl/vector.hpp>

namespace flexiblesusy {
namespace MSSMEFTHiggs_cxx_diagrams {

   class ScalarVertex
   {
   private:
      std::complex<double> val;

   public:
      ScalarVertex(std::complex<double> v) : val(v) {}

      std::complex<double> value() const { return val; }

      bool isZero() const
      {
         return (is_zero(val.real()) && is_zero(val.imag()));
      }
   };

   class ChiralVertex
   {
   private:
      std::pair<std::complex<double>, std::complex<double>> value;

   public:
      ChiralVertex(const std::complex<double>& left,
                   const std::complex<double>& right)
         : value(left, right)
      {
      }

      std::complex<double> left() const { return value.first; }
      std::complex<double> right() const { return value.second; }

      bool isZero() const
      {
         return (is_zero(value.first.real()) && is_zero(value.first.imag()) &&
                 is_zero(value.second.real()) && is_zero(value.second.imag()));
      }
   };

/** \brief A class representing a numerically evaluated
 * tree-level vertex that is proportional to a momentum.
 * It consists of a complex number as well as an index
 * corresponding to the index of the field to whose
 * momentum the vertex is proportional.
 **/
class MomentumVertex {
  std::complex<double> val;
   int ind;
public:
   /** \brief Contruct a MomentumVertex from a
    * complex number representing and a field index.
    **/
   MomentumVertex(const std::complex<double>& v, int i)
      : val(v), ind(i)
   {}

   /** \brief Retrieve the index of the field to whose
    * momentum the vertex is proportional.
    * \returns the appropriate index
    **/
   int index() const { return ind; }

   /** \brief Retrieve the numerical value of the vertex
    * \param i The index of the field to whose momentum
    * the vertex is proportional.
    * \returns the coefficient of the even permutation
    **/
   std::complex<double> value(int i) const
   {
      if (i != ind)
         throw std::invalid_argument(
            "MomentumVertex: Wrong index specified");

      return val;
   }

   bool isZero() const
   {
      return (is_zero(val.real()) && is_zero(val.imag()));
   }
};

/** \brief A class representing a numerically evaluated
 * tree-level vertex with three vector bosons.
 * It consists of one complex number as well as an \a ordering
 * encoding whether the complex number is taken to be the
 * coefficient of
 *
 * \f{equation}{
 * g[l1, l2] * (p[field1, l3] - p[field2, lIndex3]) +
 * g[l2, l3] * (p[field2, l1] - p[field3, lIndex1]) +
 * g[l1, l3] * (p[field3, l2] - p[field1, lIndex2])
 * \f}
 *
 * or its negative.
 * The former corresponds to the \a even permutation and
 * the latter to the \a odd permutation.
 **/
class TripleVectorVertex {
public:
   struct even_permutation {};
   struct odd_permutation {};
private:
   std::complex<double> val;
   bool even;
public:
   /** \brief Contruct a TripleVectorVertex from a
    * complex number representing the even coefficient.
    **/
   TripleVectorVertex(const std::complex<double>& v,
                      even_permutation)
      : val(v), even(true)
   {}

   /** \brief Contruct a TripleVectorVertex from a
    * complex number representing the odd coefficient.
    **/
   TripleVectorVertex(const std::complex<double>& v,
                      odd_permutation)
      : val(v), even(false)
   {}

   /** \brief Check whether the value in the vertex is stored
    * as proportional to the even permutation.
    * \returns true if yes and false otherwise
    **/
   bool is_even() const { return even; }

   /** \brief Retrieve the coefficient of the even permutation
    * \returns the coefficient of the even permutation
    **/
   std::complex<double> value(even_permutation) const
   { return even ? val : - val; }

   /** \brief Retrieve the coefficient of the odd permutation
    * \returns the coefficient of the odd permutation
    **/
   std::complex<double> value(odd_permutation) const
   { return even ? - val : val; }

   bool isZero() const
   {
      return (is_zero(val.real()) && is_zero(val.imag()));
   }
};

/** \brief A class representing a numerically evaluated
 * tree-level vertex with four vector bosons.
 * It consists of three complex numbers corresponding to
 * (in order) the basis expansion with respect to the basis:
 *
 * \f{equation}{
 * ( g[l1, l2] g[l3, l4], g[l1, l3] g[l2, l4], g[l1, l4] g[l2, l3] )
 * \f}
 **/
class QuadrupleVectorVertex {
   std::complex<double> part1, part2, part3;

public:
   /** \brief Contruct a QuadrupleVectorVertex from three
    * complex numbers representing the coefficients in the
    * basis expansion.
    **/
   QuadrupleVectorVertex(const std::complex<double>& p1,
                         const std::complex<double>& p2,
                         const std::complex<double>& p3)
      : part1(p1), part2(p2), part3(p3)
   {}

   /** \brief Retrieve the coefficient of \f$ g[l1, l2] g[l3, l4] \f$
    * \returns the corresponding coefficient
    **/
   std::complex<double> value1() const { return part1; }

   /** \brief Retrieve the coefficient of \f$ g[l1, l3] g[l2, l4] \f$
    * \returns the corresponding coefficient
    **/
   std::complex<double> value2() const { return part2; }

   /** \brief Retrieve the coefficient of \f$ g[l1, l4] g[l2, l3] \f$
    * \returns the corresponding coefficient
    **/
   std::complex<double> value3() const { return part3; }

   bool isZero() const
   {
      return (is_zero(part1.real()) && is_zero(part1.imag()) &&
              is_zero(part2.real()) && is_zero(part2.imag()) &&
              is_zero(part3.real()) && is_zero(part3.imag()));
   }
};

class MomentumDifferenceVertex {
   std::complex<double> val;
   int minuendIndex;
   int subtrahendIndex;
public:
   MomentumDifferenceVertex(std::complex<double> v, int mi, int si)
      : val(v), minuendIndex(mi), subtrahendIndex(si) {}

   std::complex<double> value(int mi, int si) const
   {
      if (mi == minuendIndex && si == subtrahendIndex)
         return val;
      if (mi == subtrahendIndex && si == minuendIndex)
         return -val;

      throw std::invalid_argument(
         "MomentumDifferenceVertex: Wrong index combination");
      return 0.0;
   }

   int incoming_index() const { return minuendIndex; }
   int outgoing_index() const { return subtrahendIndex; }

   bool isZero() const
   {
      return (is_zero(val.real()) && is_zero(val.imag()));
   }
};

class InverseMetricVertex {
   std::complex<double> val;
public:
   InverseMetricVertex(std::complex<double> v) : val(v) {}

   std::complex<double> value() const { return val; }

   bool isZero() const
   {
      return (is_zero(val.real()) && is_zero(val.imag()));
   }
};

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
      detail::total_number_of_field_indices<
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
         detail::total_number_of_field_indices<preceeding_fields>::value;
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
template<> struct VertexImpl<fields::Ah, fields::Ah, fields::Ah, fields::Ah>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Ah, fields::Ah, fields::hh, fields::hh>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Ah, fields::Ah, fields::hh>
{
   static ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Ah, fields::Ah, fields::Hpm, typename fields::conj<fields::Hpm>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Ah, fields::Ah, fields::Sd, typename fields::conj<fields::Sd>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Ah, fields::Ah, fields::Se, typename fields::conj<fields::Se>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Ah, fields::Ah, fields::Su, typename fields::conj<fields::Su>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Ah, fields::Ah, fields::Sv, typename fields::conj<fields::Sv>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Ah, fields::Ah, fields::VZ, fields::VZ>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Ah, fields::Ah, typename fields::conj<fields::VWm>::type, fields::VWm>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Ah, fields::hh, fields::Hpm, typename fields::conj<fields::Hpm>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Ah, fields::hh, fields::VZ>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Ah, fields::Hpm, fields::Su, typename fields::conj<fields::Sd>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Ah, fields::Hpm, fields::Sv, typename fields::conj<fields::Se>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Ah, fields::Hpm, typename fields::conj<fields::Hpm>::type>
{
   static ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Ah, fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VP>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Ah, fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VZ>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Ah, fields::Hpm, typename fields::conj<fields::VWm>::type>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Ah, fields::Sd, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Su>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Ah, fields::Sd, typename fields::conj<fields::Sd>::type>
{
   static ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Ah, fields::Se, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Sv>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Ah, fields::Se, typename fields::conj<fields::Se>::type>
{
   static ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Ah, fields::Su, typename fields::conj<fields::Su>::type>
{
   static ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Ah, typename fields::bar<fields::Fe>::type, fields::Fe>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Ah, typename fields::conj<fields::Hpm>::type, fields::VP, fields::VWm>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Ah, typename fields::conj<fields::Hpm>::type, fields::VWm, fields::VZ>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Ah, typename fields::conj<fields::Hpm>::type, fields::VWm>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Cha, fields::Fu, typename fields::conj<fields::Sd>::type>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Cha, fields::Fv, typename fields::conj<fields::Se>::type>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Chi, fields::Cha, typename fields::conj<fields::Hpm>::type>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Chi, fields::Cha, typename fields::conj<fields::VWm>::type>
{
   static ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Chi, fields::Chi, fields::Ah>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Chi, fields::Chi, fields::hh>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Chi, fields::Chi, fields::VZ>
{
   static ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Chi, fields::Fd, typename fields::conj<fields::Sd>::type>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Chi, fields::Fe, typename fields::conj<fields::Se>::type>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Chi, fields::Fu, typename fields::conj<fields::Su>::type>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Chi, fields::Fv, typename fields::conj<fields::Sv>::type>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Chi, typename fields::conj<fields::Se>::type, fields::Fe>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Glu, fields::Fd, typename fields::conj<fields::Sd>::type>
{
   static ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Glu, fields::Fu, typename fields::conj<fields::Su>::type>
{
   static ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Glu, fields::Glu, fields::VG>
{
   static ChiralVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::hh, fields::hh, fields::hh, fields::hh>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::hh, fields::hh, fields::hh>
{
   static ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::hh, fields::hh, fields::Hpm, typename fields::conj<fields::Hpm>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::hh, fields::hh, fields::Sd, typename fields::conj<fields::Sd>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::hh, fields::hh, fields::Se, typename fields::conj<fields::Se>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::hh, fields::hh, fields::Su, typename fields::conj<fields::Su>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::hh, fields::hh, fields::Sv, typename fields::conj<fields::Sv>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::hh, fields::hh, fields::VZ, fields::VZ>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::hh, fields::hh, typename fields::conj<fields::VWm>::type, fields::VWm>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::hh, fields::Hpm, fields::Su, typename fields::conj<fields::Sd>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::hh, fields::Hpm, fields::Sv, typename fields::conj<fields::Se>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::hh, fields::Hpm, typename fields::conj<fields::Hpm>::type>
{
   static ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::hh, fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VP>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::hh, fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VZ>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::hh, fields::Hpm, typename fields::conj<fields::VWm>::type>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::hh, fields::Sd, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Su>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::hh, fields::Sd, typename fields::conj<fields::Sd>::type>
{
   static ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::hh, fields::Se, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Sv>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::hh, fields::Se, typename fields::conj<fields::Se>::type>
{
   static ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::hh, fields::Su, typename fields::conj<fields::Su>::type>
{
   static ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::hh, fields::Sv, typename fields::conj<fields::Sv>::type>
{
   static ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::hh, fields::VZ, fields::VZ>
{
   static InverseMetricVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::hh, typename fields::bar<fields::Fe>::type, fields::Fe>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::hh, typename fields::conj<fields::Hpm>::type, fields::VP, fields::VWm>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::hh, typename fields::conj<fields::Hpm>::type, fields::VWm, fields::VZ>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::hh, typename fields::conj<fields::Hpm>::type, fields::VWm>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::hh, typename fields::conj<fields::VWm>::type, fields::VWm>
{
   static InverseMetricVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Hpm, fields::Hpm, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Hpm>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Hpm, fields::Sd, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Sd>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Hpm, fields::Se, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Se>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Hpm, fields::Su, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Su>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Hpm, fields::Su, typename fields::conj<fields::Sd>::type>
{
   static ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Hpm, fields::Sv, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Sv>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Hpm, fields::Sv, typename fields::conj<fields::Se>::type>
{
   static ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, fields::VP, fields::VP>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, fields::VP, fields::VZ>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, fields::VP>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, fields::VZ, fields::VZ>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, fields::VZ>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Hpm, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::VWm>::type, fields::VWm>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VP>
{
   static InverseMetricVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Hpm, typename fields::conj<fields::VWm>::type, fields::VZ>
{
   static InverseMetricVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Sd, fields::Sd, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::Sd>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Sd, fields::Se, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::Se>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Sd, fields::Su, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::Su>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Sd, fields::Sv, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::Sv>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Sd, fields::Sv, typename fields::conj<fields::Se>::type, typename fields::conj<fields::Su>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Sd, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Su>::type>
{
   static ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VG, fields::VG>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VG, fields::VP>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VG, fields::VZ>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VG>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VP, fields::VP>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VP, fields::VZ>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VP>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VZ, fields::VZ>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, fields::VZ>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Sd, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::VWm>::type, fields::VWm>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Sd, typename fields::conj<fields::Su>::type, typename fields::conj<fields::VWm>::type, fields::VG>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Sd, typename fields::conj<fields::Su>::type, typename fields::conj<fields::VWm>::type, fields::VP>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Sd, typename fields::conj<fields::Su>::type, typename fields::conj<fields::VWm>::type, fields::VZ>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Sd, typename fields::conj<fields::Su>::type, typename fields::conj<fields::VWm>::type>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Se, fields::Se, typename fields::conj<fields::Se>::type, typename fields::conj<fields::Se>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Se, fields::Su, typename fields::conj<fields::Sd>::type, typename fields::conj<fields::Sv>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Se, fields::Su, typename fields::conj<fields::Se>::type, typename fields::conj<fields::Su>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Se, fields::Sv, typename fields::conj<fields::Se>::type, typename fields::conj<fields::Sv>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Se, typename fields::conj<fields::Hpm>::type, typename fields::conj<fields::Sv>::type>
{
   static ScalarVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Se, typename fields::conj<fields::Se>::type, fields::VP, fields::VP>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Se, typename fields::conj<fields::Se>::type, fields::VP, fields::VZ>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Se, typename fields::conj<fields::Se>::type, fields::VP>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Se, typename fields::conj<fields::Se>::type, fields::VZ, fields::VZ>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Se, typename fields::conj<fields::Se>::type, fields::VZ>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Se, typename fields::conj<fields::Se>::type, typename fields::conj<fields::VWm>::type, fields::VWm>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Se, typename fields::conj<fields::Sv>::type, typename fields::conj<fields::VWm>::type, fields::VP>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Se, typename fields::conj<fields::Sv>::type, typename fields::conj<fields::VWm>::type, fields::VZ>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Se, typename fields::conj<fields::Sv>::type, typename fields::conj<fields::VWm>::type>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Su, fields::Su, typename fields::conj<fields::Su>::type, typename fields::conj<fields::Su>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Su, fields::Sv, typename fields::conj<fields::Su>::type, typename fields::conj<fields::Sv>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Su, typename fields::conj<fields::Sd>::type, fields::VG, fields::VWm>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Su, typename fields::conj<fields::Sd>::type, fields::VP, fields::VWm>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Su, typename fields::conj<fields::Sd>::type, fields::VWm, fields::VZ>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Su, typename fields::conj<fields::Sd>::type, fields::VWm>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VG, fields::VG>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VG, fields::VP>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VG, fields::VZ>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VG>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VP, fields::VP>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VP, fields::VZ>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VP>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VZ, fields::VZ>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, fields::VZ>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Su, typename fields::conj<fields::Su>::type, typename fields::conj<fields::VWm>::type, fields::VWm>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Sv, fields::Sv, typename fields::conj<fields::Sv>::type, typename fields::conj<fields::Sv>::type>
{
   static ScalarVertex evaluate(const std::array<int, 4>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Sv, typename fields::conj<fields::Se>::type, fields::VP, fields::VWm>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Sv, typename fields::conj<fields::Se>::type, fields::VWm, fields::VZ>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Sv, typename fields::conj<fields::Se>::type, fields::VWm>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Sv, typename fields::conj<fields::Sv>::type, fields::VZ, fields::VZ>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Sv, typename fields::conj<fields::Sv>::type, fields::VZ>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::Sv, typename fields::conj<fields::Sv>::type, typename fields::conj<fields::VWm>::type, fields::VWm>
{
   static InverseMetricVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VG, fields::VG, fields::VG>
{
   static TripleVectorVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VWm, fields::Ah, typename fields::conj<fields::Hpm>::type>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VWm, fields::Chi, typename fields::bar<fields::Cha>::type>
{
   static ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VWm, fields::hh, typename fields::conj<fields::Hpm>::type>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VWm, fields::hh, typename fields::conj<fields::VWm>::type>
{
   static InverseMetricVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VWm, fields::Su, typename fields::conj<fields::Sd>::type>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VWm, fields::Sv, typename fields::conj<fields::Se>::type>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VWm, typename fields::bar<fields::Cha>::type, fields::Chi>
{
   static ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VWm, typename fields::conj<fields::Hpm>::type, fields::Ah>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VWm, typename fields::conj<fields::Hpm>::type, fields::hh>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VWm, typename fields::conj<fields::Hpm>::type, fields::VP>
{
   static InverseMetricVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VWm, typename fields::conj<fields::Hpm>::type, fields::VZ>
{
   static InverseMetricVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VWm, typename fields::conj<fields::Sd>::type, fields::Su>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VWm, typename fields::conj<fields::Se>::type, fields::Sv>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VZ, fields::Ah, fields::hh>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VZ, fields::Cha, typename fields::bar<fields::Cha>::type>
{
   static ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VZ, fields::Chi, fields::Chi>
{
   static ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VZ, fields::hh, fields::Ah>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VZ, fields::hh, fields::VZ>
{
   static InverseMetricVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VZ, fields::Hpm, typename fields::conj<fields::Hpm>::type>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VZ, fields::Hpm, typename fields::conj<fields::VWm>::type>
{
   static InverseMetricVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VZ, fields::Sd, typename fields::conj<fields::Sd>::type>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VZ, fields::Se, typename fields::conj<fields::Se>::type>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VZ, fields::Su, typename fields::conj<fields::Su>::type>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VZ, fields::Sv, typename fields::conj<fields::Sv>::type>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VZ, typename fields::bar<fields::Cha>::type, fields::Cha>
{
   static ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VZ, typename fields::conj<fields::Hpm>::type, fields::Hpm>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VZ, typename fields::conj<fields::Hpm>::type, fields::VWm>
{
   static InverseMetricVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VZ, typename fields::conj<fields::Sd>::type, fields::Sd>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VZ, typename fields::conj<fields::Se>::type, fields::Se>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VZ, typename fields::conj<fields::Su>::type, fields::Su>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<fields::VZ, typename fields::conj<fields::Sv>::type, fields::Sv>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Cha>::type, fields::Cha, fields::Ah>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Cha>::type, fields::Cha, fields::hh>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Cha>::type, fields::Cha, fields::VP>
{
   static ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Cha>::type, fields::Cha, fields::VZ>
{
   static ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Cha>::type, fields::Chi, fields::Hpm>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Cha>::type, fields::Chi, fields::VWm>
{
   static ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Cha>::type, fields::Fd, typename fields::conj<fields::Su>::type>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Cha>::type, fields::Fe, typename fields::conj<fields::Sv>::type>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Cha>::type, typename fields::bar<fields::Fu>::type, fields::Sd>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Cha>::type, typename fields::bar<fields::Fv>::type, fields::Se>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fd>::type, fields::Cha, fields::Su>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fd>::type, fields::Chi, fields::Sd>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fd, fields::Ah>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fd, fields::hh>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fd, fields::VG>
{
   static ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fd, fields::VP>
{
   static ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fd, fields::VZ>
{
   static ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fu, fields::Hpm>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fu, fields::VWm>
{
   static ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fd>::type, fields::Glu, fields::Sd>
{
   static ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fe>::type, fields::Cha, fields::Sv>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fe>::type, fields::Chi, fields::Se>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fe, fields::Ah>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fe, fields::hh>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fe, fields::VP>
{
   static ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fe, fields::VZ>
{
   static ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fv, fields::Hpm>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fv, fields::VWm>
{
   static ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fe>::type, fields::Hpm, fields::Fv>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fe>::type, fields::Se, fields::Chi>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fu>::type, fields::Chi, fields::Su>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fd, typename fields::conj<fields::Hpm>::type>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fd, typename fields::conj<fields::VWm>::type>
{
   static ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fu, fields::Ah>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fu, fields::hh>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fu, fields::VG>
{
   static ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fu, fields::VP>
{
   static ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fu, fields::VZ>
{
   static ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fu>::type, fields::Glu, fields::Su>
{
   static ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fv>::type, fields::Chi, fields::Sv>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fv>::type, fields::Fe, typename fields::conj<fields::Hpm>::type>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fv>::type, fields::Fe, typename fields::conj<fields::VWm>::type>
{
   static ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fv>::type, fields::Fv, fields::VZ>
{
   static ChiralVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::Fv>::type, typename fields::conj<fields::Hpm>::type, fields::Fe>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::gG>::type, fields::gG, fields::VG>
{
   static MomentumVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::gP>::type, fields::gWmC, fields::VWm>
{
   static MomentumVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::gP>::type, fields::gWm, typename fields::conj<fields::VWm>::type>
{
   static MomentumVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gP, typename fields::conj<fields::Hpm>::type>
{
   static ScalarVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gP, typename fields::conj<fields::VWm>::type>
{
   static MomentumVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gWmC, fields::Ah>
{
   static ScalarVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gWmC, fields::hh>
{
   static ScalarVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gWmC, fields::VP>
{
   static MomentumVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gWmC, fields::VZ>
{
   static MomentumVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gZ, typename fields::conj<fields::Hpm>::type>
{
   static ScalarVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::gWmC>::type, fields::gZ, typename fields::conj<fields::VWm>::type>
{
   static MomentumVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::gWm>::type, fields::gP, fields::Hpm>
{
   static ScalarVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::gWm>::type, fields::gP, fields::VWm>
{
   static MomentumVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::gWm>::type, fields::gWm, fields::Ah>
{
   static ScalarVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::gWm>::type, fields::gWm, fields::hh>
{
   static ScalarVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::gWm>::type, fields::gWm, fields::VP>
{
   static MomentumVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::gWm>::type, fields::gWm, fields::VZ>
{
   static MomentumVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::gWm>::type, fields::gZ, fields::Hpm>
{
   static ScalarVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::gWm>::type, fields::gZ, fields::VWm>
{
   static MomentumVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::gZ>::type, fields::gWmC, fields::Hpm>
{
   static ScalarVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::gZ>::type, fields::gWmC, fields::VWm>
{
   static MomentumVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::gZ>::type, fields::gWm, typename fields::conj<fields::Hpm>::type>
{
   static ScalarVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::gZ>::type, fields::gWm, typename fields::conj<fields::VWm>::type>
{
   static MomentumVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::bar<fields::gZ>::type, fields::gZ, fields::hh>
{
   static ScalarVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::conj<fields::Hpm>::type, fields::Hpm, fields::VP>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::conj<fields::Hpm>::type, fields::VP, fields::VWm>
{
   static InverseMetricVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::conj<fields::Hpm>::type, fields::VWm, fields::VZ>
{
   static InverseMetricVertex evaluate(const std::array<int, 1>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::conj<fields::Se>::type, fields::Se, fields::VP>
{
   static MomentumDifferenceVertex evaluate(const std::array<int, 2>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::conj<fields::Sv>::type, typename fields::bar<fields::Cha>::type, fields::Fe>
{
   static ChiralVertex evaluate(const std::array<int, 3>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::conj<fields::VWm>::type, fields::VP, fields::VP, fields::VWm>
{
   static QuadrupleVectorVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::conj<fields::VWm>::type, fields::VP, fields::VWm, fields::VZ>
{
   static QuadrupleVectorVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::conj<fields::VWm>::type, fields::VP, fields::VWm>
{
   static TripleVectorVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::conj<fields::VWm>::type, fields::VWm, fields::VZ, fields::VZ>
{
   static QuadrupleVectorVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::conj<fields::VWm>::type, fields::VWm, fields::VZ>
{
   static TripleVectorVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};

template<> struct VertexImpl<typename fields::conj<fields::VWm>::type, typename fields::conj<fields::VWm>::type, fields::VWm, fields::VWm>
{
   static QuadrupleVectorVertex evaluate(const std::array<int, 0>& indices, const context_base& context);
};



ChiralVertex unit_charge(const context_base& context);
} // namespace detail

inline double unit_charge(const context_base& context)
{
   return -(detail::unit_charge(context).left().real() /
            fields::Electron::electric_charge);
}

} // namespace MSSMEFTHiggs_cxx_diagrams
} // namespace flexiblesusy

#endif
