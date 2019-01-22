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

// File generated at Tue 22 Jan 2019 18:00:13

/**
 * @file cxx_qft/CMSSMSemiAnalytic_vertices.hpp
 *
 * This file was generated at Tue 22 Jan 2019 18:00:13 with FlexibleSUSY
 * 2.3.0 and SARAH 4.14.1 .
 */

#ifndef CMSSMSemiAnalytic_CXXQFT_VERTICES_H
#define CMSSMSemiAnalytic_CXXQFT_VERTICES_H

#include "concatenate.hpp"
#include "multiindex.hpp"
#include "numerics2.hpp"
#include "wrappers.hpp"

#include <array>

#include <boost/array.hpp>
#include <boost/range/join.hpp>
#include <boost/version.hpp>

#include <boost/mpl/erase.hpp>
#include <boost/mpl/joint_view.hpp>
#include <boost/mpl/vector_c.hpp>

#include <boost/fusion/adapted/boost_array.hpp>
#include <boost/fusion/adapted/mpl.hpp>
#include <boost/fusion/include/vector.hpp>

#if BOOST_VERSION >= 105800
#include <boost/fusion/include/move.hpp>
#else
#include <boost/fusion/include/copy.hpp>
#endif

#define INPUTPARAMETER(p) context.model.get_input().p
#define MODELPARAMETER(p) context.model.get_##p()
#define DERIVEDPARAMETER(p) context.model.p()
#define PHASE(p) context.model.get_##p()

namespace flexiblesusy
{
namespace CMSSMSemiAnalytic_cxx_diagrams
{
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

   class MomentumDifferenceVertex
   {
   private:
      std::complex<double> val;
      int minuendIndex;
      int subtrahendIndex;

   public:
      MomentumDifferenceVertex(std::complex<double> v, int mi, int si)
         : val(v), minuendIndex(mi), subtrahendIndex(si)
      {
      }

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

      bool isZero() const
      {
         return (is_zero(val.real()) && is_zero(val.imag()));
      }
   };

   class InverseMetricVertex
   {
   private:
      std::complex<double> val;

   public:
      InverseMetricVertex(std::complex<double> v) : val(v) {}

      std::complex<double> value() const { return val; }

      bool isZero() const
      {
         return (is_zero(val.real()) && is_zero(val.imag()));
      }
   };

   /**
    * @class VertexData<F...>
    * @brief VertexData data for a vertex with the fields specified by F....
    */
   template <class... Fields>
   struct VertexData;

   struct context_base;

   template <class... Fields>
   class Vertex
   {
      using Data = VertexData<Fields...>;

   public:
      using index_bounds = typename boost::mpl::fold<
         boost::mpl::vector<Fields...>,
         boost::mpl::pair<boost::mpl::vector<>, boost::mpl::vector<>>,
         boost::mpl::pair<
            boost::mpl::joint_view<
               boost::mpl::first<boost::mpl::_1>,
               boost::mpl::first<meta::index_bounds<boost::mpl::_2>>>,
            boost::mpl::joint_view<
               boost::mpl::second<boost::mpl::_1>,
               boost::mpl::second<meta::index_bounds<boost::mpl::_2>>>>>::type;
      using indices_type =
         std::array<int, detail::total_number_of_field_indices<
                            boost::mpl::vector<Fields...>>::value>;
      using vertex_type = typename Data::vertex_type;

      template <int FieldIndex>
      static typename field_indices<typename boost::mpl::at_c<
         boost::mpl::vector<Fields...>, FieldIndex>::type>::type
      field_indices(const indices_type& indices)
      {
         using namespace boost::mpl;
         using fields = vector<Fields...>;

         using result_type = typename field_indices<
            typename boost::mpl::at_c<fields, FieldIndex>::type>::type;

         using preceeding_fields =
            typename erase<fields,
                           typename advance<typename begin<fields>::type,
                                            int_<FieldIndex>>::type,
                           typename end<fields>::type>::type;

         constexpr int offset =
            detail::total_number_of_field_indices<preceeding_fields>::value;
         constexpr int length = std::tuple_size<result_type>::value;

         result_type result_indices;
         std::copy(indices.begin() + offset, indices.begin() + offset + length,
                   result_indices.begin());
         return result_indices;
      }

      static vertex_type evaluate(const indices_type& indices,
                                  const context_base& context);
   };

   template<> struct VertexData<typename bar<fields::Fe>::type, fields::Ah, fields::Fe>
{
   using vertex_type = ChiralVertex;
};

template<> struct VertexData<fields::Fe, fields::Ah, typename bar<fields::Fe>::type>
{
   using vertex_type = ChiralVertex;
};

template<> struct VertexData<fields::VP, typename bar<fields::Fe>::type, fields::Fe>
{
   using vertex_type = ChiralVertex;
};

template<> struct VertexData<typename bar<fields::Fe>::type, fields::Chi, fields::Se>
{
   using vertex_type = ChiralVertex;
};

template<> struct VertexData<fields::Fe, fields::Chi, typename conj<fields::Se>::type>
{
   using vertex_type = ChiralVertex;
};

template<> struct VertexData<fields::VP, typename conj<fields::Se>::type, fields::Se>
{
   using vertex_type = MomentumDifferenceVertex;
};

template<> struct VertexData<typename bar<fields::Fe>::type, fields::Fv, fields::Hpm>
{
   using vertex_type = ChiralVertex;
};

template<> struct VertexData<fields::Fe, typename bar<fields::Fv>::type, typename conj<fields::Hpm>::type>
{
   using vertex_type = ChiralVertex;
};

template<> struct VertexData<fields::VP, typename conj<fields::Hpm>::type, fields::Hpm>
{
   using vertex_type = MomentumDifferenceVertex;
};

template<> struct VertexData<typename bar<fields::Fe>::type, fields::hh, fields::Fe>
{
   using vertex_type = ChiralVertex;
};

template<> struct VertexData<fields::Fe, fields::hh, typename bar<fields::Fe>::type>
{
   using vertex_type = ChiralVertex;
};

template<> struct VertexData<typename bar<fields::Fe>::type, fields::Sv, fields::Cha>
{
   using vertex_type = ChiralVertex;
};

template<> struct VertexData<fields::Fe, typename conj<fields::Sv>::type, typename bar<fields::Cha>::type>
{
   using vertex_type = ChiralVertex;
};

template<> struct VertexData<fields::VP, typename bar<fields::Cha>::type, fields::Cha>
{
   using vertex_type = ChiralVertex;
};

   template<> inline
Vertex<typename bar<fields::Fe>::type, fields::Ah, fields::Fe>::vertex_type
Vertex<typename bar<fields::Fe>::type, fields::Ah, fields::Fe>::evaluate(const indices_type& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,j2)))*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZEL(gt1,j2))*ZA(gt3,0);

   return vertex_type(left, right);
}

template<> inline
Vertex<fields::Fe, fields::Ah, typename bar<fields::Fe>::type>::vertex_type
Vertex<fields::Fe, fields::Ah, typename bar<fields::Fe>::type>::evaluate(const indices_type& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt3 = indices[1];
   const int gt1 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZA = MODELPARAMETER(ZA);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,j2)))*ZA(gt3,0);

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZEL(gt1,j2))*ZA(gt3,0);

   return vertex_type(left, right);
}

template<> inline
Vertex<fields::VP, typename bar<fields::Fe>::type, fields::Fe>::vertex_type
Vertex<fields::VP, typename bar<fields::Fe>::type, fields::Fe>::evaluate(const indices_type& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.5*KroneckerDelta(gt1,gt2)*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   const std::complex<double> right = 0.7745966692414834*g1*Cos(ThetaW)*KroneckerDelta(gt1,gt2);

   return vertex_type(left, right);
}

template<> inline
Vertex<typename bar<fields::Fe>::type, fields::Chi, fields::Se>::vertex_type
Vertex<typename bar<fields::Fe>::type, fields::Chi, fields::Se>::evaluate(const indices_type& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZEL = MODELPARAMETER(ZEL);

   const std::complex<double> left = -1.0954451150103321*g1*Conj(ZN(gt2,0))*SUM(j1,0,2,Conj(ZE(gt3,3 + j1))*Conj(ZER(gt1,j1))) - Conj(ZN(gt2,2))*SUM(j2,0,2,Conj(ZE(gt3,j2))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,j2)));

   const std::complex<double> right = 0.5*(SUM(j1,0,2,Conj(ZE(gt3,j1))*ZEL(gt1,j1))*(1.0954451150103321*g1*ZN(gt2,0) + 1.4142135623730951*g2*ZN(gt2,1)) - 2*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gt3,3 + j1)))*ZEL(gt1,j2))*ZN(gt2,2));

   return vertex_type(left, right);
}

template<> inline
Vertex<fields::Fe, fields::Chi, typename conj<fields::Se>::type>::vertex_type
Vertex<fields::Fe, fields::Chi, typename conj<fields::Se>::type>::evaluate(const indices_type& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const int gt3 = indices[2];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZN = MODELPARAMETER(ZN);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ZER = MODELPARAMETER(ZER);

   const std::complex<double> left = 0.5477225575051661*g1*Conj(ZN(gt1,0))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZE(gt3,j1)) + 0.7071067811865475*g2*Conj(ZN(gt1,1))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZE(gt3,j1)) - Conj(ZN(gt1,2))*SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gt3,3 + j1)));

   const std::complex<double> right = -1.0954451150103321*g1*SUM(j1,0,2,ZE(gt3,3 + j1)*ZER(gt2,j1))*ZN(gt1,0) - SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZE(gt3,j2))*ZN(gt1,2);

   return vertex_type(left, right);
}

template<> inline
Vertex<fields::VP, typename conj<fields::Se>::type, fields::Se>::vertex_type
Vertex<fields::VP, typename conj<fields::Se>::type, fields::Se>::evaluate(const indices_type& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZE = MODELPARAMETER(ZE);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*((0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*SUM(j1,0,2,Conj(ZE(gt1,j1))*ZE(gt2,j1)) + 1.5491933384829668*g1*Cos(ThetaW)*SUM(j1,0,2,Conj(ZE(gt1,3 + j1))*ZE(gt2,3 + j1)));

   return vertex_type(result, minuend_index, subtrahend_index);
}

template<> inline
Vertex<typename bar<fields::Fe>::type, fields::Fv, fields::Hpm>::vertex_type
Vertex<typename bar<fields::Fe>::type, fields::Fv, fields::Hpm>::evaluate(const indices_type& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const int gt3 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> left = SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,gt2))*ZP(gt3,0);

   const std::complex<double> right = 0;

   return vertex_type(left, right);
}

template<> inline
Vertex<fields::Fe, typename bar<fields::Fv>::type, typename conj<fields::Hpm>::type>::vertex_type
Vertex<fields::Fe, typename bar<fields::Fv>::type, typename conj<fields::Hpm>::type>::evaluate(const indices_type& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const int gt3 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZP = MODELPARAMETER(ZP);

   const std::complex<double> left = 0;

   const std::complex<double> right = SUM(j1,0,2,Conj(Ye(j1,gt1))*ZER(gt2,j1))*ZP(gt3,0);

   return vertex_type(left, right);
}

template<> inline
Vertex<fields::VP, typename conj<fields::Hpm>::type, fields::Hpm>::vertex_type
Vertex<fields::VP, typename conj<fields::Hpm>::type, fields::Hpm>::evaluate(const indices_type& indices, const context_base& context)
{
   int minuend_index = 2;
   int subtrahend_index = 1;

   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ZP = MODELPARAMETER(ZP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*(ZP(gt1,0)*ZP(gt2,0) + ZP(gt1,1)*ZP(gt2,1));

   return vertex_type(result, minuend_index, subtrahend_index);
}

template<> inline
Vertex<typename bar<fields::Fe>::type, fields::hh, fields::Fe>::vertex_type
Vertex<typename bar<fields::Fe>::type, fields::hh, fields::Fe>::evaluate(const indices_type& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,j2)))*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZEL(gt1,j2))*ZH(gt3,0);

   return vertex_type(left, right);
}

template<> inline
Vertex<fields::Fe, fields::hh, typename bar<fields::Fe>::type>::vertex_type
Vertex<fields::Fe, fields::hh, typename bar<fields::Fe>::type>::evaluate(const indices_type& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt3 = indices[1];
   const int gt1 = indices[2];
   const auto Ye = MODELPARAMETER(Ye);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZH = MODELPARAMETER(ZH);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(ZEL(gt2,j2))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,j2)))*ZH(gt3,0);

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZEL(gt1,j2))*ZH(gt3,0);

   return vertex_type(left, right);
}

template<> inline
Vertex<typename bar<fields::Fe>::type, fields::Sv, fields::Cha>::vertex_type
Vertex<typename bar<fields::Fe>::type, fields::Sv, fields::Cha>::evaluate(const indices_type& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt3 = indices[1];
   const int gt2 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto UM = MODELPARAMETER(UM);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZER = MODELPARAMETER(ZER);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto UP = MODELPARAMETER(UP);

   const std::complex<double> left = Conj(UM(gt2,1))*SUM(j2,0,2,Conj(ZV(gt3,j2))*SUM(j1,0,2,Conj(ZER(gt1,j1))*Ye(j1,j2)));

   const std::complex<double> right = -(g2*SUM(j1,0,2,Conj(ZV(gt3,j1))*ZEL(gt1,j1))*UP(gt2,0));

   return vertex_type(left, right);
}

template<> inline
Vertex<fields::Fe, typename conj<fields::Sv>::type, typename bar<fields::Cha>::type>::vertex_type
Vertex<fields::Fe, typename conj<fields::Sv>::type, typename bar<fields::Cha>::type>::evaluate(const indices_type& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt3 = indices[1];
   const int gt1 = indices[2];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ye = MODELPARAMETER(Ye);
   const auto UP = MODELPARAMETER(UP);
   const auto ZEL = MODELPARAMETER(ZEL);
   const auto ZV = MODELPARAMETER(ZV);
   const auto ZER = MODELPARAMETER(ZER);
   const auto UM = MODELPARAMETER(UM);

   const std::complex<double> left = -(g2*Conj(UP(gt1,0))*SUM(j1,0,2,Conj(ZEL(gt2,j1))*ZV(gt3,j1)));

   const std::complex<double> right = SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gt2,j1))*ZV(gt3,j2))*UM(gt1,1);

   return vertex_type(left, right);
}

template<> inline
Vertex<fields::VP, typename bar<fields::Cha>::type, fields::Cha>::vertex_type
Vertex<fields::VP, typename bar<fields::Cha>::type, fields::Cha>::evaluate(const indices_type& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto UM = MODELPARAMETER(UM);
   const auto UP = MODELPARAMETER(UP);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.5*(2*g2*Conj(UM(gt2,0))*Sin(ThetaW)*UM(gt1,0) + Conj(UM(gt2,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UM(gt1,1));

   const std::complex<double> right = 0.5*(2*g2*Conj(UP(gt1,0))*Sin(ThetaW)*UP(gt2,0) + Conj(UP(gt1,1))*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW))*UP(gt2,1));

   return vertex_type(left, right);
}

   namespace detail
   {
         static ChiralVertex unit_charge(const context_base& context)
   {
      using vertex_type = ChiralVertex;

      std::array<int, 1> electron_indices = { 0 };
      std::array<int, 0> photon_indices = {};
      std::array<int, 2> indices = concatenate(photon_indices, electron_indices, electron_indices);

         const int gt1 = indices[0];
      const int gt2 = indices[1];
      const auto g1 = MODELPARAMETER(g1);
      const auto g2 = MODELPARAMETER(g2);
      const auto ThetaW = DERIVEDPARAMETER(ThetaW);

      const std::complex<double> left = 0.5*KroneckerDelta(gt1,gt2)*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

      const std::complex<double> right = 0.7745966692414834*g1*Cos(ThetaW)*KroneckerDelta(gt1,gt2);

      return vertex_type(left, right);
   }
   }

   static double unit_charge(const context_base& context)
   {
      return -(detail::unit_charge(context).left().real() /
               fields::Electron::electric_charge);
   }

} // namespace CMSSMSemiAnalytic_cxx_diagrams
} // namespace flexiblesusy

#undef INPUTPARAMETER
#undef MODELPARAMETER
#undef DERIVEDPARAMETER
#undef PHASE

#endif
