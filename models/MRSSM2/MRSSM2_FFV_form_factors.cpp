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
 * @file MRSSM2_FFV_form_factors.cpp
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#include "MRSSM2_FFV_form_factors.hpp"
#include "MRSSM2_mass_eigenstates.hpp"
#include "concatenate.hpp"
#include "cxx_qft/MRSSM2_qft.hpp"
#include "ffv_loop_functions.hpp"
#include "wrappers.hpp"

#include <complex>
#include <valarray>

using namespace flexiblesusy;
using namespace flexiblesusy::ffv_loop_functions::one_loop;
using namespace MRSSM2_cxx_diagrams;
using namespace MRSSM2_cxx_diagrams::fields;

namespace {

static constexpr double oneOver32PiSqr = 0.5*oneOver16PiSqr;

/**
 * @class FFV_SSF
 * @brief A template that calculate contributions to the FFV form
 *        factors of a given particles in a one loop diagram
 *        specified by a vector emitters and an exchange particle.
 * @tparam Args Specifies in order the field of which to
 *              calculate the electric dipole moment,
 *              the photon emitter and the exchange particle
 *              in a one-loop diagram where the photon emitter
 *              is a scalar and the exchange particle a fermion.
 *
 * This template evaluates the contribution to the electric
 * dipole moment of a one loop diagram with fields given by
 * \a Args.
 */
template <class Fj, class Fi, class V, class S1, class S2, class F>
struct FFV_SSF {
   static std::valarray<std::complex<double>>
   value(const typename field_indices<Fj>::type& indices_in,
         const typename field_indices<Fi>::type& indices_out,
         context_base const& context,
         bool discard_SM_contributions);
};

/**
* @class FFV_FFS
* @brief A template that calculate contributions to the FFV form
*        factors of a given particle in a one loop diagram
*        specified by a vector emitters and an exchange particle.
* @tparam Args Specifies in order the field of which to
*              calculate the electric dipole moment,
*              the photon emitter and the exchange particle
*              in a one-loop diagram where the photon emitter
*              is a fermion and the exchange particle a scalar.
*
* This template evaluates the contribution to the electric
* dipole moment of a one loop diagram with fields given by
* \a Args.
*/
template <class Fj, class Fi, class V, class F1, class F2, class S>
struct FFV_FFS {
   static std::valarray<std::complex<double>>
   value(const typename field_indices<Fj>::type& indices_in,
         const typename field_indices<Fi>::type& indices_out,
         context_base const& context,
         bool discard_SM_contributions);
};

/**
 * @class FFV_VVF
 * @brief A template that calculates contributions to the FFV form
 *          factors of a given particle in a one loop diagram
 *          specified by a vector emitter and an exchange particle.
 * @tparam Args Specifies in order the field of which to
 *              calculate the electric dipole moment,
 *              the photon emitter and the exchange particle
 *              in a one-loop diagram where the photon emitter
 *              is a vector and the exchange particle a fermion.
 *
 * This template evaluates the contribution to the electric
 * dipole moment of a one loop diagram with fields given by
 * \a Args.
 */
template <class Fj, class Fi, class P, class V1, class V2, class F>
struct FFV_VVF {
    static std::valarray<std::complex<double>>
    value(const typename field_indices<Fj>::type& indices_in,
          const typename field_indices<Fi>::type& indices_out,
          context_base const& context,
          bool discard_SM_contributions);
};

/**
 * @class FFV_FFV
 * @brief A template that calculates contributions to the FFV form
 *        factors of a given particle in a one loop diagram
 *        specified by a vector emitter and an exchange particle.
 * @tparam Args Specifies in order the field of which to
 *              calculate the electric dipole moment,
 *              the photon emitter and the exchange particle
 *              in a one-loop diagram where the photon emitter
 *              is a fermion and the exchange particle a vector.
 *
 * This template evaluates the contribution to the electric
 * dipole moment of a one loop diagram with fields given by
 * \a Args.
 */
template <class Fj, class Fi, class P, class F1, class F2, class V>
struct FFV_FFV {
    static std::valarray<std::complex<double>>
    value(const typename field_indices<Fj>::type& indices_in,
          const typename field_indices<Fi>::type& indices_out,
          context_base const& context,
          bool discard_SM_contributions);
};

/**
 * @class FFV_VSF
 * @brief A template that calculates contributions to the FFV form
 *        factors of a given particle in a one loop diagram
 *        specified by a vector emitter and an exchange particle.
 * @tparam Args Specifies in order the field of which to
 *         calculate the electric dipole moment,
 *         the photon emitter and the exchange particle
 *         in a one-loop diagram where the photon is emitted
 *         from a vector-goldstone vertex and the exchange
 *         particle is a fermion.
 *
 * This template evaluates the contribution to the electric
 * dipole moment of a one loop diagram with fields given by
 * \a Args.
 */
template <class Fj, class Fi, class P, class V1, class G2, class F>
struct FFV_VSF {
    static std::valarray<std::complex<double>>
    value(const typename field_indices<Fj>::type& indices_in,
          const typename field_indices<Fi>::type& indices_out,
          context_base const& context,
          bool discard_SM_contributions);
};

/**
 * @class FFV_SVF
 * @brief A template that calculates contributions to the FFV form
 *        factors of a given particle in a one loop diagram
 *        specified by a vector emitter and an exchange particle.
 * @tparam Args Specifies in order the field of which to
 *              calculate the electric dipole moment,
 *              the photon emitter and the exchange particle
 *              in a one-loop diagram where the photon is emitted
 *              from a goldstone-vector vertex and the exchange
 *              particle is a fermion.
 *
 * This template evaluates the contribution to the electric
 * dipole moment of a one loop diagram with fields given by
 * \a Args.
 */
template <class Fj, class Fi, class P, class G1, class V2, class F>
struct FFV_SVF {
    static std::valarray<std::complex<double>>
    value(const typename field_indices<Fj>::type& indices_in,
          const typename field_indices<Fi>::type& indices_out,
          context_base const& context,
          bool discard_SM_contributions);
};

} // anonymous namespace

namespace flexiblesusy {
namespace MRSSM2_FFV_form_factors {

std::valarray<std::complex<double>> calculate_Fe_Fe_VP_form_factors (
   int generationIndex1, int generationIndex2,
   const MRSSM2_mass_eigenstates& model, bool discard_SM_contributions) {

   context_base context {model};
   std::array<int, 1> indices1 = {generationIndex1 };
   std::array<int, 1> indices2 = {generationIndex2 };

   std::valarray<std::complex<double>> val {0.0, 0.0, 0.0, 0.0};

   val += std::complex<double> {1., 0.} * FFV_FFS<Fe,Fe,VP,Fe,Fe,Ah>::value(indices1, indices2, context, discard_SM_contributions);
   val += std::complex<double> {1., 0.} * FFV_FFS<Fe,Fe,VP,Fe,Fe,hh>::value(indices1, indices2, context, discard_SM_contributions);
   val += std::complex<double> {1., 0.} * FFV_FFV<Fe,Fe,VP,Fe,Fe,VP>::value(indices1, indices2, context, discard_SM_contributions);
   val += std::complex<double> {1., 0.} * FFV_FFV<Fe,Fe,VP,Fe,Fe,VZ>::value(indices1, indices2, context, discard_SM_contributions);
   val += std::complex<double> {1., 0.} * FFV_SSF<Fe,Fe,VP,Hpm,Hpm,Fv>::value(indices1, indices2, context, discard_SM_contributions);
   val += std::complex<double> {1., 0.} * FFV_VSF<Fe,Fe,VP,VWm,Hpm,Fv>::value(indices1, indices2, context, discard_SM_contributions);
   val += std::complex<double> {1., 0.} * FFV_SSF<Fe,Fe,VP,Se,Se,Chi>::value(indices1, indices2, context, discard_SM_contributions);
   val += std::complex<double> {1., 0.} * FFV_SSF<Fe,Fe,VP,Se,Se,typename bar<Chi>::type>::value(indices1, indices2, context, discard_SM_contributions);
   val += std::complex<double> {1., 0.} * FFV_SVF<Fe,Fe,VP,Hpm,VWm,Fv>::value(indices1, indices2, context, discard_SM_contributions);
   val += std::complex<double> {1., 0.} * FFV_VVF<Fe,Fe,VP,VWm,VWm,Fv>::value(indices1, indices2, context, discard_SM_contributions);
   val += std::complex<double> {1., 0.} * FFV_FFS<Fe,Fe,VP,typename bar<Cha1>::type,typename bar<Cha1>::type,Sv>::value(indices1, indices2, context, discard_SM_contributions);

   return val;
}



template <typename Fj, typename Fi, typename V>
std::enable_if_t<
   std::is_same_v<Fj, Fe> && std::is_same_v<Fi, Fe> && std::is_same_v<V, VP>,
   std::valarray<std::complex<double>>
>
calculate_form_factors(   int generationIndex1, int generationIndex2,
   const MRSSM2_mass_eigenstates& model, bool discard_SM_contributions) {
return calculate_Fe_Fe_VP_form_factors(   generationIndex1, generationIndex2,
   model, discard_SM_contributions);
}
template std::valarray<std::complex<double>> calculate_form_factors<Fe,Fe,VP>(   int, int,
   const MRSSM2_mass_eigenstates&, bool);

} // namespace MRSSM2_FFV_form_factors
} // namespace flexiblesusy

namespace {

/**
* @defgroup FFVContributions FFV diagram massless photon contributions
* @brief Contributions to the processes Fe_I -> Fe_J gamma at the one-loop level.
*
* Diagram contributions are of the form:
* \Gamma^\mu = gamma^\mu F1 + I sigma^\mu\nu (p-pp)_\nu / (2*mj) F2
*                A1L gamma^\mu P_L + A1R gamma^\mu P_R
*              + A2L I sigma^\mu\nu (p-pp)_\nu / (2*mj) P_L
*              + A2R I sigma^\mu\nu (p-pp)_\nu / (2*mj) P_R
*/


/// emit massless vector boson from the internal scalar line
template <class Fj, class Fi, class V, class SA, class SB, class F>
std::valarray<std::complex<double>> FFV_SSF<Fj, Fi, V, SA, SB, F>::value(
   const typename field_indices<Fj>::type& indices_in,
   const typename field_indices<Fi>::type& indices_out,
   context_base const& context,
   bool discard_SM_contributions)
{
   static_assert(
      std::is_same_v<SA, SB>,
      "Internal scalars in the FFV_SSF instantiation must be of the same type."
   );

   using VertexFBarFjSBar = Vertex<typename F::lorentz_conjugate, typename SA::lorentz_conjugate, Fj>;
   using VertexFiBarFS    = Vertex<typename Fi::lorentz_conjugate, SB, F>;
   using VertexSBarSVBar  = Vertex<typename SB::lorentz_conjugate, SA, typename V::lorentz_conjugate>;

   // masses of external fermions
   const auto mj = context.mass<Fj>(indices_in);
   const auto mi = context.mass<Fi>(indices_out);

   std::valarray<std::complex<double>> res{0.0, 0.0, 0.0, 0.0};

   // loop over all possible particle generations attached to both vertices
   for (const auto& indexIn: index_range<VertexFBarFjSBar>()) {
      for (const auto& indexOut: index_range<VertexFiBarFS>()) {

         // cycle if generations of external fermions  are different then requested
         const auto jFieldIndices = VertexFBarFjSBar::template indices_of_field<2>(indexIn);
         const auto iFieldIndices = VertexFiBarFS::template indices_of_field<0>(indexOut);

         if (jFieldIndices != indices_in || iFieldIndices != indices_out) {
            continue;
         }

         // match indices of the fermion in the loop
         const auto fermionFieldIndicesIn = VertexFBarFjSBar::template indices_of_field<0>(indexIn);
         const auto fermionFieldIndicesOut = VertexFiBarFS::template indices_of_field<2>(indexOut);

         if (fermionFieldIndicesIn != fermionFieldIndicesOut) {
            continue;
         }

         // match indices of the scalar in the loop
         const auto scalarFieldIndicesIn = VertexFBarFjSBar::template indices_of_field<1>(indexIn);
         const auto scalarIndicesOut = VertexFiBarFS::template indices_of_field<1>(indexOut);

         if (scalarFieldIndicesIn != scalarIndicesOut) {
            continue;
         }

         if (discard_SM_contributions) {
            if (cxx_diagrams::isSMField<SA>(scalarFieldIndicesIn) &&
                cxx_diagrams::isSMField<F>(fermionFieldIndicesIn)) {
               continue;
            }
         }

         const auto vertexIn = VertexFBarFjSBar::evaluate(indexIn, context);
         const auto vertexOut = VertexFiBarFS::evaluate(indexOut, context);

         const auto indexEmit = concatenate(scalarFieldIndicesIn, scalarFieldIndicesIn);
         const auto vertexEmit = VertexSBarSVBar::evaluate(indexEmit, context);

         const auto mS = context.mass<SA>(scalarFieldIndicesIn);
         const auto mF = context.mass<F>(fermionFieldIndicesIn);
         const auto x = Power2(mF/mS);
         const auto LoopFunctionA = FA(x);
         const auto LoopFunctionB = FB(x);
         const auto LoopFunctionC = FC(x);

         // SARAH returns this coupling as some expression and momenta difference:
         // incoming momenta of one scalar comes with + sign,
         // incoming momenta of the other one comes with - sign.
         // + sign corresponds to the first argument of value(*,*).
         // - sign corresponds to the second argument of value(*,*).
         // Number, passed to the argument, is the position of the scalar in
         // the coupling. Example:
         //    Vertex<typename SB::lorentz_conjugate, SA, *>
         // and
         //    value(0, 1)
         // means, that we get C if vertex (in pseudolanguage) is given as
         // i C * (momenta(SBconj) - momenta(SA))
         std::complex<double> vector_boson_coupling {vertexEmit.value(0, 1)};

         // eq. 15 of hep-ph/9510309 (possibly with different sign)
         const std::complex<double> A1L =
            - 1./18. * vertexOut.right() * vertexIn.left() * LoopFunctionA;
         // eq. 16 of hep-ph/9510309 (possibly with different sign)
         const std::complex<double> A2L =
            - vertexOut.left() * vertexIn.right() * LoopFunctionB/12.
            - vertexOut.left()* vertexIn.left() * mF/mj * LoopFunctionC/3.
            - mi/mj * vertexOut.right() * vertexIn.left() * LoopFunctionB/12.;
         // eq. 15 & 16 of hep-ph/9510309 after replacement L <-> R (possibly with different sign)
         const std::complex<double> A1R =
            - 1./18. * vertexOut.left() * vertexIn.right() * LoopFunctionA;
         const std::complex<double> A2R =
            - vertexOut.right() * vertexIn.left() * LoopFunctionB/12.
            - vertexOut.right()* vertexIn.right() * mF/mj * LoopFunctionC/3.
            - mi/mj * vertexOut.left() * vertexIn.right() * LoopFunctionB/12.;

         const std::complex<double> massFactor = Power2(1.0/mS);

         res += oneOver32PiSqr * vector_boson_coupling * massFactor
            * std::valarray<std::complex<double>> {A1L, A1R, A2L, A2R};
      }
   }

   return res;
}


/// emit massless vector boson from the internal fermion line
template <class Fj, class Fi, class V, class FA, class FB, class S>
std::valarray<std::complex<double>> FFV_FFS<Fj, Fi, V, FA, FB, S>::value(
   const typename field_indices<Fj>::type& indices_in,
   const typename field_indices<Fi>::type& indices_out,
   context_base const& context,
   bool discard_SM_contributions)
{
   static_assert(
      std::is_same_v<FA, FB>,
      "Internal fermions in the FFV_FFS instantiation must be of the same type."
   );

   using VertexFBarFjSBar = Vertex<typename S::lorentz_conjugate, typename FA::lorentz_conjugate, Fj>;
   using VertexFiBarFS    = Vertex<typename Fi::lorentz_conjugate, FB, S>;
   using VertexFBarFVBar  = Vertex<typename FB::lorentz_conjugate, FA, typename V::lorentz_conjugate>;

   // masses of external fermions
   const auto mj = context.mass<Fj>(indices_in);
   const auto mi = context.mass<Fi>(indices_out);

   std::valarray<std::complex<double>> res{0.0, 0.0, 0.0, 0.0};

   // loop over all possible particle generations attached to both vertices
   for (const auto& indexIn: index_range<VertexFBarFjSBar>()) {
      for (const auto& indexOut: index_range<VertexFiBarFS>()) {

         // cycle if generations of external fermions are different then requested
         const auto jFieldIndices = VertexFBarFjSBar::template indices_of_field<2>(indexIn);
         const auto iFieldIndices = VertexFiBarFS::template indices_of_field<0>(indexOut);

         if (jFieldIndices != indices_in || iFieldIndices != indices_out) {
            continue;
         }

         const auto fermionFieldIndicesIn = VertexFBarFjSBar::template indices_of_field<1>(indexIn);
         const auto fermionFieldIndicesOut = VertexFiBarFS::template indices_of_field<1>(indexOut);

         if (fermionFieldIndicesIn != fermionFieldIndicesOut) {
            continue;
         }

         const auto scalarFieldIndicesIn = VertexFBarFjSBar::template indices_of_field<0>(indexIn);
         const auto scalarIndicesOut = VertexFiBarFS::template indices_of_field<2>(indexOut);

         if (scalarFieldIndicesIn != scalarIndicesOut) {
            continue;
         }

         if (discard_SM_contributions) {
            if (cxx_diagrams::isSMField<S>(scalarFieldIndicesIn) &&
                cxx_diagrams::isSMField<FA>(fermionFieldIndicesIn)) {
               continue;
            }
         }

         const auto vertexIn = VertexFBarFjSBar::evaluate(indexIn, context);
         const auto vertexOut = VertexFiBarFS::evaluate(indexOut, context);

         const auto indexEmit = concatenate(fermionFieldIndicesIn, fermionFieldIndicesIn);
         const auto vertexEmit = VertexFBarFVBar::evaluate(indexEmit, context);

         const auto mF = context.mass<FA>(fermionFieldIndicesIn);
         const auto mS = context.mass<S>(scalarFieldIndicesIn);
         const auto x = Power2(mF/mS);
         const auto LoopFunctionD = FD(x);
         const auto LoopFunctionE = FE(x);
         const auto LoopFunctionF = FF(x);

         const std::complex<double> vector_boson_coupling{vertexEmit.left()};

         // eq. 18 of hep-ph/9510309 (possibly with different sign)
         const std::complex<double> A1L =
            - 1./18. * vertexOut.right() * vertexIn.left() * LoopFunctionD;
         // eq. 19 of hep-ph/9510309 (possibly with different sign)
         const std::complex<double> A2L =
            - vertexOut.left() * vertexIn.right() * LoopFunctionE/12.0
            - vertexOut.left()* vertexIn.left() * mF/mj * LoopFunctionF * 2./3.
            - mi/mj * vertexOut.right() * vertexIn.left() * LoopFunctionE/12.0;
         // eq. 18 & 18 of hep-ph/9510309 after replacement L <-> R (possibly with different sign)
         const std::complex<double> A1R =
            - 1./18. * vertexOut.left() * vertexIn.right() * LoopFunctionD;
         const std::complex<double> A2R =
            - vertexOut.right() * vertexIn.left() * LoopFunctionE/12.0
            - vertexOut.right()* vertexIn.right() * mF/mj * LoopFunctionF * 2./3.
            - mi/mj * vertexOut.left() * vertexIn.right() * LoopFunctionE/12.0;

         const std::complex<double> massFactor = Power2(1.0/mS);

         res += oneOver32PiSqr * vector_boson_coupling * massFactor
            * std::valarray<std::complex<double>> {A1L, A1R, A2L, A2R};
      }
   }

   return res;
}


/// emit massless vector boson from the internal vector line in a fermion-vector loop
template <class Fj, class Fi, class P, class VA, class VB, class F>
std::valarray<std::complex<double>> FFV_VVF<Fj, Fi, P, VA, VB, F>::value(
   const typename field_indices<Fj>::type& indices_in,
   const typename field_indices<Fi>::type& indices_out,
   context_base const& context,
   bool discard_SM_contributions)
{
   static_assert(
      std::is_same_v<VA, VB>,
      "Internal vectors in the FFV_VVF instantiation must be of the same type."
   );

   // use P for external massless vector so as not to confuse with the internal vectors
   using VertexFBarFjVBar = Vertex<typename F::lorentz_conjugate, typename VA::lorentz_conjugate, Fj>;
   using VertexFiBarFV = Vertex<typename Fi::lorentz_conjugate, VB, F>;
   using VertexVBarVPBar = Vertex<typename VB::lorentz_conjugate, VA, typename P::lorentz_conjugate>;

   // masses of external fermions
   const auto mj = context.mass<Fj>(indices_in);
   const auto mi = context.mass<Fi>(indices_out);

   std::valarray<std::complex<double>> res{0.0, 0.0, 0.0, 0.0};

   // loop over all possible particle generations attached to both vertices
   for (const auto& indexIn: index_range<VertexFBarFjVBar>()) {
      for (const auto& indexOut: index_range<VertexFiBarFV>()) {

         // cycle if generations of external fermions are different then requested
         const auto jFieldIndices = VertexFBarFjVBar::template indices_of_field<2>(indexIn);
         const auto iFieldIndices = VertexFiBarFV::template indices_of_field<0>(indexOut);

         if (jFieldIndices != indices_in || iFieldIndices != indices_out) {
            continue;
         }

         // match indices of the fermion in the loop
         const auto fermionFieldIndicesIn = VertexFBarFjVBar::template indices_of_field<0>(indexIn);
         const auto fermionFieldIndicesOut = VertexFiBarFV::template indices_of_field<2>(indexOut);

         if (fermionFieldIndicesIn != fermionFieldIndicesOut) {
            continue;
         }

         // match indices of the vector in the loop
         const auto vectorFieldIndicesIn = VertexFBarFjVBar::template indices_of_field<1>(indexIn);
         const auto vectorFieldIndicesOut = VertexFiBarFV::template indices_of_field<1>(indexOut);

         if (vectorFieldIndicesIn != vectorFieldIndicesOut) {
            continue;
         }

         // ignore contribution if both vector and fermion are SM particles
         if (discard_SM_contributions) {
            if (cxx_diagrams::isSMField<VA>(vectorFieldIndicesIn) &&
                cxx_diagrams::isSMField<F>(fermionFieldIndicesIn)) {
               continue;
            }
         }

         const auto vertexIn = VertexFBarFjVBar::evaluate(indexIn, context);
         const auto vertexOut = VertexFiBarFV::evaluate(indexOut, context);

         const auto indexEmit = concatenate(vectorFieldIndicesIn, vectorFieldIndicesIn);
         const auto vertexEmit = VertexVBarVPBar::evaluate(indexEmit, context);

         const auto mV = context.mass<VA>(vectorFieldIndicesIn);
         const auto mF = context.mass<F>(fermionFieldIndicesIn);
         const auto x = Power2(mF/mV);
         const auto LoopFunctionH = FH(x);
         const auto LoopFunctionJ = FJ(x);
         const auto LoopFunctionK = FK(x);

         // triple gauge boson vertex
         // need to check that the sign and function call are correct
         const std::complex<double> vector_boson_coupling {-vertexEmit.value(cxx_diagrams::TripleVectorVertex::even_permutation{})};

         const std::complex<double> A1L =
            vertexIn.left() * vertexOut.left() * LoopFunctionH / 18.;
         const std::complex<double> A2L =
            3. * mF/mj * vertexIn.left() * vertexOut.right() * LoopFunctionK
            + mi/mj * vertexIn.left() * vertexOut.left() * LoopFunctionJ / 6.
            + vertexIn.right() * vertexOut.right() * LoopFunctionJ / 6.;
         const std::complex<double> A1R =
            vertexIn.right() * vertexOut.right() * LoopFunctionH / 18.;
         const std::complex<double> A2R =
            3. * mF/mj * vertexIn.right() * vertexOut.left() * LoopFunctionK
            + mi/mj * vertexIn.right() * vertexOut.right() * LoopFunctionJ / 6.
            + vertexIn.left() * vertexOut.left() * LoopFunctionJ / 6.;

         const std::complex<double> massFactor = Power2(1.0/mV);

         res += oneOver32PiSqr * vector_boson_coupling * massFactor * std::valarray<std::complex<double>> {A1L,A1R,A2L,A2R};
      }
   }

   return res;
}


/// emit massless vector boson from the internal vector line in a fermion-vector loop
template <class Fj, class Fi, class P, class FA, class FB, class V>
std::valarray<std::complex<double>> FFV_FFV<Fj, Fi, P, FA, FB, V>::value(
   const typename field_indices<Fj>::type& indices_in,
   const typename field_indices<Fi>::type& indices_out,
   context_base const& context,
   bool discard_SM_contributions)
{
   static_assert(
      std::is_same_v<FA, FB>,
      "Internal fermions in the FFV_FFV instantiation must be of the same type."
   );

   // use P for external massless vector so as not to confuse with the internal vectors
   using VertexFBarFjVBar = Vertex<typename V::lorentz_conjugate, typename FA::lorentz_conjugate, Fj>;
   using VertexFiBarFV = Vertex<typename Fi::lorentz_conjugate, FB, V>;
   using VertexFBarFPBar = Vertex<typename FB::lorentz_conjugate, FA, typename P::lorentz_conjugate>;

   // masses of external fermions
   const auto mj = context.mass<Fj>(indices_in);
   const auto mi = context.mass<Fi>(indices_out);

   std::valarray<std::complex<double>> res{0.0, 0.0, 0.0, 0.0};

   // loop over all possible particle generations attached to both vertices
   for (const auto& indexIn: index_range<VertexFBarFjVBar>()) {
      for (const auto& indexOut: index_range<VertexFiBarFV>()) {

         // cycle if generations of external fermions are different then requested
         const auto jFieldIndices = VertexFBarFjVBar::template indices_of_field<2>(indexIn);
         const auto iFieldIndices = VertexFiBarFV::template indices_of_field<0>(indexOut);

         if (jFieldIndices != indices_in || iFieldIndices != indices_out) {
            continue;
         }

         // match indices of the fermion in the loop
         const auto fermionFieldIndicesIn = VertexFBarFjVBar::template indices_of_field<1>(indexIn);
         const auto fermionFieldIndicesOut = VertexFiBarFV::template indices_of_field<1>(indexOut);

         if (fermionFieldIndicesIn != fermionFieldIndicesOut) {
            continue;
         }

         // match indices of the vector in the loop
         const auto vectorFieldIndicesIn = VertexFBarFjVBar::template indices_of_field<0>(indexIn);
         const auto vectorFieldIndicesOut = VertexFiBarFV::template indices_of_field<2>(indexOut);

         if (vectorFieldIndicesIn != vectorFieldIndicesOut) {
            continue;
         }

         // ignore contribution if both vector and fermion are SM particles
         if (discard_SM_contributions) {
            if (cxx_diagrams::isSMField<FA>(fermionFieldIndicesIn) &&
                cxx_diagrams::isSMField<V>(vectorFieldIndicesIn)) {
               continue;
            }
         }

         const auto vertexIn = VertexFBarFjVBar::evaluate(indexIn, context);
         const auto vertexOut = VertexFiBarFV::evaluate(indexOut, context);

         const auto indexEmit = concatenate(fermionFieldIndicesIn, fermionFieldIndicesIn);
         const auto vertexEmit = VertexFBarFPBar::evaluate(indexEmit, context);

         const auto mV = context.mass<V>(vectorFieldIndicesIn);
         const auto mF = context.mass<FA>(fermionFieldIndicesIn);
         const auto x = Power2(mF/mV);
         const auto LoopFunctionC = FC(x);
         const auto LoopFunctionL = FL(x);
         const auto LoopFunctionM = FM(x);

         const std::complex<double> vector_boson_coupling{vertexEmit.left()};

         const std::complex<double> A1L = vertexIn.left() * vertexOut.left() * LoopFunctionL / 9.;
         const std::complex<double> A2L =
            - 4. * mF/mj * vertexIn.left() * vertexOut.right() * LoopFunctionC / 3.
            + mi/mj * vertexIn.left() * vertexOut.left() * LoopFunctionM / 3.
            + vertexIn.right() * vertexOut.right() * LoopFunctionM / 3.;
         const std::complex<double> A1R = vertexIn.right() * vertexOut.right() * LoopFunctionL / 9.;
         const std::complex<double> A2R =
            - 4. * mF/mj * vertexIn.right() * vertexOut.left() * LoopFunctionC / 3.
            + mi/mj * vertexIn.right() * vertexOut.right() * LoopFunctionM / 3.
            + vertexIn.left() * vertexOut.left() * LoopFunctionM / 3.;

         std::complex<double> massFactor = 0.0;

         if (std::isinf(x)) {
            massFactor = Power2(1.0/mF);
         } else {
            massFactor = Power2(1.0/mV);
         }

         res += oneOver32PiSqr * vector_boson_coupling * massFactor * std::valarray<std::complex<double>> {A1L,A1R,A2L,A2R};
      }
   }

   return res;
}


/// emit massless vector boson from the internal vector line in a fermion-vector loop
template <class Fj, class Fi, class P, class VA, class GB, class F>
std::valarray<std::complex<double>> FFV_VSF<Fj, Fi, P, VA, GB, F>::value(
   const typename field_indices<Fj>::type& indices_in,
   const typename field_indices<Fi>::type& indices_out,
   context_base const& context,
   bool discard_SM_contributions)
{
   // do not assert that VA and GB are the same
   // static_assert(
   //     std::is_same_v<VA, GB>,
   //     "Internal vectors in the FFV_VGF instantiation must be of the same type."
   // );

   // use P for external massless vector so as not to confuse with the internal vectors
   using VertexFBarFjVBar = Vertex<typename F::lorentz_conjugate, typename VA::lorentz_conjugate, Fj>;
   using VertexFiBarFG = Vertex<typename Fi::lorentz_conjugate, GB, F>;
   using VertexGBarVPBar = Vertex<typename GB::lorentz_conjugate, VA, typename P::lorentz_conjugate>;

   // mass of incoming fermion
   const auto mj = context.mass<Fj>(indices_in);

   std::valarray<std::complex<double>> res{0.0, 0.0, 0.0, 0.0};

   // loop over all possible particle generations attached to both vertices
   for (const auto& indexIn: index_range<VertexFBarFjVBar>()) {
      for (const auto& indexOut: index_range<VertexFiBarFG>()) {

         // cycle if generations of external fermions are different then requested
         const auto jFieldIndices = VertexFBarFjVBar::template indices_of_field<2>(indexIn);
         const auto iFieldIndices = VertexFiBarFG::template indices_of_field<0>(indexOut);

         if (jFieldIndices != indices_in || iFieldIndices != indices_out) {
            continue;
         }

         // match indices of the fermion in the loop
         const auto fermionFieldIndicesIn = VertexFBarFjVBar::template indices_of_field<0>(indexIn);
         const auto fermionFieldIndicesOut = VertexFiBarFG::template indices_of_field<2>(indexOut);

         if (fermionFieldIndicesIn != fermionFieldIndicesOut) {
            continue;
         }

         // cannot match indices if scalar is not would-be-goldstone corresponding to vector
         const auto vectorFieldIndicesIn = VertexFBarFjVBar::template indices_of_field<1>(indexIn);
         const auto goldstoneFieldIndicesOut = VertexFiBarFG::template indices_of_field<1>(indexOut);

         // if (vectorFieldIndicesIn != goldstoneFieldIndicesOut) {
         //    continue;
         // }

         // ignore contribution if both vector and fermion are SM particles
         if (discard_SM_contributions) {
            if (cxx_diagrams::isSMField<VA>(vectorFieldIndicesIn) &&
                cxx_diagrams::isSMField<GB>(goldstoneFieldIndicesOut) &&
                cxx_diagrams::isSMField<F>(fermionFieldIndicesIn)) {
               continue;
            }
         }

         const auto vertexIn = VertexFBarFjVBar::evaluate(indexIn, context);
         const auto vertexOut = VertexFiBarFG::evaluate(indexOut, context);

         const auto indexEmit = concatenate(vectorFieldIndicesIn, goldstoneFieldIndicesOut);
         const auto vertexEmit = VertexGBarVPBar::evaluate(indexEmit, context);

         const auto mV = context.mass<VA>(vectorFieldIndicesIn);
         const auto mG = context.mass<GB>(goldstoneFieldIndicesOut);
         const auto mF = context.mass<F>(fermionFieldIndicesIn);
         const auto x = Power2(mF/mV);
         const auto y = Power2(mF/mG);
         const auto LoopFunctionN = FN(x,y);

         // scalar-vector-vector vertex
         // check that definition and sign is correct
         const std::complex<double> vector_boson_coupling{-vertexEmit.value()};

         const std::complex<double> A1L = 0.;
         const std::complex<double> A2L = 1./mj * vertexIn.left() * vertexOut.left() * LoopFunctionN;
         const std::complex<double> A1R = 0.;
         const std::complex<double> A2R = 1./mj * vertexIn.right() * vertexOut.right() * LoopFunctionN;

         const std::complex<double> massFactor = Power2(1.0/mV);

         res += oneOver32PiSqr * vector_boson_coupling * massFactor * std::valarray<std::complex<double>> {A1L,A1R,A2L,A2R};
      }
   }

   return res;
}


/// emit massless vector boson from the internal vector line in a fermion-vector loop
template <class Fj, class Fi, class P, class GA, class VB, class F>
std::valarray<std::complex<double>> FFV_SVF<Fj, Fi, P, GA, VB, F>::value(
   const typename field_indices<Fj>::type& indices_in,
   const typename field_indices<Fi>::type& indices_out,
   context_base const& context,
   bool discard_SM_contributions)
{
   // do not assert that GA and VB are the same
   // static_assert(
   //     std::is_same_v<GA, VB>,
   //     "Internal vectors in the FFV_VGF instantiation must be of the same type."
   // );

   // use P for external massless vector so as not to confuse with the internal vectors
   using VertexFBarFjGBar = Vertex<typename F::lorentz_conjugate, typename GA::lorentz_conjugate, Fj>;
   using VertexFiBarFV = Vertex<typename Fi::lorentz_conjugate, VB, F>;
   using VertexVBarGPBar = Vertex<typename VB::lorentz_conjugate, GA, typename P::lorentz_conjugate>;

   // mass of incoming fermion
   const auto mj = context.mass<Fj>(indices_in);

   std::valarray<std::complex<double>> res{0.0, 0.0, 0.0, 0.0};

   // loop over all possible particle generations attached to both vertices
   for (const auto& indexIn: index_range<VertexFBarFjGBar>()) {
      for (const auto& indexOut: index_range<VertexFiBarFV>()) {

         // cycle if generations of external fermions are different then requested
         const auto jFieldIndices = VertexFBarFjGBar::template indices_of_field<2>(indexIn);
         const auto iFieldIndices = VertexFiBarFV::template indices_of_field<0>(indexOut);

         if (jFieldIndices != indices_in || iFieldIndices != indices_out) {
            continue;
         }

         // match indices of the fermion in the loop
         const auto fermionFieldIndicesIn = VertexFBarFjGBar::template indices_of_field<0>(indexIn);
         const auto fermionFieldIndicesOut = VertexFiBarFV::template indices_of_field<2>(indexOut);

         if (fermionFieldIndicesIn != fermionFieldIndicesOut) {
            continue;
         }

         // cannot match indices if scalar is not would-be-goldstone corresponding to vector
         const auto goldstoneFieldIndicesIn = VertexFBarFjGBar::template indices_of_field<1>(indexIn);
         const auto vectorFieldIndicesOut = VertexFiBarFV::template indices_of_field<1>(indexOut);

         // if (goldstoneFieldIndicesIn != vectorFieldIndicesOut) {
         //    continue;
         // }

         // ignore contribution if both vector and fermion are SM particles
         if (discard_SM_contributions) {
            if (cxx_diagrams::isSMField<GA>(goldstoneFieldIndicesIn) &&
                cxx_diagrams::isSMField<VB>(vectorFieldIndicesOut) &&
                cxx_diagrams::isSMField<F>(fermionFieldIndicesIn)) {
               continue;
            }
         }

         const auto vertexIn = VertexFBarFjGBar::evaluate(indexIn, context);
         const auto vertexOut = VertexFiBarFV::evaluate(indexOut, context);

         const auto indexEmit = concatenate(goldstoneFieldIndicesIn, vectorFieldIndicesOut);
         const auto vertexEmit = VertexVBarGPBar::evaluate(indexEmit, context);

         const auto mG = context.mass<GA>(goldstoneFieldIndicesIn);
         const auto mV = context.mass<VB>(vectorFieldIndicesOut);
         const auto mF = context.mass<F>(fermionFieldIndicesIn);
         const auto x = Power2(mF/mV);
         const auto y = Power2(mF/mG);
         const auto LoopFunctionN = FN(x,y);

         // scalar-vector-vector vertex
         // check that definition and sign is correct
         const std::complex<double> vector_boson_coupling{-vertexEmit.value()};

         const std::complex<double> A1L = 0.;
         const std::complex<double> A2L = 1./mj * vertexIn.left() * vertexOut.right() * LoopFunctionN;
         const std::complex<double> A1R = 0.;
         const std::complex<double> A2R = 1./mj * vertexIn.right() * vertexOut.left() * LoopFunctionN;

         const std::complex<double> massFactor = Power2(1.0/mV);

         res += oneOver32PiSqr * vector_boson_coupling * massFactor * std::valarray<std::complex<double>> {A1L,A1R,A2L,A2R};
      }
   }

   return res;
}

} // anonymous namespace
