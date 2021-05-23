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
 * @file SplitMSSM_decays.hpp
 *
 * @brief contains class for calculating particle decays
 *
 * This file was generated with FlexibleSUSY 2.6.0 and SARAH 4.14.4 .
 */

#ifndef SplitMSSM_DECAYS_H
#define SplitMSSM_DECAYS_H

#include "SplitMSSM_decay_table.hpp"
#include "SplitMSSM_mass_eigenstates.hpp"
#include "SplitMSSM_mass_eigenstates_decoupling_scheme.hpp"
#include "cxx_qft/SplitMSSM_qft.hpp"
#include "SplitMSSM_decay_amplitudes.hpp"
#include "decays/flexibledecay_problems.hpp"
#include "lowe.h"
#include "wrappers.hpp"
#include "error.hpp"
#include "physical_input.hpp"
#include "decays/flexibledecay_settings.hpp"

namespace flexiblesusy {

template <typename Field1, typename Field2>
constexpr std::enable_if_t<!std::is_same<Field1, Field2>::value, double>
final_state_symmetry_factor(typename SplitMSSM_cxx_diagrams::field_indices<Field1>::type const& idx1,
                            typename SplitMSSM_cxx_diagrams::field_indices<Field2>::type const& idx2)
{
   return 1.;
}

template <typename Field1, typename Field2>
std::enable_if_t<std::is_same<Field1, Field2>::value, double>
final_state_symmetry_factor(typename SplitMSSM_cxx_diagrams::field_indices<Field1>::type const& idx1,
                            typename SplitMSSM_cxx_diagrams::field_indices<Field2>::type const& idx2)
{
   if (boost::range::equal(idx1, idx2)) {
      return 0.5;
   }
   else {
      return 1.;
   }
}

class SplitMSSM_decays {
public:
   SplitMSSM_decays() = default;
   SplitMSSM_decays(SplitMSSM_mass_eigenstates model_, softsusy::QedQcd const& qedqcd_,
         Physical_input const& physical_input_,
         FlexibleDecay_settings const& flexibledecay_settings_)
      : model(model_)
      , qedqcd(qedqcd_)
      , physical_input(physical_input_)
      , flexibledecay_settings(flexibledecay_settings_)
      {}
   SplitMSSM_decays(const SplitMSSM_decays&) = default;
   SplitMSSM_decays(SplitMSSM_decays&&) = default;
   ~SplitMSSM_decays() = default;
   SplitMSSM_decays& operator=(const SplitMSSM_decays&) = default;
   SplitMSSM_decays& operator=(SplitMSSM_decays&&) = default;

   const SplitMSSM_decay_table& get_decay_table() const;
   const FlexibleDecay_problems& get_problems() const;

   void clear();
   void clear_problems();
   void calculate_decays();

   const Decays_list& get_hh_decays() const { return decay_table.get_hh_decays();
      }
   void calculate_hh_decays();

double partial_width_hh_to_VGVG(SplitMSSM_mass_eigenstates_interface*) const;
double partial_width_hh_to_VPVP(SplitMSSM_mass_eigenstates_interface*) const;
double partial_width_hh_to_VPVZ(SplitMSSM_mass_eigenstates_interface*) const;
double partial_width_hh_to_VZVZ(SplitMSSM_mass_eigenstates_interface*) const;
double partial_width_hh_to_conjVWpVWp(SplitMSSM_mass_eigenstates_interface*) const;
double partial_width_hh_to_barFvFv(SplitMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_GluGlu(SplitMSSM_mass_eigenstates_interface*) const;
double partial_width_hh_to_barFdFd(SplitMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_barFuFu(SplitMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_barFeFe(SplitMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_ChiChi(SplitMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_barChaCha(SplitMSSM_mass_eigenstates_interface*, int, int) const;

private:
   SplitMSSM_mass_eigenstates model{};
   softsusy::QedQcd qedqcd{};
   Physical_input physical_input;
   FlexibleDecay_settings flexibledecay_settings {};
   bool run_to_decay_particle_scale {true};
   SplitMSSM_decay_table decay_table{};
   FlexibleDecay_problems problems{};

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   typename Decay_amplitude_type<FieldIn, FieldOut1, FieldOut2>::type
   calculate_amplitude(
      const SplitMSSM_cxx_diagrams::context_base&,
      const typename SplitMSSM_cxx_diagrams::field_indices<FieldIn>::type&,
      const typename SplitMSSM_cxx_diagrams::field_indices<FieldOut1>::type&,
      const typename SplitMSSM_cxx_diagrams::field_indices<FieldOut2>::type&) const;

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   double amplitude_squared(SplitMSSM_cxx_diagrams::context_base const& context,
                  typename SplitMSSM_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
                  typename SplitMSSM_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
                  typename SplitMSSM_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2) const;

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   double get_partial_width(
      const SplitMSSM_cxx_diagrams::context_base&,
      typename SplitMSSM_cxx_diagrams::field_indices<FieldIn>::type const&,
      typename SplitMSSM_cxx_diagrams::field_indices<FieldOut1>::type const&,
      typename SplitMSSM_cxx_diagrams::field_indices<FieldOut2>::type const&) const;
};

template<>
Decay_amplitude_SVV SplitMSSM_decays::calculate_amplitude<SplitMSSM_cxx_diagrams::fields::hh, SplitMSSM_cxx_diagrams::fields::VG, SplitMSSM_cxx_diagrams::fields::VG>(const SplitMSSM_cxx_diagrams::context_base&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::hh >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::VG >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SVV SplitMSSM_decays::calculate_amplitude<SplitMSSM_cxx_diagrams::fields::hh, SplitMSSM_cxx_diagrams::fields::VP, SplitMSSM_cxx_diagrams::fields::VP>(const SplitMSSM_cxx_diagrams::context_base&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::hh >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::VP >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SVV SplitMSSM_decays::calculate_amplitude<SplitMSSM_cxx_diagrams::fields::hh, SplitMSSM_cxx_diagrams::fields::VP, SplitMSSM_cxx_diagrams::fields::VZ>(const SplitMSSM_cxx_diagrams::context_base&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::hh >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::VP >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV SplitMSSM_decays::calculate_amplitude<SplitMSSM_cxx_diagrams::fields::hh, SplitMSSM_cxx_diagrams::fields::VZ, SplitMSSM_cxx_diagrams::fields::VZ>(const SplitMSSM_cxx_diagrams::context_base&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::hh >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::VZ >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV SplitMSSM_decays::calculate_amplitude<SplitMSSM_cxx_diagrams::fields::hh, typename SplitMSSM_cxx_diagrams::fields::conj<SplitMSSM_cxx_diagrams::fields::VWp>::type, SplitMSSM_cxx_diagrams::fields::VWp>(const SplitMSSM_cxx_diagrams::context_base&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::hh >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<typename SplitMSSM_cxx_diagrams::fields::conj<SplitMSSM_cxx_diagrams::fields::VWp>::type >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::VWp >::type&) const;

template<>
Decay_amplitude_SFF SplitMSSM_decays::calculate_amplitude<SplitMSSM_cxx_diagrams::fields::hh, typename SplitMSSM_cxx_diagrams::fields::bar<SplitMSSM_cxx_diagrams::fields::Fv>::type, SplitMSSM_cxx_diagrams::fields::Fv>(const SplitMSSM_cxx_diagrams::context_base&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::hh >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<typename SplitMSSM_cxx_diagrams::fields::bar<SplitMSSM_cxx_diagrams::fields::Fv>::type >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::Fv >::type&) const;

template<>
Decay_amplitude_SFF SplitMSSM_decays::calculate_amplitude<SplitMSSM_cxx_diagrams::fields::hh, SplitMSSM_cxx_diagrams::fields::Glu, SplitMSSM_cxx_diagrams::fields::Glu>(const SplitMSSM_cxx_diagrams::context_base&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::hh >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::Glu >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::Glu >::type&) const;

template<>
Decay_amplitude_SFF SplitMSSM_decays::calculate_amplitude<SplitMSSM_cxx_diagrams::fields::hh, typename SplitMSSM_cxx_diagrams::fields::bar<SplitMSSM_cxx_diagrams::fields::Fd>::type, SplitMSSM_cxx_diagrams::fields::Fd>(const SplitMSSM_cxx_diagrams::context_base&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::hh >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<typename SplitMSSM_cxx_diagrams::fields::bar<SplitMSSM_cxx_diagrams::fields::Fd>::type >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SFF SplitMSSM_decays::calculate_amplitude<SplitMSSM_cxx_diagrams::fields::hh, typename SplitMSSM_cxx_diagrams::fields::bar<SplitMSSM_cxx_diagrams::fields::Fu>::type, SplitMSSM_cxx_diagrams::fields::Fu>(const SplitMSSM_cxx_diagrams::context_base&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::hh >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<typename SplitMSSM_cxx_diagrams::fields::bar<SplitMSSM_cxx_diagrams::fields::Fu>::type >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::Fu >::type&) const;

template<>
Decay_amplitude_SFF SplitMSSM_decays::calculate_amplitude<SplitMSSM_cxx_diagrams::fields::hh, typename SplitMSSM_cxx_diagrams::fields::bar<SplitMSSM_cxx_diagrams::fields::Fe>::type, SplitMSSM_cxx_diagrams::fields::Fe>(const SplitMSSM_cxx_diagrams::context_base&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::hh >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<typename SplitMSSM_cxx_diagrams::fields::bar<SplitMSSM_cxx_diagrams::fields::Fe>::type >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::Fe >::type&) const;

template<>
Decay_amplitude_SFF SplitMSSM_decays::calculate_amplitude<SplitMSSM_cxx_diagrams::fields::hh, SplitMSSM_cxx_diagrams::fields::Chi, SplitMSSM_cxx_diagrams::fields::Chi>(const SplitMSSM_cxx_diagrams::context_base&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::hh >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::Chi >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF SplitMSSM_decays::calculate_amplitude<SplitMSSM_cxx_diagrams::fields::hh, typename SplitMSSM_cxx_diagrams::fields::bar<SplitMSSM_cxx_diagrams::fields::Cha>::type, SplitMSSM_cxx_diagrams::fields::Cha>(const SplitMSSM_cxx_diagrams::context_base&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::hh >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<typename SplitMSSM_cxx_diagrams::fields::bar<SplitMSSM_cxx_diagrams::fields::Cha>::type >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::Cha >::type&) const;


template<typename FieldIn, typename FieldOut1, typename FieldOut2>
double
SplitMSSM_decays::amplitude_squared(SplitMSSM_cxx_diagrams::context_base const& context,
                  typename SplitMSSM_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
                  typename SplitMSSM_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
                  typename SplitMSSM_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2) const
{

   const auto mat_elem = calculate_amplitude<FieldIn, FieldOut1, FieldOut2>(
      context, indexIn, indexOut1, indexOut2);
   return mat_elem.square();
}

// returns a squared color generator for a 3 point amplitude with FieldIn, FieldOut1 and FieldOut2
// averaged over inital state colors
// the generator is guessed from color representations of FieldIn, FieldOut1 and FieldOut2
// This is not a bulletproof solution and might fail in general but is enough for
// decays of color singlets

// 1 -> 1, 1
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
SplitMSSM_cxx_diagrams::fields::is_singlet_v<FieldIn> &&
SplitMSSM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
SplitMSSM_cxx_diagrams::fields::is_singlet_v<FieldOut2>, double>
squared_color_generator() {return 1.;}

// 1 -> 3, 3bar
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
   SplitMSSM_cxx_diagrams::fields::is_singlet_v<FieldIn>
   &&
   (
   (SplitMSSM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
   SplitMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>)
   ||
   (SplitMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
   SplitMSSM_cxx_diagrams::fields::is_triplet_v<FieldOut2>)
   ), double
>
squared_color_generator() {return 3.;}

// 1 -> 8, 8
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
SplitMSSM_cxx_diagrams::fields::is_singlet_v<FieldIn> &&
SplitMSSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
SplitMSSM_cxx_diagrams::fields::is_octet_v<FieldOut2>, double>
squared_color_generator() {return 8.;}

// 3 -> 3, 1; 3bar -> 3bar, 1
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
(
SplitMSSM_cxx_diagrams::fields::is_triplet_v<FieldIn> &&
((SplitMSSM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
SplitMSSM_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(SplitMSSM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
SplitMSSM_cxx_diagrams::fields::is_triplet_v<FieldOut2>))
)
||
(
SplitMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldIn> &&
((SplitMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
SplitMSSM_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(SplitMSSM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
SplitMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>))
)
, double>
squared_color_generator() {return 1.;}

// 3 -> 3, 8; 3bar -> 3bar, 8
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
(
SplitMSSM_cxx_diagrams::fields::is_triplet_v<FieldIn> &&
((SplitMSSM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
SplitMSSM_cxx_diagrams::fields::is_octet_v<FieldOut2>) ||
(SplitMSSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
SplitMSSM_cxx_diagrams::fields::is_triplet_v<FieldOut2>))
)
||
(
SplitMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldIn> &&
((SplitMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
SplitMSSM_cxx_diagrams::fields::is_octet_v<FieldOut2>) ||
(SplitMSSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
SplitMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>))
)
, double>
squared_color_generator() {return 4.;}

// 8 -> 8, 1
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
SplitMSSM_cxx_diagrams::fields::is_octet_v<FieldIn> &&
((SplitMSSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
SplitMSSM_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(SplitMSSM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
SplitMSSM_cxx_diagrams::fields::is_octet_v<FieldOut2>))
, double>
squared_color_generator() {return 1.;}

// 8 -> 3, 3bar
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
   SplitMSSM_cxx_diagrams::fields::is_octet_v<FieldIn>
   &&
   (
   (SplitMSSM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
   SplitMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>)
   ||
   (SplitMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
   SplitMSSM_cxx_diagrams::fields::is_triplet_v<FieldOut2>)
   ), double
>
squared_color_generator() {return 1./2.;}

// generic decay of FieldIn -> FieldOut1 FieldOut2
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
double SplitMSSM_decays::get_partial_width(
   const SplitMSSM_cxx_diagrams::context_base& context,
   typename SplitMSSM_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
   typename SplitMSSM_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
   typename SplitMSSM_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2
   ) const
{

   // on-shell masses
   const double mIn = context.physical_mass<FieldIn>(indexIn);
   const double mOut1 = context.physical_mass<FieldOut1>(indexOut1);
   const double mOut2 = context.physical_mass<FieldOut2>(indexOut2);

   // is decay kinematically allowed?
   if(mIn < mOut1 + mOut2) {
      WARNING("Called kinematically forbidden decay");
      return 0.;
   }

   // phase space without symmetry factor
   const double ps = 1./(8.*Pi) * std::sqrt(KallenLambda(1., Sqr(mOut1/mIn), Sqr(mOut2/mIn)));

   // phase space symmetry factor
   const double ps_symmetry =
      final_state_symmetry_factor<FieldOut1, FieldOut2>(indexOut1, indexOut2);

   // color factor
   constexpr double color_factor = squared_color_generator<FieldIn, FieldOut1, FieldOut2>();

   // matrix element squared
   const auto mat_elem_sq = amplitude_squared<FieldIn, FieldOut1, FieldOut2>(
      context, indexIn, indexOut1, indexOut2);

   // flux * phase space factor * symmetry factor * color factor * |matrix element|^2
   const auto result = 0.5/mIn * ps * ps_symmetry * color_factor * mat_elem_sq;

   return result;
}

template <>
double SplitMSSM_decays::get_partial_width<SplitMSSM_cxx_diagrams::fields::hh,SplitMSSM_cxx_diagrams::fields::VZ,SplitMSSM_cxx_diagrams::fields::VZ >(const SplitMSSM_cxx_diagrams::context_base&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::hh >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::VZ >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::VZ >::type&) const;
template <>
double SplitMSSM_decays::get_partial_width<SplitMSSM_cxx_diagrams::fields::hh,typename SplitMSSM_cxx_diagrams::fields::conj<SplitMSSM_cxx_diagrams::fields::VWp>::type,SplitMSSM_cxx_diagrams::fields::VWp >(const SplitMSSM_cxx_diagrams::context_base&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::hh >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<typename SplitMSSM_cxx_diagrams::fields::conj<SplitMSSM_cxx_diagrams::fields::VWp>::type >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::VWp >::type&) const;
template <>
double SplitMSSM_decays::get_partial_width<SplitMSSM_cxx_diagrams::fields::hh,SplitMSSM_cxx_diagrams::fields::VG,SplitMSSM_cxx_diagrams::fields::VG >(const SplitMSSM_cxx_diagrams::context_base&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::hh >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::VG >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::VG >::type&) const;
template <>
double SplitMSSM_decays::get_partial_width<SplitMSSM_cxx_diagrams::fields::hh,SplitMSSM_cxx_diagrams::fields::VP,SplitMSSM_cxx_diagrams::fields::VP >(const SplitMSSM_cxx_diagrams::context_base&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::hh >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::VP >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::VP >::type&) const;
template <>
double SplitMSSM_decays::get_partial_width<SplitMSSM_cxx_diagrams::fields::hh,SplitMSSM_cxx_diagrams::fields::VP,SplitMSSM_cxx_diagrams::fields::VZ >(const SplitMSSM_cxx_diagrams::context_base&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::hh >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::VP >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::VZ >::type&) const;
template <>
double SplitMSSM_decays::get_partial_width<SplitMSSM_cxx_diagrams::fields::hh,typename SplitMSSM_cxx_diagrams::fields::bar<SplitMSSM_cxx_diagrams::fields::Fu>::type,SplitMSSM_cxx_diagrams::fields::Fu >(const SplitMSSM_cxx_diagrams::context_base&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::hh >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<typename SplitMSSM_cxx_diagrams::fields::bar<SplitMSSM_cxx_diagrams::fields::Fu>::type >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::Fu >::type&) const;
template <>
double SplitMSSM_decays::get_partial_width<SplitMSSM_cxx_diagrams::fields::hh,typename SplitMSSM_cxx_diagrams::fields::bar<SplitMSSM_cxx_diagrams::fields::Fd>::type,SplitMSSM_cxx_diagrams::fields::Fd >(const SplitMSSM_cxx_diagrams::context_base&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::hh >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<typename SplitMSSM_cxx_diagrams::fields::bar<SplitMSSM_cxx_diagrams::fields::Fd>::type >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::Fd >::type&) const;
template <>
double SplitMSSM_decays::get_partial_width<SplitMSSM_cxx_diagrams::fields::hh,typename SplitMSSM_cxx_diagrams::fields::bar<SplitMSSM_cxx_diagrams::fields::Fe>::type,SplitMSSM_cxx_diagrams::fields::Fe >(const SplitMSSM_cxx_diagrams::context_base&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::hh >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<typename SplitMSSM_cxx_diagrams::fields::bar<SplitMSSM_cxx_diagrams::fields::Fe>::type >::type&, const typename SplitMSSM_cxx_diagrams::field_indices<SplitMSSM_cxx_diagrams::fields::Fe >::type&) const;

} // namespace flexiblesusy

#endif
