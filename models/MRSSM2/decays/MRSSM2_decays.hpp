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
 * @file MRSSM2_decays.hpp
 *
 * @brief contains class for calculating particle decays
 *
 * This file was generated with FlexibleSUSY 2.7.1 and SARAH 4.14.5 .
 */

#ifndef MRSSM2_DECAYS_H
#define MRSSM2_DECAYS_H

#include "MRSSM2_decay_table.hpp"
#include "MRSSM2_mass_eigenstates.hpp"
#include "MRSSM2_mass_eigenstates_decoupling_scheme.hpp"
#include "cxx_qft/MRSSM2_qft.hpp"
#include "MRSSM2_decay_amplitudes.hpp"
#include "decays/flexibledecay_problems.hpp"
#include "lowe.h"
#include "wrappers.hpp"
#include "error.hpp"
#include "physical_input.hpp"
#include "decays/flexibledecay_settings.hpp"

namespace flexiblesusy {

template <typename Field1, typename Field2>
constexpr std::enable_if_t<!std::is_same<Field1, Field2>::value, double>
final_state_symmetry_factor(typename MRSSM2_cxx_diagrams::field_indices<Field1>::type const& idx1,
                            typename MRSSM2_cxx_diagrams::field_indices<Field2>::type const& idx2)
{
   return 1.;
}

template <typename Field1, typename Field2>
std::enable_if_t<std::is_same<Field1, Field2>::value, double>
final_state_symmetry_factor(typename MRSSM2_cxx_diagrams::field_indices<Field1>::type const& idx1,
                            typename MRSSM2_cxx_diagrams::field_indices<Field2>::type const& idx2)
{
   if (boost::range::equal(idx1, idx2)) {
      return 0.5;
   }
   else {
      return 1.;
   }
}

class MRSSM2_decays {
public:
   MRSSM2_decays() = default;
   MRSSM2_decays(MRSSM2_mass_eigenstates model_, softsusy::QedQcd const& qedqcd_,
         Physical_input const& physical_input_,
         FlexibleDecay_settings const& flexibledecay_settings_)
      : model(model_)
      , qedqcd(qedqcd_)
      , physical_input(physical_input_)
      , flexibledecay_settings(flexibledecay_settings_)
      {}
   MRSSM2_decays(const MRSSM2_decays&) = default;
   MRSSM2_decays(MRSSM2_decays&&) = default;
   ~MRSSM2_decays() = default;
   MRSSM2_decays& operator=(const MRSSM2_decays&) = default;
   MRSSM2_decays& operator=(MRSSM2_decays&&) = default;

   const MRSSM2_decay_table& get_decay_table() const;
   const FlexibleDecay_problems& get_problems() const;

   void clear();
   void clear_problems();
   void calculate_decays();

   const Decays_list& get_hh_decays(int i) const { return decay_table.
      get_hh_decays(i); }
   const Decays_list& get_Ah_decays(int i) const { return decay_table.
      get_Ah_decays(i); }
   const Decays_list& get_Hpm_decays(int i) const { return decay_table.
      get_Hpm_decays(i); }
   const Decays_list& get_sigmaO_decays() const { return decay_table.
      get_sigmaO_decays(); }
   const Decays_list& get_phiO_decays() const { return decay_table.get_phiO_decays
      (); }
   void calculate_hh_decays();
   void calculate_Ah_decays();
   void calculate_Hpm_decays();
   void calculate_sigmaO_decays();
   void calculate_phiO_decays();

double partial_width_hh_to_SRdpSRum(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_SRdpHpm(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_SRdpconjSRdp(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_SRumconjSRum(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_SRumconjHpm(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_sigmaOsigmaO(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_sigmaOphiO(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_phiOphiO(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_SdconjSd(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SvconjSv(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SuconjSu(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SeconjSe(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_hhhh(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_hhAh(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_hhRh(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_hhconjRh(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_AhAh(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_AhRh(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_AhconjRh(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_RhRh(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_RhconjRh(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_HpmconjSRum(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_HpmconjHpm(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_conjSRdpconjSRum(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_conjSRdpconjHpm(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_conjRhconjRh(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_sigmaOVG(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_phiOVG(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_hhVP(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_AhVP(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_RhVP(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_conjRhVP(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_hhVZ(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_AhVZ(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_RhVZ(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_conjRhVZ(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_SRdpVWm(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_conjSRumVWm(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_conjHpmVWm(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_SRumconjVWm(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_HpmconjVWm(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_conjSRdpconjVWm(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VGVG(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VPVP(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VPVZ(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VZVZ(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_conjVWmVWm(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_GluGlu(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_barGluGlu(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_barFvFv(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_ChiChi(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barChiChi(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_Cha1Cha2(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barCha1Cha1(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barCha2Cha2(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barFeFe(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barFdFd(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barFuFu(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barGlubarGlu(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_barChibarChi(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barCha1barCha2(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SRdpSRum(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_SRdpHpm(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_SRdpconjSRdp(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_SRumconjSRum(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_SRumconjHpm(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_sigmaOsigmaO(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_sigmaOphiO(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_phiOphiO(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_SdconjSd(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SvconjSv(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SuconjSu(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SeconjSe(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_hhhh(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_hhAh(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_hhRh(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_hhconjRh(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_AhAh(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_AhRh(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_AhconjRh(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_RhRh(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_RhconjRh(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_HpmconjSRum(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_HpmconjHpm(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_conjSRdpconjSRum(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_conjSRdpconjHpm(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_conjRhconjRh(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_sigmaOVG(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_phiOVG(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_hhVP(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_AhVP(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_RhVP(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_conjRhVP(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_hhVZ(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_AhVZ(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_RhVZ(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_conjRhVZ(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_SRdpVWm(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_conjSRumVWm(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_conjHpmVWm(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_SRumconjVWm(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_HpmconjVWm(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_conjSRdpconjVWm(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VGVG(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VPVP(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VPVZ(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VZVZ(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_conjVWmVWm(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_GluGlu(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_barGluGlu(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_barFvFv(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_ChiChi(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barChiChi(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_Cha1Cha2(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barCha1Cha1(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barCha2Cha2(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barFeFe(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barFdFd(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barFuFu(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barGlubarGlu(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_barChibarChi(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barCha1barCha2(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_SRumhh(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_SRumAh(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_SRumRh(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_SRumconjRh(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_SdconjSu(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_SeconjSv(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_hhHpm(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_hhconjSRdp(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_AhHpm(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_AhconjSRdp(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_RhHpm(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_RhconjSRdp(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_HpmconjRh(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_conjSRdpconjRh(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_SRumVP(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_Hpm_to_HpmVP(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_conjSRdpVP(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_Hpm_to_SRumVZ(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_Hpm_to_HpmVZ(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_conjSRdpVZ(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_Hpm_to_hhVWm(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_AhVWm(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_RhVWm(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_conjRhVWm(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_VPVWm(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_Hpm_to_VZVWm(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_Hpm_to_ChiCha2(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_barCha1Chi(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_barChiCha2(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_barFvFe(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_barFuFd(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_barChibarCha1(MRSSM2_mass_eigenstates_interface*, int, int, int) const;
double partial_width_sigmaO_to_phiOphiO(MRSSM2_mass_eigenstates_interface*) const;
double partial_width_sigmaO_to_phiOhh(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_sigmaO_to_phiOAh(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_sigmaO_to_phiORh(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_sigmaO_to_phiOconjRh(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_sigmaO_to_SdconjSd(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_sigmaO_to_SuconjSu(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_sigmaO_to_hhVG(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_sigmaO_to_AhVG(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_sigmaO_to_RhVG(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_sigmaO_to_conjRhVG(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_sigmaO_to_phiOVP(MRSSM2_mass_eigenstates_interface*) const;
double partial_width_sigmaO_to_phiOVZ(MRSSM2_mass_eigenstates_interface*) const;
double partial_width_sigmaO_to_VGVG(MRSSM2_mass_eigenstates_interface*) const;
double partial_width_sigmaO_to_VGVP(MRSSM2_mass_eigenstates_interface*) const;
double partial_width_sigmaO_to_VGVZ(MRSSM2_mass_eigenstates_interface*) const;
double partial_width_sigmaO_to_GluGlu(MRSSM2_mass_eigenstates_interface*) const;
double partial_width_sigmaO_to_GluChi(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_sigmaO_to_barChiGlu(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_sigmaO_to_barGluChi(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_sigmaO_to_barFdFd(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_sigmaO_to_barFuFu(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_sigmaO_to_barGlubarGlu(MRSSM2_mass_eigenstates_interface*) const;
double partial_width_sigmaO_to_barGlubarChi(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_phiO_to_sigmaOsigmaO(MRSSM2_mass_eigenstates_interface*) const;
double partial_width_phiO_to_sigmaOhh(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_phiO_to_sigmaOAh(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_phiO_to_sigmaORh(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_phiO_to_sigmaOconjRh(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_phiO_to_SdconjSd(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_phiO_to_SuconjSu(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_phiO_to_hhVG(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_phiO_to_AhVG(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_phiO_to_RhVG(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_phiO_to_conjRhVG(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_phiO_to_sigmaOVP(MRSSM2_mass_eigenstates_interface*) const;
double partial_width_phiO_to_sigmaOVZ(MRSSM2_mass_eigenstates_interface*) const;
double partial_width_phiO_to_VGVG(MRSSM2_mass_eigenstates_interface*) const;
double partial_width_phiO_to_VGVP(MRSSM2_mass_eigenstates_interface*) const;
double partial_width_phiO_to_VGVZ(MRSSM2_mass_eigenstates_interface*) const;
double partial_width_phiO_to_GluGlu(MRSSM2_mass_eigenstates_interface*) const;
double partial_width_phiO_to_GluChi(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_phiO_to_barChiGlu(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_phiO_to_barGluChi(MRSSM2_mass_eigenstates_interface*, int) const;
double partial_width_phiO_to_barFdFd(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_phiO_to_barFuFu(MRSSM2_mass_eigenstates_interface*, int, int) const;
double partial_width_phiO_to_barGlubarGlu(MRSSM2_mass_eigenstates_interface*) const;
double partial_width_phiO_to_barGlubarChi(MRSSM2_mass_eigenstates_interface*, int) const;

private:
   MRSSM2_mass_eigenstates model{};
   softsusy::QedQcd qedqcd{};
   Physical_input physical_input;
   FlexibleDecay_settings flexibledecay_settings {};
   bool run_to_decay_particle_scale {true};
   MRSSM2_decay_table decay_table{};
   FlexibleDecay_problems problems{};

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   typename Decay_amplitude_type<FieldIn, FieldOut1, FieldOut2>::type
   calculate_amplitude(
      const MRSSM2_cxx_diagrams::context_base&,
      const typename MRSSM2_cxx_diagrams::field_indices<FieldIn>::type&,
      const typename MRSSM2_cxx_diagrams::field_indices<FieldOut1>::type&,
      const typename MRSSM2_cxx_diagrams::field_indices<FieldOut2>::type&) const;

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   double amplitude_squared(MRSSM2_cxx_diagrams::context_base const& context,
                  typename MRSSM2_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
                  typename MRSSM2_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
                  typename MRSSM2_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2) const;

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   double get_partial_width(
      const MRSSM2_cxx_diagrams::context_base&,
      typename MRSSM2_cxx_diagrams::field_indices<FieldIn>::type const&,
      typename MRSSM2_cxx_diagrams::field_indices<FieldOut1>::type const&,
      typename MRSSM2_cxx_diagrams::field_indices<FieldOut2>::type const&) const;
};

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::SRdp, MRSSM2_cxx_diagrams::fields::SRum>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::SRdp >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::SRum >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::SRdp, MRSSM2_cxx_diagrams::fields::Hpm>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::SRdp >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::SRdp, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::SRdp >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::SRum, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRum>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::SRum >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRum>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::SRum, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Hpm>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::SRum >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Hpm>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::sigmaO>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::phiO>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::phiO, MRSSM2_cxx_diagrams::fields::phiO>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::Sd, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Sd>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Sd >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Sd>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::Sv, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Sv>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Sv >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Sv>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::Su, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Su>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Su >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Su>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::Se, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Se>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Se >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Se>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::hh>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::Ah>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::Rh>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Rh >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::Ah>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::Rh>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Rh >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::Rh, MRSSM2_cxx_diagrams::fields::Rh>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Rh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Rh >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::Rh, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Rh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::Hpm, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRum>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRum>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::Hpm, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Hpm>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Hpm>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRum>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRum>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Hpm>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Hpm>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::VG>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::phiO, MRSSM2_cxx_diagrams::fields::VG>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::VP>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::VP>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::Rh, MRSSM2_cxx_diagrams::fields::VP>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Rh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type, MRSSM2_cxx_diagrams::fields::VP>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::VZ>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::VZ>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::Rh, MRSSM2_cxx_diagrams::fields::VZ>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Rh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type, MRSSM2_cxx_diagrams::fields::VZ>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::SRdp, MRSSM2_cxx_diagrams::fields::VWm>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::SRdp >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRum>::type, MRSSM2_cxx_diagrams::fields::VWm>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRum>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Hpm>::type, MRSSM2_cxx_diagrams::fields::VWm>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Hpm>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::SRum, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::VWm>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::SRum >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::VWm>::type >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::Hpm, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::VWm>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::VWm>::type >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::VWm>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::VWm>::type >::type&) const;

template<>
Decay_amplitude_SVV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::VG, MRSSM2_cxx_diagrams::fields::VG>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SVV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::VP, MRSSM2_cxx_diagrams::fields::VP>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SVV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::VP, MRSSM2_cxx_diagrams::fields::VZ>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::VZ, MRSSM2_cxx_diagrams::fields::VZ>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VZ >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::VWm>::type, MRSSM2_cxx_diagrams::fields::VWm>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::VWm>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::Glu, MRSSM2_cxx_diagrams::fields::Glu>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Glu >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Glu >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type, MRSSM2_cxx_diagrams::fields::Glu>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Glu >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fv>::type, MRSSM2_cxx_diagrams::fields::Fv>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fv>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Fv >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::Chi, MRSSM2_cxx_diagrams::fields::Chi>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Chi >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type, MRSSM2_cxx_diagrams::fields::Chi>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::Cha1, MRSSM2_cxx_diagrams::fields::Cha2>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Cha1 >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Cha2 >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha1>::type, MRSSM2_cxx_diagrams::fields::Cha1>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha1>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Cha1 >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha2>::type, MRSSM2_cxx_diagrams::fields::Cha2>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha2>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Cha2 >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type, MRSSM2_cxx_diagrams::fields::Fe>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Fe >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fd>::type, MRSSM2_cxx_diagrams::fields::Fd>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fd>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fu>::type, MRSSM2_cxx_diagrams::fields::Fu>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fu>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Fu >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha1>::type, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha2>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha1>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha2>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::SRdp, MRSSM2_cxx_diagrams::fields::SRum>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::SRdp >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::SRum >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::SRdp, MRSSM2_cxx_diagrams::fields::Hpm>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::SRdp >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::SRdp, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::SRdp >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::SRum, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRum>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::SRum >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRum>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::SRum, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Hpm>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::SRum >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Hpm>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::sigmaO>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::phiO>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::phiO, MRSSM2_cxx_diagrams::fields::phiO>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::Sd, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Sd>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Sd >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Sd>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::Sv, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Sv>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Sv >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Sv>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::Su, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Su>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Su >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Su>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::Se, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Se>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Se >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Se>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::hh>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::Ah>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::Rh>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Rh >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::Ah>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::Rh>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Rh >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::Rh, MRSSM2_cxx_diagrams::fields::Rh>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Rh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Rh >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::Rh, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Rh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::Hpm, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRum>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRum>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::Hpm, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Hpm>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Hpm>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRum>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRum>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Hpm>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Hpm>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::VG>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::phiO, MRSSM2_cxx_diagrams::fields::VG>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::VP>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::VP>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::Rh, MRSSM2_cxx_diagrams::fields::VP>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Rh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type, MRSSM2_cxx_diagrams::fields::VP>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::VZ>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::VZ>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::Rh, MRSSM2_cxx_diagrams::fields::VZ>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Rh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type, MRSSM2_cxx_diagrams::fields::VZ>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::SRdp, MRSSM2_cxx_diagrams::fields::VWm>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::SRdp >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRum>::type, MRSSM2_cxx_diagrams::fields::VWm>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRum>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Hpm>::type, MRSSM2_cxx_diagrams::fields::VWm>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Hpm>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::SRum, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::VWm>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::SRum >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::VWm>::type >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::Hpm, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::VWm>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::VWm>::type >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::VWm>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::VWm>::type >::type&) const;

template<>
Decay_amplitude_SVV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::VG, MRSSM2_cxx_diagrams::fields::VG>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SVV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::VP, MRSSM2_cxx_diagrams::fields::VP>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SVV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::VP, MRSSM2_cxx_diagrams::fields::VZ>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::VZ, MRSSM2_cxx_diagrams::fields::VZ>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VZ >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::VWm>::type, MRSSM2_cxx_diagrams::fields::VWm>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::VWm>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::Glu, MRSSM2_cxx_diagrams::fields::Glu>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Glu >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Glu >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type, MRSSM2_cxx_diagrams::fields::Glu>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Glu >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fv>::type, MRSSM2_cxx_diagrams::fields::Fv>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fv>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Fv >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::Chi, MRSSM2_cxx_diagrams::fields::Chi>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Chi >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type, MRSSM2_cxx_diagrams::fields::Chi>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::Cha1, MRSSM2_cxx_diagrams::fields::Cha2>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Cha1 >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Cha2 >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha1>::type, MRSSM2_cxx_diagrams::fields::Cha1>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha1>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Cha1 >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha2>::type, MRSSM2_cxx_diagrams::fields::Cha2>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha2>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Cha2 >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type, MRSSM2_cxx_diagrams::fields::Fe>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Fe >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fd>::type, MRSSM2_cxx_diagrams::fields::Fd>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fd>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fu>::type, MRSSM2_cxx_diagrams::fields::Fu>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fu>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Fu >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha1>::type, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha2>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha1>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha2>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::SRum, MRSSM2_cxx_diagrams::fields::hh>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::SRum >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::SRum, MRSSM2_cxx_diagrams::fields::Ah>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::SRum >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::SRum, MRSSM2_cxx_diagrams::fields::Rh>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::SRum >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Rh >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::SRum, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::SRum >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::Sd, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Su>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Sd >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Su>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::Se, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Sv>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Se >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Sv>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::Hpm>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::hh, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::Hpm>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::Ah, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::Rh, MRSSM2_cxx_diagrams::fields::Hpm>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Rh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::Rh, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Rh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::Hpm, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::SRum, MRSSM2_cxx_diagrams::fields::VP>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::SRum >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::VP>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type, MRSSM2_cxx_diagrams::fields::VP>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::SRum, MRSSM2_cxx_diagrams::fields::VZ>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::SRum >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::VZ>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type, MRSSM2_cxx_diagrams::fields::VZ>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::SRdp>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::VWm>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::VWm>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::Rh, MRSSM2_cxx_diagrams::fields::VWm>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Rh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type, MRSSM2_cxx_diagrams::fields::VWm>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SVV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::VP, MRSSM2_cxx_diagrams::fields::VWm>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SVV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::VZ, MRSSM2_cxx_diagrams::fields::VWm>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VZ >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, MRSSM2_cxx_diagrams::fields::Chi, MRSSM2_cxx_diagrams::fields::Cha2>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Chi >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Cha2 >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha1>::type, MRSSM2_cxx_diagrams::fields::Chi>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha1>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type, MRSSM2_cxx_diagrams::fields::Cha2>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Cha2 >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fv>::type, MRSSM2_cxx_diagrams::fields::Fe>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fv>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Fe >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fu>::type, MRSSM2_cxx_diagrams::fields::Fd>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fu>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::Hpm, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha1>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Hpm >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Cha1>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::phiO, MRSSM2_cxx_diagrams::fields::phiO>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::phiO, MRSSM2_cxx_diagrams::fields::hh>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::phiO, MRSSM2_cxx_diagrams::fields::Ah>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::phiO, MRSSM2_cxx_diagrams::fields::Rh>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Rh >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::phiO, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::Sd, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Sd>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Sd >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Sd>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::Su, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Su>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Su >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Su>::type >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::VG>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::VG>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::Rh, MRSSM2_cxx_diagrams::fields::VG>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Rh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::sigmaO, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type, MRSSM2_cxx_diagrams::fields::VG>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::phiO, MRSSM2_cxx_diagrams::fields::VP>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::phiO, MRSSM2_cxx_diagrams::fields::VZ>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::VG, MRSSM2_cxx_diagrams::fields::VG>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SVV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::VG, MRSSM2_cxx_diagrams::fields::VP>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SVV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::VG, MRSSM2_cxx_diagrams::fields::VZ>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::Glu, MRSSM2_cxx_diagrams::fields::Glu>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Glu >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Glu >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::Glu, MRSSM2_cxx_diagrams::fields::Chi>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Glu >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::sigmaO, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type, MRSSM2_cxx_diagrams::fields::Glu>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Glu >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::sigmaO, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type, MRSSM2_cxx_diagrams::fields::Chi>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::sigmaO, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fd>::type, MRSSM2_cxx_diagrams::fields::Fd>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fd>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::sigmaO, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fu>::type, MRSSM2_cxx_diagrams::fields::Fu>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fu>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Fu >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::sigmaO, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::sigmaO, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::phiO, MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::sigmaO>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::phiO, MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::hh>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::phiO, MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::Ah>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::phiO, MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::Rh>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Rh >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::phiO, MRSSM2_cxx_diagrams::fields::sigmaO, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::phiO, MRSSM2_cxx_diagrams::fields::Sd, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Sd>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Sd >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Sd>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::phiO, MRSSM2_cxx_diagrams::fields::Su, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Su>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Su >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Su>::type >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::phiO, MRSSM2_cxx_diagrams::fields::hh, MRSSM2_cxx_diagrams::fields::VG>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::phiO, MRSSM2_cxx_diagrams::fields::Ah, MRSSM2_cxx_diagrams::fields::VG>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::phiO, MRSSM2_cxx_diagrams::fields::Rh, MRSSM2_cxx_diagrams::fields::VG>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Rh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::phiO, typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type, MRSSM2_cxx_diagrams::fields::VG>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::Rh>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::phiO, MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::VP>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::phiO, MRSSM2_cxx_diagrams::fields::sigmaO, MRSSM2_cxx_diagrams::fields::VZ>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::phiO, MRSSM2_cxx_diagrams::fields::VG, MRSSM2_cxx_diagrams::fields::VG>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SVV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::phiO, MRSSM2_cxx_diagrams::fields::VG, MRSSM2_cxx_diagrams::fields::VP>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SVV MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::phiO, MRSSM2_cxx_diagrams::fields::VG, MRSSM2_cxx_diagrams::fields::VZ>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::phiO, MRSSM2_cxx_diagrams::fields::Glu, MRSSM2_cxx_diagrams::fields::Glu>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Glu >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Glu >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::phiO, MRSSM2_cxx_diagrams::fields::Glu, MRSSM2_cxx_diagrams::fields::Chi>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Glu >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::phiO, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type, MRSSM2_cxx_diagrams::fields::Glu>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Glu >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::phiO, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type, MRSSM2_cxx_diagrams::fields::Chi>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::phiO, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fd>::type, MRSSM2_cxx_diagrams::fields::Fd>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fd>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::phiO, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fu>::type, MRSSM2_cxx_diagrams::fields::Fu>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fu>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Fu >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::phiO, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type >::type&) const;

template<>
Decay_amplitude_SFF MRSSM2_decays::calculate_amplitude<MRSSM2_cxx_diagrams::fields::phiO, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type, typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type>(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::phiO >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Glu>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Chi>::type >::type&) const;


template<typename FieldIn, typename FieldOut1, typename FieldOut2>
double
MRSSM2_decays::amplitude_squared(MRSSM2_cxx_diagrams::context_base const& context,
                  typename MRSSM2_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
                  typename MRSSM2_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
                  typename MRSSM2_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2) const
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
MRSSM2_cxx_diagrams::fields::is_singlet_v<FieldIn> &&
MRSSM2_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
MRSSM2_cxx_diagrams::fields::is_singlet_v<FieldOut2>, double>
squared_color_generator() {return 1.;}

// 1 -> 3, 3bar
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
   MRSSM2_cxx_diagrams::fields::is_singlet_v<FieldIn>
   &&
   (
   (MRSSM2_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
   MRSSM2_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>)
   ||
   (MRSSM2_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
   MRSSM2_cxx_diagrams::fields::is_triplet_v<FieldOut2>)
   ), double
>
squared_color_generator() {return 3.;}

// 1 -> 8, 8
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
MRSSM2_cxx_diagrams::fields::is_singlet_v<FieldIn> &&
MRSSM2_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
MRSSM2_cxx_diagrams::fields::is_octet_v<FieldOut2>, double>
squared_color_generator() {return 8.;}

// 3 -> 3, 1; 3bar -> 3bar, 1
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
(
MRSSM2_cxx_diagrams::fields::is_triplet_v<FieldIn> &&
((MRSSM2_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
MRSSM2_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(MRSSM2_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
MRSSM2_cxx_diagrams::fields::is_triplet_v<FieldOut2>))
)
||
(
MRSSM2_cxx_diagrams::fields::is_anti_triplet_v<FieldIn> &&
((MRSSM2_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
MRSSM2_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(MRSSM2_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
MRSSM2_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>))
)
, double>
squared_color_generator() {return 1.;}

// 3 -> 3, 8; 3bar -> 3bar, 8
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
(
MRSSM2_cxx_diagrams::fields::is_triplet_v<FieldIn> &&
((MRSSM2_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
MRSSM2_cxx_diagrams::fields::is_octet_v<FieldOut2>) ||
(MRSSM2_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
MRSSM2_cxx_diagrams::fields::is_triplet_v<FieldOut2>))
)
||
(
MRSSM2_cxx_diagrams::fields::is_anti_triplet_v<FieldIn> &&
((MRSSM2_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
MRSSM2_cxx_diagrams::fields::is_octet_v<FieldOut2>) ||
(MRSSM2_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
MRSSM2_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>))
)
, double>
squared_color_generator() {return 4.;}

// 8 -> 8, 1
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
MRSSM2_cxx_diagrams::fields::is_octet_v<FieldIn> &&
((MRSSM2_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
MRSSM2_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(MRSSM2_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
MRSSM2_cxx_diagrams::fields::is_octet_v<FieldOut2>))
, double>
squared_color_generator() {return 1.;}

// 8 -> 3, 3bar
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
   MRSSM2_cxx_diagrams::fields::is_octet_v<FieldIn>
   &&
   (
   (MRSSM2_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
   MRSSM2_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>)
   ||
   (MRSSM2_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
   MRSSM2_cxx_diagrams::fields::is_triplet_v<FieldOut2>)
   ), double
>
squared_color_generator() {return 1./2.;}

// 8 -> 8, 8 with identical particles in the final state
// because of symmetry of the final state it must be proportional to d^2
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
MRSSM2_cxx_diagrams::fields::is_octet_v<FieldIn> &&
MRSSM2_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
MRSSM2_cxx_diagrams::fields::is_octet_v<FieldOut2> &&
std::is_same<FieldOut1, FieldOut2>::value
, double>
// color:   d^2 = (2 (4 - 5 Nc^2 + Nc^4) TR)/Nc = 40/3
// average: 1/8
squared_color_generator() {return 40/24.;}

// 8 -> 8, 8 with differnt particles in the final state
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
MRSSM2_cxx_diagrams::fields::is_octet_v<FieldIn> &&
MRSSM2_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
MRSSM2_cxx_diagrams::fields::is_octet_v<FieldOut2> &&
!std::is_same<FieldOut1, FieldOut2>::value
, double>
// color:   f^2 = 2 Nc (-1 + Nc^2) TR = 24
// average: 1/8
squared_color_generator() {return 3.;}

// generic decay of FieldIn -> FieldOut1 FieldOut2
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
double MRSSM2_decays::get_partial_width(
   const MRSSM2_cxx_diagrams::context_base& context,
   typename MRSSM2_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
   typename MRSSM2_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
   typename MRSSM2_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2
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
double MRSSM2_decays::get_partial_width<MRSSM2_cxx_diagrams::fields::hh,MRSSM2_cxx_diagrams::fields::VZ,MRSSM2_cxx_diagrams::fields::VZ >(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VZ >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VZ >::type&) const;
template <>
double MRSSM2_decays::get_partial_width<MRSSM2_cxx_diagrams::fields::hh,typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::VWm>::type,MRSSM2_cxx_diagrams::fields::VWm >(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::conj<MRSSM2_cxx_diagrams::fields::VWm>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VWm >::type&) const;
template <>
double MRSSM2_decays::get_partial_width<MRSSM2_cxx_diagrams::fields::hh,MRSSM2_cxx_diagrams::fields::VG,MRSSM2_cxx_diagrams::fields::VG >(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&) const;
template <>
double MRSSM2_decays::get_partial_width<MRSSM2_cxx_diagrams::fields::hh,MRSSM2_cxx_diagrams::fields::VP,MRSSM2_cxx_diagrams::fields::VP >(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&) const;
template <>
double MRSSM2_decays::get_partial_width<MRSSM2_cxx_diagrams::fields::hh,MRSSM2_cxx_diagrams::fields::VP,MRSSM2_cxx_diagrams::fields::VZ >(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VZ >::type&) const;
template <>
double MRSSM2_decays::get_partial_width<MRSSM2_cxx_diagrams::fields::hh,typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fu>::type,MRSSM2_cxx_diagrams::fields::Fu >(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fu>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Fu >::type&) const;
template <>
double MRSSM2_decays::get_partial_width<MRSSM2_cxx_diagrams::fields::hh,typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fd>::type,MRSSM2_cxx_diagrams::fields::Fd >(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fd>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Fd >::type&) const;
template <>
double MRSSM2_decays::get_partial_width<MRSSM2_cxx_diagrams::fields::hh,typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type,MRSSM2_cxx_diagrams::fields::Fe >(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::hh >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Fe >::type&) const;
template <>
double MRSSM2_decays::get_partial_width<MRSSM2_cxx_diagrams::fields::Ah,typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fd>::type,MRSSM2_cxx_diagrams::fields::Fd >(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fd>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Fd >::type&) const;
template <>
double MRSSM2_decays::get_partial_width<MRSSM2_cxx_diagrams::fields::Ah,typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fu>::type,MRSSM2_cxx_diagrams::fields::Fu >(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fu>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Fu >::type&) const;
template <>
double MRSSM2_decays::get_partial_width<MRSSM2_cxx_diagrams::fields::Ah,MRSSM2_cxx_diagrams::fields::VG,MRSSM2_cxx_diagrams::fields::VG >(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VG >::type&) const;
template <>
double MRSSM2_decays::get_partial_width<MRSSM2_cxx_diagrams::fields::Ah,MRSSM2_cxx_diagrams::fields::VP,MRSSM2_cxx_diagrams::fields::VP >(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&) const;
template <>
double MRSSM2_decays::get_partial_width<MRSSM2_cxx_diagrams::fields::Ah,MRSSM2_cxx_diagrams::fields::VP,MRSSM2_cxx_diagrams::fields::VZ >(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VP >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::VZ >::type&) const;
template <>
double MRSSM2_decays::get_partial_width<MRSSM2_cxx_diagrams::fields::Ah,typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type,MRSSM2_cxx_diagrams::fields::Fe >(const MRSSM2_cxx_diagrams::context_base&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Ah >::type&, const typename MRSSM2_cxx_diagrams::field_indices<typename MRSSM2_cxx_diagrams::fields::bar<MRSSM2_cxx_diagrams::fields::Fe>::type >::type&, const typename MRSSM2_cxx_diagrams::field_indices<MRSSM2_cxx_diagrams::fields::Fe >::type&) const;

} // namespace flexiblesusy

#endif
