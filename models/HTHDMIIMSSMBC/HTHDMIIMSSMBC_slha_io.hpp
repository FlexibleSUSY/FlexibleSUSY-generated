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

// File generated at Mon 5 Mar 2018 17:33:14

#ifndef HTHDMIIMSSMBC_SLHA_IO_H
#define HTHDMIIMSSMBC_SLHA_IO_H

#include "HTHDMIIMSSMBC_mass_eigenstates.hpp"
#include "HTHDMIIMSSMBC_model_slha.hpp"
#include "HTHDMIIMSSMBC_info.hpp"
#include "HTHDMIIMSSMBC_observables.hpp"
#include "HTHDMIIMSSMBC_physical.hpp"
#include "problems.hpp"
#include "spectrum_generator_problems.hpp"
#include "standard_model_two_scale_model.hpp"
#include "slha_io.hpp"
#include "ckm.hpp"
#include "ew_input.hpp"
#include "lowe.h"

#include <Eigen/Core>
#include <string>
#include <tuple>
#include <utility>

#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/adapted/std_tuple.hpp>

#define Pole(p) physical.p
#define PHYSICAL(p) model.get_physical().p
#define PHYSICAL_SLHA(p) model.get_physical_slha().p
#define LOCALPHYSICAL(p) physical.p
#define MODEL model
#define MODELPARAMETER(p) model.get_##p()
#define EXTRAPARAMETER(p) model.get_##p()
#define OBSERVABLES observables
#define LowEnergyConstant(p) Electroweak_constants::p
#define SCALES(p) scales.p

namespace flexiblesusy {

struct HTHDMIIMSSMBC_input_parameters;
class Spectrum_generator_settings;

template <class T>
class HTHDMIIMSSMBC;

struct HTHDMIIMSSMBC_scales {
   double HighScale{0.}, SUSYScale{0.}, LowScale{0.};
   double pole_mass_scale{0.};
};

class HTHDMIIMSSMBC_slha_io {
public:
   HTHDMIIMSSMBC_slha_io();

   void clear();

   void fill(softsusy::QedQcd& qedqcd) const { slha_io.fill(qedqcd); }
   void fill(HTHDMIIMSSMBC_input_parameters&) const;
   void fill(HTHDMIIMSSMBC_mass_eigenstates&) const;
   template <class Model> void fill(HTHDMIIMSSMBC_slha<Model>&) const;
   void fill(Physical_input&) const;
   void fill(Spectrum_generator_settings&) const;
   double get_parameter_output_scale() const;
   const SLHA_io& get_slha_io() const { return slha_io; }
   void read_from_file(const std::string&);
   void read_from_source(const std::string&);
   void read_from_stream(std::istream&);
   void set_block(const std::string& str, SLHA_io::Position position = SLHA_io::back) { slha_io.set_block(str, position); }
   void set_blocks(const std::vector<std::string>& vec, SLHA_io::Position position = SLHA_io::back) { slha_io.set_blocks(vec, position); }
   template <class Model> void set_extra(const HTHDMIIMSSMBC_slha<Model>&, const HTHDMIIMSSMBC_scales&, const HTHDMIIMSSMBC_observables&);
   void set_input(const HTHDMIIMSSMBC_input_parameters&);
   void set_modsel(const SLHA_io::Modsel&);
   void set_physical_input(const Physical_input&);
   void set_settings(const Spectrum_generator_settings&);
   void set_sminputs(const softsusy::QedQcd&);
   template <class... Ts> void set_spectrum(const std::tuple<Ts...>&);
   template <class Model> void set_spectrum(const HTHDMIIMSSMBC_slha<Model>&);
   template <class T> void set_spectrum(const HTHDMIIMSSMBC<T>&);
   void set_spectrum(const standard_model::Standard_model&);
   void set_spinfo(const Spectrum_generator_problems&);
   void set_spinfo(const Problems&);
   void set_spinfo(const std::vector<std::string>&, const std::vector<std::string>&);
   void set_print_imaginary_parts_of_majorana_mixings(bool);
   void write_to(const std::string&) const;
   void write_to_file(const std::string& file_name) const { slha_io.write_to_file(file_name); }
   void write_to_stream(std::ostream& ostr = std::cout) const { slha_io.write_to_stream(ostr); }

   static void fill_minpar_tuple(HTHDMIIMSSMBC_input_parameters&, int, double);
   static void fill_extpar_tuple(HTHDMIIMSSMBC_input_parameters&, int, double);
   static void fill_imminpar_tuple(HTHDMIIMSSMBC_input_parameters&, int, double);
   static void fill_imextpar_tuple(HTHDMIIMSSMBC_input_parameters&, int, double);

   template <class Model>
   static void fill_slhaea(SLHAea::Coll&, const HTHDMIIMSSMBC_slha<Model>&, const softsusy::QedQcd&, const HTHDMIIMSSMBC_scales&, const HTHDMIIMSSMBC_observables&);

   template <class Model>
   static SLHAea::Coll fill_slhaea(const HTHDMIIMSSMBC_slha<Model>&, const softsusy::QedQcd&, const HTHDMIIMSSMBC_scales&, const HTHDMIIMSSMBC_observables&);

private:
   SLHA_io slha_io; ///< SLHA io class
   bool print_imaginary_parts_of_majorana_mixings;

   void set_extpar(const HTHDMIIMSSMBC_input_parameters&);
   void set_imminpar(const HTHDMIIMSSMBC_input_parameters&);
   void set_imextpar(const HTHDMIIMSSMBC_input_parameters&);
   void set_minpar(const HTHDMIIMSSMBC_input_parameters&);
   void set_mass(const HTHDMIIMSSMBC_physical&, bool);
   void set_mass(const standard_model::Standard_model_physical&);
   void set_mixing_matrices(const HTHDMIIMSSMBC_physical&, bool);
   void set_mixing_matrices(const standard_model::Standard_model_physical&);
   template <class Model> void set_model_parameters(const HTHDMIIMSSMBC_slha<Model>&);
   void set_model_parameters(const standard_model::Standard_model&);
   void set_ckm(const Eigen::Matrix<std::complex<double>,3,3>&, double);
   void set_pmns(const Eigen::Matrix<std::complex<double>,3,3>&, double);
   double read_scale() const;
   void fill_drbar_parameters(HTHDMIIMSSMBC_mass_eigenstates&) const;
   void fill_physical(HTHDMIIMSSMBC_physical&) const;
};

/**
 * Reads DR-bar parameters, pole masses and mixing matrices from a
 * SLHA output file.
 */
template <class Model>
void HTHDMIIMSSMBC_slha_io::fill(HTHDMIIMSSMBC_slha<Model>& model) const
{
   fill(static_cast<HTHDMIIMSSMBC_mass_eigenstates&>(model));
   fill_physical(model.get_physical_slha());
}

template <class Model>
void HTHDMIIMSSMBC_slha_io::fill_slhaea(
   SLHAea::Coll& slhaea, const HTHDMIIMSSMBC_slha<Model>& model,
   const softsusy::QedQcd& qedqcd, const HTHDMIIMSSMBC_scales& scales,
   const HTHDMIIMSSMBC_observables& observables)
{
   HTHDMIIMSSMBC_slha_io slha_io;
   const HTHDMIIMSSMBC_input_parameters& input = model.get_input();
   const auto& problems = model.get_problems();
   const bool error = problems.have_problem();

   slha_io.set_spinfo(problems);
   slha_io.set_sminputs(qedqcd);
   slha_io.set_input(input);
   if (!error) {
      slha_io.set_spectrum(model);
      slha_io.set_extra(model, scales, observables);
   }

   slhaea = slha_io.get_slha_io().get_data();
}

template <class Model>
SLHAea::Coll HTHDMIIMSSMBC_slha_io::fill_slhaea(
   const HTHDMIIMSSMBC_slha<Model>& model, const softsusy::QedQcd& qedqcd,
   const HTHDMIIMSSMBC_scales& scales, const HTHDMIIMSSMBC_observables& observables)
{
   SLHAea::Coll slhaea;
   HTHDMIIMSSMBC_slha_io::fill_slhaea(slhaea, model, qedqcd, scales, observables);

   return slhaea;
}

/**
 * Stores the model (DR-bar) parameters in the SLHA object.
 *
 * @param model model class
 */
template <class Model>
void HTHDMIIMSSMBC_slha_io::set_model_parameters(const HTHDMIIMSSMBC_slha<Model>& model)
{
   {
      std::ostringstream block;
      block << "Block gauge Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(g1) * 0.7745966692414834), "g1 * 0.7745966692414834")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(g2)), "g2")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(g3)), "g3")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("Yu", ToMatrix(MODELPARAMETER(Yu_slha)), "Yu", model.get_scale());
   slha_io.set_block("Yd", ToMatrix(MODELPARAMETER(Yd_slha)), "Yd", model.get_scale());
   slha_io.set_block("Ye", ToMatrix(MODELPARAMETER(Ye_slha)), "Ye", model.get_scale());
   {
      std::ostringstream block;
      block << "Block HMIX Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(31, (MODELPARAMETER(Lambda1)), "Lambda1")
            << FORMAT_ELEMENT(32, (MODELPARAMETER(Lambda2)), "Lambda2")
            << FORMAT_ELEMENT(33, (MODELPARAMETER(Lambda3)), "Lambda3")
            << FORMAT_ELEMENT(34, (MODELPARAMETER(Lambda4)), "Lambda4")
            << FORMAT_ELEMENT(35, (MODELPARAMETER(Lambda5)), "Lambda5")
            << FORMAT_ELEMENT(36, (MODELPARAMETER(Lambda6)), "Lambda6")
            << FORMAT_ELEMENT(37, (MODELPARAMETER(Lambda7)), "Lambda7")
            << FORMAT_ELEMENT(20, (MODELPARAMETER(M112)), "M112")
            << FORMAT_ELEMENT(21, (MODELPARAMETER(M222)), "M222")
            << FORMAT_ELEMENT(22, (MODELPARAMETER(M122)), "M122")
            << FORMAT_ELEMENT(102, (MODELPARAMETER(v1)), "v1")
            << FORMAT_ELEMENT(103, (MODELPARAMETER(v2)), "v2")
            << FORMAT_ELEMENT(1, (MODELPARAMETER(Mu)), "Mu")
      ;
      slha_io.set_block(block);
   }


}

/**
 * Writes extra SLHA blocks
 *
 * @param model model class
 * @param scales struct of boundary condition scales
 * @param observables struct of observables
 */
template <class Model>
void HTHDMIIMSSMBC_slha_io::set_extra(
   const HTHDMIIMSSMBC_slha<Model>& model, const HTHDMIIMSSMBC_scales& scales,
   const HTHDMIIMSSMBC_observables& observables)
{
   const HTHDMIIMSSMBC_physical physical(model.get_physical_slha());

   {
      std::ostringstream block;
      block << "Block FlexibleSUSYLowEnergy Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (OBSERVABLES.a_muon), "Delta(g-2)_muon/2 FlexibleSUSY")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block EFFHIGGSCOUPLINGS" << '\n'
            << FORMAT_RANK_THREE_TENSOR(25, 22, 22, (Abs(OBSERVABLES.eff_cp_higgs_photon_photon(0))), "Abs(effective H-Photon-Photon coupling)")
            << FORMAT_RANK_THREE_TENSOR(35, 22, 22, (Abs(OBSERVABLES.eff_cp_higgs_photon_photon(1))), "Abs(effective H-Photon-Photon coupling)")
            << FORMAT_RANK_THREE_TENSOR(25, 21, 21, (Abs(OBSERVABLES.eff_cp_higgs_gluon_gluon(0))), "Abs(effective H-Gluon-Gluon coupling)")
            << FORMAT_RANK_THREE_TENSOR(35, 21, 21, (Abs(OBSERVABLES.eff_cp_higgs_gluon_gluon(1))), "Abs(effective H-Gluon-Gluon coupling)")
            << FORMAT_RANK_THREE_TENSOR(36, 22, 22, (Abs(OBSERVABLES.eff_cp_pseudoscalar_photon_photon)), "Abs(effective A-Photon-Photon coupling)")
            << FORMAT_RANK_THREE_TENSOR(36, 21, 21, (Abs(OBSERVABLES.eff_cp_pseudoscalar_gluon_gluon)), "Abs(effective A-Gluon-Gluon coupling)")
      ;
      slha_io.set_block(block);
   }

}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices of
 * all given models in the SLHA object.
 *
 * @todo Use generic lambda instead of Set_spectrum in C++14
 *
 * @param models model classes
 */
template <class... Ts>
void HTHDMIIMSSMBC_slha_io::set_spectrum(const std::tuple<Ts...>& models)
{
   Set_spectrum<HTHDMIIMSSMBC_slha_io> ss(this);
   boost::fusion::for_each(models, ss);
}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in BPMZ convention
 */
template <class T>
void HTHDMIIMSSMBC_slha_io::set_spectrum(const HTHDMIIMSSMBC<T>& model)
{
   set_spectrum(HTHDMIIMSSMBC_slha<HTHDMIIMSSMBC<T> >(model));
}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in SLHA convention
 */
template <class Model>
void HTHDMIIMSSMBC_slha_io::set_spectrum(const HTHDMIIMSSMBC_slha<Model>& model)
{
   const HTHDMIIMSSMBC_physical physical(model.get_physical_slha());
   const bool write_sm_masses = model.do_calculate_sm_pole_masses();

   set_model_parameters(model);
   set_mass(physical, write_sm_masses);
   set_mixing_matrices(physical, write_sm_masses);

   if (slha_io.get_modsel().quark_flavour_violated)
      set_ckm(model.get_ckm_matrix(), model.get_scale());

   if (slha_io.get_modsel().lepton_flavour_violated)
      set_pmns(model.get_pmns_matrix(), model.get_scale());
}

} // namespace flexiblesusy

#undef Pole
#undef PHYSICAL
#undef PHYSICAL_SLHA
#undef LOCALPHYSICAL
#undef MODEL
#undef MODELPARAMETER
#undef EXTRAPARAMETER
#undef OBSERVABLES
#undef LowEnergyConstant
#undef SCALES

#endif
