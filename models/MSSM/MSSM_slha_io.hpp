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


#ifndef MSSM_SLHA_IO_H
#define MSSM_SLHA_IO_H

#include "problems.hpp"
#include "slha_io.hpp"
#include "for_each.hpp"

#include <Eigen/Core>
#include <string>
#include <tuple>

namespace flexiblesusy {

struct MSSM_input_parameters;
class MSSM_mass_eigenstates;
struct MSSM_observables;
struct MSSM_physical;
class MSSM_slha;
class Spectrum_generator_problems;
class Spectrum_generator_settings;

namespace standard_model {
class Standard_model;
struct Standard_model_physical;
} // namespace standard_model

struct MSSM_scales {
   double HighScale{0.}, SUSYScale{0.}, LowScale{0.};
   double pole_mass_scale{0.};
};

class MSSM_slha_io {
public:
   MSSM_slha_io();

   void clear();

   void fill(softsusy::QedQcd& qedqcd) const { slha_io.fill(qedqcd); }
   void fill(MSSM_input_parameters&) const;
   void fill(MSSM_mass_eigenstates&) const;
   void fill(MSSM_slha&) const;
   void fill(Physical_input&) const;
   void fill(Spectrum_generator_settings&) const;
   double get_parameter_output_scale() const;
   const SLHA_io& get_slha_io() const { return slha_io; }
   void read_from_file(const std::string&);
   void read_from_source(const std::string&);
   void read_from_stream(std::istream&);
   void set_block(const std::string& str, SLHA_io::Position position = SLHA_io::back) { slha_io.set_block(str, position); }
   void set_blocks(const std::vector<std::string>& vec, SLHA_io::Position position = SLHA_io::back) { slha_io.set_blocks(vec, position); }
   void set_extra(const MSSM_slha&, const MSSM_scales&, const MSSM_observables&, const flexiblesusy::Spectrum_generator_settings&);
   void set_input(const MSSM_input_parameters&);
   void set_modsel(const SLHA_io::Modsel&);
   void set_physical_input(const Physical_input&);
   void set_settings(const Spectrum_generator_settings&);
   void set_sminputs(const softsusy::QedQcd&);
   template <class... Ts> void set_spectrum(const std::tuple<Ts...>&);
   void set_spectrum(const MSSM_slha&);
   void set_spectrum(const standard_model::Standard_model&);
   void set_spinfo(const Spectrum_generator_problems&);
   void set_spinfo(const Problems&);
   void set_spinfo(const std::vector<std::string>&, const std::vector<std::string>&);
   void set_print_imaginary_parts_of_majorana_mixings(bool);
   void write_to(const std::string&) const;
   void write_to_file(const std::string&) const;
   void write_to_stream() const;
   void write_to_stream(std::ostream&) const;

   static void fill_minpar_tuple(MSSM_input_parameters&, int, double);
   static void fill_extpar_tuple(MSSM_input_parameters&, int, double);
   static void fill_imminpar_tuple(MSSM_input_parameters&, int, double);
   static void fill_imextpar_tuple(MSSM_input_parameters&, int, double);

   template <class... Ts>
   void fill(const std::tuple<Ts...>&, const softsusy::QedQcd&, const MSSM_scales&, const MSSM_observables&, const Spectrum_generator_settings&);

private:
   SLHA_io slha_io; ///< SLHA io class
   bool print_imaginary_parts_of_majorana_mixings;

   void set_extpar(const MSSM_input_parameters&);
   void set_imminpar(const MSSM_input_parameters&);
   void set_imextpar(const MSSM_input_parameters&);
   void set_minpar(const MSSM_input_parameters&);
   void set_mass(const MSSM_physical&, bool);
   void set_mass(const standard_model::Standard_model_physical&);
   void set_mixing_matrices(const MSSM_physical&, bool);
   void set_mixing_matrices(const standard_model::Standard_model_physical&);
   void set_model_parameters(const MSSM_slha&);
   void set_model_parameters(const standard_model::Standard_model&);
   void set_ckm(const Eigen::Matrix<std::complex<double>,3,3>&, double);
   void set_pmns(const Eigen::Matrix<std::complex<double>,3,3>&, double);
   double read_scale() const;
   void fill_drbar_parameters(MSSM_mass_eigenstates&) const;
   void fill_physical(MSSM_physical&) const;
};

/**
 * Convenience wrapper function to fill all information at once
 */
template <class... Ts>
void MSSM_slha_io::fill(
   const std::tuple<Ts...>& models,
   const softsusy::QedQcd& qedqcd,
   const MSSM_scales& scales,
   const MSSM_observables& observables,
   const Spectrum_generator_settings& spectrum_generator_settings)
{
   const auto& model = std::get<0>(models);
   const auto& problems = model.get_problems();

   set_spinfo(problems);
   set_sminputs(qedqcd);
   set_input(model.get_input());
   if (!problems.have_problem()) {
      set_spectrum(models);
      set_extra(model, scales, observables, spectrum_generator_settings);
   }
}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices of
 * all given models in the SLHA object.
 *
 * @param models model classes
 */
template <class... Ts>
void MSSM_slha_io::set_spectrum(const std::tuple<Ts...>& models)
{
   for_each_in_tuple(models,
                     [this](const auto& model) { this->set_spectrum(model); });
}

} // namespace flexiblesusy

#endif
