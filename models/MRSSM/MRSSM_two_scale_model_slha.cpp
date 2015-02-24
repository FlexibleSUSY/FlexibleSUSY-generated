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

// File generated at Tue 24 Feb 2015 17:34:04

/**
 * @file MRSSM_two_scale_model_slha.cpp
 * @brief MRSSM model class wrapper for SLHA conversion
 */

#include "MRSSM_two_scale_model_slha.hpp"
#include "MRSSM_slha_io.hpp"

namespace flexiblesusy {

#define CLASSNAME MRSSM_slha<Two_scale>

CLASSNAME::MRSSM_slha(const MRSSM_input_parameters& input_)
   : MRSSM<Two_scale>(input_)
   , physical_slha()
{
}

/**
 * Copy constructor.  Copies from base class (two-scale model class in
 * BPMZ convention) and converts parameters to SLHA.
 *
 * @param model_ model class in BPMZ convention
 */
CLASSNAME::MRSSM_slha(const MRSSM<Two_scale>& model_)
   : MRSSM<Two_scale>(model_)
{
   convert_to_slha();
}

CLASSNAME::~MRSSM_slha()
{
}

void CLASSNAME::clear()
{
   MRSSM<Two_scale>::clear();
   physical_slha.clear();
}

void CLASSNAME::calculate_spectrum()
{
   MRSSM<Two_scale>::calculate_spectrum();
   convert_to_slha();
}

void CLASSNAME::convert_to_slha()
{
   physical_slha = get_physical();
   MRSSM_slha_io::convert_to_slha_convention(physical_slha);
}

const MRSSM_physical& CLASSNAME::get_physical_slha() const
{
   return physical_slha;
}

MRSSM_physical& CLASSNAME::get_physical_slha()
{
   return physical_slha;
}

void CLASSNAME::print(std::ostream& ostr) const
{
   MRSSM<Two_scale>::print(ostr);

   ostr << "----------------------------------------\n"
           "SLHA convention:\n"
           "----------------------------------------\n";
   physical_slha.print(ostr);
}

} // namespace flexiblesusy
