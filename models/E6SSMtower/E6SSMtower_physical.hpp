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

// File generated at Wed 12 Apr 2017 12:36:08

#ifndef E6SSMtower_PHYSICAL_H
#define E6SSMtower_PHYSICAL_H

#include <Eigen/Core>

#include <iosfwd>

namespace flexiblesusy {

struct E6SSMtower_physical {
   E6SSMtower_physical();
   void clear();
   void convert_to_hk();   ///< converts pole masses to HK convention
   void convert_to_slha(); ///< converts pole masses to SLHA convention
   Eigen::ArrayXd get() const; ///< returns array with all masses and mixings
   void set(const Eigen::ArrayXd&); ///< set all masses and mixings
   Eigen::ArrayXd get_masses() const; ///< returns array with all masses
   void set_masses(const Eigen::ArrayXd&); ///< set all masses
   void print(std::ostream&) const;

   double MVG;
   double MGlu;
   Eigen::Array<double,3,1> MFv;
   double MChaP;
   Eigen::Array<double,6,1> MSd;
   Eigen::Array<double,3,1> MSv;
   Eigen::Array<double,6,1> MSu;
   Eigen::Array<double,6,1> MSe;
   Eigen::Array<double,6,1> MSDX;
   Eigen::Array<double,3,1> Mhh;
   Eigen::Array<double,3,1> MAh;
   Eigen::Array<double,2,1> MHpm;
   Eigen::Array<double,6,1> MChi;
   Eigen::Array<double,2,1> MCha;
   Eigen::Array<double,3,1> MFe;
   Eigen::Array<double,3,1> MFd;
   Eigen::Array<double,3,1> MFu;
   Eigen::Array<double,3,1> MFDX;
   Eigen::Array<double,4,1> MSHI0;
   Eigen::Array<double,4,1> MSHIp;
   Eigen::Array<double,2,1> MChaI;
   Eigen::Array<double,4,1> MChiI;
   Eigen::Array<double,2,1> MSSI0;
   Eigen::Array<double,2,1> MFSI;
   Eigen::Array<double,2,1> MSHp0;
   Eigen::Array<double,2,1> MSHpp;
   Eigen::Array<double,2,1> MChiP;
   double MVWm;
   double MVP;
   double MVZ;
   double MVZp;

   Eigen::Matrix<double,6,6> ZD;
   Eigen::Matrix<double,3,3> ZV;
   Eigen::Matrix<double,6,6> ZU;
   Eigen::Matrix<double,6,6> ZE;
   Eigen::Matrix<double,6,6> ZDX;
   Eigen::Matrix<double,3,3> ZH;
   Eigen::Matrix<double,3,3> ZA;
   Eigen::Matrix<double,2,2> ZP;
   Eigen::Matrix<std::complex<double>,6,6> ZN;
   Eigen::Matrix<std::complex<double>,2,2> UM;
   Eigen::Matrix<std::complex<double>,2,2> UP;
   Eigen::Matrix<std::complex<double>,3,3> ZEL;
   Eigen::Matrix<std::complex<double>,3,3> ZER;
   Eigen::Matrix<std::complex<double>,3,3> ZDL;
   Eigen::Matrix<std::complex<double>,3,3> ZDR;
   Eigen::Matrix<std::complex<double>,3,3> ZUL;
   Eigen::Matrix<std::complex<double>,3,3> ZUR;
   Eigen::Matrix<std::complex<double>,3,3> ZDXL;
   Eigen::Matrix<std::complex<double>,3,3> ZDXR;
   Eigen::Matrix<double,4,4> UHI0;
   Eigen::Matrix<double,4,4> UHIp;
   Eigen::Matrix<std::complex<double>,2,2> ZMI;
   Eigen::Matrix<std::complex<double>,2,2> ZPI;
   Eigen::Matrix<std::complex<double>,4,4> ZNI;
   Eigen::Matrix<double,2,2> ZSSI;
   Eigen::Matrix<std::complex<double>,2,2> ZFSI;
   Eigen::Matrix<double,2,2> UHp0;
   Eigen::Matrix<double,2,2> UHpp;
   Eigen::Matrix<std::complex<double>,2,2> ZNp;
   Eigen::Matrix<double,3,3> ZZ;

};

std::ostream& operator<<(std::ostream&, const E6SSMtower_physical&);

} // namespace flexiblesusy

#endif
