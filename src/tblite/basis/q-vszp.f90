! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

!> @file tblite/basis/qvszp.f90
!>
!> Implements the charge dependent single-zeta q-vSZP basis set based on: 
!> Marcel Müller, Andreas Hansen, and Stefan Grimme, 
!> "An atom-in-molecule adaptive polarized valence single-ζ 
!> atomic orbital basis for electronic structure calculations"
!> J. Chem. Phys. 159, 164108 (2023). DOI: 10.1063/5.0172373
module tblite_basis_qvszp
   
   use iso_fortran_env, only: output_unit


   use mctc_env, only : wp
   use mctc_io_constants, only : pi
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type, integral_cutoff
   use tblite_basis_type, only : cgto_type
   use tblite_integral_overlap, only : overlap_cgto, overlap_grad_decontracted_cgto, msao
   use tblite_ncoord, only : ncoord_type, new_ncoord
   implicit none
   private

   public :: new_basis, add_qvszp_cgtos, scale_cgto, maxg

   !> Maximal number of primitives per CGTO
   integer, parameter :: maxg = 9
   !> Maximal number of elements
   integer, parameter :: max_elem = 103
   !> Maximal number of shells
   integer, parameter :: max_shell = 4

   !> Two over pi
   real(wp), parameter :: top = 2.0_wp / pi
   !> Double factorial, see OEIS A001147
   real(wp), parameter :: dfactorial(8) = &
      & [1.0_wp,1.0_wp,3.0_wp,15.0_wp,105.0_wp,945.0_wp,10395.0_wp,135135.0_wp]

   !> Contracted Gaussian type basis function of qvszp type
   type, public, extends(cgto_type) :: qvszp_cgto_type
      !> Contraction coefficients of the primitive Gaussian functions,
      real(wp) :: coeff0(maxg) = 0.0_wp
      !> Primitive-specific strength of the environment specific scaling
      real(wp) :: coeff1(maxg) = 0.0_wp
      !> Atom-specific parameter scaling the charge-dependence of the coefficient
      real(wp) :: k1 = 0.0_wp
      !> Atom-specific parameter scaling the CN-dependence of the coefficient
      real(wp) :: k2 = 0.0_wp
      !> Global parameter scaling the mixed charge- and CN-dependence of the coefficient
      real(wp) :: k3 = 0.0_wp
   contains
      !> Scales the coefficient of the cgto 
      procedure :: scale_cgto
   end type qvszp_cgto_type

   !> Collection of information regarding the basis set of a system
   type, public, extends(basis_type) :: qvszp_basis_type 
      !> Coordination number for modifying the self-energies
      class(ncoord_type), allocatable :: ncoord
   end type qvszp_basis_type


   integer, parameter :: nshell(max_elem) = [ &
   & 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3, & ! -20
   & 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, & ! -40
   & 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 4, 4, 4, & ! -60
   & 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 2, & ! -80
   & 3, 3, 3, 3, 3, 3, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, & ! -100
   & 4, 4, 4]

   integer, parameter :: n_prim(max_shell, max_elem) = reshape([&
   & 8, 3, 0, 0, 8, 0, 0, 0, 5, 5, 0, 0, 6, 4, 0, 0, 6, 5, 0, 0, 6, 6, 0, 0, & ! -6
   & 6, 6, 0, 0, 6, 6, 0, 0, 6, 6, 0, 0, 6, 6, 0, 0, 4, 4, 0, 0, 4, 3, 2, 0, & ! -12
   & 5, 4, 2, 0, 5, 4, 2, 0, 5, 5, 2, 0, 5, 5, 2, 0, 5, 5, 2, 0, 5, 5, 2, 0, & ! -18
   & 4, 3, 0, 0, 4, 3, 3, 0, 5, 2, 6, 0, 5, 3, 6, 0, 5, 3, 6, 0, 5, 3, 6, 0, & ! -24
   & 5, 3, 6, 0, 5, 3, 6, 0, 5, 3, 6, 0, 5, 3, 6, 0, 5, 3, 6, 0, 5, 3, 0, 0, & ! -30
   & 5, 4, 2, 0, 6, 5, 2, 0, 6, 5, 2, 0, 6, 5, 2, 0, 6, 5, 2, 0, 6, 5, 2, 0, & ! -36
   & 4, 3, 0, 0, 4, 3, 3, 0, 5, 2, 6, 0, 5, 3, 6, 0, 5, 3, 6, 0, 5, 3, 6, 0, & ! -42
   & 5, 2, 6, 0, 5, 3, 6, 0, 5, 3, 6, 0, 5, 3, 6, 0, 5, 3, 6, 0, 5, 3, 0, 0, & ! -48
   & 5, 4, 2, 0, 6, 5, 2, 0, 6, 5, 2, 0, 6, 5, 2, 0, 6, 5, 2, 0, 6, 5, 2, 0, & ! -54
   & 4, 3, 0, 0, 4, 3, 3, 0, 4, 2, 5, 0, 7, 4, 7, 7, 7, 4, 6, 7, 7, 4, 6, 7, & ! -60
   & 7, 4, 7, 7, 7, 4, 6, 7, 7, 4, 6, 7, 7, 4, 6, 7, 7, 4, 7, 7, 7, 4, 6, 7, & ! -66
   & 7, 4, 6, 7, 7, 4, 6, 7, 7, 4, 6, 7, 7, 4, 6, 7, 7, 4, 6, 7, 4, 3, 5, 0, & ! -72
   & 5, 3, 5, 0, 5, 3, 5, 0, 5, 3, 5, 0, 5, 3, 5, 0, 5, 3, 5, 0, 5, 3, 5, 0, & ! -78
   & 5, 3, 5, 0, 5, 3, 0, 0, 5, 4, 2, 0, 5, 5, 2, 0, 6, 5, 2, 0, 6, 5, 2, 0, & ! -84
   & 6, 5, 2, 0, 6, 5, 2, 0, 9, 9, 0, 0, 9, 9, 4, 0, 7, 4, 6, 7, 7, 4, 6, 7, & ! -90
   & 7, 4, 6, 7, 7, 4, 6, 7, 7, 4, 6, 7, 7, 4, 6, 7, 7, 4, 6, 7, 7, 4, 6, 7, & ! -96
   & 7, 4, 6, 7, 7, 4, 6, 7, 7, 4, 6, 7, 7, 4, 6, 7, 7, 4, 6, 7, 7, 4, 6, 7, & ! -102
   & 7, 4, 6, 7], shape(n_prim))

   integer, parameter :: ang_shell(max_shell, max_elem) = reshape([&
   & 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, & ! -6
   & 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 2, 0, & ! -12
   & 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, & ! -18
   & 0, 1, 0, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, & ! -24
   & 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 0, 0, & ! -30
   & 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, & ! -36
   & 0, 1, 0, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, & ! -42
   & 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 0, 0, & ! -48
   & 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, & ! -54
   & 0, 1, 0, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, & ! -60
   & 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, & ! -66
   & 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 0, & ! -72
   & 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, & ! -78
   & 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, & ! -84
   & 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 0, 0, 1, 2, 3, 0, 1, 2, 3, & ! -90
   & 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, & ! -96
   & 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, & ! -102
   & 0, 1, 2, 3], shape(ang_shell))

   !> Parameter for quadratic charge dependence 
   real(wp), parameter :: p_k1(max_elem) =  [&
   &    0.2272460655_wp,   -0.0596139193_wp,   -0.8162219211_wp,    0.2615605768_wp, & ! -4
   &    0.2862843025_wp,   -0.0173073504_wp,    0.2392693084_wp,    0.0118176989_wp, & ! -8
   &    0.0124716084_wp,   -0.0393471359_wp,   -0.3706940822_wp,   -0.5925288736_wp, & ! -12
   &    0.2215897181_wp,   -0.0315062343_wp,    0.3434101268_wp,    0.2263363232_wp, & ! -16
   &   -0.0492370074_wp,   -0.0538331339_wp,   -0.5261134284_wp,   -0.0304502510_wp, & ! -20
   &   -0.6498395567_wp,   -0.0343697461_wp,    0.0522662650_wp,    0.3083293812_wp, & ! -24
   &    0.5112898789_wp,    0.3426408669_wp,    0.5599083638_wp,   -0.5203072121_wp, & ! -28
   &   -0.2237774345_wp,    0.2495153629_wp,    0.0829715317_wp,   -0.0406155959_wp, & ! -32
   &    0.4647810638_wp,    0.2789378902_wp,    0.0529580888_wp,   -0.0772051837_wp, & ! -36
   &   -0.5679581136_wp,    1.6730047019_wp,   -0.5388627874_wp,   -0.0341739234_wp, & ! -40
   &    0.0361793489_wp,    0.7684617366_wp,    0.5874990182_wp,    0.2345703772_wp, & ! -44
   &    0.3357144208_wp,   -0.0370942202_wp,    0.1103011795_wp,    0.4849810383_wp, & ! -48
   &    0.0449687823_wp,   -0.0474013379_wp,    0.5613791651_wp,    0.2526477340_wp, & ! -52
   &   -0.0478303905_wp,    0.1495971686_wp,   -0.8369485378_wp,   -0.0542019762_wp, & ! -56
   &   -0.1247167330_wp,    0.0502888557_wp,   -0.0637994570_wp,   -0.0635998097_wp, & ! -60
   &   -0.0531518530_wp,   -0.0968731360_wp,   -0.1337199042_wp,   -0.0542471443_wp, & ! -64
   &   -0.0409157548_wp,   -0.0500477946_wp,   -0.0903409668_wp,   -0.0624137643_wp, & ! -68
   &   -0.0964656303_wp,   -0.0547844435_wp,   -0.0901499353_wp,   -0.1044976019_wp, & ! -72
   &    0.0748615545_wp,   -0.0409080635_wp,    0.1422175908_wp,    0.0245211400_wp, & ! -76
   &   -0.0457319531_wp,    0.3030126210_wp,   -0.2156726008_wp,    0.2730278695_wp, & ! -80
   &    0.0727943605_wp,   -0.0482926025_wp,    0.3318716965_wp,    0.3413143387_wp, & ! -84
   &   -0.0400146567_wp,   -0.0668956214_wp,    0.0876924655_wp,    0.1321612468_wp, & ! -88
   &    0.2949305821_wp,    0.4665862315_wp,   -0.1050478058_wp,    0.1694553742_wp, & ! -92
   &   -0.1161440443_wp,    0.2497086409_wp,    0.0496197475_wp,   -0.0424576824_wp, & ! -96
   &   -0.0417067664_wp,   -0.0367557847_wp,   -0.0348215574_wp,   -0.0364418652_wp, & ! -100
   &   -0.0412106965_wp,   -0.0350811514_wp,   -0.0609692849_wp]

   !> Parameter for square-root CN dependence
   real(wp), parameter :: p_k2(max_elem) = [&
   &    0.4237989919_wp,    0.9054329237_wp,    0.7273599522_wp,    1.1780944589_wp, & ! -4
   &    1.2432066425_wp,    0.9021610088_wp,    0.1729209683_wp,    0.3137107221_wp, & ! -8
   &    0.1481611609_wp,    0.0315012512_wp,    0.6948564260_wp,    2.4690486101_wp, & ! -12
   &    1.4646004177_wp,    0.7338489851_wp,    0.2672743388_wp,    0.1255489520_wp, & ! -16
   &    0.3795992335_wp,    0.3132264025_wp,    1.2745899845_wp,    2.7868396415_wp, & ! -20
   &    9.0228693106_wp,    2.0914550048_wp,    1.4881571248_wp,    1.3809899002_wp, & ! -24
   &    4.5228861834_wp,    0.3765445274_wp,    0.3755945109_wp,   -0.3805725578_wp, & ! -28
   &    2.1519012955_wp,    1.0593568196_wp,    3.6358884859_wp,    1.7008902427_wp, & ! -32
   &    0.1456480496_wp,    0.2300853535_wp,    0.1852541944_wp,    0.3599891154_wp, & ! -36
   &    1.2633917188_wp,    4.3256633463_wp,   11.3945013526_wp,    1.1180269980_wp, & ! -40
   &    2.8552731219_wp,    1.1736062932_wp,   10.8207923458_wp,    1.5285170139_wp, & ! -44
   &    5.4935663943_wp,   -0.0457626441_wp,    2.7793846592_wp,    2.2644408359_wp, & ! -48
   &    1.9338445142_wp,    0.8773559861_wp,    0.3739992316_wp,    0.4517948170_wp, & ! -52
   &    0.3041556301_wp,    0.3290807250_wp,    1.9181133121_wp,    3.7648587203_wp, & ! -56
   &   12.5691776403_wp,    0.9136690116_wp,    0.8484898306_wp,    0.5660217937_wp, & ! -60
   &    0.5064170106_wp,    0.4731908893_wp,    1.1891855140_wp,    0.5071098795_wp, & ! -64
   &    0.5824925420_wp,    0.7716650814_wp,    0.5509526594_wp,    0.3183237219_wp, & ! -68
   &    0.2877717791_wp,    0.5259043531_wp,    0.7234286343_wp,    0.9421873340_wp, & ! -72
   &    1.7311512798_wp,    1.9740750517_wp,    4.6810366190_wp,    0.6060090165_wp, & ! -76
   &    2.2323237646_wp,    0.4726592528_wp,    0.5878545519_wp,    2.5750166448_wp, & ! -80
   &    2.3103814679_wp,    0.8246952678_wp,    0.9178667857_wp,    0.4513286733_wp, & ! -84
   &    0.3927258076_wp,    2.3976443914_wp,    0.9832048840_wp,    1.3019615779_wp, & ! -88
   &    0.9237282914_wp,    0.9516418629_wp,    0.9962947453_wp,    1.3291362100_wp, & ! -92
   &    1.0607380277_wp,    1.7496211778_wp,    2.0864700431_wp,    1.1745934488_wp, & ! -96
   &    1.6410458337_wp,    1.3975001895_wp,    1.4464121202_wp,    1.3956575722_wp, & ! -100
   &    1.3080725818_wp,    1.1455579374_wp,    1.1997772501_wp]

   !> Parameter for mixed charge and CN dependence
   real(wp), parameter :: p_k3(max_elem) = [&
   &   -0.1482235383_wp,    0.0797773499_wp,    0.0129050210_wp,    0.8011212947_wp, & ! -4
   &   -0.1403501688_wp,    0.0462480322_wp,    0.1060356466_wp,    0.1795566137_wp, & ! -8
   &    0.3893098244_wp,   -0.0146687422_wp,   -0.3111690904_wp,    0.0320636218_wp, & ! -12
   &   -0.1552600500_wp,   -0.0517809794_wp,    0.0691802709_wp,    0.0491066807_wp, & ! -16
   &   -0.1037356424_wp,    0.0395605098_wp,   -0.0557489720_wp,   -0.0281758849_wp, & ! -20
   &    1.4477859307_wp,    1.4167677387_wp,   -0.2695867206_wp,   -0.2787536365_wp, & ! -24
   &   -0.4474576351_wp,   -0.2086478372_wp,   -0.2234162079_wp,    0.0102244113_wp, & ! -28
   &   -0.5878677603_wp,   -0.3191068261_wp,    0.1185759534_wp,   -0.0332403309_wp, & ! -32
   &    0.1743722400_wp,    0.0316773654_wp,    0.1235826506_wp,    0.3145302762_wp, & ! -36
   &    0.0975461822_wp,    0.1632869075_wp,    1.2804857957_wp,    0.0186250109_wp, & ! -40
   &   -0.0377836152_wp,    0.0194565795_wp,   -0.4838763981_wp,   -0.0573939752_wp, & ! -44
   &   -0.1007465431_wp,    0.0492494427_wp,   -0.1099447559_wp,   -0.0827002830_wp, & ! -48
   &   -0.0228253273_wp,    0.0412157970_wp,    0.0337923911_wp,   -0.0369295411_wp, & ! -52
   &    0.0435379390_wp,    0.1775358173_wp,    0.2804861722_wp,    0.0348907089_wp, & ! -56
   &    1.3465591240_wp,    0.2042327307_wp,    0.1936009387_wp,    0.2300270019_wp, & ! -60
   &    0.1460300571_wp,    0.1738570015_wp,    0.1110552072_wp,    0.1062222399_wp, & ! -64
   &    0.0837731401_wp,    0.1540607852_wp,    0.1449757942_wp,    0.0936359310_wp, & ! -68
   &    0.3750429911_wp,    0.0851732052_wp,    0.0677646919_wp,    0.1532929419_wp, & ! -72
   &    0.0308496571_wp,   -0.1713827638_wp,   -0.2456188248_wp,   -0.0358765818_wp, & ! -76
   &   -0.1115432202_wp,    0.8523468378_wp,    0.0755461658_wp,   -0.3092470689_wp, & ! -80
   &   -0.0776952101_wp,    0.0821136317_wp,    0.0540698456_wp,    0.0614960628_wp, & ! -84
   &   -0.0511469514_wp,   -0.1290672339_wp,   -0.1736992507_wp,    0.1358232940_wp, & ! -88
   &    0.1982890822_wp,    0.2068029681_wp,    0.0910237200_wp,   -0.0609937959_wp, & ! -92
   &   -0.0552882699_wp,    0.1101861282_wp,   -0.0640006225_wp,   -0.0888724586_wp, & ! -96
   &   -0.0580125827_wp,   -0.0426636961_wp,   -0.0631687199_wp,   -0.0696692533_wp, & ! -100
   &   -0.0827656978_wp,   -0.0629399715_wp,   -0.0776964767_wp]

   !> Empirical atomic radii for calculation of the coordination number (modified Pyykkö radii)
   real(wp), parameter :: qvszp_cov_radii(max_elem) = 1.889725949_wp * [&
   & 0.29_wp, 0.46_wp, 1.20_wp, 0.94_wp, 0.77_wp, 0.75_wp, 0.71_wp, 0.63_wp, & ! 1-8
   & 0.64_wp, 0.67_wp, 1.40_wp, 1.25_wp, 1.13_wp, 1.04_wp, 1.10_wp, 1.02_wp, & ! 9-16 
   & 0.99_wp, 0.96_wp, 1.76_wp, 1.54_wp, 1.33_wp, 1.22_wp, 1.21_wp, 1.10_wp, & ! 17-24
   & 1.07_wp, 1.04_wp, 1.00_wp, 0.99_wp, 1.01_wp, 1.09_wp, 1.12_wp, 1.09_wp, & ! 25-32
   & 1.15_wp, 1.10_wp, 1.14_wp, 1.17_wp, 1.89_wp, 1.67_wp, 1.47_wp, 1.39_wp, & ! 33-40
   & 1.32_wp, 1.24_wp, 1.15_wp, 1.13_wp, 1.13_wp, 1.08_wp, 1.15_wp, 1.23_wp, & ! 41-48
   & 1.28_wp, 1.26_wp, 1.26_wp, 1.23_wp, 1.32_wp, 1.31_wp, 2.09_wp, 1.76_wp, & ! 49-56
   & 1.62_wp, 1.47_wp, 1.58_wp, 1.57_wp, 1.56_wp, 1.55_wp, 1.51_wp, 1.52_wp, & ! 57-64
   & 1.51_wp, 1.50_wp, 1.49_wp, 1.49_wp, 1.48_wp, 1.53_wp, 1.46_wp, 1.37_wp, & ! 65-72
   & 1.31_wp, 1.23_wp, 1.18_wp, 1.16_wp, 1.11_wp, 1.12_wp, 1.13_wp, 1.32_wp, & ! 73-80
   & 1.30_wp, 1.30_wp, 1.36_wp, 1.31_wp, 1.38_wp, 1.42_wp, 2.01_wp, 1.81_wp, & ! 81-88
   & 1.67_wp, 1.58_wp, 1.52_wp, 1.53_wp, 1.54_wp, 1.55_wp, 1.49_wp, 1.49_wp, & ! 89-96
   & 1.51_wp, 1.51_wp, 1.48_wp, 1.50_wp, 1.56_wp, 1.58_wp, 1.45_wp]            ! 97-103

   !> Steepness of counting function for the coordination number
   real(wp), parameter   :: kcn = 3.75_wp

   !> CGTO exponents (Exponent of the primitive Gaussian functions)
   real(wp), protected :: exponents(maxg, max_shell, max_elem) = 0.0_wp

   !> CGTO coefficients (Contraction coefficients of the primitive Gaussian functions,
   !> might contain normalization)
   real(wp), protected :: coefficients(maxg, max_shell, max_elem) = 0.0_wp

   !> CGTO coefficients (Contraction coefficients of the primitive Gaussian functions,
   !> might contain normalization)
   real(wp), protected :: coefficients_env(maxg, max_shell, max_elem) = 0.0_wp

contains

   !> Create a new basis set
   subroutine new_basis(self, mol, nshell, cgto, acc)
      !> Instance of the basis set data
      type(qvszp_basis_type), intent(out) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Number of shells per species
      integer, intent(in) :: nshell(:)
      !> Contracted Gaussian basis functions for each shell and species
      class(qvszp_cgto_type), intent(in) :: cgto(:, :)
      !> Calculation accuracy
      real(wp), intent(in) :: acc

      integer :: iat, isp, ish, iao, ii
      real(wp) :: min_alpha

      self%nsh_id = nshell
      self%cgto = cgto
      self%intcut = integral_cutoff(acc)

      ! Make count of shells for each atom
      self%nsh_at = nshell(mol%id)

      ! Create mapping between atoms and shells
      self%nsh = sum(self%nsh_at)
      allocate(self%ish_at(mol%nat), self%sh2at(self%nsh))
      ii = 0
      do iat = 1, mol%nat
         self%ish_at(iat) = ii
         do ish = 1, self%nsh_at(iat)
            self%sh2at(ii+ish) = iat
         end do
         ii = ii + self%nsh_at(iat)
      end do

      ! Make count of spherical orbitals for each shell
      allocate(self%nao_sh(self%nsh))
      do iat = 1, mol%nat
         ii = self%ish_at(iat)
         do ish = 1, self%nsh_at(iat)
            self%nao_sh(ii+ish) = 2*cgto(ish, iat)%ang + 1
         end do
      end do

      ! Create mapping between shells and spherical orbitals, also map directly back to atoms
      self%nao = sum(self%nao_sh)
      allocate(self%iao_sh(self%nsh), self%ao2sh(self%nao), self%ao2at(self%nao))
      ii = 0
      do ish = 1, self%nsh
         self%iao_sh(ish) = ii
         do iao = 1, self%nao_sh(ish)
            self%ao2sh(ii+iao) = ish
            self%ao2at(ii+iao) = self%sh2at(ish)
         end do
         ii = ii + self%nao_sh(ish)
      end do

      ii = 0
      do iat = 1, mol%nat
         isp = mol%id(iat)
         do ish = 1, nshell(isp)
            self%iao_sh(ish+self%ish_at(iat)) = ii
            ii = ii + 2*cgto(ish, iat)%ang + 1
         end do
      end do

      min_alpha = huge(acc)
      do iat = 1, mol%nat
         isp = mol%id(iat)
         do ish = 1, nshell(isp)
            self%maxl = max(self%maxl, cgto(ish, iat)%ang)
            min_alpha = min(min_alpha, minval(cgto(ish, iat)%alpha(:cgto(ish, iat)%nprim)))
         end do
      end do

      self%min_alpha = min_alpha

      ! Allocate the coefficient derivatives 
      do iat = 1, mol%nat
         isp = mol%id(iat)
         do ish = 1, nshell(isp)
            ! Allocate memory for the gradient of the coefficients
            allocate(self%cgto(ish, iat)%dcoeffdr(3, mol%nat, maxg), &
            & self%cgto(ish, iat)%dcoeffdL(3, 3, maxg), source=0.0_wp)
         end do
      end do

      ! Create the coordination number for the environment specific scaling
      call new_ncoord(self%ncoord, mol, cn_type="erf", &
      & kcn=kcn, rcov=qvszp_cov_radii(mol%num))

   end subroutine new_basis


   !> Procedure to setup the CGTOs for each atom in the molecule
   subroutine add_qvszp_cgtos(self, mol, nsh_id, norm)
      !> Instance of the cgto data
      type(qvszp_cgto_type), allocatable, intent(inout) :: self(:,:)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Number of shells at each atom
      integer, allocatable, intent(inout) :: nsh_id(:)
      !> Include normalization in contraction coefficients
      logical, intent(in) :: norm

      integer :: iat, isp, izp, ish, il, nprim
      real(wp) :: alpha(maxg), coeff(maxg)

      !> Initialize full parameter set
      !> set up the array of CGTOs for the molecule of interest
      call setCGTOexponents()
      call setCGTOcoefficients()
      call setCGTOcoefficients_environment()

      nsh_id = nshell(mol%num)

      allocate(self(maxval(nsh_id), mol%nat))
      do iat = 1, mol%nat
         isp = mol%id(iat)
         izp = mol%num(isp)
         do ish = 1, nsh_id(isp)
            il = ang_shell(ish, izp)
            nprim = n_prim(ish, izp)

            ! Gather all necessary parameters
            self(ish, iat)%ang = il
            self(ish, iat)%nprim = nprim
            self(ish, iat)%alpha(1:nprim) = exponents(1:nprim, ish, izp)
            self(ish, iat)%coeff(1:nprim) = coefficients(1:nprim, ish, izp)
            self(ish, iat)%coeff0(1:nprim) = coefficients(1:nprim, ish, izp)
            self(ish, iat)%coeff1(1:nprim) = coefficients_env(1:nprim, ish, izp)
            self(ish, iat)%k1 = p_k1(izp)
            self(ish, iat)%k2 = p_k2(izp)
            self(ish, iat)%k3 = p_k3(izp)
         end do
      end do

      ! Normalize the CGTOs if requested
      if (norm) then
         do iat = 1, mol%nat
            isp = mol%id(iat)
            izp = mol%num(isp)
            do ish = 1, nsh_id(isp)
               call normalize_cgto(self(ish, iat), mol%nat, .false.)
            end do 
         end do 
      end if

   end subroutine add_qvszp_cgtos

   !> Procedure to scale the contraction coefficients of the CGTOs
   subroutine scale_cgto(self, qat, cn, norm, expscal, dqatdr, dqatdL, dcndr, dcndL)
      !> Instance of the cgto data
      class(qvszp_cgto_type), intent(inout) :: self
      !> Atomic charge for the charge scaling of the basis set 
      real(wp), intent(in) :: qat
      !> Coordination number
      real(wp), intent(in) :: cn
      !> Include normalization in contraction coefficients
      logical, intent(in) :: norm
      !> Exponent scaling factor
      real(wp), intent(in), optional :: expscal
      !> Derivative of the specific atomic charge w.r.t. all 3N coordinates (not spin-resolved)
      real(wp), intent(in), optional :: dqatdr(:, :)
      !> Derivative of the specific atomic charge w.r.t. the lattice vectors (not spin-resolved)
      real(wp), intent(in), optional :: dqatdL(:, :)
      !> Derivative of the specific coordination number w.r.t. all 3N coordinates
      real(wp), intent(in), optional :: dcndr(:, :)
      !> Derivative of the specific coordination number w.r.t. the lattice vectors
      real(wp), intent(in), optional :: dcndL(:, :)

      integer :: ipr
      real(wp) :: qeff, normalization
      real(wp), allocatable :: dqeffdr(:, :), dqeffdL(:, :)
      logical grad

      grad = .false.
      if(present(dqatdr) .and. present(dcndr) .and. present(dqatdL) .and. present(dcndL)) then
         allocate(dqeffdr(3, size(dqatdr, 2)), dqeffdL(3, 3))
         grad = .true.
      end if

      ! Calculate the effective charge and derivatives if requested
      qeff = qat - self%k1 * qat**2 + self%k2 * sqrt(cn) + self%k3 * cn * qat
      if(grad) then
         dqeffdr = dqatdr - 2.0_wp * self%k1 * qat * dqatdr &
            & + 0.5_wp * (self%k2 / sqrt(cn)) * dcndr &
            & + self%k3 * cn * dqatdr + self%k3 * qat * dcndr
         dqeffdL = dqatdL - 2.0_wp * self%k1 * qat * dqatdL &
            & + 0.5_wp * (self%k2 / sqrt(cn)) * dcndL &
            & + self%k3 * cn * dqatdL + self%k3 * qat * dcndL
      end if

      ! Scale the coeffients depending on the effective charge
      normalization = 1.0_wp
      do ipr = 1, maxg
         if (norm) then
            normalization = (top*self%alpha(ipr))**0.75_wp &
               & * sqrt(4*self%alpha(ipr))**self%ang / sqrt(dfactorial(self%ang+1))
         endif
         self%coeff(ipr) = (self%coeff0(ipr) + self%coeff1(ipr) * qeff) * normalization
         if(grad) then
            self%dcoeffdr(:, :, ipr) = self%coeff1(ipr) * dqeffdr * normalization
            self%dcoeffdL(:, :, ipr) = self%coeff1(ipr) * dqeffdL * normalization
         end if
      end do    
      
      ! Optionally scale the exponent 
      if(present(expscal)) then
         do ipr = 1, maxg
            self%alpha(ipr) = self%alpha(ipr) * expscal
         end do        
      end if

      ! Normalize the CGTOs if requested
      if(norm) then
         call normalize_cgto(self, size(dqatdr, 2), grad)
      end if

   end subroutine scale_cgto

   subroutine normalize_cgto(self, nat, grad)
      !> Instance of the cgto data
      type(qvszp_cgto_type), intent(inout) :: self
      !> Number of atoms in the molecule
      integer, intent(in) :: nat
      !> Flag to normalize also the coefficient gradient
      logical, intent(in) :: grad

      integer :: ipr
      real(wp) :: normalization, r2, vec(3)
      real(wp) :: overlap(msao(max_shell), msao(max_shell))
      real(wp), allocatable :: dnormalization(:, :), doverlap(:, :, :, :)

      r2 = 0.0_wp
      vec = 0.0_wp
      overlap = 0.0_wp

      ! No screening has to be applied since the gaussians are on the same center
      if(grad) then
         allocate(dnormalization(3, nat), doverlap(3, nat, msao(max_shell), msao(max_shell)), source=0.0_wp)
         ! The single center overlap has a derivative since it is not normalized at this point
         call overlap_grad_decontracted_cgto(nat, self, self, r2, vec, 100.0_wp, overlap, doverlap)
      else
         call overlap_cgto(self, self, r2, vec, 100.0_wp, overlap)
      end if

      ! Normalization constant since all diagonal entries are equivalent
      normalization = 1 / sqrt(overlap(1, 1))
      if(grad) then
         dnormalization = -0.5_wp * doverlap(:, :, 1, 1) / sqrt(overlap(1, 1))**3 
         
         do ipr = 1, self%nprim
            ! Normalization of the coefficient derivative
            self%dcoeffdr(:, :, ipr) = self%dcoeffdr(:, :, ipr) * normalization 
            ! self%dcoeffdL(:,ipr) = self%dcoeffdL(:,ipr) * normalization

            ! Derivative of the normalization constant
            self%dcoeffdr(:, :, ipr) = self%dcoeffdr(:, :, ipr) + self%coeff(ipr) * dnormalization         
            ! self%dcoeffdL(:,ipr) = self%dcoeffdL(:,ipr) + self%coeff(ipr) * dnormalization
         end do 
      end if

      do  ipr = 1, self%nprim
         ! Normalization of the coefficient
         self%coeff(ipr) = self%coeff(ipr) * normalization
      end do

   end subroutine normalize_cgto


   subroutine write_2d_matrix(matrix, name, unit, step)
      implicit none
      real(wp), intent(in) :: matrix(:, :)
      character(len=*), intent(in), optional :: name
      integer, intent(in), optional :: unit
      integer, intent(in), optional :: step
      integer :: d1, d2
      integer :: i, j, k, l, istep, iunit
  
      d1 = size(matrix, dim=1)
      d2 = size(matrix, dim=2)
  
      if (present(unit)) then
        iunit = unit
      else
        iunit = output_unit
      end if
  
      if (present(step)) then
        istep = step
      else
        istep = 6
      end if
  
      if (present(name)) write (iunit, '(/,"matrix printed:",1x,a)') name
  
      do i = 1, d2, istep
        l = min(i + istep - 1, d2)
        write (iunit, '(/,6x)', advance='no')
        do k = i, l
          write (iunit, '(6x,i7,3x)', advance='no') k
        end do
        write (iunit, '(a)')
        do j = 1, d1
          write (iunit, '(i6)', advance='no') j
          do k = i, l
            write (iunit, '(1x,f15.12)', advance='no') matrix(j, k)
          end do
          write (iunit, '(a)')
        end do
      end do
  
   end subroutine write_2d_matrix

   subroutine setCGTOexponents()
      ! set up the array of CGTO exponents during initialization of the basis set data
      exponents(:, :,  1) = reshape([&
      &  337.0050122185_wp,   53.3105310515_wp,   12.2082534814_wp,    3.4188269098_wp,    1.1087901258_wp, &
      &    0.4787177983_wp,    0.2121370877_wp,    0.0690581086_wp,    0.0000000000_wp, & ! shell 1
      &    1.4835839361_wp,    0.3372588943_wp,    0.0864167011_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :,  2) = reshape([&
      &  372.1466733263_wp,   57.4817180773_wp,   13.2069652019_wp,    3.8561112751_wp,    1.2082156387_wp, &
      &    0.4573786432_wp,    0.2389375157_wp,    0.1868021767_wp,    0.0000000000_wp, & ! shell 1
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :,  3) = reshape([&
      &    3.4910597671_wp,    0.5727038734_wp,    0.1419665372_wp,    0.0723362562_wp,    0.0304809874_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    4.7379371606_wp,    1.5009129887_wp,    0.4471972692_wp,    0.1521393085_wp,    0.0577820903_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :,  4) = reshape([&
      &    7.2931906575_wp,    1.2361123226_wp,    0.3999600873_wp,    0.1759126350_wp,    0.0652707036_wp, &
      &    0.0190742484_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.7511643918_wp,    0.6536335064_wp,    0.2120007335_wp,    0.0590415628_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :,  5) = reshape([&
      &    1.5440589040_wp,    1.2089721882_wp,    0.3018669904_wp,    0.1207968920_wp,    0.1132821864_wp, &
      &    0.0776904034_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    8.8208901866_wp,    2.3070729939_wp,    0.6785535316_wp,    0.2422072712_wp,    0.0832231508_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :,  6) = reshape([&
      &    2.2086276291_wp,    1.9179589490_wp,    0.4971530368_wp,    0.1921824673_wp,    0.1505311832_wp, &
      &    0.0781886773_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   24.5973356357_wp,    6.4881608665_wp,    2.2263586744_wp,    0.7768970031_wp,    0.2834247214_wp, &
      &    0.0986282573_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :,  7) = reshape([&
      &    2.6950730961_wp,    2.2494702128_wp,    0.6118516958_wp,    0.3037266656_wp,    0.1716282704_wp, &
      &    0.0715491343_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   27.5476401720_wp,    7.1265662782_wp,    2.3580552040_wp,    0.8450472623_wp,    0.3010630579_wp, &
      &    0.0945745258_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :,  8) = reshape([&
      &    3.4835949431_wp,    2.8598657644_wp,    0.7406810290_wp,    0.3449269192_wp,    0.2619785329_wp, &
      &    0.0946163504_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   37.6414713255_wp,    9.7130885617_wp,    3.1371771476_wp,    1.1213805162_wp,    0.3980779265_wp, &
      &    0.1230477710_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :,  9) = reshape([&
      &    6.8188387931_wp,    5.1934068072_wp,    1.3443160267_wp,    0.6789538085_wp,    0.4294176500_wp, &
      &    0.1740270721_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   54.4257988065_wp,   15.0636200476_wp,    4.8170002590_wp,    1.7248652235_wp,    0.5868819738_wp, &
      &    0.1919528054_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 10) = reshape([&
      &    8.8923160490_wp,    7.6627593978_wp,    1.7842346597_wp,    0.6991767519_wp,    0.3083126139_wp, &
      &    0.2129764427_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   73.4944041184_wp,   18.9102168118_wp,    6.1870669767_wp,    2.2209903931_wp,    0.7548834305_wp, &
      &    0.2448850968_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 11) = reshape([&
      &    2.1261969277_wp,    0.5405207606_wp,    0.0731107333_wp,    0.0267822692_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.4303461738_wp,    0.2556368092_wp,    0.0860967844_wp,    0.0263882077_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 12) = reshape([&
      &    2.6078220503_wp,    0.8030261667_wp,    0.1299956568_wp,    0.0523718547_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.8588160999_wp,    0.1754920593_wp,    0.0569256804_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.4747794463_wp,    0.1150613841_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 13) = reshape([&
      &    3.1961316172_wp,    1.0470244669_wp,    0.3770937988_wp,    0.1632089175_wp,    0.0640732091_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.9449487491_wp,    0.3983949738_wp,    0.1527350144_wp,    0.0592354825_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.5970683998_wp,    0.1724426664_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 14) = reshape([&
      &    4.4809894623_wp,    1.3118331904_wp,    0.3049298502_wp,    0.2010748901_wp,    0.0769532886_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.2036915405_wp,    0.4884596995_wp,    0.1837801098_wp,    0.0650477655_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    1.0625303376_wp,    0.2885597553_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 15) = reshape([&
      &    8.0175344409_wp,    1.4444963516_wp,    0.4552363098_wp,    0.1935391086_wp,    0.0674143274_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.1707461912_wp,    0.7735984606_wp,    0.3310762816_wp,    0.1608853026_wp,    0.0671391937_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    1.6034741474_wp,    0.3628987891_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 16) = reshape([&
      &    7.5364413985_wp,    1.9294689410_wp,    0.5221222113_wp,    0.2305237994_wp,    0.0862615041_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.6434259779_wp,    1.1519870166_wp,    0.4274171111_wp,    0.1941236769_wp,    0.0861259080_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    1.5683376734_wp,    0.3588872912_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 17) = reshape([&
      &   16.4401524867_wp,    2.1283197520_wp,    0.7039922357_wp,    0.2933317022_wp,    0.1112209830_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.6783569387_wp,    1.6724835438_wp,    0.6520547879_wp,    0.2616677727_wp,    0.1031425929_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    2.0726630450_wp,    0.4763394499_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 18) = reshape([&
      &   16.9311039744_wp,    2.2930114551_wp,    1.7775061724_wp,    0.4735254551_wp,    0.1728644199_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    3.2415547019_wp,    1.1204952462_wp,    0.4404486462_wp,    0.2349439290_wp,    0.1254954683_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    1.7670574130_wp,    0.3722739347_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 19) = reshape([&
      &    0.9096291775_wp,    0.3075153536_wp,    0.0452932389_wp,    0.0202391205_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.2097555495_wp,    0.0817195861_wp,    0.0219522293_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 20) = reshape([&
      &    0.8069471761_wp,    0.4946552308_wp,    0.0500345462_wp,    0.0387154055_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.3538365975_wp,    0.0802283356_wp,    0.0237643919_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    1.5523774589_wp,    0.4588094612_wp,    0.1224154789_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 21) = reshape([&
      &    0.7404248520_wp,    0.4899005688_wp,    0.1013651184_wp,    0.0385135596_wp,    0.0120083001_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.5405272830_wp,    0.0879253830_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   37.4491137329_wp,    9.6603347463_wp,    3.1174995835_wp,    1.1112292201_wp,    0.3775979697_wp, &
      &    0.1266399593_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 22) = reshape([&
      &    0.9697757290_wp,    0.5340801100_wp,    0.2416035765_wp,    0.1075356765_wp,    0.0465975840_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.5796825299_wp,    0.1234237064_wp,    0.0730447552_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   41.1761976464_wp,   11.1358285314_wp,    3.6875504496_wp,    1.3338185901_wp,    0.4511802208_wp, &
      &    0.1456642297_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 23) = reshape([&
      &    0.9028729898_wp,    0.5698447744_wp,    0.3038012890_wp,    0.1134239408_wp,    0.0360491166_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.7631683919_wp,    0.1248613589_wp,    0.0561105580_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   41.1703945214_wp,   11.7905049594_wp,    4.0451713506_wp,    1.4981538118_wp,    0.5275160351_wp, &
      &    0.1797703925_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 24) = reshape([&
      &    1.0951382808_wp,    0.2908102872_wp,    0.4709988662_wp,    0.1117239766_wp,    0.0310474525_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.7793776741_wp,    0.1252213107_wp,    0.0228203514_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   49.6314851957_wp,   13.8226667153_wp,    4.7919721935_wp,    1.8111120193_wp,    0.6461675807_wp, &
      &    0.2161466076_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 25) = reshape([&
      &    1.2710215875_wp,    0.7088872015_wp,    0.2589864707_wp,    0.1630347400_wp,    0.0629368683_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.8372023834_wp,    0.1582365482_wp,    0.0509274840_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   49.7708386060_wp,   14.4934220080_wp,    5.1530869746_wp,    1.9465732365_wp,    0.6956418205_wp, &
      &    0.2360335316_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 26) = reshape([&
      &    1.6572154176_wp,    0.8131349274_wp,    0.4125700590_wp,    0.1043183705_wp,    0.0185408395_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.8817530766_wp,    0.1524026741_wp,    0.0409311933_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   53.6873317631_wp,   15.6563029879_wp,    5.5728284406_wp,    2.0953165211_wp,    0.7410790629_wp, &
      &    0.2456877765_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 27) = reshape([&
      &    1.4974029220_wp,    0.8306274545_wp,    0.4052311256_wp,    0.1178195973_wp,    0.0280679884_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.9955620249_wp,    0.1636273297_wp,    0.0568760750_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   59.7428956025_wp,   18.0706494974_wp,    6.5261362709_wp,    2.4801217833_wp,    0.8932838719_wp, &
      &    0.2962976877_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 28) = reshape([&
      &    1.4207925823_wp,    0.9553898154_wp,    0.2026327115_wp,    0.1060225662_wp,    0.0270998585_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.0980015291_wp,    0.1794387155_wp,    0.0630422263_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   64.1765855187_wp,   18.3142207840_wp,    6.4425915978_wp,    2.3837181928_wp,    0.8223463933_wp, &
      &    0.2421737828_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 29) = reshape([&
      &    1.4953072401_wp,    1.0020794930_wp,    0.3227739980_wp,    0.1276219035_wp,    0.0558545578_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.0189632325_wp,    0.1888184159_wp,    0.0733905246_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   67.7505543157_wp,   20.2133476456_wp,    7.3284136186_wp,    2.7771207677_wp,    0.9898811925_wp, &
      &    0.3059680593_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 30) = reshape([&
      &    1.8204663224_wp,    1.1327458207_wp,    0.3150460049_wp,    0.1486172772_wp,    0.0512931422_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.2237530210_wp,    0.1836401118_wp,    0.0474059648_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 31) = reshape([&
      &    3.3018831782_wp,    1.7261110515_wp,    0.2149517181_wp,    0.1849832576_wp,    0.0701476710_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.0552674413_wp,    0.7750459739_wp,    0.1912102468_wp,    0.0666015979_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.3363760028_wp,    0.1087662159_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 32) = reshape([&
      &    3.2260486960_wp,    1.9113375882_wp,    0.2844558436_wp,    0.1934849783_wp,    0.0715603147_wp, &
      &    0.0597175231_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.2039608566_wp,    0.8028078617_wp,    0.2290148936_wp,    0.0791035546_wp,    0.0731817815_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.3848760903_wp,    0.1387639644_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 33) = reshape([&
      &    3.4199269206_wp,    1.8840896804_wp,    0.3644420076_wp,    0.2215550465_wp,    0.1279642877_wp, &
      &    0.0706925317_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.3421339764_wp,    0.9607964455_wp,    0.3334829227_wp,    0.1702027767_wp,    0.0659750113_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.2861382957_wp,    0.0802841315_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 34) = reshape([&
      &    3.6934659349_wp,    2.1866375691_wp,    0.3572912768_wp,    0.1902160148_wp,    0.1289984689_wp, &
      &    0.0691751289_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.4641392362_wp,    1.0310897408_wp,    0.3677824918_wp,    0.1851325719_wp,    0.0801214756_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.3149667465_wp,    0.1107870749_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 35) = reshape([&
      &    4.6018981817_wp,    2.2695913083_wp,    0.4295141419_wp,    0.2434407790_wp,    0.1423478157_wp, &
      &    0.0766627777_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.7679775563_wp,    1.3186715921_wp,    0.4687569925_wp,    0.2087338113_wp,    0.0901306557_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.4524437537_wp,    0.1809741952_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 36) = reshape([&
      &    4.1549129477_wp,    2.8168550350_wp,    0.4621616822_wp,    0.2202312656_wp,    0.1674883460_wp, &
      &    0.0684611374_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.0787378157_wp,    1.4841675497_wp,    0.5630822122_wp,    0.2519355994_wp,    0.1037780649_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.4883089480_wp,    0.1737894756_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 37) = reshape([&
      &    0.5570613598_wp,    0.3483923768_wp,    0.0308045042_wp,    0.0166736600_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1636920577_wp,    0.0749579066_wp,    0.0224242436_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 38) = reshape([&
      &    0.6201677265_wp,    0.3917110150_wp,    0.0570874630_wp,    0.0300150324_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1941351247_wp,    0.1013068251_wp,    0.0329240872_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.8240801848_wp,    0.2797884665_wp,    0.0861235400_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 39) = reshape([&
      &    0.6636926668_wp,    0.3749945899_wp,    0.1016698645_wp,    0.0568591814_wp,    0.0330378428_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.3558004721_wp,    0.0716976747_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    8.5068117647_wp,    2.2308087109_wp,    0.9283189428_wp,    0.3765930358_wp,    0.1560030182_wp, &
      &    0.0655774736_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 40) = reshape([&
      &    0.5097879863_wp,    0.3407163126_wp,    0.3224715162_wp,    0.1564767409_wp,    0.0484822988_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.3215997408_wp,    0.1367038485_wp,    0.0588868984_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    9.9445666726_wp,    2.2820649381_wp,    0.9351962665_wp,    0.3779862649_wp,    0.1556926700_wp, &
      &    0.0633356617_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 41) = reshape([&
      &    0.7801941256_wp,    0.4283293090_wp,    0.1798252664_wp,    0.1109596022_wp,    0.0487692913_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.4258370219_wp,    0.0915136842_wp,    0.0674062698_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   10.2134004528_wp,    1.8790153326_wp,    0.8975842658_wp,    0.3861929737_wp,    0.1528975160_wp, &
      &    0.0667442529_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 42) = reshape([&
      &    0.7209507628_wp,    0.4631058108_wp,    0.1294540850_wp,    0.0456107839_wp,    0.0361701269_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.3271246927_wp,    0.2102661767_wp,    0.0615197268_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   11.4413105195_wp,    3.2668359146_wp,    1.3099653506_wp,    0.5560299472_wp,    0.2424854011_wp, &
      &    0.1100552828_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 43) = reshape([&
      &    0.9915409334_wp,    0.4560288487_wp,    0.1658845151_wp,    0.1583503624_wp,    0.0564879250_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.5406335343_wp,    0.0922334765_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   14.6904983733_wp,    3.9428852607_wp,    1.6156628715_wp,    0.7638891071_wp,    0.3491413249_wp, &
      &    0.1425047750_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 44) = reshape([&
      &    1.0490528210_wp,    0.6000998177_wp,    0.1169604882_wp,    0.0830887561_wp,    0.0390944365_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.5281721682_wp,    0.1349946825_wp,    0.0418914420_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   14.1570444463_wp,    3.8624522400_wp,    1.7205843720_wp,    0.8092993808_wp,    0.3541739276_wp, &
      &    0.1358271050_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 45) = reshape([&
      &    1.1993100378_wp,    0.6357123830_wp,    0.1552525471_wp,    0.0650680783_wp,    0.0543454893_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.4911273047_wp,    0.1589247853_wp,    0.0602135387_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   17.1172085044_wp,    5.5168983213_wp,    1.8745132663_wp,    0.8471905715_wp,    0.3630990879_wp, &
      &    0.1390332678_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 46) = reshape([&
      &    1.2232518777_wp,    0.6585080596_wp,    0.1207094402_wp,    0.0560533164_wp,    0.0188479414_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.5613490923_wp,    0.1345989577_wp,    0.0515572793_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   18.4164924558_wp,    3.6729952878_wp,    1.8586760589_wp,    0.8584679173_wp,    0.3670597614_wp, &
      &    0.1358465633_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 47) = reshape([&
      &    1.6542562745_wp,    0.6764920403_wp,    0.1337395602_wp,    0.0631344237_wp,    0.0441222219_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.5184959249_wp,    0.1644034514_wp,    0.0588852857_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   17.6013877351_wp,    4.5305317828_wp,    2.3199306701_wp,    1.0812950139_wp,    0.4479984557_wp, &
      &    0.1697934801_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 48) = reshape([&
      &    1.3361916176_wp,    0.8680090265_wp,    0.3623183322_wp,    0.1276117136_wp,    0.0489971823_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.7687446034_wp,    0.1555020535_wp,    0.0395657892_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 49) = reshape([&
      &    1.3312962151_wp,    0.9531669045_wp,    0.3193324693_wp,    0.1437739173_wp,    0.0612304330_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.4427843003_wp,    1.2093296930_wp,    0.1921867827_wp,    0.0611238738_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1036681014_wp,    0.0802156566_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 50) = reshape([&
      &    2.2727327462_wp,    1.2011451235_wp,    0.4213517133_wp,    0.1808863445_wp,    0.1335372526_wp, &
      &    0.0488179294_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.3108181035_wp,    1.6014828842_wp,    0.2180851901_wp,    0.0805628099_wp,    0.0593288795_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1929955327_wp,    0.0890528038_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 51) = reshape([&
      &    1.9506458255_wp,    1.4010102528_wp,    0.3381446226_wp,    0.2232807839_wp,    0.1464752207_wp, &
      &    0.0511828620_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.0686051625_wp,    1.3877670528_wp,    0.2784128300_wp,    0.1325725818_wp,    0.0559347514_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1857722673_wp,    0.1565279907_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 52) = reshape([&
      &    2.2679141921_wp,    1.5930708311_wp,    0.3355932222_wp,    0.2189037043_wp,    0.1804550515_wp, &
      &    0.0784602088_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.3208698414_wp,    1.5520446211_wp,    0.3596100840_wp,    0.1658344506_wp,    0.0652014390_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1881948376_wp,    0.0412966287_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 53) = reshape([&
      &    2.5829786025_wp,    1.6352048003_wp,    0.3339789084_wp,    0.2719730234_wp,    0.2080233318_wp, &
      &    0.0995937345_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.9371458179_wp,    0.6647456229_wp,    0.2148231135_wp,    0.1887970512_wp,    0.0827058985_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.2474029037_wp,    0.1240023810_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 54) = reshape([&
      &    2.6656085130_wp,    1.7642851869_wp,    0.3876272308_wp,    0.3362948782_wp,    0.1802314033_wp, &
      &    0.0715665673_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.0366864486_wp,    0.7726965515_wp,    0.2753071369_wp,    0.1853956128_wp,    0.0774713138_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.2057564060_wp,    0.0658635473_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 55) = reshape([&
      &    0.3494861553_wp,    0.2922119502_wp,    0.0231876896_wp,    0.0198269471_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1233530403_wp,    0.0983126528_wp,    0.0214548993_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 56) = reshape([&
      &    0.4501221558_wp,    0.3348673432_wp,    0.0354643131_wp,    0.0291592698_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1321386724_wp,    0.1018410312_wp,    0.0252662227_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.3277409002_wp,    0.1145675228_wp,    0.0315162358_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 57) = reshape([&
      &    0.3577274426_wp,    0.2071073072_wp,    0.1269988628_wp,    0.0359388367_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1864590486_wp,    0.0640135382_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    3.1520136611_wp,    0.5662056822_wp,    0.3231855058_wp,    0.1355847988_wp,    0.0452017436_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 58) = reshape([&
      &    7.6267298509_wp,    3.8419733114_wp,    1.3375407638_wp,    0.4638722575_wp,    0.2851468066_wp, &
      &    0.0493823263_wp,    0.0166395057_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.7032462712_wp,    0.8268294026_wp,    0.3241783612_wp,    0.0323216787_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   14.7744185697_wp,    8.4637344375_wp,    4.0672731372_wp,    2.0016722206_wp,    0.4958243935_wp, &
      &    0.1577075244_wp,    0.0415553782_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &   33.2454226169_wp,   11.5647490556_wp,    4.4433226937_wp,    1.7498818631_wp,    0.6520481705_wp, &
      &    0.2413152779_wp,    0.0815497049_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 59) = reshape([&
      &    6.8152630044_wp,    3.6176788571_wp,    1.8259715275_wp,    0.6252111136_wp,    0.3111609509_wp, &
      &    0.0497151345_wp,    0.0116559097_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    3.6132621171_wp,    1.0845784007_wp,    0.2067284822_wp,    0.0440968278_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   11.5286372250_wp,    4.5324225317_wp,    1.8127849130_wp,    0.5338776853_wp,    0.1774779233_wp, &
      &    0.0534240170_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &   38.3851186825_wp,   13.7951061925_wp,    5.1605007759_wp,    1.9208967128_wp,    0.6947016522_wp, &
      &    0.2249687743_wp,    0.0815235727_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 60) = reshape([&
      &    7.6309027915_wp,    3.8166599570_wp,    1.8462568112_wp,    0.4865447110_wp,    0.2650200404_wp, &
      &    0.0553944916_wp,    0.0207998752_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.6596927326_wp,    1.2301282056_wp,    0.1865917883_wp,    0.0422426598_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    9.9485076355_wp,    4.8490000571_wp,    2.1138354422_wp,    0.5501122146_wp,    0.1650167860_wp, &
      &    0.0596752195_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &   38.6966072595_wp,   14.0298565647_wp,    5.4681679722_wp,    2.0966146329_wp,    0.7752501074_wp, &
      &    0.2750390721_wp,    0.0913863665_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 61) = reshape([&
      &    8.3733147950_wp,    3.9784105549_wp,    1.8364187167_wp,    0.5466430105_wp,    0.2916681484_wp, &
      &    0.0600809486_wp,    0.0162889793_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    3.4068105153_wp,    1.0173993223_wp,    0.2354323114_wp,    0.0438389171_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   11.0410552020_wp,    8.1952063080_wp,    5.1526805955_wp,    1.6731743259_wp,    0.5468765162_wp, &
      &    0.1902246725_wp,    0.0463467561_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &   40.9924007653_wp,   14.6848383748_wp,    5.6541300606_wp,    2.2121225791_wp,    0.8401058478_wp, &
      &    0.3123376871_wp,    0.0852891569_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 62) = reshape([&
      &    9.9693675744_wp,    3.1845090194_wp,    2.5756556242_wp,    0.4183887364_wp,    0.5267128792_wp, &
      &    0.0781447212_wp,    0.0218912075_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    3.9012494139_wp,    1.0643296540_wp,    0.2572184216_wp,    0.1438422298_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   11.3316952182_wp,    5.2638412227_wp,    3.0188848745_wp,    0.6197113013_wp,    0.2158338157_wp, &
      &    0.0566889345_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &   45.5091386140_wp,   16.3495667350_wp,    6.2992633282_wp,    2.4682328550_wp,    0.9619766433_wp, &
      &    0.3578351734_wp,    0.0536195282_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 63) = reshape([&
      &    9.3751438182_wp,    4.6452527676_wp,    2.1863117142_wp,    0.6142749070_wp,    0.3321365030_wp, &
      &    0.0562282121_wp,    0.0157522377_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.1484660488_wp,    0.8772074068_wp,    0.2194329134_wp,    0.0519605782_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   13.4474016699_wp,    4.9155146590_wp,    0.9526227516_wp,    0.5683894250_wp,    0.1704700455_wp, &
      &    0.0419348364_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &   46.2647739651_wp,   16.7240709127_wp,    6.5316581338_wp,    2.5835877894_wp,    0.9667385613_wp, &
      &    0.3234077422_wp,    0.0577771911_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 64) = reshape([&
      &    9.3354451210_wp,    4.4082279966_wp,    1.9558022299_wp,    0.6419704707_wp,    0.3909632299_wp, &
      &    0.0507236709_wp,    0.0167501620_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.8330520628_wp,    1.0860843853_wp,    0.1951418653_wp,    0.0356533156_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   16.3899702127_wp,    5.1480753750_wp,    1.0766426775_wp,    0.6356283506_wp,    0.1742210148_wp, &
      &    0.0457407054_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &   50.2449377090_wp,   18.3921842814_wp,    7.2336473746_wp,    2.8799215848_wp,    1.1110827857_wp, &
      &    0.3887170489_wp,    0.0886370135_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 65) = reshape([&
      &    8.7654853413_wp,    3.7787856154_wp,    2.2842716814_wp,    0.5326878958_wp,    0.3386853603_wp, &
      &    0.0545687743_wp,    0.0137134336_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    3.4731633241_wp,    1.0980371362_wp,    0.2233108724_wp,    0.0561449534_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   18.4257376205_wp,    7.8527478855_wp,    5.6417234048_wp,    1.4415445816_wp,    0.6989522406_wp, &
      &    0.1912358618_wp,    0.0429641659_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &   50.8923675244_wp,   18.4849032097_wp,    7.1344626226_wp,    2.7457406431_wp,    1.0187331257_wp, &
      &    0.3564467897_wp,    0.0783136299_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 66) = reshape([&
      &    9.9130613806_wp,    4.6415489128_wp,    2.4166391806_wp,    0.5485967997_wp,    0.3365521830_wp, &
      &    0.0582870988_wp,    0.0153131631_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    4.0458013975_wp,    0.9222651741_wp,    0.1873489341_wp,    0.0509551057_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   28.7701535729_wp,    6.6923660666_wp,    2.5046639455_wp,    0.6599077084_wp,    0.2246200449_wp, &
      &    0.0623359401_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &   57.4711675112_wp,   21.2081424907_wp,    8.3832151002_wp,    3.3403307445_wp,    1.2593703704_wp, &
      &    0.4354047088_wp,    0.1435592511_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 67) = reshape([&
      &   10.1209579895_wp,    5.3611014863_wp,    2.8221168551_wp,    0.5863857548_wp,    0.2835396060_wp, &
      &    0.0632609409_wp,    0.0186268562_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    3.7419041641_wp,    1.0183065818_wp,    0.1826312099_wp,    0.0626130022_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   29.9246007212_wp,    6.6378039866_wp,    1.4193965055_wp,    0.6884600064_wp,    0.1904244132_wp, &
      &    0.0493114280_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &   62.5869637712_wp,   23.1237099515_wp,    9.1276994404_wp,    3.6243016344_wp,    1.3781838795_wp, &
      &    0.4770224959_wp,    0.1310517979_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 68) = reshape([&
      &   11.7049012316_wp,    4.8491277634_wp,    2.7268784786_wp,    0.6369420260_wp,    0.4451189398_wp, &
      &    0.0611796172_wp,    0.0164493522_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.5540988596_wp,    1.1021463158_wp,    0.2111850582_wp,    0.0571794137_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   39.1634638404_wp,    8.4947810489_wp,    3.3164897433_wp,    0.6718831616_wp,    0.1821191783_wp, &
      &    0.0441426248_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &   63.3429996524_wp,   23.5007381172_wp,    9.4321413203_wp,    3.7404886251_wp,    1.4504431970_wp, &
      &    0.5196305601_wp,    0.1432044608_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 69) = reshape([&
      &   10.9934977441_wp,    4.3800324475_wp,    1.8149011030_wp,    0.8631818942_wp,    0.3761304679_wp, &
      &    0.0574203655_wp,    0.0179114763_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    3.4248117423_wp,    1.0680724040_wp,    0.2986593943_wp,    0.0714662266_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   34.1701660375_wp,    7.1942403540_wp,    0.9993688215_wp,    0.7083266494_wp,    0.2155872462_wp, &
      &    0.0574971860_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &   64.1651369178_wp,   23.8030320258_wp,    9.4885687792_wp,    3.7736675912_wp,    1.4118515885_wp, &
      &    0.4768339257_wp,    0.1370681146_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 70) = reshape([&
      &   10.8955562918_wp,    5.6915531474_wp,    2.4235945804_wp,    0.7413749984_wp,    0.4020683255_wp, &
      &    0.0686723158_wp,    0.0198263641_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    3.2204435325_wp,    0.9620125488_wp,    0.2176773444_wp,    0.0556719203_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   36.8844642802_wp,    7.8918735723_wp,    1.3313671398_wp,    0.6937669257_wp,    0.1920421078_wp, &
      &    0.0524593556_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &   65.1088419639_wp,   23.9891386911_wp,    9.6442223448_wp,    3.8271285786_wp,    1.4362862069_wp, &
      &    0.4817960703_wp,    0.1280618563_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 71) = reshape([&
      &   11.9119675403_wp,    5.7445069952_wp,    1.8261910330_wp,    0.9353464146_wp,    0.3483936487_wp, &
      &    0.0737007282_wp,    0.0155748421_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    3.2282451982_wp,    1.2962669393_wp,    0.2260711155_wp,    0.0499798305_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   37.7168422826_wp,    8.3636916419_wp,    2.0431527549_wp,    0.6925441676_wp,    0.1871723502_wp, &
      &    0.0519659254_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &   70.8824405084_wp,   26.0470687242_wp,   10.4160090170_wp,    4.1497991197_wp,    1.5699613186_wp, &
      &    0.5272796899_wp,    0.0463966297_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 72) = reshape([&
      &    0.6556753520_wp,    0.4655928403_wp,    0.1119422838_wp,    0.0351664588_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.3560492293_wp,    0.1437655681_wp,    0.0549054073_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    6.5551599573_wp,    0.9888874488_wp,    0.3987652642_wp,    0.1633444220_wp,    0.0624393638_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 73) = reshape([&
      &    0.6607358220_wp,    0.4402949513_wp,    0.2577729208_wp,    0.0964194714_wp,    0.0392973632_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.3741485044_wp,    0.1539067361_wp,    0.0639480493_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    6.0816668691_wp,    1.0727726060_wp,    0.4399462963_wp,    0.1800986802_wp,    0.0698297015_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 74) = reshape([&
      &    1.1003356644_wp,    0.5869660475_wp,    0.1434823407_wp,    0.0659932944_wp,    0.0380847147_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.4536963553_wp,    0.1487348778_wp,    0.0640512178_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    7.4179724676_wp,    1.1042084986_wp,    0.4470530779_wp,    0.1805880204_wp,    0.0828295617_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 75) = reshape([&
      &    1.1077648679_wp,    0.6126121951_wp,    0.1439787951_wp,    0.0621195584_wp,    0.0488926380_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.4171001651_wp,    0.2050183678_wp,    0.0849021459_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    6.9912387455_wp,    1.6542427116_wp,    1.1006002843_wp,    0.4418243983_wp,    0.1586051479_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 76) = reshape([&
      &    1.0031084376_wp,    0.6627665611_wp,    0.1414598806_wp,    0.0572714415_wp,    0.0249572501_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.5080063820_wp,    0.1473602271_wp,    0.0568106458_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    6.8202761038_wp,    1.2949105939_wp,    0.5531647435_wp,    0.2426876005_wp,    0.1030452281_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 77) = reshape([&
      &    1.0955792053_wp,    0.6326640163_wp,    0.2648085670_wp,    0.0861084204_wp,    0.0557504035_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.3493745241_wp,    0.2965819635_wp,    0.0810046021_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    8.4175798553_wp,    1.3302072836_wp,    0.5761463457_wp,    0.2685646358_wp,    0.1155698754_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 78) = reshape([&
      &    1.3191897603_wp,    0.7488858401_wp,    0.1253479012_wp,    0.0904705929_wp,    0.0357746981_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.5072903541_wp,    0.1573879397_wp,    0.0528464911_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    8.4852796212_wp,    1.4233998509_wp,    0.6125301378_wp,    0.2610340357_wp,    0.1014231424_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 79) = reshape([&
      &    1.2921326063_wp,    0.7424125197_wp,    0.1492726529_wp,    0.0556519435_wp,    0.0470371144_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.5023436706_wp,    0.1802689959_wp,    0.0559685199_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    8.8896798355_wp,    1.4953607053_wp,    0.6615743075_wp,    0.3054506088_wp,    0.1332194553_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 80) = reshape([&
      &    1.4311065622_wp,    0.8909389260_wp,    0.3089310173_wp,    0.1384973599_wp,    0.0601406748_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.6000012770_wp,    0.1468846154_wp,    0.0265839764_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 81) = reshape([&
      &    1.4214809893_wp,    0.9734689493_wp,    0.2554521187_wp,    0.1699892882_wp,    0.0631270117_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.1797953179_wp,    1.0632755773_wp,    0.1788776404_wp,    0.0601538861_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.7893481417_wp,    0.0966320789_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 82) = reshape([&
      &    1.4297680962_wp,    1.0378214188_wp,    0.5202036242_wp,    0.1862514692_wp,    0.0627056312_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.2748970817_wp,    1.0116674642_wp,    0.2209187666_wp,    0.0762582312_wp,    0.0589215102_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.2013719395_wp,    0.0842980022_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 83) = reshape([&
      &    1.5614257341_wp,    1.0078141617_wp,    0.4021503281_wp,    0.1876698799_wp,    0.1710607382_wp, &
      &    0.0552223895_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.3684575253_wp,    1.1902687785_wp,    0.2882332519_wp,    0.1466745227_wp,    0.0616442984_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1243462490_wp,    0.1122497776_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 84) = reshape([&
      &    1.7991270466_wp,    1.3549882454_wp,    0.6942960354_wp,    0.2532694700_wp,    0.1465403361_wp, &
      &    0.0728071673_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.6591185448_wp,    1.2416040624_wp,    0.2513493327_wp,    0.1432553755_wp,    0.0761738436_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1597043344_wp,    0.0425127067_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 85) = reshape([&
      &    1.9477903997_wp,    1.3024571193_wp,    0.5595178331_wp,    0.3871054525_wp,    0.2308105647_wp, &
      &    0.1060220278_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.9258620584_wp,    1.2235357817_wp,    0.3307997893_wp,    0.1459441139_wp,    0.0662087939_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1786859802_wp,    0.1695330234_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 86) = reshape([&
      &    1.9463155734_wp,    1.4587807631_wp,    0.3693756217_wp,    0.2131516658_wp,    0.1544345131_wp, &
      &    0.1039965012_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.7773177967_wp,    1.4633849188_wp,    0.3659823508_wp,    0.1563626118_wp,    0.0566094420_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1624420471_wp,    0.0278514465_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 87) = reshape([&
      &   10.1728435892_wp,    6.4450043340_wp,    2.8736065747_wp,    1.3484757605_wp,    0.5988150795_wp, &
      &    0.2099404232_wp,    0.0429268794_wp,    0.0228034252_wp,    0.0091671830_wp, & ! shell 1
      &    5.4917500437_wp,    2.5993866274_wp,    1.4386273644_wp,    0.4269289293_wp,    0.1593272308_wp, &
      &    0.0704211411_wp,    0.0334030112_wp,    0.0128930232_wp,    0.0056891250_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 88) = reshape([&
      &   10.7446244528_wp,    6.7151384330_wp,    2.8154390397_wp,    1.4811064759_wp,    0.5515954870_wp, &
      &    0.2390849042_wp,    0.0525426006_wp,    0.0193822644_wp,    0.0103482565_wp, & ! shell 1
      &    5.2505537935_wp,    2.5547948088_wp,    1.1785584062_wp,    0.3234258277_wp,    0.1642752210_wp, &
      &    0.0686152294_wp,    0.0362512648_wp,    0.0149048474_wp,    0.0069597827_wp, & ! shell 2
      &    1.1403098280_wp,    0.4528863207_wp,    0.2456968223_wp,    0.0831410947_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 89) = reshape([&
      &    6.0016980622_wp,    2.8633121428_wp,    1.3534218353_wp,    0.4170886032_wp,    0.3753346761_wp, &
      &    0.0509603723_wp,    0.0171013138_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.2805718199_wp,    0.9749910066_wp,    0.3013480027_wp,    0.0271505530_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   11.1399896296_wp,    4.7662866030_wp,    2.1195753327_wp,    0.3819263428_wp,    0.1545764845_wp, &
      &    0.0507141469_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    7.6065211097_wp,    3.0604092782_wp,    1.9904936030_wp,    0.7801257005_wp,    0.2772370057_wp, &
      &    0.1912080117_wp,    0.0509806616_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 90) = reshape([&
      &    6.0726048657_wp,    3.0624841920_wp,    1.4227710764_wp,    0.4470441200_wp,    0.3687904549_wp, &
      &    0.0529461033_wp,    0.0162034039_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.1991166978_wp,    0.9986774230_wp,    0.3112299473_wp,    0.0274867620_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   11.2682601861_wp,    4.8303051423_wp,    2.3261890625_wp,    0.4075567776_wp,    0.1584398338_wp, &
      &    0.0470337317_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    8.1615276289_wp,    3.2394158997_wp,    2.3100427369_wp,    0.8903941566_wp,    0.3259365786_wp, &
      &    0.1955627117_wp,    0.0541752166_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 91) = reshape([&
      &    6.4399751769_wp,    3.1231123241_wp,    1.4469389038_wp,    0.4521711961_wp,    0.3912424021_wp, &
      &    0.0536239335_wp,    0.0152588360_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.1107736823_wp,    1.0520097208_wp,    0.3236308765_wp,    0.0267541213_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   10.2403680235_wp,    5.2504253132_wp,    2.5032980250_wp,    0.3971381099_wp,    0.1572975645_wp, &
      &    0.0490846349_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    8.0403341724_wp,    3.8430158094_wp,    2.4513526225_wp,    1.0204806363_wp,    0.4124579602_wp, &
      &    0.2017241524_wp,    0.0551909080_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 92) = reshape([&
      &    6.5230465662_wp,    3.2449573555_wp,    1.4024699312_wp,    0.4580789044_wp,    0.3779013570_wp, &
      &    0.0548418725_wp,    0.0166269112_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.1973408211_wp,    1.0926066815_wp,    0.3174102614_wp,    0.0298676070_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    8.6141505335_wp,    5.8134832544_wp,    2.7066956465_wp,    0.4122261241_wp,    0.1459709719_wp, &
      &    0.0401514465_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &   10.3452055112_wp,    4.6936355030_wp,    2.6553895791_wp,    1.1142653084_wp,    0.4349684397_wp, &
      &    0.1739785777_wp,    0.0585613195_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 93) = reshape([&
      &    6.6954006093_wp,    3.2428431471_wp,    1.4869401220_wp,    0.4679064994_wp,    0.4162408120_wp, &
      &    0.0563355612_wp,    0.0179742097_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.0847507016_wp,    1.1851609119_wp,    0.3046167128_wp,    0.0289692345_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    9.0641057249_wp,    5.9477714169_wp,    2.8388440382_wp,    0.4154086812_wp,    0.1465748848_wp, &
      &    0.0396111051_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &   12.1964726953_wp,    5.1557768789_wp,    2.7370716549_wp,    1.1266697717_wp,    0.4479703666_wp, &
      &    0.1990949718_wp,    0.0619953667_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 94) = reshape([&
      &    6.6303117355_wp,    3.5466077014_wp,    1.4430320263_wp,    0.4600369969_wp,    0.3943355086_wp, &
      &    0.0642228196_wp,    0.0161546839_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.9850558282_wp,    1.0920239440_wp,    0.2700788005_wp,    0.0316845111_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    9.2284044015_wp,    7.0601875208_wp,    2.9686897673_wp,    0.4330935943_wp,    0.1503418894_wp, &
      &    0.0395026063_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &   17.5531078430_wp,    6.1671571925_wp,    2.9224649551_wp,    1.1746242807_wp,    0.4718807083_wp, &
      &    0.2199544748_wp,    0.0589542836_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 95) = reshape([&
      &    7.0424500779_wp,    3.5775826451_wp,    1.5177274070_wp,    0.4741165153_wp,    0.4505584314_wp, &
      &    0.0568807786_wp,    0.0164432909_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.0349259246_wp,    1.0944107493_wp,    0.2806318835_wp,    0.0318669819_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    9.2533032610_wp,    7.1898312518_wp,    3.1576326210_wp,    0.4335400810_wp,    0.1422294620_wp, &
      &    0.0402329381_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &   22.0064234839_wp,    6.8863125525_wp,    3.0868737987_wp,    1.2777797360_wp,    0.5355416438_wp, &
      &    0.2421878809_wp,    0.0586301830_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 96) = reshape([&
      &    6.9341153400_wp,    3.6225082743_wp,    1.5047460716_wp,    0.4855796144_wp,    0.4501347952_wp, &
      &    0.0567110582_wp,    0.0157310729_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.0822959819_wp,    1.1238926455_wp,    0.2864063044_wp,    0.0328786335_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    9.4083679284_wp,    7.4920567095_wp,    3.2707904809_wp,    0.4488055445_wp,    0.1463427263_wp, &
      &    0.0393433753_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &   25.3267286214_wp,    7.7009138138_wp,    3.2367439792_wp,    1.3548370814_wp,    0.5751211772_wp, &
      &    0.2605448497_wp,    0.0643172804_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 97) = reshape([&
      &    7.4950989526_wp,    3.4243303793_wp,    1.5864713404_wp,    0.4882423579_wp,    0.4561501937_wp, &
      &    0.0539403961_wp,    0.0172116334_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.2243578919_wp,    1.1804603178_wp,    0.2868002591_wp,    0.0345143164_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    8.1889421318_wp,    8.1901765294_wp,    3.6408812774_wp,    0.4526306901_wp,    0.1473981013_wp, &
      &    0.0415959901_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &   26.9160948078_wp,    8.2487969208_wp,    3.4083409786_wp,    1.4232647472_wp,    0.5894958818_wp, &
      &    0.2413126753_wp,    0.0616604244_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 98) = reshape([&
      &    7.9459576053_wp,    3.8911918217_wp,    1.6964057794_wp,    0.5191278857_wp,    0.4846313963_wp, &
      &    0.0635981414_wp,    0.0212160316_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.0592884537_wp,    1.0948281728_wp,    0.3056870537_wp,    0.0359858600_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   10.2965357068_wp,   10.2961961727_wp,    3.7087473515_wp,    0.4405616856_wp,    0.1341671984_wp, &
      &    0.0491374299_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &   28.8778862843_wp,    9.0703146481_wp,    3.6261090658_wp,    1.5600373766_wp,    0.6607578887_wp, &
      &    0.2523710264_wp,    0.0612669222_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 99) = reshape([&
      &    7.6882429650_wp,    3.7813249201_wp,    1.8453477108_wp,    0.5286304942_wp,    0.4720639279_wp, &
      &    0.0637464013_wp,    0.0216822117_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.1640877745_wp,    1.0965302889_wp,    0.3123338101_wp,    0.0367928373_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   11.7491049968_wp,   10.6700703538_wp,    3.8351743478_wp,    0.4642600544_wp,    0.1440255055_wp, &
      &    0.0572976116_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &   33.0035392828_wp,   10.5555511750_wp,    4.0641153179_wp,    1.6832040660_wp,    0.6839240184_wp, &
      &    0.2559324600_wp,    0.0579029571_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 100) = reshape([&
      &    8.2369399121_wp,    4.2137141108_wp,    1.8327379111_wp,    0.5060988337_wp,    0.5048567637_wp, &
      &    0.0700665318_wp,    0.0210969004_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.8961806248_wp,    1.0335975210_wp,    0.3614338102_wp,    0.0351670299_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   11.2512607251_wp,   10.2242689301_wp,    4.1358305685_wp,    0.4635625882_wp,    0.1594717863_wp, &
      &    0.0612982855_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &   33.0321960994_wp,   10.7776628002_wp,    4.3009148057_wp,    1.8506523048_wp,    0.7553841824_wp, &
      &    0.2811112724_wp,    0.0663325373_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 101) = reshape([&
      &    8.0478911868_wp,    4.0957879345_wp,    1.6325931418_wp,    0.5636929285_wp,    0.4890392567_wp, &
      &    0.0636170692_wp,    0.0224854933_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.0762900087_wp,    0.9889693711_wp,    0.3247942138_wp,    0.0363462016_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   14.3137277144_wp,   11.3750666030_wp,    3.7778349185_wp,    0.4651190133_wp,    0.1716670765_wp, &
      &    0.0465668554_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &   34.8161829262_wp,   11.2183524455_wp,    4.3970350375_wp,    1.8038936592_wp,    0.7082227529_wp, &
      &    0.2507268825_wp,    0.0624128345_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 102) = reshape([&
      &    7.3156183110_wp,    4.0209428247_wp,    1.6967889710_wp,    0.5497342641_wp,    0.4944136800_wp, &
      &    0.0730857100_wp,    0.0232675455_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.3699233533_wp,    1.1483874489_wp,    0.3174529183_wp,    0.0378796116_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   11.8573676894_wp,   14.5440595400_wp,    4.4806031954_wp,    0.4952398557_wp,    0.1724814384_wp, &
      &    0.0494814122_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &   33.2484770982_wp,   10.7236714702_wp,    4.1467833485_wp,    1.7013913095_wp,    0.6701263677_wp, &
      &    0.2493630187_wp,    0.0652077156_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      exponents(:, :, 103) = reshape([&
      &    8.5722777227_wp,    4.3221758001_wp,    1.7369169223_wp,    0.6462505715_wp,    0.4942253720_wp, &
      &    0.0696099985_wp,    0.0195701020_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.0176396835_wp,    1.3587235670_wp,    0.2866424405_wp,    0.0335984325_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   12.9394430468_wp,    9.6274373969_wp,    4.1938099787_wp,    0.5072944338_wp,    0.1520666553_wp, &
      &    0.0543134141_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &   35.0059226210_wp,   10.9879237590_wp,    4.3700229591_wp,    1.8188377809_wp,    0.7245877475_wp, &
      &    0.2777464988_wp,    0.0517780716_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))

   end subroutine setCGTOexponents

   subroutine setCGTOcoefficients()
      ! set up the array of CGTO coefficients during initialization of the basis set data
      coefficients(:, :,  1) = reshape([&
      &    0.0008797951_wp,    0.0069174165_wp,    0.0342454085_wp,    0.1433864428_wp,    0.3867191850_wp, &
      &    0.5050461002_wp,    0.6518542832_wp,    0.3855926056_wp,    0.0000000000_wp, & ! shell 1
      &    0.2137741444_wp,    0.5182428239_wp,    0.8280851349_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :,  2) = reshape([&
      &    0.0024473520_wp,    0.0184503394_wp,    0.0909195777_wp,    0.2932289516_wp,    0.6018599857_wp, &
      &    0.6332182440_wp,    0.0273485662_wp,    0.3761054064_wp,    0.0000000000_wp, & ! shell 1
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :,  3) = reshape([&
      &    0.0091601628_wp,   -0.1692728979_wp,    0.0077909725_wp,    0.4086745672_wp,    0.8967648389_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0021234303_wp,    0.0358000079_wp,    0.0984515864_wp,    0.3765269279_wp,    0.9204610845_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :,  4) = reshape([&
      &    0.0111724591_wp,   -0.2500794922_wp,    0.0911334125_wp,    0.6150053486_wp,    0.7419423986_wp, &
      &    0.0178891935_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0488242092_wp,    0.2075209455_wp,    0.4973082342_wp,    0.8409731113_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :,  5) = reshape([&
      &   -0.5411167832_wp,    0.4187078610_wp,    0.5213985993_wp,    0.4785957335_wp,    0.1145214505_wp, &
      &    0.1336069456_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0212573319_wp,    0.1117723723_wp,    0.4404252790_wp,    0.5315229425_wp,    0.7145376112_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :,  6) = reshape([&
      &   -0.6278730641_wp,    0.5339759085_wp,    0.3846239012_wp,    0.3939820565_wp,    0.0312304588_wp, &
      &    0.1285005712_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0076955434_wp,    0.0650844755_wp,    0.2184596616_wp,    0.5573292664_wp,    0.5354300006_wp, &
      &    0.5921815341_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :,  7) = reshape([&
      &   -0.6295987180_wp,    0.5951877076_wp,    0.3977529686_wp,    0.2569623688_wp,    0.1371739153_wp, &
      &    0.0793932043_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0165125409_wp,    0.1143319210_wp,    0.3448816460_wp,    0.6526610499_wp,    0.6161896255_wp, &
      &    0.2491106129_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :,  8) = reshape([&
      &   -0.6120497222_wp,    0.6027699174_wp,    0.4391296589_wp,    0.2150047976_wp,    0.1281369272_wp, &
      &    0.0811330500_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0176755874_wp,    0.1167992881_wp,    0.3497875054_wp,    0.6129920371_wp,    0.6438595392_wp, &
      &    0.2708871669_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :,  9) = reshape([&
      &   -0.4831574684_wp,    0.3050162101_wp,    0.5473143900_wp,    0.3420574984_wp,    0.4805329668_wp, &
      &    0.1614176769_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0171530702_wp,    0.1002706538_wp,    0.3030250851_wp,    0.5779777167_wp,    0.6760861482_wp, &
      &    0.3266139085_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 10) = reshape([&
      &   -0.4457054459_wp,    0.2875953102_wp,    0.5036922891_wp,    0.6607219681_wp,    0.1547577125_wp, &
      &    0.0665297101_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0141913108_wp,    0.1019938420_wp,    0.3062604555_wp,    0.5839956635_wp,    0.6735526759_wp, &
      &    0.3176102235_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 11) = reshape([&
      &    0.0170781613_wp,   -0.2348190977_wp,    0.3659683683_wp,    0.9003529759_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0987902537_wp,    0.0831657441_wp,    0.2597219446_wp,    0.9570101652_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 12) = reshape([&
      &    0.0304936108_wp,   -0.2883216386_wp,    0.5475865395_wp,    0.7849138514_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0131634221_wp,    0.0883963755_wp,    0.9959983961_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0272443312_wp,    0.9996288043_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 13) = reshape([&
      &    0.0417237738_wp,   -0.3588648490_wp,    0.1065464072_wp,    0.6954331606_wp,    0.6119605617_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1233746959_wp,    0.2691220282_wp,    0.5531575486_wp,    0.7786968247_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.2458922077_wp,    0.9692971795_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 14) = reshape([&
      &    0.0414345795_wp,   -0.3963751030_wp,    0.4606821646_wp,    0.4086599922_wp,    0.6796608767_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1031873781_wp,    0.3415806338_wp,    0.7379107291_wp,    0.5728549480_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1780066147_wp,    0.9840292908_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 15) = reshape([&
      &    0.0284680112_wp,   -0.3899003269_wp,    0.4655149132_wp,    0.7553699769_wp,    0.2447026174_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1884168735_wp,    0.3742550489_wp,    0.6274297912_wp,    0.5523288210_wp,    0.3545376859_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.2348501761_wp,    0.9720315812_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 16) = reshape([&
      &    0.0406614047_wp,   -0.4025216708_wp,    0.5121117218_wp,    0.7362830574_wp,    0.1787506602_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.2566234456_wp,    0.3838628533_wp,    0.6347486028_wp,    0.5691389704_wp,    0.2448851988_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.2650215302_wp,    0.9642424947_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 17) = reshape([&
      &    0.0191286317_wp,   -0.3938152784_wp,    0.4816874711_wp,    0.7512864570_wp,    0.2192930953_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1454702085_wp,    0.1702339017_wp,    0.6008784632_wp,    0.6939222987_wp,    0.3275297741_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.2628355192_wp,    0.9648406552_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 18) = reshape([&
      &    0.0091546612_wp,   -0.6019941338_wp,    0.4419591263_wp,    0.6206929237_wp,    0.2386036045_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0816454155_wp,    0.4980643204_wp,    0.7556730066_wp,    0.2654866387_wp,    0.3220886689_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.3151363069_wp,    0.9490464204_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 19) = reshape([&
      &    0.0562318304_wp,   -0.3467397064_wp,    0.4194511813_wp,    0.8370604899_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0509016713_wp,    0.2393287825_wp,    0.9696034002_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 20) = reshape([&
      &    0.1471893494_wp,   -0.3559717629_wp,    0.9224088917_wp,   -0.0279505978_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0375806482_wp,    0.5278867414_wp,    0.8484829304_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.2639113036_wp,    0.6059555246_wp,    0.7504456849_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 21) = reshape([&
      &    0.1618472261_wp,   -0.4630145068_wp,    0.4380445203_wp,    0.7532900185_wp,    0.0097050631_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0856452797_wp,    0.9963256928_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0112965030_wp,    0.0851812013_wp,    0.2821877935_wp,    0.5230350798_wp,    0.5698808591_wp, &
      &    0.5609427002_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 22) = reshape([&
      &    0.0764328206_wp,   -0.3905229863_wp,    0.0197075952_wp,    0.4003435540_wp,    0.8252190440_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1193097356_wp,    0.3989070958_wp,    0.9091965222_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0127453115_wp,    0.0909051855_wp,    0.3043841753_wp,    0.5522076698_wp,    0.6055762230_wp, &
      &    0.4767265522_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 23) = reshape([&
      &    0.0986138931_wp,   -0.4174831051_wp,    0.0248020947_wp,    0.5951632066_wp,    0.6790793552_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1109066879_wp,    0.5671380496_wp,    0.8161213998_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0183547569_wp,    0.1099001008_wp,    0.3485793408_wp,    0.5934605985_wp,    0.6013053385_wp, &
      &    0.3902741625_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 24) = reshape([&
      &    0.0034415715_wp,    0.1272233059_wp,   -0.3879090551_wp,    0.5508882893_wp,    0.7279086781_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0761040328_wp,    0.4072806121_wp,    0.9101267380_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0169670556_wp,    0.1080759127_wp,    0.3376774166_wp,    0.5796710186_wp,    0.6285907811_wp, &
      &    0.3779693354_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 25) = reshape([&
      &    0.0525208479_wp,   -0.3305291098_wp,    0.0040389215_wp,    0.2910486993_wp,    0.8962513096_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0962486721_wp,    0.4190716793_wp,    0.9028372615_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0216733320_wp,    0.1238086560_wp,    0.3730430949_wp,    0.6203116725_wp,    0.5841462452_wp, &
      &    0.3450030812_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 26) = reshape([&
      &    0.0489400913_wp,   -0.3058988563_wp,    0.0033110529_wp,    0.8735574857_wp,    0.3753892823_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1328305457_wp,    0.6296558018_wp,    0.7654342672_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0227451276_wp,    0.1270304738_wp,    0.3748468573_wp,    0.6008341549_wp,    0.5945006501_wp, &
      &    0.3583337087_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 27) = reshape([&
      &    0.0589229501_wp,   -0.3302113218_wp,    0.0033842491_wp,    0.8345179531_wp,    0.4371005625_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0894253730_wp,    0.4062430642_wp,    0.9093787305_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0224379944_wp,    0.1216614429_wp,    0.3648798796_wp,    0.5962902663_wp,    0.5958693832_wp, &
      &    0.3754135049_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 28) = reshape([&
      &    0.1137751160_wp,   -0.3579580997_wp,    0.0145224213_wp,    0.8880433245_wp,    0.2647439800_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1505671349_wp,    0.7225177901_wp,    0.6747574237_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0244185080_wp,    0.1440804529_wp,    0.4115395674_wp,    0.6194290011_wp,    0.5814412107_wp, &
      &    0.2958269344_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 29) = reshape([&
      &    0.0988321253_wp,   -0.3566702806_wp,    0.0255824319_wp,    0.6199513440_wp,    0.6913930808_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1033026693_wp,    0.5422685471_wp,    0.8338305471_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0255800251_wp,    0.1394024928_wp,    0.3953184275_wp,    0.6034710176_wp,    0.5715106996_wp, &
      &    0.3644642635_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 30) = reshape([&
      &    0.1577424595_wp,   -0.4118626434_wp,    0.0176025730_wp,    0.7188377332_wp,    0.5370744289_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0524465131_wp,    0.4882709861_wp,    0.8711146924_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 31) = reshape([&
      &    0.1712628331_wp,   -0.4835181846_wp,    0.2251571753_wp,    0.5791940948_wp,    0.5922141961_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.3052156060_wp,    0.3005052051_wp,    0.4292243617_wp,    0.7951770261_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0140287134_wp,    0.9999015928_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 32) = reshape([&
      &    0.2372792977_wp,   -0.5611582741_wp,    0.2335290062_wp,    0.7037364758_wp,    0.1352620950_wp, &
      &    0.2464209158_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.3456179981_wp,    0.3238424031_wp,    0.6687411234_wp,    0.4566695252_wp,    0.3462839180_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.2820184263_wp,    0.9594089885_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 33) = reshape([&
      &    0.1810417528_wp,   -0.5452131708_wp,    0.3469962539_wp,    0.6874476641_wp,    0.0951799837_wp, &
      &    0.2606080615_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.4436672132_wp,    0.4025510111_wp,    0.4798742943_wp,    0.5518884119_wp,    0.3259630805_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.9232101079_wp,    0.3842955851_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 34) = reshape([&
      &    0.2312152963_wp,   -0.5881334406_wp,    0.6750007099_wp,    0.3376099979_wp,    0.1480183673_wp, &
      &    0.0955124950_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.3673706348_wp,    0.3693860165_wp,    0.5745986713_wp,    0.5576643771_wp,    0.2957018716_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.9559071832_wp,    0.2936689585_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 35) = reshape([&
      &    0.1401068844_wp,   -0.5307574947_wp,    0.6949480335_wp,    0.3833236240_wp,    0.2520879086_wp, &
      &    0.0723080839_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.3987787976_wp,    0.3516729392_wp,    0.5520983695_wp,    0.5772111044_wp,    0.2816315781_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.7209557045_wp,    0.6929811485_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 36) = reshape([&
      &    0.2904161981_wp,   -0.6120999197_wp,    0.6672857422_wp,    0.2866877459_wp,    0.1126491136_wp, &
      &    0.0290201986_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.3495470302_wp,    0.2943973288_wp,    0.5538464789_wp,    0.6319851534_wp,    0.2915406151_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.6005861208_wp,    0.7995600737_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 37) = reshape([&
      &    0.2603571353_wp,   -0.5301059633_wp,    0.5927627842_wp,    0.5475711018_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1610247858_wp,    0.1558798218_wp,    0.9745627222_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 38) = reshape([&
      &    0.2447165478_wp,   -0.5508472131_wp,    0.6094201822_wp,    0.5150613561_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1779667291_wp,    0.4484713057_wp,    0.8759002976_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1340774449_wp,    0.4035484803_wp,    0.9050811360_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 39) = reshape([&
      &    0.1384194123_wp,   -0.4891277340_wp,    0.3963786370_wp,    0.0282728699_wp,    0.7639887087_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1322882753_wp,    0.9912112854_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0113967902_wp,    0.0365703004_wp,    0.3006091191_wp,    0.5396601122_wp,    0.4861701613_wp, &
      &    0.6169055206_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 40) = reshape([&
      &    0.1189309415_wp,   -0.5469519868_wp,    0.0567130094_wp,    0.2775273461_wp,    0.7787561634_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1483428583_wp,    0.2162757064_wp,    0.9649970027_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0121875392_wp,    0.0687087357_wp,    0.4397634943_wp,    0.6762830027_wp,    0.4688148993_wp, &
      &    0.3529766749_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 41) = reshape([&
      &    0.0848701921_wp,   -0.4546407590_wp,    0.0681861191_wp,    0.3545036838_wp,    0.8098003594_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1499020311_wp,    0.5978665514_wp,    0.7874547402_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0153736775_wp,    0.1702474680_wp,    0.5076541403_wp,    0.7131624631_wp,    0.4101877874_wp, &
      &    0.1902945188_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 42) = reshape([&
      &    0.1587183711_wp,   -0.5780942295_wp,    0.6227009570_wp,    0.2841147887_wp,    0.4148949811_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.3961700779_wp,    0.4997813982_wp,    0.7702388093_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0179072480_wp,    0.0575847757_wp,    0.4877527996_wp,    0.6855583923_wp,    0.4696057199_wp, &
      &    0.2606543479_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 43) = reshape([&
      &    0.0457356544_wp,   -0.4312126915_wp,    0.3213055532_wp,    0.0503524486_wp,    0.8403518531_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0984630388_wp,    0.9951407086_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0152185716_wp,    0.0415236730_wp,    0.4443520315_wp,    0.6568013459_wp,    0.5256080046_wp, &
      &    0.3048666417_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 44) = reshape([&
      &    0.1158592386_wp,   -0.4438037092_wp,    0.5821157386_wp,    0.0435440404_wp,    0.6699702143_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1403563065_wp,    0.4757122924_wp,    0.8683305374_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0186687577_wp,    0.0527359954_wp,    0.4527321090_wp,    0.6348823665_wp,    0.5592273616_wp, &
      &    0.2758499021_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 45) = reshape([&
      &    0.1074080596_wp,   -0.4776741101_wp,    0.4879137086_wp,    0.1190235605_wp,    0.7127864745_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1693022496_wp,    0.3890084347_wp,    0.9055435859_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0184646658_wp,    0.0447801115_wp,    0.5096982829_wp,    0.6734074589_wp,    0.4917718145_wp, &
      &    0.2062627812_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 46) = reshape([&
      &    0.0978718293_wp,   -0.4603250007_wp,    0.6933366340_wp,    0.5173540324_wp,    0.1736407667_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.2040021873_wp,    0.7777688036_wp,    0.5945240077_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0154996026_wp,    0.0968231881_wp,    0.5121189480_wp,    0.6442696082_wp,    0.5135620659_wp, &
      &    0.2220132701_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 47) = reshape([&
      &    0.0498175129_wp,   -0.3774884095_wp,    0.5911548953_wp,    0.0667235644_wp,    0.7078874006_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.2356807166_wp,    0.6272503186_wp,    0.7423015813_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0194789342_wp,    0.0485863911_wp,    0.4528706162_wp,    0.6748473286_wp,    0.5405740718_wp, &
      &    0.2110187070_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 48) = reshape([&
      &    0.2314826607_wp,   -0.5292876368_wp,    0.0335698089_wp,    0.6422604123_wp,    0.5026380467_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0996010230_wp,    0.4852484727_wp,    0.8686849579_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 49) = reshape([&
      &    0.3700762240_wp,   -0.6985898373_wp,    0.1653350518_wp,    0.4449435170_wp,    0.3869178401_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.2480828373_wp,   -0.3625619255_wp,    0.3952502150_wp,    0.8067099997_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.2689821861_wp,    0.9631451519_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 50) = reshape([&
      &    0.2172795496_wp,   -0.6339959249_wp,    0.1215496176_wp,    0.5518275967_wp,    0.3616294061_wp, &
      &    0.3174506751_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1261549945_wp,   -0.2781739271_wp,    0.7364820624_wp,    0.4904828341_wp,    0.3517455684_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.5636094751_wp,    0.8260413788_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 51) = reshape([&
      &    0.4187333702_wp,   -0.7353974710_wp,    0.1536594472_wp,    0.2537511775_wp,    0.4197196119_wp, &
      &    0.1403121018_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1698396681_wp,   -0.3952911874_wp,    0.6034059180_wp,    0.5790399146_wp,    0.3398726816_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.9649954000_wp,    0.2622668067_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 52) = reshape([&
      &    0.3315713584_wp,   -0.6332758814_wp,    0.0974014827_wp,    0.6555907562_wp,    0.0347170932_wp, &
      &    0.2202964533_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1945171972_wp,   -0.4349644621_wp,    0.4625248850_wp,    0.6761709294_wp,    0.3191121772_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.9567206813_wp,    0.2910077970_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 53) = reshape([&
      &    0.3331637936_wp,   -0.7324861637_wp,    0.3658877200_wp,    0.2516379618_wp,    0.2798428301_wp, &
      &    0.2774137877_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.5436185166_wp,    0.6081942493_wp,    0.4877318153_wp,    0.1044121258_wp,    0.2929068928_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.8884690195_wp,    0.4589365985_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 54) = reshape([&
      &    0.3605569215_wp,   -0.7243674614_wp,    0.4382595924_wp,    0.0819559243_wp,    0.3758081028_wp, &
      &    0.0725982999_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.5823050794_wp,    0.6402162549_wp,    0.3185444982_wp,    0.3361960216_wp,    0.1911689805_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.5969336668_wp,    0.8022905942_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 55) = reshape([&
      &    0.5411165388_wp,   -0.7071909468_wp,    0.3889459502_wp,    0.2362094495_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.2409527693_wp,    0.1746738794_wp,    0.9546888492_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 56) = reshape([&
      &    0.3589947574_wp,   -0.5869148088_wp,    0.7225617684_wp,    0.0675149042_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.3210782681_wp,    0.6228060848_wp,    0.7134573053_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0382293382_wp,    0.8052274854_wp,    0.5917323842_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 57) = reshape([&
      &    0.1125351619_wp,   -0.6498934146_wp,    0.4677954202_wp,    0.5883381951_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.2841108709_wp,    0.9587914335_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0211639960_wp,    0.2417075911_wp,    0.3257269735_wp,    0.6986702453_wp,    0.5889748320_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 72) = reshape([&
      &    0.2107171529_wp,   -0.5534548310_wp,    0.5710553379_wp,    0.5684908377_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1976775682_wp,    0.3293759382_wp,    0.9232741036_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0148383090_wp,    0.3728906326_wp,    0.6463339435_wp,    0.5024996632_wp,    0.4364389105_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 73) = reshape([&
      &    0.2036816385_wp,   -0.7299149457_wp,    0.4310202489_wp,    0.3109789845_wp,    0.3784859025_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.2234018561_wp,    0.4055731707_wp,    0.8863419283_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0221840631_wp,    0.3984719577_wp,    0.6539247798_wp,    0.5150778153_wp,    0.3844544092_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 74) = reshape([&
      &    0.1388448140_wp,   -0.5437354269_wp,    0.5789600957_wp,    0.2508987918_wp,    0.5356574530_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1919177025_wp,    0.3886080406_wp,    0.9011944220_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0174373071_wp,    0.4855666476_wp,    0.7127397284_wp,    0.4414771791_wp,    0.2470241911_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 75) = reshape([&
      &    0.1296388571_wp,   -0.5083456833_wp,    0.5850318872_wp,    0.0781271639_wp,    0.6135244658_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.2239090087_wp,    0.2635198362_wp,    0.9383080794_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0243136225_wp,    0.0428773652_wp,    0.5417981842_wp,    0.7192672447_wp,    0.4320645062_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 76) = reshape([&
      &    0.2032560741_wp,   -0.5689765294_wp,    0.6674772081_wp,    0.3281236331_wp,    0.2859400907_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.3026180522_wp,    0.7667961264_wp,    0.5660795130_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0302778289_wp,    0.5292335161_wp,    0.6839314002_wp,    0.4357467889_wp,    0.2477048936_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 77) = reshape([&
      &    0.1606361546_wp,   -0.6437784939_wp,    0.4355343170_wp,    0.4803948704_wp,    0.3731968701_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.5758104822_wp,    0.6147568227_wp,    0.5389956749_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0219640049_wp,    0.5913224677_wp,    0.6563731283_wp,    0.4128253415_wp,    0.2204651346_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 78) = reshape([&
      &    0.1339882829_wp,   -0.4839138714_wp,    0.7579465194_wp,    0.1048846445_wp,    0.4029774066_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.2541777513_wp,    0.6404779603_wp,    0.7246941790_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0260539307_wp,    0.5904130377_wp,    0.6655962369_wp,    0.4225252599_wp,    0.1708440569_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 79) = reshape([&
      &    0.1221444769_wp,   -0.4885939937_wp,    0.7345265275_wp,    0.0602220356_wp,    0.4507779087_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.3129794168_wp,    0.6976403354_wp,    0.6444702065_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0272142229_wp,    0.6351737108_wp,    0.6348454866_wp,    0.3988885966_wp,    0.1835016046_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 80) = reshape([&
      &    0.2500087888_wp,   -0.5834572637_wp,    0.0347702407_wp,    0.6403978885_wp,    0.4309928093_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1442167476_wp,    0.6219844395_wp,    0.7696342552_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 81) = reshape([&
      &    0.3450301959_wp,   -0.6656616443_wp,    0.0357889952_wp,    0.5749495194_wp,    0.3255778512_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.4174483273_wp,   -0.5202101027_wp,    0.3707016436_wp,    0.6462960889_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0037544065_wp,    0.9999929522_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 82) = reshape([&
      &    0.3085798951_wp,   -0.6458378848_wp,    0.0682777796_wp,    0.6529758632_wp,    0.2379759270_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.3587852971_wp,   -0.5800636971_wp,    0.5451357598_wp,    0.4375998987_wp,    0.2147848923_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0413177536_wp,    0.9991460570_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 83) = reshape([&
      &    0.0823877467_wp,   -0.4796066120_wp,    0.1525140361_wp,    0.8206600476_wp,    0.1100906299_wp, &
      &    0.2330801689_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.4816856782_wp,   -0.6540503362_wp,    0.2961889083_wp,    0.3316365116_wp,    0.3774737341_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1346102166_wp,    0.9908986273_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 84) = reshape([&
      &    0.4010728273_wp,   -0.7295820208_wp,    0.0734539808_wp,    0.5150561810_wp,    0.1527558631_wp, &
      &    0.1133046846_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.2296779161_wp,   -0.4469864296_wp,    0.7053358031_wp,    0.2495638920_wp,    0.4332094820_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.8124765082_wp,    0.5829939310_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 85) = reshape([&
      &    0.3342127970_wp,   -0.7542585250_wp,    0.0631922385_wp,    0.2998744672_wp,    0.4275126700_wp, &
      &    0.2066660246_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1087412670_wp,   -0.3320658662_wp,    0.6044866700_wp,    0.6556670498_wp,    0.2874094343_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.3736695688_wp,    0.9275618865_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 86) = reshape([&
      &    0.3939436080_wp,   -0.7320742471_wp,    0.5004670761_wp,    0.0427547759_wp,    0.2287297851_wp, &
      &    0.0652928071_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.2950890441_wp,   -0.4731990988_wp,    0.5394535997_wp,    0.6120232142_wp,    0.1530440069_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0175178810_wp,    0.9998465501_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 87) = reshape([&
      &    0.0178989118_wp,   -0.0423135281_wp,   -0.0150291546_wp,    0.2957721758_wp,   -0.3065576989_wp, &
      &   -0.4889117998_wp,    0.6486340148_wp,    0.3920393113_wp,    0.0524295837_wp, & ! shell 1
      &   -0.0221671923_wp,    0.0549264751_wp,    0.0353287876_wp,   -0.2729886071_wp,   -0.2844370069_wp, &
      &    0.1902142931_wp,    0.6432490965_wp,    0.4924666212_wp,    0.3838517080_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 88) = reshape([&
      &    0.0162362907_wp,   -0.0424443307_wp,   -0.0150306095_wp,    0.2644472873_wp,   -0.2894595949_wp, &
      &   -0.4719372908_wp,    0.7308575055_wp,    0.2892687206_wp,    0.0586159110_wp, & ! shell 1
      &   -0.0233286348_wp,    0.0430738191_wp,    0.0402972925_wp,   -0.2749420671_wp,   -0.2627989463_wp, &
      &    0.2222870690_wp,    0.6769680758_wp,    0.5275055017_wp,    0.2556574888_wp, & ! shell 2
      &   -0.0352668026_wp,    0.0487468206_wp,    0.6005959367_wp,    0.7972857210_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 89) = reshape([&
      &   -0.1253340463_wp,    0.1637248631_wp,    0.2566758480_wp,   -0.4865790976_wp,   -0.1916033991_wp, &
      &    0.7616252719_wp,    0.1950869540_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1418819586_wp,   -0.0760306377_wp,   -0.2230446704_wp,    0.9614259862_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0248986904_wp,   -0.0302084258_wp,   -0.3261056700_wp,    0.4749458656_wp,    0.7563511507_wp, &
      &    0.3073791791_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0172280161_wp,    0.1042108137_wp,    0.3721203780_wp,    0.6509265151_wp,    0.5902743270_wp, &
      &    0.2263085933_wp,    0.1643929363_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 90) = reshape([&
      &   -0.1271959559_wp,    0.1624899021_wp,    0.2503222995_wp,   -0.4742204353_wp,   -0.1883289044_wp, &
      &    0.7718515856_wp,    0.1965942500_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1347030695_wp,   -0.0792796506_wp,   -0.2063222462_wp,    0.9659197435_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0233708730_wp,   -0.0304109981_wp,   -0.3170288992_wp,    0.4739281137_wp,    0.7819880014_wp, &
      &    0.2488143067_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0195789916_wp,    0.1096308449_wp,    0.4042263904_wp,    0.6852199716_wp,    0.5306655023_wp, &
      &    0.2255237594_wp,    0.1490151518_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 91) = reshape([&
      &   -0.1201086826_wp,    0.1582956653_wp,    0.2407406610_wp,   -0.4725522158_wp,   -0.1832968678_wp, &
      &    0.7699991523_wp,    0.2296917236_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1285352812_wp,   -0.0759160038_wp,   -0.1987257572_wp,    0.9686193862_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0238798731_wp,   -0.0376296914_wp,   -0.3119634720_wp,    0.5360535297_wp,    0.7398003262_wp, &
      &    0.2569720622_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0258489656_wp,    0.1258415409_wp,    0.4369515250_wp,    0.7044959938_wp,    0.4434869084_wp, &
      &    0.2669452435_wp,    0.1682679296_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 92) = reshape([&
      &   -0.1154449394_wp,    0.1561529001_wp,    0.2341880440_wp,   -0.4712991452_wp,   -0.1864964316_wp, &
      &    0.7752364378_wp,    0.2225968552_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1367713544_wp,   -0.0871627823_wp,   -0.2088745728_wp,    0.9644001549_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0274479126_wp,   -0.0450605004_wp,   -0.2965094925_wp,    0.5567861393_wp,    0.7294289332_wp, &
      &    0.2592699572_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0281587436_wp,    0.1152420189_wp,    0.4762275503_wp,    0.6816945635_wp,    0.4850543969_wp, &
      &    0.2141075582_wp,    0.1153533280_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 93) = reshape([&
      &   -0.1171539572_wp,    0.1635394203_wp,    0.2268364114_wp,   -0.4645611425_wp,   -0.1905889424_wp, &
      &    0.7838510446_wp,    0.2037434430_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1312043245_wp,   -0.0876499466_wp,   -0.2027194646_wp,    0.9664407539_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0247805768_wp,   -0.0437949937_wp,   -0.2847188278_wp,    0.5594652805_wp,    0.7262433205_wp, &
      &    0.2756308219_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0272770922_wp,    0.1247156575_wp,    0.5303855607_wp,    0.6885869104_wp,    0.4274681297_wp, &
      &    0.1895222940_wp,    0.0979463473_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 94) = reshape([&
      &   -0.1184660502_wp,    0.1505210117_wp,    0.2357347072_wp,   -0.4663553720_wp,   -0.2087050340_wp, &
      &    0.7755352212_wp,    0.2126931223_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1342213033_wp,   -0.1066045101_wp,   -0.2073565314_wp,    0.9631320725_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0204732522_wp,   -0.0378899898_wp,   -0.2819129338_wp,    0.5455427721_wp,    0.7447448711_wp, &
      &    0.2576983762_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0217368663_wp,    0.1262957006_wp,    0.5706234577_wp,    0.6824438383_wp,    0.3732167204_wp, &
      &    0.2055556090_wp,    0.1034038347_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 95) = reshape([&
      &   -0.1127966078_wp,    0.1524067130_wp,    0.2230550580_wp,   -0.4460934313_wp,   -0.1968449638_wp, &
      &    0.7930260695_wp,    0.2183069493_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1318341657_wp,   -0.1036842585_wp,   -0.2076055772_wp,    0.9637267515_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0167462522_wp,   -0.0382076921_wp,   -0.2579374391_wp,    0.5452679752_wp,    0.7144572232_wp, &
      &    0.3520819850_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0229441505_wp,    0.1502684738_wp,    0.5725531245_wp,    0.6645469271_wp,    0.3876731132_wp, &
      &    0.2151971918_wp,    0.1041776318_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 96) = reshape([&
      &   -0.1103081023_wp,    0.1491186625_wp,    0.2227492234_wp,   -0.4436675497_wp,   -0.1927149771_wp, &
      &    0.7951945061_wp,    0.2228548204_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1288081152_wp,   -0.1009165960_wp,   -0.2077541237_wp,    0.9643974980_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0171831016_wp,   -0.0408839830_wp,   -0.2603318234_wp,    0.5646992312_wp,    0.7404758971_wp, &
      &    0.2511390188_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0249812993_wp,    0.1728181057_wp,    0.6061573470_wp,    0.6491470857_wp,    0.3762714905_wp, &
      &    0.1888298305_wp,    0.0587726933_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 97) = reshape([&
      &   -0.1118690072_wp,    0.1390448349_wp,    0.2121985768_wp,   -0.4327759621_wp,   -0.1861908769_wp, &
      &    0.8096215451_wp,    0.2137159350_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1313120712_wp,   -0.1039668197_wp,   -0.2139007161_wp,    0.9623900062_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0170344616_wp,   -0.0403887212_wp,   -0.2665912529_wp,    0.6157543596_wp,    0.6960972003_wp, &
      &    0.2516007524_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0266580037_wp,    0.1765737619_wp,    0.6157055252_wp,    0.6158124477_wp,    0.4033622159_wp, &
      &    0.1969820438_wp,    0.0910482851_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 98) = reshape([&
      &   -0.1251102536_wp,    0.1954843664_wp,    0.1804607524_wp,   -0.4150186872_wp,   -0.2294540519_wp, &
      &    0.8011996173_wp,    0.2162329863_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1501288098_wp,   -0.1119263303_wp,   -0.2231557617_wp,    0.9566270658_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0198355540_wp,   -0.0457969396_wp,   -0.2435963860_wp,    0.6548291670_wp,    0.6237374312_wp, &
      &    0.3468722687_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0419142658_wp,    0.2289988744_wp,    0.6118367966_wp,    0.5925530078_wp,    0.4244717774_wp, &
      &    0.1895547241_wp,    0.0650545450_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 99) = reshape([&
      &   -0.1202131436_wp,    0.1886627771_wp,    0.1633458020_wp,   -0.4007009752_wp,   -0.2254773337_wp, &
      &    0.8196301302_wp,    0.2001960432_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1608773420_wp,   -0.1157811318_wp,   -0.2246561397_wp,    0.9540664700_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0225197562_wp,   -0.0483209189_wp,   -0.2591300806_wp,    0.6868952957_wp,    0.6093959601_wp, &
      &    0.2946539789_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0433487378_wp,    0.2273558396_wp,    0.6010934208_wp,    0.6133426625_wp,    0.4132144075_wp, &
      &    0.1838531634_wp,    0.0661782045_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 100) = reshape([&
      &   -0.1167654967_wp,    0.1713728238_wp,    0.1726437930_wp,   -0.4365525385_wp,   -0.1973156516_wp, &
      &    0.8050598261_wp,    0.2226171284_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1561597368_wp,   -0.0982804490_wp,   -0.2060440453_wp,    0.9609895636_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0239733827_wp,   -0.0496517477_wp,   -0.2569215824_wp,    0.7084997881_wp,    0.6169511589_wp, &
      &    0.2198876967_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0496045267_wp,    0.2410740708_wp,    0.5743988457_wp,    0.6112091951_wp,    0.4375409407_wp, &
      &    0.1965266667_wp,    0.0764667489_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 101) = reshape([&
      &   -0.1158350861_wp,    0.1802826113_wp,    0.1828745644_wp,   -0.4015689892_wp,   -0.2295658366_wp, &
      &    0.8190809311_wp,    0.1891708445_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1501592210_wp,   -0.1072366072_wp,   -0.2223289995_wp,    0.9573517297_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0211398066_wp,   -0.0497702692_wp,   -0.2318346713_wp,    0.6422054989_wp,    0.6560811164_wp, &
      &    0.3169516997_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0462415661_wp,    0.2303851769_wp,    0.6154894119_wp,    0.6196769472_wp,    0.3960366263_wp, &
      &    0.1506456720_wp,    0.0491784994_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 102) = reshape([&
      &   -0.1220124651_wp,    0.1755904250_wp,    0.1891745041_wp,   -0.3934160539_wp,   -0.2687997970_wp, &
      &    0.7998942687_wp,    0.2272302896_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1523117064_wp,   -0.1224604975_wp,   -0.2007160315_wp,    0.9599571060_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0198549835_wp,   -0.0444862062_wp,   -0.2457650683_wp,    0.5920108562_wp,    0.7327400857_wp, &
      &    0.2232518784_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0626911809_wp,    0.2884228542_wp,    0.6510278131_wp,    0.5831744311_wp,    0.3664395745_wp, &
      &    0.1112809278_wp,    0.0478647664_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients(:, :, 103) = reshape([&
      &   -0.1147823174_wp,    0.1861682623_wp,    0.1725037079_wp,   -0.3837435952_wp,   -0.2376802820_wp, &
      &    0.8295832754_wp,    0.1744975486_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1181853768_wp,   -0.0995774642_wp,   -0.1874940767_wp,    0.9700322245_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0207984776_wp,   -0.0548218015_wp,   -0.2360952679_wp,    0.7071312041_wp,    0.6206594542_wp, &
      &    0.2357293362_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0459821308_wp,    0.2410759170_wp,    0.6490170772_wp,    0.6182994441_wp,    0.3564383574_wp, &
      &    0.0916093632_wp,    0.0284622344_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))

   end subroutine setCGTOcoefficients

   subroutine setCGTOcoefficients_environment()
      ! set up the array of CGTO coeffients of environment dependence during initialization of the basis set data
      coefficients_env(:, :,  1) = reshape([&
      &    0.0005296397_wp,    0.0011495460_wp,    0.0162334608_wp,    0.0420068739_wp,    0.1975533619_wp, &
      &    0.0002878381_wp,    0.0001587756_wp,   -0.8711386555_wp,    0.0000000000_wp, & ! shell 1
      &    0.1201006344_wp,    0.0102506589_wp,   -1.9886029516_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :,  2) = reshape([&
      &    0.0001463287_wp,    0.0003609613_wp,    0.0003145524_wp,    0.0024513496_wp,    0.0386169855_wp, &
      &    0.0009688548_wp,   -0.4350280799_wp,    0.0010449002_wp,    0.0000000000_wp, & ! shell 1
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :,  3) = reshape([&
      &   -0.0000308921_wp,   -0.1287949117_wp,    0.1502773795_wp,    0.0022913851_wp,   -0.5840919456_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0011614373_wp,    0.0035213134_wp,    0.0677562968_wp,   -0.0011215077_wp,   -0.5167236620_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :,  4) = reshape([&
      &    0.0003121112_wp,   -0.0011918362_wp,    0.0006086080_wp,    0.0004522728_wp,   -0.2965863324_wp, &
      &    0.0004637573_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0158814159_wp,    0.0908916665_wp,    0.0644513740_wp,   -0.3034802415_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :,  5) = reshape([&
      &    0.0009785608_wp,    0.0019415360_wp,    0.0053340156_wp,    0.0130127277_wp,   -0.3922652287_wp, &
      &    0.0130618843_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0009965729_wp,    0.0120956207_wp,    0.0018039597_wp,    0.0020547870_wp,   -0.3601246642_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :,  6) = reshape([&
      &    0.0005347158_wp,    0.0003901904_wp,    0.0003412173_wp,    0.0002301530_wp,   -0.2286195285_wp, &
      &    0.0001344651_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0008374230_wp,    0.0005178058_wp,    0.0032892179_wp,    0.0008561450_wp,    0.0005591193_wp, &
      &   -0.3153886518_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :,  7) = reshape([&
      &    0.0003605272_wp,    0.0003930935_wp,    0.0003773702_wp,    0.0001870998_wp,    0.0003770533_wp, &
      &   -0.1031510244_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0014424263_wp,    0.0021086746_wp,    0.0268142354_wp,    0.0203970898_wp,    0.0008521880_wp, &
      &   -0.3492299852_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :,  8) = reshape([&
      &   -0.0011087353_wp,    0.0038143623_wp,    0.0089804075_wp,    0.0012623104_wp,    0.0010086400_wp, &
      &   -0.0648728802_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0011757148_wp,    0.0061148138_wp,    0.0344819432_wp,    0.0365926222_wp,    0.0012780482_wp, &
      &   -0.2713024898_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :,  9) = reshape([&
      &   -0.0036420322_wp,    0.0068572925_wp,    0.0069031964_wp,    0.0028950497_wp,    0.0032597791_wp, &
      &   -0.1066346371_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0000625381_wp,    0.0005578566_wp,    0.0020823764_wp,    0.0029238286_wp,    0.0018456845_wp, &
      &   -0.2378848787_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 10) = reshape([&
      &   -0.0007130439_wp,    0.0013624745_wp,    0.0011207564_wp,    0.0014306716_wp,    0.0013718809_wp, &
      &   -0.0501177337_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0004590389_wp,    0.0012217305_wp,    0.0010442677_wp,    0.0012027281_wp,    0.0012419841_wp, &
      &   -0.1166055260_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 11) = reshape([&
      &    0.0000427255_wp,    0.0000397207_wp,    0.0014142536_wp,   -0.7676018943_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0014099993_wp,    0.0009304416_wp,    0.0054038469_wp,   -0.8107170305_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 12) = reshape([&
      &    0.0007888795_wp,   -0.0024058118_wp,    0.0023034910_wp,   -0.1507191623_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0244469178_wp,    0.1670794669_wp,   -0.1816535280_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1014213309_wp,   -0.1210816486_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 13) = reshape([&
      &    0.0009206259_wp,   -0.0008776361_wp,    0.0005959114_wp,    0.0015292523_wp,   -0.2189931705_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0005224137_wp,    0.0113022684_wp,    0.0035424305_wp,   -0.3303802887_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0158475649_wp,   -0.2504303108_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 14) = reshape([&
      &    0.0004084319_wp,   -0.0005843188_wp,    0.0031367318_wp,    0.0023741221_wp,   -0.2607070110_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0196185543_wp,    0.0716970759_wp,    0.0005283833_wp,   -0.2985639246_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1221119198_wp,   -0.0991085712_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 15) = reshape([&
      &    0.0004882180_wp,   -0.0004801835_wp,    0.0067676303_wp,    0.0018923361_wp,   -0.0774915115_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0023506293_wp,    0.0113380048_wp,    0.0017086639_wp,    0.0013268832_wp,   -0.3756159050_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0307186620_wp,   -0.0422614891_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 16) = reshape([&
      &    0.0006110785_wp,   -0.0006525088_wp,    0.0297470805_wp,    0.0019610002_wp,   -0.0885706253_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0029888373_wp,    0.0081919372_wp,    0.0029902686_wp,    0.0036224785_wp,   -0.2660290621_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0242455324_wp,   -0.0356259808_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 17) = reshape([&
      &    0.0003530395_wp,   -0.0001206536_wp,    0.0215370869_wp,    0.0014785576_wp,   -0.0870009673_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0009181752_wp,    0.0029776780_wp,    0.0008041102_wp,    0.0008543416_wp,   -0.2324067935_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0680146527_wp,   -0.0228949342_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 18) = reshape([&
      &    0.0004137385_wp,   -0.0005413185_wp,    0.0105596441_wp,    0.0010584302_wp,   -0.0550725585_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0008487849_wp,    0.0012454450_wp,    0.0011239510_wp,    0.0013147683_wp,   -0.1862978145_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0119716188_wp,   -0.0129172210_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 19) = reshape([&
      &   -0.0009205683_wp,   -0.0020460109_wp,   -0.0028254187_wp,   -0.3889172057_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1039105928_wp,    0.0286464904_wp,   -0.4243630718_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 20) = reshape([&
      &   -0.0000215642_wp,   -0.0007477592_wp,    0.0012975488_wp,   -0.1123073043_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0586897624_wp,    0.0597825501_wp,   -0.2224256517_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0189661376_wp,    0.0018102425_wp,    0.0008131504_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 21) = reshape([&
      &    0.0000544327_wp,   -0.0061515235_wp,    0.0018515458_wp,   -0.0508925783_wp,    0.0004446958_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0003424361_wp,   -0.0470334866_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000238501_wp,    0.0000527551_wp,    0.0006151565_wp,    0.0009012069_wp,    0.0136747019_wp, &
      &   -0.0105689799_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 22) = reshape([&
      &    0.0003771154_wp,   -0.0006700780_wp,    0.0001974796_wp,    0.0009151255_wp,   -0.1231539296_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0039050719_wp,    0.0148118595_wp,   -0.1396217238_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0001450231_wp,    0.0001353201_wp,    0.0000058656_wp,    0.0001917173_wp,    0.0182088541_wp, &
      &   -0.0167708810_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 23) = reshape([&
      &    0.0009310154_wp,   -0.0008877679_wp,    0.0002608874_wp,    0.0010264592_wp,   -0.2799869060_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0031944621_wp,    0.0296057709_wp,   -0.3776604222_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0002087508_wp,    0.0008359028_wp,    0.0006552667_wp,    0.0004175356_wp,    0.0240272887_wp, &
      &   -0.0068312796_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 24) = reshape([&
      &    0.0004305312_wp,    0.0010128571_wp,   -0.0009567035_wp,    0.0015951083_wp,   -0.3441085685_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0061058620_wp,    0.0017313234_wp,   -0.4019474311_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0002769565_wp,    0.0010072349_wp,    0.0002240819_wp,    0.0002411219_wp,    0.0010566318_wp, &
      &   -0.0013964950_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 25) = reshape([&
      &    0.0015024156_wp,   -0.0000898389_wp,    0.0003275977_wp,    0.0046530206_wp,   -0.1348924682_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0063030483_wp,    0.0103494254_wp,   -0.1374123062_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0001116721_wp,    0.0003996802_wp,    0.0006458858_wp,    0.0000780123_wp,    0.0079831970_wp, &
      &   -0.0009749076_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 26) = reshape([&
      &    0.0013400518_wp,   -0.0004648697_wp,    0.0004110043_wp,    0.0009202207_wp,   -0.5446634320_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0040216111_wp,    0.0460961733_wp,   -1.1552547006_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0001381274_wp,    0.0017329564_wp,    0.0006867330_wp,    0.0002501617_wp,    0.0006281418_wp, &
      &   -0.0656833803_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 27) = reshape([&
      &    0.0009895233_wp,   -0.0004962635_wp,    0.0004134322_wp,    0.0009706105_wp,   -0.6428988110_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0023440373_wp,    0.0094926023_wp,   -1.4164075983_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0002158402_wp,    0.0018556607_wp,    0.0006510962_wp,    0.0001651919_wp,    0.0003237446_wp, &
      &   -0.0545985192_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 28) = reshape([&
      &    0.0032358915_wp,   -0.0002021844_wp,    0.0003473983_wp,    0.0013866109_wp,   -0.1820396405_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0030068434_wp,    0.0196579775_wp,   -0.9933409226_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0002290243_wp,    0.0022505250_wp,    0.0095130421_wp,    0.0003829713_wp,    0.0002461333_wp, &
      &   -0.1410873163_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 29) = reshape([&
      &    0.0011724002_wp,   -0.0002247600_wp,    0.0003779710_wp,    0.0007705077_wp,   -0.2663780052_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0034426994_wp,    0.0132749369_wp,   -0.2857756649_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0001829206_wp,    0.0004245731_wp,    0.0002275513_wp,    0.0001487139_wp,   -0.0000479724_wp, &
      &   -0.0091410272_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 30) = reshape([&
      &    0.0013135971_wp,   -0.0015858861_wp,    0.0034039288_wp,    0.0031375738_wp,   -0.3100080159_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0040580946_wp,    0.0034360129_wp,   -0.5993163106_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 31) = reshape([&
      &    0.0022051937_wp,   -0.0035095060_wp,    0.0033226895_wp,    0.0019603321_wp,   -0.0898320152_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0043921896_wp,    0.0033877106_wp,    0.0032614110_wp,   -0.1396958883_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0606106588_wp,   -0.1385195222_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 32) = reshape([&
      &    0.0031723530_wp,   -0.0033190404_wp,    0.0051513746_wp,    0.0049578395_wp,   -0.0925287793_wp, &
      &    0.0026184370_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0025674462_wp,    0.0030357076_wp,    0.0051065750_wp,    0.0047336515_wp,   -0.2303793879_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0358224576_wp,   -0.2399984054_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 33) = reshape([&
      &    0.0025276521_wp,   -0.0023390676_wp,    0.0175474970_wp,    0.0027901575_wp,    0.0020374128_wp, &
      &   -0.0893688529_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0036488646_wp,    0.0061018385_wp,    0.0036031838_wp,    0.0024996846_wp,   -0.3009938385_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1959264592_wp,   -0.3300298590_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 34) = reshape([&
      &    0.0050234116_wp,   -0.0018652523_wp,    0.0163801844_wp,    0.0020474025_wp,    0.0027112856_wp, &
      &   -0.0748000427_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0017847325_wp,    0.0032226151_wp,    0.0023297686_wp,    0.0018711005_wp,   -0.2849597396_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0313016999_wp,   -0.0689272916_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 35) = reshape([&
      &    0.0045456875_wp,   -0.0009199939_wp,    0.0026451754_wp,    0.0040367671_wp,    0.0025345510_wp, &
      &   -0.0642103484_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0036374844_wp,    0.0070464130_wp,    0.0041444968_wp,    0.0034793975_wp,   -0.2022546674_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0433391346_wp,   -0.2129075600_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 36) = reshape([&
      &    0.0034638346_wp,   -0.0008422280_wp,    0.0050583771_wp,    0.0019188680_wp,    0.0018829865_wp, &
      &    0.0036110515_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0033198339_wp,    0.0051964196_wp,    0.0037553061_wp,    0.0033068070_wp,   -0.1327384530_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0170715113_wp,   -0.0167830443_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 37) = reshape([&
      &    0.0022496224_wp,   -0.0087561463_wp,    0.0036820185_wp,   -0.3405486473_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0006726147_wp,    0.0011974231_wp,   -0.3531512625_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 38) = reshape([&
      &    0.0008477202_wp,    0.0000504519_wp,    0.0003282471_wp,   -0.1097988678_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0429249498_wp,    0.0012453216_wp,   -0.1577332687_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0152016671_wp,    0.0447623889_wp,   -0.0720897872_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 39) = reshape([&
      &    0.0003038564_wp,   -0.0047463892_wp,    0.0003872298_wp,    0.0000735049_wp,   -0.0419966631_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0006719365_wp,   -0.0421837851_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0000087130_wp,    0.0000939915_wp,    0.0002571519_wp,    0.0041686176_wp,    0.0007523279_wp, &
      &   -0.0291093562_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 40) = reshape([&
      &    0.0010441424_wp,   -0.0006631136_wp,    0.0363951001_wp,    0.0009832517_wp,   -0.2762037480_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0017448344_wp,    0.0039768741_wp,   -0.3413606974_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0004103659_wp,    0.0004333441_wp,    0.0014956117_wp,    0.0005993425_wp,    0.0009943705_wp, &
      &   -0.1268008050_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 41) = reshape([&
      &    0.0005308735_wp,   -0.0011540738_wp,    0.0158160467_wp,    0.0008275638_wp,   -0.1578543101_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0012070677_wp,    0.0068440652_wp,   -0.1894033822_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0002036216_wp,    0.0001147754_wp,    0.0005343970_wp,    0.0002915832_wp,    0.0023186355_wp, &
      &   -0.0356150757_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 42) = reshape([&
      &    0.0001502833_wp,   -0.0001130311_wp,    0.0001746929_wp,    0.0001830460_wp,   -0.2833206445_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0089842860_wp,    0.0024281721_wp,   -0.2958282766_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0003647661_wp,    0.0001780550_wp,    0.0008172559_wp,    0.0010180277_wp,    0.0223313926_wp, &
      &   -0.0459800843_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 43) = reshape([&
      &    0.0000710099_wp,   -0.0005104263_wp,   -0.0000135646_wp,    0.0000582057_wp,   -0.0479220970_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0009986806_wp,   -0.0485849880_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0000619642_wp,    0.0002048716_wp,    0.0002900506_wp,    0.0001231811_wp,    0.0067818808_wp, &
      &   -0.0005949755_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 44) = reshape([&
      &    0.0000574979_wp,   -0.0001348331_wp,   -0.0000189219_wp,    0.0000568562_wp,   -0.2534784357_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0202789199_wp,    0.0208750307_wp,   -0.3311602261_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0000798897_wp,    0.0002223137_wp,    0.0009430127_wp,    0.0004107685_wp,    0.0125471178_wp, &
      &   -0.0218849787_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 45) = reshape([&
      &    0.0000906490_wp,   -0.0000624158_wp,    0.0000122100_wp,    0.0001182668_wp,   -0.0858560101_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0016298275_wp,    0.0068406418_wp,   -0.0994133122_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0000398754_wp,    0.0001990758_wp,    0.0015327130_wp,    0.0002060990_wp,    0.0103828112_wp, &
      &   -0.0020362664_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 46) = reshape([&
      &    0.0000527210_wp,   -0.0001033510_wp,    0.0000135561_wp,    0.0000989938_wp,   -0.3085709591_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0023376635_wp,    0.0073169814_wp,   -0.6626483400_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0001404756_wp,    0.0003277090_wp,    0.0259307113_wp,    0.0012055677_wp,    0.0005430567_wp, &
      &   -0.1245502992_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 47) = reshape([&
      &    0.0000739444_wp,   -0.0004923076_wp,    0.0000234551_wp,    0.0001557065_wp,   -0.2010614295_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0018066595_wp,    0.0133579525_wp,   -0.2001656975_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0000539089_wp,    0.0002180303_wp,    0.0006100929_wp,    0.0005500497_wp,    0.0011373162_wp, &
      &   -0.0032119982_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 48) = reshape([&
      &    0.0015274143_wp,   -0.0025939651_wp,    0.0057270477_wp,    0.0029348872_wp,   -0.1669149182_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0084534061_wp,    0.0033232862_wp,   -0.3245067452_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 49) = reshape([&
      &    0.0044484068_wp,   -0.0036391498_wp,    0.0062809043_wp,    0.0037431621_wp,   -0.1061104066_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0059695649_wp,   -0.0078006601_wp,    0.0031684847_wp,   -0.2265699246_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1588047763_wp,   -0.3792968058_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 50) = reshape([&
      &    0.0027773983_wp,   -0.0030385425_wp,    0.0040284764_wp,    0.0031375116_wp,    0.0027420115_wp, &
      &   -0.1042329224_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0020592548_wp,   -0.0026780332_wp,    0.0060137181_wp,    0.0028744514_wp,   -0.2763614607_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0881269712_wp,   -0.4062612750_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 51) = reshape([&
      &    0.0020890584_wp,   -0.0046038674_wp,    0.0026195413_wp,    0.0026763730_wp,    0.0027990623_wp, &
      &    0.0015960339_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0027896760_wp,   -0.0043782723_wp,    0.0043537546_wp,    0.0037499263_wp,   -0.3331974704_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0473395676_wp,   -0.0657653866_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 52) = reshape([&
      &    0.0009824975_wp,   -0.0008284961_wp,    0.0159773517_wp,    0.0009284976_wp,    0.0012947587_wp, &
      &   -0.0370891572_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0039153061_wp,   -0.0033856183_wp,    0.0045407903_wp,    0.0036147269_wp,   -0.2391790969_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0773960085_wp,   -0.2627191151_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 53) = reshape([&
      &    0.0027514173_wp,   -0.0014007840_wp,    0.0157651118_wp,    0.0021833961_wp,    0.0023100528_wp, &
      &   -0.0428177862_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0015704289_wp,    0.0078364870_wp,    0.0053931164_wp,    0.0039582799_wp,   -0.1579812994_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0843796435_wp,   -0.2597986847_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 54) = reshape([&
      &    0.0046390594_wp,   -0.0008206744_wp,    0.0022257414_wp,    0.0051542292_wp,    0.0021118672_wp, &
      &   -0.0349658491_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0027642590_wp,    0.0050038498_wp,    0.0046415929_wp,    0.0034092191_wp,   -0.1044984354_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0747489802_wp,   -0.4863492333_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 55) = reshape([&
      &    0.0051671294_wp,   -0.0041640485_wp,    0.0052322062_wp,   -0.1377106953_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0306044799_wp,    0.0013541063_wp,   -0.1860792254_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 56) = reshape([&
      &    0.0023161273_wp,    0.0010648131_wp,    0.0010688523_wp,   -0.0969864033_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0880410491_wp,    0.0003151140_wp,   -0.1380726676_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1339023617_wp,    0.0020952535_wp,   -0.1035128370_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 57) = reshape([&
      &    0.0016242310_wp,   -0.0048713040_wp,    0.0000300141_wp,   -0.0304020227_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0003696334_wp,   -0.0425083027_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0000429983_wp,    0.0002485657_wp,    0.0218650813_wp,    0.0002021293_wp,   -0.0272314524_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 72) = reshape([&
      &    0.0000377931_wp,   -0.0002660731_wp,    0.0000104570_wp,   -0.2024776748_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0006692014_wp,    0.0036449384_wp,   -0.3176740300_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0001247106_wp,    0.0000989576_wp,    0.0011692253_wp,    0.0008537226_wp,   -0.1586467347_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 73) = reshape([&
      &    0.0000567576_wp,   -0.0002590049_wp,    0.0000141659_wp,    0.0001889176_wp,   -0.1158724906_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0019257132_wp,    0.0002977642_wp,   -0.2373074242_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0002044768_wp,    0.0106440802_wp,    0.0226329492_wp,    0.0008169736_wp,   -0.1074306689_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 74) = reshape([&
      &    0.0002368718_wp,   -0.0000915458_wp,    0.0000292389_wp,    0.0001740168_wp,   -0.1618978397_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0075560208_wp,    0.0092933570_wp,   -0.2322286402_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0001035733_wp,    0.0004385680_wp,    0.0009745325_wp,    0.0034613098_wp,   -0.0599180721_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 75) = reshape([&
      &    0.0001291738_wp,   -0.0002662872_wp,    0.0000389531_wp,    0.0001972034_wp,   -0.0788493755_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0008030493_wp,    0.0025279821_wp,   -0.1138125356_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0000530504_wp,    0.0004635277_wp,    0.0004892240_wp,    0.0065054947_wp,   -0.0053078879_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 76) = reshape([&
      &    0.0005138004_wp,   -0.0002540934_wp,    0.0001039399_wp,    0.0005981915_wp,   -0.2774005988_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0021118645_wp,    0.0050262974_wp,   -0.4001459505_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0000944970_wp,    0.0009683771_wp,    0.0006650624_wp,    0.0021549810_wp,   -0.1179017105_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 77) = reshape([&
      &    0.0000987264_wp,   -0.0003188333_wp,    0.0000651020_wp,    0.0001663179_wp,   -0.1268776164_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0015583031_wp,    0.0041455288_wp,   -0.1230518758_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0003017284_wp,    0.0005098442_wp,    0.0008375335_wp,    0.0007627978_wp,   -0.0212219720_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 78) = reshape([&
      &    0.0001058543_wp,   -0.0003113459_wp,    0.0000540767_wp,    0.0002933045_wp,   -0.1981044765_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0036503203_wp,    0.0142484689_wp,   -0.3383232454_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0000860260_wp,    0.0090920803_wp,    0.0008308612_wp,    0.0007091927_wp,   -0.0584087019_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 79) = reshape([&
      &    0.0001153518_wp,   -0.0005256949_wp,    0.0000355296_wp,    0.0002034069_wp,   -0.2169436228_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0018941036_wp,    0.0312998581_wp,   -0.3748966486_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0000790090_wp,    0.0007288050_wp,    0.0009277617_wp,    0.0010303549_wp,   -0.0530588957_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 80) = reshape([&
      &    0.0017329090_wp,   -0.0014857323_wp,    0.0099163652_wp,    0.0047081178_wp,   -0.1085629997_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0108825715_wp,    0.0060640597_wp,   -0.2668197606_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 81) = reshape([&
      &    0.0026092265_wp,   -0.0025600695_wp,    0.0025637936_wp,    0.0059866087_wp,   -0.0691110477_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0029634480_wp,   -0.0066872166_wp,    0.0105986403_wp,   -0.1621737081_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0652185597_wp,   -0.0423637324_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 82) = reshape([&
      &    0.0040015496_wp,   -0.0029532982_wp,    0.0025395021_wp,    0.0068030466_wp,   -0.0503648134_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0028881979_wp,   -0.0030550812_wp,    0.0052453040_wp,    0.0045906228_wp,   -0.2091033489_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.3530470896_wp,   -0.0646290719_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 83) = reshape([&
      &    0.0010270510_wp,   -0.0020117714_wp,    0.0027635110_wp,    0.0016214486_wp,    0.0026618006_wp, &
      &   -0.0315735925_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0024847300_wp,   -0.0015425374_wp,    0.0021954904_wp,    0.0024918401_wp,   -0.1737159730_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0916347992_wp,   -0.0454515167_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 84) = reshape([&
      &    0.0036393083_wp,   -0.0010515856_wp,    0.0030341602_wp,    0.0020879467_wp,    0.0017121894_wp, &
      &   -0.0373302361_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0059932703_wp,   -0.0011957727_wp,    0.0030153740_wp,    0.0026025586_wp,   -0.2574077776_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1912306581_wp,   -0.5034986267_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 85) = reshape([&
      &    0.0024424753_wp,   -0.0010107286_wp,    0.0012850641_wp,    0.0017319439_wp,    0.0015802107_wp, &
      &   -0.0171406864_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0028869978_wp,   -0.0018425066_wp,    0.0095753108_wp,    0.0027478767_wp,   -0.2131708907_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0585392112_wp,   -0.0651771937_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 86) = reshape([&
      &    0.0034128544_wp,   -0.0013780733_wp,    0.0014229514_wp,    0.0013167755_wp,    0.0023195507_wp, &
      &   -0.0169981823_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0023457428_wp,   -0.0034361839_wp,    0.0004255498_wp,    0.0292559868_wp,   -0.0152046243_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.3234416731_wp,   -0.3178862109_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 87) = reshape([&
      &    0.0009304264_wp,    0.0010535475_wp,    0.0010454718_wp,    0.0010997782_wp,    0.0009678037_wp, &
      &    0.0013803830_wp,    0.0010294515_wp,    0.0012367014_wp,   -0.0163723329_wp, & ! shell 1
      &    0.0009398286_wp,    0.0010622010_wp,    0.0009911964_wp,    0.0009130391_wp,    0.0009894133_wp, &
      &    0.0010413611_wp,    0.0011368879_wp,    0.0010914697_wp,   -0.0537729895_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 88) = reshape([&
      &    0.0009455425_wp,    0.0010007950_wp,    0.0009154790_wp,    0.0009643503_wp,    0.0009731908_wp, &
      &    0.0010010192_wp,    0.0010583884_wp,    0.0011155418_wp,   -0.0409708908_wp, & ! shell 1
      &    0.0011629143_wp,    0.0009318714_wp,    0.0011141663_wp,    0.0010330026_wp,    0.0010631764_wp, &
      &    0.0010287207_wp,    0.0010917686_wp,    0.0010361639_wp,   -0.0603527956_wp, & ! shell 2
      &    0.0024267785_wp,    0.0019022832_wp,    0.0059331522_wp,   -0.0629346867_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 89) = reshape([&
      &    0.0000000000_wp,    0.0000000000_wp,   -0.0014667733_wp,   -0.0031449361_wp,   -0.0029689855_wp, &
      &    0.0025166777_wp,   -0.0172796851_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0042112884_wp,   -0.0056318168_wp,    0.0214339323_wp,   -0.0862954332_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0024846817_wp,    0.0040508006_wp,   -0.0127650772_wp,    0.0291475468_wp,   -0.0237214788_wp, &
      &   -0.0549435330_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0028544067_wp,    0.0085975338_wp,    0.0065431414_wp,    0.0147507634_wp,    0.0061260925_wp, &
      &    0.0050391708_wp,   -0.0291826442_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 90) = reshape([&
      &    0.0000000000_wp,    0.0000000000_wp,   -0.0014913992_wp,   -0.0031676519_wp,   -0.0028554980_wp, &
      &    0.0024328467_wp,   -0.0160793529_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0041311214_wp,   -0.0055999643_wp,    0.0201921025_wp,   -0.0804639292_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0024046843_wp,    0.0039058467_wp,   -0.0126205404_wp,    0.0280131612_wp,   -0.0211469399_wp, &
      &   -0.0536548170_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0031484387_wp,    0.0090189678_wp,    0.0076106863_wp,    0.0175016436_wp,    0.0062878564_wp, &
      &    0.0054117700_wp,   -0.0326634834_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 91) = reshape([&
      &    0.0000000000_wp,    0.0000000000_wp,   -0.0014102550_wp,   -0.0029535604_wp,   -0.0025802606_wp, &
      &    0.0022876463_wp,   -0.0022148423_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0037312724_wp,   -0.0055562147_wp,    0.0194341264_wp,   -0.0768074333_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0020068778_wp,    0.0036977546_wp,   -0.0135713755_wp,    0.0310699764_wp,   -0.0198560884_wp, &
      &   -0.0372020661_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0034464042_wp,    0.0103106527_wp,    0.0085963701_wp,    0.0200242327_wp,    0.0061949951_wp, &
      &    0.0048436388_wp,   -0.0347479523_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 92) = reshape([&
      &    0.0000000000_wp,    0.0000000000_wp,   -0.0013915966_wp,   -0.0029256452_wp,   -0.0026388357_wp, &
      &    0.0024064817_wp,   -0.0033121902_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0041974543_wp,   -0.0057287065_wp,    0.0190971325_wp,   -0.0770532522_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0022188558_wp,    0.0038424159_wp,   -0.0124503420_wp,    0.0266368616_wp,   -0.0189432236_wp, &
      &   -0.0714747912_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0032195433_wp,    0.0112324115_wp,    0.0085056388_wp,    0.0223695246_wp,    0.0068980873_wp, &
      &    0.0035634033_wp,   -0.0351623640_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 93) = reshape([&
      &    0.0000000000_wp,    0.0000000000_wp,   -0.0014164283_wp,   -0.0030603986_wp,   -0.0026155415_wp, &
      &    0.0023470443_wp,   -0.0026160991_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0042261790_wp,   -0.0057258215_wp,    0.0187175882_wp,   -0.0773216649_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0021996597_wp,    0.0034373701_wp,   -0.0118018201_wp,    0.0263309936_wp,   -0.0193268981_wp, &
      &   -0.0757983381_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0028774967_wp,    0.0108541000_wp,    0.0074447383_wp,    0.0231706341_wp,    0.0059016031_wp, &
      &    0.0031742967_wp,   -0.0348975801_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 94) = reshape([&
      &    0.0000000000_wp,    0.0000000000_wp,   -0.0014047973_wp,   -0.0031059468_wp,   -0.0020091020_wp, &
      &    0.0024629364_wp,   -0.0022558870_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0040206873_wp,   -0.0041179699_wp,    0.0160288841_wp,   -0.0734668140_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0013069996_wp,    0.0028585360_wp,   -0.0109805730_wp,    0.0271897663_wp,   -0.0187538368_wp, &
      &   -0.0629529030_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0031425062_wp,    0.0100299713_wp,    0.0083139444_wp,    0.0260494552_wp,    0.0036467217_wp, &
      &    0.0042986190_wp,   -0.0345773316_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 95) = reshape([&
      &    0.0000000000_wp,    0.0000000000_wp,   -0.0014184874_wp,   -0.0030716119_wp,   -0.0017739291_wp, &
      &    0.0027140729_wp,   -0.0019764160_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0039728940_wp,   -0.0041425257_wp,    0.0160594006_wp,   -0.0712718903_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0011960569_wp,    0.0028327232_wp,   -0.0104360657_wp,    0.0246362171_wp,   -0.0151785174_wp, &
      &   -0.0921557953_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0024592991_wp,    0.0100850751_wp,    0.0085311201_wp,    0.0272548854_wp,    0.0034309356_wp, &
      &    0.0042269344_wp,   -0.0253685366_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 96) = reshape([&
      &    0.0000000000_wp,    0.0000000000_wp,   -0.0014093246_wp,   -0.0029565878_wp,   -0.0017052181_wp, &
      &    0.0024659479_wp,   -0.0020104778_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0037567203_wp,   -0.0042578908_wp,    0.0158140411_wp,   -0.0628937088_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0012246119_wp,    0.0030098111_wp,   -0.0108496223_wp,    0.0250328710_wp,   -0.0173849023_wp, &
      &   -0.0598396979_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0027198266_wp,    0.0102817166_wp,    0.0087414956_wp,    0.0300939433_wp,    0.0036357818_wp, &
      &    0.0044430367_wp,   -0.0193591011_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 97) = reshape([&
      &    0.0000000000_wp,    0.0000000000_wp,   -0.0012388839_wp,   -0.0028623116_wp,   -0.0017886171_wp, &
      &    0.0025827198_wp,   -0.0019235148_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0037342860_wp,   -0.0044834469_wp,    0.0156046793_wp,   -0.0587540239_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0012025829_wp,    0.0033308802_wp,   -0.0109419612_wp,    0.0259138154_wp,   -0.0155277365_wp, &
      &   -0.0582129879_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0029710089_wp,    0.0113588875_wp,    0.0089485759_wp,    0.0332767568_wp,    0.0040856877_wp, &
      &    0.0046490484_wp,   -0.0210401698_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 98) = reshape([&
      &    0.0000000000_wp,    0.0000000000_wp,   -0.0012421066_wp,   -0.0022115446_wp,   -0.0019706343_wp, &
      &    0.0030343387_wp,   -0.0017100783_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0032932289_wp,   -0.0037758236_wp,    0.0130746925_wp,   -0.0550502825_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0010385421_wp,    0.0037108017_wp,   -0.0123469412_wp,    0.0307280659_wp,   -0.0170126343_wp, &
      &   -0.0816066870_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0025910198_wp,    0.0105373077_wp,    0.0090590856_wp,    0.0271302122_wp,    0.0040288776_wp, &
      &    0.0050215275_wp,   -0.0188689028_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 99) = reshape([&
      &    0.0000000000_wp,    0.0000000000_wp,   -0.0012697954_wp,   -0.0021143573_wp,   -0.0021219603_wp, &
      &    0.0031382349_wp,   -0.0015616319_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0032952386_wp,   -0.0036529062_wp,    0.0136591047_wp,   -0.0498992658_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0010407730_wp,    0.0040689762_wp,   -0.0124374033_wp,    0.0270375638_wp,   -0.0197289799_wp, &
      &   -0.0849531509_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0022448476_wp,    0.0098938096_wp,    0.0088601232_wp,    0.0275787736_wp,    0.0038309816_wp, &
      &    0.0049632789_wp,   -0.0183362858_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 100) = reshape([&
      &    0.0000000000_wp,    0.0000000000_wp,   -0.0011210053_wp,   -0.0020442206_wp,   -0.0016673869_wp, &
      &    0.0033065768_wp,   -0.0016540399_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0028515473_wp,   -0.0034207502_wp,    0.0103543622_wp,   -0.0456630973_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0009971514_wp,    0.0041661174_wp,   -0.0120780890_wp,    0.0278136947_wp,   -0.0240153764_wp, &
      &   -0.0861647316_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0025271940_wp,    0.0118392491_wp,    0.0110356052_wp,    0.0288519193_wp,    0.0038129826_wp, &
      &    0.0047063670_wp,   -0.0214652195_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 101) = reshape([&
      &    0.0000000000_wp,    0.0000000000_wp,   -0.0011475709_wp,   -0.0022305508_wp,   -0.0019439730_wp, &
      &    0.0031818860_wp,   -0.0017806995_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0034640252_wp,   -0.0042780170_wp,    0.0132773665_wp,   -0.0574102413_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0010021946_wp,    0.0036565929_wp,   -0.0117010491_wp,    0.0228562757_wp,   -0.0179239590_wp, &
      &   -0.0742666928_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0026047654_wp,    0.0100820361_wp,    0.0086374156_wp,    0.0310271992_wp,    0.0033744769_wp, &
      &    0.0054640456_wp,   -0.0255931339_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 102) = reshape([&
      &    0.0000000000_wp,    0.0000000000_wp,   -0.0014053586_wp,   -0.0024203428_wp,   -0.0020890092_wp, &
      &    0.0030983945_wp,   -0.0018982354_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0039714497_wp,   -0.0033945344_wp,    0.0142221416_wp,   -0.0558054540_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0009440856_wp,    0.0030978789_wp,   -0.0100447411_wp,    0.0202792496_wp,   -0.0166578372_wp, &
      &   -0.0561793632_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0035369282_wp,    0.0120537343_wp,    0.0079842773_wp,    0.0334275732_wp,    0.0031432516_wp, &
      &    0.0079965938_wp,   -0.0128164284_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
      coefficients_env(:, :, 103) = reshape([&
      &    0.0000000000_wp,    0.0000000000_wp,   -0.0012568206_wp,   -0.0024693353_wp,   -0.0019239256_wp, &
      &    0.0037339322_wp,   -0.0014456821_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0031706633_wp,   -0.0022476270_wp,    0.0100366496_wp,   -0.0491509118_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0008403351_wp,    0.0033944871_wp,   -0.0110707406_wp,    0.0194169554_wp,   -0.0100290389_wp, &
      &   -0.0655883236_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0023593427_wp,    0.0080435407_wp,    0.0128656232_wp,    0.0287143618_wp,    0.0033036594_wp, &
      &    0.0043979774_wp,   -0.0182597575_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/maxg, max_shell/))
      
   end subroutine setCGTOcoefficients_environment

end module tblite_basis_qvszp