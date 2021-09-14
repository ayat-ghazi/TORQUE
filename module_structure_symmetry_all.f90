module module_structure_symmetry_all
implicit double precision (a-h,o-z)

integer, parameter :: nsgmax=700
integer nsg
dimension ntrigonal(nsgmax),nunique(nsgmax),norigin(nsgmax)
dimension ssym(3,3,48,48,nsgmax),ttau(3,48,48,nsgmax),nspace(nsgmax)
dimension nshift(nsgmax),vshift(3,4,nsgmax),multi(48,nsgmax)
dimension nsites(nsgmax),nunique_alt(nsgmax),ncell_alt(nsgmax)
dimension number_sg(nsgmax)
dimension xorigin(3,nsgmax)
dimension qqall(3,3,nsgmax),ppall(3,3,nsgmax)
character (len=100) :: cwyck(48,nsgmax),cschoen(nsgmax)
character (len=100) :: cpg(48,nsgmax)
character (len=100) :: cspace(nsgmax),cspace_alt(nsgmax)
character (len=100) :: cspace_alt_ext(nsgmax)
character (len=100) :: csgpath
integer cspace_id(2,nsgmax),cspace_alt_id(2,nsgmax)
integer cspace_alt_ext_id(2,nsgmax)
!!!!! az inja male modul_structure_symmetry_local hast@@@
integer :: ntrigonall,nuniquel,noriginl,nsitesl,nsgtot
integer :: nshiftl,multil(48),number_sgl
integer :: nuniquel_alt,ncelll_alt
dimension ssyml(3,3,48,48),ttaul(3,48,48)
dimension vshiftl(3,4)
dimension xoriginl(3)
character (len=100) :: cwyckl(48),cspacel
character (len=100) :: cpgl(48)
character (len=100) :: cspacel_alt
character (len=100) :: cspacel_alt_ext,cschoenl

! integer, ALLOCATABLE,DIMENSION (:) :: n_oxidation,lO_elm,lO_wykf
 integer, ALLOCATABLE,DIMENSION (:) :: lO_elm,lO_wykf,lH_ip
 integer,allocatable,dimension(:,:)::n_wykf_oxid
 real,allocatable,dimension(:,:,:):: xyz_wykf_cartz
 real,allocatable,dimension(:,:,:):: xyz_H_cartz
 real,allocatable,dimension(:,:,:):: xyz_H_lattice
 real,allocatable,dimension(:,:,:):: force_H
 real,allocatable,dimension(:,:,:):: torque_H
 real,allocatable,dimension(:,:,:):: delta_force,delta_torque,ox_shift_cz,ox_shift_latt
 real,allocatable,dimension(:) :: delta_charge
 real,allocatable,dimension(:,:):: H_lattice_POSCAR
 real :: delta_distance,delta_q,delta_angle
 integer:: iwater,model
 real :: w_angle
 real:: h_oxidation,h_charge,h_oxidation_init
 real,allocatable,dimension(:):: n_oxidation !,oxidation_init !,atom_charge !####
dimension qq(3,3),pp(3,3)

end module module_structure_symmetry_all
