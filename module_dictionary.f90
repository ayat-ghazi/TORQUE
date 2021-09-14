module module_dictionary
implicit double precision (a-h,o-z)

 CHARACTER (len=50):: c_name
 CHARACTER (len=15):: c_SG
 real,dimension(6):: acell
 INTEGER:: list_SG,num_SG,natypes,iaPos_max,iwykf_max,list_max,nfold_max
 real:: dis_cutoff,d_O_H
 real, dimension(3,3):: T_metric, aL_matrix
 integer :: mainctrl,l_adrs_ary

 integer,allocatable,dimension(:):: iaPos,iwykf
 real,allocatable,dimension(:,:):: occ
 real,allocatable,dimension(:,:,:)::r_wykf,ave_hydrogen_coord
 real,allocatable,dimension(:,:):: hydrogen_coords 
 real,allocatable,dimension(:,:,:,:) ::R_H_kol
 Character(len=5),allocatable,dimension(:,:)::cAMC
 character (len=100),allocatable,dimension (:) :: c_line
 character (len=50),allocatable,dimension (:) :: c_address
 real:: TTmax,TTmin,eps
 integer:: Seed
 integer, allocatable, dimension (:) :: iseed,iseed_first,iseed_last ! seed array

end module module_dictionary


