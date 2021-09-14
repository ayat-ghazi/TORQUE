module module_input
implicit double precision (a-h,o-z)
 integer,dimension(2):: Nwatom 
 real, dimension(3,3) :: platt
 integer,allocatable,dimension(:):: iaP,lwykf
 real,allocatable,dimension(:,:):: occu
 real,allocatable,dimension(:,:,:)::posatm,Hloc
 Character(len=5),allocatable,dimension(:,:)::cSpcs
 real,allocatable,dimension(:):: achg
 real,dimension(5):: pramt
 integer, ALLOCATABLE,DIMENSION (:) :: ip_Ox,iwyk_Ox,ip_Hy
 integer, allocatable,dimension (:,:) :: iwpos	
 real,dimension(1000,5):: sym


end module module_input


