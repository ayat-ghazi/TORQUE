subroutine  arries_assigning(sym2)
use module_dictionary
use module_structure_symmetry_all
implicit double precision (a-h,o-z)

 integer:: i,j,ii,jj,k ! ,ierror,stats,nline,n_column,ll,ip,iw,kntrl,nt1,nt2
 integer:: np,icnt_empthy, icnt_comment,il
 real,dimension(3)::syml
 real:: taul 
 integer:: intg_cnt,it_cnt,ii_cnt,isyz
 logical :: fexist
 real,dimension(1000,5),intent(in):: sym2



acell=0.
do j=1,3
acell(1)=acell(1)+aL_matrix(j,1)**2
acell(2)=acell(2)+aL_matrix(j,2)**2
acell(3)=acell(3)+aL_matrix(j,3)**2 
enddo

acell(1)=sqrt(acell(1))
acell(2)=sqrt(acell(2))
acell(3)=sqrt(acell(3))

do j=1,3
acell(6)= acell(6)+ aL_matrix(j,1)*aL_matrix(j,2) !gama
acell(5)= acell(5)+ aL_matrix(j,3)*aL_matrix(j,1) !beta
acell(4)= acell(4)+ aL_matrix(j,2)*aL_matrix(j,3) !alpha
enddo

acell(4)= acos( acell(4)/( acell(2)*acell(3)) )*180/3.14159265
acell(5)= acos( acell(5)/( acell(3)*acell(1)) )*180/3.14159265
acell(6)= acos( acell(6)/( acell(1)*acell(2)) )*180/3.14159265

if (acell(4)> 180 .or. acell(5)> 180 .or. acell(6)> 180) then
print *, "Please check the lattice angles"
!read *, ieror
stop
endif 



!read(c_line(nn_MT+i),*) (T_metric(i,j),j=1,3 ) !  !metric tensor
T_metric(1,1)=acell(1)**2
T_metric(2,2)=acell(2)**2
T_metric(3,3)=acell(3)**2

T_metric(1,2)=acell(1)*acell(2)*cos(acell(6)*3.14159265/180.)
T_metric(1,3)=acell(1)*acell(3)*cos(acell(5)*3.14159265/180.)
T_metric(2,3)=acell(2)*acell(3)*cos(acell(4)*3.14159265/180.)

T_metric(2,1)=T_metric(1,2)
T_metric(3,1)=T_metric(1,3)
T_metric(3,2)=T_metric(2,3)
!--------------------------------------------

do i=1,6
 if ((sym2(i,1) -1.) < 0.0001 .and. (sym2(i,1) -1.) > -0.0001 ) then 
	nshiftl=i-1
	exit
 endif
enddo



do i=1,nshiftl
vshiftl(1,i)=sym2(i,2)
vshiftl(2,i)=sym2(i,3)
vshiftl(3,i)=sym2(i,4)
enddo

it_cnt=1
ii_cnt=0
do i=nshiftl+1,1000
	if (sym2(i,1) < 0.0001 .and. sym2(i,1) > -0.0001) exit
	intg_cnt=int(sym2(i,1))
	if (intg_cnt==it_cnt) then
	ii_cnt=ii_cnt+1
	multil(1)=intg_cnt
	ssyml(ii_cnt,1,intg_cnt,1)=sym2(i,2)
	ssyml(ii_cnt,2,intg_cnt,1)=sym2(i,3)
	ssyml(ii_cnt,3,intg_cnt,1)=sym2(i,4)
	ttaul(ii_cnt,intg_cnt,1)= sym2(i,5)
!	print '(I5,3f12.6,10x,f12.6)', intg_cnt,(ssyml(ii_cnt,j,it_cnt,1),j=1,3),ttaul(ii_cnt,it_cnt,1)
	else
	it_cnt=intg_cnt
	ii_cnt=1
	ssyml(ii_cnt,1,intg_cnt,1)=sym2(i,2)
	ssyml(ii_cnt,2,intg_cnt,1)=sym2(i,3)
	ssyml(ii_cnt,3,intg_cnt,1)=sym2(i,4)
	ttaul(ii_cnt,intg_cnt,1)= sym2(i,5)
!	print '(I5,3f12.6,10x,f12.6)', intg_cnt,(ssyml(ii_cnt,j,it_cnt,1),j=1,3),ttaul(ii_cnt,it_cnt,1)
	endif



enddo


!print *, intg_cnt, multil(1)


return
END

subroutine  column_finder(c_string,n_column)
use module_dictionary
implicit double precision (a-h,o-z)
 Character (len=100), intent (in):: c_string
 integer, intent (out):: n_column
 integer:: n_on,n_off,i,ltr
 n_column=0
 n_off=1

  ltr=len_trim(c_string)

 do i=1,ltr
 if (c_string(i:i)==" ") then
  n_off=1
 endif
 if (c_string(i:i)/=" " .and. n_off==1) then
  n_off=0
  n_column=n_column+1
 endif
 enddo

! write(*,*) n_column

return
end
      subroutine q_multi(qtmp,bb,aa)
      implicit real*8 (a-h,o-z)
      real,dimension(4),intent(out)::qtmp
      real,dimension(4),intent(in)::aa,bb 
!      dimension qtmp(4),bb(4),aa(4)
!c
      qtmp(1)=bb(1)*aa(1)-bb(2)*aa(2)-bb(3)*aa(3)-bb(4)*aa(4)
      qtmp(2)=bb(1)*aa(2)+bb(2)*aa(1)+bb(3)*aa(4)-bb(4)*aa(3)
      qtmp(3)=bb(1)*aa(3)+bb(3)*aa(1)+bb(4)*aa(2)-bb(2)*aa(4)
      qtmp(4)=bb(1)*aa(4)+bb(4)*aa(1)+bb(2)*aa(3)-bb(3)*aa(2)
!c
      return
      end
      subroutine qxqml2(hklrot,hklc,rot)
      implicit real*8 (a-h,o-z)
      real,dimension(4),intent(in):: rot
      real,dimension(4),intent(out)::hklrot
      real,dimension(3),intent(in):: hklc
      real,dimension(4):: quat(4),quatm1(4),vv(4)!,rot(4),hklc(3),hklrot(4)
      real,dimension(4):: qtmp

      vv(1)=0.d0
      do ii=1, 3
       vv(ii+1)=hklc(ii)
      enddo

      sina2=sin(rot(1)/2.d0)
      cosa2=cos(rot(1)/2.d0)

      quat(1)=cosa2
      quat(2)=sina2*rot(2)
      quat(3)=sina2*rot(3)
      quat(4)=sina2*rot(4)

      quatm1(1)=cosa2
      quatm1(2)=-(sina2*rot(2))
      quatm1(3)=-(sina2*rot(3))
      quatm1(4)=-(sina2*rot(4))
      call q_multi(qtmp,quat,quatm1)
      call q_multi(qtmp,vv,quatm1)
      call q_multi(hklrot,quat,qtmp)

      return
      end
subroutine  read_input_file_2nd(c_filename)
use module_input
implicit double precision (a-h,o-z)

 character(len=100),intent(in):: c_filename
 character(len=100)::c_temp
 character(len=50):: c_temp2,c_file,c_file_adres_total
 character(len=30):: cct
 character(len=3) :: c_fldr
 character(len=1) :: cPF
 integer:: i,j,ii,jj,k,ierror,stats,nline,n_column,ll,ip,iw,kntrl,nt1,nt2
 integer::nn_name,nn_SG,nn_SG_ext,nn_MT,nn_ATB3,nn_wykf,nn_cut_off,nn_coords,nn_LM
 integer:: nn_lattice,nn_ATB1,nn_ATB2,nn_O,nn_H,nn_model,nn_symm,nn_pramt
 real:: x,y,z,oc
 integer,dimension(2000):: NI1_w,ijkf
 integer:: np,nw,icnt_empthy, icnt_comment,il,iwt
 real,dimension(3)::syml
 real:: taul 
! integer:: nshiftl
 integer:: intg_cnt,it_cnt,ii_cnt

 logical :: fexist

 character (len=100),allocatable,dimension (:) :: c_stFile




 c_file="./"//trim(c_filename)
inquire (file=c_file , exist= fexist)
if (fexist) then
 c_file_adres_total=c_file

else
 print *, "The input file name has not found"
 stop
endif


 open (unit=100, file=c_file_adres_total,status='old', action='read' , IOstat= ierror)
nline=0



 do        
	read (100,'(A100)',IOstat=stats) c_temp
	if (stats/=0) exit
	nline=nline+1
 end do  
 REWIND (UNIT=100)



ALLOCATE (c_stFile(nline),STAT= stats) 

 c_stFile(:)= ""

 il=0
 reading: DO i=1,nline
  read (100,'(A100)',IOstat=stats) c_temp
 
	icnt_empthy=1
	icnt_comment=0
	 do j=1,100
	 	 if (c_temp(j:j) /= " " ) then
			
			icnt_empthy=0
			if (index(c_temp(j:j),"#") /=0 .or. index(c_temp(j:j),"!") /=0 .or. index(c_temp(j:j),"c") /=0 ) then
			 icnt_comment=1
			!print *, j, c_temp(j:j) 
			endif
			exit
		endif
	 enddo

	if (icnt_empthy ==0 .and. icnt_comment==0 )  then
	 il=il+1
	c_stFile(il)=c_temp

	endif


  END DO reading

	nline=il



  pramt(1)=0.00001
  pramt(2)=0.5
  pramt(3)=0.005
  pramt(4)=12345.
  pramt(5)=1.
do i=1,nline
  if (index(c_stFile(i),"LATTICE_MATRIX_")/=0 .and.index(c_stFile(i),"END_LATTICE_")==0 ) nn_lattice=i
  if (index(c_stFile(i),"ATOMIC_SPECIES_CHARGES")/=0 .and.index(c_stFile(i),"END_")==0) nn_ATB1=i
  if (index(c_stFile(i),"ATOMIC_POSITIONS_OCCUPANCY")/=0 .and.index(c_stFile(i),"END_")==0) nn_ATB2=i
  if (index(c_stFile(i),"WATER_IDENTIFICATION_NUMBER_")/=0 .and. index(c_stFile(i),"END_")==0) nn_O=i
  if (index(c_stFile(i),"HYDROGEN_IDENTIFICATION_NUMBER_")/=0 .and. index(c_stFile(i),"END_")==0) nn_H=i
  if (index(c_stFile(i),"WATER_MODEL_")/=0 .and. index(c_stFile(i),"END_")==0) nn_model=i
  if (index(c_stFile(i),"SYMMETRY_OPERATORS_")/=0 .and. index(c_stFile(i),"END_")==0) nn_symm=i
  if (index(c_stFile(i),"PARAMETERS_")/=0 .and. index(c_stFile(i),"END_")==0) nn_pramt=i

  if (index(c_stFile(i),"Threshold")/=0)  then
	do j=1,100
	  if (c_stFile(i)(j:j)=="=") then
	  	cct=c_stFile(i)(j+1:j+30)
		exit
	  endif
	enddo
	read (cct,*) pramt(1)
  endif

  if (index(c_stFile(i),"Max_Angle")/=0)  then
	do j=1,100
	  if (c_stFile(i)(j:j)=="=") then
	  	cct=c_stFile(i)(j+1:j+30)
		exit
	  endif
	enddo
	read (cct,*) pramt(2)
  endif

  if (index(c_stFile(i),"Min_Angle")/=0)  then
	do j=1,100
	  if (c_stFile(i)(j:j)=="=") then
	  	cct=c_stFile(i)(j+1:j+30)
		exit
	  endif
	enddo
	read (cct,*) pramt(3)
  endif

  if (index(c_stFile(i),"Seed")/=0)  then
	do j=1,100
	  if (c_stFile(i)(j:j)=="=") then
	  	cct=c_stFile(i)(j+1:j+30)
		exit
	  endif
	enddo
	read (cct,*) pramt(4)
  endif


  if (index(c_stFile(i),"Iteration")/=0)  then
	do j=1,100
	  if (c_stFile(i)(j:j)=="=") then
	  	cct=c_stFile(i)(j+1:j+30)
		exit
	  endif
	enddo
	read (cct,*) pramt(5)
  endif

enddo

 close (100)


 
Nwatom(1)=nn_ATB2-nn_ATB1-2   !max of main positions(does not consider wyckoff posisions)


iaP_max=50  !this is for setting maximum numberof atom species in one position
lwykf_max=200 !max of wyckoff position is 192

list_max=nn_O-nn_ATB2-2 !max of all positions


ALLOCATE (iaP(Nwatom(1)),STAT= nst)
ALLOCATE (lwykf(Nwatom(1)),STAT= nst)
ALLOCATE (occu(Nwatom(1),lwykf_max),STAT= nst)
ALLOCATE (posatm(Nwatom(1),lwykf_max,3),STAT= nst)
ALLOCATE (cSpcs(Nwatom(1),1),STAT= nst)
ALLOCATE (achg(Nwatom(1)),STAT= ierr)
ALLOCATE (ip_Hy(Nwatom(1)),STAT= ierr)
ALLOCATE (iwpos(Nwatom(1),3),STAT=ierr)


iapos=1






!*******************************Lattice
do i=1,3
read(c_stFile(nn_lattice+i),*) (platt(j,i),j=1,3 ) !lattice matrix **** ein bayad taghir kone too kole barname, jaye i,j avaz beshe
enddo
!!


!do i=1,3
!print '(3f12.6)',  (T_metric(i,j),j=1,3)
!enddo



!Nwatom(1)==nn_ATB2-nn_ATB1-2
do ip=1,Nwatom(1)
read(c_stFile(ip+nn_ATB1),*)  k,cSpcs(ip,1),achg(ip)

if (k/=ip) then

Print *, "There is a mismatch in the sequence of the ATOMIC_SPECIES_CHARGES_ section. Please look at the input file"
print *, "At line :", k, "IN1 = ", ip
stop
endif

enddo


lwykf(:)=0
!ip, lwykf(ip), posatm(ip,iw,123),occu(ip,1)


icnt=0
do j=nn_ATB2+1,nn_O-2
 call column_finder(c_stFile(j),n_column)
 if (n_column==5) then
	read(c_stFile(j),*) np,x,y,z,oc
 else if (n_column==4) then
	read(c_stFile(j),*) np,x,y,z	
 	oc=1.0
 else
	print *, "There is a problem in  ATOMIC_POSITIONS_OCCUPANCY_ section of input file"
	print *, "Please look at TORQUE_USER_Manual file"
	stop
 endif

!print *, np,x,y,z,oc
 if (j==nn_ATB2+1) then
   ip=np
   lwykf(ip)=1
   posatm(ip,lwykf(ip),1)=x
   posatm(ip,lwykf(ip),2)=y
   posatm(ip,lwykf(ip),3)=z
   occu(ip,lwykf(ip))=oc
 !  icnt=icnt+1
 else
	if (np==ip) then
		lwykf(ip)=lwykf(ip)+1
		posatm(ip,lwykf(ip),1)=x
		posatm(ip,lwykf(ip),2)=y
	   	posatm(ip,lwykf(ip),3)=z
		occu(ip,lwykf(ip))=oc

	else
		ip=np
		lwykf(ip)=1
		posatm(ip,lwykf(ip),1)=x
		posatm(ip,lwykf(ip),2)=y
	   	posatm(ip,lwykf(ip),3)=z
		occu(ip,lwykf(ip))=oc
	endif

 endif

! print '(2 I5,3f12.6,f5.2)', ip,lwykf(ip),posatm(ip,lwykf(ip),1),posatm(ip,lwykf(ip),2),posatm(ip,lwykf(ip),3),occu(ip,lwykf(ip))

enddo

!-----------------
iwpos=0
do i=1, Nwatom(1)
  iwpos(i,1)=i
enddo



!***************************************water positions
 iwt=0
 call column_finder(c_stFile(nn_O+1),n_column)
! print *, " n col= ", n_column
 NI1_w=0
 if (index(c_stFile(nn_O+1),"END")==0) then
 read(c_stFile(nn_O+1),*)  (NI1_w(j),j=1,n_column) 
 else
 Print *, "Please specify WATER_IDENTIFICATION_NUMBER_ in the input file."
 Print *, "For more information, please read the TORQUE_User_Manual."
 Stop
 endif
 


 do j=1,n_column
   iwt=iwt+ lwykf(NI1_w(j)) 
   iwpos(j,2)=NI1_w(j)
!   print *, j,NI1_w(j), lwykf(NI1_w(j))
 enddo
 ! Nwater was set


!===================================================
Nwatom(2)=iwt


!Nwater=iwt 
ALLOCATE (ip_Ox(Nwatom(2)),STAT= ierr) 
ALLOCATE (iwyk_Ox(Nwatom(2)),STAT= ierr) 
ALLOCATE (Hloc(Nwatom(2),2,3),STAT= ierr)
!print *, "size:  ",size(iwyk_Ox),ierr 

iwat=0
 do j=1,n_column
 do jj=1,lwykf(NI1_w(j))
   iwat=iwat+1
   ip_Ox(iwat)=NI1_w(j)
   iwyk_Ox(iwat)=jj
 enddo
enddo

do iwat=1,Nwatom(2)

!print *, ip_Ox(iwat),iwyk_Ox(iwat)
enddo
!----------------------------------------------------
!




!*******************Hydrogen
ip_Hy=0

 call column_finder(c_stFile(nn_H+1),n_column)
 NI1_w(:)=0


 if (index(c_stFile(nn_H+1),"END")==0) then
 read(c_stFile(nn_H+1),*)  (NI1_w(j),j=1,n_column)  
 do j=1,n_column
 ip_Hy(NI1_w(j))=NI1_w(j)
 iwpos(j,3)=NI1_w(j)
 enddo
 endif


!-----------------------
!do i=1,Nwatom(1)
!print *, (iwpos(i,j),j=1,3)
!enddo

Hloc=0.


!*************** water model
call column_finder(c_stFile(nn_model+1),n_column)
!print *, "K ",nn_model+1,n_column,trim(c_stFile(nn_model+1))

 read (c_stFile(nn_model+1),*)  model !water model


!*****************Symetry
sym=0.
icnt=1
isym=5
do i=nn_symm+1,nn_pramt-2	!!***
read (c_stFile(i),*) intg_cnt,syml(1),syml(2),syml(3),taul
	sym(icnt,1)=intg_cnt
	sym(icnt,2)=syml(1)
	sym(icnt,3)=syml(2)
	sym(icnt,4)=syml(3)
	sym(icnt,5)=taul
icnt=icnt+1
enddo		!!**


if (icnt==1) then
sym=0.
	do i=2,4
		sym(i,1)=1
		sym(i,i)=1.
	enddo
endif
!**********************

 DEALLOCATE(c_stFile)




return
END




subroutine  tabbe2_calculate_forece_on_an_atom(ipH,iwH,N,nH,nTIP,irn)
use module_dictionary
use module_structure_symmetry_all
implicit double precision (a-h,o-z)

integer,intent(in):: ipH,iwH,N,nH,irn,nTIP
!ipH : shomare molecool AAb ro neshoon mide
!iwH: shomare H oon molecool AAb

real,dimension(3) :: r1,r2,r_atom_cz,r_atom_l,R_H_elm,R12,H_transferd,F !in F bayad taghir kone hala, moteghayer global beshe ehtemalan
real::dis,d1,d2,dis_max,occy
integer :: ix,iy,iz,ip,iw,ierr,icont,jtip,Nm
!real:: start,finish



Nm=N

dis_max=45000.
F=0.
occy=1.

!write(*,*) "#############", ipH,lO_elm(ipH)

do ip=1,natypes
icont=0
do iwat=1,iwater
 if (lO_wykf(iwat)==1) then
	if (ip==lO_elm(iwat)) icont=1
 endif
enddo


!read (*,*) lll
if (icont==1) goto 155 !Nirooye hasel az Ox haye Aab ro mohasebe nemikonim 
!*************************
!************************* Nirooye hasel az H ha mohasebe mishe!!
!************************ 

if (lH_ip(ip)/=0 ) goto 155
!print *, "force ip=", ip

do iw=1,iwykf(ip)
!write(*,*) ">>>  ",ip, occ(ip,iw),occ(lO_elm(ipH),1),lO_elm(ipH),ipH


do ix=-1*Nm , N
do iy=-1*Nm , N
do iz=-1*N , N
		


	r_atom_l(1)=r_wykf(ip,iw,1)+ix
	r_atom_l(2)=r_wykf(ip,iw,2)+iy
	r_atom_l(3)=r_wykf(ip,iw,3)+iz


!write(*,*) cAMC(ip,1),n_oxidation(ip), r_atom_l(:)
	call tabbe_tabdil_lattic_cartezian(r_atom_l,r_atom_cz)
!write(*,*) ip, iw, r_atom_cz(:)

	do i=1,3
 	!R_H_elm(i)=r_atom_cz(i)-xyz_H_cartz(ipH,iwH,i) !bordar vasel H be atom digar==> jahat bordar F niroo
	 R_H_elm(i)=xyz_H_cartz(ipH,iwH,i)-r_atom_cz(i)
!	**	R_H_elm(i)=H_transferd(i)-r_atom_cz(i)
	enddo





!	call tabbe_direct_distance_cartz(H_transferd(:),r_atom_cz(:),dis)
	call tabbe_direct_distance_cartz(xyz_H_cartz(ipH,iwH,:),r_atom_cz(:),dis)
	!dis ==> fasele ta atom mored nazar




if (dis==0) then
write(*,*) "dis==0, elements", ip,iw ,ipH,iwH
read (*,*) i
endif

!if (dis>40) goto 45    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!in khat 1.3.2016 ezafe shod baraye inke bekhaym niro ro motegaren konim !**************************************
!if (abs(R_H_elm(1))> N*a_cell(1) .or. abs(R_H_elm(2))> N*a_cell(2) .or. abs(R_H_elm(3))> N*a_cell(3)  ) goto 45!***
!***************************************************************************************************************

	!baar e atom H ra +1 dar nazar gereftim
	!......mohasebe niroo
	if (dis <dis_max ) then
	do i=1,3
	F(i)=F(i)+ ( ( occ(ip,iw)*n_oxidation(ip)*occ(lO_elm(ipH),1)*h_oxidation )/(dis**3) )*R_H_elm(i)
	enddo
	endif

45 continue

enddo
enddo
enddo

!write(*,*) "!!!!!!!!!!!!!!!!!!!"
55 continue
enddo !!iw 


155 continue
enddo !ip

!write (*,'(A10,3f10.5)') "Force...>" , force_H(ipH,iwH,:)

!if (irn==250) write(*,*) "######"
!print '(A15,3f12.6)' , "Force...>" , F(:)







!goto 35
do jtip=1,nTIP     !tamame molecoolhaye aab ro mesle model dar nazar gereftim
do iwat=1,iwater
!do iw=1,lO_wykf(iwat)

!write (*,*) ox_shift_latt(iwat,jtip,:)
!read (*,*) iiiii

do ix=-1*Nm , N
do iy=-1*Nm , N
do iz=-1*N , N
		

!ipH is the associated number of water molecule, which its wyckoff position is 1 ==> just this water would be sent to this subroutine
if (iwat /= ipH .or. ix /= 0 .or. iy /= 0 .or. iz /=0 ) then
!if (iwat /= ipH ) then

	r_atom_l(1)=ox_shift_latt(iwat,jtip,1)+ix
	r_atom_l(2)=ox_shift_latt(iwat,jtip,2)+iy
	r_atom_l(3)=ox_shift_latt(iwat,jtip,3)+iz

	call tabbe_tabdil_lattic_cartezian(r_atom_l,r_atom_cz)

	do i=1,3
	 R_H_elm(i)=xyz_H_cartz(ipH,iwH,i)-r_atom_cz(i)
	enddo

	call tabbe_direct_distance_cartz(xyz_H_cartz(ipH,iwH,:),r_atom_cz(:),dis)

if (dis==0) then
write(*,*) "dis==0, elements", ip,iw ,ipH,iwH
!read (*,*) i
endif
	if (dis > 0.0 ) then
	do i=1,3
	F(i)=F(i)+ ( ( (occ(lO_elm(ipH),1)**2)  *delta_q*h_oxidation )/(dis**3) )*R_H_elm(i) !*** IN dastoore asli ast
	enddo
	endif

!************************************************************************************
!************************** In tike baraye Grimselite Ezafe shode
!	do i=1,3
!	if (ix==0 .and. iy==0 .and. iz==0 .and. iwat==3) then
!	F(i)=F(i)+ ((delta_q*h_oxidation )/(dis**3) )*R_H_elm(i)
!	elseif (ix==0 .and. iy==0 .and. iz==0 .and. iwat==2) then
!	F(i)=F(i)
!	elseif (ix==0 .and. iy==0 .and. iz==0 .and. iwat==4) then
!	F(i)=F(i)
!	else
!	F(i)=F(i)+ 0.5*((delta_q*h_oxidation )/(dis**3) )*R_H_elm(i)
!	endif		
!	enddo
	
!*************************************************************************************

endif !

enddo
enddo
enddo
!enddo

enddo
enddo !jtip
35 continue









!write(*,*) "~~~~~~~~~~~~>  ",irn
!if (irn<1) goto 75   !dar irn<n ghadam aval nirooye hasel az H ha ro hesab nemikone
!.................................niroo ye hasel az H ha

do inpH=1,iwater  !***************************************************************



do inwH=1,2   !har AAb 2 ta H dare     !iwater*2

!if (ipH==inpH) then

!write(*,*) "*****", xyz_H_lattice(inpH,inwH,:)
! goto 65   !nirooye hasel az atom H yek molekool water ro dar nazar nemigire
!if (iwH==inwH) goto 65
!if (irn==250) write(*,*) inpH
!endif


!if (irn==700 .and. ipH==1 .and. iwH==2) write(373,'(2I5,3f12.7)')  inpH,inwH, xyz_H_lattice(inpH,inwH,:)


do ix=-1*Nm , N
do iy=-1*Nm , N
do iz=-1*N , N

if (ipH /=inpH .or. ix /= 0 .or. iy /= 0 .or. iz /=0 ) then
!if (iwH /=inwH ) then
	r_atom_l(1)=xyz_H_lattice(inpH,inwH,1)+ix
	r_atom_l(2)=xyz_H_lattice(inpH,inwH,2)+iy
	r_atom_l(3)=xyz_H_lattice(inpH,inwH,3)+iz

	call tabbe_tabdil_lattic_cartezian(r_atom_l,r_atom_cz)

	do i=1,3
 	!R_H_elm(i)=r_atom_cz(i)-xyz_H_cartz(ipH,iwH,i) !bordar vasel H be atom digar==> jahat bordar F niroo
	R_H_elm(i)=xyz_H_cartz(ipH,iwH,i)-r_atom_cz(i)
!	R_H_elm(i)=H_transferd(i)-r_atom_cz(i)
	enddo

!	call tabbe_direct_distance_cartz(H_transferd(:),r_atom_cz(:),dis)
	call tabbe_direct_distance_cartz(xyz_H_cartz(ipH,iwH,:),r_atom_cz(:),dis)
	!dis ==> fasele ta atom mored nazar

!if (dis>40) goto 60    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	!bar atom H ra +1 dar nazar gereftim
	!......mohasebe niroo

!in khat 1.3.2016 ezafe shod baraye inke bekhaym niro ro motegaren konim !**************************************
!if (abs(R_H_elm(1))> N*a_cell(1) .or. abs(R_H_elm(2))> N*a_cell(2) .or. abs(R_H_elm(3))> N*a_cell(3)  ) goto 60!***
!***************************************************************************************************************


if (irn==700 .and. ipH==1 .and. iwH==2) then
!write(373,'(f12.7)') dis
! write (373,*) "555555" 
! write(*,*) "!!!!", irn
! read (*,*) i
endif







	if (dis > 0. ) then
	do i=1,3
	F(i)=F(i)+  ( ( (occ(lO_elm(ipH),1)**2)*(h_oxidation**2) ) /(dis**3) ) * R_H_elm(i)   ! in Dastoor asli ast
	enddo
	endif

!************************************************************************************
!************************** In tike baraye Grimselite Ezafe shode
!	do i=1,3
!	if (ix==0 .and. iy==0 .and. iz==0 .and. iwat==3) then
!	F(i)=F(i)+ (h_oxidation**2/(dis**3) )*R_H_elm(i)
!	elseif (ix==0 .and. iy==0 .and. iz==0 .and. iwat==2) then
!	F(i)=F(i)
!	elseif (ix==0 .and. iy==0 .and. iz==0 .and. iwat==4) then
!	F(i)=F(i)
!	else
!	F(i)=F(i)+ 0.5*(h_oxidation**2/(dis**3) )*R_H_elm(i)
!	endif		
!	enddo
	
!*************************************************************************************






!endif !iwH

endif

60 continue
enddo
enddo
enddo

65 continue
enddo
enddo



75 continue

do i=1,3
force_H(ipH,iwH,i)=F(i)
enddo

!write (*,'(A10,3f10.5)') "Force...>" , force_H(ipH,iwH,:)
!read (*,*) i

! close (373)
!CALL CPU_TIME(finish)
!write(*,'(A10,f12.8)') "Time:", finish-start
!read (*,*) mmmm

return
end subroutine tabbe2_calculate_forece_on_an_atom
!niroo be yek poit
!mokhtasate cartezi migire
subroutine  tabbe2_calculate_forece_on_a_point(point_cz, N,FF)
use module_dictionary
use module_structure_symmetry_all
implicit double precision (a-h,o-z)

integer,intent(in):: N
real,dimension(3),intent(inout) :: point_cz 
real,dimension(3),intent(out) :: FF
!ipH : shomare molecool AAb ro neshoon mide
!iwH: shomare H oon molecool AAb

real,dimension(3) :: r1,r2,r_atom_cz,r_atom_l,R_H_elm,F !in F bayad taghir kone hala, moteghayer global beshe ehtemalan
real::dis
integer :: ix,iy,iz,ip,iw

!-1.14819   0.14076   5.65789
!point_cz(1)=-1.14819
!point_cz(2)=0.14076
!point_cz(3)=5.65789


F=0.

!write(*,*) "#############", ipH,lO_elm(ipH)


do ip=1,natypes
do iw=1,iwykf(ip)
 
!if (ip==1 .and. iw==1 )  write(*,*) "++++", xyz_H_lattice(ip,iw,:)


! ip hydrogenha ro mostasna mikonim
!if (ip==9 .or. ip==10 ) goto 55   !Haidingerite 10412
!if (ip==6 .or. ip==7 ) goto 55   !BaCl2_2H2O  20023
!if (ip==8 .or. ip==9 ) goto 55   !BaCl2_2H2O   20023

!if (ip==5 .or. ip==6 ) goto 55   !BaCl2_H2O   20024

!if (ip==9 .or. ip==10 ) goto 55  	!Apophyllite 00006
!if (ip==11 ) goto 55                    !Apophylite   
!if (ip==35 .or. ip==36 ) goto 55   !blue
!if (ip==37 .or. ip==38 ) goto 55  !blue
!if (ip==13 .or. ip==14 ) goto 55   !adolf
!if (ip==9 .or. ip==10 ) goto 55   !Hydroxyapophyllite
!if (ip==20 .or. ip==22 ) goto 55   !Kernite
!if (ip==18 .or. ip==25 ) goto 55   !Kernite 00363

!if (ip==4 .or. ip==5 ) goto 55   ! NaNiF3(3H2O) 30005 20025

!if (ip==20 .and. iw==1 ) then
!write(*,*) "==1===", r_wykf(ip,iw,:)
! goto 55   !Kernite
!endif

!if (ip==22 .and. iw==1 ) then
!write(*,*) "==2===", r_wykf(ip,iw,:)
! goto 55   !Kernite
!endif
!if (irn==250) write(*,*) ip,iw

!nemikhahim oxygenhayi ke vojood nadaran hesab shavand
!*****************   agar atomi nabayad hazf beshe inja bayad moshakhas beshe
!******************  agar nabayad chizi hazf beshe bere be khat 225
goto 225
if (ip==3) then
if (iw==11 .or. iw==12) goto 55  
if (iw==15 .or. iw==16) goto 55  
if (iw==19 .or. iw==20) goto 55  
if (iw==23 .or. iw==24) goto 55  
if (iw==27 .or. iw==28) goto 55  
if (iw==31 .or. iw==32) goto 55  

endif
!225 continue

if (ip==8) then
if (iw==2 .or. iw==4) goto 55    !for grimselite
endif
225 continue
!!******************************************************************************

!write(*,*) "ff", r_wykf(ip,iw,:)


do ix=-1*N , N
do iy=-1*N , N
do iz=-1*N , N
		


	r_atom_l(1)=r_wykf(ip,iw,1)+ix
	r_atom_l(2)=r_wykf(ip,iw,2)+iy
	r_atom_l(3)=r_wykf(ip,iw,3)+iz



	call tabbe_tabdil_lattic_cartezian(r_atom_l,r_atom_cz)

	do i=1,3
 	!R_H_elm(i)=r_atom_cz(i)-point_cz(i) !bordar vasel H be atom digar==> jahat bordar F niroo
	R_H_elm(i)=point_cz(i)-r_atom_cz(i)
	!R_H_elm(i)=xyz_H_cartz(ipH,iwH,i)-r_atom_cz(i)
	enddo


	! R_H_elm() is an array to keep the components of vector R=r|H-r|elm

	call tabbe_direct_distance_cartz(point_cz(:),r_atom_cz(:),dis)
	!dis ==> fasele ta atom mored nazar

	if (ix==-1 .and. iy==-1 .and. iz==-1 .and. ip==24) then
!	write(*,*) "+++++++++",point_cz(:)
!	write(*,*) "ff",dis
!	write(*,*) "---------",r_atom_cz(:)
!	write(*,*) iwykf(ip)
!	write(*,*) "ff",ip,iw
	endif


if (dis==0) then
write(*,*) ip,iw
read (*,*) i
endif



	!baar e atom H ra +1 dar nazar gereftim
	!......mohasebe niroo
	do i=1,3
	F(i)=F(i)+ (n_oxidation(ip)/(dis**3) )*R_H_elm(i)
	enddo

45 continue

enddo
enddo
enddo

!write(*,*) "fffffffffffff", F(:)

!write(*,*) "!!!!!!!!!!!!!!!!!!!"
55 continue


enddo !!iw 
enddo !ip



!write (*,'(A10,3f10.5)') "Force...>" , force_H(ipH,iwH,:)


do inpH=1,iwater  !***************************************************************



do inwH=1,2   !har AAb 2 ta H dare     !iwater*2

if (inpH==1) then  !farz kardim mokhtasate aab shomare yek mad nazare
 goto 65  
endif

!write(*,'(2I5,3f12.7)')  inpH,inwH, xyz_H_lattice(inpH,inwH,:)


do ix=-1*N , N
do iy=-1*N , N
do iz=-1*N , N


	r_atom_l(1)=xyz_H_lattice(inpH,inwH,1)+ix
	r_atom_l(2)=xyz_H_lattice(inpH,inwH,2)+iy
	r_atom_l(3)=xyz_H_lattice(inpH,inwH,3)+iz

	call tabbe_tabdil_lattic_cartezian(r_atom_l,r_atom_cz)

	do i=1,3
	!R_H_elm(i)=xyz_H_cartz(ipH,iwH,i)-r_atom_cz(i)
	R_H_elm(i)=point_cz(i)-r_atom_cz(i)
	!R_H_elm(i)=xyz_H_cartz(ipH,iwH,i)-r_atom_cz(i)
	enddo

	call tabbe_direct_distance_cartz(point_cz(:),r_atom_cz(:),dis)



	do i=1,3
	!F(i)=F(i)+ (1/(dis**3) )*R_H_elm(i)
	F(i)=F(i)+ (h_oxidation/(dis**3) )*R_H_elm(i)
	enddo


! write(*,'(f12.7)') dis

60 continue
enddo
enddo
enddo

65 continue
enddo
enddo









75 continue

do i=1,3
FF(i)=F(i)
enddo

!write (*,'(A10,3f10.5)') "........>" , FF(:)
!read (*,*) i



return
end subroutine tabbe2_calculate_forece_on_a_point
!niroo be yek poit
!mokhtasate cartezi migire
subroutine  tabbe2_calculate_forece_on_TIPxP(point_cz,nTIP_indx,nTIP,Nwat,N)
use module_dictionary
use module_structure_symmetry_all
implicit double precision (a-h,o-z)

integer,intent(in):: nTIP,Nwat,N,nTIP_indx
real,dimension(3),intent(in) :: point_cz 
!real,dimension(3),intent(out) :: FF
real,dimension(3):: FF
!ipH : shomare molecool AAb ro neshoon mide
!iwH: shomare H oon molecool AAb

real,dimension(3) :: r1,r2,r_atom_cz,r_atom_l,R_H_elm,F !in F bayad taghir kone hala, moteghayer global beshe ehtemalan
real::dis,dis_max
integer :: ix,iy,iz,ip,iw

!-1.14819   0.14076   5.65789
!point_cz(1)=-1.14819
!point_cz(2)=0.14076
!point_cz(3)=5.65789

dis_max=45000.
F=0.

!write(*,*) "#############", nTIP_indx,nTIP,Nwat



do ip=1,natypes

icont=0
do iwat=1,iwater
 if (lO_wykf(iwat)==1) then
	if (ip==lO_elm(iwat))  then
	icont=1
!	write(*,*) iwat
	endif
 endif
enddo

!write(*,*) ip, icont
!read (*,*) lll
if (icont==1) goto 155


if (lH_ip(ip)/=0 ) goto 155
 

do iw=1,iwykf(ip)
!nemikhahim oxygenhayi ke vojood nadaran hesab shavand
!*****************   agar atomi nabayad hazf beshe inja bayad moshakhas beshe
!******************  agar nabayad chizi hazf beshe bere be khat 225
goto 225
if (ip==3) then
if (iw==11 .or. iw==12) goto 55  
if (iw==15 .or. iw==16) goto 55  
if (iw==19 .or. iw==20) goto 55  
if (iw==23 .or. iw==24) goto 55  
if (iw==27 .or. iw==28) goto 55  
if (iw==31 .or. iw==32) goto 55  

endif
!225 continue

if (ip==lO_elm(ipH)) then     !for  Grimselite 06225/ 40000  !2 ta AMC4000 vojood dare ke ba ham motefavetan, too Output va Dic
if (iw==2 .or. iw==4) goto 55  
				!ein baraye vaghtiye ke partial occupancy darim. 
endif
225 continue

!!******************************************************************************

!write(*,*) "ff", r_wykf(ip,iw,:)


do ix=-1*N , N
do iy=-1*N , N
do iz=-1*N , N
		


	r_atom_l(1)=r_wykf(ip,iw,1)+ix
	r_atom_l(2)=r_wykf(ip,iw,2)+iy
	r_atom_l(3)=r_wykf(ip,iw,3)+iz



	call tabbe_tabdil_lattic_cartezian(r_atom_l,r_atom_cz)

	do i=1,3
 	!R_H_elm(i)=r_atom_cz(i)-point_cz(i) !bordar vasel H be atom digar==> jahat bordar F niroo
	R_H_elm(i)=point_cz(i)-r_atom_cz(i)
	!R_H_elm(i)=xyz_H_cartz(ipH,iwH,i)-r_atom_cz(i)
	enddo


	! R_H_elm() is an array to keep the components of vector R=r|H-r|elm

	call tabbe_direct_distance_cartz(point_cz(:),r_atom_cz(:),dis)
	!dis ==> fasele ta atom mored nazar

	if (ix==-1 .and. iy==-1 .and. iz==-1 .and. ip==24) then
!	write(*,*) "+++++++++",point_cz(:)
!	write(*,*) "ff",dis
!	write(*,*) "---------",r_atom_cz(:)
!	write(*,*) iwykf(ip)
!	write(*,*) "ff",ip,iw
	endif


if (dis==0.) then
!write(*,*) "Dis =0 force=?",ip,iw
!read (*,*) i
endif



	!baar e atom H ra +1 dar nazar gereftim
	!......mohasebe niroo
	if (dis > 0.0 ) then
	do i=1,3
	F(i)=F(i)+ ( (occ(ip,iw)*occ(lO_elm(Nwat),1)*n_oxidation(ip)*delta_q ) /(dis**3) )*R_H_elm(i)
	enddo
	endif

45 continue

enddo
enddo
enddo

!write(*,*) "fffffffffffff", F(:)

!write(*,*) "!!!!!!!!!!!!!!!!!!!"
55 continue


enddo !!iw 
155 continue
enddo !ip



!write(*,'(A10,3f12.6)') "1-----  ",F(:)












! mikhaym nirooye hasel az bar O jabeja shode ro, hesab konim! tavajoh shavad, in faghar baraye Oxygen haye khaej az unit sell hast!
! in momkene mansha eshtebah bashe! choon baraye grimselite, O ha tooye unit cell nistand!!

!Goto 35
do jtip=1,nTIP
do iwat=1,iwater
!do iw=1,lO_wykf(iwat)

do ix=-1*N , N
do iy=-1*N , N
do iz=-1*N , N
		
if (iwat /= Nwat .or. ix /= 0 .or. iy /= 0 .or. iz /=0 ) then !mikhaym nirooye hasel az khodeshoon ro hesab nakonim
!if (iwat /= Nwat ) then

	r_atom_l(1)=ox_shift_latt(iwat,jtip,1)+ix
	r_atom_l(2)=ox_shift_latt(iwat,jtip,2)+iy
	r_atom_l(3)=ox_shift_latt(iwat,jtip,3)+iz

	call tabbe_tabdil_lattic_cartezian(r_atom_l,r_atom_cz)

	do i=1,3
	! R_H_elm(i)=xyz_H_cartz(ipH,iwH,i)-r_atom_cz(i)
	R_H_elm(i)=point_cz(i)-r_atom_cz(i)
	enddo

!	call tabbe_direct_distance_cartz(xyz_H_cartz(ipH,iwH,:),r_atom_cz(:),dis)
	call tabbe_direct_distance_cartz(point_cz(:),r_atom_cz(:),dis)
!if (dis < 0.1) write(*,*) "dis==0, elements", iwat,lO_elm(iwat),lO_wykf(iwat)

!if (dis > 0.1) then   !04/05/2017 in activated

!write(*,*) "dis==0, elements", iwat,lO_elm(iwat),lO_wykf(iwat)
!read (*,*) i
!endif
	!baar e atom H ra +1 dar nazar gereftim
	!......mohasebe niroo
	if (dis > 0.0 ) then
	do i=1,3
	F(i)=F(i)+ (occ(lO_elm(Nwat),1) *occ(lO_elm(iwat),1) *(delta_q**2 )/(dis**3) )*R_H_elm(i)
	enddo
	endif

!endif    !04/05/2017 in activated
endif

enddo
enddo
enddo
!enddo

enddo
enddo !jtip
35 continue

!write(*,'(A10,3f12.6)') "2-----  ",F(:)














!write (*,'(A10,3f10.5)') "Force...>" , force_H(ipH,iwH,:)


do inpH=1,iwater  !***************************************************************



do inwH=1,2   !har AAb 2 ta H dare     !iwater*2

!if (inpH==Nwat) then  
! goto 65  
!endif

!write(*,'(2I5,3f12.7)')  inpH,inwH, xyz_H_lattice(inpH,inwH,:)


do ix=-1*N , N
do iy=-1*N , N
do iz=-1*N , N

if (inpH/=Nwat .or. ix /= 0 .or. iy /= 0 .or. iz /=0 ) then
!if (inpH/=Nwat ) then

	r_atom_l(1)=xyz_H_lattice(inpH,inwH,1)+ix
	r_atom_l(2)=xyz_H_lattice(inpH,inwH,2)+iy
	r_atom_l(3)=xyz_H_lattice(inpH,inwH,3)+iz

	call tabbe_tabdil_lattic_cartezian(r_atom_l,r_atom_cz)

	do i=1,3
	!R_H_elm(i)=xyz_H_cartz(ipH,iwH,i)-r_atom_cz(i)
	R_H_elm(i)=point_cz(i)-r_atom_cz(i)
	!R_H_elm(i)=xyz_H_cartz(ipH,iwH,i)-r_atom_cz(i)
	enddo

	call tabbe_direct_distance_cartz(point_cz(:),r_atom_cz(:),dis)


	if (dis > 0.0 ) then
	do i=1,3
	F(i)=F(i)+ ( occ(lO_elm(Nwat),1) *occ(lO_elm(inpH),1) *(h_oxidation*delta_q) /(dis**3) )*R_H_elm(i)
	enddo
	endif

! write(*,'(f12.7)') dis
endif

60 continue
enddo
enddo
enddo

65 continue
enddo
enddo






!write(*,'(A10,3f12.6)') "3-----  ",F(:)


75 continue

do i=1,3
!FF(i)=F(i)
delta_force(Nwat,nTIP_indx,i)=F(i)
enddo

!write (*,'(A10,3f10.5)') "........>" , FF(:)
!read (*,*) i



return
end subroutine tabbe2_calculate_forece_on_TIPxP
! mokhtasat kartezi ro migire, ye noghte markazi ro ham migire, be markaziyat in noghte, mokhtasat koravi
!noghte aval ro bedast miyare

subroutine tabbe2_cartz_spherical(r_origin,r_xyz,R,thet,phi) !,r_lattic,r_cartezian)
use module_dictionary
use module_structure_symmetry_all
implicit double precision (a-h,o-z)

real, dimension(3),intent(in) :: r_origin,r_xyz
!real,intent(in):: r
real,intent(out):: thet,phi,R !r_lattic, r_cartezian

real, dimension(3)::r_sub  !r_orig_cartz,r_prime
!real:: cs_tht,si_tht,cs_phi,si_phi,rad_tht,rad_phi
real::Rxy,rad_tht,rad_phi

R=0.
Rxy=0.
r_sub=0.



do i=1,3
r_sub(i)=r_xyz(i)-r_origin(i)
R=R+(r_sub(i))**2
enddo

Rxy=R-(r_xyz(3)-r_origin(3))**2 !moalefe Z ro kam kardin



if (R==0.) goto 900

R=sqrt(R)
Rxy=sqrt(Rxy)



thet=acos(r_sub(3)/R)
thet=thet*180/3.141593

if (Rxy>=0.0001) then


phi=acos(abs(r_sub(1)/Rxy))
phi=phi*180/3.141593

if (abs(r_sub(1)/Rxy) >= 1.0000000000000000) phi=0.00000

!write(*,*) "tabbe_cartz_spherical**********",thet,phi
!write(*,'(f21.16)') abs(r_sub(1)/Rxy)
!write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~",acos(1.00000) , acos(1.000001),acos(0.9999999)
if (r_sub(1)>=0) then
	if (r_sub(2) >=0) then
	 phi=phi
	else 
	phi=360-phi
	endif
else
	if (r_sub(2) >=0) then
	 phi=180-phi
	else 
	phi=phi+180
	endif
endif


else !Rxy==0
phi=0.
endif




!write(*,*) "r=",r_sub
!write(*,*) "R==", R
!write(*,*) "Rxy==", Rxy





!write(*,*) "@@@@@@@@@@==",thet,phi

!phi=0

!write(*,*) "tabbe_cartz_spherical*****",thet,phi


900 continue

return
end subroutine tabbe2_cartz_spherical
!mokhtase ye O ro migire,cartezi, ye Theta , Phi, Xi, migire, mokhtasate H haye molecoole AAb ro mide__be Cartezi
!Theta,Phi, Xi darvaghe mokhtasate markaze 2 ta H ast
!t==Theta      f==Phi

subroutine  tabbe2_finding_Hs_by_midelpoint(rO,tht,phi,Xi,rH1,rH2)
use module_dictionary
use module_structure_symmetry_all
implicit double precision (a-h,o-z)

!integer,intent(in):: iwater

real, dimension(3),intent(in) :: rO
real, dimension(3),intent(out) :: rH1,rH2
real, intent(in):: tht,phi,Xi
real :: xd ,xi_cs,xi_si,t,f,x ,rHH, rHH_t    !t==tht, f==phi, x==Xi
real:: RR_H1,tht_H1,Phi_H1,RR_H2,tht_H2,Phi_H2,PI,gama

real, dimension(3) :: r_xi1,r_xi2,r_m
real, dimension(3) :: r_hat,tht_hat,phi_hat,r_mines,r_tmp
integer :: i,itht,iphi
!real:: w_angle=54.
!real, dimension(3) :: RR,point_lat,point_cz,FF,Trq,r_hat,tht_hat,phi_hat,F_tan_cz
!real:: thet,phi,trq_mag,F_mag,t,f,tresh, F_tht,F_phi,F_cz_mag
! open (unit=643, file='./Results/Torque_contour_plot',status='replace', action='write' , IOstat= ierr)

!rHH, rHH_t



rH1=0.
rH2=0.

PI=3.141593
t=tht*3.141593/180.0
f=phi*3.141593/180.0
x=Xi*3.141593/180.0


r_hat(1)= sin(t)*cos(f)
r_hat(2)= sin(t)*sin(f)
r_hat(3)= cos(t)

tht_hat(1)= cos(t)*cos(f)
tht_hat(2)= cos(t)*sin(f)
tht_hat(3)= -1*sin(t)

phi_hat(1)= -1*sin(f)
phi_hat(2)= cos(f)
phi_hat(3)= 0.


xd=d_O_H*tan(w_angle*3.141593/180.0)   !mikhaym az midpoint be in andaze jabeja beshim ta mokhtasate H ro peyda konim
			! 52 nesfe zaviye beyne do H ast==>104


!mokhtasate noghte vasat ...mokhtasate Oxygen bayad behesh ezafe beshe
r_m(1)=d_O_H*r_hat(1)+rO(1)
r_m(2)=d_O_H*r_hat(2)+rO(2)
r_m(3)=d_O_H*r_hat(3)+rO(3)

!write(*,*) rO
!write(*,'(A10,3f12.5)') "rmmm",r_m


!**************************************************
!*** Hydrogen 1
!***************************************************
xi_cs=xd*cos(x)
xi_si=xd*sin(x)

!write(*,*)  "Hydrogen11 ", xd,xi_cs,xi_si
!............................

r_xi1(1)=xi_cs*tht_hat(1)+xi_si*phi_hat(1)
r_xi1(2)=xi_cs*tht_hat(2)+xi_si*phi_hat(2)
r_xi1(3)=xi_cs*tht_hat(3)+xi_si*phi_hat(3)
!...........................


r_xi1(1)=r_xi1(1)+r_m(1)
r_xi1(2)=r_xi1(2)+r_m(2)
r_xi1(3)=r_xi1(3)+r_m(3)




 call tabbe2_cartz_spherical(rO,r_xi1,RR_H1,tht_H1,Phi_H1)
!write(*,*) ,RR_H1,tht_H1,Phi_H1


rH1(1)=d_O_H*sin(tht_H1*3.141593/180.0)*cos(Phi_H1*3.141593/180.0)+rO(1)
rH1(2)=d_O_H*sin(tht_H1*3.141593/180.0)*sin(Phi_H1*3.141593/180.0)+rO(2)
rH1(3)=d_O_H*cos(tht_H1*3.141593/180.0)+rO(3)

!**************************************************
!*** Hydrogen 2
!***************************************************
xi_cs=xd*cos(x+PI)
xi_si=xd*sin(x+PI)

!write(*,*) "Hydrogen22 ", xd,xi_cs,xi_si
!............................

r_xi2(1)=xi_cs*tht_hat(1)+xi_si*phi_hat(1)
r_xi2(2)=xi_cs*tht_hat(2)+xi_si*phi_hat(2)
r_xi2(3)=xi_cs*tht_hat(3)+xi_si*phi_hat(3)
!...........................


r_xi2(1)=r_xi2(1)+r_m(1)
r_xi2(2)=r_xi2(2)+r_m(2)
r_xi2(3)=r_xi2(3)+r_m(3)


!write(*,'(A10,3f12.5)') ">>",r_xi1
!write(*,'(A10,3f12.5)') "<<",r_xi2


 call tabbe2_cartz_spherical(rO,r_xi2,RR_H2,tht_H2,Phi_H2)
!write(*,*) ,RR_H2,tht_H2,Phi_H2,"~~~"


rH2(1)=d_O_H*sin(tht_H2*PI/180.0)*cos(Phi_H2*PI/180.0)+rO(1)
rH2(2)=d_O_H*sin(tht_H2*PI/180.0)*sin(Phi_H2*PI/180.0)+rO(2)
rH2(3)=d_O_H*cos(tht_H2*PI/180.0)+rO(3)


 call tabbe_zaviyeyeab_cartz(rO,rH1,rH2,gama)


!goto 111
!*************************************************
	rHH=0.
	do i=1,3
	r_mines(i)=rH1(i)-rH2(i)  !r1-(r2-R)
	rHH=rHH+r_mines(i)**2
	enddo
rHH=sqrt(rHH)
rHH_t=0.
do i=1,3
rHH_t=rHH_t+ r_mines(i)*phi_hat(i) !zarb dakheli bordar vasl konande H ha ba phi^
enddo
gama=acos(abs(rHH_t/rHH) )

if (abs(rHH_t/rHH)>1.0000000) gama=0.

gama=gama*180.0/3.141593

if (rHH_t/rHH <0 ) gama = 180. - gama

if (gama > 90.0) then
r_tmp(:)=rH1(:)
rH1(:)=rH2(:)
rH2(:)=r_tmp(:)
endif
!############################


!write(*,*) gama
!write(*,*) "&&^&&"


111 continue


return
end subroutine tabbe2_finding_Hs_by_midelpoint
!barname ei baraye tabdil mokhtasate tamam atomha be mokhtasate cartezian
!arraye "n_wykf_oxid" ham inja meghdar dehi mishavad (tariif nakardam felan)
subroutine  tabbe2_make_xyz_cartz_array()
use module_dictionary
use module_structure_symmetry_all
implicit double precision (a-h,o-z)

integer :: ip,iw,ierr,j

ALLOCATE (xyz_wykf_cartz(natypes,iwykf_max,3),STAT= ierr) 
!ALLOCATE (n_wykf_oxid(natypes,iwykf_max),STAT= ierr)
!ALLOCATE (r_wykf(natypes,iwykf_max,3),STAT= nst)


do ip=1,natypes
do iw=1,iwykf(ip)
call tabbe_tabdil_lattic_cartezian(r_wykf(ip,iw,:),xyz_wykf_cartz(ip,iw,:))

!write(*,*)  (r_wykf(ip,iw,j),j=1,3)
!write(*,*)  (xyz_wykf_cartz(ip,iw,j),j=1,3)
!write(*,*) n_oxidation(ip)
!write(*,*) "....."
enddo
enddo





return
end subroutine tabbe2_make_xyz_cartz_array
!Tabe ei ke mokhtasate molecul AAb ro migire, mokhtasate theta phi noghte vasat 2 ta h ro mide va Xi (jahatgiri)
!mokhtasat bayad kartezi bashad
!ma'koos in barname "tabbe2_finding_Hs_by_midelpoint.f90" ast.
subroutine tabbe2_median(rO,rH1,rH2,tht,phi,xi)
use module_dictionary
use module_structure_symmetry_all
implicit double precision (a-h,o-z)

real, dimension(3),intent(in) :: rO,rH1,rH2
real,intent(out)::tht,phi,xi

real,dimension(3) ::tansor,r_plus,r_mines,r_hat,tht_hat,phi_hat
integer:: ii,i,ierr,mosbat
real:: rr,rHH,rHH_t,tt,ff,rHH_p

! open (unit=621, file='./Results/Theta_Phi_path',status='old', action='write' , IOstat= ierr)
	rHH=0.
	do i=1,3
	r_plus(i)=rH1(i)+rH2(i)
	r_mines(i)=rH1(i)-rH2(i)  !r1-(r2-R)
	r_plus(i)=r_plus(i)/2
	rHH=rHH+r_mines(i)**2
	enddo
	rHH=sqrt(rHH) !fasele beyne do noghte
	
	
	
	call tabbe2_cartz_spherical(rO,r_plus,rr,tht,phi)

	tt=tht*3.141593/180
	ff=phi*3.141593/180	

 
! unit vector in spherical coordinate, at midel point
r_hat(1)= sin(tt)*cos(ff)
r_hat(2)= sin(tt)*sin(ff)
r_hat(3)= cos(tt)

tht_hat(1)= cos(tt)*cos(ff)
tht_hat(2)= cos(tt)*sin(ff)
tht_hat(3)= -1*sin(tt)

phi_hat(1)= -1*sin(ff)
phi_hat(2)= cos(ff)
phi_hat(3)= 0.
!........................................

! zaviye Xi, zaviye ba bordar Theta^ ast
rHH_t=0.
rHH_p=0.
do i=1,3
rHH_t=rHH_t+ r_mines(i)*tht_hat(i) !zarb dakheli bordar vasl konande H ha ba theta^
rHH_p=rHH_p+ r_mines(i)*phi_hat(i)
enddo

xi=acos(abs(rHH_t/rHH) )


if (abs(rHH_t/rHH)>1.0000000) then
 if (rHH_t>= 0) then
	 xi=0.
 else
	xi=3.141593
 endif 
endif
xi=xi*180.0/3.141593


if (rHH_p/rHH>=0) then
 mosbat=1
else
 mosbat=-1
endif



!if (rHH_t/rHH <0 ) xi = 180. - xi   !baraye bargashtan be gabl az 5/1/2016 in khat ro faal, baghi khotoot ro pak kon
if (rHH_t/rHH <0 .and. mosbat==1 ) xi = 180. - xi 
if (rHH_t/rHH >0 .and. mosbat==-1 ) xi = 180. - xi




!write(*,*) "zzzzzz", tht,phi,xi
! write( 621,'(2f10.4)') tht,phi

return
end subroutine tabbe2_median



!tabbe ei baraye tashih makan, mikhahim motmaen shavim fasele H bishtar az 1 nemishavad
subroutine  tabbe2_modify_positions()
use module_dictionary
use module_structure_symmetry_all
implicit double precision (a-h,o-z)
real,dimension(3):: RR
!integer,intent(in) :: iwater
real:: rr_mag,bozorg,koochik,tht
real:: rd,Theta,phi,si,cs
integer:: icntrl

do iwat=1,iwater
do iH=1,2

call tabbe2_cartz_spherical(xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),:),xyz_H_cartz(iwat,iH,:),rd,theta,phi)


theta=theta*3.141593/180.0
phi=phi*3.141593/180.0



goto 111

rd=0.
do i=1,3
RR(i)=xyz_H_cartz(iwat,iH,i)- xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),i)
rd= rd+ RR(i)**2
enddo

rd=sqrt(rd)
theta=Acos(RR(3)/rd)

!if (iwat==1 .and. iH==1) write(*,*) theta
si=sin(theta)



phi=Acos( abs(RR(1))/(rd*si) )

if (RR(1)<0 .and. RR(2) >=0 ) phi=3.141593-phi  !Area2
if (RR(1)<0 .and. RR(2) < 0 ) phi=3.141593+phi  !Area3
if (RR(1)>=0 .and. RR(2) < 0 ) phi=2*3.141593-phi  !Area4


 
!write(*,*) "phi=  ",phi


!phi=Asin(-1*RR(2)/(rd*si) )
!write(*,*) phi*180/3.141593
!write(*,*) ""

!write(*,*) xyz_H_cartz(iwat,iH,:)
!write(*,*) xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),:)
!write(*,*) "~~~~~~~~~~~~", iwat,iH

111 continue

si=sin(theta)
RR(1)= si* cos(phi)*d_O_H
RR(2)= si* sin(phi) *d_O_H
RR(3)= cos(theta)*d_O_H
!write(*,*) "~~>",RR(:)

!xyz_H_cartz(iwat,iH,1)

do i=1,3
xyz_H_cartz(iwat,iH,i)=RR(i)+ xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),i)
enddo



enddo  !iw
!call tabbe_zaviyeyeab_cartz(xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),:),xyz_H_cartz(iwat,1,:),xyz_H_cartz(iwat,2,:),tht)
!write(*,*) "water zaviye   " , tht
enddo







return
end subroutine tabbe2_modify_positions



subroutine tabbe2_move_by_symmetry(r_ref,r_comper,H_ini,H_trans,is)
use module_dictionary
use module_structure_symmetry_all
implicit double precision (a-h,o-z)

real, dimension(3),intent(in) :: r_ref,r_comper,H_ini
real, dimension(3),intent(out) :: H_trans
integer, intent(in):: is
real, dimension(3) :: r_l,h_tmp,h_tmp2,rr_mog
real:: d_x,d_y,d_z,dis,sss
integer:: istats,mlt,nal,it_rf,mlt_rf,icntr=0
integer:: ii,ix,iy,iz,ix_rf=0,iy_rf=0,iz_rf=0
nal=multil(is)*nshiftl
allocate (hydrogen_coords(nal,3),STAT=istats)



r_l=r_comper

do it=1,multil(is)    !nshiftl considered at following	

d_x= ssyml(1,1,it,is) * r_l(1)+ssyml(1,2,it,is) * r_l(2)+ssyml(1,3,it,is) * r_l(3)
d_y= ssyml(2,1,it,is) * r_l(1)+ssyml(2,2,it,is) * r_l(2)+ssyml(2,3,it,is) * r_l(3)
d_z= ssyml(3,1,it,is) * r_l(1)+ssyml(3,2,it,is) * r_l(2)+ssyml(3,3,it,is) * r_l(3)


d_x=d_x +ttaul(1,it,is)
d_y=d_y +ttaul(2,it,is)
d_z=d_z +ttaul(3,it,is)


icntr=0
do ii=1,nshiftl
mlt=ii-1   
hydrogen_coords(it+(mlt*multil(is)),1)=d_x+vshiftl(1,ii)
hydrogen_coords(it+(mlt*multil(is)),2)=d_y+vshiftl(2,ii)
hydrogen_coords(it+(mlt*multil(is)),3)=d_z+vshiftl(3,ii)


do ix=-1,1
do iy=-1,1
do iz=-1,1
rr_mog(1)=hydrogen_coords(it+(mlt*multil(is)),1)+ix
rr_mog(2)=hydrogen_coords(it+(mlt*multil(is)),2)+iy
rr_mog(3)=hydrogen_coords(it+(mlt*multil(is)),3)+iz

	if (r_ref(1)<=rr_mog(1)+0.001 .and. r_ref(1)>=rr_mog(1)-0.001) then
	if (r_ref(2)<=rr_mog(2)+0.001 .and. r_ref(2)>=rr_mog(2)-0.001) then	
	if (r_ref(3)<=rr_mog(3)+0.001 .and. r_ref(3)>=rr_mog(3)-0.001) then	
		icntr=1
		it_rf=it
		mlt_rf=mlt
		ix_rf=ix
		iy_rf=iy
		iz_rf=iz
		exit
	endif
	endif
	endif

enddo
enddo
enddo

110 continue

enddo !ii

if (icntr==1) exit
enddo



if (icntr==0) then
write(*,*) "There is a Problem at Symmetry operators"
stop
endif



r_l=H_ini

d_x= ssyml(1,1,it_rf,is) * r_l(1)+ssyml(1,2,it_rf,is) * r_l(2)+ssyml(1,3,it_rf,is) * r_l(3)
d_y= ssyml(2,1,it_rf,is) * r_l(1)+ssyml(2,2,it_rf,is) * r_l(2)+ssyml(2,3,it_rf,is) * r_l(3)
d_z= ssyml(3,1,it_rf,is) * r_l(1)+ssyml(3,2,it_rf,is) * r_l(2)+ssyml(3,3,it_rf,is) * r_l(3)



d_x=d_x +ttaul(1,it_rf,is)
d_y=d_y +ttaul(2,it_rf,is)
d_z=d_z +ttaul(3,it_rf,is)


ii=mlt_rf+1
h_tmp(1)=d_x+vshiftl(1,ii)+ix_rf
h_tmp(2)=d_y+vshiftl(2,ii)+iy_rf
h_tmp(3)=d_z+vshiftl(3,ii)+iz_rf


call tabbe_direct_distance(r_ref,h_tmp,dis)



if (dis>=d_O_H-0.05*d_O_H .and. dis<=d_O_H+0.05*d_O_H ) then
	H_trans=h_tmp
else
 do ix=-1,1
 do iy=-1,1
 do iz=-1,1
	h_tmp2(1)=h_tmp(1)+ix
	h_tmp2(2)=h_tmp(2)+iy
	h_tmp2(3)=h_tmp(3)+iz
	
	call tabbe_direct_distance(r_ref,h_tmp2,dis)
	if (dis>=d_O_H-0.07*d_O_H .and. dis<=d_O_H+0.07*d_O_H ) then
	H_trans=h_tmp2
	exit
	endif
 enddo
 enddo
 enddo
endif


return
end subroutine tabbe2_move_by_symmetry
!in barname asliye ke charkhesh ro baraye ma anjam mide
!dar tarikh 10/4/2016 taghiresh dadam, goto 111 ro negah kon
!**NOKTE: shive charkhesh ro be Quatronium taghir dadam
!dar format ghadim, del_phi jabejaei bood dar vaghe, ke 0.01 dar nazar migereftim
!dar shive jadid, del_phi vaghean zaviye charkheshe, ke bayad 0.5 bashe, moadele 0.01 jabejaei
subroutine  tabbe2_rotation(rO,rH,rH_prim,Trq,del_phi)
use module_dictionary
use module_structure_symmetry_all
implicit double precision (a-h,o-z)

real, dimension(3),intent(in)::rO,rH,Trq
real,intent(in):: del_phi
real, dimension(3),intent(out)::rH_prim
real, dimension(3) :: r_prime,del_r
real, dimension(3) :: RR,torq,rot_axis
real:: torq_mag
real:: RR_mag,rp_mag
integer:: icnt
real, dimension (4) :: rot4vectors,quat
!real,dimension(3):: quat

do i=1,3
RR(i)=rH(i)-rO(i)
enddo


icnt=0
!do kk=1,300

!write(*,'("RO--",3f12.6)') rO
!write(*,'("RH--",3f12.6)') rH
!write(*,'("Torque--",3f12.6)') Trq
!write(*,'("del--",f12.6)') del_phi

!....Cross product torq=r x F
!torq(1)= RR(2)*F(3)-RR(3)*F(2)
!torq(2)= RR(3)*F(1)-RR(1)*F(3)
!torq(3)= RR(1)*F(2)-RR(2)*F(1)
!........................
torq(1)=Trq(1)
torq(2)=Trq(2)
torq(3)=Trq(3)

torq_mag=torq(1)**2+torq(2)**2+torq(3)**2
torq_mag=sqrt(torq_mag)
if (torq_mag==0.) then
!write(*,*) "torq====0."
!read(*,*) icnt
icnt=1
goto 100
endif


do i=1,3
rot_axis(i)=torq(i)/torq_mag    !rotation axis
enddo


goto 111
!........finding delta_r = (r * n ) del_phi
del_r(1)= (RR(2)*rot_axis(3)-RR(3)*rot_axis(2)) * del_phi
del_r(2)= (RR(3)*rot_axis(1)-RR(1)*rot_axis(3)) * del_phi
del_r(3)= (RR(1)*rot_axis(2)-RR(2)*rot_axis(1)) * del_phi
!....................................

do i=1,3
r_prime(i)=RR(i)-del_r(i)
enddo

rp_mag=r_prime(1)**2+r_prime(2)**2+r_prime(3)**2
RR_mag=RR(1)**2+RR(2)**2+RR(3)**2

rp_mag=sqrt(rp_mag)
RR_mag= sqrt(RR_mag)



!write(*,*) "......."
!write(*,'(3f10.5)') RR
!write(*,'(3f10.5)') r_prime
!write(*,*) RR_mag,rp_mag
!write(*,'(3f10.5)') torq
!write(*,'(3f10.5)') rot_axis
!write(*,'(f10.5)') torq_mag

do i=1,3
rH_prim(i)=r_prime(i)+rO(i)
enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
111 continue

!goto 100
rot4vectors(1)=del_phi*3.141593/180.0
do ii=1,3
rot4vectors(1+ii)=rot_axis(ii)
enddo
!rot4vectors(1)=180.*3.141593/180.0
!rot4vectors(2)=0.
!rot4vectors(3)=0.
!rot4vectors(4)=0.

!RR(1)=0.
!RR(2)=1.
!RR(3)=0.

call qxqml2(quat,RR,rot4vectors)
do jj=2,4
!rH_prim(jj-1)=quat(jj)+rO(jj-1)
quat(jj)=quat(jj)+rO(jj-1)
enddo


!write(*,'("my----    ",3f10.5)') rH_prim
!write(*,'("quat--    ",3f10.5)') (quat(k),k=2,4)
!write(*,*) ""


!	call tabbe_zaviyeyeab_cartz(xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),:),xyz_H_cartz(iwat,2,:),rH_prim,dis1)!@@@@@@

!write(*,'("axis--",3f12.6)') rot_axis(2),rot_axis(3),rot_axis(4)
!write(*,'("r_Rot--",3f12.6)') rH_prim


do jj=2,4
rH_prim(jj-1)=quat(jj)       !+rO(jj-1)
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
100 continue
!RR(i)=rH(i)
if (icnt==1) then
rH_prim(:)=rH(:)
endif
!enddo !KK

return
end subroutine tabbe2_rotation
!barname'i baraye mokhtasat dehi ebtedayi be H
!bayad joori taghir kone ke masalan mokhtasat O ha ro az jayi bekhone ya ...
!felan mokhtasat O behesh dade mishe

!04/16/2016
!in barname ro mikham joori taghir bedam ke faghat baraye molecule aval, Thet, Phi , Xi ro begire
!mokhtasat bede


!r_wykf(lO_elm(iwat),lO_wykf(iwat),:)  Mokhtasate latise molecoolhaye Oxygen
! xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),:) Mokhtasate cartezi Molecoolhaye O


subroutine  tabbe2_set_initial_H()
use module_dictionary
use module_structure_symmetry_all
implicit double precision (a-h,o-z)

real :: tht1,tht2,phi1,phi2,dis,dis1,dis2,xrand
real :: tt,pph,gama,psi
integer :: ierr,itht,iphi,icntrl,ict
!integer:: ndim
real,dimension(3) :: r_l,r_H,H1_lattice,H2_lattice,rr1,rr2
integer:: zp,zb,ict2
real:: thet,phi,xi
!integer, allocatable, dimension (:) :: iseed,iseed_first,iseed_last ! seed array

zp=104
zb=106


!integer:: ndim
!integer, allocatable, dimension (:) :: iseed,iseed_first,iseed_last ! seed array
!call random_seed(size=ndim)             ! seeding   
!allocate(iseed(ndim))
!allocate(iseed_first(ndim))
!allocate(iseed_last(ndim))
!iseed= Seed
!deallocate(iseed)
!deallocate(iseed_first)
!deallocate(iseed_last)



do iwat=1,iwater
call random_seed(put=iseed)
call random_seed(get=iseed)
call random_number(xrand)

call random_number (tht1)
tht1=tht1*180.

call random_number (phi1)
phi1=phi1*360.

call random_number (psi)
psi=psi*180.


call random_seed(get=iseed_last) ! Final seed if desired
iseed = iseed_last




!print *, "Random orientation", tht1,phi1,psi


!call tabbe2_finding_Hs_by_midelpoint(xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),:),tht1,phi1,psi
!,xyz_H_cartz(iwat,1,:),xyz_H_cartz(iwat,2,:))
call tabbe2_finding_Hs_by_midelpoint(xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),:),tht1,phi1,psi,rr1,rr2)
xyz_H_cartz(iwat,1,:)=rr1
xyz_H_cartz(iwat,2,:)=rr2
call tabbe_zaviyeyeab_cartz(xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),:),xyz_H_cartz(iwat,1,:),xyz_H_cartz(iwat,2,:),gama)
call tabbe_direct_distance_cartz(xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),:),xyz_H_cartz(iwat,1,:),dis1)
call tabbe_direct_distance_cartz(xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),:),xyz_H_cartz(iwat,2,:),dis2)
!write(*,*) iwat,gama,dis1,dis2

enddo !iwater





ict=0

if (ict /= 0 ) then  !start from a known orientation
thet=50
phi=100
xi=90
!call tabbe2_finding_Hs_by_midelpoint(xyz_wykf_cartz(lO_elm(ict),lO_wykf(ict),:)
!,thet,phi,xi,xyz_H_cartz(ict,1,:),xyz_H_cartz(ict,2,:))
call tabbe2_finding_Hs_by_midelpoint(xyz_wykf_cartz(lO_elm(ict),lO_wykf(ict),:),thet,phi,xi,rr1,rr2)
xyz_H_cartz(iwat,1,:)=rr1
xyz_H_cartz(iwat,2,:)=rr2
write(*,*) "Initial Position for iwat== " , ict, " has been set"
write(*,'(3f12.6)') xyz_H_cartz(ict,1,:)
write(*,'(3f12.6)') xyz_H_cartz(ict,2,:)
write(*,*) "##$$##"
call tabbe2_median(xyz_wykf_cartz(lO_elm(ict),lO_wykf(ict),:),xyz_H_cartz(ict,1,:),xyz_H_cartz(ict,2,:),thet,phi,xi)
write(*,*) "Test the preset angles value: ",thet,phi,xi




endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*) "Random initial positions for hydrogen atoms:"
do i=1,iwater
do j=1,2
call tabbe_tabdil_cartezian_lattic(xyz_H_cartz(i,j,:),xyz_H_lattice(i,j,:))
write(*,'(3f8.4)') xyz_H_lattice(i,j,:)
enddo
enddo
write(*,*) ""
write(*,*) "-----------------------------------------------------"


!deallocate(iseed)
!deallocate(iseed_first)
!deallocate(iseed_last)


return
end subroutine tabbe2_set_initial_H




subroutine tabbe2_spherical_cartz(r_center,R,thet,phi,r_xyz) 
use module_dictionary
use module_structure_symmetry_all
implicit double precision (a-h,o-z)

real, dimension(3),intent(out) :: r_xyz
real, dimension(3),intent(in) ::r_center
!real,intent(in):: r
real,intent(in):: thet,phi,R !r_lattic, r_cartezian

real, dimension(3)::r_sub  !r_orig_cartz,r_prime
real::Rxy,rad_tht,rad_phi
real :: xx,yy,zz,tt,ph

!R=0.
!Rxy=0.
!r_sub=0.


tt=thet*3.141593/180.
ph=phi*3.141593/180.


xx= R*sin(tt)*cos(ph)
yy= R*sin(tt)*sin(ph)
zz= R*cos(tt)


r_xyz(1)=xx + r_center(1)
r_xyz(2)=yy + r_center(2)
r_xyz(3)=zz + r_center(3)



return
end subroutine tabbe2_spherical_cartz
! In tabe, shekle molekool Aab ro, ba asas modeli ke darim, misaze
! masalan, bar asas model TIP4P zaviye ro 104 migire va jebejayi bar Ox ro 0.15A
! be har hal tabe asli dar shekl dehi model hast

!*** In 2 mokhtasat ke jabejayi bar OX hastand ro dar akhar bayad moshakhas konim
!ox_shift_cz(iwat,nTIP,:)
!ox_shift_latt(iwat,nTIP,:)

!r_wykf(lO_elm(iwat),lO_wykf(iwat),:)  Mokhtasate latise molecoolhaye Oxygen
! xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),:) Mokhtasate cartezi Molecoolhaye Oxygen
!xyz_H_cartz(iwat,1,:) mokhtasat H, moolecool iwat
!xyz_H_lattice(iwat,iH,:)

! DELTA_DISTANCE: fasee bar Ox ta makan haghighi Ox
!delta_angle: zaviye ei ke barhaye majazi ba ham misazan, dar TIP5P model

subroutine  tabbe2_TIPxP(nTIP)
use module_dictionary
use module_structure_symmetry_all
implicit double precision (a-h,o-z)

!integer, intent(in):: iwater
integer,intent(in) :: nTIP
integer:: ierr
real,dimension(3) :: r_plus,r_plus_hat,r_vert,r_vert_hat,roh1,roh2,r_plus_inv
real:: r_plus_mag,dis,mag_test,r_vert_mag,ang,dd,tht
integer:: i,j,k,mm,jtip
!real:: tm,sm,tt



! mokhtase aval, shomare molecool Water ro neshoon mide

!write(*,*) "********TABBE TIPXP*******"

ang= 3.141593* delta_angle/(2. * 180.)
dd=delta_distance
!write(*,*) dd,ang


if (nTIP==1) then  !TIP4P model
do iwat=1,iwater
	r_plus_mag=0.
	do i=1,3
	r_plus(i)=xyz_H_cartz(iwat,1,i)+xyz_H_cartz(iwat,2,i)
	r_plus(i)=r_plus(i)/2
	r_plus(i)=r_plus(i)-xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),i) 
	r_plus_mag=r_plus_mag+r_plus(i)**2
	enddo
	r_plus_mag=sqrt(r_plus_mag) 


	do j=1,3
	r_plus_hat(j)=r_plus(j)/r_plus_mag
	enddo

	do k=1,3
	ox_shift_cz(iwat,nTIP,k)=delta_distance* r_plus_hat(k)+ xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),k)
	enddo

!	call tabbe_direct_distance_cartz(xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),:),ox_shift_cz(iwat,nTIP,:),dis)
!	WRITE(*,'(3F12.6)') xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),:)
!	WRITE(*,'(3F12.6)') ox_shift_cz(iwat,nTIP,:)
!	write(*,*) dis
!	read (*,*) mm

	call tabbe_tabdil_cartezian_lattic(ox_shift_cz(iwat,nTIP,:),ox_shift_latt(iwat,nTIP,:))

enddo

else if (nTIP==2) then

 do iwat=1,iwater

	r_plus_mag=0.
	do i=1,3
	r_plus(i)=xyz_H_cartz(iwat,1,i)+xyz_H_cartz(iwat,2,i)
	r_plus(i)=r_plus(i)/2
	r_plus(i)=r_plus(i)-xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),i) 
	r_plus_mag=r_plus_mag+r_plus(i)**2
	enddo
	r_plus_mag=sqrt(r_plus_mag) 

	do j=1,3
	r_plus_hat(j)=r_plus(j)/r_plus_mag
	enddo 
!Same as nTIP=1, but now we whant to change the direction of r_plus_hat
	r_plus_inv(:)=-1*r_plus_hat(:)

!need to find perpendicular direction to H-O-H surface
    roh1(:)=xyz_H_cartz(iwat,1,:)-xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),:)
    roh2(:)=xyz_H_cartz(iwat,2,:)-xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),:)

    r_vert(1)=roh1(2)*roh2(3)-roh1(3)*roh2(2)
    r_vert(2)=roh1(3)*roh2(1)-roh1(1)*roh2(3)
    r_vert(3)=roh1(1)*roh2(2)-roh1(2)*roh2(1)

!	r_vert(:)-
    r_vert_mag=sqrt(r_vert(1)**2 + r_vert(2)**2 + r_vert(3)**2 )
    r_vert_hat(:)= r_vert(:)/r_vert_mag

  do k=1,3
 	ox_shift_cz(iwat,1,k)=dd*cos(ang)*r_plus_inv(k)+ dd*sin(ang)* r_vert_hat(k)+xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),k)
	ox_shift_cz(iwat,2,k)=dd*cos(ang)*r_plus_inv(k)- dd*sin(ang)* r_vert_hat(k)+xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),k)
  enddo

  do jtip=1,nTIP
	call tabbe_tabdil_cartezian_lattic(ox_shift_cz(iwat,jtip,:),ox_shift_latt(iwat,jtip,:))
  enddo

 enddo
endif

!write(*,'(3f12.6)') xyz_wykf_cartz(lO_elm(5),lO_wykf(5),:)
!write(*,'(3f12.6)') ox_shift_cz(1,2,:)
!write(*,'(3f12.6)') ox_shift_cz(1,2,:)
!call tabbe_zaviyeyeab_cartz(xyz_wykf_cartz(lO_elm(5),lO_wykf(5),:),ox_shift_cz(5,1,:),ox_shift_cz(5,2,:),tht)
!xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),i)
!write(*,'(A15,f16.10)') "z-------", tht
!	call tabbe_tabdil_cartezian_lattic(xyz_wykf_cartz(lO_elm(5),lO_wykf(5),:),roh1)
!	write(*,'(A7,3f12.6)') " O  0",roh1
 !	write(*,'(A7,3f12.6)') " He  0",xyz_H_lattice(5,1,:)!
!	write(*,'(A7,3f12.6)') " He  0",xyz_H_lattice(5,2,:)
!	write(*,*) ""
!	write(*,'(A7,3f12.6)') " c  0",ox_shift_latt(5,1,:)
!	write(*,'(A7,3f12.6)') " c  0",ox_shift_latt(5,2,:)
!	read (*,*) mm

!sm=1.
!do i=5,90,5
!tm=real(i)*3.1415/180.
!tt=cos(tm)*2.*3.1415/0.0981
!write(*,*) tt
!sm=sm+tt
!enddo
!write(*,*) "ineee" ,sm
!read(*,*) mm


return
end subroutine tabbe2_TIPxP
!tabbe2_Torque
subroutine  tabbe2_Torque(iwat,iH,nTIP)
use module_dictionary
use module_structure_symmetry_all
implicit double precision (a-h,o-z)

!integer,intent(in):: iwater

integer,intent(in) :: iwat, iH
integer:: i,jtip
real, dimension(3) :: RR
real,dimension (nTIP,3)::dRR
!xyz_H_cartz(iwat,iH,i)
!-2.31577   1.06908   5.51181



!do iwat=1,iwater
!do iH=1,2

!write(*,*) "***",lO_elm(iwat),lO_wykf(iwat)
!write(*,*) xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),:)

do i=1,3
RR(i)=xyz_H_cartz(iwat,iH,i)-xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),i) !riidemal shod, vali "xyz_wykf_cartz" !mokhtasate oxygen ast.
!do jtip=1,nTIP
!dRR(jtip,i)=  ox_shift_cz(iwat,nTIP,i)-xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),i)	
!enddo
enddo


! T= RR*F 
torque_H(iwat,iH,1)= RR(2)*force_H(iwat,iH,3)-RR(3)*force_H(iwat,iH,2)
torque_H(iwat,iH,2)= RR(3)*force_H(iwat,iH,1)-RR(1)*force_H(iwat,iH,3)
torque_H(iwat,iH,3)= RR(1)*force_H(iwat,iH,2)-RR(2)*force_H(iwat,iH,1)
!torq(2)= RR(3)*F(1)-RR(1)*F(3)
!torq(3)= RR(1)*F(2)-RR(2)*F(1)


!delta_torque(iwat,nTIP,:)
!delta_force(iwat,jtip,:)

!do jtip=1,nTIP
!delta_torque(iwat,jtip,1)= dRR(jtip,2)*delta_force(iwat,jtip,3)-dRR(jtip,3)*delta_force(iwat,jtip,2)
!delta_torque(iwat,jtip,2)= dRR(jtip,3)*delta_force(iwat,jtip,1)-dRR(jtip,1)*delta_force(iwat,jtip,3)
!delta_torque(iwat,jtip,3)= dRR(jtip,1)*delta_force(iwat,jtip,2)-dRR(jtip,2)*delta_force(iwat,jtip,1)
!enddo


!enddo
!enddo



!write(*,'(A15,3f12.7)') "Torq..>  " ,torque_H(1,1,1),torque_H(1,1,2),torque_H(1,1,3)




return
end subroutine tabbe2_Torque
!tabbe2_Torque be yek point
! mokhtasate cartezi
subroutine  tabbe2_Torque_on_point(r_O,point_cz,FF,Trq)
use module_dictionary
use module_structure_symmetry_all
implicit double precision (a-h,o-z)

!integer,intent(in):: iwater

integer :: iwat,j , k, iH
real, dimension(3) :: RR
real, dimension(3),intent(in)::r_O,point_cz,FF
real, dimension(3),intent(out)::Trq

!do iwat=1,iwater
!do iH=1,2

!write(*,*) "***",lO_elm(iwat),lO_wykf(iwat)
!write(*,*) xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),:)





do i=1,3
RR(i)=point_cz(i)-r_O(i) 
								
enddo


! T= RR*F 
Trq(1)= RR(2)*FF(3)-RR(3)*FF(2)
Trq(2)= RR(3)*FF(1)-RR(1)*FF(3)
Trq(3)= RR(1)*FF(2)-RR(2)*FF(1)
!torq(2)= RR(3)*F(1)-RR(1)*F(3)
!torq(3)= RR(1)*F(2)-RR(2)*F(1)



!enddo
!enddo








return
end subroutine tabbe2_Torque_on_point
!tabbe2_Torque
subroutine  tabbe2_Torque_TIPxP(iwat,iH,nTIP)
use module_dictionary
use module_structure_symmetry_all
implicit double precision (a-h,o-z)

!integer,intent(in):: iwater

integer,intent(in) :: iwat, iH
integer:: i,jtip
real, dimension(3) :: RR
real,dimension (nTIP,3)::dRR
!xyz_H_cartz(iwat,iH,i)
!-2.31577   1.06908   5.51181



!do iwat=1,iwater
!do iH=1,2

!write(*,*) "***",lO_elm(iwat),lO_wykf(iwat)
!write(*,*) xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),:)


!RR(i)=xyz_H_cartz(iwat,iH,i)-xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),i) !riidemal shod, vali "xyz_wykf_cartz" !mokhtasate oxygen ast.
do jtip=1,nTIP
do i=1,3
!dRR(jtip,i)=  ox_shift_cz(iwat,nTIP,i)-xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),i)	
dRR(jtip,i)=  ox_shift_cz(iwat,jtip,i)-xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),i)	
enddo
enddo




!delta_torque(iwat,jtip,:)
!delta_force(iwat,jtip,:)

do jtip=1,nTIP
delta_torque(iwat,jtip,1)= dRR(jtip,2)*delta_force(iwat,jtip,3)-dRR(jtip,3)*delta_force(iwat,jtip,2)
delta_torque(iwat,jtip,2)= dRR(jtip,3)*delta_force(iwat,jtip,1)-dRR(jtip,1)*delta_force(iwat,jtip,3)
delta_torque(iwat,jtip,3)= dRR(jtip,1)*delta_force(iwat,jtip,2)-dRR(jtip,2)*delta_force(iwat,jtip,1)
enddo


!enddo
!enddo



!write(*,'(A15,3f12.7)') "Torq..>  " ,torque_H(1,1,1),torque_H(1,1,2),torque_H(1,1,3)




return
end subroutine tabbe2_Torque_TIPxP
!tabe'e peyda kardan fasele do noghte dar fazaye dekarti
subroutine tabbe_direct_distance_cartz(r1,r2,dis)
use module_dictionary
use module_structure_symmetry_all
implicit double precision (a-h,o-z)


real, dimension(3),intent(in) :: r1,r2
real,intent(out)::dis

real,dimension(3) ::r_subtract
integer:: ii

	dis=0.

	r_subtract(1)=r1(1)-r2(1) !r1-(r2-R)
	r_subtract(2)=r1(2)-r2(2)
	r_subtract(3)=r1(3)-r2(3)


	do ii=1,3
	dis=dis+r_subtract(ii)*r_subtract(ii)
	enddo


	dis=sqrt(dis)





return
end subroutine tabbe_direct_distance_cartz

!mokhtasate do noghte dar fazaye lattice ro migire, fasele ro mide
subroutine tabbe_direct_distance(r1,r2,dis)
use module_dictionary
use module_structure_symmetry_all
implicit double precision (a-h,o-z)

real, dimension(3),intent(in) :: r1,r2
real,intent(out)::dis

real,dimension(3) ::tansor,r_subtract
integer:: ii

	r_subtract(1)=r1(1)-r2(1) !r1-(r2-R)
	r_subtract(2)=r1(2)-r2(2)
	r_subtract(3)=r1(3)-r2(3)

	!amin_dis=0.
	tansor(1)=0.0
	tansor(2)=0.0
	tansor(3)=0.0
	dis=0.

!
	do ii=1,3
	tansor(1)=tansor(1)+T_metric(1,ii)*r_subtract(ii)
	tansor(2)=tansor(2)+T_metric(2,ii)*r_subtract(ii)
	tansor(3)=tansor(3)+T_metric(3,ii)*r_subtract(ii)
	enddo

	do ii=1,3
	dis=dis+r_subtract(ii)*tansor(ii)
	enddo

	dis=sqrt(dis)


!dis^2 = rTr

return
end subroutine tabbe_direct_distance
!mokhtasate do noghte dar fazaye lattice ro migire, fasele ro mide
!r2 is constant
!r1 is a lattice point which can varies by symmetry
subroutine tabbe_direct_nearest_distance(r1,r2,dis,i,j,k)
use module_dictionary
use module_structure_symmetry_all
implicit double precision (a-h,o-z)

real, dimension(3),intent(in) :: r1,r2
real,intent(out)::dis
integer,intent(out) :: i,j,k
real:: tmp_dis,dis2

real,dimension(3) ::tansor,r_subtract,r1p
integer:: ii,ix,iy,iz


dis=1000.0

ix=0
iy=0
ix=0



r1p=r1

!write(*,*) r1,r1p

do ix= -1,1
	r1p(1)=r1(1)+ix
do iy= -1,1
	r1p(2)=r1(2)+iy
do iz= -1,1
	r1p(3)=r1(3)+iz
	
	
	

!write(*,*) ix,iy,iz
!write(*,'(3f12.6,A2)') r1p(:),"!!"

	r_subtract(1)=r1p(1)-r2(1) !r1-(r2-R)
	r_subtract(2)=r1p(2)-r2(2)
	r_subtract(3)=r1p(3)-r2(3)

	!amin_dis=0.
	tansor(1)=0.0
	tansor(2)=0.0
	tansor(3)=0.0
	tmp_dis=0.

!
	do ii=1,3
	tansor(1)=tansor(1)+T_metric(1,ii)*r_subtract(ii)
	tansor(2)=tansor(2)+T_metric(2,ii)*r_subtract(ii)
	tansor(3)=tansor(3)+T_metric(3,ii)*r_subtract(ii)
	enddo

	do ii=1,3
	tmp_dis=tmp_dis+r_subtract(ii)*tansor(ii)
	enddo

	tmp_dis=sqrt(tmp_dis)
	
	!write(*,'(2f12.6)') tmp_dis,dis
	if (dis>tmp_dis) then
		 dis=tmp_dis
		i=ix
		j=iy
		k=iz	
	endif

enddo
enddo
enddo

!write(*,*) dis

return
end subroutine tabbe_direct_nearest_distance
subroutine tabbe_tabdil_cartezian_lattic(r_cartezian,r_lattic)
use module_dictionary
use module_structure_symmetry_all
implicit double precision (a-h,o-z)

real,dimension(3),intent(out) :: r_lattic
real,dimension(3),intent(in):: r_cartezian
real,dimension(3,3):: invers_al

r_lattic=0.

call inverse(aL_matrix,invers_al,3)

!do i=1,3
!write(*,*) (invers_al(i,j),j=1,3)
!enddo

do k=1,3
r_lattic(1)= r_lattic(1)+ invers_al(1,k) *r_cartezian(k)    !===   invers_al_matrix(3*3) x r_cartezian(3*1)
r_lattic(2)= r_lattic(2)+ invers_al(2,k) *r_cartezian(k)
r_lattic(3)= r_lattic(3)+ invers_al(3,k) *r_cartezian(k)
!write(*,*) aL_matrix(k,1)

enddo






return
end subroutine tabbe_tabdil_cartezian_lattic

 subroutine inverse(aa,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer:: n
real,dimension(3,3), intent(in) :: aa
real,dimension(3,3), intent(out) :: c
double precision:: a(n,n)
double precision:: L(n,n), U(n,n), b(n), d(n), x(n)
double precision:: coeff
integer:: i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
a=aa
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse
subroutine tabbe_tabdil_lattic_cartezian(r_lattic,r_cartezian)
use module_dictionary
use module_structure_symmetry_all
implicit double precision (a-h,o-z)

real,dimension(3),intent(in) :: r_lattic
real,dimension(3),intent(out):: r_cartezian

integer :: i,j,k

r_cartezian=0.
do k=1,3
r_cartezian(1)= r_cartezian(1)+ aL_matrix(1,k) *r_lattic(k)    !===   al_matrix(3*3) x r_lattice(3*1)
r_cartezian(2)= r_cartezian(2)+ aL_matrix(2,k) *r_lattic(k)
r_cartezian(3)= r_cartezian(3)+ aL_matrix(3,k) *r_lattic(k)
!write(*,*) aL_matrix(k,1)

enddo

!write(*,*) r_cartezian

return
end subroutine tabbe_tabdil_lattic_cartezian

subroutine tabbe_zaviyeyeab_cartz(r_cen_cartezian,r_1_cartezian,r_2_cartezian,theta)
use module_dictionary
use module_structure_symmetry_all
!use module_structure
!use module_structure_symmetry_local
implicit double precision (a-h,o-z)

!real,dimension(3),intent(in) :: r_center,r_1,r_2
real,dimension(3),intent(in) :: r_1_cartezian,r_2_cartezian,r_cen_cartezian
real, intent(out) :: theta
real,dimension(3) :: R_hat_1,R_hat_2!,r_1_cartezian,r_2_cartezian,r_cen_cartezian
real :: s_bar1,s_bar2, dot_pro
integer :: i,j ,k

!write(*,*) "r_cen" ,r_cen_cartezian
!write(*,*) "r_1" ,r_1_cartezian
!write(*,*) "r_2" ,r_2_cartezian





 s_bar1=0.
 s_bar2=0.
  do i=1,3
   R_hat_1(i)=r_1_cartezian(i) - r_cen_cartezian(i)
   R_hat_2(i)=r_2_cartezian(i) - r_cen_cartezian(i)

   s_bar1=s_bar1+R_hat_1(i)**2
   s_bar2=s_bar2+R_hat_2(i)**2
!    write(*,*) ">>>>>>>>>>>>>", s_bar1,s_bar2
  enddo
!write(*,*) "R_hat1",R_hat_1
!write(*,*) "R_hat2",R_hat_2


 s_bar1=sqrt(s_bar1)
 s_bar2=sqrt(s_bar2)

  do k=1,3
   R_hat_1(k)=R_hat_1(k)/s_bar1
   R_hat_2(k)=R_hat_2(k)/s_bar2
  enddo

!write(*,*) "R_hat1",R_hat_1
!write(*,*) "R_hat2",R_hat_2


dot_pro=0.
do k=1,3
dot_pro=dot_pro+ R_hat_1(k)*R_hat_2(k)

enddo


theta=acos(dot_pro)

theta=theta*180/3.141593

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!------------------------

!write(*,*) "()()()()()()()()",theta
return
end subroutine tabbe_zaviyeyeab_cartz




subroutine  TORQUE(Nw,llwykf,platt2,pos,iwpos2,cSpcs2,achg2,occu2,pramt2,imodel,Hloc2,sym2)
use module_dictionary
use module_structure_symmetry_all

implicit double precision (a-h,o-z)

 integer,dimension(2),intent(inout):: Nw
 integer,dimension(Nw(1)),intent(inout):: llwykf
 real,dimension(3,3),intent(inout)::platt2
 real,dimension(Nw(1),200,3),intent(inout):: pos
 integer,dimension(Nw(1),3),intent(inout):: iwpos2
 Character(len=5),dimension(Nw(1),1):: cSpcs2
 real,dimension(Nw(1)),intent(inout):: achg2
 real,dimension(Nw(1),200),intent(inout):: occu2
 real,dimension(5),intent(inout)::pramt2
 integer,intent(inout):: imodel
 real,dimension(Nw(2),2,3),intent(inout):: Hloc2
 real,dimension(1000,5),intent(in):: sym2
! integer,intent(in):: N_Run

real,dimension(3):: r_l,rH_prim 
integer:: ierr
integer:: icont,ishm
real:: del_phi2 
real,allocatable,dimension(:,:):: torq_magnitud,delta_torq_mag,H_pos_extend, net_torque,torque_pstep
real,allocatable,dimension(:) :: torque_net_mag,torque_net_mag_pre,H_dis_extend,work,DelTht,trq_max
real,allocatable,dimension(:,:,:):: water_H_LatPos
real,allocatable,dimension(:,:,:):: hydrogen_coord
integer :: irun,nH=2 ,NN ,N_Run
real,dimension(3) :: shifted_H,r_ref,r_comper,H_ini,H_trans,H_pos_trans,rH1,rH2,rHc1,rHc2,O_lattice 
integer :: iwat_mar,mSteps
integer :: iTIP   
real:: start,finish,var_angle,dis1,dis2,disHh
real:: sum_ntrq,sum_ntrq_rate,DelRot
!integer:: T1, Ph1, X1
real :: hh1,hh2,Total_Err,atom
integer :: ix,iy,iz
integer:: ndim







!Hloc
!sym
CALL CPU_TIME(start)
 mSteps=12500
 NN=4
 !N_Run=1

 natypes=Nw(1)
 iwater=Nw(2)


iaPos_max=50  
iwykf_max=200 

ALLOCATE (iaPos(natypes),STAT= nst)
ALLOCATE (iwykf(natypes),STAT= nst)
ALLOCATE (occ(natypes,iwykf_max),STAT= nst)
ALLOCATE (r_wykf(natypes,iwykf_max,3),STAT= nst)
ALLOCATE (cAMC(natypes,1),STAT= nst)
ALLOCATE (n_oxidation(natypes),STAT= ierr)
ALLOCATE (lH_ip(natypes),STAT= ierr)
ALLOCATE (lO_elm(iwater),STAT= ierr)
ALLOCATE (lO_wykf(iwater),STAT= ierr)

call random_seed(size=ndim)             ! seeding   
allocate(iseed(ndim))
allocate(iseed_first(ndim))
allocate(iseed_last(ndim))



lH_ip=0
iaPos=1
	do i=1,3
	do j=1,3
	aL_matrix(i,j)=	platt2(i,j)
	enddo
	enddo

	do ip=1,natypes
	iwykf(ip)=llwykf(ip)
	cAMC(ip,1)=cSpcs2(ip,1)
	n_oxidation(ip)=achg2(ip)
	enddo

	do ip=1,natypes
	do iw=1,iwykf(ip)
		occ(ip,iw)=occu2(ip,iw)
		do i=1,3
			r_wykf(ip,iw,i)=pos(ip,iw,i)
		enddo
	enddo
	enddo


	iwat=0
	do ip=1,natypes
	  if (iwpos2(ip,2)/=0) then
		do jj=1,iwykf(iwpos2(ip,2))
		iwat=iwat+1
		lO_elm(iwat)=iwpos2(ip,2)
		lO_wykf(iwat)=jj
		enddo
	  endif
	  if (iwpos2(ip,3)/=0) then
		lH_ip(iwpos2(ip,3))=iwpos2(ip,3)
	  endif	
	enddo


!------------------------
	eps=pramt2(1)
	TTmax=pramt2(2)
	TTmin=pramt2(3)
	Seed=NINT(pramt2(4))
	N_Run=NINT(pramt2(5))
!	print *, eps,TTmax,TTmin,Seed



iseed= Seed



 	call arries_assigning(sym2)





	
imodel=model
iTIP=imodel
iTIP=0
print *, iTIP !ppt

if (iTIP==1) then 
h_oxidation=0.5564 !TIP4P/2005
else if (iTIP==2) then
h_oxidation=0.241 !TIP5P
else if (iTIP==0) then
h_oxidation=0.417 !TIP3P
endif


w_angle=104.52/2.
iwrite=2
d_O_H=0.9572
!--------------------------------------

if (iTIP==1) then 
delta_distance=0.1546 !TIP4P/2005
else if (iTIP==2) then
delta_distance=0.7 !TIP5P
else if (iTIP==0) then
delta_distance=0. !TIP3P
endif

delta_angle=109.47


if (iTIP==1 .or. iTIP==0) then 
delta_q=-2*h_oxidation
else if (iTIP==2) then
delta_q=-1*h_oxidation
endif

var_angle=w_angle
if (iTIP==0) iTIP=1 

atom=0
do ip=1,natypes
do j=1,iapos(ip)
atom = atom+ occ(ip,j)*iwykf(ip)
enddo
enddo



write (5,'(A17,3f12.5,10x,3f12.5)') "Lattice:  " , (acell(j),j=1,6)
write (5,'(A25,I5)') "Number of Atoms:  ", NINT(atom)
write (5,*) ""


ALLOCATE (torque_net_mag(iwater),STAT= ierr)
ALLOCATE (torque_net_mag_pre(iwater),STAT= ierr)
ALLOCATE (force_H(iwater,2,3),STAT= ierr)
ALLOCATE (torque_H(iwater,2,3),STAT= ierr)
ALLOCATE (hydrogen_coord(iwater,2,3),STAT= ierr)
ALLOCATE (delta_force(iwater,iTIP,3),STAT= ierr)
ALLOCATE (delta_torque(iwater,iTIP,3),STAT= ierr)
ALLOCATE (delta_charge(iTIP),STAT= ierr)
ALLOCATE (ox_shift_cz(iwater,iTIP,3),STAT= ierr)
ALLOCATE (ox_shift_latt(iwater,iTIP,3),STAT= ierr)
ALLOCATE (H_pos_extend(2*27*iwater,3),STAT= ierr)
ALLOCATE (H_dis_extend(2*27*iwater),STAT= ierr)
ALLOCATE (H_lattice_POSCAR(2*iwater,3),STAT= ierr)

ALLOCATE (torq_magnitud(iwater,2),STAT= ierr)
ALLOCATE (delta_torq_mag(iwater,iTIP),STAT= ierr)
ALLOCATE (net_torque(iwater,3),STAT= ierr)
ALLOCATE (torque_pstep(iwater,3),STAT= ierr)
ALLOCATE (work(iwater),STAT= ierr)
ALLOCATE (DelTht(iwater),STAT= ierr)
ALLOCATE (trq_max(iwater),STAT= ierr)
ALLOCATE (xyz_H_cartz(iwater,2,3),STAT= ierr)
ALLOCATE (xyz_H_lattice(iwater,2,3),STAT= ierr) 
ALLOCATE (water_H_LatPos(iwater,2,3),STAT= ierr) 
!--------------------------------------------------


call tabbe2_make_xyz_cartz_array() 



DO irun=1,N_Run  !"ijk" is a variable as a counter

print *, "Iteration:   " , irun, " out of ",N_Run
!T1=70
!Ph1=120
!X1=50

xyz_H_cartz=0.
xyz_H_lattice=0.

do i=1,iwater
torque_net_mag(i)=0.0  !*** shayad lazem nabashe
torque_net_mag_pre(i)=0.0
work(i)=0.0
enddo




call tabbe2_set_initial_H()

!-----------------------------------------
do iwat=1,iwater
 if (lO_wykf(iwat)==1) iwat_mar=iwat
do iH=1,2

  call tabbe_tabdil_cartezian_lattic(xyz_wykf_cartz(lO_elm(iwat_mar),1,:),r_comper)  
  call tabbe_tabdil_cartezian_lattic(xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),:),r_ref)
 	 H_ini(:)=xyz_H_lattice(iwat_mar,iH,:)

  call tabbe2_move_by_symmetry(r_ref,r_comper,H_ini,H_trans,1)
	xyz_H_lattice(iwat,iH,:)=H_trans
  call tabbe_tabdil_lattic_cartezian(H_trans,xyz_H_cartz(iwat,iH,:))
enddo
enddo
!------------------------------------------
call tabbe2_TIPxP(iTIP)  







DO iii=1,mSteps                  !<<<<<<<<<<<<<<<<<<main Do>>>>>>>>>>>>>>>>>>>>>>
  

if (mod(iii,2500)==0) TTmin=TTmin/2.

!****************           Force calculation 
do iwat=1,iwater
 if (lO_wykf(iwat)==1) then   

do iittpp=1,iTIP
call tabbe2_calculate_forece_on_TIPxP(ox_shift_cz(iwat,iittpp,:),iittpp,iTIP,iwat,NN) 
enddo




do iH=1,2
call tabbe2_calculate_forece_on_an_atom(iwat,iH,NN,iwater*2,iTIP,iii) 
					
enddo
endif

enddo
!*************************************         Calculating the Torques   
do iwat=1,iwater
 if (lO_wykf(iwat)==1) then 

call tabbe2_Torque_TIPxP(iwat,iH,iTIP)

do iH=1,2
call tabbe2_Torque(iwat,iH,iTIP)   
				
enddo

endif
enddo


net_torque=0.

!...........Magnitude of Torques & net_torque
do iwat=1,iwater
 if (lO_wykf(iwat)==1) then

do jtip=1,iTIP
net_torque(iwat,1)=net_torque(iwat,1)+ delta_torque(iwat,jtip,1) !torque on each water== net_torque(iwat,:)
net_torque(iwat,2)=net_torque(iwat,2)+ delta_torque(iwat,jtip,2)
net_torque(iwat,3)=net_torque(iwat,3)+ delta_torque(iwat,jtip,3)
delta_torq_mag(iwat,jtip)=sqrt(delta_torque(iwat,jtip,1)**2+delta_torque(iwat,jtip,2)**2+delta_torque(iwat,jtip,3)**2)
enddo


do iH=1,2
torq_magnitud(iwat,iH)=torque_H(iwat,iH,1)**2+torque_H(iwat,iH,2)**2+torque_H(iwat,iH,3)**2
torq_magnitud(iwat,iH)=sqrt(torq_magnitud(iwat,iH)) !torque on each H

net_torque(iwat,1)=net_torque(iwat,1)+ torque_H(iwat,iH,1) !torque on each water== net_torque(iwat,:)
net_torque(iwat,2)=net_torque(iwat,2)+ torque_H(iwat,iH,2)
net_torque(iwat,3)=net_torque(iwat,3)+ torque_H(iwat,iH,3)
enddo



torque_net_mag(iwat)= sqrt(net_torque(iwat,1)**2 + net_torque(iwat,2)**2 + net_torque(iwat,3)**2) !total net torque

endif
enddo

if (iii==1) torque_net_mag_pre(:)=torque_net_mag(:)

!***********************************		Comparing TORQUE ==> if T_net < 10^-5 exit the loop
mtrq=0
sum_ntrq=0.
do iwat=1,iwater
 if (lO_wykf(iwat)==1) then 

if (torque_net_mag(iwat) > eps ) mtrq=mtrq+1

sum_ntrq=sum_ntrq+ torque_net_mag(iwat)

endif
enddo

if (mtrq ==0 ) then 
print *, ""
print *, "The net torque reached to the threshold"
  !write (5,'(A20, I5,10x,A20,f12.6)') "Number of steps: ", iii,"Net Torque magnitude: ", sum_ntrq
	write (5,'(A80,I5)') "The net torque magnitude has riched to the threshold at step number: ",iii
exit
endif




!***********************************  angle of rotation
if (iii==1) then
trq_max=torque_net_mag(1)
endif
	do iwat=1,iwater
	 if (lO_wykf(iwat)==1) then
		DelTht(iwat)=( (TTmax-TTmin)* torque_net_mag(iwat)/trq_max(iwat) ) +TTmin
		if (trq_max(iwat) < torque_net_mag(iwat)) then
		trq_max(iwat)=torque_net_mag(iwat)
		endif
	endif
	enddo




!****************     Rotation   
do iwat=1,iwater
 if (lO_wykf(iwat)==1) then  

del_phi2= DelTht(iwat)
 call tabbe2_rotation(xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),:),xyz_H_cartz(iwat,1,:),rH_prim,net_torque(iwat,:),del_phi2)

	do k=1,3
	xyz_H_cartz(iwat,1,k)=rH_prim(k) 
	enddo
 call tabbe2_rotation(xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),:),xyz_H_cartz(iwat,2,:),rH_prim,net_torque(iwat,:),del_phi2)

	do k=1,3
	xyz_H_cartz(iwat,2,k)=rH_prim(k) 
	enddo

work(iwat)=work(iwat)+torque_net_mag(iwat)*(del_phi2*3.14159265/180.0)*14.400304


endif
enddo !iwater  ***********charkhandan  
!******************************************************************************

if (mod (iii,100)==0) then
!lw=1
	!write (*,'(3f10.5)') net_torque(lw,1),net_torque(lw,2),net_torque(lw,3)
	print *, "Steps: " , iii, "    Net_Torque: ", sum_ntrq
	write (5,'(A20, I5,10x,A20,f12.6)') "Number of steps: ", iii,"Net Torque magnitude: ", sum_ntrq
endif


!!*************update H_lattice**********************************************************
do iwat=1,iwater
 if (lO_wykf(iwat)==1) iwat_mar=iwat

call tabbe_tabdil_cartezian_lattic(xyz_H_cartz(iwat_mar,1,:),xyz_H_lattice(iwat_mar,1,:))
call tabbe_tabdil_cartezian_lattic(xyz_H_cartz(iwat_mar,2,:),xyz_H_lattice(iwat_mar,2,:))

do iH=1,2
  call tabbe_tabdil_cartezian_lattic(xyz_wykf_cartz(lO_elm(iwat_mar),1,:),r_comper)  
  call tabbe_tabdil_cartezian_lattic(xyz_wykf_cartz(lO_elm(iwat),lO_wykf(iwat),:),r_ref)
  H_ini(:)=xyz_H_lattice(iwat_mar,iH,:)

 call tabbe2_move_by_symmetry(r_ref,r_comper,H_ini,H_trans,1)
	xyz_H_lattice(iwat,iH,:)=H_trans

	call tabbe_tabdil_lattic_cartezian(H_trans,xyz_H_cartz(iwat,iH,:))
enddo
enddo
!***********************************************************

call tabbe2_TIPxP(iTIP) 

ENDDO                !<<<<<<<<<<<<<<<<<<<iii>>>>>>>>>>>>>>>>>>>>>>>

!*********************************8Cartezian coordination
print *, ""
print *, "Hloc has found:"

do iwat=1,iwater
do iH=1,2
	call tabbe_tabdil_cartezian_lattic(xyz_H_cartz(iwat,iH,:),xyz_H_lattice(iwat,iH,:))
	write(*,'(3f10.5)') xyz_H_lattice(iwat,iH,:)
enddo
enddo


	Total_Err=0.
do iwat=1,iwater
do iH=1,2
	icnt=0
	dis1= acell(1)+acell(2)+acell(3)
do ip =1, natypes
  if (lH_ip(ip) /= 0 ) then
	icnt=1
do iw=1,iwykf(ip)
	do ix=-1,1
	do iy=-1,1
	do iz=-1,1

	rH1(1)=r_wykf(ip,iw,1)+ix
	rH1(2)=r_wykf(ip,iw,2)+iy
	rH1(3)=r_wykf(ip,iw,3)+iz
	
	
	call tabbe_direct_distance(rH1, xyz_H_lattice(iwat,iH,:),dis2)
	
		if (dis2 < dis1) then
			dis1=dis2
			water_H_LatPos(iwat,iH,:)=rH1(:)
		endif
	enddo
	enddo
	enddo
    enddo !iw
   endif
  enddo !ip
	if (icnt==1)	Total_Err=Total_Err+dis1

 enddo  !iH
enddo   !iwat

Total_Err=Total_Err/real(2*iwater)

atom=0
do ip=1,natypes
do j=1,iapos(ip)
atom = atom+ occ(ip,j)*iwykf(ip)
enddo
enddo




	
k=0
do iwat=1,iwater
!do iH=1,2
   if (lO_wykf(iwat)==1) then
	k=k+1
	call tabbe_tabdil_cartezian_lattic(xyz_wykf_cartz(lO_elm(iwat),1,:),O_lattice)
	write (5,*) "Water number:  ", k
	write (5,*) "     Oxygen Position:" 
	write (5,'(A5,I5,A5,3f12.5)') "IN1: ",lO_elm(iwat)," O ", O_lattice
	write (5,*) "     H positions (Experiment) :"
	write (5,'(5x,A10,3f12.5)') "      H ", water_H_LatPos(iwat,1,:)
	write (5,'(5x,A10,3f12.5)') "      H ", water_H_LatPos(iwat,2,:)
	write (5,*) "     H positions (Torque Method) :"
	write (5,'(5x,A10,3f12.5)') "      H ", xyz_H_lattice(iwat,1,:)
	write (5,'(5x,A10,3f12.5)') "      H ", xyz_H_lattice(iwat,2,:)
	write (5,*) ""	
   endif
!enddo
enddo

write (5,*) ""
write (5,*) ""
write (5,'(A24,f10.3,A10)') "Mean Absolute Error: ", Total_Err, " Angstrom"




	icnt=0
do ip =1, natypes
  if (lH_ip(ip) /= 0 ) then
	icnt=1
	exit
	!write(*,'(A10,f12.6)') "MAE = ",Total_Err
  endif
enddo

	if (icnt==1) then
	write(*,'(A10,f12.6)') "MAE = ",Total_Err
	else
	write(*,'(A30)') "MAE =  N/A"
	endif
!------------------------------------

H_pos_trans=0.00

icont=0
ishm=0
do iwat=1,iwater
do iH=1,2

	call tabbe_tabdil_cartezian_lattic(xyz_H_cartz(iwat,iH,:),r_l)
	do i=-1,1
	do j=-1,1
	do k=-1,1
	icont=icont+1
	H_pos_extend(icont,1)=r_l(1)+i
	H_pos_extend(icont,2)=r_l(2)+j	
	H_pos_extend(icont,3)=r_l(3)+k
		
	call tabbe_direct_distance(H_pos_extend(icont,:),H_pos_trans,disHh) !distance from (0,0,0)
	H_dis_extend(icont)=disHh
	enddo
	enddo
	enddo

		ishm=ishm+1
		H_lattice_POSCAR(ishm,:)=r_l(:)
enddo
enddo
!H_lattice_POSCAR



CALL CPU_TIME(finish)

write(*,'(A10,f12.8)') "Time:", finish-start
write(5,'(A25,f10.2,A10)') "Time of calculation:  ", finish-start, "Seconds"


 write (5,*) "------------------------------------------------------"
 write (5,*) ""
 write (5,*) ""
 write (5,*) ""
 write (5,*) ""
 write (5,*) ""


ENDDO 

call zz_POSCAR()


DEALLOCATE (n_oxidation)
DEALLOCATE (xyz_wykf_cartz)
DEALLOCATE (force_H)
DEALLOCATE (torque_H)
DEALLOCATE (lO_elm)
DEALLOCATE (lO_wykf)
DEALLOCATE (torq_magnitud)
DEALLOCATE (hydrogen_coord)
DEALLOCATE (torque_net_mag)
DEALLOCATE (xyz_H_lattice)
DEALLOCATE (xyz_H_cartz)
DEALLOCATE (delta_force)
DEALLOCATE (delta_torque)
DEALLOCATE (ox_shift_cz)
DEALLOCATE (ox_shift_latt)
DEALLOCATE (delta_charge)
DEALLOCATE (H_pos_extend)
DEALLOCATE (H_dis_extend)
DEALLOCATE (H_lattice_POSCAR)
DEALLOCATE (net_torque)
DEALLOCATE (work)
DEALLOCATE (torque_pstep)
DEALLOCATE (DelTht)
DEALLOCATE (trq_max)
DEALLOCATE (lH_ip)
DEALLOCATE (water_H_LatPos)
deallocate(iseed)
deallocate(iseed_first)
deallocate(iseed_last)



return
end subroutine TORQUE





subroutine  zz_POSCAR()
use module_dictionary
use module_structure_symmetry_all
implicit double precision (a-h,o-z)




 integer:: ierror,istat,ierr,j,i,ip,iw



 open (unit=20, file='OUTPUT_Structure',status='replace', action='write' , IOstat= ierror)


write(20,*) ""
write(20,*) "   1.0000000000000000"

do i=1,3
write (20,'(3f21.16)') (aL_matrix(j,i),j=1,3 ) !lattice matrix
enddo


write(20,*)  "Direct"



do ip=1,natypes
do iw=1,iwykf(ip)
write(20,'(3(f10.5,"00000000000"))')  (r_wykf(ip,iw,j),j=1,3)
enddo
enddo

do i=1,iwater*2
!print*, ">>>>  ",H_lattice_POSCAR(i,:)
if ( H_lattice_POSCAR(i,1) < 0.00000000 ) H_lattice_POSCAR(i,1)=H_lattice_POSCAR(i,1)+1
if ( H_lattice_POSCAR(i,1) > 1.00000000 ) H_lattice_POSCAR(i,1)=H_lattice_POSCAR(i,1)-1	
if ( H_lattice_POSCAR(i,2) < 0.00000000 ) H_lattice_POSCAR(i,2)=H_lattice_POSCAR(i,2)+1
if ( H_lattice_POSCAR(i,2) > 1.00000000 ) H_lattice_POSCAR(i,2)=H_lattice_POSCAR(i,2)-1
if ( H_lattice_POSCAR(i,3) < 0.00000000 ) H_lattice_POSCAR(i,3)=H_lattice_POSCAR(i,3)+1
if ( H_lattice_POSCAR(i,3) > 1.00000000 ) H_lattice_POSCAR(i,3)=H_lattice_POSCAR(i,3)-1	

write(20,'(3(f10.5,"00000000000"))')  (H_lattice_POSCAR(i,j),j=1,3)
enddo


 close(20)
end subroutine zz_POSCAR

