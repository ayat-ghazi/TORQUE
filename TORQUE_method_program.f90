program  TORQUE_method_program
use module_dictionary
use module_input
use module_structure_symmetry_all
implicit double precision (a-h,o-z)
integer:: ii,ierror
 character(len=128)::c_input,c_abs_adrs,PREFIX,c_OUT
 logical :: fexist




CALL system("clear")

IF (iargc()==0) then
	write (*,*) "Please insert the name of an input file as the argument."
	Stop
	
ELSE
  DO ii = 1, iargc()
    CALL getarg(ii, c_input)
    CALL GETCWD(PREFIX)
	c_abs_adrs= PREFIX(:LNBLNK(PREFIX)) // '/' // c_input(:LNBLNK(c_input))

	INQUIRE (file=c_abs_adrs , exist= fexist)
	IF (fexist .eqv. .false.) THEN
	PRINT *, "The input file does not exist."
	STOP
	ENDIF

    PRINT *, c_abs_adrs
  END DO
ENDIF


 open (unit=5, file='./OUTPUT_Results',status='replace', action='write' , IOstat= ierror)

 call read_input_file_2nd(c_input)

print *, Hloc
 call TORQUE(Nwatom,lwykf,platt,posatm,iwpos,cSpcs,achg,occu,pramt,model,Hloc,sym)




 DEALLOCATE(iaP)
 DEALLOCATE(lwykf)
 DEALLOCATE(occu)
 DEALLOCATE(posatm)
 DEALLOCATE(cSpcs)
 DEALLOCATE(achg)
 DEALLOCATE(ip_Ox)
 DEALLOCATE(iwyk_Ox)
 DEALLOCATE(ip_Hy)
 DEALLOCATE (Hloc)





 close (5)
END
