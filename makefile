# Generates the geometry program!
 FC = ifort
 FCFLAGS =
# FC = g95
#FC = gfortran
#FFLAGS = -Wall -fcheck=all -pedantic-errors -ffixed-line-length-132

OBJECTS = TORQUE_method_program.o sub_tools.o module_dictionary.o module_structure_symmetry_all.o module_input.o




TQ:  $(OBJECTS)
	$(FC) -132 -o $@ $^ 

TORQUE_method_program.o: TORQUE_method_program.f90 module_dictionary.mod module_structure_symmetry_all.mod module_input.mod
	$(FC) $(FCFLAGS) -c $<

sub_tools.o: sub_tools.f90 module_dictionary.mod module_structure_symmetry_all.mod module_input.mod
	$(FC) $(FCFLAGS) -c $<


module_dictionary.o: module_dictionary.f90
	$(FC) $(FCFLAGS) -c $<


module_structure_symmetry_all.o: module_structure_symmetry_all.f90
	$(FC) $(FCFLAGS) -c $<


module_input.o: module_input.f90
	$(FC) $(FCFLAGS) -c $<

 
%.o %.mod: %.f90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f $(EXEC) *.mod *.o 
