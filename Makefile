options = -Wall -ffree-line-length-0 -cpp -ggdb -fmax-errors=1 -fcheck=all -O0 #-g

compiler = gfortran
linker = gfortran

OBJECTS = \
PhaseEnvelope.o \
EOS.o \
GaussElimination.o \
cubic_solver.o \

compile:	
	$(compiler) $(options) -c GaussElimination.f90 -o GaussElimination.o #compilar
	$(compiler) $(options) -c cubic_solver.f90 -o cubic_solver.o #compilar
	$(compiler) $(options) -c EOS.f90 -o EOS.o #compilar
	$(compiler) $(options) -c PhaseEnvelope.f90 -o PhaseEnvelope.o #compilar
	$(linker) $(options) $(OBJECTS) -o Trace_The_Envelope.elf #linkar
		

clean:
	rm -f *.o *.elf *.mod
