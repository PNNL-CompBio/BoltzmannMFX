FC=gfortran
FCFLAGS=-c -fbacktrace -g -O0

ascii-tools: \
	 ascii-reader.mod \
	 ascii-tools.o
	${FC} -o ascii-tools ascii-tools.o ascii-reader.o

# Link command:

ascii-reader.mod: ascii-reader.f08
	${FC} ${FCFLAGS} ascii-reader.f08

ascii-tools.o: ascii-tools.f08 ascii_reader.mod
	${FC} ${FCFLAGS} ascii-tools.f08 -o ascii-tools.o

.PHONY : clean
clean :
	rm -f ascii-tools *.o *.mod
