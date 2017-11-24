all: use_matrices

#For debugging
OPT=-g -Wall
#For optimistaion
#OPT=-O

#All objects (except main) come from cpp and hpp
%.o:	%.cpp %.hpp
	g++ ${OPT} -c -o $@ $<
#use_vectors relies on objects which rely on headers
use_matrices:	use_matrices.cpp Matrix.o Vector.o Exception.o
		g++ ${OPT} -o use_matrices use_matrices.cpp Matrix.o Vector.o Exception.o

clean:
	rm -f *.o *~ use_matrices
