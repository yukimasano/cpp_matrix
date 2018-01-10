all: vary_N SRDD vary_kappa

#For debugging
#OPT=-g -Wall
#For optimistaion
OPT=-O

#All objects (except main) come from cpp and hpp
%.o:	%.cpp %.hpp
	g++ ${OPT} -c -std=c++11 -o $@ $<
#use_vectors relies on objects which rely on headers
vary_N:	vary_N.cpp Matrix.o Vector.o Exception.o
		g++ ${OPT} -o vary_N vary_N.cpp Matrix.o Vector.o Exception.o
SRDD:	SRDD.cpp Matrix.o Vector.o Exception.o
				g++ ${OPT} -o SRDD SRDD.cpp Matrix.o Vector.o Exception.o
vary_kappa:	vary_kappa.cpp Matrix.o Vector.o Exception.o
				g++ ${OPT} -o vary_kappa vary_kappa.cpp Matrix.o Vector.o Exception.o
clean:
	rm -f *.o *~ use_matrices
