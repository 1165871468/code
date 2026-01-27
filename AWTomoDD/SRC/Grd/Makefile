fc=mpicc
#fc=icc

all: Gradient

Gradient: SEM2D.o Readmesh.o Meshcount.o Matrix.o File.o Solve.o Rfile.o
	$(fc) -o $@ $^ -lm
Matrix.o: Matrix.c File.c
	$(fc) -c -O3 $^ -lm
Meshcount.o: Meshcount.c File.c
	$(fc) -c -O3 $^ -lm
Readmesh.o: Readmesh.c File.c
	$(fc) -c -O3 $^ -lm
Solve.o: Solve.c File.c
	$(fc) -c -O3 $^ -lm
.c.o:
	$(fc) -c -O3 $< -lm
clean:
	rm *.o

