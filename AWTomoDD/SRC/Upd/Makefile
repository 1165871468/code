fc=mpicc

all: Update

Update: File.o Kernel.o Meshcount.o Matrix.o Readmesh.o Rfile.o SEM2D.o Solve.o SourceTime.o
	$(fc) -o $@ $^ -lm
Kernel.o: Kernel.c File.c
	$(fc) -c -O3 $^ -lm
Matrix.o: Matrix.c File.c
	$(fc) -c -O3 $^ -lm
Meshcount.o: Meshcount.c File.c
	$(fc) -c -O3 $^ -lm
Readmesh.o: Readmesh.c File.c
	$(fc) -c -O3 $^ -lm
Solve.o: Solve.c File.c
	$(fc) -c -O3 $^ -lm
SourceTime.o: Common.c SourceTime.c File.c
	$(fc) -c -O3 $^ -lm
.c.o:
	$(fc) -c -O3 $<
clean:
	rm *.o

