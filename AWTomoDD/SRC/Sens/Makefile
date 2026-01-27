fc = mpicc
flag = -g

all: Sensitivity

Sensitivity: Common.o Interpolation.o Rfile.o Sensitivity.o
	$(fc) -o $@ $^ -lm
.c.o:
	$(fc) -c -O3 $<
clean:
	rm *.o
