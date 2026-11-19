fc=mpicc
CFLAGS = -O3
LDFLAGS = -lm

all: Update

Update: File.o Kernel.o Meshcount.o Matrix.o Readmesh.o Rfile.o SEM2D.o Solve.o SourceTime.o Python.o
	$(fc) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(fc) -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o Update