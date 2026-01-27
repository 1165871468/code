fc=gcc

all: GetKernel

GetKernel: Common.o GetKernel.o Rfile.o
	$(fc) -o $@ $^ -lm
Common.o: Common.c
	$(fc) -c -O3 $^ -lm
GetKernel.o: GetKernel.c
	$(fc) -c -O3 $^ -lm
Rfile.o: Rfile.c
	$(fc) -c -O3 $^ -lm
.c.o:
	$(fc) -c -O3 $<
clean:
	rm *.o

