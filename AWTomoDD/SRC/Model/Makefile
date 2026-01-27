fc=gcc

all: Initial modelUpdate

Initial: Common.o Initial.o Rfile.o Interpolation.o
	$(fc) -o $@ $^ -lm
modelUpdate: Common.o modelUpdate.o Rfile.o
	$(fc) -o $@ $^ -lm

Common.o: Common.c
	$(fc) -c -O3 $^ -lm
Initial.o: Initial.c
	$(fc) -c -O3 $^ -lm
Interpolation.o: Interpolation.c
	$(fc) -c -O3 $^ -lm
modelUpdate.o: modelUpdate.c
	$(fc) -c -O3 $^ -lm
Rfile.o: Rfile.c
	$(fc) -c -O3 $^ -lm
.c.o:
	$(fc) -c -O3 $<
clean:
	rm *.o

