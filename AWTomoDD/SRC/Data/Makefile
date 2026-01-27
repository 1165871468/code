fc=gcc

all: Get_trval

Get_trval: Common.o Get_trval.o Rfile.o
	$(fc) -o $@ $^ -lm
Common.o: Common.c
	$(fc) -c -O3 $^ -lm
Get_trval.o: Get_trval.c
	$(fc) -c -O3 $^ -lm
Rfile.o: Rfile.c
	$(fc) -c -O3 $^ -lm
.c.o:
	$(fc) -c -O3 $<
clean:
	rm *.o

