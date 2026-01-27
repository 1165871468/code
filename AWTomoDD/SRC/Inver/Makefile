fc = gcc
flag = -g

all: Inversion

Inversion: Inversion.o
	$(fc) -o $@ $^ -lm
Inversion.o: Inversion.c
	$(fc) -c -O3 $^ -lm
.c.o:
	$(fc) -c -O3 $<
clean:
	rm *.o
