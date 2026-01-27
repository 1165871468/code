fc=gcc

all: timeDiff Error

timeDiff: Common.o timeDiff.o Rfile.o
	$(fc) -o $@ $^ -lm
Error: Common.o Error.o Rfile.o
	$(fc) -o $@ $^ -lm
Common.o: Common.c
	$(fc) -c -O3 $^ -lm
Error.o: Error.c
	$(fc) -c -O3 $^ -lm
timeDiff.o: timeDiff.c Common.c
	$(fc) -c -O3 $^ -lm
Rfile.o: Rfile.c
	$(fc) -c -O3 $^ -lm
.c.o:
	$(fc) -c -O3 $<
clean:
	rm *.o

