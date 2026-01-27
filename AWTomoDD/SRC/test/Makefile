CC=mpicc
CFLAGS=-Wall -I/usr/local/include
LDFLAGS=-L/usr/local/lib -lgsl -lgslcblas -lm

all: gsl_example

gsl_example: cc1.o
	$(CC) -o $@ $^ -lgsl -lgslcblas -lm
cc1.o: cc1.c
	$(CC) -c -O3 $^ -lgsl -lgslcblas -lm
clean:
	rm -f *.o gsl_example
