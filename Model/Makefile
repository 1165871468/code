
fc=gcc
CFLAGS = -O3
LDFLAGS = -lm

all: Initial

Initial: Common.o Initial.o Rfile.o Interpolation.o
	$(fc) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(fc) -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o Initial
