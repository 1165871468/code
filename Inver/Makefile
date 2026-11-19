CC = gcc
# CFLAGS = -O3 -Wall
CFLAGS = -O3
LDLIBS = -lm

TARGETS = Inversion_tt_3d Inversion_cs Inversion_cr

all: $(TARGETS)

# 明确强制使用自定义规则，并将 LDLIBS 置于末尾
%: %.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o $(TARGETS)
