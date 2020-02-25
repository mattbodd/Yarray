CC = gcc
CFLAGS = -Wall -g
DEPS = yarray.h yutilities.h
OBJ = yoperations.o yutilities.o yarray.o

# store binaries
#binaries = y_ops

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

y_ops: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

# Use a phony pseudo target to distinguish `clean` target from being a dependency for anyone else
#.PHONY: clean
#clean:
#	rm -f $(binaries) *.o
