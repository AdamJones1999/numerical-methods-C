# Author: Adam Jones

CC = gcc # c compiler used
CFLAGS = -g  # compiler flags
LIBS = # syntax: -l<lib_name> (omit <>)
LIB_DIRS = # syntax: -L<dir> (omit <>)
CPPFLAGS = -Wall -Wextra -pedantic

LFLAGS = $(#linker flags)

objects = rk4.o
c-review: $(objects)
	$(CC) -o rk4.c $(objects) $(CFLAGS)
c-review.o: c-review.c
	$(CC) -O2 -c rk4.c

.PHONY: clean 
clean:
	rm rk4 $(objects)
