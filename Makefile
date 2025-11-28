# Author: Adam Jones

CC = gcc # c compiler used
CFLAGS = -g  # compiler flags
LIBS = # syntax: -l<lib_name> (omit <>)
LIB_DIRS = # syntax: -L<dir> (omit <>)
CPPFLAGS = -Wall -Wextra -pedantic

LFLAGS = $(#linker flags)

objects = c-review.o
c-review: $(objects)
	$(CC) -o c-review $(objects) $(CFLAGS)
c-review.o: c-review.c
	$(CC) -O2 -c c-review.c

.PHONY: clean 
clean:
	rm c-review $(objects)