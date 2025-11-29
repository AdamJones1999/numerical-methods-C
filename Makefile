# Author: Adam Jones

CC = gcc # c compiler used
CFLAGS = -g  # compiler flags (-g: enable debugging)
LIBS = # syntax: -l<lib_name> (omit <>)
LIB_DIRS = # syntax: -L<dir> (omit <>)
CPPFLAGS = -O2 -Wall -Wextra -pedantic

LFLAGS = $(LIBS) $(LIB_DIRS) # linker flags

objects = rk4.o
rk4: $(objects)
	$(CC) -o rk4 $(objects)
rk4.o: rk4.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) rk4.c

.PHONY: clean 
clean:
	rm rk4 $(objects)
