# Author: Adam Jones

CC = gcc # c compiler used
CFLAGS = -g  # compiler flags (-g: enable debugging)
LIBS = # syntax: -l<lib_name> (omit <>)
LIB_DIRS = # syntax: -L<dir> (omit <>)
CPPFLAGS = -O2 -Wall -Wextra -pedantic

LFLAGS = $(LIBS) $(LIB_DIRS) # linker flags
SRC_DIR = .

objects = euler_method.o euler_test.o
euler_method: $(objects)
	$(CC) -o euler_test $(objects)
euler_test.o: $(SRC_DIR)/euler_method.c $(SRC_DIR)/nummethods.h
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(SRC_DIR)/euler_test.c
euler_method.o: $(SRC_DIR)/euler_method.c $(SRC_DIR)/nummethods.h
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(SRC_DIR)/euler_method.c

.PHONY: clean 
clean:
	rm euler_test $(objects)
