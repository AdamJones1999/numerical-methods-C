# Author: Adam Jones

CC = gcc # c compiler used
CFLAGS = -g  # compiler flags (-g: enable debugging)
LIBS = # syntax: -l<lib_name> (omit <>)
LIB_DIRS = # syntax: -L<dir> (omit <>)
CPPFLAGS = -O2 -Wall -Wextra -pedantic

LFLAGS = $(LIBS) $(LIB_DIRS) # linker flags
# directories to search for header files, libraries, 
# and parts of compiler. All during preprocessing.
# makes nummethods.h, etc. considered system headers by compiler.
INCLUDE_DIRS = -I. 
SRC_DIR = .
BUILD_DIR = ./build
OBJ_DIR = $(BUILD_DIR)/obj
BIN_DIR = $(BUILD_DIR)/bin
LIB_DIR = $(BUILD_DIR)/lib

TARGET = $(BIN_DIR)/euler_test

all: make-dir $(TARGET)
make-dir:
	mkdir -p $(BIN_DIR) $(OBJ_DIR) $(LIB_DIR)

OBJS = $(OBJ_DIR)/nummethods.o $(OBJ_DIR)/euler_test.o

$(TARGET): $(OBJS)
	$(CC) -o $(TARGET) $(OBJS)

$(OBJ_DIR)/euler_test.o: $(SRC_DIR)/euler_test.c $(SRC_DIR)/nummethods.c $(SRC_DIR)/nummethods.h
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(SRC_DIR)/euler_test.c \
-o $(OBJ_DIR)/euler_test.o $(INCLUDE_DIRS)

$(OBJ_DIR)/nummethods.o: $(SRC_DIR)/nummethods.c $(SRC_DIR)/nummethods.h
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(SRC_DIR)/nummethods.c \
-o $(OBJ_DIR)/nummethods.o $(INCLUDE_DIRS)

.PHONY: clean 
clean:
	rm -r $(BUILD_DIR)
