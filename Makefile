# Author: Adam Jones

CC = gcc # c compiler used
CFLAGS = -g  # compiler flags (-g: enable debugging)
LIBS = # syntax: -l<lib_name> (omit <>)
LIB_DIRS = # syntax: -L<dir> (omit <>)
CPPFLAGS = -O2 -Wall -Wextra -pedantic

LFLAGS = $(LIBS) $(LIB_DIRS) # linker flags
SRC_DIR = .
BUILD_DIR = ./build
OBJ_DIR = $(BUILD_DIR)/obj
BIN_DIR = $(BUILD_DIR)/bin
LIB_DIR = $(BUILD_DIR)/lib

TARGET = euler_test

all: make-dir $(TARGET)
make-dir:
	mkdir -p $(BIN_DIR) $(OBJ_DIR) $(LIB_DIR)

OBJS = $(OBJ_DIR)/euler_method.o $(OBJ_DIR)/euler_test.o
$(TARGET): $(OBJS)
	$(CC) -o $(BIN_DIR)/euler_test $(OBJS)

$(OBJ_DIR)/euler_test.o: $(SRC_DIR)/euler_method.c $(SRC_DIR)/nummethods.h
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(SRC_DIR)/euler_test.c \
-o $(OBJ_DIR)/euler_test.o

$(OBJ_DIR)/euler_method.o: $(SRC_DIR)/euler_method.c $(SRC_DIR)/nummethods.h
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(SRC_DIR)/euler_method.c \
-o $(OBJ_DIR)/euler_method.o

.PHONY: clean 
clean:
	rm -r $(BUILD_DIR)
