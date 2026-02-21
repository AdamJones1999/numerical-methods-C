# Author: Adam Jones

CC = gcc # c compiler used
CFLAGS = -g # compiler flags (-g: enable debugging)
LIBS = -lm # syntax: -l<lib_name> (omit <>)
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
DATA_DIR = ./data

TARGET = $(BIN_DIR)/nummethods-test

all: make-dir $(TARGET)
make-dir:
	mkdir -p $(BIN_DIR) $(OBJ_DIR) $(LIB_DIR) $(DATA_DIR)

OBJS = $(OBJ_DIR)/nummethods.o $(OBJ_DIR)/nummethods-test.o $(OBJ_DIR)/odesystems.o

$(TARGET): $(OBJS)
	$(CC) -o $(TARGET) $(OBJS) $(LFLAGS)

$(OBJ_DIR)/nummethods-test.o: $(SRC_DIR)/nummethods-test.c \
$(SRC_DIR)/nummethods.c $(SRC_DIR)/nummethods.h \
$(SRC_DIR)/odesystems.c $(SRC_DIR)/odesystems.h
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(SRC_DIR)/nummethods-test.c \
-o $(OBJ_DIR)/nummethods-test.o $(INCLUDE_DIRS)

$(OBJ_DIR)/odesystems.o: $(SRC_DIR)/odesystems.c $(SRC_DIR)/odesystems.h
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(SRC_DIR)/odesystems.c \
-o $(OBJ_DIR)/odesystems.o $(INCLUDE_DIRS)

$(OBJ_DIR)/nummethods.o: $(SRC_DIR)/nummethods.c $(SRC_DIR)/nummethods.h
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(SRC_DIR)/nummethods.c \
-o $(OBJ_DIR)/nummethods.o $(INCLUDE_DIRS)

.PHONY: clean clean_data
clean:
	rm -r $(BUILD_DIR) $(DATA_DIR)
#clean_data:
#	rm -r $(DATA_DIR)
