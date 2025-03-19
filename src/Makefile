# Compiler and flags
CC = mpicc
CFLAGS = -Wall -Wextra -O2 -lm -fopenmp

# Source files
SRCS = parallel-main.c dictionary.c util.c
HEADERS = dictionary.h util.h

# Object files
OBJS = $(SRCS:.c=.o)

# Executable name
TARGET = mitm
SEQ_TARGET = seq
SEQ_SRC = sequential.c

# Default target
all: $(TARGET) $(SEQ_TARGET)

# Rule to build the parallel executable
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS)

# Rule to build the sequential executable
$(SEQ_TARGET): $(SEQ_SRC)
	$(CC) -Wall -Wextra -O2 -lm -o $@ $(SEQ_SRC)

# Rule to build object files
%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
