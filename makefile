# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2
# FINAL DEBUGGING STEP: Removed -mwindows to make the console visible
LDFLAGS = -lSDL3 -lOpenGL32 -lShell32 -lgdi32

# Project structure
ROOT_DIR = .
BIN_DIR = $(ROOT_DIR)/bin
# FINAL FIX: Pointing to the 64-bit (x86_64) library to match your compiler
SDL_DIR = $(BIN_DIR)/SDL3-devel-3.2.22-mingw/SDL3-3.2.22/x86_64-w64-mingw32

# Paths
INCLUDE_PATHS = -I"$(SDL_DIR)/include"
LIB_PATHS = -L"$(SDL_DIR)/lib"

# Source files (only main.cpp for this test)
SRCS = main.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Executable name
TARGET = minimal_test

.PHONY: all clean run

# Default target
all: $(TARGET)

# Link the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIB_PATHS) $(LDFLAGS)

# Compile source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATHS) -c $< -o $@

# Clean up build files
clean:
	rm -f $(OBJS) $(TARGET)

# Command to run the demo
run: all
	./$(TARGET)

