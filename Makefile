# Compiler
CC = clang

# Library paths
LIBOMP_PATH = /opt/homebrew/Cellar/libomp/20.1.1
GLFW_PATH = /opt/homebrew/Cellar/glfw/3.4
GLEW_PATH = /opt/homebrew/Cellar/glew/2.2.01

# Compilation flags
CFLAGS = -I$(LIBOMP_PATH)/include -Xpreprocessor -fopenmp

# Linker flags
LDFLAGS = -L$(LIBOMP_PATH)/lib -lomp \
          -L$(GLFW_PATH)/lib -lglfw \
          -framework OpenGL \
          -framework GLUT \
          -L$(GLEW_PATH)/lib -lglew

# Source files
SRCS = main.c fem.c formation_mtrx.c

TARGET = a.out

# Set DYLD_LIBRARY_PATH for GLEW
export DYLD_LIBRARY_PATH := $(GLEW_PATH)/lib:$(DYLD_LIBRARY_PATH)

.PHONY: all clean

all:
	$(CC) $(CFLAGS) $(SRCS) $(LDFLAGS) -o $(TARGET)

clean:
	rm -rf $(TARGET)