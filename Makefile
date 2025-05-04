# компилятор который будет использоваться
CC=clang

# флаги компиляции
CFLAGS := ${CFLAGS} -I$(CURDIR)/dependencies/libomp/20.1.1/include -Xpreprocessor -fopenmp

# флаги линковки
LDFLAGS := ${LDFLAGS} -L$(CURDIR)/dependencies/libomp/20.1.1/lib -lomp -L$(CURDIR)/dependencies/glfw/3.4/lib -lglfw -framework OpenGL -framework GLUT -L$(CURDIR)/dependencies/glew/2.2.01/lib -lglew

# Исходные файлы
SRCS = main.c fem.c formation_mtrx.c

# Исполняемый файл
TARGET = a.out

# путь к glew библиотеке
DYLD_LIBRARY_PATH := $(CURDIR)/dependencies/glew/2.2.01/lib:
# export DYLD_LIBRARY_PATH="$(DYLD_LIBRARY_PATH)";

all: compilation
	
compilation: $(SRCS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(SRCS) -o $(TARGET)

clean:
	rm -rf $(TARGET)