# flags
CC=clang

OS=$(shell uname -s)

# Default paths
LIBOMP_PATH := /usr/lib/libomp
GLFW_PATH := /usr/lib/libglfw
GLEW_PATH := /usr/lib/libglew

CFLAGS=-Wall -Werror -Wextra -std=c99 \
	-I$(LIBOMP_PATH)/include -Xpreprocessor

ifeq ($(OS), Darwin)
	OPENMP_FLAG=-Xpreprocessor -fopenmp
else
	OPENMP_FLAG=-fopenmp
endif

MAC_FLAGS=-lcheck -lm
LINUX_FLAGS=-lcheck -lm -lrt -lpthread -lsubunit
GCOV_FLAGS=-fprofile-arcs -ftest-coverage

LDFLAGS = -L$(LIBOMP_PATH)/lib -lomp \
	-L$(GLFW_PATH)/lib -lglfw \
	-framework OpenGL \
	-framework GLUT \
	-L$(GLEW_PATH)/lib -lglew \
	$(OPENMP_FLAG)

OBJ=$(patsubst %.c, %.o, $(FILES))

FILES=$(wildcard *.c)
TEST_FILES=$(wildcard tests/test_*.c)

MY_LIB=fem-01.a

all: check_deps $(MY_LIB)

.PHONY: all test memory_test $(MY_LIB) gcov_report clean rebuild check_deps

check_deps:
ifeq ($(OS), Darwin)
	@echo "Checking dependencies for MacOS..."
	@brew update
	@if ! brew list libomp; then \
		echo "Installing libomp..."; \
		brew install libomp; \
	fi
	@if ! brew list glfw; then \
		echo "Installing glfw..."; \
		brew install glfw; \
	fi
	@if ! brew list glew; then \
		echo "Installing glew..."; \
		brew install glew; \
	fi
	$(eval LIBOMP_PATH := $(shell brew --prefix libomp))
	$(eval GLFW_PATH := $(shell brew --prefix glfw))
	$(eval GLEW_PATH := $(shell brew --prefix glew))
else
	@echo "Checking dependencies for Linux..."
	@if ! dpkg -l | grep -q libomp-dev; then \
		echo "Installing libomp-dev..."; \
		sudo apt update && sudo apt install -y libomp-dev; \
	fi
	@if ! dpkg -l | grep -q libglfw3-dev; then \
		echo "Installing libglfw3-dev..."; \
		sudo apt update && sudo apt install -y libglfw3-dev; \
	fi
	@if ! dpkg -l | grep -q libglew-dev; then \
		echo "Installing libglew-dev..."; \
		sudo apt update && sudo apt install -y libglew-dev; \
	fi
endif

$(MY_LIB): $(FILES)
	$(CC) -c $(CFLAGS) $(FILES)
	ar rc $@ $(OBJ)
	ranlib $@
	rm -rf *.0

test: $(MY_LIB)
ifeq ($(OS), Darwin)
	$(CC) $(CFLAGS) $(GCOV_FLAGS) $(TEST_FILES) $< -o test $(MAC_FLAGS) $(LDFLAGS)
else
	$(CC) $(CFLAGS) $(GCOV_FLAGS) $(TEST_FILES) $< -o test $(LINUX_FLAGS) $(LDFLAGS)
endif
	./test

gcov_report: clean
ifeq ($(OS), Darwin)
	$(CC) $(CFLAGS) $(GCOV_FLAGS) $(TEST_FILES) -o test $(MAC_FLAGS) $(LDFLAGS)
else
	$(CC) $(CFLAGS) $(GCOV_FLAGS) $(TEST_FILES) -o test $(LINUX_FLAGS) $(LDFLAGS)
endif
	./test
	rm test_*
	lcov -t "fem-01_test" -o postreport.info -c -d .
	genheml -o report postreport.info
	$(OPEN) report/index.html

style_check:
	clang-format -i -style=google tests/*.h tests/*.c *.c *.h
	clang-format -n -style=google tests/*.h tests/*.c *.c *.h

clean_obj:
	rm -rf obj/*.o *.o

clean_lib:
	rm -rf *.a

clean_test:
	rm -rf *.gcda
	rm -rf *.gcno
	rm -rf *.info
	rm -rf test
	rm -rf report
	rm -rf valgrind.txt

clean: clean_lib clean_obj clean_test
	rm -rf unot_test
	rm -rf RESULT_VALGRIND.txt

rebuild: clean all