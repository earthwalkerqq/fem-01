# fem-01

clang draw_triangl.c -lglew -lglfw -framework OpenGL -framework GLUT

clang -I/opt/homebrew/Cellar/libomp/20.1.1/include -L/opt/homebrew/Cellar/libomp/20.1.1/lib -Xpreprocessor -fopenmp -lomp main.c fem.c formation_mtrx.c draw.c -framework OpenGL -framework GLUT

export DYLD_LIBRARY_PATH="/opt/homebrew/Cellar/glew/2.2.01/lib:$DYLD_LIBRARY_PATH"

clang -I/opt/homebrew/Cellar/libomp/20.1.1/include -L/opt/homebrew/Cellar/libomp/20.1.1/lib -Xpreprocessor -fopenmp -lomp main.c fem.c formation_mtrx.c -L/opt/homebrew/Cellar/glfw/3.4/lib -lglfw -framework OpenGL -framework GLUT -L/opt/homebrew/Cellar/glew/2.2.01/lib -lglew