FLAGS_G++ = -g -O2 -pedantic -std=c++17 -Wall -Wextra -Werror

OUTPUT_NAME = solver

# Compile solution
solver: main.cpp hpde.hpp
	g++ $(FLAGS_G++) main.cpp -o $(OUTPUT_NAME)

clean:
	rm -rf *.o $(OUTPUT_NAME)
