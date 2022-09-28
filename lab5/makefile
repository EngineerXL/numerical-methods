FLAGS_G++ = -g -pedantic -std=c++17 -Wall -Wextra -Werror

OUTPUT_NAME = solver

# Compile solution
compile: main.cpp ppde.hpp
	g++ $(FLAGS_G++) main.cpp -o $(OUTPUT_NAME)

# Run solution
run: compile
	./$(OUTPUT_NAME)

# Run solution with valgrind
vrun: compile
	valgrind ./$(OUTPUT_NAME)

clean:
	rm -rf *.o $(OUTPUT_NAME)
