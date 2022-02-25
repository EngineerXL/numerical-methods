FLAGS_G++ = -g -pedantic -std=c++17 -Wall -Wextra -Werror
# FLAGS_G++ = -g -O2 -pedantic -std=c++17 -Wall -Wextra

# Compile solution
compile: main.cpp
	g++ $(FLAGS_G++) main.cpp -o solution

# Compile solution for answer generation
correct: main.cpp
	g++ $(FLAGS_G++) main.cpp -o correct

# Compile with sanitizer
sanitize: main.cpp
	g++ $(FLAGS_G++) main.cpp -fsanitize=address -o solution

# Compile with gprof
gprof: main.cpp
	g++ $(FLAGS_G++) main.cpp -pg -o solution

# Clean tests and create new template file
new: clean
	bash scripts/generate_tests.sh
	cat template.cpp > main.cpp

# Check solution
check: compile
	bash scripts/checker.sh

# Check solution with valgrind
vcheck: compile
	bash scripts/vchecker.sh

# Run solution
run: compile
	./solution

# Run solution with valgrind
vrun: compile
	valgrind ./solution

# Generate tests
test: test_gen.py
	bash scripts/generate_tests.sh
	python3 test_gen.py

# Generate answers based on tests and curent solution
ans:
	bash scripts/gen_ans.sh

# Generate tests and run check
test_check: test check

# Generate tests, answers and run check
test_ans_check: test ans check

# Generate tests and run check with valgrind
test_vcheck: test vcheck

# Generate tests, answers and run check with valgrind
test_ans_vcheck: test ans vcheck

clean:
	rm -rf *.o solution correct
