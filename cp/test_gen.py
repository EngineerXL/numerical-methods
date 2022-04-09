from math import *
from random import *

TEST_FILE_COUNT = 10
MIN_A = 1
MAX_A = 10 ** 10000

for test in range(1, TEST_FILE_COUNT + 1):
    testFile = open("tests/" + str(test) + ".in", "w")
    ansFile = open("tests/" + str(test) + ".out", "w")

    a = randint(MIN_A, MAX_A)
    b = randint(MIN_A, MAX_A)

    # .in
    testFile.write(" ".join(str(elem) for elem in [a]) + "\n")
    testFile.write(" ".join(str(elem) for elem in [b]) + "\n")

    # .out
    ansFile.write(" ".join(str(elem) for elem in [a * b]) + "\n")

    testFile.close()
    ansFile.close()
