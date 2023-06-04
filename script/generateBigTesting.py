import sys
import random

with open('../testFile/30Test.txt', 'w')as fp:
    for i in range(30):
        x = random.randint(0, 20)
        y = random.randint(0, 20)
        fp.write(str(x) + " " + str(y) + "\n")
