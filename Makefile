CXXPROC = clang++
CPPFLAGS = -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` `python3-config --includes --ldflags`# for debug
PROG_TARGET = flute_test.py
LIB_TARGET = ./src/repo/_flute.so
 

all: $(LIB_TARGET)

$(LIB_TARGET): ./src/repo/flute.cpp ./src/repo/flute.h
	$(CXXPROC) -o $@  $(CPPFLAGS) ${LDFLAGS} $<

.PHONY: test
test: _matrix.so
	env PYTHONPATH=".:$PYTHONPATH" python3 -m pytest -v

.PHONY: clean
clean:
	rm -rf *.o *.so .pytest_cache/ __pycache__/