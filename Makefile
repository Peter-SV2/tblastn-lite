# The command-line tool. The GUI needs ImGui and is built with CMake:
#   cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build
CXX      ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall
LDFLAGS  ?= -pthread
# a truly standalone binary on mingw/windows: make LDFLAGS=-static
BIN      := tblastn_lite

$(BIN): tblastn_lite.cpp tblastn_core.h
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

test: $(BIN)
	./$(BIN) --selftest
	./$(BIN) -q test/ecoli_5.tsv -a thrA -d test/ecoli_50kb.fna -e 1e-10

clean:
	rm -f $(BIN) $(BIN).exe

.PHONY: test clean
