CXX      ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall
LDFLAGS  ?= -pthread
# a truly standalone binary on mingw/windows: make LDFLAGS=-static
BIN      := tblastn_lite

$(BIN): tblastn_lite.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

test: $(BIN)
	./$(BIN) --selftest
	./$(BIN) -q test/thrA.faa -d test/ecoli_50kb.fna -e 1e-10

clean:
	rm -f $(BIN) $(BIN).exe

.PHONY: test clean
