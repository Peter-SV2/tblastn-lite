CXX      ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall
LDFLAGS  ?= -pthread
# a truly standalone binary on mingw/windows: make LDFLAGS=-static
BIN      := tblastn_lite
GUI      := tblastn_gui

$(BIN): tblastn_lite.cpp tblastn_core.h
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

# Windows only (Win32 API, no toolkit)
gui: tblastn_gui.cpp tblastn_core.h
	$(CXX) $(CXXFLAGS) -mwindows -o $(GUI) $< $(LDFLAGS) -lcomdlg32 -lshell32

test: $(BIN)
	./$(BIN) --selftest
	./$(BIN) -q test/thrA.faa -d test/ecoli_50kb.fna -e 1e-10

clean:
	rm -f $(BIN) $(BIN).exe $(GUI) $(GUI).exe

.PHONY: gui test clean
