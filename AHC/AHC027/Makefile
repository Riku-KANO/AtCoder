SRC_DIR := ./src
BIN_DIR := ./bin
SRC_NAME := $(SRC_DIR)/solve.cpp
EXE_NAME := $(BIN_DIR)/solve
CXX := g++
CXXFLAGS := -O2 -DLOCAL -DDEBUG -std=c++17

.PHONY: build clean

build: $(EXE_NAME)

$(EXE_NAME): $(SRC_NAME)
	$(CXX) $< -o $@ $(CXXFLAGS)

clean:
	rm -f $(EXE_NAME)