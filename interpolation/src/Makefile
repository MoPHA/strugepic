CXX=g++
LIB_FLAGS=-fPIC -shared
CXX_FLAGS=-O3 -std=c++14 -Wall -Werror
LIB_NAME=interpol
LIB_SRC=./poly8_range2.cpp
LIB_DIR=../lib

all: $(LIB_DIR)/lib$(LIB_NAME).so

$(LIB_DIR)/lib$(LIB_NAME).so: $(LIB_SRC) 
	$(CXX) $(CXX_FLAGS) $(LIB_FLAGS) $(LIB_SRC) -o $(LIB_DIR)/lib$(LIB_NAME).so


.PHONY: clean

clean:
	rm -f $(LIB_DIR)/lib$(LIB_NAME).so
