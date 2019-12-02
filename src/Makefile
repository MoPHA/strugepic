CXX=g++
LIB_FLAGS=-fvisibility=hidden -fPIC -shared
CXX_FLAGS=-O3 -std=c++11 -Wall -Werror
LIB_NAME=interpol
LIB_SRC=internal_interpolation.hpp interpolation.cpp
LIB_HEADDER=interpolation.hpp
LIB_DIR=../lib

all: $(LIB_DIR)/lib$(LIB_NAME).so

$(LIB_DIR)/lib$(LIB_NAME).so: $(LIB_SRC) 
	$(CXX) $(CXX_FLAGS) $(LIB_FLAGS) $(LIB_SRC) -o $(LIB_DIR)/lib$(LIB_NAME).so


.PHONY: clean

clean:
	rm -f $(LIB_DIR)/lib$(LIB_NAME).so
