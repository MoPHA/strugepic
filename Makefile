
CXX=g++
LIB_FLAGS=-fvisibility=hidden -fPIC -shared
CXX_FLAGS=-O3 -std=c++11 -Wall -Werror
LIB_NAME=interpol
LIB_SRC=internal_interpolation.hpp interpolation.cpp
LIB_HEADDER=interpolation.hpp
LINK_FLAGS=-Wl,-rpath=$(PWD) -L$(PWD) -l$(LIB_NAME)
TEST_SRC=test_interpolation.cpp
TEST_MAIN=test_main
CLEAN_C=rm -f *.out
MAXIMA_C=maxima  --very-quiet --batch-string
ERR_FILE=err.out

all: lib$(LIB_NAME).so $(TEST_MAIN)

lib$(LIB_NAME).so: $(LIB_SRC) 
	$(CXX) $(CXX_FLAGS) $(LIB_FLAGS) $(LIB_SRC) -o lib$(LIB_NAME).so


$(TEST_MAIN): lib$(LIB_NAME).so $(TEST_SRC) 
	$(CXX) $(CXX_FLAGS) $(TEST_SRC)  $(LINK_FLAGS) -o $(TEST_MAIN) 

run-test: $(TEST_MAIN) lib$(LIB_NAME).so
	$(CLEAN_C)
	$(MAXIMA_C) "load(\"generate_testdata.mc\")$$" >/dev/null
	@./$(TEST_MAIN) 2>$(ERR_FILE)  && $(CLEAN_C) || echo "Test failed, see $(ERR_FILE) for details"

run-test-symbolic:
	@$(MAXIMA_C) "load(\"verify_symbolic.mc\")$$"


.PHONY: clean

clean:
	rm -f $(TEST_MAIN) lib$(LIB_NAME).so
	$(CLEAN_C)
