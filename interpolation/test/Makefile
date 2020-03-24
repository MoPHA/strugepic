CXX=g++
CXX_FLAGS=-O3 -std=c++14 -Wall -Werror
LIB_NAME=interpol
LIB_DIR=../lib
INC_DIR=../include
LINK_FLAGS=-Wl,-rpath=$(LIB_DIR) -L$(LIB_DIR) -l$(LIB_NAME) -I$(INC_DIR)
TEST_SRC=test_interpolation.cpp
TEST_MAIN=test_main
CLEAN_C=rm -f *.out
MAXIMA_C=maxima  --very-quiet --batch-string
ERR_FILE=err.out

all: $(TEST_MAIN)

$(TEST_MAIN): $(TEST_SRC) 
	$(CXX) $(CXX_FLAGS) $(TEST_SRC) $(LINK_FLAGS) -o $(TEST_MAIN) 

run-test: $(TEST_MAIN) 
	$(CLEAN_C)
	$(MAXIMA_C) "load(\"generate_testdata.mc\")$$" >/dev/null
	@./$(TEST_MAIN) 2>$(ERR_FILE)  && $(CLEAN_C) || echo "Test failed, see $(ERR_FILE) for details"

gen-data: 
	$(MAXIMA_C) "load(\"generate_testdata.mc\")$$" >/dev/null

run-test-symbolic:
	@$(MAXIMA_C) "load(\"verify_symbolic.mc\")$$"


.PHONY: clean

clean:
	rm -f $(TEST_MAIN)
	$(CLEAN_C)
