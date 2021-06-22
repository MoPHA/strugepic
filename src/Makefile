SRCS    := $(wildcard *.cpp)
OBJS    := $(patsubst %.cpp,%.o,$(SRCS))


all:  $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS)  $(AMREX_CFLAGS) -I ../include $< -o $@  

.PHONY: clean
clean:
	rm -f *.o
