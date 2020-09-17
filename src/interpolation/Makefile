SRCS    := $(wildcard *.cpp)
HDS     := $(wildcard *.hpp)
OBJS    := $(patsubst %.cpp,%.o,$(SRCS))



all: $(OBJS)

%.o:  %.cpp $(HDS)
	$(CXX) -fPIC $(CXXFLAGS) $(INTERPOLATION_NAME)  $< -o $@ -c

.PHONY: clean
clean:
	rm -f *.o
