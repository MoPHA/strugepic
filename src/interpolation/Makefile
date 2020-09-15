SRCS    := $(wildcard *.cpp)
OBJS    := $(patsubst %.cpp,%.o,$(SRCS))

all: $(OBJS)

%.o : %.cpp %.hpp
	$(CXX) $(CXXFLAGS) $< -o $@ -c

.PHONY: clean
clean:
	rm -f *.o
