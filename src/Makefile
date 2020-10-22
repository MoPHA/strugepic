SRCS    := $(wildcard *.cpp)
OBJS    := $(patsubst %.cpp,%.o,$(SRCS))
CXXFLAGS+=-fPIC

all: $(OBJS)

%.o: %.cpp
	$(CXX) $(NVB)  $(CXXFLAGS) $(NVE) -I ../include $< -o $@ -c 

.PHONY: clean
clean:
	rm -f *.o
