SRCS    := $(wildcard *.cpp)
OBJS    := $(patsubst %.cpp,%.o,$(SRCS))

#CXXFLAGS+=-fPIC

all: $(DOBJS) $(OBJS)

%.o: %.cpp
	$(CXX) $(NVB)  $(CXXFLAGS) $(NVE) -I ../include $< -o $@  

.PHONY: clean
clean:
	rm -f *.o
