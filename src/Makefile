SRCS    := $(wildcard *.cpp)
OBJS    := $(patsubst %.cpp,%.o,$(SRCS))

all: $(OBJS)

%.o: %.cpp
	$(MPICXX) $(CXXFLAGS) -fPIC -I ../include $< -o $@ -c 



.PHONY: clean
clean:
	rm -f *.o
