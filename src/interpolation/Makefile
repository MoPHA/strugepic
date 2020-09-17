SRCS    := $(wildcard *.cpp)
HDS     := $(wildcard *.hpp)
OBJS    := $(patsubst %.cpp,%.o,$(SRCS))


ifeq ($(INTERPOLATION_FUNC),P8R2)
CXXFLAGS+=-DINTERPOLATION_P8R2=1
endif
ifeq ($(INTERPOLATION_FUNC),PWL)
CXXFLAGS+=-DINTERPOLATION_PWL=1
endif

all: $(OBJS)

%.o:  %.cpp $(HDS)
	$(CXX) -fPIC $(CXXFLAGS)  $< -o $@ -c

.PHONY: clean
clean:
	rm -f *.o
