SRCS    := $(wildcard *.cpp)
HDS     := $(wildcard *.hpp)
OBJS    := $(patsubst %.cpp,%.o,$(SRCS))


ifeq ($(INTERPOLATION_FUNC),P8R2)
CXXFLAGS+=-DINTERPOLATION_P8R2=1
endif
ifeq ($(INTERPOLATION_FUNC),PWL)
CXXFLAGS+=-DINTERPOLATION_PWL=1
endif

ifeq ($(BUILD),GPU)
CXXFLAGS+=-DGPUBUILD
endif

all: $(DOBJS) $(OBJS)

%.o:  %.cpp $(HDS)
	$(CXX) $(NVB) $(CXXFLAGS) $(NVE)  $< -o $@ 

.PHONY: clean
clean:
	rm -f *.o
