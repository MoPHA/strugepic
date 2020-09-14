DEPS=src/propagators.o src/amrex_util.o 
CXXFLAGS+=-fPIC -I ./include 

lib/lib$(LNAME).so: $(DEPS) 
	mkdir -p lib
	$(MPICXX) -shared  $(DEPS) $(CXXFLAGS) -o lib/lib$(LNAME).so

src/%.o : src/%.cpp
	cd src
	$(MPICXX) $(CXXFLAGS) $<  -o $@ -c


.PHONY: clean
clean:
	rm -rf lib
	cd src && rm -f *.o
.PHONY: tests
tests:
	cd tests && $(MAKE)

.PHONY: install
install:
	mkdir $(PREFIX)/lib
	mkdir $(PREFIX)/include
	cp include/*.hpp $(PREFIX)/include
	cp lib/lib$(LNAME).so $(PREFIX)/lib


