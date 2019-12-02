

all: libs tests 

libs:
	cd src && $(MAKE) all

tests:
	cd test && $(MAKE) all
run-test:
	cd test && $(MAKE) run-test

.PHONY: clean

clean:
	cd test && $(MAKE) clean
	cd src && $(MAKE) clean
