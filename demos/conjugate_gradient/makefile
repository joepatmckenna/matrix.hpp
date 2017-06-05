COMPILER = g++
OPTIONS = -std=c++11 -O3 -Wall
#PREPROCESS = sed -e 's/[0:1]/[{0:1}]/g' $(PREPROCESS) main.cpp | 

all: main

main: main.o
	$(COMPILER) $(OPTIONS) main.o -o main

main.o: main.cpp
	$(COMPILER) $(OPTIONS) -c main.cpp

clean:
	\rm *.o main
