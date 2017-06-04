COMPILER = g++
OPTIONS = -std=c++11 -O3 -Wall
#PREPROCESS = sed -e 's/[0:1]/[{0:1}]/g' $(PREPROCESS) main.cpp | 

all: conj_grad

conj_grad: conj_grad.o
	$(COMPILER) $(OPTIONS) conj_grad.o -o conj_grad

conj_grad.o: conj_grad.cpp
	$(COMPILER) $(OPTIONS) -c conj_grad.cpp

clean:
	\rm *.o conj_grad
