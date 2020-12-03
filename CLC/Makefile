main:
	g++ -o lemon_gc main.cpp mpc/mpc.cpp  -I . -O2 -lemon

gc: simple.cpp Graph.hpp loader.hpp anchor.hpp
	g++ -o simple_gc simple.cpp -std=c++17 -O2 -lemon

gen: gen.cpp Graph.hpp loader.hpp anchor.hpp
	g++ -o gen gen.cpp -std=c++17 -O2 -lemon

brute: brute.cpp Graph.hpp loader.hpp anchor.hpp
	g++ -o brute brute.cpp -std=c++17 -O2 -lemon

all: main gc