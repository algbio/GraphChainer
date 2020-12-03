main:
	g++ -o lemon_gc main.cpp mpc/mpc.cpp  -I . -O2 -lemon

simple: simple.cpp Graph.hpp loader.hpp anchor.hpp
	g++ -o simple simple.cpp -std=c++17 -O2 -lemon

gen: gen.cpp Graph.hpp loader.hpp anchor.hpp
	g++ -o gen gen.cpp -std=c++17 -O2 -lemon

brute: brute.cpp Graph.hpp loader.hpp anchor.hpp
	g++ -o brute brute.cpp -std=c++17 -O2 -lemon


genref: genref.cpp Graph.hpp loader.hpp anchor.hpp
	g++ -o genref genref.cpp -std=c++17 -O2 -lemon

genread: genread.cpp Graph.hpp loader.hpp anchor.hpp
	g++ -o genread genread.cpp -std=c++17 -O2 -lemon

evaluate: evaluate.cpp Graph.hpp loader.hpp anchor.hpp
	g++ -o evaluate evaluate.cpp -std=c++17 -O2 -lemon

all: main gc