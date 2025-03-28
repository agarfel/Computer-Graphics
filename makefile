all: RayTracer

clean:
	rm -f Vector.o RayTracer.o RayTracer

RayTracer: RayTracer.o Vector.o
	g++ -g -O3 -Wall -o RayTracer RayTracer.o -lm

RayTracer.o: Vector.cpp
	g++ -c -O3 -Wall -o RayTracer.o Vector.cpp
