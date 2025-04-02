all: RayTracer

clean:
	rm -f RayTracer

RayTracer: Vector.cpp
	g++ -fopenmp -o RayTracer Vector.cpp
