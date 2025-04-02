#include <vector>

class Vector {
public:
	explicit Vector(double, double, double);
	double norm2() const;
	double norm() const;
	void normalize();
	double operator[](int) const;
	double& operator[](int);
	double data[3];
};

Vector operator+(const Vector&, const Vector&);
Vector operator-(const Vector&, const Vector&);
Vector operator*(const double, const Vector&);
Vector operator*(const Vector&, const double);
Vector operator/(const Vector&, const double);
double dot(const Vector&, const Vector&);
Vector cross(const Vector&, const Vector&);

class Ray{
    public:
        explicit Ray(Vector, Vector);
        Vector O;
        Vector u;
    };

class Light{
	public:
		explicit Light(Vector, double);
		Vector P; // Poisition
		double I; // Intensity
	};

class Sphere {
public:
	explicit Sphere(Vector, double, Vector, bool mirror = false);
    Vector C;
    double R;
    Vector albedo;
	bool Mirror;
	struct Intersection intersect(Ray&) const;
};


struct Intersection{
	bool intersects;
	Vector point;
};


class Scene {
	public:
		explicit Scene();
        std::vector<Sphere> arr;
        struct Intersection intersect(Ray&) const;
        Vector get_colour(Ray&, int) const;
		std::vector<Light> lights;
};

Vector random_vector(Vector&);

int main();