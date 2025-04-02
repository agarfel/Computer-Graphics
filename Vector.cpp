#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include "Vector.h"
#include <iostream>
#include <random>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

std::default_random_engine engine;
std::uniform_real_distribution<double> unif(0.0, 1.0);

Vector::Vector(double x = 0, double y = 0, double z = 0) {
	data[0] = x;
	data[1] = y;
	data[2] = z;
}
double Vector::norm2() const {
	return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
}
double Vector::norm() const {
	return sqrt(norm2());
}
void Vector::normalize() {
	double n = norm();
	data[0] /= n;
	data[1] /= n;
	data[2] /= n;
}
double Vector::operator[](int i) const { return data[i]; };
double& Vector::operator[](int i) { return data[i]; };

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}

Vector operator*(const Vector& a, const Vector& b) {
	return Vector(a[0]*b[0], a[0]*b[1], a[0]*b[2]);
}
Vector operator*(const Vector& a, const double b) {
	return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator/(const Vector& a, const double b) {
	return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

Sphere::Sphere(Vector c, double r, Vector a, bool reflect) {
	this->C = c;
	this->R = r;
	this->albedo = a;
	this->Mirror = reflect;
}
struct Intersection Sphere::intersect(Ray& r) const {
    Vector OC = r.O - this->C;
    double b = 2 * dot(r.u, OC);
    double c = dot(OC, OC) - this->R * this->R;
    double delta = b * b - 4 * c;

    struct Intersection result;
    result.intersects = false;
    if (delta < 0) {
        return result;
    }

    double sqrtDelta = sqrt(delta);
    double t1 = (-b - sqrtDelta) / 2.0;
    double t2 = (-b + sqrtDelta) / 2.0;

    double t;
    if (t1 > 0 && t2 > 0) {
        t = std::min(t1, t2);
    } else if (t1 > 0) {
        t = t1;
    } else if (t2 > 0) {
        t = t2;
    } else {
        return result;
    }

    result.intersects = true;
    result.point = r.O + t * r.u;
    return result;
}


Light::Light(Vector p, double i) {
	this->P = p;
	this->I = i;
}

Ray::Ray(Vector o, Vector dir) {
	this->O = o;
	this->u = dir;
}

Scene::Scene(){
    this->arr.emplace_back(Vector(0, 0, 0), 10, Vector(0.8, 0.8, 0.8), true);
	this->arr.emplace_back(Vector(1000,0,0), 940, Vector(0.6, 0.5, 0.1));
	this->arr.emplace_back(Vector(-1000,0,0), 940, Vector(0.9, 0.2, 0.9));
	this->arr.emplace_back(Vector(0,0,-1000), 940, Vector(0.4, 0.8, 0.7));
	this->arr.emplace_back(Vector(0,1000,0), 940, Vector(0.2, 0.5, 0.9));
	this->arr.emplace_back(Vector(0,-1000,0), 990, Vector(0.3, 0.4, 0.7));
	this->arr.emplace_back(Vector(0,0,1000), 940,  Vector(0.9, 0.4, 0.3));
	this->lights.emplace_back(Light(Vector(-10,20,40), 1*pow(10,10)));
}



Vector Scene::get_colour(Ray& ray, int max_reflection) const {

	double dist = 100000;
	struct Intersection intersection;
	Sphere sphere = this->arr[0];
	intersection.intersects = false;
	for (auto s : this->arr) {
		struct Intersection tmp = s.intersect(ray);
		if (tmp.intersects && ((ray.O - tmp.point).norm() < dist)) {
			intersection = tmp;
			sphere = s;
			dist = (ray.O - tmp.point).norm();
		}
	}

	if (! intersection.intersects){
		return Vector(100, 100, 100);
	}

	Vector normal = intersection.point - sphere.C;
	normal.normalize();

	intersection.point = intersection.point + pow(10, -12)*normal;
	Vector LP = this->lights[0].P - intersection.point;

	dist = LP.norm();
	Ray r = Ray(intersection.point, LP/LP.norm());



	if (sphere.Mirror && max_reflection > 0){
		Ray mirror_ray = Ray(intersection.point, ray.u - 2*dot(ray.u, normal)*normal);
		return this->get_colour(mirror_ray, max_reflection);
	}

	// Check for shadow
	bool shadow = false;
	for (auto s : this->arr) {
		struct Intersection tmp = s.intersect(r);
		if (tmp.intersects && ((intersection.point - tmp.point).norm() < dist)) {
			shadow = true;
		}
	} 
	// Direct Light
	Vector direct_light = Vector();
	if (not shadow){
		Vector material = sphere.albedo/M_PI;
		double attenuation = this->lights[0].I / (4 * M_PI * LP.norm2());
		double angle = dot(normal, LP/LP.norm());
		direct_light =  angle*material*attenuation;
	}
	if (max_reflection == 0){return direct_light; }


	// Indirect Light
	Vector random_dir = random_vector(normal);
	Ray indirect_ray = Ray(intersection.point, random_dir);
	Vector indirect = sphere.albedo * get_colour(indirect_ray, max_reflection -1);

		
	return direct_light + indirect;
}

Vector random_vector(Vector& normal){

	// Sample Random Direction
		// Compute x, y, z
		double r1 = unif(engine);
		double r2 = unif(engine);

		double x = cos(2*M_PI*r1)*sqrt(1-r2);
		double y = sin(2*M_PI*r1)*sqrt(1-r2);
		double z = sqrt(r2);


		Vector T1 = Vector();

		if (abs(normal[2]) < abs(normal[0]) and abs(normal[2]) < abs(normal[1])){
			T1 = Vector(- normal[1], normal[0], 0);
			T1.normalize();
			Vector T2 = cross(normal,T1);
			return x*T1 +y*T2 + z*normal;
		} else if (abs(normal[0]) < abs(normal[1]) and abs(normal[0]) < abs(normal[2])){
			T1 = Vector(0, - normal[2], normal[1]);
			T1.normalize();
			Vector T2 = cross(normal,T1);
			return x*T1 +y*T2 + z*normal;
		}
		T1 = Vector(- normal[2], 0, normal[0]);
		T1.normalize();
		Vector T2 = cross(normal,T1);
		return x*T1 +y*T2 + z*normal;
}

int main() {
	int W = 512;
	int H = 512;
	double alpha = 60*M_PI / 180;
	Vector camera = Vector(0,0,55);
	int max_reflection = 7;
	std::vector<unsigned char> image(W * H * 3, 0);
	Scene scene = Scene();

	#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < H; i++) {
		// int tid = omp_get_thread_num();
		for (int j = 0; j < W; j++) {
			Vector colour = Vector();
			int N = 100;
			for (int n = 0; n < N; n++) {
				double r1 = unif(engine);
				double r2 = unif(engine);
				double x1 = sqrt(-2*log(r1))*cos(2.*M_PI*r2) * 0.25; // Standard deviation is 0.25
				double x2 = sqrt(-2*log(r1))*sin(2.*M_PI*r2) * 0.25;

				Vector dir = Vector(j - W/2 +0.5 +x1, H/2 - i - 0.5 +x2, (-W/(2.*tan(alpha/2))));
				dir.normalize();
				Ray r = Ray(camera,  dir);
				colour = colour + scene.get_colour(r, max_reflection);
			}
			colour = colour/N;
			

            image[(i * W + j) * 3 + 0] = std::max(0., std::min(255., std::pow(colour.data[0], 1.0/2.2)));
            image[(i * W + j) * 3 + 1] = std::max(0., std::min(255., std::pow(colour.data[1], 1.0/2.2)));
            image[(i * W + j) * 3 + 2] = std::max(0., std::min(255., std::pow(colour.data[2], 1.0/2.2)));

		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}