#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <iostream>
#include <random>
#include <string>
#include <stdio.h>
#include <algorithm>
#include "Vector.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

std::default_random_engine engine;
std::uniform_real_distribution<double> unif(0.0, 1.0);

void TriangleMesh::readOBJ(const char* obj) {
		char matfile[255];
		char grp[255];
	
		FILE* f;
		f = fopen(obj, "r");
		int curGroup = -1;
		while (!feof(f)) {
			char line[255];
			if (!fgets(line, 255, f)) break;
	
			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());
	
			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				curGroup++;
			}
	
			if (line[0] == 'v' && line[1] == ' ') {
				Vector vec;
	
				Vector col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
					col[0] = std::min(1., std::max(0., col[0]));
					col[1] = std::min(1., std::max(0., col[1]));
					col[2] = std::min(1., std::max(0., col[2]));
	
					vertices.push_back(vec);
					vertexcolors.push_back(col);
	
				} else {
					sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector vec;
				sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector vec;
				sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;
				t.group = curGroup;
	
				char* consumedline = line + 1;
				int offset;
	
				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
					if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
					if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
					if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
					indices.push_back(t);
				} else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
						if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
						if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
						if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
						indices.push_back(t);
					} else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
							indices.push_back(t);
						} else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
							if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
							if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
							if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
							indices.push_back(t);
						}
					}
				}
	
				consumedline = consumedline + offset;
	
				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.group = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
						if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
						if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
						if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					} else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
							if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
							if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
							if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						} else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
								if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
								if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
								if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;                             
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							} else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
									if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
									if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								} else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}
	
			}
	
		}
		fclose(f);
	
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
Vector operator-(const Vector& a) {
	return Vector(-a[0], -a[1], -a[2]);
}
Vector operator*(const double a, const Vector& b) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const Vector& b) {
	return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
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
Sphere::Sphere(Vector c, double r, Vector a, bool reflect, bool transparent, bool inner) {
	this->C = c;
	this->R = r;
	this->albedo = a;
	this->Mirror = reflect;
	this->Transparent = transparent;
	this->Inner = inner;
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
	result.normal = result.point - this->C;
	result.normal.normalize();
    return result;
}

void TriangleMesh::computeBoundingBox(){
    boundingBox.min = Vector(1E9, 1E9, 1E9);
    boundingBox.max = Vector(-1E9, -1E9, -1E9);

    for (const Vector& v : vertices) {
        if (v[0] < boundingBox.min[0]) boundingBox.min[0] = v[0];
        if (v[1] < boundingBox.min[1]) boundingBox.min[1] = v[1];
        if (v[2] < boundingBox.min[2]) boundingBox.min[2] = v[2];

        if (v[0] > boundingBox.max[0]) boundingBox.max[0] = v[0];
        if (v[1] > boundingBox.max[1]) boundingBox.max[1] = v[1];
        if (v[2] > boundingBox.max[2]) boundingBox.max[2] = v[2];
    }
}


bool BoundingBox::Intersect(Ray& r) const{
	double tx1 = (min[0] - r.O[0])/r.u[0];
	double tx2 = (max[0] - r.O[0])/r.u[0];
	double txMin = std::min(tx1, tx2);
	double txMax = std::max(tx1, tx2);

	double ty1 = (min[1] - r.O[1])/r.u[1];
	double ty2 = (max[1] - r.O[1])/r.u[1];
	double tyMin = std::min(ty1, ty2);
	double tyMax = std::max(ty1, ty2);

	double tz1 = (min[2] - r.O[2])/r.u[2];
	double tz2 = (max[2] - r.O[2])/r.u[2];
	double tzMin = std::min(tz1, tz2);
	double tzMax = std::max(tz1, tz2);

	if (std::min(txMax, std::min(tyMax, tzMax)) > std::max(txMin, std::max(tyMin, tzMin))){
		return true;
	}
	return false;
}

void TriangleMesh::scale_translate(double s, const Vector& t){
	for (int i = 0; i < vertices.size(); i++){
		vertices[i][0] = vertices[i][0]*s + t[0];
		vertices[i][1] = vertices[i][1]*s + t[1];
		vertices[i][2] = vertices[i][2]*s + t[2];
	}
}
struct Intersection TriangleMesh::intersect(Ray& r) const {
	struct Intersection intersection;
	intersection.intersects = false;
	if	(! boundingBox.Intersect(r)){
		return intersection;
	}
	double t_prev = 1E9;

	for (int i= 0; i < indices.size(); i++){
		const Vector& A = vertices[indices[i].vtxi];
		const Vector& B = vertices[indices[i].vtxj];
		const Vector& C = vertices[indices[i].vtxk];

		const Vector e1 = B - A;
		const Vector e2 = C - A;
		Vector N = cross(e1, e2);		
		const Vector AOu = cross((A-r.O),r.u);
		double invuN = 1/dot(r.u, N);
		double t = dot(A - r.O, N)*invuN;
		if (t < 0) continue;
		double beta = dot(e2, AOu)*invuN;
		if (beta < 0 or beta > 1) continue;
		double gamma = -dot(e1, AOu)*invuN;
		if (gamma < 0 or gamma > 1) continue;
		double alpha = 1 - beta - gamma;
		if (alpha < 0 or alpha > 1) continue;

		if (t > 0.000001 and t < t_prev){
			Vector P = A + beta*e1 + gamma*e2;
			t_prev = t;
			intersection.point = P;
			intersection.intersects = true;
			intersection.normal = alpha*normals[indices[i].ni] + beta*normals[indices[i].nj] + gamma*normals[indices[i].nk];
			intersection.normal.normalize();
		}
	}
	return intersection;
};

Light::Light(Vector p, double i) {
	this->P = p;
	this->I = i;
}
Ray::Ray(Vector o, Vector dir) {
	this->O = o;
	this->u = dir;
}
Scene::Scene(){
	this->arr.push_back(new Sphere(Vector(-20, 0, 0), 10, Vector(0.8, 0.8, 0.8), true));
    this->arr.push_back(new Sphere(Vector(0, 0, 0), 10, Vector(0.8, 0.8, 0.8), false, true));
	this->arr.push_back(new Sphere(Vector(20, 0, 0), 10, Vector(0.8, 0.8, 0.8), false, true));
    this->arr.push_back(new Sphere(Vector(20, 0, 0), 9.5, Vector(0.8, 0.8, 0.8), false, true, true));
	this->arr.push_back(new Sphere(Vector(1000,0,0), 940, Vector(0.6, 0.5, 0.1)));
	this->arr.push_back(new Sphere( Vector(-1000,0,0), 940, Vector(0.9, 0.2, 0.9)));
	this->arr.push_back(new Sphere( Vector(0,0,-1000), 940, Vector(0.4, 0.8, 0.7)));
	this->arr.push_back(new Sphere( Vector(0,1000,0), 940, Vector(0.2, 0.5, 0.9)));
	this->arr.push_back(new Sphere( Vector(0,-1000,0), 990, Vector(0.3, 0.4, 0.7)));
	this->arr.push_back(new Sphere( Vector(0,0,1000), 940,  Vector(0.9, 0.4, 0.3)));
	this->lights.emplace_back(Light(Vector(-10,20,40), 8*pow(10,9)));
}
Vector Scene::get_colour(Ray& ray, int max_reflection, double n1) const {

	double dist = 100000;
	struct Intersection intersection;
	const Geometry* object = this->arr[0];
	intersection.intersects = false;
	for (auto s : this->arr) {
		struct Intersection tmp = s->intersect(ray);
		if (tmp.intersects && ((ray.O - tmp.point).norm() < dist)) {
			intersection = tmp;
			object = s;
			dist = (ray.O - tmp.point).norm();
		}
	}

	if (! intersection.intersects){
		return Vector(0, 0, 0);
	}

	Vector normal = intersection.normal;
	normal.normalize();
	if (object->Inner){normal = -normal;}
	Vector LP = this->lights[0].P - intersection.point + normal*0.00001;

	dist = LP.norm();
	Ray r = Ray(intersection.point + (LP/LP.norm()*0.00001), LP/LP.norm());

	if (object->Transparent && (max_reflection > 0)){
		double n2 = 1.5;
		double iN = dot(ray.u, normal);
		Vector refraction_dir = Vector();
		if (iN > 0){ //Ray is exiting the sphere
			normal = -normal;
			std::swap(n1,n2);
			iN = dot(ray.u, normal);
		}
		double delta = 1 - pow(n1/n2, 2)*(1 - pow(iN, 2));

		if (delta < 0){ // Total internal reflection
			refraction_dir = ray.u - 2*dot(ray.u, normal)*normal;
			refraction_dir.normalize();
			Ray mirror_ray = Ray(intersection.point + 0.00001*refraction_dir, refraction_dir);
			return this->get_colour(mirror_ray, max_reflection -1, n1);
		}
		Vector tn = - sqrt(delta)*normal;
		Vector tt = (n1/n2)*(ray.u - iN*normal);
		refraction_dir = tt + tn;
		double k0 = pow((n1 - n2) / (n1 + n2), 2);
		double R = k0 + (1 - k0)*pow(1 - abs(iN), 5);
		double T = 1 - R;
		if (unif(engine) > T) {
			refraction_dir =  ray.u - 2*dot(ray.u, normal)*normal;
		}
		refraction_dir.normalize();
		Ray refraction_ray = Ray(intersection.point + 0.00001*refraction_dir, refraction_dir);
		return this->get_colour(refraction_ray, max_reflection -1, n1);
	}
	

	if (object->Mirror && max_reflection > 0){
		Ray mirror_ray = Ray(intersection.point, ray.u - 2*dot(ray.u, normal)*normal);
		// std::cout << "Mirror" << std::endl;
		return this->get_colour(mirror_ray, max_reflection -1, n1);
	}

	// Check for shadow
	bool shadow = false;
	for (auto s : this->arr) {
		struct Intersection tmp = s->intersect(r);
		if (tmp.intersects && ((intersection.point - tmp.point).norm() < dist)) {
			shadow = true;
		}
	} 
	// Direct Light
	Vector direct_light = Vector();
	if (not shadow){
		Vector material = object->albedo/M_PI;
		double attenuation = this->lights[0].I / (4 * M_PI * LP.norm2());
		double angle = dot(normal, LP/LP.norm());
		if (angle >= 0 ){
			direct_light =  angle*material*attenuation;
		}
	}
	if (max_reflection == 0){return direct_light; }


	// Indirect Light
	Vector random_dir = random_vector(normal);
	Ray indirect_ray = Ray(intersection.point + 0.00001*random_dir, random_dir);
	// std::cout << "Indirect" << std::endl;
	Vector indirect = object->albedo * get_colour(indirect_ray, max_reflection -1, n1);

		
	return direct_light; // + indirect;
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
	int W = 256;
	int H = 256;

	Scene scene = Scene();

	// TriangleMesh* mesh = new TriangleMesh();
	// mesh->readOBJ("cat.obj");
	// mesh->albedo = Vector(0.8, 0.8, 0.8);
	// mesh->scale_translate(0.6, Vector(0,-10,0));
	// mesh->computeBoundingBox();
	// scene.arr.push_back(mesh);


	double alpha = 60*M_PI / 180;
	Vector camera = Vector(0,0,55);

	int max_reflection = 5;
	int N = 100;

	std::vector<unsigned char> image(W * H * 3, 0);

	#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < H; i++) {
		// int tid = omp_get_thread_num();
		for (int j = 0; j < W; j++) {
			Vector colour = Vector();
			for (int n = 0; n < N; n++) {
				double r1 = unif(engine);
				double r2 = unif(engine);
				double x1 = sqrt(-2*log(r1))*cos(2.*M_PI*r2) * 0.25; // Standard deviation is 0.25
				double x2 = sqrt(-2*log(r1))*sin(2.*M_PI*r2) * 0.25;

				Vector dir = Vector(j - W/2 +0.5 +x1, H/2 - i - 0.5 +x2, (-W/(2.*tan(alpha/2))));
				dir.normalize();
				Ray r = Ray(camera,  dir);
				colour = colour + scene.get_colour(r, max_reflection, 1.0);
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