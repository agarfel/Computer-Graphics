#include <vector>

class Vector {
public:
	explicit Vector(double x = 0, double y = 0, double z = 0){
		data[0] = x;
		data[1] = y;
		data[2] = z;
	};
	double norm2() const;
	double norm() const;
	void normalize();
	double operator[](int) const;
	double& operator[](int);
	double data[3];
};


struct Intersection{
	bool intersects;
	Vector point;
	Vector normal;
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

class Geometry{
public:
	Geometry() : albedo(0, 0, 0), Transparent(false), Mirror(false) {};
	Geometry(const Vector& albedo, bool Mirror, bool Transparent) : albedo(albedo),Transparent(Transparent), Mirror(Mirror){};
	virtual Intersection intersect(Ray&) const = 0;
    Vector albedo;
	bool Mirror;
	bool Transparent;
};

class Sphere : public Geometry {
public:
	explicit Sphere(Vector, double, Vector, bool mirror = false, bool transparent = false);
    Vector C;
    double R;
	Intersection intersect(Ray&) const override;
};

class TriangleIndices {
	public:
		TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
		};
		int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
		int uvi, uvj, uvk;  // indices within the uv coordinates array
		int ni, nj, nk;  // indices within the normals array
		int group;       // face group
	};


class TriangleMesh : public Geometry {
public:
	~TriangleMesh() {};
	TriangleMesh() {};
	void scale_translate(double, const Vector&);
	Intersection intersect(Ray&) const override;
	void readOBJ(const char* obj);
	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
};


class Scene {
	public:
		explicit Scene();
        std::vector<const Geometry*> arr;
        struct Intersection intersect(Ray&) const;
        Vector get_colour(Ray&, int, double) const;
		std::vector<Light> lights;
};

Vector random_vector(Vector&);

int main();



/*
Bounding box implementation:
	1) Class BoundingBox
		public:
			BoundingBox(const Vector& m, const Vector& M)
			Vector Min, Max;
	2) add BoundingBox to triabgleMesh class

	3) add function:
	void computeBoundingBox(){
		bounding_box.m = Vector(1E9, 1E9, 1E9);
		bounding_box.m = Vector(-1E9, -1E9, -1E9);

		// Progressively try to improve bound
		for (int i = 0; i < vertoces.size(); i++){
			for (int j = 0; j < vertoces.size(); i++){
				for (int k = 0; k < vertoces.size(); i++){

			boundingbox.m[i] = std::min(bounding_box.,[0], vertives[i][o])
			boundingbox.M[j] = std::max(bounding_box.,[1], vertive[j][1])

		}
	}

	4) Compute bounding box after scaling mesh

	5) Wrtie bounding box intersect function:

	bool intersect(const Ray& r){
		double tx1 = (m[0] - r.O[0])/r.u[0];
		double tx2 = (M[0] - r.O[0])/r.u[0];
		double txMin = std::min(tx1, tx2);
		double txMax = std::max(tx1, tx2);

		double ty1 = (m[1] - r.O[1])/r.u[1];
		double ty2 = (M[1] - r.O[1])/r.u[1];
		double tyMin = std::min(ty1, ty2);
		double tyMax = std::max(ty1, ty2);

		double tz1 = (m[2] - r.O[2])/r.u[2];
		double tz2 = (M[2] - r.O[2])/r.u[2];
		double tzMin = std::min(tz1, tz2);
		double tzMax = std::max(tz1, tz2);

		if (std::min(txmax, tyMax, tzMax) > std::max(txmin, tyMin, tzMin))){
		return true
		}
		return false
	}

	6)In intersect, start by checking if hits bounding box:ADJ_OFFSET_SINGLESHOT
	if	(!boundingbox.)

*/