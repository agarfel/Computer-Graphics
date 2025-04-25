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

class BoundingBox{
	public:
		BoundingBox(): min(0,0,0), max(0,0,0){};
		Vector min, max;
		bool Intersect(Ray&) const;
};

class BVHNode{
public:
	BoundingBox BBox;
	BVHNode* left_child = nullptr;
	BVHNode* right_child = nullptr;
	int first, last;
	void construct(int first, int last);
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
	BoundingBox computeBoundingBox(int, int);
	void build_bvh(BVHNode*, int, int);

	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
	BoundingBox boundingBox;
	BVHNode root;
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
Implement:
	1) Construct BVH
	2) Intersection Routine

construct (BVHNode *c, int first, int last) {
	c->first = first;
	c->last = last;
	c->BBox = computeBBox(first, last) [compute bounding box for all triangles between first and last]
	axis = longest axis(c->BBox)  [Check which axis to split along]
	M = Middle(BBox, axis) [find middle point along axis]
	pivot = first [keep track of pivot, progressively swap with current value is value is smaller than axis]
	for (i = first; i < last; i++){
		if (bary(triangle[i]) < M) {
			swap:: Triangle[i], Triangle[pivot];
			pivot++;
		} 
	}
	construct(c->left, first, pivot);
	construct(c->right, pivot +1, last);	[rescursive calls]
}
BVHNode{
	BoundingBox BBox(of current node)
	BVHNode* left_child
	BVHNode* right_child
	int first, last [array index of first and last elements]
}

Depth first traversal of bounding box tree

intersect(ray){
	if (!intersect(BBox)){return false}
	list<BVHNode> l;
	l->push(root);
	while (!l->empty()) {
		c = l->back;
		l->pop(back)
		if (!leaf){
			if(l->left->BBox.intersect(ray)){
				l->push(c->left)
			}
			if(l->right->BBox.intersect(ray)){
				l->push(c->right)
			}
		} else {
			for (i = c->first; i < c->last; i++) {
				as before [Moller algorithm]
			} 
		}
	}
}


std::vector<std::vector<double<< textures;
std::vector<int> texturesW
std::vector<int> texturesh

void load_texture (const char* filename){
	int W,H,C;
	unsigned char* tex = stbi_load(filename, &W, &H, &C, 3);
	testuresW.pushback(W);
	testuresH.pushback(H);
	std::vector<double< current_tex(W*H*3);
	for (i=0; i< W*H*3; i++){
		current_tex[i] = std::pow(tex[i]/255., 2.2);
	}
}



after mesh.readOBJ:
mesh.load_texture("texture.png");


when getting intersection point:
	return albedo of texture?


	// for mesh intersection:
	Vector UV = alpha* (same as normal with uv instead of n)
	UV[0] = UV[0]*std::min(texturesW[indices[i].group], texturesW[indices[i].group] -1.); //add also >0
	UV[1] = (1 - UV[1])*std::min(texturesH[indices[i].group], texturesH[indices[i].group] -1.);
	int U = UV[0]
	int V = UV[1]
	int texHeight = texturesH[indices[i].group]
	VEctor albedo = Vector(
		textures[indices[i].group][3*(V*texHeight + U) + 1]
		textures[indices[i].group][3*(V*texHeight + U) + 2]
		textures[indices[i].group][3*(V*texHeight + U) + 3]
		)
*/