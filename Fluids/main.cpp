#include "svg.cpp"
std::default_random_engine engine;
std::uniform_real_distribution<double> unif(0.0, 1.0);

bool inside(Vector P0, Vector Pi, Vector X){
    return (X - P0).norm2() <= (X - Pi).norm2();
};

class VornoiDiagram{
public:
    VornoiDiagram(){};
    Polygon clip_by_bisector(const Polygon& V, const Vector P0, const Vector Pi){
        Polygon result;
        for (int i = 0; i < V.vertices.size(); i++){
            const Vector& A = V.vertices[(i==0)? V.vertices.size()-1: i-1];

            const Vector& B = V.vertices[i];

            Vector M = (P0 + Pi)/2;
            double t = dot(M - A, Pi - P0)/dot(B - A, Pi - P0);
            Vector P = A + t*(B - A);


            if (inside(P0, Pi, B)){
                if (! inside(P0, Pi, A)){
                    result.vertices.push_back(P);
                }
                result.vertices.push_back(B);
            }
            else {
                if (inside(P0, Pi, A)){
                    result.vertices.push_back(P);
                }
            }
        }
        return result;
    
    
    }
    void compute(){
        Polygon square;
        square.vertices.push_back(Vector(0,0));
        square.vertices.push_back(Vector(0,1));
        square.vertices.push_back(Vector(1,1));
        square.vertices.push_back(Vector(1,0));
        cells.resize(points.size());
        #pragma omp parallel for schedule(dynamic, 1)

        for (int i = 0; i < points.size(); i++){
            Polygon V = square;
            for (int j = 0; j < points.size(); j++){
                if(i==j) continue;
                V = clip_by_bisector(V, points[i], points[j]);
            }
            cells[i] = V;
        }
    }
    std::vector<Vector> points;
    std::vector<Polygon> cells;
};



int main() {
    int N = 100;
    VornoiDiagram Vor;
    for (int i=0; i<N; i++){
        Vor.points.push_back(Vector(unif(engine), unif(engine), 0));
    }
    Vor.compute();
    save_svg(Vor.cells, "test.svg");
    return 0;
}
