#include "svg.cpp"
#include "lbfgs.c"


std::default_random_engine engine;
std::uniform_real_distribution<double> unif(0.0, 1.0);

bool inside(Vector P0, Vector Pi, Vector X, double w0, double wi){
    return ((X - P0).norm2() - w0) <= ((X - Pi).norm2() - wi) + 1e-10;
};

class VornoiDiagram{
public:
    VornoiDiagram(){};
    Polygon clip_by_bisector(const Polygon& V, const Vector P0, const Vector Pi, double w0, double wi){
        Polygon result;
        for (int i = 0; i < V.vertices.size(); i++){
            const Vector& A = V.vertices[(i==0)? V.vertices.size()-1: i-1];
            const Vector& B = V.vertices[i];

            Vector D = Pi - P0;
            Vector M = (P0 + Pi) / 2 + ((w0 - wi) / (2 * D.norm2())) * D;
            
            double t = dot(M - A, Pi - P0)/dot(B - A, Pi - P0);
            Vector P = A + t*(B - A);


            if (inside(P0, Pi, B, w0, wi)){
                if (! inside(P0, Pi, A, w0, wi)){
                    result.vertices.push_back(P);
                }
                result.vertices.push_back(B);
            }
            else {
                if (inside(P0, Pi, A, w0, wi)){
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
        if (weights.empty()){
            weights = std::vector<double>(points.size(), 1.);
        }

        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < points.size(); i++){
            Polygon V = square;
            for (int j = 0; j < points.size(); j++){
                if(i==j) continue;
                double w0 = weights[i];
                double wi = weights[j];
                V = clip_by_bisector(V, points[i], points[j], w0, wi);
            }
            cells[i] = V;
        }
    }
    std::vector<Vector> points;
    std::vector<Polygon> cells;
    std::vector<double> weights;

};
static lbfgsfloatval_t evaluate(void* , const lbfgsfloatval_t*, lbfgsfloatval_t*, const int, const lbfgsfloatval_t);

class OptimalTransport{
public:
    OptimalTransport(){};
    VornoiDiagram vor;


    void optimize(){
        int i, ret = 0;
        int N = vor.points.size();
        lbfgsfloatval_t fx;
        std::vector<double> weights(N, 0);
        lbfgs_parameter_t param;

        lbfgs_parameter_init(&param);

        ret = lbfgs(N, &weights[0], &fx, evaluate, NULL, (void*) this, &param);
        memcpy(&vor.weights[0], &weights[0], N * sizeof(double));

        vor.compute();

    }
};

static lbfgsfloatval_t evaluate(void* instance, const lbfgsfloatval_t* x, lbfgsfloatval_t* g, const int n, const lbfgsfloatval_t step){
    OptimalTransport* ot = (OptimalTransport*)(instance);
    int N = ot->vor.points.size();
    memcpy(&ot->vor.weights[0], x, N*sizeof(x[0]));
    ot->vor.compute(); //compute vornoi diagram

    int i;
    lbfgsfloatval_t fx = 0.;
    for (int i = 0; i < n; i++){
        double current_area = ot->vor.cells[i].area();
        g[i] = current_area - (1.0 / n);
        fx += ot->vor.cells[i].integral_sqr_d(ot->vor.points[i]) - x[i]*(current_area - 1./n);

    }
    return -fx;
}

int main() {
    int N = 50;
    VornoiDiagram Vor;
    for (int i=0; i<N; i++){
        Vor.points.push_back(Vector(unif(engine), unif(engine), 0));
    }
    Vor.compute();
    OptimalTransport ot;
    ot.vor = Vor;
    ot.optimize();
    save_svg(ot.vor.cells, "test.svg");
    return 0;
}
