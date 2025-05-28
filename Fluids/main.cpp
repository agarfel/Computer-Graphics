#include "svg.cpp"
#include "lbfgs.c"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <sstream>

#define VOL_FLUID 0.6

std::default_random_engine engine;
std::uniform_real_distribution<double> unif(0.0, 1.0);

int sgn(double x){
    if (x > 0){return 1;}
    if (x < 0){return -1;}
    return 0;
}
void save_frame(const std::vector<Polygon> &cells, std::string filename, int frameid = 0) {
    int W = 1000, H = 1000;
    std::vector<unsigned char> image(W*H * 3, 255);
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < cells.size(); i++) {

        double bminx = 1E9, bminy = 1E9, bmaxx = -1E9, bmaxy = -1E9;
        for (int j = 0; j < cells[i].vertices.size(); j++) {
            bminx = std::min(bminx, cells[i].vertices[j][0]);
            bminy = std::min(bminy, cells[i].vertices[j][1]);
            bmaxx = std::max(bmaxx, cells[i].vertices[j][0]);
            bmaxy = std::max(bmaxy, cells[i].vertices[j][1]);
        }
        bminx = std::min(W-1., std::max(0., W * bminx));
        bminy = std::min(H-1., std::max(0., H * bminy));
        bmaxx = std::max(W-1., std::max(0., W * bmaxx));
        bmaxy = std::max(H-1., std::max(0., H * bmaxy));

        for (int y = bminy; y < bmaxy; y++) {
            for (int x = bminx; x < bmaxx; x++) {
                int prevSign = 0;
                bool isInside = true;
                double mindistEdge = 1E9;
                for (int j = 0; j < cells[i].vertices.size(); j++) {
                    double x0 = cells[i].vertices[j][0] * W;
                    double y0 = cells[i].vertices[j][1] * H;
                    double x1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][0] * W;
                    double y1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][1] * H;
                    double det = (x - x0)*(y1-y0) - (y - y0)*(x1-x0);
                    int sign = sgn(det);
                    if (prevSign == 0) prevSign = sign; else
                        if (sign == 0) sign = prevSign; else
                        if (sign != prevSign) {
                            isInside = false;
                            break;
                        }
                    prevSign = sign;
                    double edgeLen = sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0));
                    double distEdge = std::abs(det)/ edgeLen;
                    double dotp = (x - x0)*(x1 - x0) + (y - y0)*(y1 - y0);
                    if (dotp<0 || dotp>edgeLen*edgeLen) distEdge = 1E9;
                    mindistEdge = std::min(mindistEdge, distEdge);
                }
                if (isInside) {
                    //if (i < N) {   // the N first particles may represent fluid, displayed in blue
                    //  image[((H - y - 1)*W + x) * 3] = 0;
                    //  image[((H - y - 1)*W + x) * 3 + 1] = 0;
                    //  image[((H - y - 1)*W + x) * 3 + 2] = 255;
                    //}
                    if (mindistEdge <= 2) {
                        image[((H - y - 1)*W + x) * 3] = 0;
                        image[((H - y - 1)*W + x) * 3 + 1] = 0;
                        image[((H - y - 1)*W + x) * 3 + 2] = 0;
                    }

                }
                
            }
        }
    }
    std::ostringstream os;
    os << filename << frameid << ".png";
    stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
};

bool inside(Vector P0, Vector Pi, Vector X, double w0, double wi){
    return ((X - P0).norm2() - w0) <= ((X - Pi).norm2() - wi) + 1e-10;
};

class VornoiDiagram{
public:
    VornoiDiagram(){
        n_disk = 100;
        unit_disk.resize(n_disk);
        for (int i = 0; i < n_disk; i++){
            double theta = i*2.*M_PI / ((double)n_disk);
            unit_disk[i] = Vector(-cos(theta), sin(theta));
        }
    };
    std::vector<Vector> unit_disk;
    int n_disk;
    Polygon clip_by_edge(const Polygon& V, const Vector u, const Vector v){
        Polygon result;
        Vector N(-(v[1]-u[1]), v[0]-u[0], 0);

        for (int i = 0; i < V.vertices.size(); i++){
            const Vector& A = V.vertices[(i==0)? V.vertices.size()-1: i-1];
            const Vector& B = V.vertices[i];
            double tmp = dot(B - A, N);
            Vector P = A + (dot(u - A, N)/dot(B - A, N))*(B - A);

            if ( dot(u - B, N) >= 0 ){  //B INSIDE
                if ( dot(u - A, N) < 0 ){   // A OUTSIDE
                    result.vertices.push_back(P);
                }
                result.vertices.push_back(B);
            }
            else {
                if ( dot(u - A, N) >= 0 ){  // A INSIDE
                    result.vertices.push_back(P);
                }
            }
        }
        return result;
    }
    Polygon clip_by_bisector(const Polygon& V, const Vector P0, const Vector Pi, double w0, double wi){
        Polygon result;
        for (int i = 0; i < V.vertices.size(); i++){
            const Vector& A = V.vertices[(i==0)? V.vertices.size()-1: i-1];
            const Vector& B = V.vertices[i];

            Vector D = Pi - P0;
            Vector M = (P0 + Pi) / 2 + ((w0 - wi) / (2 * D.norm2())) * D;
            double tmp = dot(B - A, Pi - P0);
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
            Polygon V_prev;

            for (int j = 0; j < points.size(); j++){
                if(i==j) continue;
                double w0 = weights[i];
                double wi = weights[j];
                V_prev = V;
                V = clip_by_bisector(V, points[i], points[j], w0, wi);

            }

            // Clip V by disk
            double radius = sqrt(weights[i] - weights[weights.size()-1]);

            for (int j=0; j < n_disk; j++){
                Vector u = unit_disk[j]*radius + points[i];
                Vector v = unit_disk[(j+1)%n_disk]*radius + points[i];
                V = clip_by_edge(V, u,v);
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
        int N = vor.weights.size(); //N is number particles + 1
        lbfgsfloatval_t fx;
        std::vector<double> weights(N, 0);
        memcpy(&weights[0], &vor.weights[0], N*sizeof(weights[0]));

        lbfgs_parameter_t param;
        lbfgs_parameter_init(&param);

        ret = lbfgs(N, &weights[0], &fx, evaluate, NULL, (void*) this, &param);
        memcpy(&vor.weights[0], &weights[0], (N-1) * sizeof(weights[0]));

        vor.compute();

    }
};

static lbfgsfloatval_t evaluate(void* instance, const lbfgsfloatval_t* x, lbfgsfloatval_t* g, const int n, const lbfgsfloatval_t step){
    OptimalTransport* ot = (OptimalTransport*)(instance);
    int N = ot->vor.weights.size();
    memcpy(&ot->vor.weights[0], x, N*sizeof(x[0]));
    ot->vor.compute(); //compute vornoi diagram

    int i;
    lbfgsfloatval_t fx = 0.;
    double sum_areas = 0;
    for (int i = 0; i < N-1; i++){ //only particles
        double current_area = ot->vor.cells[i].area();
        g[i] = current_area - (VOL_FLUID / (N-1)); // 1/n --> VolFluids/ (N-1)
        fx += ot->vor.cells[i].integral_sqr_d(ot->vor.points[i]) - x[i]*(current_area - (VOL_FLUID / (N-1)));
        sum_areas += current_area;
    }
    double est_air_vol = 1 - sum_areas;
    // std::cout << est_air_vol << std::endl;
    g[N-1] = - ((1-VOL_FLUID) - est_air_vol);   //desired vol air - estimated vol air
    fx += x[N-1]*((1-VOL_FLUID) - est_air_vol); //
    return -fx;
}


class Fluid {
public:

    Fluid(int N = 1000):N(N){
        particles.resize(N);
        velocities.resize(N, Vector(0,0,0));
        for (int i = 0; i < N; i++){
            particles[i] = Vector(unif(engine), unif(engine), 0);
        }
        fluid_volume = VOL_FLUID;
        ot.vor.points = particles;
        ot.vor.weights.resize(N+1);
        std::fill(ot.vor.weights.begin(), ot.vor.weights.end(), 1.); //we want wi > wair. wi = 1, wair = 0.9
        ot.vor.weights[N] = 0.9;
    };

    void time_step(double dt){
        ot.vor.points = particles;
        ot.optimize();
        double epsilon2 = 0.004*0.004;
        Vector g(0, -9.81,0);
        double m_i = 200;
        for (int i = 0; i < particles.size(); i++){
            Vector center_cell = ot.vor.cells[i].centroid();
            Vector spring_force = (center_cell - particles[i])/epsilon2;
            Vector all_forces = g*m_i + spring_force;
            velocities[i] = velocities[i] + dt/m_i * all_forces;
            particles[i] = particles[i] + dt*velocities[i];
        }
    }

        void run_simulation(){
            double dt = 0.005;
            for (int i = 0; i < 150; i++){
                time_step(dt);
                save_frame(ot.vor.cells, "test/test", i);
            }
        }


    OptimalTransport ot;
    std::vector<Vector> particles;
    std::vector<Vector> velocities;
    int N;
    double fluid_volume;
};




int main() {

    Fluid fluid(1000);
    fluid.run_simulation();

    // int N = 50;
    // VornoiDiagram Vor;
    // for (int i=0; i<N; i++){
    //     Vor.points.push_back(Vector(unif(engine), unif(engine), 0));
    // }
    // Vor.compute();
    // OptimalTransport ot;
    // ot.vor = Vor;
    // ot.optimize();
    // save_svg(ot.vor.cells, "test.svg", "lightblue");
    return 0;
}