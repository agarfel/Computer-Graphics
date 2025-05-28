#include "obj_reader.cpp"
#include <map>
#include <set>
#include <random>


void TriangleMesh::write_obj(const char* filename) {
    FILE* f = fopen(filename, "w+");
    for (int i = 0; i < vertices.size(); i++) {
        fprintf(f, "v %3.5f %3.5f %3.5f\n", vertices[i][0], vertices[i][1], vertices[i][2]);
    }
    for (int i = 0; i < indices.size(); i++) {
        fprintf(f, "f %u %u %u\n", indices[i].vtxi+1, indices[i].vtxj + 1, indices[i].vtxk + 1);
    }
    fclose(f);

}
/*
Copy vector and triangle mesh and indices classes

*/
void TriangleMesh::tutte(){
    // FInd boundary veritces
    std::vector<std::set<int>> vertex_to_neighbors_set(vertices.size());

    for (int i = 0; i < indices.size(); i++){ // AGAIN NOT EFFICIENT - TRY TO IMPROVE
        vertex_to_neighbors_set[indices[i].vtxi].insert(indices[i].vtxj);
        vertex_to_neighbors_set[indices[i].vtxi].insert(indices[i].vtxk);
        vertex_to_neighbors_set[indices[i].vtxj].insert(indices[i].vtxi);
        vertex_to_neighbors_set[indices[i].vtxj].insert(indices[i].vtxk);
        vertex_to_neighbors_set[indices[i].vtxk].insert(indices[i].vtxi);
        vertex_to_neighbors_set[indices[i].vtxk].insert(indices[i].vtxj);
    }
    std::vector<std::vector<int>> vertex_to_neighbors(vertices.size());
    for (int i = 0; i < vertices.size(); i++){
        for (auto it : vertex_to_neighbors_set[i]){
            vertex_to_neighbors[i].push_back(it);
        }
    }

    std::vector<bool> boundary_vertex(vertices.size(), false);


    std::map<std::pair<int,int>, std::vector<int>> edge_to_triangles;
    for (int i = 0; i<indices.size(); i++){

        int minIJ = std::min(indices[i].vtxi, indices[i].vtxj);
        int maxIJ = std::max(indices[i].vtxi, indices[i].vtxj);
    
        std::pair<int,int> edge1(minIJ, maxIJ); //enforce an order [edge(i,j) <==> edge(j,i)]
        edge_to_triangles[edge1].push_back(i);

        int minIK = std::min(indices[i].vtxi, indices[i].vtxk);
        int maxIK = std::max(indices[i].vtxi, indices[i].vtxk);
    
        std::pair<int,int> edge2(minIK, maxIK); //enforce an order [edge(i,j) <==> edge(j,i)]
        edge_to_triangles[edge2].push_back(i);

        int minJK = std::min(indices[i].vtxj, indices[i].vtxk);
        int maxJK = std::max(indices[i].vtxj, indices[i].vtxk);
    
        std::pair<int,int> edge3(minJK, maxJK); //enforce an order [edge(i,j) <==> edge(j,i)]
        edge_to_triangles[edge3].push_back(i);
    }

    std::vector<std::pair<int,int>> edges_on_the_boundary;
    for (auto it = edge_to_triangles.begin(); it != edge_to_triangles.end(); it++){
        if (it->second.size() == 1){ //edge of boundary
            edges_on_the_boundary.push_back(it->first);
            boundary_vertex[it->first.first] = true;
            boundary_vertex[it->first.second] = true;
        }
    }

    // We have al edges in the boundary. Now we want to "order" them
    std::vector<int> boundary_chain;

    int starting_vertex = edges_on_the_boundary[0].first;
    boundary_chain.push_back(starting_vertex);
    int current_vertex = edges_on_the_boundary[0].second;
    int current_edge = 0;
    // THIS IS VERY UGLY. TRY TO IMPROVE
    while(current_vertex != starting_vertex){
        boundary_chain.push_back(current_vertex);
        for (int i = 0; i < edges_on_the_boundary.size(); i++){
            if (i != current_edge) {
                if (edges_on_the_boundary[i].first == current_vertex) {
                    current_vertex = edges_on_the_boundary[i].second;
                    current_edge = i;
                    break;
                } else if (edges_on_the_boundary[i].second == current_vertex) {
                    current_vertex = edges_on_the_boundary[i].first;
                    current_edge = i;
                    break;
                }
            }
            
        }
    }

    // Start parametrization
    TriangleMesh parametrization;
    parametrization.indices = indices; //indices of original mesh
    parametrization.vertices = vertices; //better if you initialise them close to unit circle

    // Initialise close to unit circel:
    // for (int i =0; i < vertices.size(); i++){
    //     parametrization.vertices[i] = Vector(rand()((double)RAND_MAX)*2.-1, (rand()/(double)RAND_M));
    // }

    for (int i = 0; i < edges_on_the_boundary.size(); i++){
        double angle = i*2.*M_PI/((double) boundary_chain.size());
        double x = cos(angle);
        double y = sin(angle);
        parametrization.vertices[boundary_chain[i]] = Vector(x,y,0); // Something on the circle
    }

    // Iterate
    for (int iter=0; iter < 1000; iter++){
        std::vector<Vector> new_vertices(vertices.size());
        for (int i= 0; i<vertices.size(); i++){
            if (!boundary_vertex[i]){ //if not a boundary vertex
                // new position is average of neighbours
                Vector avg(0,0,0);
                int num_neighbours = vertex_to_neighbors[i].size();
                for (int j=0; j<num_neighbours; j++){
                    avg = avg + parametrization.vertices[vertex_to_neighbors[i][j]];
                }
                avg = avg / num_neighbours;
                new_vertices[i] = avg;
            } else {
                new_vertices[i] = parametrization.vertices[i];
            }
        }
        parametrization.vertices = new_vertices;
    }

    parametrization.write_obj("tutte_result.obj");

}


int main() {
	TriangleMesh* mesh = new TriangleMesh();
	mesh->readOBJ("goethe.obj");
    mesh->tutte();
    return 0;
}