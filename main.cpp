#include <OpenMesh/Core/IO/importer/ImporterT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"

// These are required for functions like get_feature_edge
struct simpleTraits : public OpenMesh::DefaultTraits {
    VertexAttributes(OpenMesh::Attributes::Status);
    FaceAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Color);
    EdgeAttributes(OpenMesh::Attributes::Status);
};

typedef OpenMesh::TriMesh_ArrayKernelT<simpleTraits> simpleMesh;


// This function unpacks openmesh data to a simpler cpp containers suitable for polyscope
std::pair<std::vector<std::array<double, 3>>, std::vector<std::array<unsigned int, 3>>>
unpack_open_mesh(simpleMesh &reference_mesh) {
    std::vector<std::array<double, 3> > vertexPositions;
    std::vector<std::array<unsigned int, 3> > meshConnectivity;
    for (auto v: reference_mesh.vertices()) {
        auto current_point = reference_mesh.point(v);
        std::array<double, 3> current_array_point;
        current_array_point[0] = current_point[0];
        current_array_point[1] = current_point[1];
        current_array_point[2] = current_point[2];
        vertexPositions.push_back(current_array_point);
    }
    for (auto f: reference_mesh.faces()) {
        std::array<unsigned int, 3> current_triangle;
        unsigned int temp_index = 0;
        for (auto fv_it = reference_mesh.fv_iter(f); fv_it.is_valid(); ++fv_it) {
            current_triangle[temp_index] = fv_it->idx();
            temp_index++;
        }
        meshConnectivity.push_back(current_triangle);
        temp_index = 0;
    }
    return std::make_pair(vertexPositions, meshConnectivity);
}

// This function ensures that only the vertices that are part of feature network are present
// If extra vertices are present, polyscope results look weird
std::pair<std::vector<std::array<double, 3>>, std::vector<std::array<unsigned int, 2> > >
get_feature_edges(simpleMesh &reference_mesh, double feature_angle = 45.0) {
    std::vector<std::array<unsigned int, 2> > feature_edges_global;
    std::set<unsigned int> unique_vertex_ids;
    reference_mesh.find_feature_edges(feature_angle);
    std::map<unsigned int, std::array<double, 3> > global_id_to_vertex_positions;
    for (auto e: reference_mesh.edges()) {
        if (reference_mesh.status(e).feature()) {
            std::array<unsigned int, 2> current_feature_edge;
            current_feature_edge[0] = e.v0().idx();
            current_feature_edge[1] = e.v1().idx();
            feature_edges_global.push_back(current_feature_edge);
            unique_vertex_ids.insert(current_feature_edge[0]);
            unique_vertex_ids.insert(current_feature_edge[1]);
            auto vertex_one = reference_mesh.point(e.v0());
            auto vertex_two = reference_mesh.point(e.v1());
            std::array<double, 3> vertex_one_array = std::array<double, 3>{
                    vertex_one[0], vertex_one[1], vertex_one[2]
            };
            std::array<double, 3> vertex_two_array = std::array<double, 3>{
                    vertex_two[0], vertex_two[1], vertex_two[2]
            };
            global_id_to_vertex_positions[current_feature_edge[0]] = vertex_one_array;
            global_id_to_vertex_positions[current_feature_edge[1]] = vertex_two_array;
        }
    }
    std::map<unsigned int, unsigned int> global_id_to_local_id;
    std::map<unsigned int, unsigned int> local_id_to_global_id;
    unsigned int new_id = 0;
    for (auto v: unique_vertex_ids) {
        global_id_to_local_id[v] = new_id;
        local_id_to_global_id[new_id] = v;
        new_id++;
    }
    std::vector<std::array<double, 3> > new_positions;
    std::vector<std::array<unsigned int, 2> > new_feature_edges;
    for (unsigned int i = 0; i < new_id; ++i) {
        new_positions.push_back(global_id_to_vertex_positions[local_id_to_global_id[i]]);
    }
    for (unsigned int i = 0; i < feature_edges_global.size(); ++i) {
        std::array<unsigned int, 2> current_feature_edge;
        current_feature_edge[0] = global_id_to_local_id[feature_edges_global.at(i)[0]];
        current_feature_edge[1] = global_id_to_local_id[feature_edges_global.at(i)[1]];
        new_feature_edges.push_back(current_feature_edge);
    }
    return std::make_pair(new_positions, new_feature_edges);
}

// This is quite useless, try to do some useful computation to visualize
void initialize_random_scalar_field(simpleMesh &reference_mesh, const std::string field_name) {
    OpenMesh::VPropHandleT<double> field_handle;
    reference_mesh.add_property(field_handle, field_name);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distribution(0, 1); // define the range
    for (auto vertex_id: reference_mesh.vertices()) {
        if (reference_mesh.status(vertex_id).deleted()) {
            continue;
        }
        reference_mesh.property(field_handle, vertex_id) = distribution(gen);
    }
}

// This can be extended to other entities like edge or face
template<typename T>
std::vector<T> get_vertex_scalar_quantity(simpleMesh &reference_mesh, const std::string identifier) {
    std::vector<T> scalar_quantity;
    OpenMesh::VPropHandleT<T> scalar_handle;
    reference_mesh.get_property_handle(scalar_handle, identifier);
    for (auto v: reference_mesh.vertices()) {
        scalar_quantity.push_back(reference_mesh.property(scalar_handle, v));
    }
    return scalar_quantity;
}


int main() {
    polyscope::init();

    // Let's read a mesh and draw it in polyscope
    simpleMesh input_mesh;
    OpenMesh::IO::read_mesh(input_mesh, "../data/Bot_Shoulder.stl");
    auto unpackOpenMesh = unpack_open_mesh(input_mesh);
    polyscope::registerSurfaceMesh("openMesh", unpackOpenMesh.first,
                                   unpackOpenMesh.second);

    // Now detect the feature edges in openmesh and draw the same in polyscope
    auto featureData = get_feature_edges(input_mesh, 45.0);
    if (featureData.first.size() > 0) {
        polyscope::registerCurveNetwork("feature_edges", featureData.first, featureData.second);
    }

    // Here is an example for visualizing a scalar quantity
    initialize_random_scalar_field(input_mesh, "RANDOM_SCALAR");
    auto random_scalar = get_vertex_scalar_quantity<double>(input_mesh, "RANDOM_SCALAR");
    polyscope::getSurfaceMesh("openMesh")->addVertexScalarQuantity("VERTEX_SCALAR", random_scalar);

    // Finally show it all
    polyscope::show();
}