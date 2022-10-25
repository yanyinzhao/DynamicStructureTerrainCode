#include "algorithms.h"

int main(int argc, char **argv)
{

    // std::string input_file = "../input/EP_600996.off";
    // std::string input_file = "../input/SCpre_500000.off";
    // std::string input_file = "../input/VSpost.off";
    // std::string pre_input_file = "../exp_input/AUpre_2504322.off";
    // std::string post_input_file = "../exp_input/AUpost_2504322.off";
    std::string pre_input_file = "../exp_input/AUpre_2504322.off";
    std::string post_input_file = "../exp_input/AUpost_2504322.off";

    // std::vector<double> points;
    // std::vector<unsigned> faces;
    // geodesic::read_mesh_from_file(&input_file[0], points, faces);
    // geodesic::Mesh mesh;
    // mesh.initialize_mesh_data(points, faces);

    std::vector<double> pre_points;
    std::vector<unsigned> pre_faces;
    geodesic::read_mesh_from_file(&pre_input_file[0], pre_points, pre_faces);
    geodesic::Mesh pre_mesh;
    pre_mesh.initialize_mesh_data(pre_points, pre_faces);

    std::vector<double> post_points;
    std::vector<unsigned> post_faces;
    geodesic::read_mesh_from_file(&post_input_file[0], post_points, post_faces);
    geodesic::Mesh post_mesh;
    post_mesh.initialize_mesh_data(post_points, post_faces);

    std::vector<int> changed_face_index_list;
    changed_face_index_list.clear();
    assert(pre_mesh.faces().size() == post_mesh.faces().size());
    int face_difference_error = 100;
    for (int i = 0; i < pre_mesh.faces().size(); i++)
    {
        if (abs(pre_mesh.faces()[i].adjacent_vertices()[0]->x() - post_mesh.faces()[i].adjacent_vertices()[0]->x()) > face_difference_error ||
            abs(pre_mesh.faces()[i].adjacent_vertices()[0]->y() - post_mesh.faces()[i].adjacent_vertices()[0]->y()) > face_difference_error ||
            abs(pre_mesh.faces()[i].adjacent_vertices()[0]->z() - post_mesh.faces()[i].adjacent_vertices()[0]->z()) > face_difference_error ||
            abs(pre_mesh.faces()[i].adjacent_vertices()[1]->x() - post_mesh.faces()[i].adjacent_vertices()[1]->x()) > face_difference_error ||
            abs(pre_mesh.faces()[i].adjacent_vertices()[1]->y() - post_mesh.faces()[i].adjacent_vertices()[1]->y()) > face_difference_error ||
            abs(pre_mesh.faces()[i].adjacent_vertices()[1]->z() - post_mesh.faces()[i].adjacent_vertices()[1]->z()) > face_difference_error ||
            abs(pre_mesh.faces()[i].adjacent_vertices()[2]->x() - post_mesh.faces()[i].adjacent_vertices()[2]->x()) > face_difference_error ||
            abs(pre_mesh.faces()[i].adjacent_vertices()[2]->y() - post_mesh.faces()[i].adjacent_vertices()[2]->y()) > face_difference_error ||
            abs(pre_mesh.faces()[i].adjacent_vertices()[2]->z() - post_mesh.faces()[i].adjacent_vertices()[2]->z()) > face_difference_error)
        {
            changed_face_index_list.push_back(i);
            // std::cout << i << std::endl;
        }
    }
    std::cout << changed_face_index_list.size() << std::endl;

    // auto start_time = std::chrono::high_resolution_clock::now();
    // geodesic::GeodesicAlgorithmExact algorithm(&mesh);
    // double const distance_limit = geodesic::GEODESIC_INF;
    // geodesic::SurfacePoint source(&mesh.vertices()[0]);
    // geodesic::SurfacePoint destination(&mesh.vertices()[1]);
    // // geodesic::SurfacePoint source(&mesh.vertices()[216336]);
    // // geodesic::SurfacePoint destination(&mesh.vertices()[230408]);
    // std::vector<geodesic::SurfacePoint> one_source_poi_list(1, source);
    // std::vector<geodesic::SurfacePoint> one_destination_poi_list(1, destination);
    // algorithm.propagate(one_source_poi_list, distance_limit, &one_destination_poi_list);
    // double d;
    // algorithm.best_source(destination, d);
    // std::cout << d << std::endl;
    // auto stop_time = std::chrono::high_resolution_clock::now();
    // auto duration_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time);
    // double time = duration_time.count();
    // std::ofstream ofs("../output/output.txt", std::ios_base::app);
    // ofs << "time: " << time << " ms\n";
    // ofs.close();

    // auto start_time2 = std::chrono::high_resolution_clock::now();
    // geodesic::SurfacePoint source2(&mesh.vertices()[0]);
    // geodesic::SurfacePoint destination2(&mesh.vertices()[1]);
    // // geodesic::SurfacePoint source2(&mesh.vertices()[216336]);
    // // geodesic::SurfacePoint destination2(&mesh.vertices()[230408]);
    // std::vector<geodesic::SurfacePoint> one_source_poi_list2(1, source2);
    // std::vector<geodesic::SurfacePoint> one_destination_poi_list2(1, destination2);
    // algorithm.propagate(one_source_poi_list2);
    // algorithm.best_source(destination2, d);
    // std::cout << d << std::endl;
    // auto stop_time2 = std::chrono::high_resolution_clock::now();
    // auto duration_time2 = std::chrono::duration_cast<std::chrono::milliseconds>(stop_time2 - start_time2);
    // double time2 = duration_time2.count();
    // std::ofstream ofs2("../output/output.txt", std::ios_base::app);
    // ofs2 << "time: " << time2 << " ms\n";
    // ofs2.close();

    // std::vector<double> eps = {0.05, 0.1, 0.25, 0.5, 0.75, 1};
    // // std::vector<double> subdivision_level = {10, 20, 30, 40, 50, 60};
    // for (int i = 0; i < eps.size(); i++)
    // {
    //     auto start_time = std::chrono::high_resolution_clock::now();

    //     double subdivision_level = floor((std::pow(1.0 / (eps[i] + 0.5), 1.5) + 1.0) * (11.1 * eps[i] + 27.8));
    //     std::cout << "subdivision_level: " << subdivision_level << std::endl;
    //     geodesic::GeodesicAlgorithmSubdivision algorithm(&mesh, subdivision_level);
    //     double const distance_limit = geodesic::GEODESIC_INF;
    //     geodesic::SurfacePoint source(&mesh.vertices()[0]);
    //     geodesic::SurfacePoint destination(&mesh.vertices()[1]);
    //     std::vector<geodesic::SurfacePoint> one_source_poi_list(1, source);
    //     std::vector<geodesic::SurfacePoint> one_destination_poi_list(1, destination);
    //     algorithm.propagate(one_source_poi_list, distance_limit, &one_destination_poi_list);
    //     double d;
    //     algorithm.best_source(destination, d);

    //     auto stop_time = std::chrono::high_resolution_clock::now();
    //     auto duration_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time);
    //     double time = duration_time.count();

    //     std::ofstream ofs("../output/output.txt", std::ios_base::app);
    //     ofs << "subdivision_level: " << subdivision_level << ", time: " << time << " ms\n";
    //     ofs.close();
    // }
}