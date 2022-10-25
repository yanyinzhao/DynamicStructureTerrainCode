#include "algorithms.h"

void pre_CG_preprocessing_post_CG_updating(
    int poi_num, geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
    geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list,
    std::vector<std::vector<double>> &pre_pairwise_distance_poi_to_poi,
    std::vector<std::vector<double>> &post_pairwise_distance_poi_to_poi,
    std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> &pre_pairwise_path_poi_to_poi,
    std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> &post_pairwise_path_poi_to_poi,
    double &pre_preprocessing_time, double &pre_memory_usage,
    double &post_updating_time, double &post_memory_usage)
{
    std::vector<std::vector<std::vector<int>>> pre_face_sequence_index_list;
    std::vector<std::vector<double>> pairwise_distance_poi_to_poi;
    std::vector<std::vector<bool>> pairwise_distance_poi_to_poi_changed;
    std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> pairwise_path_poi_to_poi;
    std::vector<std::vector<double>> pairwise_distance_poi_to_vertex(poi_num, std::vector<double>(pre_mesh->vertices().size(), 0));

    pre_complete_graph_preprocessing(poi_num, pre_mesh, pre_poi_list, pre_face_sequence_index_list,
                                     pairwise_distance_poi_to_poi, pairwise_distance_poi_to_poi_changed,
                                     pairwise_path_poi_to_poi, pairwise_distance_poi_to_vertex,
                                     pre_preprocessing_time, pre_memory_usage);

    pre_pairwise_distance_poi_to_poi = pairwise_distance_poi_to_poi;
    pre_pairwise_path_poi_to_poi = pairwise_path_poi_to_poi;

    post_complete_graph_updating(poi_num, pre_mesh, pre_poi_list, post_mesh, post_poi_list,
                                 pre_face_sequence_index_list, pairwise_distance_poi_to_poi,
                                 pairwise_distance_poi_to_poi_changed,
                                 pairwise_path_poi_to_poi, pairwise_distance_poi_to_vertex,
                                 post_updating_time, post_memory_usage);

    post_pairwise_distance_poi_to_poi = pairwise_distance_poi_to_poi;
    post_pairwise_path_poi_to_poi = pairwise_path_poi_to_poi;
}

void UE_N1_and_UE_N2_and_UE_exp(
    int poi_num, geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
    geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list, double epsilon,
    int source_poi_index, int destination_poi_index, double pre_exact_distance,
    double post_exact_distance, int &pre_MST_weight, int &post_MST_weight,
    std::string write_file_header)
{
    std::cout << "== UE_N1 ==" << std::endl;

    std::vector<std::vector<double>> pre_pairwise_distance_poi_to_poi;
    std::vector<std::vector<double>> post_pairwise_distance_poi_to_poi;
    std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> pre_pairwise_path_poi_to_poi;
    std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> post_pairwise_path_poi_to_poi;
    double pre_preprocessing_time = 0;
    double pre_memory_usage = 0;
    double post_updating_time = 0;
    double post_memory_usage = 0;

    pre_CG_preprocessing_post_CG_updating(
        poi_num, pre_mesh, pre_poi_list, post_mesh, post_poi_list,
        pre_pairwise_distance_poi_to_poi, post_pairwise_distance_poi_to_poi,
        pre_pairwise_path_poi_to_poi, post_pairwise_path_poi_to_poi,
        pre_preprocessing_time, pre_memory_usage, post_updating_time, post_memory_usage);

    double pre_hash_mapping_time_UE_N1 = 0;
    double pre_query_time_UE_N1 = 0;
    double pre_hash_mapping_memory_usage_UE_N1 = 0;
    double pre_index_size_UE_N1 = 0;
    int pre_index_edge_num_UE_N1 = 0;
    double pre_index_weight_UE_N1 = 0;
    double pre_approximate_distance_UE_N1 = 0;
    std::vector<geodesic::SurfacePoint> pre_approximate_path_UE_N1;
    pre_approximate_path_UE_N1.clear();

    std::unordered_map<int, double> pairwise_distance_poi_to_poi_pre_map;
    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> pairwise_path_poi_to_poi_pre_map;
    get_pairwise_distance_and_path_poi_to_poi_map(pre_pairwise_distance_poi_to_poi, pairwise_distance_poi_to_poi_pre_map,
                                                  pre_pairwise_path_poi_to_poi, pairwise_path_poi_to_poi_pre_map,
                                                  pre_index_size_UE_N1, pre_index_edge_num_UE_N1, pre_index_weight_UE_N1,
                                                  pre_hash_mapping_time_UE_N1, pre_hash_mapping_memory_usage_UE_N1);
    complete_graph_query(poi_num, pairwise_distance_poi_to_poi_pre_map, pairwise_path_poi_to_poi_pre_map, source_poi_index,
                         destination_poi_index, pre_approximate_distance_UE_N1, pre_approximate_path_UE_N1, pre_query_time_UE_N1);
    calculate_MST_weight(pre_pairwise_distance_poi_to_poi, pre_MST_weight);

    // for (int i = 0; i < pre_pairwise_distance_poi_to_poi.size(); i++)
    // {
    //     for (int j = 0; j < pre_pairwise_distance_poi_to_poi[i].size(); j++)
    //     {
    //         std::cout << pre_pairwise_distance_poi_to_poi[i][j] << std::endl;
    //     }
    // }

    std::cout << "Pre terrain preprocessing time (CG): " << pre_preprocessing_time << " ms" << std::endl;
    std::cout << "Pre terrain preprocessing time (hash mapping): " << pre_hash_mapping_time_UE_N1 << " ms" << std::endl;
    std::cout << "Pre terrain query time: " << pre_query_time_UE_N1 << " ms" << std::endl;
    std::cout << "Pre terrain memory usage (CG): " << pre_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain memory usage (hash mapping): " << pre_hash_mapping_memory_usage_UE_N1 / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index size: " << pre_index_size_UE_N1 / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index edge number: " << pre_index_edge_num_UE_N1 << " edges" << std::endl;
    std::cout << "Pre terrain index/MST weight: " << pre_index_weight_UE_N1 / pre_MST_weight << std::endl;
    std::cout << "Pre terrain approximate distance: " << pre_approximate_distance_UE_N1 << ", pre terrain exact distance: " << pre_exact_distance << ", distance error: " << pre_approximate_distance_UE_N1 / pre_exact_distance - 1 << std::endl;
    std::cout << std::endl;

    std::ofstream ofs1("../output/output.txt", std::ios_base::app);
    ofs1 << "== UE_N1 ==\n";
    ofs1 << write_file_header << "\t"
         << pre_preprocessing_time << "\t"
         << pre_hash_mapping_time_UE_N1 << "\t"
         << pre_query_time_UE_N1 << "\t"
         << pre_memory_usage / 1e6 << "\t"
         << pre_hash_mapping_memory_usage_UE_N1 / 1e6 << "\t"
         << pre_index_size_UE_N1 / 1e6 << "\t"
         << pre_index_edge_num_UE_N1 << "\t"
         << pre_index_weight_UE_N1 / pre_MST_weight << "\t"
         << pre_approximate_distance_UE_N1 / pre_exact_distance - 1 << "\t";
    ofs1.close();

    double post_hash_mapping_time_UE_N1 = 0;
    double post_query_time_UE_N1 = 0;
    double post_hash_mapping_memory_usage_UE_N1 = 0;
    double post_index_size_UE_N1 = 0;
    int post_index_edge_num_UE_N1 = 0;
    double post_index_weight_UE_N1 = 0;
    double post_approximate_distance_UE_N1 = 0;
    std::vector<geodesic::SurfacePoint> post_approximate_path_UE_N1;
    post_approximate_path_UE_N1.clear();

    std::unordered_map<int, double> pairwise_distance_poi_to_poi_post_map;
    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> pairwise_path_poi_to_poi_post_map;
    get_pairwise_distance_and_path_poi_to_poi_map(post_pairwise_distance_poi_to_poi, pairwise_distance_poi_to_poi_post_map,
                                                  post_pairwise_path_poi_to_poi, pairwise_path_poi_to_poi_post_map,
                                                  post_index_size_UE_N1, post_index_edge_num_UE_N1, post_index_weight_UE_N1,
                                                  post_hash_mapping_time_UE_N1, post_hash_mapping_memory_usage_UE_N1);
    complete_graph_query(poi_num, pairwise_distance_poi_to_poi_post_map, pairwise_path_poi_to_poi_post_map, source_poi_index,
                         destination_poi_index, post_approximate_distance_UE_N1, post_approximate_path_UE_N1, post_query_time_UE_N1);
    calculate_MST_weight(post_pairwise_distance_poi_to_poi, post_MST_weight);

    std::cout << "Post terrain updating time (CG): " << post_updating_time << " ms" << std::endl;
    std::cout << "Post terrain preprocessing time (hash mapping): " << post_hash_mapping_time_UE_N1 << " ms" << std::endl;
    std::cout << "Post terrain query time: " << post_query_time_UE_N1 << " ms" << std::endl;
    std::cout << "Post terrain memory usage (CG): " << post_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain memory usage (hash mapping): " << post_hash_mapping_memory_usage_UE_N1 / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index size: " << post_index_size_UE_N1 / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index edge number: " << post_index_edge_num_UE_N1 << " edges" << std::endl;
    std::cout << "Post terrain index/MST weight: " << post_index_weight_UE_N1 / post_MST_weight << std::endl;
    std::cout << "Post terrain approximate distance: " << post_approximate_distance_UE_N1 << ", post terrain exact distance: " << post_exact_distance << ", distance error: " << post_approximate_distance_UE_N1 / post_exact_distance - 1 << std::endl;

    std::ofstream ofs2("../output/output.txt", std::ios_base::app);
    ofs2 << post_updating_time << "\t"
         << post_hash_mapping_time_UE_N1 << "\t"
         << post_query_time_UE_N1 << "\t"
         << post_memory_usage / 1e6 << "\t"
         << post_hash_mapping_memory_usage_UE_N1 / 1e6 << "\t"
         << post_index_size_UE_N1 / 1e6 << "\t"
         << post_index_edge_num_UE_N1 << "\t"
         << post_index_weight_UE_N1 / post_MST_weight << "\t"
         << post_approximate_distance_UE_N1 / post_exact_distance - 1 << "\n\n";
    ofs2.close();

    std::cout << std::endl;

    std::cout << "== UE_N2 ==" << std::endl;

    double pre_GS_time_UE_N2 = 0;
    double pre_query_time_UE_N2 = 0;
    double pre_GS_memory_usage_UE_N2 = 0;
    double pre_index_size_UE_N2 = 0;
    int pre_index_edge_num_UE_N2 = 0;
    double pre_index_weight_UE_N2 = 0;
    double pre_approximate_distance_UE_N2 = 0;
    std::vector<geodesic::SurfacePoint> pre_approximate_path_UE_N2;
    pre_approximate_path_UE_N2.clear();

    std::unordered_map<int, double> pairwise_distance_poi_to_poi_pre_greedy_spanner_map;
    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> pairwise_path_poi_to_poi_pre_greedy_spanner_map;
    greedy_spanner(pre_pairwise_distance_poi_to_poi, epsilon, pairwise_distance_poi_to_poi_pre_greedy_spanner_map,
                   pre_pairwise_path_poi_to_poi, pairwise_path_poi_to_poi_pre_greedy_spanner_map,
                   pre_index_size_UE_N2, pre_index_edge_num_UE_N2, pre_index_weight_UE_N2, pre_GS_time_UE_N2, pre_GS_memory_usage_UE_N2);
    complete_graph_query(poi_num, pairwise_distance_poi_to_poi_pre_greedy_spanner_map, pairwise_path_poi_to_poi_pre_greedy_spanner_map,
                         source_poi_index, destination_poi_index, pre_approximate_distance_UE_N2, pre_approximate_path_UE_N2, pre_query_time_UE_N2);

    // for (int i = 0; i < pre_pairwise_distance_poi_to_poi.size(); i++)
    // {
    //     for (int j = 0; j < pre_pairwise_distance_poi_to_poi[i].size(); j++)
    //     {
    //         std::cout << pre_pairwise_distance_poi_to_poi[i][j] << std::endl;
    //     }
    // }

    std::cout << "Pre terrain preprocessing time (CG): " << pre_preprocessing_time << " ms" << std::endl;
    std::cout << "Pre terrain preprocessing time (GS): " << pre_GS_time_UE_N2 << " ms" << std::endl;
    std::cout << "Pre terrain query time: " << pre_query_time_UE_N2 << " ms" << std::endl;
    std::cout << "Pre terrain memory usage (CG): " << pre_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain memory usage (GS): " << pre_GS_memory_usage_UE_N2 / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index size: " << pre_index_size_UE_N2 / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index edge number: " << pre_index_edge_num_UE_N2 << " edges" << std::endl;
    std::cout << "Pre terrain index/MST weight: " << pre_index_weight_UE_N2 / pre_MST_weight << std::endl;
    std::cout << "Pre terrain approximate distance: " << pre_approximate_distance_UE_N2 << ", pre terrain exact distance: " << pre_exact_distance << ", distance error: " << pre_approximate_distance_UE_N2 / pre_exact_distance - 1 << std::endl;
    std::cout << std::endl;

    std::ofstream ofs3("../output/output.txt", std::ios_base::app);
    ofs3 << "== UE_N2 ==\n";
    ofs3 << write_file_header << "\t"
         << pre_preprocessing_time << "\t"
         << pre_GS_time_UE_N2 << "\t"
         << pre_query_time_UE_N2 << "\t"
         << pre_memory_usage / 1e6 << "\t"
         << pre_GS_memory_usage_UE_N2 / 1e6 << "\t"
         << pre_index_size_UE_N2 / 1e6 << "\t"
         << pre_index_edge_num_UE_N2 << "\t"
         << pre_index_weight_UE_N2 / pre_MST_weight << "\t"
         << pre_approximate_distance_UE_N2 / pre_exact_distance - 1 << "\t";
    ofs3.close();

    double post_GS_time_UE_N2 = 0;
    double post_query_time_UE_N2 = 0;
    double post_GS_memory_usage_UE_N2 = 0;
    double post_index_size_UE_N2 = 0;
    int post_index_edge_num_UE_N2 = 0;
    double post_index_weight_UE_N2 = 0;
    double post_approximate_distance_UE_N2 = 0;
    std::vector<geodesic::SurfacePoint> post_approximate_path_UE_N2;
    post_approximate_path_UE_N2.clear();

    std::unordered_map<int, double> pairwise_distance_poi_to_poi_post_greedy_spanner_map;
    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> pairwise_path_poi_to_poi_post_greedy_spanner_map;
    greedy_spanner(post_pairwise_distance_poi_to_poi, epsilon, pairwise_distance_poi_to_poi_post_greedy_spanner_map,
                   post_pairwise_path_poi_to_poi, pairwise_path_poi_to_poi_post_greedy_spanner_map,
                   post_index_size_UE_N2, post_index_edge_num_UE_N2, post_index_weight_UE_N2, post_GS_time_UE_N2, post_GS_memory_usage_UE_N2);
    complete_graph_query(poi_num, pairwise_distance_poi_to_poi_post_greedy_spanner_map, pairwise_path_poi_to_poi_post_greedy_spanner_map,
                         source_poi_index, destination_poi_index, post_approximate_distance_UE_N2, post_approximate_path_UE_N2, post_query_time_UE_N2);

    std::cout << "Post terrain updating time (CG): " << post_updating_time << " ms" << std::endl;
    std::cout << "Post terrain updating time (GS): " << post_GS_time_UE_N2 << " ms" << std::endl;
    std::cout << "Post terrain query time: " << post_query_time_UE_N2 << " ms" << std::endl;
    std::cout << "Post terrain memory usage (CG): " << post_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain memory usage (GS): " << post_GS_memory_usage_UE_N2 / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index size: " << post_index_size_UE_N2 / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index edge number: " << post_index_edge_num_UE_N2 << " edges" << std::endl;
    std::cout << "Post terrain index/MST weight: " << post_index_weight_UE_N2 / post_MST_weight << std::endl;
    std::cout << "Post terrain approximate distance: " << post_approximate_distance_UE_N2 << ", post terrain exact distance: " << post_exact_distance << ", distance error: " << post_approximate_distance_UE_N2 / post_exact_distance - 1 << std::endl;

    std::ofstream ofs4("../output/output.txt", std::ios_base::app);
    ofs4 << post_updating_time << "\t"
         << post_GS_time_UE_N2 << "\t"
         << post_query_time_UE_N2 << "\t"
         << post_memory_usage / 1e6 << "\t"
         << post_GS_memory_usage_UE_N2 / 1e6 << "\t"
         << post_index_size_UE_N2 / 1e6 << "\t"
         << post_index_edge_num_UE_N2 << "\t"
         << post_index_weight_UE_N2 / post_MST_weight << "\t"
         << post_approximate_distance_UE_N2 / post_exact_distance - 1 << "\n\n";
    ofs4.close();

    std::cout << std::endl;

    std::cout << "== UE ==" << std::endl;

    double pre_HGS_time_UE = 0;
    double pre_query_time_UE = 0;
    double pre_HGS_memory_usage_UE = 0;
    double pre_index_size_UE = 0;
    int pre_index_edge_num_UE = 0;
    double pre_index_weight_UE = 0;
    double pre_approximate_distance_UE = 0;
    std::vector<geodesic::SurfacePoint> pre_approximate_path_UE;
    pre_approximate_path_UE.clear();

    std::unordered_map<int, double> pairwise_distance_poi_to_poi_pre_hierarchy_greedy_spanner_map;
    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> pairwise_path_poi_to_poi_pre_hierarchy_greedy_spanner_map;
    hierarchy_greedy_spanner(pre_pairwise_distance_poi_to_poi, epsilon, pairwise_distance_poi_to_poi_pre_hierarchy_greedy_spanner_map,
                             pre_pairwise_path_poi_to_poi, pairwise_path_poi_to_poi_pre_hierarchy_greedy_spanner_map,
                             pre_index_size_UE, pre_index_edge_num_UE, pre_index_weight_UE, pre_HGS_time_UE, pre_HGS_memory_usage_UE);
    complete_graph_query(poi_num, pairwise_distance_poi_to_poi_pre_hierarchy_greedy_spanner_map, pairwise_path_poi_to_poi_pre_hierarchy_greedy_spanner_map,
                         source_poi_index, destination_poi_index, pre_approximate_distance_UE, pre_approximate_path_UE, pre_query_time_UE);

    // for (int i = 0; i < pre_pairwise_distance_poi_to_poi.size(); i++)
    // {
    //     for (int j = 0; j < pre_pairwise_distance_poi_to_poi[i].size(); j++)
    //     {
    //         std::cout << pre_pairwise_distance_poi_to_poi[i][j] << std::endl;
    //     }
    // }

    std::cout << "Pre terrain preprocessing time (CG): " << pre_preprocessing_time << " ms" << std::endl;
    std::cout << "Pre terrain preprocessing time (HGS): " << pre_HGS_time_UE << " ms" << std::endl;
    std::cout << "Pre terrain query time: " << pre_query_time_UE << " ms" << std::endl;
    std::cout << "Pre terrain memory usage (CG): " << pre_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain memory usage (HGS): " << pre_HGS_memory_usage_UE / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index size: " << pre_index_size_UE / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index edge number: " << pre_index_edge_num_UE << " edges" << std::endl;
    std::cout << "Pre terrain index/MST weight: " << pre_index_weight_UE / pre_MST_weight << std::endl;
    std::cout << "Pre terrain approximate distance: " << pre_approximate_distance_UE << ", pre terrain exact distance: " << pre_exact_distance << ", distance error: " << pre_approximate_distance_UE / pre_exact_distance - 1 << std::endl;
    std::cout << std::endl;

    std::ofstream ofs5("../output/output.txt", std::ios_base::app);
    ofs5 << "== UE ==\n";
    ofs5 << write_file_header << "\t"
         << pre_preprocessing_time << "\t"
         << pre_HGS_time_UE << "\t"
         << pre_query_time_UE << "\t"
         << pre_memory_usage / 1e6 << "\t"
         << pre_HGS_memory_usage_UE / 1e6 << "\t"
         << pre_index_size_UE / 1e6 << "\t"
         << pre_index_edge_num_UE << "\t"
         << pre_index_weight_UE / pre_MST_weight << "\t"
         << pre_approximate_distance_UE / pre_exact_distance - 1 << "\t";
    ofs5.close();

    double post_HGS_time_UE = 0;
    double post_query_time_UE = 0;
    double post_HGS_memory_usage_UE = 0;
    double post_index_size_UE = 0;
    int post_index_edge_num_UE = 0;
    double post_index_weight_UE = 0;
    double post_approximate_distance_UE = 0;
    std::vector<geodesic::SurfacePoint> post_approximate_path_UE;
    post_approximate_path_UE.clear();

    std::unordered_map<int, double> pairwise_distance_poi_to_poi_post_hierarchy_greedy_spanner_map;
    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> pairwise_path_poi_to_poi_post_hierarchy_greedy_spanner_map;
    hierarchy_greedy_spanner(post_pairwise_distance_poi_to_poi, epsilon, pairwise_distance_poi_to_poi_post_hierarchy_greedy_spanner_map,
                             post_pairwise_path_poi_to_poi, pairwise_path_poi_to_poi_post_hierarchy_greedy_spanner_map,
                             post_index_size_UE, post_index_edge_num_UE, post_index_weight_UE, post_HGS_time_UE, post_HGS_memory_usage_UE);
    complete_graph_query(poi_num, pairwise_distance_poi_to_poi_post_hierarchy_greedy_spanner_map, pairwise_path_poi_to_poi_post_hierarchy_greedy_spanner_map,
                         source_poi_index, destination_poi_index, post_approximate_distance_UE, post_approximate_path_UE, post_query_time_UE);

    std::cout << "Post terrain updating time (CG): " << post_updating_time << " ms" << std::endl;
    std::cout << "Post terrain updating time (HGS): " << post_HGS_time_UE << " ms" << std::endl;
    std::cout << "Post terrain query time: " << post_query_time_UE << " ms" << std::endl;
    std::cout << "Post terrain memory usage (CG): " << post_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain memory usage (HGS): " << post_HGS_memory_usage_UE / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index size: " << post_index_size_UE / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index edge number: " << post_index_edge_num_UE << " edges" << std::endl;
    std::cout << "Post terrain index/MST weight: " << post_index_weight_UE / post_MST_weight << std::endl;
    std::cout << "Post terrain approximate distance: " << post_approximate_distance_UE << ", post terrain exact distance: " << post_exact_distance << ", distance error: " << post_approximate_distance_UE / post_exact_distance - 1 << std::endl;

    std::ofstream ofs6("../output/output.txt", std::ios_base::app);
    ofs6 << post_updating_time << "\t"
         << post_HGS_time_UE << "\t"
         << post_query_time_UE << "\t"
         << post_memory_usage / 1e6 << "\t"
         << post_HGS_memory_usage_UE / 1e6 << "\t"
         << post_index_size_UE / 1e6 << "\t"
         << post_index_edge_num_UE << "\t"
         << post_index_weight_UE / post_MST_weight << "\t"
         << post_approximate_distance_UE / post_exact_distance - 1 << "\n\n";
    ofs6.close();

    std::cout << std::endl;
}

void WSPD_oracle_exp(int poi_num, geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
                     double epsilon, int source_poi_index, int destination_poi_index,
                     double pre_exact_distance, int pre_MST_weight, std::string write_file_header)
{
    std::cout << "== WSPD_oracle ==" << std::endl;

    double pre_preprocessing_time = 0;
    double pre_query_time = 0;
    double pre_memory_usage = 0;
    double pre_index_size = 0;
    int pre_index_edge_num = 0;
    double pre_index_weight = 0;
    double pre_approximate_distance = 0;
    std::vector<geodesic::SurfacePoint> pre_approximate_path;
    pre_approximate_path.clear();

    pre_or_post_WSPD_oracle_preprocessing_and_query(poi_num, pre_mesh, pre_poi_list, epsilon,
                                                    source_poi_index, destination_poi_index,
                                                    pre_preprocessing_time, pre_query_time,
                                                    pre_memory_usage, pre_index_size, pre_index_edge_num, pre_index_weight,
                                                    pre_approximate_distance, pre_approximate_path);

    std::cout << "Pre terrain preprocessing time: " << pre_preprocessing_time << " ms" << std::endl;
    std::cout << "Pre terrain query time: " << pre_query_time << " ms" << std::endl;
    std::cout << "Pre terrain memory usage: " << pre_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index size: " << pre_index_size / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index edge number: " << pre_index_edge_num << " edges" << std::endl;
    std::cout << "Pre terrain index/MST weight: " << pre_index_weight / pre_MST_weight << std::endl;
    std::cout << "Pre terrain approximate distance: " << pre_approximate_distance << ", pre terrain exact distance: " << pre_exact_distance << ", distance error: " << abs(pre_approximate_distance / pre_exact_distance - 1) << std::endl;
    std::cout << std::endl;

    std::cout << "Post terrain preprocessing time: " << pre_preprocessing_time << " ms" << std::endl;
    std::cout << "Post terrain query time: " << pre_query_time << " ms" << std::endl;
    std::cout << "Post terrain memory usage: " << pre_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index size: " << pre_index_size / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index edge number: " << pre_index_edge_num << " edges" << std::endl;
    std::cout << "Post terrain index/MST weight: " << pre_index_weight / pre_MST_weight << std::endl;
    std::cout << "Post terrain approximate distance: " << pre_approximate_distance << ", post terrain exact distance: " << pre_exact_distance << ", distance error: " << abs(pre_approximate_distance / pre_exact_distance - 1) << std::endl;

    std::ofstream ofs("../output/output.txt", std::ios_base::app);
    ofs << "== WSPD_oracle ==\n";
    ofs << write_file_header << "\t"
        << pre_preprocessing_time << "\t"
        << 0 << "\t"
        << pre_query_time << "\t"
        << pre_memory_usage / 1e6 << "\t"
        << 0 << "\t"
        << pre_index_size / 1e6 << "\t"
        << pre_index_edge_num << "\t"
        << pre_index_weight / pre_MST_weight << "\t"
        << abs(pre_approximate_distance / pre_exact_distance - 1) << "\t"
        << pre_preprocessing_time << "\t"
        << 0 << "\t"
        << pre_query_time << "\t"
        << pre_memory_usage / 1e6 << "\t"
        << 0 << "\t"
        << pre_index_size / 1e6 << "\t"
        << pre_index_edge_num << "\t"
        << pre_index_weight / pre_MST_weight << "\t"
        << abs(pre_approximate_distance / pre_exact_distance - 1) << "\n\n";
    ofs.close();

    std::cout << std::endl;
}

void SP_oracle_exp(geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list, double epsilon,
                   int source_poi_index, int destination_poi_index, double pre_exact_distance,
                   std::string write_file_header)
{
    std::cout << "== SP_oracle ==" << std::endl;

    double pre_preprocessing_time = 0;
    double pre_query_time = 0;
    double pre_memory_usage = 0;
    double pre_index_size = 0;
    double pre_approximate_distance = 0;
    std::vector<geodesic::SurfacePoint> pre_approximate_path;
    pre_approximate_path.clear();

    pre_or_post_SP_oracle_preprocessing_and_query(pre_mesh, pre_poi_list,
                                                  epsilon, source_poi_index, destination_poi_index,
                                                  pre_preprocessing_time, pre_query_time,
                                                  pre_memory_usage, pre_index_size, pre_approximate_distance,
                                                  pre_approximate_path);

    std::cout << "Pre terrain preprocessing time: " << pre_preprocessing_time << " ms" << std::endl;
    std::cout << "Pre terrain query time: " << pre_query_time << " ms" << std::endl;
    std::cout << "Pre terrain memory usage: " << pre_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index size: " << pre_index_size / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain approximate distance: " << pre_approximate_distance << ", pre terrain exact distance: " << pre_exact_distance << ", distance error: " << pre_approximate_distance / pre_exact_distance - 1 << std::endl;
    std::cout << std::endl;

    std::cout << "Post terrain preprocessing time: " << pre_preprocessing_time << " ms" << std::endl;
    std::cout << "Post terrain query time: " << pre_query_time << " ms" << std::endl;
    std::cout << "Post terrain memory usage: " << pre_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index size: " << pre_index_size / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain approximate distance: " << pre_approximate_distance << ", post terrain exact distance: " << pre_exact_distance << ", distance error: " << pre_approximate_distance / pre_exact_distance - 1 << std::endl;

    std::ofstream ofs("../output/output.txt", std::ios_base::app);
    ofs << "== SP_oracle ==\n";
    ofs << write_file_header << "\t"
        << pre_preprocessing_time << "\t"
        << 0 << "\t"
        << pre_query_time << "\t"
        << pre_memory_usage / 1e6 << "\t"
        << 0 << "\t"
        << pre_index_size / 1e6 << "\t"
        << 0 << "\t"
        << 0 << "\t"
        << pre_approximate_distance / pre_exact_distance - 1 << "\t"
        << pre_preprocessing_time << "\t"
        << 0 << "\t"
        << pre_query_time << "\t"
        << pre_memory_usage / 1e6 << "\t"
        << 0 << "\t"
        << pre_index_size / 1e6 << "\t"
        << 0 << "\t"
        << 0 << "\t"
        << pre_approximate_distance / pre_exact_distance - 1 << "\n\n";
    ofs.close();

    std::cout << std::endl;
}

void KF_exp(geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
            double epsilon, int source_poi_index, int destination_poi_index,
            double pre_exact_distance, std::string write_file_header)
{
    std::cout << "== KF ==" << std::endl;

    double pre_query_time = 0;
    double pre_memory_usage = 0;
    double pre_approximate_distance = 0;
    std::vector<geodesic::SurfacePoint> pre_approximate_path;
    pre_approximate_path.clear();

    pre_or_post_KF_query(pre_mesh, pre_poi_list, epsilon, source_poi_index,
                         destination_poi_index, pre_query_time, pre_memory_usage,
                         pre_approximate_distance, pre_approximate_path);

    std::cout << "Pre terrain query time: " << pre_query_time << " ms" << std::endl;
    std::cout << "Pre terrain memory usage: " << pre_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain approximate distance: " << pre_approximate_distance << ", pre terrain exact distance: " << pre_exact_distance << ", distance error: " << pre_approximate_distance / pre_exact_distance - 1 << std::endl;
    std::cout << std::endl;

    std::cout << "Post terrain query time: " << pre_query_time << " ms" << std::endl;
    std::cout << "Post terrain memory usage: " << pre_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain approximate distance: " << pre_approximate_distance << ", post terrain exact distance: " << pre_exact_distance << ", distance error: " << pre_approximate_distance / pre_exact_distance - 1 << std::endl;

    std::ofstream ofs("../output/output.txt", std::ios_base::app);
    ofs << "== KF ==\n";
    ofs << write_file_header << "\t"
        << 0 << "\t"
        << 0 << "\t"
        << pre_query_time << "\t"
        << pre_memory_usage / 1e6 << "\t"
        << 0 << "\t"
        << 0 << "\t"
        << 0 << "\t"
        << 0 << "\t"
        << pre_approximate_distance / pre_exact_distance - 1 << "\t"
        << 0 << "\t"
        << 0 << "\t"
        << pre_query_time << "\t"
        << pre_memory_usage / 1e6 << "\t"
        << 0 << "\t"
        << 0 << "\t"
        << 0 << "\t"
        << 0 << "\t"
        << pre_approximate_distance / pre_exact_distance - 1 << "\n\n";
    ofs.close();

    std::cout << std::endl;
}

struct input_poi_struct
{
    std::string pre_poi;
    std::string post_poi;
    int source_poi_index;
    int destination_poi_index;

    input_poi_struct() {}

    input_poi_struct(std::string _pre_poi, std::string _post_poi,
                     int _source_poi_index, int _destination_poi_index)
    {
        pre_poi = _pre_poi;
        post_poi = _post_poi;
        source_poi_index = _source_poi_index;
        destination_poi_index = _destination_poi_index;
    }
};

void run_algorithms(std::string input_folder, std::string input_poi_folder,
                    std::map<int, std::pair<std::string, std::string>> input_file,
                    std::map<int, input_poi_struct> input_poi,
                    double epsilon, int input_file_index, int input_poi_index,
                    bool run_WSPD_oracle, bool run_SP_oracle, bool run_KF)
{
    std::string pre_input_file = input_folder + input_file[input_file_index].first;
    std::string post_input_file = input_folder + input_file[input_file_index].second;
    std::string pre_input_poi = input_poi_folder + input_poi[input_poi_index].pre_poi;
    std::string post_input_poi = input_poi_folder + input_poi[input_poi_index].post_poi;
    int source_poi_index = input_poi[input_poi_index].source_poi_index;
    int destination_poi_index = input_poi[input_poi_index].destination_poi_index;
    assert(source_poi_index != destination_poi_index);

    std::cout << "============" << std::endl;

    std::vector<double> pre_points;
    std::vector<unsigned> pre_faces;
    geodesic::read_mesh_from_file(&pre_input_file[0], pre_points, pre_faces);
    geodesic::Mesh pre_mesh;
    pre_mesh.initialize_mesh_data(pre_points, pre_faces);
    std::vector<int> pre_poi_list;
    int poi_num;
    std::ifstream pre_input(&pre_input_poi[0], std::ios::in);
    pre_input >> poi_num;
    pre_poi_list.resize(poi_num);
    for (int i = 0; i < poi_num; i++)
    {
        pre_input >> pre_poi_list[i];
    }
    std::vector<double> post_points;
    std::vector<unsigned> post_faces;
    geodesic::read_mesh_from_file(&post_input_file[0], post_points, post_faces);
    geodesic::Mesh post_mesh;
    post_mesh.initialize_mesh_data(post_points, post_faces);
    std::vector<int> post_poi_list;
    std::ifstream post_input(&post_input_poi[0], std::ios::in);
    post_input >> poi_num;
    post_poi_list.resize(poi_num);
    for (int i = 0; i < poi_num; i++)
    {
        post_input >> post_poi_list[i];
    }
    assert(pre_mesh.faces().size() == post_mesh.faces().size());

    std::string write_file_header = input_file[input_file_index].first + "\t" +
                                    input_file[input_file_index].second + "\t" +
                                    std::to_string(pre_mesh.faces().size()) + "\t" +
                                    std::to_string(poi_num) + "\t" +
                                    std::to_string(epsilon);

    std::cout << "pre_dataset: " << input_file[input_file_index].first << "\tpost_dataset: " << input_file[input_file_index].second << "\tdatasize: " << pre_mesh.faces().size() << "\tpoi_num: " << poi_num << "\tepsilon: " << epsilon << std::endl;
    std::cout << std::endl;

    std::ofstream ofs("../output/output.txt", std::ofstream::app);
    ofs << "# pre_dataset\tpost_dataset\tdatasize\tpoi_num\tepsilon\tpre_preprocessing_time1\tpre_preprocessing_time2\tpre_query_time\tpre_memory_usage1\tpre_memory_usage2\tpre_index_size\tpre_index_edge_num\tpre_index/MST_weight\tpre_distance_error\tpost_updating_time1\tpost_updating_time2\tpost_query_time\tpost_memory_usage1\tpost_memory_usage2\tpost_index_size\tpost_index_edge_num\tpost_index/MST_weight\tpost_distance_error\n\n";
    ofs.close();

    assert(source_poi_index <= poi_num - 1 && destination_poi_index <= poi_num - 1);

    double pre_exact_distance = 0;
    double post_exact_distance = 0;
    int pre_MST_weight = 0;
    int post_MST_weight = 0;

    exact_distance(&pre_mesh, pre_poi_list, &post_mesh, post_poi_list, source_poi_index,
                   destination_poi_index, pre_exact_distance, post_exact_distance);

    UE_N1_and_UE_N2_and_UE_exp(
        poi_num, &pre_mesh, pre_poi_list, &post_mesh, post_poi_list, epsilon, source_poi_index,
        destination_poi_index, pre_exact_distance, post_exact_distance, pre_MST_weight,
        post_MST_weight, write_file_header);

    if (run_WSPD_oracle)
    {
        WSPD_oracle_exp(poi_num, &pre_mesh, pre_poi_list, epsilon, source_poi_index, destination_poi_index,
                        pre_exact_distance, pre_MST_weight, write_file_header);
    }
    if (run_SP_oracle)
    {
        SP_oracle_exp(&pre_mesh, pre_poi_list, epsilon, source_poi_index, destination_poi_index,
                      pre_exact_distance, write_file_header);
    }
    if (run_KF)
    {
        KF_exp(&pre_mesh, pre_poi_list, epsilon, source_poi_index,
               destination_poi_index, pre_exact_distance, write_file_header);
    }
}

int main(int argc, char **argv)
{
    int option = std::stoi(argv[1]);
    std::string input_folder = "../exp_input/";
    std::string input_poi_folder = "../exp_input_poi/";
    // std::cout.precision(30);

    // run small dataset when varying epsilon
    if (option == 1)
    {
        std::map<int, std::pair<std::string, std::string>> input_file{
            {0, std::make_pair("SCpre_10082.off", "SCpost_10082.off")},
            {1, std::make_pair("AUpre_10082.off", "AUpost_10082.off")},
            {2, std::make_pair("VSpre_10082.off", "VSpost_10082.off")},
        };

        std::map<int, input_poi_struct> input_poi{
            {0, input_poi_struct("SCpre_50_poi_on_10082.txt", "SCpost_50_poi_on_10082.txt", 0, 1)},
            {1, input_poi_struct("AUpre_50_poi_on_10082.txt", "AUpost_50_poi_on_10082.txt", 0, 1)},
            {2, input_poi_struct("VSpre_50_poi_on_10082.txt", "VSpost_50_poi_on_10082.txt", 0, 1)},
        };

        std::vector<double> epsilon_list = {0.05, 0.1, 0.25, 0.5, 0.75, 1};

        for (int i = 0; i < input_file.size(); i++)
        {
            for (int j = 0; j < epsilon_list.size(); j++)
            {
                //********** change here **********
                double epsilon = epsilon_list[j];
                int input_file_index = i;
                int input_poi_index = i;
                //********** change here **********

                run_algorithms(&input_folder[0], &input_poi_folder[0], input_file, input_poi, epsilon,
                               input_file_index, input_poi_index, true, true, true);
            }
        }
    }

    // run small dataset when varying POI number
    if (option == 2)
    {
        std::map<int, std::pair<std::string, std::string>> input_file{
            {0, std::make_pair("SCpre_10082.off", "SCpost_10082.off")},
            {1, std::make_pair("AUpre_10082.off", "AUpost_10082.off")},
            {2, std::make_pair("VSpre_10082.off", "VSpost_10082.off")},
        };

        std::map<int, input_poi_struct> input_poi{
            {0, input_poi_struct("SCpre_100_poi_on_10082.txt", "SCpost_100_poi_on_10082.txt", 0, 1)},
            {1, input_poi_struct("SCpre_150_poi_on_10082.txt", "SCpost_150_poi_on_10082.txt", 0, 1)},
            {2, input_poi_struct("SCpre_200_poi_on_10082.txt", "SCpost_200_poi_on_10082.txt", 0, 1)},
            {3, input_poi_struct("SCpre_250_poi_on_10082.txt", "SCpost_250_poi_on_10082.txt", 0, 1)},
            {4, input_poi_struct("AUpre_100_poi_on_10082.txt", "AUpost_100_poi_on_10082.txt", 0, 1)},
            {5, input_poi_struct("AUpre_150_poi_on_10082.txt", "AUpost_150_poi_on_10082.txt", 0, 1)},
            {6, input_poi_struct("AUpre_200_poi_on_10082.txt", "AUpost_200_poi_on_10082.txt", 0, 1)},
            {7, input_poi_struct("AUpre_250_poi_on_10082.txt", "AUpost_250_poi_on_10082.txt", 0, 1)},
            {8, input_poi_struct("VSpre_100_poi_on_10082.txt", "VSpost_100_poi_on_10082.txt", 0, 1)},
            {9, input_poi_struct("VSpre_150_poi_on_10082.txt", "VSpost_150_poi_on_10082.txt", 0, 1)},
            {10, input_poi_struct("VSpre_200_poi_on_10082.txt", "VSpost_200_poi_on_10082.txt", 0, 1)},
            {11, input_poi_struct("VSpre_250_poi_on_10082.txt", "VSpost_250_poi_on_10082.txt", 0, 1)},
        };

        for (int i = 0; i < input_file.size(); i++)
        {
            for (int j = 0; j < input_poi.size() / input_file.size(); j++)
            {
                //********** change here **********
                double epsilon = 0.1;
                int input_file_index = i;
                int input_poi_index = i * input_poi.size() / input_file.size() + j;
                //********** change here **********

                run_algorithms(&input_folder[0], &input_poi_folder[0], input_file, input_poi, epsilon,
                               input_file_index, input_poi_index, true, false, false);
            }
        }
    }

    // run small dataset when varying dataset size
    if (option == 3)
    {
        std::map<int, std::pair<std::string, std::string>> input_file{
            {0, std::make_pair("SCpre_20000.off", "SCpost_20000.off")},
            {1, std::make_pair("SCpre_30258.off", "SCpost_30258.off")},
            {2, std::make_pair("SCpre_40328.off", "SCpost_40328.off")},
            {3, std::make_pair("SCpre_50562.off", "SCpost_50562.off")},
        };

        std::map<int, input_poi_struct> input_poi{
            {0, input_poi_struct("SCpre_50_poi_on_20000.txt", "SCpost_50_poi_on_20000.txt", 0, 1)},
            {1, input_poi_struct("SCpre_50_poi_on_30258.txt", "SCpost_50_poi_on_30258.txt", 0, 1)},
            {2, input_poi_struct("SCpre_50_poi_on_40328.txt", "SCpost_50_poi_on_40328.txt", 0, 1)},
            {3, input_poi_struct("SCpre_50_poi_on_50562.txt", "SCpost_50_poi_on_50562.txt", 0, 1)},
        };

        for (int i = 0; i < input_file.size(); i++)
        {
            //********** change here **********
            double epsilon = 0.1;
            int input_file_index = i;
            int input_poi_index = i;
            //********** change here **********

            run_algorithms(&input_folder[0], &input_poi_folder[0], input_file, input_poi, epsilon,
                           input_file_index, input_poi_index, true, true, true);
        }
    }

    // run large dataset when varying dataset size
    if (option == 4)
    {
        std::map<int, std::pair<std::string, std::string>> input_file{
            {0, std::make_pair("SCpre_500000.off", "SCpost_500000.off")},
            {1, std::make_pair("AUpre_500000.off", "AUpost_500000.off")},
            {2, std::make_pair("VSpre_500000.off", "VSpost_500000.off")},
        };

        std::map<int, input_poi_struct> input_poi{
            {0, input_poi_struct("SCpre_500_poi_on_500000.txt", "SCpost_500_poi_on_500000.txt", 0, 1)},
            {1, input_poi_struct("AUpre_500_poi_on_500000.txt", "AUpost_500_poi_on_500000.txt", 0, 1)},
            {2, input_poi_struct("VSpre_500_poi_on_500000.txt", "VSpost_500_poi_on_500000.txt", 0, 1)},
        };

        std::vector<double> epsilon_list = {0.05, 0.1, 0.25, 0.5, 0.75, 1};

        for (int i = 0; i < input_file.size(); i++)
        {
            for (int j = 0; j < epsilon_list.size(); j++)
            {
                //********** change here **********
                double epsilon = epsilon_list[j];
                int input_file_index = i;
                int input_poi_index = i;
                //********** change here **********

                run_algorithms(&input_folder[0], &input_poi_folder[0], input_file, input_poi, epsilon,
                               input_file_index, input_poi_index, false, false, true);
            }
        }
    }

    // run large dataset when varying POI number
    if (option == 5)
    {
        std::map<int, std::pair<std::string, std::string>> input_file{
            {0, std::make_pair("SCpre_500000.off", "SCpost_500000.off")},
            {1, std::make_pair("AUpre_500000.off", "AUpost_500000.off")},
            {2, std::make_pair("VSpre_500000.off", "VSpost_500000.off")},
        };

        std::map<int, input_poi_struct> input_poi{
            {0, input_poi_struct("SCpre_1000_poi_on_500000.txt", "SCpost_1000_poi_on_500000.txt", 0, 1)},
            {1, input_poi_struct("SCpre_1500_poi_on_500000.txt", "SCpost_1500_poi_on_500000.txt", 0, 1)},
            {2, input_poi_struct("SCpre_2000_poi_on_500000.txt", "SCpost_2000_poi_on_500000.txt", 0, 1)},
            {3, input_poi_struct("SCpre_2500_poi_on_500000.txt", "SCpost_2500_poi_on_500000.txt", 0, 1)},
            {4, input_poi_struct("AUpre_1000_poi_on_500000.txt", "AUpost_1000_poi_on_500000.txt", 0, 1)},
            {5, input_poi_struct("AUpre_1500_poi_on_500000.txt", "AUpost_1500_poi_on_500000.txt", 0, 1)},
            {6, input_poi_struct("AUpre_2000_poi_on_500000.txt", "AUpost_2000_poi_on_500000.txt", 0, 1)},
            {7, input_poi_struct("AUpre_2500_poi_on_500000.txt", "AUpost_2500_poi_on_500000.txt", 0, 1)},
            {8, input_poi_struct("VSpre_1000_poi_on_500000.txt", "VSpost_1000_poi_on_500000.txt", 0, 1)},
            {9, input_poi_struct("VSpre_1500_poi_on_500000.txt", "VSpost_1500_poi_on_500000.txt", 0, 1)},
            {10, input_poi_struct("VSpre_2000_poi_on_500000.txt", "VSpost_2000_poi_on_500000.txt", 0, 1)},
            {11, input_poi_struct("VSpre_2500_poi_on_500000.txt", "VSpost_2500_poi_on_500000.txt", 0, 1)},
        };

        for (int i = 0; i < input_file.size(); i++)
        {
            for (int j = 0; j < input_poi.size() / input_file.size(); j++)
            {
                //********** change here **********
                double epsilon = 0.25;
                int input_file_index = i;
                int input_poi_index = i * input_poi.size() / input_file.size() + j;
                //********** change here **********

                run_algorithms(&input_folder[0], &input_poi_folder[0], input_file, input_poi, epsilon,
                               input_file_index, input_poi_index, false, false, false);
            }
        }
    }

    // run large dataset when varying dataset size
    if (option == 6)
    {
        std::map<int, std::pair<std::string, std::string>> input_file{
            {0, std::make_pair("SCpre_1002528.off", "SCpost_1002528.off")},
            {1, std::make_pair("SCpre_1503378.off", "SCpost_1503378.off")},
            {2, std::make_pair("SCpre_2000000.off", "SCpost_2000000.off")},
            {3, std::make_pair("SCpre_2504322.off", "SCpost_2504322.off")},
            {4, std::make_pair("AUpre_1002528.off", "AUpost_1002528.off")},
            {5, std::make_pair("AUpre_1503378.off", "AUpost_1503378.off")},
            {6, std::make_pair("AUpre_2000000.off", "AUpost_2000000.off")},
            {7, std::make_pair("AUpre_2504322.off", "AUpost_2504322.off")},
            {8, std::make_pair("VSpre_1002528.off", "VSpost_1002528.off")},
            {9, std::make_pair("VSpre_1503378.off", "VSpost_1503378.off")},
            {10, std::make_pair("VSpre_2000000.off", "VSpost_2000000.off")},
            {11, std::make_pair("VSpre_2504322.off", "VSpost_2504322.off")},
        };

        std::map<int, input_poi_struct> input_poi{
            {0, input_poi_struct("SCpre_500_poi_on_1002528.txt", "SCpost_500_poi_on_1002528.txt", 0, 1)},
            {1, input_poi_struct("SCpre_500_poi_on_1503378.txt", "SCpost_500_poi_on_1503378.txt", 0, 1)},
            {2, input_poi_struct("SCpre_500_poi_on_2000000.txt", "SCpost_500_poi_on_2000000.txt", 0, 1)},
            {3, input_poi_struct("SCpre_500_poi_on_2504322.txt", "SCpost_500_poi_on_2504322.txt", 0, 1)},
            {4, input_poi_struct("AUpre_500_poi_on_1002528.txt", "AUpost_500_poi_on_1002528.txt", 0, 1)},
            {5, input_poi_struct("AUpre_500_poi_on_1503378.txt", "AUpost_500_poi_on_1503378.txt", 0, 1)},
            {6, input_poi_struct("AUpre_500_poi_on_2000000.txt", "AUpost_500_poi_on_2000000.txt", 0, 1)},
            {7, input_poi_struct("AUpre_500_poi_on_2504322.txt", "AUpost_500_poi_on_2504322.txt", 0, 1)},
            {8, input_poi_struct("VSpre_500_poi_on_1002528.txt", "VSpost_500_poi_on_1002528.txt", 0, 1)},
            {9, input_poi_struct("VSpre_500_poi_on_1503378.txt", "VSpost_500_poi_on_1503378.txt", 0, 1)},
            {10, input_poi_struct("VSpre_500_poi_on_2000000.txt", "VSpost_500_poi_on_2000000.txt", 0, 1)},
            {11, input_poi_struct("VSpre_500_poi_on_2504322.txt", "VSpost_500_poi_on_2504322.txt", 0, 1)},
        };

        for (int i = 0; i < input_file.size(); i++)
        {
            //********** change here **********
            double epsilon = 0.25;
            int input_file_index = i;
            int input_poi_index = i;
            //********** change here **********

            run_algorithms(&input_folder[0], &input_poi_folder[0], input_file, input_poi, epsilon,
                           input_file_index, input_poi_index, false, false, true);
        }
    }
}
