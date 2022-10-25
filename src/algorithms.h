#include "helper.h"
#include <sstream>
#include <unistd.h>
#include <sys/resource.h>
#include <sys/times.h>
#include <iostream>
#include <sstream>
#include <chrono>
#include <unordered_map>
#include <set>

void pre_or_post_exact_distance(geodesic::Mesh *mesh, std::vector<int> &poi_list,
                                int source_poi_index, int destination_poi_index,
                                double &exact_distance)
{
    geodesic::GeodesicAlgorithmExact algorithm(mesh);
    double const distance_limit = geodesic::GEODESIC_INF;
    geodesic::SurfacePoint source(&mesh->vertices()[poi_list[source_poi_index]]);
    geodesic::SurfacePoint destination(&mesh->vertices()[poi_list[destination_poi_index]]);
    std::vector<geodesic::SurfacePoint> one_source_poi_list(1, source);
    std::vector<geodesic::SurfacePoint> one_destination_poi_list(1, destination);
    algorithm.propagate(one_source_poi_list, distance_limit, &one_destination_poi_list);
    algorithm.best_source(destination, exact_distance);
    exact_distance = round(exact_distance * 1000000000.0) / 1000000000.0;
    // std::cout << "Exact distance: " << exact_distance << std::endl;
}

void calculate_MST_weight(std::vector<std::vector<double>> pairwise_distance_poi_to_poi, int &MST_weight)
{
    int V = pairwise_distance_poi_to_poi.size();
    int E = V * (V - 1) / 2;
    Graph MST_graph(V, E);
    for (int i = 0; i < pairwise_distance_poi_to_poi.size(); i++)
    {
        for (int j = 0; j < pairwise_distance_poi_to_poi[i].size(); j++)
        {
            if (j == 0)
            {
                continue;
            }
            MST_graph.add_edge_Kruskal(i, i + j, pairwise_distance_poi_to_poi[i][j]);
        }
    }
    MST_weight = MST_graph.MST_Kruskal();
    // std::cout << "MST weight: " << MST_weight << std::endl;
}

void pre_complete_graph_preprocessing(int poi_num, geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
                                      std::vector<std::vector<std::vector<int>>> &pre_face_sequence_index_list,
                                      std::vector<std::vector<double>> &pairwise_distance_poi_to_poi,
                                      std::vector<std::vector<bool>> &pairwise_distance_poi_to_poi_changed,
                                      std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> &pairwise_path_poi_to_poi,
                                      std::vector<std::vector<double>> &pairwise_distance_poi_to_vertex,
                                      double &preprocessing_time, double &memory_usage)
{
    auto start_preprocessing_time = std::chrono::high_resolution_clock::now();

    geodesic::GeodesicAlgorithmExact pre_algorithm(pre_mesh);

    // std::vector<std::vector<std::vector<int>>> pre_face_sequence_index_list;
    int pre_face_sequence_index_list_size = 0;
    std::vector<std::vector<int>> one_poi_to_other_poi_pre_face_sequence_index_list;
    std::vector<int> one_poi_pre_face_sequence_index_list;
    pre_face_sequence_index_list.clear();

    double const distance_limit = geodesic::GEODESIC_INF;
    // std::vector<std::vector<double>> pairwise_distance_poi_to_poi(poi_num, std::vector<double>(poi_num, 0));
    // std::vector<std::vector<double>> pairwise_distance_poi_to_poi;
    pairwise_distance_poi_to_poi.clear();
    // std::vector<std::vector<bool>> pairwise_distance_poi_to_poi_changed(poi_num, std::vector<bool>(poi_num, false));
    // std::vector<std::vector<bool>> pairwise_distance_poi_to_poi_changed;
    pairwise_distance_poi_to_poi_changed.clear();
    pairwise_path_poi_to_poi.clear();
    int pairwise_path_poi_to_poi_size = 0;
    // std::vector<std::vector<double>> pairwise_distance_poi_to_vertex(poi_num, std::vector<double>(pre_mesh->vertices().size(), 0));
    std::vector<geodesic::SurfacePoint> one_source_poi_list;
    std::vector<geodesic::SurfacePoint> destinations_poi_list;
    std::vector<geodesic::SurfacePoint> destinations_vertex_list;

    // calculate the pairwise geodesic distance on pre terrain
    for (int i = 0; i < poi_num; i++)
    {
        std::vector<double> current_poi_to_other_poi_distance;
        current_poi_to_other_poi_distance.clear();
        std::vector<bool> current_poi_to_other_poi_distance_changed;
        current_poi_to_other_poi_distance_changed.clear();
        std::vector<std::vector<geodesic::SurfacePoint>> current_poi_to_other_poi_path;
        current_poi_to_other_poi_path.clear();
        one_source_poi_list.clear();
        destinations_poi_list.clear();
        destinations_vertex_list.clear();
        one_poi_to_other_poi_pre_face_sequence_index_list.clear();
        one_poi_pre_face_sequence_index_list.clear();
        one_source_poi_list.push_back(geodesic::SurfacePoint(&pre_mesh->vertices()[pre_poi_list[i]]));

        for (int j = i; j < poi_num; j++)
        {
            destinations_poi_list.push_back(geodesic::SurfacePoint(&pre_mesh->vertices()[pre_poi_list[j]]));
            // std::cout << poi_list[j] << std::endl;
        }

        for (int j = 0; j < pre_mesh->vertices().size(); j++)
        {
            destinations_vertex_list.push_back(geodesic::SurfacePoint(&pre_mesh->vertices()[j]));
        }

        // auto st3 = std::chrono::high_resolution_clock::now();
        pre_algorithm.propagate(one_source_poi_list, distance_limit, &destinations_vertex_list);
        // auto st4 = std::chrono::high_resolution_clock::now();
        // auto d3 = std::chrono::duration_cast<std::chrono::milliseconds>(st4 - st3);
        // std::cout << "before update propagate: " << d3.count() << " ms" << std::endl;

        for (int j = 0; j < i; j++)
        {
            one_poi_to_other_poi_pre_face_sequence_index_list.push_back(one_poi_pre_face_sequence_index_list);
            // one_poi_to_other_poi_pre_face_sequence_index_list.push_back(pre_face_sequence_index_list[j][i]);
        }

        for (int j = i; j < poi_num; j++)
        {
            std::vector<geodesic::SurfacePoint> path;
            pre_algorithm.trace_back(destinations_poi_list[j - i], path);
            current_poi_to_other_poi_distance.push_back(length(path));
            current_poi_to_other_poi_distance_changed.push_back(false);
            current_poi_to_other_poi_path.push_back(path);
            // pairwise_distance_poi_to_poi[i][j] = length(path);
            // std::cout << path.size() << std::endl;
            get_face_sequence(pre_mesh, path, one_poi_pre_face_sequence_index_list);
            // for (int k = 0; k < one_poi_pre_face_sequence_index_list.size(); k++)
            // {
            //     std::cout << one_poi_pre_face_sequence_index_list[k] << " ";
            // }
            // std::cout << std::endl;
            one_poi_to_other_poi_pre_face_sequence_index_list.push_back(one_poi_pre_face_sequence_index_list);
            pre_face_sequence_index_list_size += one_poi_pre_face_sequence_index_list.size();
            pairwise_path_poi_to_poi_size += path.size();

            // std::cout << "i " << i << ", j " << j << " - " << distance << std::endl;
        }
        pairwise_distance_poi_to_poi.push_back(current_poi_to_other_poi_distance);
        pairwise_distance_poi_to_poi_changed.push_back(current_poi_to_other_poi_distance_changed);
        pairwise_path_poi_to_poi.push_back(current_poi_to_other_poi_path);
        pre_face_sequence_index_list.push_back(one_poi_to_other_poi_pre_face_sequence_index_list);
        // for (int j = 0; j < poi_num; j++)
        // {
        //     // std::cout << "i " << i << ", j " << j << ": ";
        //     for (int k = 0; k < pre_face_sequence_index_list[i][j].size(); k++)
        //     {
        //         std::cout << pre_face_sequence_index_list[i][j][k] << " ";
        //     }
        //     std::cout << std::endl;
        // }

        // for (int j = 0; j < pairwise_distance_poi_to_poi[i].size(); j++)
        // {
        //     std::cout << "i " << i << ", j " << j << " - " << pairwise_distance_poi_to_poi[i][j] << std::endl;
        // }
        for (int j = 0; j < pre_mesh->vertices().size(); j++)
        {
            double distance;
            pre_algorithm.best_source(destinations_vertex_list[j], distance);
            pairwise_distance_poi_to_vertex[i][j] = distance;
            // std::cout << "i " << i << ", j " << j << " - " << distance << std::endl;
        }
    }
    // std::cout << "pre_algorithm.get_memory(): " << pre_algorithm.get_memory() << " , pre_face_sequence_index_list_size * sizeof(int): " << pre_face_sequence_index_list_size * sizeof(int) << std::endl;
    memory_usage += pre_algorithm.get_memory() + pre_face_sequence_index_list_size * sizeof(int); // + 0.5 * poi_num * (poi_num + 1) * sizeof(double) + pairwise_path_poi_to_poi_size * sizeof(geodesic::SurfacePoint);

    auto stop_preprocessing_time = std::chrono::high_resolution_clock::now();
    auto duration_preprocessing_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_preprocessing_time - start_preprocessing_time);
    preprocessing_time = duration_preprocessing_time.count();
}

void post_complete_graph_updating(int poi_num, geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
                                  geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list,
                                  std::vector<std::vector<std::vector<int>>> &pre_face_sequence_index_list,
                                  std::vector<std::vector<double>> &pairwise_distance_poi_to_poi,
                                  std::vector<std::vector<bool>> &pairwise_distance_poi_to_poi_changed,
                                  std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> &pairwise_path_poi_to_poi,
                                  std::vector<std::vector<double>> &pairwise_distance_poi_to_vertex,
                                  double &updating_time, double &memory_usage)
{
    auto start_updating_time = std::chrono::high_resolution_clock::now();

    // auto st1 = std::chrono::high_resolution_clock::now();

    geodesic::GeodesicAlgorithmExact post_algorithm(post_mesh);

    std::vector<geodesic::SurfacePoint> one_source_poi_list;
    std::vector<geodesic::SurfacePoint> destinations_poi_list;
    double const distance_limit = geodesic::GEODESIC_INF;

    int changed_pairwise_path_poi_to_poi_size = 0;
    int changed_pairwise_distance_poi_to_poi_size = 0;

    // compare the pre and post terrain to detect changed face
    std::vector<int> changed_face_index_list;
    changed_face_index_list.clear();
    std::unordered_map<int, int> changed_face_index_unordered_map;
    changed_face_index_unordered_map.clear();
    assert(pre_mesh->faces().size() == post_mesh->faces().size());
    for (int i = 0; i < pre_mesh->faces().size(); i++)
    {
        if (pre_mesh->faces()[i].adjacent_vertices()[0]->x() != post_mesh->faces()[i].adjacent_vertices()[0]->x() ||
            pre_mesh->faces()[i].adjacent_vertices()[0]->y() != post_mesh->faces()[i].adjacent_vertices()[0]->y() ||
            pre_mesh->faces()[i].adjacent_vertices()[0]->z() != post_mesh->faces()[i].adjacent_vertices()[0]->z() ||
            pre_mesh->faces()[i].adjacent_vertices()[1]->x() != post_mesh->faces()[i].adjacent_vertices()[1]->x() ||
            pre_mesh->faces()[i].adjacent_vertices()[1]->y() != post_mesh->faces()[i].adjacent_vertices()[1]->y() ||
            pre_mesh->faces()[i].adjacent_vertices()[1]->z() != post_mesh->faces()[i].adjacent_vertices()[1]->z() ||
            pre_mesh->faces()[i].adjacent_vertices()[2]->x() != post_mesh->faces()[i].adjacent_vertices()[2]->x() ||
            pre_mesh->faces()[i].adjacent_vertices()[2]->y() != post_mesh->faces()[i].adjacent_vertices()[2]->y() ||
            pre_mesh->faces()[i].adjacent_vertices()[2]->z() != post_mesh->faces()[i].adjacent_vertices()[2]->z())
        {
            // int face_difference_error = 100;
            // for (int i = 0; i < pre_mesh->faces().size(); i++)
            // {
            //     if (abs(pre_mesh->faces()[i].adjacent_vertices()[0]->x() - post_mesh->faces()[i].adjacent_vertices()[0]->x()) > face_difference_error ||
            //         abs(pre_mesh->faces()[i].adjacent_vertices()[0]->y() - post_mesh->faces()[i].adjacent_vertices()[0]->y()) > face_difference_error ||
            //         abs(pre_mesh->faces()[i].adjacent_vertices()[0]->z() - post_mesh->faces()[i].adjacent_vertices()[0]->z()) > face_difference_error ||
            //         abs(pre_mesh->faces()[i].adjacent_vertices()[1]->x() - post_mesh->faces()[i].adjacent_vertices()[1]->x()) > face_difference_error ||
            //         abs(pre_mesh->faces()[i].adjacent_vertices()[1]->y() - post_mesh->faces()[i].adjacent_vertices()[1]->y()) > face_difference_error ||
            //         abs(pre_mesh->faces()[i].adjacent_vertices()[1]->z() - post_mesh->faces()[i].adjacent_vertices()[1]->z()) > face_difference_error ||
            //         abs(pre_mesh->faces()[i].adjacent_vertices()[2]->x() - post_mesh->faces()[i].adjacent_vertices()[2]->x()) > face_difference_error ||
            //         abs(pre_mesh->faces()[i].adjacent_vertices()[2]->y() - post_mesh->faces()[i].adjacent_vertices()[2]->y()) > face_difference_error ||
            //         abs(pre_mesh->faces()[i].adjacent_vertices()[2]->z() - post_mesh->faces()[i].adjacent_vertices()[2]->z()) > face_difference_error)
            //     {
            changed_face_index_list.push_back(i);
            changed_face_index_unordered_map[i] = i;
            // std::cout << i << std::endl;
        }
    }
    assert(changed_face_index_list.size() == changed_face_index_unordered_map.size());
    // std::cout << "changed_face_index_list.size(): " << changed_face_index_list.size() << std::endl;

    // compare the pre and post terrain to detect changed vertex
    std::vector<int> changed_vertex_index_list;
    changed_vertex_index_list.clear();
    assert(pre_mesh->vertices().size() == post_mesh->vertices().size());
    for (int i = 0; i < pre_mesh->vertices().size(); i++)
    {
        if (pre_mesh->vertices()[i].x() != post_mesh->vertices()[i].x() ||
            pre_mesh->vertices()[i].y() != post_mesh->vertices()[i].y() ||
            pre_mesh->vertices()[i].z() != post_mesh->vertices()[i].z())
        {
            changed_vertex_index_list.push_back(i);
        }
    }

    // stores the changed poi index such that (1) the index of poi changed, (2) this poi is in the changed area, (3) the path that connects this poi passes the changed area
    std::vector<int> changed_poi_index_list(poi_num, 0);

    // stores only the poi in the changed area (used in calculation the 2D distance between these pois and other pois)
    std::vector<int> poi_in_the_changed_area_index_list;
    poi_in_the_changed_area_index_list.clear();

    // for the pre and post terrain, the index of poi may change, but the actual vertex on the terrain that these two pois stand for are very close
    // if the poi for the pre and post list changed, then the value in this list becomes 1, which indicates changes
    for (int i = 0; i < poi_num; i++)
    {
        if (pre_poi_list[i] != post_poi_list[i])
        {
            // std::cout << pre_poi_list[i] << " " << post_poi_list[i] << std::endl;
            changed_poi_index_list[i] = 1;
        }
    }

    // auto st2 = std::chrono::high_resolution_clock::now();
    // auto d = std::chrono::duration_cast<std::chrono::milliseconds>(st2 - st1);
    // std::cout << "after terrain preprocessing 1: " << d.count() << " ms" << std::endl;

    // st1 = std::chrono::high_resolution_clock::now();

    // for the pre and post terrain, even if the index of poi not changes, if the poi is on the changed area, we also indicate it by 3 in the following list
    for (int i = 0; i < poi_num; i++)
    {
        if (changed_poi_index_list[i] == 1)
        {
            continue;
        }
        for (int j = 0; j < changed_face_index_list.size(); j++)
        {
            if (pre_poi_list[i] == post_poi_list[i] &&
                (pre_mesh->faces()[changed_face_index_list[j]].adjacent_vertices()[0]->id() == pre_poi_list[i] ||
                 pre_mesh->faces()[changed_face_index_list[j]].adjacent_vertices()[1]->id() == pre_poi_list[i] ||
                 pre_mesh->faces()[changed_face_index_list[j]].adjacent_vertices()[2]->id() == pre_poi_list[i]))
            {
                // std::cout << pre_poi_list[i] << std::endl;
                changed_poi_index_list[i] = 3;
                poi_in_the_changed_area_index_list.push_back(i);
                break;
            }
        }
    }

    // st2 = std::chrono::high_resolution_clock::now();
    // d = std::chrono::duration_cast<std::chrono::milliseconds>(st2 - st1);
    // std::cout << "after terrain preprocessing 2: " << d.count() << " ms" << std::endl;

    // st1 = std::chrono::high_resolution_clock::now();

    // // if the path passes the changed area ((1) the two endpoints, i.e., the two pois, in the changed area, or (2) the center part of
    // // the path passes the changed area), we also indicate the poi by 4
    // if two pois are not in the changed area, but the path passes the changed area, we also indicate the poi by 4
    // note that for a poi, if one of the path that connects this poi passes the changed area, then we need to run SSAD
    // for this poi again involves all the path connects to this poi
    for (int i = 0; i < poi_num; i++)
    {
        bool break_loop = false;
        if (changed_poi_index_list[i] == 1 || changed_poi_index_list[i] == 3)
        {
            continue;
        }
        for (int j = i; j < poi_num; j++)
        {
            if (changed_poi_index_list[j] == 1 || changed_poi_index_list[j] == 3)
            {
                continue;
            }
            for (int k = 0; k < pre_face_sequence_index_list[i][j].size(); k++)
            {
                if (changed_face_index_unordered_map.count(pre_face_sequence_index_list[i][j][k]) != 0)
                {
                    changed_poi_index_list[i] = 4;
                    break_loop = true;
                    break;
                }
                // for (int m = 0; m < changed_face_index_list.size(); m++)
                // {
                //     // std::cout << changed_face_index_list[m] << " " << pre_face_sequence_index_list[i][j][k] << std::endl;

                //     if (changed_face_index_list[m] == pre_face_sequence_index_list[i][j][k])
                //     {
                //         changed_poi_index_list[i] = 4;
                //         break_loop = true;
                //         break;
                //     }
                // }
                // if (break_loop)
                // {
                //     break;
                // }
            }
            if (break_loop)
            {
                break;
            }
        }
    }

    // st2 = std::chrono::high_resolution_clock::now();
    // d = std::chrono::duration_cast<std::chrono::milliseconds>(st2 - st1);
    // std::cout << "after terrain preprocessing 3: " << d.count() << " ms" << std::endl;

    // st1 = std::chrono::high_resolution_clock::now();

    // update the pairwise geodesic distance on post terrain for changed poi
    for (int i = 0; i < poi_num; i++)
    {
        one_source_poi_list.clear();
        destinations_poi_list.clear();

        if (changed_poi_index_list[i] == 0)
        {
            continue;
        }

        one_source_poi_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_poi_list[i]]));

        for (int j = 0; j < poi_num; j++)
        {
            if ((changed_poi_index_list[j] == 1 || changed_poi_index_list[j] == 3 || changed_poi_index_list[j] == 4) && j < i)
            {
                continue;
            }
            destinations_poi_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_poi_list[j]]));
            // std::cout << "i " << i << ", j " << j << ", j's x " << post_mesh->vertices()[post_poi_list[j]].x() << ", j's y " << post_mesh->vertices()[post_poi_list[j]].y() << ", j's z " << post_mesh->vertices()[post_poi_list[j]].z() << std::endl;
        }

        // auto st3 = std::chrono::high_resolution_clock::now();
        post_algorithm.propagate(one_source_poi_list, distance_limit, &destinations_poi_list);
        // auto st4 = std::chrono::high_resolution_clock::now();
        // auto d3 = std::chrono::duration_cast<std::chrono::milliseconds>(st4 - st3);
        // std::cout << "first propagate: " << d3.count() << " ms" << std::endl;

        int index = 0;
        for (int j = 0; j < poi_num; j++)
        {
            if ((changed_poi_index_list[j] == 1 || changed_poi_index_list[j] == 3 || changed_poi_index_list[j] == 4) && j < i)
            {
                index++;
                continue;
            }
            // double distance;
            // post_algorithm.best_source(destinations_poi_list[j - index], distance);
            std::vector<geodesic::SurfacePoint> path;
            post_algorithm.trace_back(destinations_poi_list[j - index], path);
            changed_pairwise_path_poi_to_poi_size += path.size();
            changed_pairwise_distance_poi_to_poi_size++;
            // std::cout << "i " << i << ", j " << j << ", j's x " << destinations_poi_list[j - index].x() << ", j's y " << destinations_poi_list[j - index].y() << ", j's z " << destinations_poi_list[j - index].z() << std::endl;

            if (i <= j)
            {
                pairwise_distance_poi_to_poi[i][j - i] = length(path);
                pairwise_distance_poi_to_poi_changed[i][j - i] = true;
                pairwise_path_poi_to_poi[i][j - i] = path;
            }
            else
            {
                pairwise_distance_poi_to_poi[j][i - j] = length(path);
                pairwise_distance_poi_to_poi_changed[j][i - j] = true;
                std::reverse(path.begin(), path.end());
                pairwise_path_poi_to_poi[j][i - j] = path;
            }
        }

        // for (int j = 0; j < pairwise_distance_poi_to_poi[i].size(); j++)
        // {
        //     std::cout << "i " << i << ", j " << j << " - " << pairwise_distance_poi_to_poi[i][j] << std::endl;
        // }
    }

    // st2 = std::chrono::high_resolution_clock::now();
    // d = std::chrono::duration_cast<std::chrono::milliseconds>(st2 - st1);
    // std::cout << d.count() << " ms" << std::endl;
    // for (int i = 0; i < poi_num; i++)
    // {
    //     for (int j = 0; j < pairwise_distance_poi_to_poi[i].size(); j++)
    //     {
    //         std::cout << "i " << i << ", j " << j << " - " << pairwise_distance_poi_to_poi[i][j] << std::endl;
    //     }
    // }

    // if two pois are not in the changed area, and the path doesn't pass the changed area, but one of pois is close to the
    // changed area, so we may need to update the new path with these two pois as endpoints on the new terrain, the following
    // is to check whether the original path is too close to the changed area or not, if so, we directly update the path
    // we also indicate this type of poi in changed_poi_index_list, but indicate it as 2 for clarification

    // st1 = std::chrono::high_resolution_clock::now();

    std::vector<double> euclidean_distance_of_poi_to_changed_area(poi_num, 0);
    std::vector<std::pair<double, int>> euclidean_distance_of_poi_to_changed_area_and_original_index;

    for (int i = 0; i < poi_num; i++)
    {
        if (changed_poi_index_list[i] == 1 || changed_poi_index_list[i] == 3 || changed_poi_index_list[i] == 4)
        {
            continue;
        }
        if (poi_in_the_changed_area_index_list.size() > 0)
        {
            for (int j = 0; j < poi_in_the_changed_area_index_list.size(); j++)
            {
                euclidean_distance_of_poi_to_changed_area[i] += euclidean_distance(post_mesh->vertices()[post_poi_list[i]].x(), post_mesh->vertices()[post_poi_list[i]].y(),
                                                                                   post_mesh->vertices()[poi_in_the_changed_area_index_list[j]].x(), post_mesh->vertices()[poi_in_the_changed_area_index_list[j]].y());
            }
        }
        else
        {
            euclidean_distance_of_poi_to_changed_area[i] = euclidean_distance(post_mesh->vertices()[post_poi_list[i]].x(), post_mesh->vertices()[post_poi_list[i]].y(),
                                                                              post_mesh->vertices()[changed_vertex_index_list[0]].x(), post_mesh->vertices()[changed_vertex_index_list[0]].y());
        }
        // std::cout << "euclidean_distance_of_poi_to_changed_area: " << euclidean_distance_of_poi_to_changed_area[i] << std::endl;
    }

    // auto st3 = std::chrono::high_resolution_clock::now();

    sort_min_to_max_and_get_original_index(euclidean_distance_of_poi_to_changed_area, euclidean_distance_of_poi_to_changed_area_and_original_index);

    // auto st4 = std::chrono::high_resolution_clock::now();
    // auto d3 = std::chrono::duration_cast<std::chrono::milliseconds>(st4 - st3);
    // std::cout << "sorting: " << d3.count() << " ms" << std::endl;

    assert(euclidean_distance_of_poi_to_changed_area.size() == euclidean_distance_of_poi_to_changed_area_and_original_index.size());

    for (int i = 0; i < euclidean_distance_of_poi_to_changed_area_and_original_index.size(); i++)
    {
        if (euclidean_distance_of_poi_to_changed_area_and_original_index[i].first == 0)
        {
            continue;
        }
        int current_poi_index = euclidean_distance_of_poi_to_changed_area_and_original_index[i].second;
        // std::cout << "euclidean_distance_of_poi_to_changed_area_and_original_index: " << euclidean_distance_of_poi_to_changed_area_and_original_index[i].first << ", index: " << current_poi_index << std::endl;

        assert(pairwise_distance_poi_to_poi[current_poi_index].size() == pairwise_distance_poi_to_poi_changed[current_poi_index].size());

        double max_distance = 0;
        for (int j = 0; j < euclidean_distance_of_poi_to_changed_area_and_original_index.size(); j++)
        {
            if (euclidean_distance_of_poi_to_changed_area_and_original_index[j].first == 0)
            {
                continue;
            }
            int checking_poi_index = euclidean_distance_of_poi_to_changed_area_and_original_index[j].second;

            if (current_poi_index <= checking_poi_index)
            {
                if (pairwise_distance_poi_to_poi_changed[current_poi_index][checking_poi_index - current_poi_index])
                {
                    continue;
                }
                max_distance = std::max(max_distance, pairwise_distance_poi_to_poi[current_poi_index][checking_poi_index - current_poi_index]);
            }
            else
            {
                if (pairwise_distance_poi_to_poi_changed[checking_poi_index][current_poi_index - checking_poi_index])
                {
                    continue;
                }
                max_distance = std::max(max_distance, pairwise_distance_poi_to_poi[checking_poi_index][current_poi_index - checking_poi_index]);
            }
        }
        // std::cout << max_distance << std::endl;

        // if the current poi is too close to the changed area, we need to run SSAD for this poi to update it path on the new terrain
        for (int k = 0; k < changed_vertex_index_list.size(); k++)
        {
            if (pairwise_distance_poi_to_vertex[current_poi_index][changed_vertex_index_list[k]] < max_distance / 2)
            {
                changed_poi_index_list[current_poi_index] = 2;

                std::vector<int> destinations_poi_index_list;
                destinations_poi_index_list.clear();
                one_source_poi_list.clear();
                destinations_poi_list.clear();
                one_source_poi_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_poi_list[current_poi_index]]));

                for (int j = 0; j < euclidean_distance_of_poi_to_changed_area_and_original_index.size(); j++)
                {
                    int the_other_poi_index = euclidean_distance_of_poi_to_changed_area_and_original_index[j].second;

                    if ((current_poi_index <= the_other_poi_index && !pairwise_distance_poi_to_poi_changed[current_poi_index][the_other_poi_index - current_poi_index]) ||
                        (current_poi_index > the_other_poi_index && !pairwise_distance_poi_to_poi_changed[the_other_poi_index][current_poi_index - the_other_poi_index]))
                    {
                        destinations_poi_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_poi_list[the_other_poi_index]]));
                        destinations_poi_index_list.push_back(the_other_poi_index);
                    }
                }

                // auto st3 = std::chrono::high_resolution_clock::now();
                post_algorithm.propagate(one_source_poi_list, distance_limit, &destinations_poi_list);
                // auto st4 = std::chrono::high_resolution_clock::now();
                // auto d3 = std::chrono::duration_cast<std::chrono::milliseconds>(st4 - st3);
                // std::cout << "second propagate: " << d3.count() << " ms" << std::endl;

                assert(destinations_poi_list.size() == destinations_poi_index_list.size());
                for (int j = 0; j < destinations_poi_index_list.size(); j++)
                {
                    int the_other_poi_index = destinations_poi_index_list[j];

                    // double distance;
                    // post_algorithm.best_source(destinations_poi_list[j], distance);
                    std::vector<geodesic::SurfacePoint> path;
                    post_algorithm.trace_back(destinations_poi_list[j], path);
                    changed_pairwise_path_poi_to_poi_size += path.size();
                    changed_pairwise_distance_poi_to_poi_size++;
                    // std::cout << "current_poi_index " << current_poi_index << ", the_other_poi_index " << the_other_poi_index << std::endl;

                    if (current_poi_index <= the_other_poi_index)
                    {
                        pairwise_distance_poi_to_poi[current_poi_index][the_other_poi_index - current_poi_index] = length(path);
                        pairwise_distance_poi_to_poi_changed[current_poi_index][the_other_poi_index - current_poi_index] = true;
                        pairwise_path_poi_to_poi[current_poi_index][the_other_poi_index - current_poi_index] = path;
                    }
                    else
                    {
                        pairwise_distance_poi_to_poi[the_other_poi_index][current_poi_index - the_other_poi_index] = length(path);
                        pairwise_distance_poi_to_poi_changed[the_other_poi_index][current_poi_index - the_other_poi_index] = true;
                        std::reverse(path.begin(), path.end());
                        pairwise_path_poi_to_poi[the_other_poi_index][current_poi_index - the_other_poi_index] = path;
                    }
                }
                break;
            }
        }
    }
    // st2 = std::chrono::high_resolution_clock::now();
    // d = std::chrono::duration_cast<std::chrono::milliseconds>(st2 - st1);
    // std::cout << d.count() << " ms" << std::endl;

    // st1 = std::chrono::high_resolution_clock::now();

    // std::cout << "post_algorithm.get_memory(): " << post_algorithm.get_memory() << " , changed_face_index_list.size() * sizeof(int): " << changed_face_index_list.size() * sizeof(int) << std::endl;
    memory_usage += post_algorithm.get_memory() + changed_face_index_list.size() * sizeof(int); // + changed_pairwise_distance_poi_to_poi_size * sizeof(double) + changed_pairwise_path_poi_to_poi_size * sizeof(geodesic::SurfacePoint);
    // std::cout << "changed_pairwise_distance_poi_to_poi_size: " << changed_pairwise_distance_poi_to_poi_size << std::endl;

    // st2 = std::chrono::high_resolution_clock::now();
    // d = std::chrono::duration_cast<std::chrono::milliseconds>(st2 - st1);
    // std::cout << "memory count time: " << d.count() << " ms" << std::endl;

    auto stop_updating_time = std::chrono::high_resolution_clock::now();
    auto duration_updating_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_updating_time - start_updating_time);
    updating_time = duration_updating_time.count();

    // for (int i = 0; i < poi_num; i++)
    // {
    //     for (int j = 0; j < pairwise_distance_poi_to_poi[i].size(); j++)
    //     {
    //         std::cout << "i: " << i << " , j: " << j << " , distance: " << pairwise_distance_poi_to_poi[i][j] << " , changed: " << pairwise_distance_poi_to_poi_changed[i][j] << std::endl;
    //     }
    // }

    for (int i = 0; i < changed_poi_index_list.size(); i++)
    {
        std::cout << changed_poi_index_list[i] << " ";
    }
    std::cout << std::endl;
}

void get_pairwise_distance_and_path_poi_to_poi_map(std::vector<std::vector<double>> pairwise_distance_poi_to_poi,
                                                   std::unordered_map<int, double> &pairwise_distance_poi_to_poi_map,
                                                   std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> pairwise_path_poi_to_poi,
                                                   std::unordered_map<int, std::vector<geodesic::SurfacePoint>> &pairwise_path_poi_to_poi_map,
                                                   double &complete_graph_size, int &complete_graph_edge_num,
                                                   double &complete_graph_weight, double &hash_mapping_time,
                                                   double &hash_mapping_memory_usage)
{
    auto start_hash_mapping_time = std::chrono::high_resolution_clock::now();

    int pairwise_path_poi_to_poi_size = 0;

    pairwise_distance_poi_to_poi_map.clear();
    pairwise_path_poi_to_poi_map.clear();
    for (int i = 0; i < pairwise_distance_poi_to_poi.size(); i++)
    {
        for (int j = 0; j < pairwise_distance_poi_to_poi[i].size(); j++)
        {
            int i_j;
            hash_function_two_keys_to_one_key(pairwise_distance_poi_to_poi.size(), i, i + j, i_j);
            pairwise_distance_poi_to_poi_map[i_j] = pairwise_distance_poi_to_poi[i][j];
            pairwise_path_poi_to_poi_map[i_j] = pairwise_path_poi_to_poi[i][j];
            // std::cout << "i: " << i << ", j: " << i + j << ", distance: " << pairwise_distance_poi_to_poi[i][j] << std::endl;

            if (j == 0)
            {
                continue;
            }
            complete_graph_edge_num++;
            complete_graph_weight += pairwise_distance_poi_to_poi[i][j];
            pairwise_path_poi_to_poi_size += pairwise_path_poi_to_poi[i][j].size();
        }
    }
    complete_graph_size = complete_graph_edge_num * sizeof(double) + pairwise_path_poi_to_poi_size * sizeof(geodesic::SurfacePoint);
    hash_mapping_memory_usage += complete_graph_edge_num * sizeof(double) + pairwise_path_poi_to_poi_size * sizeof(geodesic::SurfacePoint);
    // std::cout << "Complete graph size: " << complete_graph_size << std::endl;

    auto stop_hash_mapping_time = std::chrono::high_resolution_clock::now();
    auto duration_hash_mapping_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_hash_mapping_time - start_hash_mapping_time);
    hash_mapping_time = duration_hash_mapping_time.count();
    hash_mapping_time /= 1000;
}

void greedy_spanner(std::vector<std::vector<double>> pairwise_distance_poi_to_poi, double epsilon,
                    std::unordered_map<int, double> &pairwise_distance_poi_to_poi_greedy_spanner_map,
                    std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> pairwise_path_poi_to_poi,
                    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> &pairwise_path_poi_to_poi_greedy_spanner_map,
                    double &greedy_spanner_size, int &greedy_spanner_edge_num, double &greedy_spanner_weight,
                    double &GS_time, double &GS_memory_usage)
{
    auto start_GS_time = std::chrono::high_resolution_clock::now();

    std::vector<std::pair<double, std::pair<int, std::pair<int, int>>>> min_to_max_pairwise_distance_poi_to_poi;
    min_to_max_pairwise_distance_poi_to_poi.clear();

    int pairwise_path_poi_to_poi_size = 0;

    // for (int i = 0; i < pairwise_distance_poi_to_poi.size(); i++)
    // {
    //     for (int j = 0; j < pairwise_distance_poi_to_poi[i].size(); j++)
    //     {
    //         std::cout << "i: " << i << " , j: " << j << " , distance: " << pairwise_distance_poi_to_poi[i][j] << std::endl;
    //     }
    // }

    for (int i = 0; i < pairwise_distance_poi_to_poi.size(); i++)
    {
        for (int j = 0; j < pairwise_distance_poi_to_poi[i].size(); j++)
        {
            if (j == 0)
            {
                continue;
            }
            min_to_max_pairwise_distance_poi_to_poi.push_back(std::make_pair(pairwise_distance_poi_to_poi[i][j], std::make_pair(pairwise_path_poi_to_poi[i][j].size(), std::make_pair(i, i + j))));
        }
    }
    std::sort(min_to_max_pairwise_distance_poi_to_poi.begin(), min_to_max_pairwise_distance_poi_to_poi.end());
    // for (int i = 0; i < min_to_max_pairwise_distance_poi_to_poi.size(); i++)
    // {
    //     std::cout << "i: " << min_to_max_pairwise_distance_poi_to_poi[i].second.second.first << ", j: " << min_to_max_pairwise_distance_poi_to_poi[i].second.second.second << ", distance: " << min_to_max_pairwise_distance_poi_to_poi[i].first << std::endl;
    // }

    Graph graph(pairwise_distance_poi_to_poi.size());
    for (int i = 0; i < min_to_max_pairwise_distance_poi_to_poi.size(); i++)
    {
        std::vector<double> current_poi_to_other_poi_distance_greedy_spanner(pairwise_distance_poi_to_poi.size(), INF);
        std::vector<std::vector<int>> current_poi_to_other_poi_path_index_greedy_spanner(pairwise_distance_poi_to_poi.size());
        // graph.shortest_path_Dijkstra(min_to_max_pairwise_distance_poi_to_poi[i].second.second.first, min_to_max_pairwise_distance_poi_to_poi[i].second.second.second, current_poi_to_other_poi_distance_greedy_spanner);
        graph.shortest_distance_Dijkstra(min_to_max_pairwise_distance_poi_to_poi[i].second.second.first, current_poi_to_other_poi_distance_greedy_spanner);
        double distance_on_graph = current_poi_to_other_poi_distance_greedy_spanner[min_to_max_pairwise_distance_poi_to_poi[i].second.second.second];
        // std::cout << "i: " << min_to_max_pairwise_distance_poi_to_poi[i].second.second.first << ", j: " << min_to_max_pairwise_distance_poi_to_poi[i].second.second.second << ", real distance: " << min_to_max_pairwise_distance_poi_to_poi[i].first << ", distance_on_graph: " << distance_on_graph << std::endl;

        if (distance_on_graph > (1 + epsilon) * min_to_max_pairwise_distance_poi_to_poi[i].first)
        {
            // std::cout << "^^^ added" << std::endl;
            graph.add_edge_and_geo_path_size_Dijkstra(min_to_max_pairwise_distance_poi_to_poi[i].second.second.first, min_to_max_pairwise_distance_poi_to_poi[i].second.second.second, min_to_max_pairwise_distance_poi_to_poi[i].first, min_to_max_pairwise_distance_poi_to_poi[i].second.first);
        }
    }

    pairwise_distance_poi_to_poi_greedy_spanner_map.clear();
    pairwise_path_poi_to_poi_greedy_spanner_map.clear();

    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> pairwise_path_poi_to_poi_map;
    pairwise_path_poi_to_poi_map.clear();
    for (int i = 0; i < pairwise_distance_poi_to_poi.size(); i++)
    {
        for (int j = 0; j < pairwise_distance_poi_to_poi[i].size(); j++)
        {
            int i_j;
            hash_function_two_keys_to_one_key(pairwise_distance_poi_to_poi.size(), i, i + j, i_j);
            pairwise_path_poi_to_poi_map[i_j] = pairwise_path_poi_to_poi[i][j];
            // std::cout << "i: " << i << ", j: " << i + j << ", distance: " << pairwise_distance_poi_to_poi[i][j] << std::endl;

            if (j == 0)
            {
                continue;
            }
        }
    }

    for (int i = 0; i < pairwise_distance_poi_to_poi.size(); i++)
    {
        std::vector<double> current_poi_to_other_poi_distance_greedy_spanner(pairwise_distance_poi_to_poi.size(), INF);
        std::vector<std::vector<int>> current_poi_to_other_poi_path_index_greedy_spanner(pairwise_distance_poi_to_poi.size());
        graph.shortest_path_Dijkstra(i, current_poi_to_other_poi_distance_greedy_spanner, current_poi_to_other_poi_path_index_greedy_spanner, INF);

        // for (int j = 0; j < current_poi_to_other_poi_distance_greedy_spanner.size(); j++)
        // {
        //     std::cout << "$$   i: " << i << ", j: " << j << ", approximate distance: " << current_poi_to_other_poi_distance_greedy_spanner[j] << std::endl;
        // }

        for (int j = i; j < pairwise_distance_poi_to_poi.size(); j++)
        {
            int i_j;
            hash_function_two_keys_to_one_key(pairwise_distance_poi_to_poi.size(), i, j, i_j);
            // pairwise_distance_poi_to_poi_greedy_spanner_map.insert(std::pair<int, double>(i_j, current_poi_to_other_poi_distance_greedy_spanner[j]));
            pairwise_distance_poi_to_poi_greedy_spanner_map[i_j] = current_poi_to_other_poi_distance_greedy_spanner[j];

            assert(current_poi_to_other_poi_distance_greedy_spanner[j] < INF);
            std::vector<geodesic::SurfacePoint> current_poi_to_other_poi_path_greedy_spanner;
            current_poi_to_other_poi_path_greedy_spanner.clear();

            // std::cout << "dist: " << current_poi_to_other_poi_distance_greedy_spanner[j] << ", src: " << i << ", dest: " << j << ", path: ";
            // for (int k = 0; k < current_poi_to_other_poi_path_index_greedy_spanner[j].size(); k++)
            // {
            //     std::cout << current_poi_to_other_poi_path_index_greedy_spanner[j][k] << " ";
            // }
            // std::cout << std::endl;

            for (int k = 0; k < current_poi_to_other_poi_path_index_greedy_spanner[j].size() - 1; k++)
            {
                int dest_index = current_poi_to_other_poi_path_index_greedy_spanner[j][k];
                int src_index = current_poi_to_other_poi_path_index_greedy_spanner[j][k + 1];
                int src_dest_index;
                bool reverse_path = false;
                if (src_index > dest_index)
                {
                    int temp = dest_index;
                    dest_index = src_index;
                    src_index = temp;
                    reverse_path = true;
                }
                // std::cout << "(" << dest_index << ", " << src_index << ") ";

                hash_function_two_keys_to_one_key(pairwise_distance_poi_to_poi.size(), src_index, dest_index, src_dest_index);
                // std::cout << length(pairwise_path_poi_to_poi_map[src_dest_index]) << ", ";
                if (reverse_path)
                {
                    for (int m = pairwise_path_poi_to_poi_map[src_dest_index].size() - 1; m >= 0; m--)
                    {
                        current_poi_to_other_poi_path_greedy_spanner.push_back(pairwise_path_poi_to_poi_map[src_dest_index][m]);
                    }
                }
                else
                {
                    for (int m = 0; m < pairwise_path_poi_to_poi_map[src_dest_index].size(); m++)
                    {
                        current_poi_to_other_poi_path_greedy_spanner.push_back(pairwise_path_poi_to_poi_map[src_dest_index][m]);
                    }
                }
            }
            // std::cout << std::endl;
            pairwise_path_poi_to_poi_greedy_spanner_map[i_j] = current_poi_to_other_poi_path_greedy_spanner;
            pairwise_path_poi_to_poi_size += current_poi_to_other_poi_path_greedy_spanner.size();
            // std::cout << length(pairwise_path_poi_to_poi_greedy_spanner_map[i_j]) << std::endl;
            // std::cout << "i: " << i << ", j: " << j << ", approximate distance: " << current_poi_to_other_poi_distance_greedy_spanner[j] << ", real distance: " << pairwise_distance_poi_to_poi[i][j - i] << ", ratio: " << current_poi_to_other_poi_distance_greedy_spanner[j] / pairwise_distance_poi_to_poi[i][j - i] << std::endl;
        }
    }
    greedy_spanner_edge_num = graph.get_edge_num_Dijkstra();
    greedy_spanner_weight = graph.get_total_weight_Dijkstra();
    greedy_spanner_size = greedy_spanner_edge_num * sizeof(double) + graph.get_total_geo_path_size_Dijkstra() * sizeof(geodesic::SurfacePoint);
    GS_memory_usage += 0.5 * pairwise_distance_poi_to_poi.size() * (pairwise_distance_poi_to_poi.size() - 1) * sizeof(double) + pairwise_path_poi_to_poi_size * sizeof(geodesic::SurfacePoint);
    // std::cout << "Complete graph size: " << min_to_max_pairwise_distance_poi_to_poi.size() << ", Greedy spanner size: " << graph.get_edge_num_Dijkstra() << ", Greedy spanner total weight: " << graph.get_total_weight_Dijkstra() << std::endl;

    auto stop_GS_time = std::chrono::high_resolution_clock::now();
    auto duration_GS_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_GS_time - start_GS_time);
    GS_time = duration_GS_time.count();
    GS_time /= 1000;
    // std::cout << "Total GS time: " << GS_time << " ms" << std::endl;
}

void hierarchy_greedy_spanner(std::vector<std::vector<double>> pairwise_distance_poi_to_poi, double epsilon,
                              std::unordered_map<int, double> &pairwise_distance_poi_to_poi_hierarchy_greedy_spanner_map,
                              std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> pairwise_path_poi_to_poi,
                              std::unordered_map<int, std::vector<geodesic::SurfacePoint>> &pairwise_path_poi_to_poi_hierarchy_greedy_spanner_map,
                              double &hierarchy_greedy_spanner_size, int &hierarchy_greedy_spanner_edge_num,
                              double &hierarchy_greedy_spanner_weight, double &HGS_time, double &HGS_memory_usage)
{
    auto start_HGS_time = std::chrono::high_resolution_clock::now();

    std::vector<std::pair<double, std::pair<int, std::pair<int, int>>>> min_to_max_pairwise_distance_poi_to_poi;
    min_to_max_pairwise_distance_poi_to_poi.clear();

    int pairwise_path_poi_to_poi_size = 0;

    // for (int i = 0; i < pairwise_distance_poi_to_poi.size(); i++)
    // {
    //     for (int j = 0; j < pairwise_distance_poi_to_poi[i].size(); j++)
    //     {
    //         std::cout << "i: " << i << " , j: " << j << " , distance: " << pairwise_distance_poi_to_poi[i][j] << std::endl;
    //     }
    // }

    for (int i = 0; i < pairwise_distance_poi_to_poi.size(); i++)
    {
        for (int j = 0; j < pairwise_distance_poi_to_poi[i].size(); j++)
        {
            if (j == 0)
            {
                continue;
            }
            min_to_max_pairwise_distance_poi_to_poi.push_back(std::make_pair(pairwise_distance_poi_to_poi[i][j], std::make_pair(pairwise_path_poi_to_poi[i][j].size(), std::make_pair(i, i + j))));
        }
    }
    std::sort(min_to_max_pairwise_distance_poi_to_poi.begin(), min_to_max_pairwise_distance_poi_to_poi.end());

    // for (int i = 0; i < min_to_max_pairwise_distance_poi_to_poi.size(); i++)
    // {
    //     std::cout << "i: " << min_to_max_pairwise_distance_poi_to_poi[i].second.second.first << ", j: " << min_to_max_pairwise_distance_poi_to_poi[i].second.second.second << ", distance: " << min_to_max_pairwise_distance_poi_to_poi[i].first << std::endl;
    // }
    double max_pairwise_distance = min_to_max_pairwise_distance_poi_to_poi[min_to_max_pairwise_distance_poi_to_poi.size() - 1].first;
    // std::cout << "max_pairwise_distance: " << max_pairwise_distance << std::endl;

    Graph graph(pairwise_distance_poi_to_poi.size());
    std::vector<std::pair<double, std::pair<int, std::pair<int, int>>>> added_edge_in_greedy_spanner;

    std::vector<std::vector<std::pair<double, std::pair<int, std::pair<int, int>>>>> distance_interval;
    std::vector<std::pair<double, std::pair<int, std::pair<int, int>>>> one_distance_interval;
    distance_interval.clear();
    one_distance_interval.clear();

    int interval_zero_item_num = 0;
    for (int i = 0; i < min_to_max_pairwise_distance_poi_to_poi.size(); i++)
    {
        if (min_to_max_pairwise_distance_poi_to_poi[i].first > 0 &&
            min_to_max_pairwise_distance_poi_to_poi[i].first <= max_pairwise_distance / pairwise_distance_poi_to_poi.size())
        {
            // std::cout << "==very small edge== i: " << min_to_max_pairwise_distance_poi_to_poi[i].second.second.first << ", j: " << min_to_max_pairwise_distance_poi_to_poi[i].second.second.second << ", real distance: " << min_to_max_pairwise_distance_poi_to_poi[i].first << std::endl;
            // std::cout << "^^^ added" << std::endl;
            graph.add_edge_and_geo_path_size_Dijkstra(min_to_max_pairwise_distance_poi_to_poi[i].second.second.first, min_to_max_pairwise_distance_poi_to_poi[i].second.second.second, min_to_max_pairwise_distance_poi_to_poi[i].first, min_to_max_pairwise_distance_poi_to_poi[i].second.first);
            added_edge_in_greedy_spanner.push_back(min_to_max_pairwise_distance_poi_to_poi[i]);
            one_distance_interval.push_back(min_to_max_pairwise_distance_poi_to_poi[i]);
            interval_zero_item_num++;
        }
        else
        {
            break;
        }
    }
    distance_interval.push_back(one_distance_interval);

    int interval_index = 1;
    one_distance_interval.clear();

    // std::cout << "interval_zero_item_num: " << interval_zero_item_num << std::endl;

    for (int i = interval_zero_item_num; i < min_to_max_pairwise_distance_poi_to_poi.size(); i++)
    {
        if (min_to_max_pairwise_distance_poi_to_poi[i].first > pow(2, (interval_index - 1)) * max_pairwise_distance / pairwise_distance_poi_to_poi.size() &&
            min_to_max_pairwise_distance_poi_to_poi[i].first <= pow(2, interval_index) * max_pairwise_distance / pairwise_distance_poi_to_poi.size())
        {
            one_distance_interval.push_back(min_to_max_pairwise_distance_poi_to_poi[i]);
            // std::cout << "i: " << i << ", distance: " << min_to_max_pairwise_distance_poi_to_poi[i].first << std::endl;
        }
        else
        {
            distance_interval.push_back(one_distance_interval);
            interval_index++;
            one_distance_interval.clear();
            one_distance_interval.push_back(min_to_max_pairwise_distance_poi_to_poi[i]);
        }
    }
    distance_interval.push_back(one_distance_interval);

    // for (int i = 0; i < distance_interval.size(); i++)
    // {
    //     for (int j = 0; j < distance_interval[i].size(); j++)
    //     {
    //         std::cout << "i: " << i << ", interval distance: " << distance_interval[i][j].first << std::endl;
    //     }
    // }

    double W = max_pairwise_distance / pairwise_distance_poi_to_poi.size();
    double delta = 0.5 * ((sqrt(epsilon + 1) - 1) / (sqrt(epsilon + 1) + 3)) * (log(epsilon / 4) / log(0.875));
    // std::cout << "delta: " << delta << std::endl;
    // double delta = 0.17;
    // double delta = 0.49;

    for (int k = 1; k <= ceil(log2(pairwise_distance_poi_to_poi.size())); k++)
    {
        // std::cout << "k: " << k << std::endl;
        // std::cout << "W: " << W << ", 2 * W: " << 2 * W << ", delta * W: " << delta * W << std::endl;

        // hierarchy graph
        std::vector<int> unprocessed_poi;
        unprocessed_poi.clear();

        for (int i = 0; i < pairwise_distance_poi_to_poi.size(); i++)
        {
            unprocessed_poi.push_back(i);
        }

        // std::unordered_map<int, std::vector<std::pair<int, double>>> centers;
        std::vector<int> centers;
        std::unordered_map<int, int> centers_unordered_map;
        std::unordered_map<int, std::pair<int, double>> non_centers;
        centers.clear();
        centers_unordered_map.clear();
        non_centers.clear();

        int count = 0;
        auto ite = unprocessed_poi.begin();
        while (ite != unprocessed_poi.end())
        {
            // std::vector<std::pair<int, double>> center_coveres_poi;
            // std::vector<std::pair<int, double>> non_center_covered_by_poi;
            // center_coveres_poi.clear();
            // non_center_covered_by_poi.clear();

            int current_index = *ite;
            std::vector<double> current_poi_to_other_poi_distance_greedy_spanner(pairwise_distance_poi_to_poi.size(), INF);
            graph.shortest_distance_Dijkstra(current_index, current_poi_to_other_poi_distance_greedy_spanner, delta * W);
            unprocessed_poi.erase(ite);
            centers.push_back(current_index);
            centers_unordered_map[current_index] = count;
            count++;

            auto ite2 = unprocessed_poi.begin();
            while (ite2 != unprocessed_poi.end())
            {
                if (current_poi_to_other_poi_distance_greedy_spanner[*ite2] <= delta * W &&
                    current_poi_to_other_poi_distance_greedy_spanner[*ite2] > 0)
                {
                    // centers[current_index].push_back(std::make_pair(*ite2, current_poi_to_other_poi_distance_greedy_spanner[*ite2]));
                    non_centers[*ite2] = std::make_pair(current_index, current_poi_to_other_poi_distance_greedy_spanner[*ite2]);
                    unprocessed_poi.erase(ite2);
                }
                else
                {
                    ite2++;
                }
            }
        }

        // for (int i = 0; i < centers.size(); i++)
        // {
        //     std::cout << "center: " << centers[i] << std::endl;
        // }

        // std::cout << "center size: " << centers.size() << std::endl;

        // calculate inter-hierarchy edges of first type
        std::unordered_map<int, double> pairwise_distance_center_to_center;
        pairwise_distance_center_to_center.clear();
        std::unordered_map<int, double> potential_second_type_inter_edge;
        potential_second_type_inter_edge.clear();
        Graph hierarchy_graph(centers.size());

        for (int i = 0; i < centers.size(); i++)
        {
            std::vector<double> current_center_to_other_center_and_non_center_distance_greedy_spanner(pairwise_distance_poi_to_poi.size(), INF);
            graph.shortest_distance_Dijkstra(centers[i], current_center_to_other_center_and_non_center_distance_greedy_spanner, W + 2 * W * delta);
            for (int j = i; j < centers.size(); j++)
            {
                // std::cout << "first type inter hierarchy edge -- i: " << centers[i] << " , j: " << centers[j] << ", first type inter hierarchy length: " << current_center_to_other_center_and_non_center_distance_greedy_spanner[centers[j]] << std::endl;

                if (current_center_to_other_center_and_non_center_distance_greedy_spanner[centers[j]] < W &&
                    current_center_to_other_center_and_non_center_distance_greedy_spanner[centers[j]] > 0)
                {
                    // std::cout << "^^^ added first type inter-hierarchy edge" << std::endl;
                    int i_j;
                    hash_function_two_keys_to_one_key(centers.size(), i, j, i_j);
                    pairwise_distance_center_to_center[i_j] = current_center_to_other_center_and_non_center_distance_greedy_spanner[centers[j]];
                    hierarchy_graph.add_edge_Dijkstra(i, j, current_center_to_other_center_and_non_center_distance_greedy_spanner[centers[j]]);
                }
                else if (current_center_to_other_center_and_non_center_distance_greedy_spanner[centers[j]] < W + 2 * W * delta &&
                         current_center_to_other_center_and_non_center_distance_greedy_spanner[centers[j]] >= W)
                {
                    // std::cout << "^^^ potential second type inter-hierarchy edge" << std::endl;
                    int i_j;
                    hash_function_two_keys_to_one_key(centers.size(), i, j, i_j);
                    potential_second_type_inter_edge[i_j] = current_center_to_other_center_and_non_center_distance_greedy_spanner[centers[j]];
                    // hierarchy_graph.add_edge_Dijkstra(i, j, current_center_to_other_center_and_non_center_distance_greedy_spanner[centers[j]]);
                }
            }
        }

        // calculate inter-hierarchy edges of second type
        for (int i = 0; i < added_edge_in_greedy_spanner.size(); i++)
        {
            int one_endpoint_index = added_edge_in_greedy_spanner[i].second.second.first;
            int another_endpoint_index = added_edge_in_greedy_spanner[i].second.second.second;

            // if both two endpoints of an edge is the center
            if (centers_unordered_map.count(one_endpoint_index) != 0 && centers_unordered_map.count(another_endpoint_index) != 0)
            {
                int i_j;
                if (one_endpoint_index > another_endpoint_index)
                {
                    int temp = another_endpoint_index;
                    another_endpoint_index = one_endpoint_index;
                    one_endpoint_index = temp;
                }
                hash_function_two_keys_to_one_key(centers.size(), centers_unordered_map[one_endpoint_index], centers_unordered_map[another_endpoint_index], i_j);
                if (potential_second_type_inter_edge.count(i_j) != 0)
                {
                    if (pairwise_distance_center_to_center.count(i_j) == 0)
                    {
                        // std::cout << "second type inter hierarchy edge ==both center== center i: " << one_endpoint_index << ", center j: " << another_endpoint_index << ", original i: " << one_endpoint_index << ", original j: " << another_endpoint_index << ", original edge distance: " << added_edge_in_greedy_spanner[i].first << ", second type inter hierarchy distance: " << added_edge_in_greedy_spanner[i].first << std::endl;
                        // std::cout << "^^^ added inter-hierarchy edge second type" << std::endl;
                        pairwise_distance_center_to_center[i_j] = added_edge_in_greedy_spanner[i].first;
                        hierarchy_graph.add_edge_Dijkstra(centers_unordered_map[one_endpoint_index], centers_unordered_map[another_endpoint_index], added_edge_in_greedy_spanner[i].first);
                    }
                    else if (pairwise_distance_center_to_center.count(i_j) != 0 && pairwise_distance_center_to_center[i_j] > added_edge_in_greedy_spanner[i].first)
                    {
                        // std::cout << "second type inter hierarchy edge ==both center== center i: " << one_endpoint_index << ", center j: " << another_endpoint_index << ", original i: " << one_endpoint_index << ", original j: " << another_endpoint_index << ", original edge distance: " << added_edge_in_greedy_spanner[i].first << ", second type inter hierarchy distance: " << added_edge_in_greedy_spanner[i].first << std::endl;
                        // std::cout << "^^^ updated inter-hierarchy edge second type" << std::endl;
                        pairwise_distance_center_to_center[i_j] = added_edge_in_greedy_spanner[i].first;
                        hierarchy_graph.update_edge_Dijkstra(centers_unordered_map[one_endpoint_index], centers_unordered_map[another_endpoint_index], added_edge_in_greedy_spanner[i].first);
                    }
                }
            }
            // if one endpoint of an edge is the non center and another enpoint is the center
            else if (centers_unordered_map.count(one_endpoint_index) == 0 && centers_unordered_map.count(another_endpoint_index) != 0)
            {
                int one_non_center_endpoint_center_index = non_centers[one_endpoint_index].first;
                double one_non_center_endpoint_to_center_distance = non_centers[one_endpoint_index].second;

                int i_j;
                if (one_non_center_endpoint_center_index > another_endpoint_index)
                {
                    int temp = another_endpoint_index;
                    another_endpoint_index = one_non_center_endpoint_center_index;
                    one_non_center_endpoint_center_index = temp;
                }
                hash_function_two_keys_to_one_key(centers.size(), centers_unordered_map[one_non_center_endpoint_center_index], centers_unordered_map[another_endpoint_index], i_j);
                if (potential_second_type_inter_edge.count(i_j) != 0)
                {
                    if (pairwise_distance_center_to_center.count(i_j) == 0)
                    {
                        // std::cout << "second type inter hierarchy edge ==one non-center, another center== center i: " << one_endpoint_index << ", center j: " << another_endpoint_index << ", original i: " << one_non_center_endpoint_center_index << ", original j: " << another_endpoint_index << ", original edge distance: " << added_edge_in_greedy_spanner[i].first << ", second type inter hierarchy distance: " << added_edge_in_greedy_spanner[i].first + one_non_center_endpoint_to_center_distance << std::endl;
                        // std::cout << "^^^ added inter-hierarchy edge second type" << std::endl;
                        pairwise_distance_center_to_center[i_j] = added_edge_in_greedy_spanner[i].first + one_non_center_endpoint_to_center_distance;
                        hierarchy_graph.add_edge_Dijkstra(centers_unordered_map[one_non_center_endpoint_center_index], centers_unordered_map[another_endpoint_index], added_edge_in_greedy_spanner[i].first + one_non_center_endpoint_to_center_distance);
                    }
                    else if (pairwise_distance_center_to_center.count(i_j) != 0 && pairwise_distance_center_to_center[i_j] > added_edge_in_greedy_spanner[i].first + one_non_center_endpoint_to_center_distance)
                    {
                        // std::cout << "second type inter hierarchy edge ==one non-center, another center== center i: " << one_endpoint_index << ", center j: " << another_endpoint_index << ", original i: " << one_non_center_endpoint_center_index << ", original j: " << another_endpoint_index << ", original edge distance: " << added_edge_in_greedy_spanner[i].first << ", second type inter hierarchy distance: " << added_edge_in_greedy_spanner[i].first + one_non_center_endpoint_to_center_distance << std::endl;
                        // std::cout << "^^^ updated inter-hierarchy edge second type" << std::endl;
                        pairwise_distance_center_to_center[i_j] = added_edge_in_greedy_spanner[i].first + one_non_center_endpoint_to_center_distance;
                        hierarchy_graph.update_edge_Dijkstra(centers_unordered_map[one_non_center_endpoint_center_index], centers_unordered_map[another_endpoint_index], added_edge_in_greedy_spanner[i].first + one_non_center_endpoint_to_center_distance);
                    }
                }
            }
            // if one endpoint of an edge is the center and another enpoint is the non center
            else if (centers_unordered_map.count(one_endpoint_index) != 0 && centers_unordered_map.count(another_endpoint_index) == 0)
            {

                int another_non_center_endpoint_center_index = non_centers[another_endpoint_index].first;
                double another_non_center_endpoint_to_center_distance = non_centers[another_endpoint_index].second;

                int i_j;
                if (one_endpoint_index > another_non_center_endpoint_center_index)
                {
                    int temp = another_non_center_endpoint_center_index;
                    another_non_center_endpoint_center_index = one_endpoint_index;
                    one_endpoint_index = temp;
                }
                hash_function_two_keys_to_one_key(centers.size(), centers_unordered_map[one_endpoint_index], centers_unordered_map[another_non_center_endpoint_center_index], i_j);
                if (potential_second_type_inter_edge.count(i_j) != 0)
                {
                    if (pairwise_distance_center_to_center.count(i_j) == 0)
                    {
                        // std::cout << "second type inter hierarchy edge ==one center, another non-center== center i: " << one_endpoint_index << ", center j: " << another_non_center_endpoint_center_index << ", original i: " << one_endpoint_index << ", original j: " << another_endpoint_index << ", original edge distance: " << added_edge_in_greedy_spanner[i].first << ", second type inter hierarchy distance: " << added_edge_in_greedy_spanner[i].first + another_non_center_endpoint_to_center_distance << std::endl;
                        // std::cout << "^^^ added inter-hierarchy edge second type" << std::endl;
                        pairwise_distance_center_to_center[i_j] = added_edge_in_greedy_spanner[i].first + another_non_center_endpoint_to_center_distance;
                        hierarchy_graph.add_edge_Dijkstra(centers_unordered_map[one_endpoint_index], centers_unordered_map[another_non_center_endpoint_center_index], added_edge_in_greedy_spanner[i].first + another_non_center_endpoint_to_center_distance);
                    }
                    else if (pairwise_distance_center_to_center.count(i_j) != 0 && pairwise_distance_center_to_center[i_j] > added_edge_in_greedy_spanner[i].first + another_non_center_endpoint_to_center_distance)
                    {
                        // std::cout << "second type inter hierarchy edge ==one center, another non-center== center i: " << one_endpoint_index << ", center j: " << another_non_center_endpoint_center_index << ", original i: " << one_endpoint_index << ", original j: " << another_endpoint_index << ", original edge distance: " << added_edge_in_greedy_spanner[i].first << ", second type inter hierarchy distance: " << added_edge_in_greedy_spanner[i].first + another_non_center_endpoint_to_center_distance << std::endl;
                        // std::cout << "^^^ updated inter-hierarchy edge second type" << std::endl;
                        pairwise_distance_center_to_center[i_j] = added_edge_in_greedy_spanner[i].first + another_non_center_endpoint_to_center_distance;
                        hierarchy_graph.update_edge_Dijkstra(centers_unordered_map[one_endpoint_index], centers_unordered_map[another_non_center_endpoint_center_index], added_edge_in_greedy_spanner[i].first + another_non_center_endpoint_to_center_distance);
                    }
                }
            }
            // if both two endpoints of an edge is the non center
            else
            {
                int one_non_center_endpoint_center_index = non_centers[one_endpoint_index].first;
                double one_non_center_endpoint_to_center_distance = non_centers[one_endpoint_index].second;

                int another_non_center_endpoint_center_index = non_centers[another_endpoint_index].first;
                double another_non_center_endpoint_to_center_distance = non_centers[another_endpoint_index].second;

                int i_j;
                if (one_non_center_endpoint_center_index > another_non_center_endpoint_center_index)
                {
                    int temp = another_non_center_endpoint_center_index;
                    another_non_center_endpoint_center_index = one_non_center_endpoint_center_index;
                    one_non_center_endpoint_center_index = temp;
                }
                hash_function_two_keys_to_one_key(centers.size(), centers_unordered_map[one_non_center_endpoint_center_index], centers_unordered_map[another_non_center_endpoint_center_index], i_j);
                if (potential_second_type_inter_edge.count(i_j) != 0)
                {
                    if (pairwise_distance_center_to_center.count(i_j) == 0)
                    {
                        // std::cout << "second type inter hierarchy edge ==both non-center== center i: " << one_non_center_endpoint_center_index << ", center j: " << another_non_center_endpoint_center_index << ", original i: " << one_endpoint_index << ", original j: " << another_endpoint_index << ", original edge distance: " << added_edge_in_greedy_spanner[i].first << ", second type inter hierarchy distance: " << added_edge_in_greedy_spanner[i].first + one_non_center_endpoint_to_center_distance + another_non_center_endpoint_to_center_distance << std::endl;
                        // std::cout << "^^^ added inter-hierarchy edge second type" << std::endl;
                        pairwise_distance_center_to_center[i_j] = added_edge_in_greedy_spanner[i].first + one_non_center_endpoint_to_center_distance + another_non_center_endpoint_to_center_distance;
                        hierarchy_graph.add_edge_Dijkstra(centers_unordered_map[one_non_center_endpoint_center_index], centers_unordered_map[another_non_center_endpoint_center_index], added_edge_in_greedy_spanner[i].first + one_non_center_endpoint_to_center_distance + another_non_center_endpoint_to_center_distance);
                    }
                    else if (pairwise_distance_center_to_center.count(i_j) != 0 && pairwise_distance_center_to_center[i_j] > added_edge_in_greedy_spanner[i].first + one_non_center_endpoint_to_center_distance + another_non_center_endpoint_to_center_distance)
                    {
                        // std::cout << "second type inter hierarchy edge ==both non-center== center i: " << one_non_center_endpoint_center_index << ", center j: " << another_non_center_endpoint_center_index << ", original i: " << one_endpoint_index << ", original j: " << another_endpoint_index << ", original edge distance: " << added_edge_in_greedy_spanner[i].first << ", second type inter hierarchy distance: " << added_edge_in_greedy_spanner[i].first + one_non_center_endpoint_to_center_distance + another_non_center_endpoint_to_center_distance << std::endl;
                        // std::cout << "^^^ updated inter-hierarchy edge second type" << std::endl;
                        pairwise_distance_center_to_center[i_j] = added_edge_in_greedy_spanner[i].first + one_non_center_endpoint_to_center_distance + another_non_center_endpoint_to_center_distance;
                        hierarchy_graph.update_edge_Dijkstra(centers_unordered_map[one_non_center_endpoint_center_index], centers_unordered_map[another_non_center_endpoint_center_index], added_edge_in_greedy_spanner[i].first + one_non_center_endpoint_to_center_distance + another_non_center_endpoint_to_center_distance);
                    }
                }
            }
        }

        // for each edge in E_i
        for (int i = 0; i < distance_interval[k].size(); i++)
        {
            int one_endpoint_index = distance_interval[k][i].second.second.first;
            int another_endpoint_index = distance_interval[k][i].second.second.second;
            // std::cout << "one_endpoint_index: " << one_endpoint_index << ", another_endpoint_index: " << another_endpoint_index << std::endl;

            // if both two endpoints of an edge is the center
            if (centers_unordered_map.count(one_endpoint_index) != 0 && centers_unordered_map.count(another_endpoint_index) != 0)
            {
                int i_j;
                if (one_endpoint_index > another_endpoint_index)
                {
                    int temp = another_endpoint_index;
                    another_endpoint_index = one_endpoint_index;
                    one_endpoint_index = temp;
                }
                hash_function_two_keys_to_one_key(centers.size(), centers_unordered_map[one_endpoint_index], centers_unordered_map[another_endpoint_index], i_j);
                // std::cout << "==both center== center i: " << one_endpoint_index << ", center j: " << another_endpoint_index << ", original i: " << one_endpoint_index << ", original j: " << another_endpoint_index << ", real distance: " << distance_interval[k][i].first;

                if (pairwise_distance_center_to_center.count(i_j) == 0)
                {
                    std::vector<double> current_center_to_other_center_distance_hierarchy_graph(centers.size(), INF);
                    hierarchy_graph.shortest_distance_Dijkstra(centers_unordered_map[one_endpoint_index], centers_unordered_map[another_endpoint_index], current_center_to_other_center_distance_hierarchy_graph);
                    double distance_on_hierarchy_graph = current_center_to_other_center_distance_hierarchy_graph[centers_unordered_map[another_endpoint_index]];
                    // std::cout << ", distance on hierarchy graph: " << distance_on_hierarchy_graph << std::endl;

                    if (distance_on_hierarchy_graph > (1 + epsilon) * distance_interval[k][i].first)
                    {
                        // std::cout << "^^^ added, the direct hierarchy graph doesn't exist" << std::endl;
                        graph.add_edge_and_geo_path_size_Dijkstra(one_endpoint_index, another_endpoint_index, distance_interval[k][i].first, distance_interval[k][i].second.first);
                        added_edge_in_greedy_spanner.push_back(distance_interval[k][i]);
                        hierarchy_graph.add_edge_Dijkstra(centers_unordered_map[one_endpoint_index], centers_unordered_map[another_endpoint_index], distance_interval[k][i].first);
                        pairwise_distance_center_to_center[i_j] = distance_interval[k][i].first;
                    }
                }
                else if (pairwise_distance_center_to_center.count(i_j) != 0 && pairwise_distance_center_to_center[i_j] > (1 + epsilon) * distance_interval[k][i].first)
                {
                    // std::cout << ", distance on hierarchy graph: " << pairwise_distance_center_to_center[i_j] << std::endl;
                    // std::cout << "^^^ added, the direct hierarchy graph exists, but longer than (1+e)*edge length" << std::endl;
                    graph.add_edge_and_geo_path_size_Dijkstra(one_endpoint_index, another_endpoint_index, distance_interval[k][i].first, distance_interval[k][i].second.first);
                    added_edge_in_greedy_spanner.push_back(distance_interval[k][i]);
                    hierarchy_graph.update_edge_Dijkstra(centers_unordered_map[one_endpoint_index], centers_unordered_map[another_endpoint_index], distance_interval[k][i].first);
                    pairwise_distance_center_to_center[i_j] = distance_interval[k][i].first;
                }
                else if (pairwise_distance_center_to_center.count(i_j) != 0 && pairwise_distance_center_to_center[i_j] <= (1 + epsilon) * distance_interval[k][i].first)
                {
                    // std::cout << ", distance on hierarchy graph: " << pairwise_distance_center_to_center[i_j] << std::endl;
                    if (pairwise_distance_center_to_center[i_j] > distance_interval[k][i].first)
                    {
                        hierarchy_graph.update_edge_Dijkstra(centers_unordered_map[one_endpoint_index], centers_unordered_map[another_endpoint_index], distance_interval[k][i].first);
                        pairwise_distance_center_to_center[i_j] = distance_interval[k][i].first;
                    }
                }
            }
            // if one endpoint of an edge is the non center and another enpoint is the center
            else if (centers_unordered_map.count(one_endpoint_index) == 0 && centers_unordered_map.count(another_endpoint_index) != 0)
            {
                int one_non_center_endpoint_center_index = non_centers[one_endpoint_index].first;
                double one_non_center_endpoint_to_center_distance = non_centers[one_endpoint_index].second;

                int i_j;
                if (one_non_center_endpoint_center_index > another_endpoint_index)
                {
                    int temp = another_endpoint_index;
                    another_endpoint_index = one_non_center_endpoint_center_index;
                    one_non_center_endpoint_center_index = temp;
                }
                hash_function_two_keys_to_one_key(centers.size(), centers_unordered_map[one_non_center_endpoint_center_index], centers_unordered_map[another_endpoint_index], i_j);
                // std::cout << "==one non-center, another center== center i: " << one_non_center_endpoint_center_index << ", center j: " << another_endpoint_index << ", original i: " << one_endpoint_index << ", original j: " << another_endpoint_index << ", real distance: " << distance_interval[k][i].first;
                if (pairwise_distance_center_to_center.count(i_j) == 0)
                {
                    std::vector<double> current_center_to_other_center_distance_hierarchy_graph(centers.size(), INF);
                    hierarchy_graph.shortest_distance_Dijkstra(centers_unordered_map[one_non_center_endpoint_center_index], centers_unordered_map[another_endpoint_index], current_center_to_other_center_distance_hierarchy_graph);
                    double distance_on_hierarchy_graph = current_center_to_other_center_distance_hierarchy_graph[centers_unordered_map[another_endpoint_index]];
                    // std::cout << ", distance on hierarchy graph: " << distance_on_hierarchy_graph << std::endl;

                    if (distance_on_hierarchy_graph > (1 + epsilon) * distance_interval[k][i].first)
                    {
                        // std::cout << "^^^ added, the direct hierarchy graph doesn't exist" << std::endl;
                        graph.add_edge_and_geo_path_size_Dijkstra(one_endpoint_index, another_endpoint_index, distance_interval[k][i].first, distance_interval[k][i].second.first);
                        added_edge_in_greedy_spanner.push_back(distance_interval[k][i]);
                        hierarchy_graph.add_edge_Dijkstra(centers_unordered_map[one_non_center_endpoint_center_index], centers_unordered_map[another_endpoint_index], distance_interval[k][i].first + one_non_center_endpoint_to_center_distance);
                        pairwise_distance_center_to_center[i_j] = distance_interval[k][i].first + one_non_center_endpoint_to_center_distance;
                    }
                }
                else if (pairwise_distance_center_to_center.count(i_j) != 0 && pairwise_distance_center_to_center[i_j] > (1 + epsilon) * distance_interval[k][i].first)
                {
                    // std::cout << ", distance on hierarchy graph: " << pairwise_distance_center_to_center[i_j] << std::endl;
                    // std::cout << "^^^ added, the direct hierarchy graph exists, but longer than (1+e)*edge length" << std::endl;
                    graph.add_edge_and_geo_path_size_Dijkstra(one_endpoint_index, another_endpoint_index, distance_interval[k][i].first, distance_interval[k][i].second.first);
                    added_edge_in_greedy_spanner.push_back(distance_interval[k][i]);
                    hierarchy_graph.update_edge_Dijkstra(centers_unordered_map[one_non_center_endpoint_center_index], centers_unordered_map[another_endpoint_index], distance_interval[k][i].first + one_non_center_endpoint_to_center_distance);
                    pairwise_distance_center_to_center[i_j] = distance_interval[k][i].first + one_non_center_endpoint_to_center_distance;
                }
                else if (pairwise_distance_center_to_center.count(i_j) != 0 && pairwise_distance_center_to_center[i_j] <= (1 + epsilon) * distance_interval[k][i].first)
                {
                    // std::cout << ", distance on hierarchy graph: " << pairwise_distance_center_to_center[i_j] << std::endl;
                    if (pairwise_distance_center_to_center[i_j] > distance_interval[k][i].first + one_non_center_endpoint_to_center_distance)
                    {
                        hierarchy_graph.update_edge_Dijkstra(centers_unordered_map[one_non_center_endpoint_center_index], centers_unordered_map[another_endpoint_index], distance_interval[k][i].first + one_non_center_endpoint_to_center_distance);
                        pairwise_distance_center_to_center[i_j] = distance_interval[k][i].first + one_non_center_endpoint_to_center_distance;
                    }
                }
            }
            // if one endpoint of an edge is the center and another enpoint is the non center
            else if (centers_unordered_map.count(one_endpoint_index) != 0 && centers_unordered_map.count(another_endpoint_index) == 0)
            {
                int another_non_center_endpoint_center_index = non_centers[another_endpoint_index].first;
                double another_non_center_endpoint_to_center_distance = non_centers[another_endpoint_index].second;

                int i_j;
                if (one_endpoint_index > another_non_center_endpoint_center_index)
                {
                    int temp = another_non_center_endpoint_center_index;
                    another_non_center_endpoint_center_index = one_endpoint_index;
                    one_endpoint_index = temp;
                }
                hash_function_two_keys_to_one_key(centers.size(), centers_unordered_map[one_endpoint_index], centers_unordered_map[another_non_center_endpoint_center_index], i_j);
                // std::cout << "==one center, another non-center== center i: " << one_endpoint_index << ", center j: " << another_non_center_endpoint_center_index << ", original i: " << one_endpoint_index << ", original j: " << another_endpoint_index << ", real distance: " << distance_interval[k][i].first;

                if (pairwise_distance_center_to_center.count(i_j) == 0)
                {
                    std::vector<double> current_center_to_other_center_distance_hierarchy_graph(centers.size(), INF);
                    hierarchy_graph.shortest_distance_Dijkstra(centers_unordered_map[one_endpoint_index], centers_unordered_map[another_non_center_endpoint_center_index], current_center_to_other_center_distance_hierarchy_graph);
                    double distance_on_hierarchy_graph = current_center_to_other_center_distance_hierarchy_graph[centers_unordered_map[another_non_center_endpoint_center_index]];
                    // std::cout << ", distance on hierarchy graph: " << distance_on_hierarchy_graph << std::endl;

                    if (distance_on_hierarchy_graph > (1 + epsilon) * distance_interval[k][i].first)
                    {
                        // std::cout << "^^^ added, the direct hierarchy graph doesn't exist" << std::endl;
                        graph.add_edge_and_geo_path_size_Dijkstra(one_endpoint_index, another_endpoint_index, distance_interval[k][i].first, distance_interval[k][i].second.first);
                        added_edge_in_greedy_spanner.push_back(distance_interval[k][i]);
                        hierarchy_graph.add_edge_Dijkstra(centers_unordered_map[one_endpoint_index], centers_unordered_map[another_non_center_endpoint_center_index], distance_interval[k][i].first + another_non_center_endpoint_to_center_distance);
                        pairwise_distance_center_to_center[i_j] = distance_interval[k][i].first + another_non_center_endpoint_to_center_distance;
                    }
                }
                else if (pairwise_distance_center_to_center.count(i_j) != 0 && pairwise_distance_center_to_center[i_j] > (1 + epsilon) * distance_interval[k][i].first)
                {
                    // std::cout << ", distance on hierarchy graph: " << pairwise_distance_center_to_center[i_j] << std::endl;
                    // std::cout << "^^^ added, the direct hierarchy graph exists, but longer than (1+e)*edge length" << std::endl;
                    graph.add_edge_and_geo_path_size_Dijkstra(one_endpoint_index, another_endpoint_index, distance_interval[k][i].first, distance_interval[k][i].second.first);
                    added_edge_in_greedy_spanner.push_back(distance_interval[k][i]);
                    hierarchy_graph.update_edge_Dijkstra(centers_unordered_map[one_endpoint_index], centers_unordered_map[another_non_center_endpoint_center_index], distance_interval[k][i].first + another_non_center_endpoint_to_center_distance);
                    pairwise_distance_center_to_center[i_j] = distance_interval[k][i].first + another_non_center_endpoint_to_center_distance;
                }
                else if (pairwise_distance_center_to_center.count(i_j) != 0 && pairwise_distance_center_to_center[i_j] <= (1 + epsilon) * distance_interval[k][i].first)
                {
                    // std::cout << ", distance on hierarchy graph: " << pairwise_distance_center_to_center[i_j] << std::endl;
                    if (pairwise_distance_center_to_center[i_j] > distance_interval[k][i].first + another_non_center_endpoint_to_center_distance)
                    {
                        hierarchy_graph.update_edge_Dijkstra(centers_unordered_map[one_endpoint_index], centers_unordered_map[another_non_center_endpoint_center_index], distance_interval[k][i].first + another_non_center_endpoint_to_center_distance);
                        pairwise_distance_center_to_center[i_j] = distance_interval[k][i].first + another_non_center_endpoint_to_center_distance;
                    }
                }
            }
            // if both two endpoints of an edge is the non center
            else
            {
                int one_non_center_endpoint_center_index = non_centers[one_endpoint_index].first;
                double one_non_center_endpoint_to_center_distance = non_centers[one_endpoint_index].second;

                int another_non_center_endpoint_center_index = non_centers[another_endpoint_index].first;
                double another_non_center_endpoint_to_center_distance = non_centers[another_endpoint_index].second;

                int i_j;
                if (one_non_center_endpoint_center_index > another_non_center_endpoint_center_index)
                {
                    int temp = another_non_center_endpoint_center_index;
                    another_non_center_endpoint_center_index = one_non_center_endpoint_center_index;
                    one_non_center_endpoint_center_index = temp;
                }
                hash_function_two_keys_to_one_key(centers.size(), centers_unordered_map[one_non_center_endpoint_center_index], centers_unordered_map[another_non_center_endpoint_center_index], i_j);
                // std::cout << "==both non-center== center i: " << one_non_center_endpoint_center_index << ", center j: " << another_non_center_endpoint_center_index << ", original i: " << one_endpoint_index << ", original j: " << another_endpoint_index << ", real distance: " << distance_interval[k][i].first;

                if (pairwise_distance_center_to_center.count(i_j) == 0)
                {
                    std::vector<double> current_center_to_other_center_distance_hierarchy_graph(centers.size(), INF);
                    hierarchy_graph.shortest_distance_Dijkstra(centers_unordered_map[one_non_center_endpoint_center_index], centers_unordered_map[another_non_center_endpoint_center_index], current_center_to_other_center_distance_hierarchy_graph);
                    double distance_on_hierarchy_graph = current_center_to_other_center_distance_hierarchy_graph[centers_unordered_map[another_non_center_endpoint_center_index]];
                    // std::cout << ", distance on hierarchy graph: " << distance_on_hierarchy_graph << std::endl;

                    if (distance_on_hierarchy_graph > (1 + epsilon) * distance_interval[k][i].first)
                    {
                        // std::cout << "^^^ added, the direct hierarchy graph doesn't exist" << std::endl;
                        graph.add_edge_and_geo_path_size_Dijkstra(one_endpoint_index, another_endpoint_index, distance_interval[k][i].first, distance_interval[k][i].second.first);
                        added_edge_in_greedy_spanner.push_back(distance_interval[k][i]);
                        hierarchy_graph.add_edge_Dijkstra(centers_unordered_map[one_non_center_endpoint_center_index], centers_unordered_map[another_non_center_endpoint_center_index], distance_interval[k][i].first + one_non_center_endpoint_to_center_distance + another_non_center_endpoint_to_center_distance);
                        pairwise_distance_center_to_center[i_j] = distance_interval[k][i].first + one_non_center_endpoint_to_center_distance + another_non_center_endpoint_to_center_distance;
                    }
                }
                else if (pairwise_distance_center_to_center.count(i_j) != 0 && pairwise_distance_center_to_center[i_j] > (1 + epsilon) * distance_interval[k][i].first)
                {
                    // std::cout << ", distance on hierarchy graph: " << pairwise_distance_center_to_center[i_j] << std::endl;
                    // std::cout << "^^^ added, the direct hierarchy graph exists, but longer than (1+e)*edge length" << std::endl;
                    graph.add_edge_and_geo_path_size_Dijkstra(one_endpoint_index, another_endpoint_index, distance_interval[k][i].first, distance_interval[k][i].second.first);
                    added_edge_in_greedy_spanner.push_back(distance_interval[k][i]);
                    hierarchy_graph.add_edge_Dijkstra(centers_unordered_map[one_non_center_endpoint_center_index], centers_unordered_map[another_non_center_endpoint_center_index], distance_interval[k][i].first + one_non_center_endpoint_to_center_distance + another_non_center_endpoint_to_center_distance);
                    pairwise_distance_center_to_center[i_j] = distance_interval[k][i].first + one_non_center_endpoint_to_center_distance + another_non_center_endpoint_to_center_distance;
                }
                else if (pairwise_distance_center_to_center.count(i_j) != 0 && pairwise_distance_center_to_center[i_j] <= (1 + epsilon) * distance_interval[k][i].first)
                {
                    // std::cout << ", distance on hierarchy graph: " << pairwise_distance_center_to_center[i_j] << std::endl;
                    if (pairwise_distance_center_to_center[i_j] > distance_interval[k][i].first + one_non_center_endpoint_to_center_distance + another_non_center_endpoint_to_center_distance)
                    {
                        hierarchy_graph.add_edge_Dijkstra(centers_unordered_map[one_non_center_endpoint_center_index], centers_unordered_map[another_non_center_endpoint_center_index], distance_interval[k][i].first + one_non_center_endpoint_to_center_distance + another_non_center_endpoint_to_center_distance);
                        pairwise_distance_center_to_center[i_j] = distance_interval[k][i].first + one_non_center_endpoint_to_center_distance + another_non_center_endpoint_to_center_distance;
                    }
                }
            }
        }
        W = 2 * W;
    }

    pairwise_distance_poi_to_poi_hierarchy_greedy_spanner_map.clear();
    pairwise_path_poi_to_poi_hierarchy_greedy_spanner_map.clear();

    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> pairwise_path_poi_to_poi_map;
    pairwise_path_poi_to_poi_map.clear();
    for (int i = 0; i < pairwise_distance_poi_to_poi.size(); i++)
    {
        for (int j = 0; j < pairwise_distance_poi_to_poi[i].size(); j++)
        {
            int i_j;
            hash_function_two_keys_to_one_key(pairwise_distance_poi_to_poi.size(), i, i + j, i_j);
            pairwise_path_poi_to_poi_map[i_j] = pairwise_path_poi_to_poi[i][j];
            // std::cout << "i: " << i << ", j: " << i + j << ", distance: " << pairwise_distance_poi_to_poi[i][j] << std::endl;

            if (j == 0)
            {
                continue;
            }
        }
    }

    for (int i = 0; i < pairwise_distance_poi_to_poi.size(); i++)
    {
        std::vector<double> current_poi_to_other_poi_distance_hierarchy_greedy_spanner(pairwise_distance_poi_to_poi.size(), INF);
        std::vector<std::vector<int>> current_poi_to_other_poi_path_index_hierarchy_greedy_spanner(pairwise_distance_poi_to_poi.size());
        graph.shortest_path_Dijkstra(i, current_poi_to_other_poi_distance_hierarchy_greedy_spanner, current_poi_to_other_poi_path_index_hierarchy_greedy_spanner, INF);

        // for (int j = 0; j < current_poi_to_other_poi_distance_hierarchy_greedy_spanner.size(); j++)
        // {
        //     std::cout << "$$   i: " << i << ", j: " << j << ", approximate distance: " << current_poi_to_other_poi_distance_hierarchy_greedy_spanner[j] << std::endl;
        // }

        for (int j = i; j < pairwise_distance_poi_to_poi.size(); j++)
        {
            int i_j;
            hash_function_two_keys_to_one_key(pairwise_distance_poi_to_poi.size(), i, j, i_j);
            // pairwise_distance_poi_to_poi_hierarchy_greedy_spanner_map.insert(std::pair<int, double>(i_j, current_poi_to_other_poi_distance_greedy_spanner[j]));
            pairwise_distance_poi_to_poi_hierarchy_greedy_spanner_map[i_j] = current_poi_to_other_poi_distance_hierarchy_greedy_spanner[j];

            // if (current_poi_to_other_poi_distance_hierarchy_greedy_spanner[j] >= INF)
            // {
            //     std::cout << "j: " << j << ", approximate distance: " << current_poi_to_other_poi_distance_hierarchy_greedy_spanner[j] << std::endl;
            // }

            assert(current_poi_to_other_poi_distance_hierarchy_greedy_spanner[j] < INF);
            std::vector<geodesic::SurfacePoint> current_poi_to_other_poi_path_hierarchy_greedy_spanner;
            current_poi_to_other_poi_path_hierarchy_greedy_spanner.clear();
            for (int k = 0; k < current_poi_to_other_poi_path_index_hierarchy_greedy_spanner[j].size() - 1; k++)
            {
                int dest_index = current_poi_to_other_poi_path_index_hierarchy_greedy_spanner[j][k];
                int src_index = current_poi_to_other_poi_path_index_hierarchy_greedy_spanner[j][k + 1];
                int src_dest_index;
                bool reverse_path = false;
                if (src_index > dest_index)
                {
                    int temp = dest_index;
                    dest_index = src_index;
                    src_index = temp;
                    reverse_path = true;
                }
                hash_function_two_keys_to_one_key(pairwise_distance_poi_to_poi.size(), src_index, dest_index, src_dest_index);
                if (reverse_path)
                {
                    for (int m = pairwise_path_poi_to_poi_map[src_dest_index].size() - 1; m >= 0; m--)
                    {
                        current_poi_to_other_poi_path_hierarchy_greedy_spanner.push_back(pairwise_path_poi_to_poi_map[src_dest_index][m]);
                    }
                }
                else
                {
                    for (int m = 0; m < pairwise_path_poi_to_poi_map[src_dest_index].size(); m++)
                    {
                        current_poi_to_other_poi_path_hierarchy_greedy_spanner.push_back(pairwise_path_poi_to_poi_map[src_dest_index][m]);
                    }
                }
            }
            pairwise_path_poi_to_poi_hierarchy_greedy_spanner_map[i_j] = current_poi_to_other_poi_path_hierarchy_greedy_spanner;
            pairwise_path_poi_to_poi_size += current_poi_to_other_poi_path_hierarchy_greedy_spanner.size();
            // std::cout << "i: " << i << ", j: " << j << ", approximate distance: " << current_poi_to_other_poi_distance_greedy_spanner[j] << ", real distance: " << pairwise_distance_poi_to_poi[i][j - i] << ", ratio: " << current_poi_to_other_poi_distance_greedy_spanner[j] / pairwise_distance_poi_to_poi[i][j - i] << std::endl;
        }
    }
    hierarchy_greedy_spanner_edge_num = graph.get_edge_num_Dijkstra();
    hierarchy_greedy_spanner_weight = graph.get_total_weight_Dijkstra();
    hierarchy_greedy_spanner_size = hierarchy_greedy_spanner_edge_num * sizeof(double) + graph.get_total_geo_path_size_Dijkstra() * sizeof(geodesic::SurfacePoint);
    HGS_memory_usage += 0.5 * pairwise_distance_poi_to_poi.size() * (pairwise_distance_poi_to_poi.size() - 1) * sizeof(double) + pairwise_path_poi_to_poi_size * sizeof(geodesic::SurfacePoint);
    // std::cout << "Complete graph size: " << min_to_max_pairwise_distance_poi_to_poi.size() << ", Hierarchy greedy spanner size: " << graph.get_edge_num_Dijkstra() << ", Hierarchy greedy spanner total weight: " << graph.get_total_weight_Dijkstra() << std::endl;

    auto stop_HGS_time = std::chrono::high_resolution_clock::now();
    auto duration_HGS_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_HGS_time - start_HGS_time);
    HGS_time = duration_HGS_time.count();
    HGS_time /= 1000;
    // std::cout << "Total HGS time: " << HGS_time << " ms" << std::endl;
}

void complete_graph_query(int poi_num, std::unordered_map<int, double> &pairwise_distance_poi_to_poi_map,
                          std::unordered_map<int, std::vector<geodesic::SurfacePoint>> &pairwise_path_poi_to_poi_map,
                          int source_poi_index, int destination_poi_index, double &approximate_distance,
                          std::vector<geodesic::SurfacePoint> &approximate_path, double &query_time)
{
    auto start_query_time = std::chrono::high_resolution_clock::now();

    int x_y;
    if (source_poi_index > destination_poi_index)
    {
        int temp = destination_poi_index;
        destination_poi_index = source_poi_index;
        source_poi_index = temp;
    }
    hash_function_two_keys_to_one_key(poi_num, source_poi_index, destination_poi_index, x_y);
    approximate_distance = pairwise_distance_poi_to_poi_map[x_y];
    approximate_distance = round(approximate_distance * 1000000000.0) / 1000000000.0;

    approximate_path = pairwise_path_poi_to_poi_map[x_y];
    // double dist = 0;
    // for (int i = 0; i < pairwise_path_poi_to_poi_map[x_y].size() - 1; i++)
    // {
    //     dist += pairwise_path_poi_to_poi_map[x_y][i].distance(&pairwise_path_poi_to_poi_map[x_y][i + 1]);
    // }
    // std::cout << "path dist: " << dist << std::endl;

    auto stop_query_time = std::chrono::high_resolution_clock::now();
    auto duration_query_time = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_query_time - start_query_time);
    query_time = duration_query_time.count();
    query_time /= 1000000;
}

void pre_or_post_WSPD_oracle_preprocessing_and_query(
    int poi_num, geodesic::Mesh *mesh, std::vector<int> &poi_list, double epsilon,
    int source_poi_index, int destination_poi_index, double &preprocessing_time,
    double &query_time, double &memory_usage, double &WSPD_oracle_size, int &WSPD_oracle_edge_num,
    double &WSPD_oracle_weight, double &approximate_distance, std::vector<geodesic::SurfacePoint> &approximate_path)
{
    auto start_preprocessing_time = std::chrono::high_resolution_clock::now();

    geodesic::GeodesicAlgorithmExact algorithm(mesh);

    int geo_tree_node_id = 1; // the root node has 0 id
    std::unordered_map<int, GeoPair *> geopairs;
    geopairs.clear();
    std::vector<GeoNode *> all_poi;
    all_poi.clear();
    std::vector<std::pair<int, GeoNode *>> pois;
    pois.clear();
    std::unordered_map<int, int> poi_unordered_map;
    poi_unordered_map.clear();

    for (int i = 0; i < poi_num; i++)
    {
        GeoNode *n = new GeoNode(poi_list[i], 0);
        all_poi.push_back(n);
        std::pair<int, GeoNode *> m(poi_list[i], n);
        pois.push_back(m);
        poi_unordered_map[poi_list[i]] = i;
    }
    // std::cout << poi_list[0] << std::endl;
    // std::cout << "aa " << all_poi[0]->index << std::endl;
    // std::cout << "aa " << pois[0].first << std::endl;
    // printf("file read finished\n");

    double radius = 0;
    double distance;
    stx::btree<int, GeoNode *> pois_B_tree(pois.begin(), pois.end());
    // for (stx::btree<int, GeoNode *>::iterator ite = pois_B_tree.begin(); ite != pois_B_tree.end(); ite++)
    // {
    //     std::cout << ite.key() << " ";
    // }

    double const distance_limit = 0;
    std::vector<geodesic::SurfacePoint> one_source_poi_list;
    std::vector<geodesic::SurfacePoint> destinations_poi_list;
    one_source_poi_list.clear();
    destinations_poi_list.clear();
    one_source_poi_list.push_back(geodesic::SurfacePoint(&mesh->vertices()[poi_list[0]]));
    for (int i = 0; i < poi_num; i++)
    {
        destinations_poi_list.push_back(geodesic::SurfacePoint(&mesh->vertices()[poi_list[i]]));
    }

    algorithm.propagate(one_source_poi_list, distance_limit, &destinations_poi_list);
    for (int i = 0; i < poi_num; i++)
    {
        geodesic::SurfacePoint p(&mesh->vertices()[poi_list[i]]);
        algorithm.best_source(p, distance);
        radius = std::max(distance, radius);
    }
    // std::cout << "radius stopdis: " << radius << " " << algorithm.distance_stopped() << std::endl;
    GeoNode root_geo(0, poi_list[0], radius);
    // std::cout << root_geo.radius << std::endl;

    stx::btree<int, GeoNode *> pois_as_center_each_parent_layer;
    pois_as_center_each_parent_layer.clear();
    build_geo_tree(geo_tree_node_id, mesh, root_geo, poi_num, pois_B_tree, pois_as_center_each_parent_layer, algorithm);

    std::vector<GeoNode *> partition_tree_to_compressed_partition_tree_to_be_removed_nodes;
    partition_tree_to_compressed_partition_tree_to_be_removed_nodes.clear();
    std::unordered_map<int, GeoNode *> geo_node_in_partition_tree_unordered_map;
    geo_node_in_partition_tree_unordered_map.clear();
    partition_tree_to_compressed_partition_tree(root_geo, partition_tree_to_compressed_partition_tree_to_be_removed_nodes, geo_node_in_partition_tree_unordered_map);

    // print_geo_tree(root_geo);

    std::unordered_map<int, int> geo_pair_unordered_map;
    geo_pair_unordered_map.clear();
    std::unordered_map<int, double> pairwise_distance_unordered_map;
    pairwise_distance_unordered_map.clear();
    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> pairwise_path_unordered_map;
    pairwise_path_unordered_map.clear();
    int pairwise_path_poi_to_poi_size = 0;
    generate_geo_pair(geo_tree_node_id, WSPD_oracle_edge_num, WSPD_oracle_weight, mesh, root_geo, root_geo, algorithm, epsilon, geopairs, poi_unordered_map, geo_pair_unordered_map, pairwise_distance_unordered_map, pairwise_path_unordered_map, pairwise_path_poi_to_poi_size);
    WSPD_oracle_size = WSPD_oracle_edge_num * sizeof(double) + pairwise_path_poi_to_poi_size * sizeof(geodesic::SurfacePoint);
    // std::cout << "WSPD_oracle_edge_num: " << WSPD_oracle_edge_num << std::endl;
    // std::cout << "WSPD_oracle_weight: " << WSPD_oracle_weight << std::endl;

    memory_usage += algorithm.get_memory() + (geo_tree_node_id + 1) * sizeof(GeoNode) + WSPD_oracle_edge_num * sizeof(double) + pairwise_path_poi_to_poi_size * sizeof(geodesic::SurfacePoint);

    auto stop_preprocessing_time = std::chrono::high_resolution_clock::now();
    auto duration_preprocessing_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_preprocessing_time - start_preprocessing_time);
    preprocessing_time = duration_preprocessing_time.count();

    auto start_query_time = std::chrono::high_resolution_clock::now();

    approximate_distance = distance_query_geo(geo_tree_node_id, *geo_node_in_partition_tree_unordered_map[all_poi[source_poi_index]->index], *geo_node_in_partition_tree_unordered_map[all_poi[destination_poi_index]->index], geopairs, poi_unordered_map, approximate_path);
    approximate_distance = round(approximate_distance * 1000000000.0) / 1000000000.0;
    // std::cout << "i: " << all_poi[source_poi_index]->index << ", j: " << all_poi[destination_poi_index]->index << ", approximate_distance:" << approximate_distance << std::endl;

    // std::cout << "dist: " << length(approximate_path) << std::endl;

    auto stop_query_time = std::chrono::high_resolution_clock::now();
    auto duration_query_time = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_query_time - start_query_time);
    query_time = duration_query_time.count();
    query_time /= 1000000;

    // delete_geo_tree(root_geo);
}

void pre_or_post_SP_oracle_preprocessing_and_query(geodesic::Mesh *mesh,
                                                   std::vector<int> &poi_list, double epsilon,
                                                   int source_poi_index, int destination_poi_index,
                                                   double &preprocessing_time, double &query_time,
                                                   double &memory_usage, double &SP_oracle_size,
                                                   double &approximate_distance,
                                                   std::vector<geodesic::SurfacePoint> &approximate_path)
{
    auto start_preprocessing_time = std::chrono::high_resolution_clock::now();

    std::vector<geodesic::GeodesicAlgorithmSubdivision *> landmarks;
    // int landmark_size = floor(((double)mesh->vertices().size()) / (500.0 * (epsilon + 0.05)) + 100);
    int landmark_size = 500;
    // int landmark_size = 3;
    assert(landmark_size <= mesh->vertices().size());
    // std::cout << "landmark_size: " << landmark_size << std::endl;
    landmarks.resize(landmark_size);
    std::vector<int> x;
    x.resize(landmark_size);
    for (int i = 0; i < landmark_size; i++)
    {
        x[i] = rand() * rand() % mesh->vertices().size();
        for (int j = 0; j < i; j++)
        {
            if (i == 0)
                break;
            if (x[i] == x[j])
            {
                i--;
                break;
            }
        }
    }
    double subdivision_level = floor((std::pow(1.0 / (epsilon + 0.5), 1.5) + 1.0) * (11.1 * epsilon + 27.8));
    // double subdivision_level = floor((std::pow(1.0 / (epsilon + 0.05), 1.5) + 1.0));
    std::cout << "subdivision_level: " << subdivision_level << std::endl;
    for (int i = 0; i < landmark_size; i++)
    {
        landmarks[i] = new geodesic::GeodesicAlgorithmSubdivision(mesh, subdivision_level);
        geodesic::SurfacePoint p(&mesh->vertices()[x[i]]);
        std::vector<geodesic::SurfacePoint> sources;
        sources.clear();
        sources.push_back(p);
        landmarks[i]->propagate(sources);
        memory_usage += landmarks[i]->get_memory();
        SP_oracle_size += landmarks[i]->get_memory();
    }
    geodesic::GeodesicAlgorithmExact algorithm(mesh);

    auto stop_preprocessing_time = std::chrono::high_resolution_clock::now();
    auto duration_preprocessing_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_preprocessing_time - start_preprocessing_time);
    preprocessing_time = duration_preprocessing_time.count();

    auto start_query_time = std::chrono::high_resolution_clock::now();

    geodesic::SurfacePoint source(&mesh->vertices()[poi_list[source_poi_index]]);
    geodesic::SurfacePoint destination(&mesh->vertices()[poi_list[destination_poi_index]]);
    std::vector<geodesic::SurfacePoint> one_source_poi_list(1, source);
    std::vector<geodesic::SurfacePoint> one_destination_poi_list(1, destination);
    algorithm.propagate_landmark(one_source_poi_list, &landmarks, &one_destination_poi_list);
    // double dist;
    // algorithm.best_source(destination, dist);
    // std::cout << "dist: " << dist << std::endl;
    algorithm.trace_back(destination, approximate_path);
    approximate_distance = length(approximate_path);
    approximate_distance = round(approximate_distance * 1000000000.0) / 1000000000.0;

    memory_usage += algorithm.get_memory() + approximate_path.size() * sizeof(geodesic::SurfacePoint) + sizeof(double);
    // std::cout << "approximate_distance: " << approximate_distance << std::endl;

    auto stop_query_time = std::chrono::high_resolution_clock::now();
    auto duration_query_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_query_time - start_query_time);
    query_time = duration_query_time.count();
}

void pre_or_post_KF_query(geodesic::Mesh *mesh, std::vector<int> &poi_list,
                          double epsilon, int source_poi_index, int destination_poi_index,
                          double &query_time, double &memory_usage, double &approximate_distance,
                          std::vector<geodesic::SurfacePoint> &approximate_path)
{
    auto start_query_time = std::chrono::high_resolution_clock::now();

    double max_edge_length = -1e100;
    double min_edge_length = 1e100;
    std::vector<double> edge_length_list;
    std::vector<double> edge_length_without_outliers_list;
    for (unsigned i = 0; i < mesh->edges().size(); ++i)
    {
        geodesic::Edge &e = mesh->edges()[i];
        double edge_length = e.length();
        edge_length_list.push_back(edge_length);
    }
    remove_outliers(edge_length_list, edge_length_without_outliers_list);
    for (int i = 0; i < edge_length_without_outliers_list.size(); i++)
    {
        max_edge_length = std::max(max_edge_length, edge_length_without_outliers_list[i]);
        min_edge_length = std::min(min_edge_length, edge_length_without_outliers_list[i]);
    }
    // std::cout << "max_edge_length: " << max_edge_length << " , min_edge_length: " << min_edge_length << std::endl;
    double subdivision_level = floor(max_edge_length / min_edge_length * (1.0 / epsilon + 1.0) * 2.0);
    std::cout << "subdivision_level: " << subdivision_level << std::endl;
    geodesic::GeodesicAlgorithmSubdivision algorithm(mesh, subdivision_level);
    double const distance_limit = geodesic::GEODESIC_INF;
    geodesic::SurfacePoint source(&mesh->vertices()[poi_list[source_poi_index]]);
    geodesic::SurfacePoint destination(&mesh->vertices()[poi_list[destination_poi_index]]);
    std::vector<geodesic::SurfacePoint> one_source_poi_list(1, source);
    std::vector<geodesic::SurfacePoint> one_destination_poi_list(1, destination);
    algorithm.propagate(one_source_poi_list, distance_limit, &one_destination_poi_list);
    algorithm.trace_back(destination, approximate_path);
    approximate_distance = length(approximate_path);
    // algorithm.best_source(destination, approximate_distance);
    approximate_distance = round(approximate_distance * 1000000000.0) / 1000000000.0;
    memory_usage += algorithm.get_memory() + approximate_path.size() * sizeof(geodesic::SurfacePoint) + sizeof(double);

    auto stop_query_time = std::chrono::high_resolution_clock::now();
    auto duration_query_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_query_time - start_query_time);
    query_time = duration_query_time.count();
    query_time /= 1000;
}

void exact_distance(geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
                    geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list,
                    int source_poi_index, int destination_poi_index,
                    double &pre_exact_distance, double &post_exact_distance)
{
    pre_or_post_exact_distance(pre_mesh, pre_poi_list, source_poi_index,
                               destination_poi_index, pre_exact_distance);
    pre_or_post_exact_distance(post_mesh, post_poi_list, source_poi_index,
                               destination_poi_index, post_exact_distance);
}

void UE_N1(int poi_num, geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
           geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list,
           int source_poi_index, int destination_poi_index, double pre_exact_distance,
           double post_exact_distance, int &pre_MST_weight, int &post_MST_weight,
           std::string write_file_header)
{
    std::vector<std::vector<std::vector<int>>> pre_face_sequence_index_list;
    std::vector<std::vector<double>> pairwise_distance_poi_to_poi;
    std::vector<std::vector<bool>> pairwise_distance_poi_to_poi_changed;
    std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> pairwise_path_poi_to_poi;
    std::vector<std::vector<double>> pairwise_distance_poi_to_vertex(poi_num, std::vector<double>(pre_mesh->vertices().size(), 0));

    double pre_preprocessing_time = 0;
    double pre_hash_mapping_time = 0;
    double pre_query_time = 0;
    double pre_memory_usage = 0;
    double pre_hash_mapping_memory_usage = 0;
    double pre_index_size = 0;
    int pre_index_edge_num = 0;
    double pre_index_weight = 0;
    double pre_approximate_distance = 0;
    std::vector<geodesic::SurfacePoint> pre_approximate_path;
    pre_approximate_path.clear();

    pre_complete_graph_preprocessing(poi_num, pre_mesh, pre_poi_list, pre_face_sequence_index_list,
                                     pairwise_distance_poi_to_poi, pairwise_distance_poi_to_poi_changed,
                                     pairwise_path_poi_to_poi, pairwise_distance_poi_to_vertex,
                                     pre_preprocessing_time, pre_memory_usage);
    std::unordered_map<int, double> pairwise_distance_poi_to_poi_pre_map;
    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> pairwise_path_poi_to_poi_pre_map;
    get_pairwise_distance_and_path_poi_to_poi_map(pairwise_distance_poi_to_poi, pairwise_distance_poi_to_poi_pre_map,
                                                  pairwise_path_poi_to_poi, pairwise_path_poi_to_poi_pre_map,
                                                  pre_index_size, pre_index_edge_num, pre_index_weight, pre_hash_mapping_time,
                                                  pre_hash_mapping_memory_usage);
    complete_graph_query(poi_num, pairwise_distance_poi_to_poi_pre_map, pairwise_path_poi_to_poi_pre_map, source_poi_index,
                         destination_poi_index, pre_approximate_distance, pre_approximate_path, pre_query_time);
    calculate_MST_weight(pairwise_distance_poi_to_poi, pre_MST_weight);

    std::cout << "Pre terrain preprocessing time (CG): " << pre_preprocessing_time << " ms" << std::endl;
    std::cout << "Pre terrain preprocessing time (hash mapping): " << pre_hash_mapping_time << " ms" << std::endl;
    std::cout << "Pre terrain query time: " << pre_query_time << " ms" << std::endl;
    std::cout << "Pre terrain memory usage (CG): " << pre_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain memory usage (hash mapping): " << pre_hash_mapping_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index size: " << pre_index_size / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index edge number: " << pre_index_edge_num << " edges" << std::endl;
    std::cout << "Pre terrain index/MST weight: " << pre_index_weight / pre_MST_weight << std::endl;
    std::cout << "Pre terrain approximate distance: " << pre_approximate_distance << ", pre terrain exact distance: " << pre_exact_distance << ", distance error: " << pre_approximate_distance / pre_exact_distance - 1 << std::endl;
    std::cout << std::endl;

    std::ofstream ofs1("../output/output.txt", std::ios_base::app);
    ofs1 << "== UE_N1 ==\n";
    ofs1 << write_file_header << "\t"
         << pre_preprocessing_time << "\t"
         << pre_hash_mapping_time << "\t"
         << pre_query_time << "\t"
         << pre_memory_usage / 1e6 << "\t"
         << pre_hash_mapping_memory_usage / 1e6 << "\t"
         << pre_index_size / 1e6 << "\t"
         << pre_index_edge_num << "\t"
         << pre_index_weight / pre_MST_weight << "\t"
         << pre_approximate_distance / pre_exact_distance - 1 << "\t";
    ofs1.close();

    double post_updating_time = 0;
    double post_hash_mapping_time = 0;
    double post_query_time = 0;
    double post_memory_usage = 0;
    double post_hash_mapping_memory_usage = 0;
    double post_index_size = 0;
    int post_index_edge_num = 0;
    double post_index_weight = 0;
    double post_approximate_distance = 0;
    std::vector<geodesic::SurfacePoint> post_approximate_path;
    post_approximate_path.clear();

    post_complete_graph_updating(poi_num, pre_mesh, pre_poi_list, post_mesh, post_poi_list,
                                 pre_face_sequence_index_list, pairwise_distance_poi_to_poi,
                                 pairwise_distance_poi_to_poi_changed,
                                 pairwise_path_poi_to_poi, pairwise_distance_poi_to_vertex,
                                 post_updating_time, post_memory_usage);
    std::unordered_map<int, double> pairwise_distance_poi_to_poi_post_map;
    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> pairwise_path_poi_to_poi_post_map;
    get_pairwise_distance_and_path_poi_to_poi_map(pairwise_distance_poi_to_poi, pairwise_distance_poi_to_poi_post_map,
                                                  pairwise_path_poi_to_poi, pairwise_path_poi_to_poi_post_map,
                                                  post_index_size, post_index_edge_num, post_index_weight, post_hash_mapping_time,
                                                  post_hash_mapping_memory_usage);
    complete_graph_query(poi_num, pairwise_distance_poi_to_poi_post_map, pairwise_path_poi_to_poi_post_map, source_poi_index,
                         destination_poi_index, post_approximate_distance, post_approximate_path, post_query_time);
    calculate_MST_weight(pairwise_distance_poi_to_poi, post_MST_weight);

    std::cout << "Post terrain updating time (CG): " << post_updating_time << " ms" << std::endl;
    std::cout << "Post terrain preprocessing time (hash mapping): " << post_hash_mapping_time << " ms" << std::endl;
    std::cout << "Post terrain query time: " << post_query_time << " ms" << std::endl;
    std::cout << "Post terrain memory usage (CG): " << post_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain memory usage (hash mapping): " << post_hash_mapping_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index size: " << post_index_size / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index edge number: " << post_index_edge_num << " edges" << std::endl;
    std::cout << "Post terrain index/MST weight: " << post_index_weight / post_MST_weight << std::endl;
    std::cout << "Post terrain approximate distance: " << post_approximate_distance << ", post terrain exact distance: " << post_exact_distance << ", distance error: " << post_approximate_distance / post_exact_distance - 1 << std::endl;

    std::ofstream ofs2("../output/output.txt", std::ios_base::app);
    ofs2 << post_updating_time << "\t"
         << post_hash_mapping_time << "\t"
         << post_query_time << "\t"
         << post_memory_usage / 1e6 << "\t"
         << post_hash_mapping_memory_usage / 1e6 << "\t"
         << post_index_size / 1e6 << "\t"
         << post_index_edge_num << "\t"
         << post_index_weight / post_MST_weight << "\t"
         << post_approximate_distance / post_exact_distance - 1 << "\n\n";
    ofs2.close();
}

void UE_N2(int poi_num, geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
           geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list, double epsilon,
           int source_poi_index, int destination_poi_index, double pre_exact_distance,
           double post_exact_distance, int pre_MST_weight, int post_MST_weight,
           std::string write_file_header)
{
    std::vector<std::vector<std::vector<int>>> pre_face_sequence_index_list;
    std::vector<std::vector<double>> pairwise_distance_poi_to_poi;
    std::vector<std::vector<bool>> pairwise_distance_poi_to_poi_changed;
    std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> pairwise_path_poi_to_poi;
    std::vector<std::vector<double>> pairwise_distance_poi_to_vertex(poi_num, std::vector<double>(pre_mesh->vertices().size(), 0));

    double pre_preprocessing_time = 0;
    double pre_GS_time = 0;
    double pre_query_time = 0;
    double pre_memory_usage = 0;
    double pre_GS_memory_usage = 0;
    double pre_index_size = 0;
    int pre_index_edge_num = 0;
    double pre_index_weight = 0;
    double pre_approximate_distance = 0;
    std::vector<geodesic::SurfacePoint> pre_approximate_path;
    pre_approximate_path.clear();

    pre_complete_graph_preprocessing(poi_num, pre_mesh, pre_poi_list, pre_face_sequence_index_list,
                                     pairwise_distance_poi_to_poi, pairwise_distance_poi_to_poi_changed,
                                     pairwise_path_poi_to_poi, pairwise_distance_poi_to_vertex,
                                     pre_preprocessing_time, pre_memory_usage);
    std::unordered_map<int, double> pairwise_distance_poi_to_poi_pre_greedy_spanner_map;
    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> pairwise_path_poi_to_poi_pre_greedy_spanner_map;
    greedy_spanner(pairwise_distance_poi_to_poi, epsilon, pairwise_distance_poi_to_poi_pre_greedy_spanner_map,
                   pairwise_path_poi_to_poi, pairwise_path_poi_to_poi_pre_greedy_spanner_map,
                   pre_index_size, pre_index_edge_num, pre_index_weight, pre_GS_time, pre_GS_memory_usage);
    complete_graph_query(poi_num, pairwise_distance_poi_to_poi_pre_greedy_spanner_map, pairwise_path_poi_to_poi_pre_greedy_spanner_map,
                         source_poi_index, destination_poi_index, pre_approximate_distance, pre_approximate_path, pre_query_time);

    std::cout << "Pre terrain preprocessing time (CG): " << pre_preprocessing_time << " ms" << std::endl;
    std::cout << "Pre terrain preprocessing time (GS): " << pre_GS_time << " ms" << std::endl;
    std::cout << "Pre terrain query time: " << pre_query_time << " ms" << std::endl;
    std::cout << "Pre terrain memory usage (CG): " << pre_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain memory usage (GS): " << pre_GS_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index size: " << pre_index_size / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index edge number: " << pre_index_edge_num << " edges" << std::endl;
    std::cout << "Pre terrain index/MST weight: " << pre_index_weight / pre_MST_weight << std::endl;
    std::cout << "Pre terrain approximate distance: " << pre_approximate_distance << ", pre terrain exact distance: " << pre_exact_distance << ", distance error: " << pre_approximate_distance / pre_exact_distance - 1 << std::endl;
    std::cout << std::endl;

    std::ofstream ofs1("../output/output.txt", std::ios_base::app);
    ofs1 << "== UE_N2 ==\n";
    ofs1 << write_file_header << "\t"
         << pre_preprocessing_time << "\t"
         << pre_GS_time << "\t"
         << pre_query_time << "\t"
         << pre_memory_usage / 1e6 << "\t"
         << pre_GS_memory_usage / 1e6 << "\t"
         << pre_index_size / 1e6 << "\t"
         << pre_index_edge_num << "\t"
         << pre_index_weight / pre_MST_weight << "\t"
         << pre_approximate_distance / pre_exact_distance - 1 << "\t";
    ofs1.close();

    double post_updating_time = 0;
    double post_GS_time = 0;
    double post_query_time = 0;
    double post_memory_usage = 0;
    double post_GS_memory_usage = 0;
    double post_index_size = 0;
    int post_index_edge_num = 0;
    double post_index_weight = 0;
    double post_approximate_distance = 0;
    std::vector<geodesic::SurfacePoint> post_approximate_path;
    post_approximate_path.clear();

    post_complete_graph_updating(poi_num, pre_mesh, pre_poi_list, post_mesh, post_poi_list,
                                 pre_face_sequence_index_list, pairwise_distance_poi_to_poi,
                                 pairwise_distance_poi_to_poi_changed,
                                 pairwise_path_poi_to_poi, pairwise_distance_poi_to_vertex,
                                 post_updating_time, post_memory_usage);
    std::unordered_map<int, double> pairwise_distance_poi_to_poi_post_greedy_spanner_map;
    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> pairwise_path_poi_to_poi_post_greedy_spanner_map;
    greedy_spanner(pairwise_distance_poi_to_poi, epsilon, pairwise_distance_poi_to_poi_post_greedy_spanner_map,
                   pairwise_path_poi_to_poi, pairwise_path_poi_to_poi_post_greedy_spanner_map,
                   post_index_size, post_index_edge_num, post_index_weight, post_GS_time, post_GS_memory_usage);
    complete_graph_query(poi_num, pairwise_distance_poi_to_poi_post_greedy_spanner_map, pairwise_path_poi_to_poi_post_greedy_spanner_map,
                         source_poi_index, destination_poi_index, post_approximate_distance, post_approximate_path, post_query_time);

    std::cout << "Post terrain updating time (CG): " << post_updating_time << " ms" << std::endl;
    std::cout << "Post terrain updating time (GS): " << post_GS_time << " ms" << std::endl;
    std::cout << "Post terrain query time: " << post_query_time << " ms" << std::endl;
    std::cout << "Post terrain memory usage (CG): " << post_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain memory usage (GS): " << post_GS_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index size: " << post_index_size / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index edge number: " << post_index_edge_num << " edges" << std::endl;
    std::cout << "Post terrain index/MST weight: " << post_index_weight / post_MST_weight << std::endl;
    std::cout << "Post terrain approximate distance: " << post_approximate_distance << ", post terrain exact distance: " << post_exact_distance << ", distance error: " << post_approximate_distance / post_exact_distance - 1 << std::endl;

    std::ofstream ofs2("../output/output.txt", std::ios_base::app);
    ofs2 << post_updating_time << "\t"
         << post_GS_time << "\t"
         << post_query_time << "\t"
         << post_memory_usage / 1e6 << "\t"
         << post_GS_memory_usage / 1e6 << "\t"
         << post_index_size / 1e6 << "\t"
         << post_index_edge_num << "\t"
         << post_index_weight / post_MST_weight << "\t"
         << post_approximate_distance / post_exact_distance - 1 << "\n\n";
    ofs2.close();
}

void UE(int poi_num, geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
        geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list, double epsilon,
        int source_poi_index, int destination_poi_index, double pre_exact_distance,
        double post_exact_distance, int pre_MST_weight, int post_MST_weight,
        std::string write_file_header)
{
    std::vector<std::vector<std::vector<int>>> pre_face_sequence_index_list;
    std::vector<std::vector<double>> pairwise_distance_poi_to_poi;
    std::vector<std::vector<bool>> pairwise_distance_poi_to_poi_changed;
    std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> pairwise_path_poi_to_poi;
    std::vector<std::vector<double>> pairwise_distance_poi_to_vertex(poi_num, std::vector<double>(pre_mesh->vertices().size(), 0));

    double pre_preprocessing_time = 0;
    double pre_HGS_time = 0;
    double pre_query_time = 0;
    double pre_memory_usage = 0;
    double pre_HGS_memory_usage = 0;
    double pre_index_size = 0;
    int pre_index_edge_num = 0;
    double pre_index_weight = 0;
    double pre_approximate_distance = 0;
    std::vector<geodesic::SurfacePoint> pre_approximate_path;
    pre_approximate_path.clear();

    pre_complete_graph_preprocessing(poi_num, pre_mesh, pre_poi_list, pre_face_sequence_index_list,
                                     pairwise_distance_poi_to_poi, pairwise_distance_poi_to_poi_changed,
                                     pairwise_path_poi_to_poi, pairwise_distance_poi_to_vertex,
                                     pre_preprocessing_time, pre_memory_usage);
    std::unordered_map<int, double> pairwise_distance_poi_to_poi_pre_hierarchy_greedy_spanner_map;
    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> pairwise_path_poi_to_poi_pre_hierarchy_greedy_spanner_map;
    hierarchy_greedy_spanner(pairwise_distance_poi_to_poi, epsilon, pairwise_distance_poi_to_poi_pre_hierarchy_greedy_spanner_map,
                             pairwise_path_poi_to_poi, pairwise_path_poi_to_poi_pre_hierarchy_greedy_spanner_map,
                             pre_index_size, pre_index_edge_num, pre_index_weight, pre_HGS_time, pre_HGS_memory_usage);
    complete_graph_query(poi_num, pairwise_distance_poi_to_poi_pre_hierarchy_greedy_spanner_map, pairwise_path_poi_to_poi_pre_hierarchy_greedy_spanner_map,
                         source_poi_index, destination_poi_index, pre_approximate_distance, pre_approximate_path, pre_query_time);

    std::cout << "Pre terrain preprocessing time (CG): " << pre_preprocessing_time << " ms" << std::endl;
    std::cout << "Pre terrain preprocessing time (HGS): " << pre_HGS_time << " ms" << std::endl;
    std::cout << "Pre terrain query time: " << pre_query_time << " ms" << std::endl;
    std::cout << "Pre terrain memory usage (CG): " << pre_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain memory usage (HGS): " << pre_HGS_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index size: " << pre_index_size / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index edge number: " << pre_index_edge_num << " edges" << std::endl;
    std::cout << "Pre terrain index/MST weight: " << pre_index_weight / pre_MST_weight << std::endl;
    std::cout << "Pre terrain approximate distance: " << pre_approximate_distance << ", pre terrain exact distance: " << pre_exact_distance << ", distance error: " << pre_approximate_distance / pre_exact_distance - 1 << std::endl;
    std::cout << std::endl;

    std::ofstream ofs1("../output/output.txt", std::ios_base::app);
    ofs1 << "== UE ==\n";
    ofs1 << write_file_header << "\t"
         << pre_preprocessing_time << "\t"
         << pre_HGS_time << "\t"
         << pre_query_time << "\t"
         << pre_memory_usage / 1e6 << "\t"
         << pre_HGS_memory_usage / 1e6 << "\t"
         << pre_index_size / 1e6 << "\t"
         << pre_index_edge_num << "\t"
         << pre_index_weight / pre_MST_weight << "\t"
         << pre_approximate_distance / pre_exact_distance - 1 << "\t";
    ofs1.close();

    double post_updating_time = 0;
    double post_HGS_time = 0;
    double post_query_time = 0;
    double post_memory_usage = 0;
    double post_HGS_memory_usage = 0;
    double post_index_size = 0;
    int post_index_edge_num = 0;
    double post_index_weight = 0;
    double post_approximate_distance = 0;
    std::vector<geodesic::SurfacePoint> post_approximate_path;
    post_approximate_path.clear();

    post_complete_graph_updating(poi_num, pre_mesh, pre_poi_list, post_mesh, post_poi_list,
                                 pre_face_sequence_index_list, pairwise_distance_poi_to_poi,
                                 pairwise_distance_poi_to_poi_changed,
                                 pairwise_path_poi_to_poi, pairwise_distance_poi_to_vertex,
                                 post_updating_time, post_memory_usage);
    std::unordered_map<int, double> pairwise_distance_poi_to_poi_post_hierarchy_greedy_spanner_map;
    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> pairwise_path_poi_to_poi_post_hierarchy_greedy_spanner_map;
    hierarchy_greedy_spanner(pairwise_distance_poi_to_poi, epsilon, pairwise_distance_poi_to_poi_post_hierarchy_greedy_spanner_map,
                             pairwise_path_poi_to_poi, pairwise_path_poi_to_poi_post_hierarchy_greedy_spanner_map,
                             post_index_size, post_index_edge_num, post_index_weight, post_HGS_time, post_HGS_memory_usage);
    complete_graph_query(poi_num, pairwise_distance_poi_to_poi_post_hierarchy_greedy_spanner_map, pairwise_path_poi_to_poi_post_hierarchy_greedy_spanner_map,
                         source_poi_index, destination_poi_index, post_approximate_distance, post_approximate_path, post_query_time);

    std::cout << "Post terrain updating time (CG): " << post_updating_time << " ms" << std::endl;
    std::cout << "Post terrain updating time (HGS): " << post_HGS_time << " ms" << std::endl;
    std::cout << "Post terrain query time: " << post_query_time << " ms" << std::endl;
    std::cout << "Post terrain memory usage (CG): " << post_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain memory usage (HGS): " << post_HGS_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index size: " << post_index_size / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index edge number: " << post_index_edge_num << " edges" << std::endl;
    std::cout << "Post terrain index/MST weight: " << post_index_weight / post_MST_weight << std::endl;
    std::cout << "Post terrain approximate distance: " << post_approximate_distance << ", post terrain exact distance: " << post_exact_distance << ", distance error: " << post_approximate_distance / post_exact_distance - 1 << std::endl;

    std::ofstream ofs2("../output/output.txt", std::ios_base::app);
    ofs2 << post_updating_time << "\t"
         << post_HGS_time << "\t"
         << post_query_time << "\t"
         << post_memory_usage / 1e6 << "\t"
         << post_HGS_memory_usage / 1e6 << "\t"
         << post_index_size / 1e6 << "\t"
         << post_index_edge_num << "\t"
         << post_index_weight / post_MST_weight << "\t"
         << post_approximate_distance / post_exact_distance - 1 << "\n\n";
    ofs2.close();
}

void WSPD_oracle(int poi_num, geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
                 geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list, double epsilon,
                 int source_poi_index, int destination_poi_index, double pre_exact_distance,
                 double post_exact_distance, int pre_MST_weight, int post_MST_weight,
                 std::string write_file_header)
{
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
                                                    pre_memory_usage, pre_index_size, pre_index_edge_num,
                                                    pre_index_weight, pre_approximate_distance,
                                                    pre_approximate_path);

    std::cout << "Pre terrain preprocessing time: " << pre_preprocessing_time << " ms" << std::endl;
    std::cout << "Pre terrain query time: " << pre_query_time << " ms" << std::endl;
    std::cout << "Pre terrain memory usage: " << pre_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index size: " << pre_index_size / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index edge number: " << pre_index_edge_num << " edges" << std::endl;
    std::cout << "Pre terrain index/MST weight: " << pre_index_weight / pre_MST_weight << std::endl;
    std::cout << "Pre terrain approximate distance: " << pre_approximate_distance << ", pre terrain exact distance: " << pre_exact_distance << ", distance error: " << abs(pre_approximate_distance / pre_exact_distance - 1) << std::endl;
    std::cout << std::endl;

    std::ofstream ofs1("../output/output.txt", std::ios_base::app);
    ofs1 << "== WSPD_oracle ==\n";
    ofs1 << write_file_header << "\t"
         << pre_preprocessing_time << "\t"
         << 0 << "\t"
         << pre_query_time << "\t"
         << pre_memory_usage / 1e6 << "\t"
         << 0 << "\t"
         << pre_index_size / 1e6 << "\t"
         << pre_index_edge_num << "\t"
         << pre_index_weight / pre_MST_weight << "\t"
         << abs(pre_approximate_distance / pre_exact_distance - 1) << "\t";
    ofs1.close();

    double post_preprocessing_time = 0;
    double post_query_time = 0;
    double post_memory_usage = 0;
    double post_index_size = 0;
    int post_index_edge_num = 0;
    double post_index_weight = 0;
    double post_approximate_distance = 0;
    std::vector<geodesic::SurfacePoint> post_approximate_path;
    post_approximate_path.clear();

    pre_or_post_WSPD_oracle_preprocessing_and_query(poi_num, post_mesh, post_poi_list, epsilon,
                                                    source_poi_index, destination_poi_index,
                                                    post_preprocessing_time, post_query_time,
                                                    post_memory_usage, post_index_size, post_index_edge_num,
                                                    post_index_weight, post_approximate_distance,
                                                    post_approximate_path);

    std::cout << "Post terrain preprocessing time: " << post_preprocessing_time << " ms" << std::endl;
    std::cout << "Post terrain query time: " << post_query_time << " ms" << std::endl;
    std::cout << "Post terrain memory usage: " << post_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index size: " << pre_index_size / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index edge number: " << post_index_edge_num << " edges" << std::endl;
    std::cout << "Post terrain index/MST weight: " << post_index_weight / post_MST_weight << std::endl;
    std::cout << "Post terrain approximate distance: " << post_approximate_distance << ", post terrain exact distance: " << post_exact_distance << ", distance error: " << abs(post_approximate_distance / post_exact_distance - 1) << std::endl;

    std::ofstream ofs2("../output/output.txt", std::ios_base::app);
    ofs2 << post_preprocessing_time << "\t"
         << 0 << "\t"
         << post_query_time << "\t"
         << post_memory_usage / 1e6 << "\t"
         << 0 << "\t"
         << post_index_size / 1e6 << "\t"
         << post_index_edge_num << "\t"
         << post_index_weight / post_MST_weight << "\t"
         << abs(post_approximate_distance / post_exact_distance - 1) << "\n\n";
    ofs2.close();
}

void SP_oracle(geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
               geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list, double epsilon,
               int source_poi_index, int destination_poi_index, double pre_exact_distance,
               double post_exact_distance, std::string write_file_header)
{
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

    std::ofstream ofs1("../output/output.txt", std::ios_base::app);
    ofs1 << "== SP_oracle ==\n";
    ofs1 << write_file_header << "\t"
         << pre_preprocessing_time << "\t"
         << 0 << "\t"
         << pre_query_time << "\t"
         << pre_memory_usage / 1e6 << "\t"
         << 0 << "\t"
         << pre_index_size / 1e6 << "\t"
         << 0 << "\t"
         << 0 << "\t"
         << pre_approximate_distance / pre_exact_distance - 1 << "\t";
    ofs1.close();

    double post_preprocessing_time = 0;
    double post_query_time = 0;
    double post_memory_usage = 0;
    double post_index_size = 0;
    double post_approximate_distance = 0;
    std::vector<geodesic::SurfacePoint> post_approximate_path;
    post_approximate_path.clear();

    pre_or_post_SP_oracle_preprocessing_and_query(post_mesh, post_poi_list,
                                                  epsilon, source_poi_index, destination_poi_index,
                                                  post_preprocessing_time, post_query_time,
                                                  post_memory_usage, post_index_size, post_approximate_distance,
                                                  post_approximate_path);

    std::cout << "Post terrain preprocessing time: " << post_preprocessing_time << " ms" << std::endl;
    std::cout << "Post terrain query time: " << post_query_time << " ms" << std::endl;
    std::cout << "Post terrain memory usage: " << post_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index size: " << post_index_size / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain approximate distance: " << post_approximate_distance << ", post terrain exact distance: " << post_exact_distance << ", distance error: " << post_approximate_distance / post_exact_distance - 1 << std::endl;

    std::ofstream ofs2("../output/output.txt", std::ios_base::app);
    ofs2 << post_preprocessing_time << "\t"
         << 0 << "\t"
         << post_query_time << "\t"
         << post_memory_usage / 1e6 << "\t"
         << 0 << "\t"
         << post_index_size / 1e6 << "\t"
         << 0 << "\t"
         << 0 << "\t"
         << post_approximate_distance / post_exact_distance - 1 << "\n\n";
    ofs2.close();
}

void KF(geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
        geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list, double epsilon,
        int source_poi_index, int destination_poi_index, double pre_exact_distance,
        double post_exact_distance, std::string write_file_header)
{
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

    std::ofstream ofs1("../output/output.txt", std::ios_base::app);
    ofs1 << "== KF ==\n";
    ofs1 << write_file_header << "\t"
         << 0 << "\t"
         << 0 << "\t"
         << pre_query_time << "\t"
         << pre_memory_usage / 1e6 << "\t"
         << 0 << "\t"
         << 0 << "\t"
         << 0 << "\t"
         << 0 << "\t"
         << pre_approximate_distance / pre_exact_distance - 1 << "\t";
    ofs1.close();

    double post_query_time = 0;
    double post_memory_usage = 0;
    double post_approximate_distance = 0;
    std::vector<geodesic::SurfacePoint> post_approximate_path;
    post_approximate_path.clear();

    pre_or_post_KF_query(post_mesh, post_poi_list, epsilon, source_poi_index,
                         destination_poi_index, post_query_time, post_memory_usage,
                         post_approximate_distance, post_approximate_path);

    std::cout << "Post terrain query time: " << post_query_time << " ms" << std::endl;
    std::cout << "Post terrain memory usage: " << post_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain approximate distance: " << post_approximate_distance << ", post terrain exact distance: " << post_exact_distance << ", distance error: " << post_approximate_distance / post_exact_distance - 1 << std::endl;

    std::ofstream ofs2("../output/output.txt", std::ios_base::app);
    ofs2 << 0 << "\t"
         << 0 << "\t"
         << post_query_time << "\t"
         << post_memory_usage / 1e6 << "\t"
         << 0 << "\t"
         << 0 << "\t"
         << 0 << "\t"
         << 0 << "\t"
         << post_approximate_distance / post_exact_distance - 1 << "\n\n";
    ofs2.close();
}