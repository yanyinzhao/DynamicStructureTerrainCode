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

void exact_distance(geodesic::Mesh *mesh, std::vector<int> &poi_list,
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
}

void pre_complete_graph_construction(int poi_num, geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
                                     std::vector<std::vector<std::vector<int>>> &pre_face_sequence_index_list,
                                     std::vector<std::vector<double>> &pairwise_distance_poi_to_poi,
                                     std::vector<std::vector<bool>> &pairwise_distance_poi_to_poi_changed,
                                     std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> &pairwise_path_poi_to_poi,
                                     std::vector<std::vector<double>> &pairwise_distance_poi_to_vertex,
                                     double &construction_time, double &memory_usage)
{
    auto start_construction_time = std::chrono::high_resolution_clock::now();

    geodesic::GeodesicAlgorithmExact pre_algorithm(pre_mesh);

    int pre_face_sequence_index_list_size = 0;
    std::vector<std::vector<int>> one_poi_to_other_poi_pre_face_sequence_index_list;
    std::vector<int> one_poi_pre_face_sequence_index_list;
    pre_face_sequence_index_list.clear();

    double const distance_limit = geodesic::GEODESIC_INF;
    pairwise_distance_poi_to_poi.clear();
    pairwise_distance_poi_to_poi_changed.clear();
    pairwise_path_poi_to_poi.clear();
    int pairwise_path_poi_to_poi_size = 0;
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
        }

        for (int j = 0; j < pre_mesh->vertices().size(); j++)
        {
            destinations_vertex_list.push_back(geodesic::SurfacePoint(&pre_mesh->vertices()[j]));
        }

        pre_algorithm.propagate(one_source_poi_list, distance_limit, &destinations_vertex_list);

        for (int j = 0; j < i; j++)
        {
            one_poi_to_other_poi_pre_face_sequence_index_list.push_back(one_poi_pre_face_sequence_index_list);
        }

        for (int j = i; j < poi_num; j++)
        {
            std::vector<geodesic::SurfacePoint> path;
            pre_algorithm.trace_back(destinations_poi_list[j - i], path);
            current_poi_to_other_poi_distance.push_back(length(path));
            current_poi_to_other_poi_distance_changed.push_back(false);
            current_poi_to_other_poi_path.push_back(path);
            get_face_sequence(path, one_poi_pre_face_sequence_index_list);
            one_poi_to_other_poi_pre_face_sequence_index_list.push_back(one_poi_pre_face_sequence_index_list);
            pre_face_sequence_index_list_size += one_poi_pre_face_sequence_index_list.size();
            pairwise_path_poi_to_poi_size += path.size();
        }
        pairwise_distance_poi_to_poi.push_back(current_poi_to_other_poi_distance);
        pairwise_distance_poi_to_poi_changed.push_back(current_poi_to_other_poi_distance_changed);
        pairwise_path_poi_to_poi.push_back(current_poi_to_other_poi_path);
        pre_face_sequence_index_list.push_back(one_poi_to_other_poi_pre_face_sequence_index_list);

        for (int j = 0; j < pre_mesh->vertices().size(); j++)
        {
            double distance;
            pre_algorithm.best_source(destinations_vertex_list[j], distance);
            pairwise_distance_poi_to_vertex[i][j] = distance;
        }
    }
    memory_usage += pre_algorithm.get_memory() + pre_face_sequence_index_list_size * sizeof(int); // + 0.5 * poi_num * (poi_num + 1) * sizeof(double) + pairwise_path_poi_to_poi_size * sizeof(geodesic::SurfacePoint);

    auto stop_construction_time = std::chrono::high_resolution_clock::now();
    auto duration_construction_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_construction_time - start_construction_time);
    construction_time = duration_construction_time.count();
}

void post_complete_graph_update(int poi_num, geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
                                geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list,
                                std::vector<std::vector<std::vector<int>>> &pre_face_sequence_index_list,
                                std::vector<std::vector<double>> &pairwise_distance_poi_to_poi,
                                std::vector<std::vector<bool>> &pairwise_distance_poi_to_poi_changed,
                                std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> &pairwise_path_poi_to_poi,
                                std::vector<std::vector<double>> &pairwise_distance_poi_to_vertex,
                                double &update_time, double &memory_usage)
{
    auto start_update_time = std::chrono::high_resolution_clock::now();

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
            changed_face_index_list.push_back(i);
            changed_face_index_unordered_map[i] = i;
        }
    }
    assert(changed_face_index_list.size() == changed_face_index_unordered_map.size());

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

    // stores the changed poi index such that (1) the index of poi changed, (2) this poi is in the changed face, (3) the path that connects this poi passes the changed face
    std::vector<int> changed_poi_index_list(poi_num, 0);

    // stores only the poi in the changed face (used in calculation the 2D distance between these pois and other pois)
    std::vector<int> poi_in_the_changed_face_index_list;
    poi_in_the_changed_face_index_list.clear();

    // for the pre and post terrain, the index of poi may change, but the actual vertex on the terrain that these two pois stand for are very close
    // if the poi for the pre and post list changed, then the value in this list becomes 1, which indicates changes
    for (int i = 0; i < poi_num; i++)
    {
        if (pre_poi_list[i] != post_poi_list[i])
        {
            changed_poi_index_list[i] = 1;
        }
    }

    // for the pre and post terrain, even if the index of poi not changes, if the poi is on the changed face, we also indicate it by 3 in the following list
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
                changed_poi_index_list[i] = 3;
                poi_in_the_changed_face_index_list.push_back(i);
                break;
            }
        }
    }

    // if two pois are not in the changed face, but the path passes the changed face, we also indicate the poi by 4
    // note that for a poi, if one of the path that connects this poi passes the changed face, then we need to run SSAD
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
            }
            if (break_loop)
            {
                break;
            }
        }
    }

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
        }

        post_algorithm.propagate(one_source_poi_list, distance_limit, &destinations_poi_list);

        int index = 0;
        for (int j = 0; j < poi_num; j++)
        {
            if ((changed_poi_index_list[j] == 1 || changed_poi_index_list[j] == 3 || changed_poi_index_list[j] == 4) && j < i)
            {
                index++;
                continue;
            }
            std::vector<geodesic::SurfacePoint> path;
            post_algorithm.trace_back(destinations_poi_list[j - index], path);
            changed_pairwise_path_poi_to_poi_size += path.size();
            changed_pairwise_distance_poi_to_poi_size++;

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
    }

    // if two pois are not in the changed face, and the path doesn't pass the changed face, but one of pois is close to the
    // changed face, so we may need to update the new path with these two pois as endpoints on the new terrain, the following
    // is to check whether the original path is too close to the changed face or not, if so, we directly update the path
    // we also indicate this type of poi in changed_poi_index_list, but indicate it as 2 for clarification

    std::vector<double> euclidean_distance_of_poi_to_changed_face(poi_num, 0);
    std::vector<std::pair<double, int>> euclidean_distance_of_poi_to_changed_face_and_original_index;

    for (int i = 0; i < poi_num; i++)
    {
        if (changed_poi_index_list[i] == 1 || changed_poi_index_list[i] == 3 || changed_poi_index_list[i] == 4)
        {
            continue;
        }
        if (poi_in_the_changed_face_index_list.size() > 0)
        {
            for (int j = 0; j < poi_in_the_changed_face_index_list.size(); j++)
            {
                euclidean_distance_of_poi_to_changed_face[i] += euclidean_distance(post_mesh->vertices()[post_poi_list[i]].x(), post_mesh->vertices()[post_poi_list[i]].y(),
                                                                                   post_mesh->vertices()[poi_in_the_changed_face_index_list[j]].x(), post_mesh->vertices()[poi_in_the_changed_face_index_list[j]].y());
            }
        }
        else
        {
            euclidean_distance_of_poi_to_changed_face[i] = euclidean_distance(post_mesh->vertices()[post_poi_list[i]].x(), post_mesh->vertices()[post_poi_list[i]].y(),
                                                                              post_mesh->vertices()[changed_vertex_index_list[0]].x(), post_mesh->vertices()[changed_vertex_index_list[0]].y());
        }
    }

    sort_min_to_max_and_get_original_index(euclidean_distance_of_poi_to_changed_face, euclidean_distance_of_poi_to_changed_face_and_original_index);

    assert(euclidean_distance_of_poi_to_changed_face.size() == euclidean_distance_of_poi_to_changed_face_and_original_index.size());

    for (int i = 0; i < euclidean_distance_of_poi_to_changed_face_and_original_index.size(); i++)
    {
        if (euclidean_distance_of_poi_to_changed_face_and_original_index[i].first == 0)
        {
            continue;
        }
        int current_poi_index = euclidean_distance_of_poi_to_changed_face_and_original_index[i].second;

        assert(pairwise_distance_poi_to_poi[current_poi_index].size() == pairwise_distance_poi_to_poi_changed[current_poi_index].size());

        double max_distance = 0;
        for (int j = 0; j < euclidean_distance_of_poi_to_changed_face_and_original_index.size(); j++)
        {
            if (euclidean_distance_of_poi_to_changed_face_and_original_index[j].first == 0)
            {
                continue;
            }
            int checking_poi_index = euclidean_distance_of_poi_to_changed_face_and_original_index[j].second;

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

        // if the current poi is too close to the changed face, we need to run SSAD for this poi to update it path on the new terrain
        for (int k = 0; k < changed_vertex_index_list.size(); k++)
        {
            if (pairwise_distance_poi_to_vertex[current_poi_index][changed_vertex_index_list[k]] < 0)
            {
                changed_poi_index_list[current_poi_index] = 2;

                std::vector<int> destinations_poi_index_list;
                destinations_poi_index_list.clear();
                one_source_poi_list.clear();
                destinations_poi_list.clear();
                one_source_poi_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_poi_list[current_poi_index]]));

                for (int j = 0; j < euclidean_distance_of_poi_to_changed_face_and_original_index.size(); j++)
                {
                    int the_other_poi_index = euclidean_distance_of_poi_to_changed_face_and_original_index[j].second;

                    if ((current_poi_index <= the_other_poi_index && !pairwise_distance_poi_to_poi_changed[current_poi_index][the_other_poi_index - current_poi_index]) ||
                        (current_poi_index > the_other_poi_index && !pairwise_distance_poi_to_poi_changed[the_other_poi_index][current_poi_index - the_other_poi_index]))
                    {
                        destinations_poi_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_poi_list[the_other_poi_index]]));
                        destinations_poi_index_list.push_back(the_other_poi_index);
                    }
                }

                post_algorithm.propagate(one_source_poi_list, distance_limit, &destinations_poi_list);

                assert(destinations_poi_list.size() == destinations_poi_index_list.size());
                for (int j = 0; j < destinations_poi_index_list.size(); j++)
                {
                    int the_other_poi_index = destinations_poi_index_list[j];

                    std::vector<geodesic::SurfacePoint> path;
                    post_algorithm.trace_back(destinations_poi_list[j], path);
                    changed_pairwise_path_poi_to_poi_size += path.size();
                    changed_pairwise_distance_poi_to_poi_size++;

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
    memory_usage += post_algorithm.get_memory() + changed_face_index_list.size() * sizeof(int);

    auto stop_update_time = std::chrono::high_resolution_clock::now();
    auto duration_update_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_update_time - start_update_time);
    update_time = duration_update_time.count();
}

void post_complete_graph_update_RanUpdSeq_NoDistAppr(int poi_num, geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
                                                     geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list,
                                                     std::vector<std::vector<std::vector<int>>> &pre_face_sequence_index_list,
                                                     std::vector<std::vector<double>> &pairwise_distance_poi_to_poi,
                                                     std::vector<std::vector<bool>> &pairwise_distance_poi_to_poi_changed,
                                                     std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> &pairwise_path_poi_to_poi,
                                                     std::vector<std::vector<double>> &pairwise_distance_poi_to_vertex,
                                                     double &update_time, double &memory_usage)
{
    auto start_update_time = std::chrono::high_resolution_clock::now();

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
            changed_face_index_list.push_back(i);
            changed_face_index_unordered_map[i] = i;
        }
    }
    assert(changed_face_index_list.size() == changed_face_index_unordered_map.size());

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

    // stores only the poi in the changed face (used in calculation the 2D distance between these pois and other pois)
    std::vector<int> poi_in_the_changed_face_index_list;
    poi_in_the_changed_face_index_list.clear();

    // update the pairwise geodesic distance on post terrain for changed poi
    for (int i = 0; i < poi_num; i++)
    {
        one_source_poi_list.clear();
        destinations_poi_list.clear();
        one_source_poi_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_poi_list[i]]));

        for (int j = 0; j < poi_num; j++)
        {
            destinations_poi_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_poi_list[j]]));
        }
        post_algorithm.propagate(one_source_poi_list, distance_limit, &destinations_poi_list);

        for (int j = 0; j < poi_num; j++)
        {
            std::vector<geodesic::SurfacePoint> path;
            post_algorithm.trace_back(destinations_poi_list[j], path);
            changed_pairwise_path_poi_to_poi_size += path.size();
            changed_pairwise_distance_poi_to_poi_size++;

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
    }

    memory_usage += post_algorithm.get_memory() + changed_face_index_list.size() * sizeof(int);

    auto stop_update_time = std::chrono::high_resolution_clock::now();
    auto duration_update_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_update_time - start_update_time);
    update_time = duration_update_time.count();
}

void post_complete_graph_update_FullRad(int poi_num, geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
                                        geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list,
                                        std::vector<std::vector<std::vector<int>>> &pre_face_sequence_index_list,
                                        std::vector<std::vector<double>> &pairwise_distance_poi_to_poi,
                                        std::vector<std::vector<bool>> &pairwise_distance_poi_to_poi_changed,
                                        std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> &pairwise_path_poi_to_poi,
                                        std::vector<std::vector<double>> &pairwise_distance_poi_to_vertex,
                                        double &update_time, double &memory_usage)
{
    auto start_update_time = std::chrono::high_resolution_clock::now();

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

            changed_face_index_list.push_back(i);
            changed_face_index_unordered_map[i] = i;
        }
    }
    assert(changed_face_index_list.size() == changed_face_index_unordered_map.size());

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

    // stores the changed poi index such that (1) the index of poi changed, (2) this poi is in the changed face, (3) the path that connects this poi passes the changed face
    std::vector<int> changed_poi_index_list(poi_num, 0);

    // stores only the poi in the changed face (used in calculation the 2D distance between these pois and other pois)
    std::vector<int> poi_in_the_changed_face_index_list;
    poi_in_the_changed_face_index_list.clear();

    // for the pre and post terrain, the index of poi may change, but the actual vertex on the terrain that these two pois stand for are very close
    // if the poi for the pre and post list changed, then the value in this list becomes 1, which indicates changes
    for (int i = 0; i < poi_num; i++)
    {
        if (pre_poi_list[i] != post_poi_list[i])
        {
            changed_poi_index_list[i] = 1;
        }
    }

    // for the pre and post terrain, even if the index of poi not changes, if the poi is on the changed face, we also indicate it by 3 in the following list
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
                changed_poi_index_list[i] = 3;
                poi_in_the_changed_face_index_list.push_back(i);
                break;
            }
        }
    }

    // if two pois are not in the changed face, but the path passes the changed face, we also indicate the poi by 4
    // note that for a poi, if one of the path that connects this poi passes the changed face, then we need to run SSAD
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
            }
            if (break_loop)
            {
                break;
            }
        }
    }

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
        }
        post_algorithm.propagate(one_source_poi_list, distance_limit, &destinations_poi_list);

        int index = 0;
        for (int j = 0; j < poi_num; j++)
        {
            if ((changed_poi_index_list[j] == 1 || changed_poi_index_list[j] == 3 || changed_poi_index_list[j] == 4) && j < i)
            {
                index++;
                continue;
            }
            std::vector<geodesic::SurfacePoint> path;
            post_algorithm.trace_back(destinations_poi_list[j - index], path);
            changed_pairwise_path_poi_to_poi_size += path.size();
            changed_pairwise_distance_poi_to_poi_size++;

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
    }

    // if two pois are not in the changed face, and the path doesn't pass the changed face, but one of pois is close to the
    // changed face, so we may need to update the new path with these two pois as endpoints on the new terrain, the following
    // is to check whether the original path is too close to the changed face or not, if so, we directly update the path
    // we also indicate this type of poi in changed_poi_index_list, but indicate it as 2 for clarification

    std::vector<double> euclidean_distance_of_poi_to_changed_face(poi_num, 0);
    std::vector<std::pair<double, int>> euclidean_distance_of_poi_to_changed_face_and_original_index;

    for (int i = 0; i < poi_num; i++)
    {
        if (changed_poi_index_list[i] == 1 || changed_poi_index_list[i] == 3 || changed_poi_index_list[i] == 4)
        {
            continue;
        }
        if (poi_in_the_changed_face_index_list.size() > 0)
        {
            for (int j = 0; j < poi_in_the_changed_face_index_list.size(); j++)
            {
                euclidean_distance_of_poi_to_changed_face[i] += euclidean_distance(post_mesh->vertices()[post_poi_list[i]].x(), post_mesh->vertices()[post_poi_list[i]].y(),
                                                                                   post_mesh->vertices()[poi_in_the_changed_face_index_list[j]].x(), post_mesh->vertices()[poi_in_the_changed_face_index_list[j]].y());
            }
        }
        else
        {
            euclidean_distance_of_poi_to_changed_face[i] = euclidean_distance(post_mesh->vertices()[post_poi_list[i]].x(), post_mesh->vertices()[post_poi_list[i]].y(),
                                                                              post_mesh->vertices()[changed_vertex_index_list[0]].x(), post_mesh->vertices()[changed_vertex_index_list[0]].y());
        }
    }

    sort_min_to_max_and_get_original_index(euclidean_distance_of_poi_to_changed_face, euclidean_distance_of_poi_to_changed_face_and_original_index);

    assert(euclidean_distance_of_poi_to_changed_face.size() == euclidean_distance_of_poi_to_changed_face_and_original_index.size());

    for (int i = 0; i < euclidean_distance_of_poi_to_changed_face_and_original_index.size(); i++)
    {
        if (euclidean_distance_of_poi_to_changed_face_and_original_index[i].first == 0)
        {
            continue;
        }
        int current_poi_index = euclidean_distance_of_poi_to_changed_face_and_original_index[i].second;
        assert(pairwise_distance_poi_to_poi[current_poi_index].size() == pairwise_distance_poi_to_poi_changed[current_poi_index].size());

        double max_distance = 0;
        for (int j = 0; j < euclidean_distance_of_poi_to_changed_face_and_original_index.size(); j++)
        {
            if (euclidean_distance_of_poi_to_changed_face_and_original_index[j].first == 0)
            {
                continue;
            }
            int checking_poi_index = euclidean_distance_of_poi_to_changed_face_and_original_index[j].second;

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

        // if the current poi is too close to the changed face, we need to run SSAD for this poi to update it path on the new terrain
        for (int k = 0; k < changed_vertex_index_list.size(); k++)
        {
            if (pairwise_distance_poi_to_vertex[current_poi_index][changed_vertex_index_list[k]] < max_distance * 0.8)
            {
                changed_poi_index_list[current_poi_index] = 2;
                std::vector<int> destinations_poi_index_list;
                destinations_poi_index_list.clear();
                one_source_poi_list.clear();
                destinations_poi_list.clear();
                one_source_poi_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_poi_list[current_poi_index]]));

                for (int j = 0; j < euclidean_distance_of_poi_to_changed_face_and_original_index.size(); j++)
                {
                    int the_other_poi_index = euclidean_distance_of_poi_to_changed_face_and_original_index[j].second;

                    if ((current_poi_index <= the_other_poi_index && !pairwise_distance_poi_to_poi_changed[current_poi_index][the_other_poi_index - current_poi_index]) ||
                        (current_poi_index > the_other_poi_index && !pairwise_distance_poi_to_poi_changed[the_other_poi_index][current_poi_index - the_other_poi_index]))
                    {
                        destinations_poi_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_poi_list[the_other_poi_index]]));
                        destinations_poi_index_list.push_back(the_other_poi_index);
                    }
                }
                post_algorithm.propagate(one_source_poi_list, distance_limit, &destinations_poi_list);

                assert(destinations_poi_list.size() == destinations_poi_index_list.size());
                for (int j = 0; j < destinations_poi_index_list.size(); j++)
                {
                    int the_other_poi_index = destinations_poi_index_list[j];

                    std::vector<geodesic::SurfacePoint> path;
                    post_algorithm.trace_back(destinations_poi_list[j], path);
                    changed_pairwise_path_poi_to_poi_size += path.size();
                    changed_pairwise_distance_poi_to_poi_size++;

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
    memory_usage += post_algorithm.get_memory() + changed_face_index_list.size() * sizeof(int);

    auto stop_update_time = std::chrono::high_resolution_clock::now();
    auto duration_update_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_update_time - start_update_time);
    update_time = duration_update_time.count();
}

void post_complete_graph_update_NoEffIntChe(int poi_num, geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
                                            geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list,
                                            std::vector<std::vector<std::vector<int>>> &pre_face_sequence_index_list,
                                            std::vector<std::vector<double>> &pairwise_distance_poi_to_poi,
                                            std::vector<std::vector<bool>> &pairwise_distance_poi_to_poi_changed,
                                            std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> &pairwise_path_poi_to_poi,
                                            std::vector<std::vector<double>> &pairwise_distance_poi_to_vertex,
                                            double &update_time, double &memory_usage)
{
    auto start_update_time = std::chrono::high_resolution_clock::now();

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
            changed_face_index_list.push_back(i);
            changed_face_index_unordered_map[i] = i;
        }
    }
    assert(changed_face_index_list.size() == changed_face_index_unordered_map.size());

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

    // stores the changed poi index such that (1) the index of poi changed, (2) this poi is in the changed face, (3) the path that connects this poi passes the changed face
    std::vector<int> changed_poi_index_list(poi_num, 0);

    // stores only the poi in the changed face (used in calculation the 2D distance between these pois and other pois)
    std::vector<int> poi_in_the_changed_face_index_list;
    poi_in_the_changed_face_index_list.clear();

    // for the pre and post terrain, the index of poi may change, but the actual vertex on the terrain that these two pois stand for are very close
    // if the poi for the pre and post list changed, then the value in this list becomes 1, which indicates changes
    for (int i = 0; i < poi_num; i++)
    {
        if (pre_poi_list[i] != post_poi_list[i])
        {
            changed_poi_index_list[i] = 1;
        }
    }

    // for the pre and post terrain, even if the index of poi not changes, if the poi is on the changed face, we also indicate it by 3 in the following list
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
                changed_poi_index_list[i] = 3;
                poi_in_the_changed_face_index_list.push_back(i);
                break;
            }
        }
    }

    // if two pois are not in the changed face, but the path passes the changed face, we also indicate the poi by 4
    // note that for a poi, if one of the path that connects this poi passes the changed face, then we need to run SSAD
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
            }
            if (break_loop)
            {
                break;
            }
        }
    }

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
        }
        post_algorithm.propagate(one_source_poi_list, distance_limit, &destinations_poi_list);

        int index = 0;
        for (int j = 0; j < poi_num; j++)
        {
            if ((changed_poi_index_list[j] == 1 || changed_poi_index_list[j] == 3 || changed_poi_index_list[j] == 4) && j < i)
            {
                index++;
                continue;
            }
            std::vector<geodesic::SurfacePoint> path;
            post_algorithm.trace_back(destinations_poi_list[j - index], path);
            changed_pairwise_path_poi_to_poi_size += path.size();
            changed_pairwise_distance_poi_to_poi_size++;

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
    }

    // if two pois are not in the changed face, and the path doesn't pass the changed face, but one of pois is close to the
    // changed face, so we may need to update the new path with these two pois as endpoints on the new terrain, the following
    // is to check whether the original path is too close to the changed face or not, if so, we directly update the path
    // we also indicate this type of poi in changed_poi_index_list, but indicate it as 2 for clarification

    std::vector<double> euclidean_distance_of_poi_to_changed_face(poi_num, 0);
    std::vector<std::pair<double, int>> euclidean_distance_of_poi_to_changed_face_and_original_index;

    for (int i = 0; i < poi_num; i++)
    {
        if (changed_poi_index_list[i] == 1 || changed_poi_index_list[i] == 3 || changed_poi_index_list[i] == 4)
        {
            continue;
        }
        if (poi_in_the_changed_face_index_list.size() > 0)
        {
            for (int j = 0; j < poi_in_the_changed_face_index_list.size(); j++)
            {
                euclidean_distance_of_poi_to_changed_face[i] += euclidean_distance(post_mesh->vertices()[post_poi_list[i]].x(), post_mesh->vertices()[post_poi_list[i]].y(),
                                                                                   post_mesh->vertices()[poi_in_the_changed_face_index_list[j]].x(), post_mesh->vertices()[poi_in_the_changed_face_index_list[j]].y());
            }
        }
        else
        {
            euclidean_distance_of_poi_to_changed_face[i] = euclidean_distance(post_mesh->vertices()[post_poi_list[i]].x(), post_mesh->vertices()[post_poi_list[i]].y(),
                                                                              post_mesh->vertices()[changed_vertex_index_list[0]].x(), post_mesh->vertices()[changed_vertex_index_list[0]].y());
        }
    }
    sort_min_to_max_and_get_original_index(euclidean_distance_of_poi_to_changed_face, euclidean_distance_of_poi_to_changed_face_and_original_index);

    assert(euclidean_distance_of_poi_to_changed_face.size() == euclidean_distance_of_poi_to_changed_face_and_original_index.size());

    for (int i = 0; i < euclidean_distance_of_poi_to_changed_face_and_original_index.size(); i++)
    {
        if (euclidean_distance_of_poi_to_changed_face_and_original_index[i].first == 0)
        {
            continue;
        }
        int current_poi_index = euclidean_distance_of_poi_to_changed_face_and_original_index[i].second;

        assert(pairwise_distance_poi_to_poi[current_poi_index].size() == pairwise_distance_poi_to_poi_changed[current_poi_index].size());

        double max_distance = 0;
        for (int j = 0; j < euclidean_distance_of_poi_to_changed_face_and_original_index.size(); j++)
        {
            if (euclidean_distance_of_poi_to_changed_face_and_original_index[j].first == 0)
            {
                continue;
            }
            int checking_poi_index = euclidean_distance_of_poi_to_changed_face_and_original_index[j].second;

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

        // if the current poi is too close to the changed face, we need to run SSAD for this poi to update it path on the new terrain
        for (int k = 0; k < changed_vertex_index_list.size(); k++)
        {
            for (int t = 0; t < poi_num / 2; t++)
            {
                for (int m = 0; m < euclidean_distance_of_poi_to_changed_face_and_original_index.size(); m++)
                {
                    int the_other_poi_index = euclidean_distance_of_poi_to_changed_face_and_original_index[m].second;
                    if ((current_poi_index <= the_other_poi_index && !pairwise_distance_poi_to_poi_changed[current_poi_index][the_other_poi_index - current_poi_index]) ||
                        (current_poi_index > the_other_poi_index && !pairwise_distance_poi_to_poi_changed[the_other_poi_index][current_poi_index - the_other_poi_index]))
                    {
                        if (pairwise_distance_poi_to_vertex[current_poi_index][changed_vertex_index_list[k]] < 0 ||
                            pairwise_distance_poi_to_vertex[the_other_poi_index][changed_vertex_index_list[k]] < 0)
                        {
                            int counter = 0;
                            for (int n = 0; n < changed_vertex_index_list.size(); n++)
                            {
                                if (pairwise_distance_poi_to_vertex[current_poi_index][changed_vertex_index_list[n]] < 0)
                                {
                                    counter++;
                                }
                                else if (pairwise_distance_poi_to_vertex[the_other_poi_index][changed_vertex_index_list[n]] < 0)
                                {
                                    counter++;
                                }
                            }
                            changed_poi_index_list[current_poi_index] = 2;
                            break;
                        }
                    }
                }
            }

            if (changed_poi_index_list[current_poi_index] == 2)
            {
                std::vector<int> destinations_poi_index_list;
                destinations_poi_index_list.clear();
                one_source_poi_list.clear();
                destinations_poi_list.clear();
                one_source_poi_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_poi_list[current_poi_index]]));

                for (int j = 0; j < euclidean_distance_of_poi_to_changed_face_and_original_index.size(); j++)
                {
                    int the_other_poi_index = euclidean_distance_of_poi_to_changed_face_and_original_index[j].second;

                    if ((current_poi_index <= the_other_poi_index && !pairwise_distance_poi_to_poi_changed[current_poi_index][the_other_poi_index - current_poi_index]) ||
                        (current_poi_index > the_other_poi_index && !pairwise_distance_poi_to_poi_changed[the_other_poi_index][current_poi_index - the_other_poi_index]))
                    {
                        destinations_poi_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_poi_list[the_other_poi_index]]));
                        destinations_poi_index_list.push_back(the_other_poi_index);
                    }
                }
                post_algorithm.propagate(one_source_poi_list, distance_limit, &destinations_poi_list);

                assert(destinations_poi_list.size() == destinations_poi_index_list.size());
                for (int j = 0; j < destinations_poi_index_list.size(); j++)
                {
                    int the_other_poi_index = destinations_poi_index_list[j];

                    std::vector<geodesic::SurfacePoint> path;
                    post_algorithm.trace_back(destinations_poi_list[j], path);
                    changed_pairwise_path_poi_to_poi_size += path.size();
                    changed_pairwise_distance_poi_to_poi_size++;

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
    memory_usage += post_algorithm.get_memory() + changed_face_index_list.size() * sizeof(int);

    auto stop_update_time = std::chrono::high_resolution_clock::now();
    auto duration_update_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_update_time - start_update_time);
    update_time = duration_update_time.count();
}

void get_pairwise_distance_and_path_poi_to_poi_map(std::vector<std::vector<double>> pairwise_distance_poi_to_poi,
                                                   std::unordered_map<int, double> &pairwise_distance_poi_to_poi_map,
                                                   std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> pairwise_path_poi_to_poi,
                                                   std::unordered_map<int, std::vector<geodesic::SurfacePoint>> &pairwise_path_poi_to_poi_map,
                                                   double &complete_graph_size, int &complete_graph_edge_num,
                                                   double &complete_graph_weight, double &hash_mapping_time)
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

    auto stop_hash_mapping_time = std::chrono::high_resolution_clock::now();
    auto duration_hash_mapping_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_hash_mapping_time - start_hash_mapping_time);
    hash_mapping_time = duration_hash_mapping_time.count();
    hash_mapping_time /= 1000;
}

void greedy_spanner(double epsilon, Graph &graph, std::vector<std::vector<double>> pairwise_distance_poi_to_poi,
                    std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> pairwise_path_poi_to_poi,
                    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> &path_poi_to_poi_map,
                    double &greedy_spanner_size, int &greedy_spanner_edge_num, double &greedy_spanner_weight,
                    double &GreSpan_time, double &complete_graph_size)
{
    auto start_GreSpan_time = std::chrono::high_resolution_clock::now();

    std::vector<std::pair<double, std::pair<int, std::pair<int, int>>>> min_to_max_pairwise_distance_poi_to_poi;
    min_to_max_pairwise_distance_poi_to_poi.clear();

    int pairwise_path_poi_to_poi_size = 0;
    int complete_graph_edge_num = 0;

    for (int i = 0; i < pairwise_distance_poi_to_poi.size(); i++)
    {
        for (int j = 0; j < pairwise_distance_poi_to_poi[i].size(); j++)
        {
            if (j == 0)
            {
                continue;
            }
            if (pairwise_distance_poi_to_poi[i][j] >= 0)
            {
                complete_graph_edge_num++;
                pairwise_path_poi_to_poi_size += pairwise_path_poi_to_poi[i][j].size();
                min_to_max_pairwise_distance_poi_to_poi.push_back(std::make_pair(pairwise_distance_poi_to_poi[i][j], std::make_pair(pairwise_path_poi_to_poi[i][j].size(), std::make_pair(i, i + j))));
            }
        }
    }
    std::sort(min_to_max_pairwise_distance_poi_to_poi.begin(), min_to_max_pairwise_distance_poi_to_poi.end());

    for (int i = 0; i < min_to_max_pairwise_distance_poi_to_poi.size(); i++)
    {
        std::vector<double> current_poi_to_other_poi_distance_greedy_spanner(pairwise_distance_poi_to_poi.size(), INF);
        std::vector<std::vector<int>> current_poi_to_other_poi_path_index_greedy_spanner(pairwise_distance_poi_to_poi.size());
        graph.shortest_distance_Dijkstra(min_to_max_pairwise_distance_poi_to_poi[i].second.second.first, current_poi_to_other_poi_distance_greedy_spanner);
        double distance_on_graph = current_poi_to_other_poi_distance_greedy_spanner[min_to_max_pairwise_distance_poi_to_poi[i].second.second.second];

        if (distance_on_graph > (1 + epsilon) * min_to_max_pairwise_distance_poi_to_poi[i].first)
        {
            graph.add_edge_and_geo_path_size_Dijkstra(min_to_max_pairwise_distance_poi_to_poi[i].second.second.first, min_to_max_pairwise_distance_poi_to_poi[i].second.second.second, min_to_max_pairwise_distance_poi_to_poi[i].first, min_to_max_pairwise_distance_poi_to_poi[i].second.first);
        }
    }

    for (int i = 0; i < pairwise_distance_poi_to_poi.size(); i++)
    {
        for (int j = 0; j < pairwise_distance_poi_to_poi[i].size(); j++)
        {
            if (pairwise_distance_poi_to_poi[i][j] >= 0)
            {
                int i_j;
                hash_function_two_keys_to_one_key(pairwise_distance_poi_to_poi.size(), i, i + j, i_j);
                path_poi_to_poi_map[i_j] = pairwise_path_poi_to_poi[i][j];
            }
            if (j == 0)
            {
                continue;
            }
        }
    }

    greedy_spanner_edge_num = graph.get_edge_num_Dijkstra();
    greedy_spanner_weight = graph.get_total_weight_Dijkstra();
    greedy_spanner_size = greedy_spanner_edge_num * sizeof(double) + graph.get_total_geo_path_size_Dijkstra() * sizeof(geodesic::SurfacePoint);
    complete_graph_size = complete_graph_edge_num * sizeof(double) + pairwise_path_poi_to_poi_size * sizeof(geodesic::SurfacePoint);

    auto stop_GreSpan_time = std::chrono::high_resolution_clock::now();
    auto duration_GreSpan_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_GreSpan_time - start_GreSpan_time);
    GreSpan_time = duration_GreSpan_time.count();
    GreSpan_time /= 1000;
}

void hierarchy_greedy_spanner(double epsilon, Graph &graph, std::vector<std::vector<double>> pairwise_distance_poi_to_poi,
                              std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> pairwise_path_poi_to_poi,
                              std::unordered_map<int, std::vector<geodesic::SurfacePoint>> &path_poi_to_poi_map,
                              double &hierarchy_greedy_spanner_size, int &hierarchy_greedy_spanner_edge_num,
                              double &hierarchy_greedy_spanner_weight, double &HieGreSpan_time, double &complete_graph_size)
{
    auto start_HieGreSpan_time = std::chrono::high_resolution_clock::now();

    std::vector<std::pair<double, std::pair<int, std::pair<int, int>>>> min_to_max_pairwise_distance_poi_to_poi;
    min_to_max_pairwise_distance_poi_to_poi.clear();

    int pairwise_path_poi_to_poi_size = 0;
    int complete_graph_edge_num = 0;

    for (int i = 0; i < pairwise_distance_poi_to_poi.size(); i++)
    {
        for (int j = 0; j < pairwise_distance_poi_to_poi[i].size(); j++)
        {
            if (j == 0)
            {
                continue;
            }
            if (pairwise_distance_poi_to_poi[i][j] >= 0)
            {
                complete_graph_edge_num++;
                pairwise_path_poi_to_poi_size += pairwise_path_poi_to_poi[i][j].size();
                min_to_max_pairwise_distance_poi_to_poi.push_back(std::make_pair(pairwise_distance_poi_to_poi[i][j], std::make_pair(pairwise_path_poi_to_poi[i][j].size(), std::make_pair(i, i + j))));
            }
        }
    }
    std::sort(min_to_max_pairwise_distance_poi_to_poi.begin(), min_to_max_pairwise_distance_poi_to_poi.end());

    double max_pairwise_distance = min_to_max_pairwise_distance_poi_to_poi[min_to_max_pairwise_distance_poi_to_poi.size() - 1].first;

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

    for (int i = interval_zero_item_num; i < min_to_max_pairwise_distance_poi_to_poi.size(); i++)
    {
        if (min_to_max_pairwise_distance_poi_to_poi[i].first > pow(2, (interval_index - 1)) * max_pairwise_distance / pairwise_distance_poi_to_poi.size() &&
            min_to_max_pairwise_distance_poi_to_poi[i].first <= pow(2, interval_index) * max_pairwise_distance / pairwise_distance_poi_to_poi.size())
        {
            one_distance_interval.push_back(min_to_max_pairwise_distance_poi_to_poi[i]);
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

    double W = max_pairwise_distance / pairwise_distance_poi_to_poi.size();
    double delta = 0.5 * ((sqrt(epsilon + 1) - 1) / (sqrt(epsilon + 1) + 3)) * (log(epsilon / 4) / log(0.875));

    for (int k = 1; k <= ceil(log2(pairwise_distance_poi_to_poi.size())); k++)
    {
        // hierarchy graph
        std::vector<int> unprocessed_poi;
        unprocessed_poi.clear();

        for (int i = 0; i < pairwise_distance_poi_to_poi.size(); i++)
        {
            unprocessed_poi.push_back(i);
        }

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
                    non_centers[*ite2] = std::make_pair(current_index, current_poi_to_other_poi_distance_greedy_spanner[*ite2]);
                    unprocessed_poi.erase(ite2);
                }
                else
                {
                    ite2++;
                }
            }
        }

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
                if (current_center_to_other_center_and_non_center_distance_greedy_spanner[centers[j]] < W &&
                    current_center_to_other_center_and_non_center_distance_greedy_spanner[centers[j]] > 0)
                {
                    int i_j;
                    hash_function_two_keys_to_one_key(centers.size(), i, j, i_j);
                    pairwise_distance_center_to_center[i_j] = current_center_to_other_center_and_non_center_distance_greedy_spanner[centers[j]];
                    hierarchy_graph.add_edge_Dijkstra(i, j, current_center_to_other_center_and_non_center_distance_greedy_spanner[centers[j]]);
                }
                else if (current_center_to_other_center_and_non_center_distance_greedy_spanner[centers[j]] < W + 2 * W * delta &&
                         current_center_to_other_center_and_non_center_distance_greedy_spanner[centers[j]] >= W)
                {
                    int i_j;
                    hash_function_two_keys_to_one_key(centers.size(), i, j, i_j);
                    potential_second_type_inter_edge[i_j] = current_center_to_other_center_and_non_center_distance_greedy_spanner[centers[j]];
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
                        pairwise_distance_center_to_center[i_j] = added_edge_in_greedy_spanner[i].first;
                        hierarchy_graph.add_edge_Dijkstra(centers_unordered_map[one_endpoint_index], centers_unordered_map[another_endpoint_index], added_edge_in_greedy_spanner[i].first);
                    }
                    else if (pairwise_distance_center_to_center.count(i_j) != 0 && pairwise_distance_center_to_center[i_j] > added_edge_in_greedy_spanner[i].first)
                    {
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
                        pairwise_distance_center_to_center[i_j] = added_edge_in_greedy_spanner[i].first + one_non_center_endpoint_to_center_distance;
                        hierarchy_graph.add_edge_Dijkstra(centers_unordered_map[one_non_center_endpoint_center_index], centers_unordered_map[another_endpoint_index], added_edge_in_greedy_spanner[i].first + one_non_center_endpoint_to_center_distance);
                    }
                    else if (pairwise_distance_center_to_center.count(i_j) != 0 && pairwise_distance_center_to_center[i_j] > added_edge_in_greedy_spanner[i].first + one_non_center_endpoint_to_center_distance)
                    {
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
                        pairwise_distance_center_to_center[i_j] = added_edge_in_greedy_spanner[i].first + another_non_center_endpoint_to_center_distance;
                        hierarchy_graph.add_edge_Dijkstra(centers_unordered_map[one_endpoint_index], centers_unordered_map[another_non_center_endpoint_center_index], added_edge_in_greedy_spanner[i].first + another_non_center_endpoint_to_center_distance);
                    }
                    else if (pairwise_distance_center_to_center.count(i_j) != 0 && pairwise_distance_center_to_center[i_j] > added_edge_in_greedy_spanner[i].first + another_non_center_endpoint_to_center_distance)
                    {
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
                        pairwise_distance_center_to_center[i_j] = added_edge_in_greedy_spanner[i].first + one_non_center_endpoint_to_center_distance + another_non_center_endpoint_to_center_distance;
                        hierarchy_graph.add_edge_Dijkstra(centers_unordered_map[one_non_center_endpoint_center_index], centers_unordered_map[another_non_center_endpoint_center_index], added_edge_in_greedy_spanner[i].first + one_non_center_endpoint_to_center_distance + another_non_center_endpoint_to_center_distance);
                    }
                    else if (pairwise_distance_center_to_center.count(i_j) != 0 && pairwise_distance_center_to_center[i_j] > added_edge_in_greedy_spanner[i].first + one_non_center_endpoint_to_center_distance + another_non_center_endpoint_to_center_distance)
                    {
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

                if (pairwise_distance_center_to_center.count(i_j) == 0)
                {
                    std::vector<double> current_center_to_other_center_distance_hierarchy_graph(centers.size(), INF);
                    hierarchy_graph.shortest_distance_Dijkstra(centers_unordered_map[one_endpoint_index], centers_unordered_map[another_endpoint_index], current_center_to_other_center_distance_hierarchy_graph);
                    double distance_on_hierarchy_graph = current_center_to_other_center_distance_hierarchy_graph[centers_unordered_map[another_endpoint_index]];

                    if (distance_on_hierarchy_graph > (1 + epsilon) * distance_interval[k][i].first)
                    {
                        graph.add_edge_and_geo_path_size_Dijkstra(one_endpoint_index, another_endpoint_index, distance_interval[k][i].first, distance_interval[k][i].second.first);
                        added_edge_in_greedy_spanner.push_back(distance_interval[k][i]);
                        hierarchy_graph.add_edge_Dijkstra(centers_unordered_map[one_endpoint_index], centers_unordered_map[another_endpoint_index], distance_interval[k][i].first);
                        pairwise_distance_center_to_center[i_j] = distance_interval[k][i].first;
                    }
                }
                else if (pairwise_distance_center_to_center.count(i_j) != 0 && pairwise_distance_center_to_center[i_j] > (1 + epsilon) * distance_interval[k][i].first)
                {
                    graph.add_edge_and_geo_path_size_Dijkstra(one_endpoint_index, another_endpoint_index, distance_interval[k][i].first, distance_interval[k][i].second.first);
                    added_edge_in_greedy_spanner.push_back(distance_interval[k][i]);
                    hierarchy_graph.update_edge_Dijkstra(centers_unordered_map[one_endpoint_index], centers_unordered_map[another_endpoint_index], distance_interval[k][i].first);
                    pairwise_distance_center_to_center[i_j] = distance_interval[k][i].first;
                }
                else if (pairwise_distance_center_to_center.count(i_j) != 0 && pairwise_distance_center_to_center[i_j] <= (1 + epsilon) * distance_interval[k][i].first)
                {
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
                if (pairwise_distance_center_to_center.count(i_j) == 0)
                {
                    std::vector<double> current_center_to_other_center_distance_hierarchy_graph(centers.size(), INF);
                    hierarchy_graph.shortest_distance_Dijkstra(centers_unordered_map[one_non_center_endpoint_center_index], centers_unordered_map[another_endpoint_index], current_center_to_other_center_distance_hierarchy_graph);
                    double distance_on_hierarchy_graph = current_center_to_other_center_distance_hierarchy_graph[centers_unordered_map[another_endpoint_index]];

                    if (distance_on_hierarchy_graph > (1 + epsilon) * distance_interval[k][i].first)
                    {
                        graph.add_edge_and_geo_path_size_Dijkstra(one_endpoint_index, another_endpoint_index, distance_interval[k][i].first, distance_interval[k][i].second.first);
                        added_edge_in_greedy_spanner.push_back(distance_interval[k][i]);
                        hierarchy_graph.add_edge_Dijkstra(centers_unordered_map[one_non_center_endpoint_center_index], centers_unordered_map[another_endpoint_index], distance_interval[k][i].first + one_non_center_endpoint_to_center_distance);
                        pairwise_distance_center_to_center[i_j] = distance_interval[k][i].first + one_non_center_endpoint_to_center_distance;
                    }
                }
                else if (pairwise_distance_center_to_center.count(i_j) != 0 && pairwise_distance_center_to_center[i_j] > (1 + epsilon) * distance_interval[k][i].first)
                {
                    graph.add_edge_and_geo_path_size_Dijkstra(one_endpoint_index, another_endpoint_index, distance_interval[k][i].first, distance_interval[k][i].second.first);
                    added_edge_in_greedy_spanner.push_back(distance_interval[k][i]);
                    hierarchy_graph.update_edge_Dijkstra(centers_unordered_map[one_non_center_endpoint_center_index], centers_unordered_map[another_endpoint_index], distance_interval[k][i].first + one_non_center_endpoint_to_center_distance);
                    pairwise_distance_center_to_center[i_j] = distance_interval[k][i].first + one_non_center_endpoint_to_center_distance;
                }
                else if (pairwise_distance_center_to_center.count(i_j) != 0 && pairwise_distance_center_to_center[i_j] <= (1 + epsilon) * distance_interval[k][i].first)
                {
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

                if (pairwise_distance_center_to_center.count(i_j) == 0)
                {
                    std::vector<double> current_center_to_other_center_distance_hierarchy_graph(centers.size(), INF);
                    hierarchy_graph.shortest_distance_Dijkstra(centers_unordered_map[one_endpoint_index], centers_unordered_map[another_non_center_endpoint_center_index], current_center_to_other_center_distance_hierarchy_graph);
                    double distance_on_hierarchy_graph = current_center_to_other_center_distance_hierarchy_graph[centers_unordered_map[another_non_center_endpoint_center_index]];

                    if (distance_on_hierarchy_graph > (1 + epsilon) * distance_interval[k][i].first)
                    {
                        graph.add_edge_and_geo_path_size_Dijkstra(one_endpoint_index, another_endpoint_index, distance_interval[k][i].first, distance_interval[k][i].second.first);
                        added_edge_in_greedy_spanner.push_back(distance_interval[k][i]);
                        hierarchy_graph.add_edge_Dijkstra(centers_unordered_map[one_endpoint_index], centers_unordered_map[another_non_center_endpoint_center_index], distance_interval[k][i].first + another_non_center_endpoint_to_center_distance);
                        pairwise_distance_center_to_center[i_j] = distance_interval[k][i].first + another_non_center_endpoint_to_center_distance;
                    }
                }
                else if (pairwise_distance_center_to_center.count(i_j) != 0 && pairwise_distance_center_to_center[i_j] > (1 + epsilon) * distance_interval[k][i].first)
                {
                    graph.add_edge_and_geo_path_size_Dijkstra(one_endpoint_index, another_endpoint_index, distance_interval[k][i].first, distance_interval[k][i].second.first);
                    added_edge_in_greedy_spanner.push_back(distance_interval[k][i]);
                    hierarchy_graph.update_edge_Dijkstra(centers_unordered_map[one_endpoint_index], centers_unordered_map[another_non_center_endpoint_center_index], distance_interval[k][i].first + another_non_center_endpoint_to_center_distance);
                    pairwise_distance_center_to_center[i_j] = distance_interval[k][i].first + another_non_center_endpoint_to_center_distance;
                }
                else if (pairwise_distance_center_to_center.count(i_j) != 0 && pairwise_distance_center_to_center[i_j] <= (1 + epsilon) * distance_interval[k][i].first)
                {
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

                if (pairwise_distance_center_to_center.count(i_j) == 0)
                {
                    std::vector<double> current_center_to_other_center_distance_hierarchy_graph(centers.size(), INF);
                    hierarchy_graph.shortest_distance_Dijkstra(centers_unordered_map[one_non_center_endpoint_center_index], centers_unordered_map[another_non_center_endpoint_center_index], current_center_to_other_center_distance_hierarchy_graph);
                    double distance_on_hierarchy_graph = current_center_to_other_center_distance_hierarchy_graph[centers_unordered_map[another_non_center_endpoint_center_index]];

                    if (distance_on_hierarchy_graph > (1 + epsilon) * distance_interval[k][i].first)
                    {
                        graph.add_edge_and_geo_path_size_Dijkstra(one_endpoint_index, another_endpoint_index, distance_interval[k][i].first, distance_interval[k][i].second.first);
                        added_edge_in_greedy_spanner.push_back(distance_interval[k][i]);
                        hierarchy_graph.add_edge_Dijkstra(centers_unordered_map[one_non_center_endpoint_center_index], centers_unordered_map[another_non_center_endpoint_center_index], distance_interval[k][i].first + one_non_center_endpoint_to_center_distance + another_non_center_endpoint_to_center_distance);
                        pairwise_distance_center_to_center[i_j] = distance_interval[k][i].first + one_non_center_endpoint_to_center_distance + another_non_center_endpoint_to_center_distance;
                    }
                }
                else if (pairwise_distance_center_to_center.count(i_j) != 0 && pairwise_distance_center_to_center[i_j] > (1 + epsilon) * distance_interval[k][i].first)
                {
                    graph.add_edge_and_geo_path_size_Dijkstra(one_endpoint_index, another_endpoint_index, distance_interval[k][i].first, distance_interval[k][i].second.first);
                    added_edge_in_greedy_spanner.push_back(distance_interval[k][i]);
                    hierarchy_graph.add_edge_Dijkstra(centers_unordered_map[one_non_center_endpoint_center_index], centers_unordered_map[another_non_center_endpoint_center_index], distance_interval[k][i].first + one_non_center_endpoint_to_center_distance + another_non_center_endpoint_to_center_distance);
                    pairwise_distance_center_to_center[i_j] = distance_interval[k][i].first + one_non_center_endpoint_to_center_distance + another_non_center_endpoint_to_center_distance;
                }
                else if (pairwise_distance_center_to_center.count(i_j) != 0 && pairwise_distance_center_to_center[i_j] <= (1 + epsilon) * distance_interval[k][i].first)
                {
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

    for (int i = 0; i < pairwise_distance_poi_to_poi.size(); i++)
    {
        for (int j = 0; j < pairwise_distance_poi_to_poi[i].size(); j++)
        {
            if (pairwise_distance_poi_to_poi[i][j] >= 0)
            {
                int i_j;
                hash_function_two_keys_to_one_key(pairwise_distance_poi_to_poi.size(), i, i + j, i_j);
                path_poi_to_poi_map[i_j] = pairwise_path_poi_to_poi[i][j];
            }
            if (j == 0)
            {
                continue;
            }
        }
    }

    hierarchy_greedy_spanner_edge_num = graph.get_edge_num_Dijkstra();
    hierarchy_greedy_spanner_weight = graph.get_total_weight_Dijkstra();
    hierarchy_greedy_spanner_size = hierarchy_greedy_spanner_edge_num * sizeof(double) + graph.get_total_geo_path_size_Dijkstra() * sizeof(geodesic::SurfacePoint);
    complete_graph_size = complete_graph_edge_num * sizeof(double) + pairwise_path_poi_to_poi_size * sizeof(geodesic::SurfacePoint);

    auto stop_HieGreSpan_time = std::chrono::high_resolution_clock::now();
    auto duration_HieGreSpan_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_HieGreSpan_time - start_HieGreSpan_time);
    HieGreSpan_time = duration_HieGreSpan_time.count();
    HieGreSpan_time /= 1000;
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

    auto stop_query_time = std::chrono::high_resolution_clock::now();
    auto duration_query_time = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_query_time - start_query_time);
    query_time = duration_query_time.count();
    query_time /= 1000000;
}

void spanner_query(int poi_num, Graph graph,
                   std::unordered_map<int, std::vector<geodesic::SurfacePoint>> &path_poi_to_poi_map,
                   int source_poi_index, int destination_poi_index, double &approximate_distance,
                   std::vector<geodesic::SurfacePoint> &approximate_path, double &query_time)
{
    auto start_query_time = std::chrono::high_resolution_clock::now();

    if (source_poi_index > destination_poi_index)
    {
        int temp1 = destination_poi_index;
        destination_poi_index = source_poi_index;
        source_poi_index = temp1;
    }

    std::vector<double> current_poi_to_other_poi_distance_greedy_spanner(poi_num, INF);
    std::vector<std::vector<int>> current_poi_to_other_poi_path_index_greedy_spanner(poi_num);

    graph.shortest_path_Dijkstra(source_poi_index, destination_poi_index, current_poi_to_other_poi_distance_greedy_spanner, current_poi_to_other_poi_path_index_greedy_spanner);
    approximate_distance = current_poi_to_other_poi_distance_greedy_spanner[destination_poi_index];
    assert(current_poi_to_other_poi_distance_greedy_spanner[destination_poi_index] < INF);

    for (int k = 0; k < current_poi_to_other_poi_path_index_greedy_spanner[destination_poi_index].size() - 1; k++)
    {
        int dest_index = current_poi_to_other_poi_path_index_greedy_spanner[destination_poi_index][k];
        int src_index = current_poi_to_other_poi_path_index_greedy_spanner[destination_poi_index][k + 1];
        int src_dest_index;
        bool reverse_path = false;
        if (src_index > dest_index)
        {
            int temp2 = dest_index;
            dest_index = src_index;
            src_index = temp2;
            reverse_path = true;
        }

        hash_function_two_keys_to_one_key(poi_num, src_index, dest_index, src_dest_index);
        if (reverse_path)
        {
            for (int m = path_poi_to_poi_map[src_dest_index].size() - 1; m >= 0; m--)
            {
                approximate_path.push_back(path_poi_to_poi_map[src_dest_index][m]);
            }
        }
        else
        {
            for (int m = 0; m < path_poi_to_poi_map[src_dest_index].size(); m++)
            {
                approximate_path.push_back(path_poi_to_poi_map[src_dest_index][m]);
            }
        }
    }

    auto stop_query_time = std::chrono::high_resolution_clock::now();
    auto duration_query_time = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_query_time - start_query_time);
    query_time = duration_query_time.count();
    query_time /= 1000000;
}

void pre_or_post_WSPD_Oracle_and_pre_WSPD_Oracle_Adapt_construction(
    int poi_num, geodesic::Mesh *mesh, std::vector<int> &poi_list, double epsilon,
    int &geo_tree_node_id, std::vector<GeoNode *> &all_poi, std::unordered_map<int, GeoNode *> &geo_node_in_partition_tree_unordered_map,
    std::unordered_map<int, GeoPair *> &geopairs, std::unordered_map<int, int> &poi_unordered_map, std::vector<std::vector<double>> &pairwise_distance_poi_to_vertex,
    bool run_WSPD_Oracle_Adapt, double &construction_time, double &memory_usage, double &WSPD_Oracle_size, int &WSPD_Oracle_edge_num, double &WSPD_Oracle_weight)
{
    auto start_construction_time = std::chrono::high_resolution_clock::now();

    geodesic::GeodesicAlgorithmExact algorithm(mesh);

    std::unordered_map<int, double> pre_distance_poi_to_poi_map;
    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> pre_path_poi_to_poi_map;
    double pre_memory_usage = 0;

    pre_compute_WSPD_Oracle(poi_num, mesh, poi_list, pre_distance_poi_to_poi_map, pre_path_poi_to_poi_map, pre_memory_usage);

    geopairs.clear();
    all_poi.clear();
    std::vector<std::pair<int, GeoNode *>> pois;
    pois.clear();
    poi_unordered_map.clear();

    for (int i = 0; i < poi_num; i++)
    {
        GeoNode *n = new GeoNode(poi_list[i], 0);
        all_poi.push_back(n);
        std::pair<int, GeoNode *> m(poi_list[i], n);
        pois.push_back(m);
        poi_unordered_map[poi_list[i]] = i;
    }

    double radius = 0;
    stx::btree<int, GeoNode *> pois_B_tree(pois.begin(), pois.end());

    for (int i = 0; i < poi_num; i++)
    {
        int x_in_poi_list = 0;
        int y_in_poi_list = i;
        int x_y_in_poi_list;
        if (x_in_poi_list <= y_in_poi_list)
        {
            hash_function_two_keys_to_one_key(poi_num, x_in_poi_list, y_in_poi_list, x_y_in_poi_list);
        }
        else
        {
            hash_function_two_keys_to_one_key(poi_num, y_in_poi_list, x_in_poi_list, x_y_in_poi_list);
        }
        radius = std::max(pre_distance_poi_to_poi_map[x_y_in_poi_list], radius);
    }
    GeoNode *root_geo = new GeoNode(0, poi_list[0], radius);

    stx::btree<int, GeoNode *> pois_as_center_each_parent_layer;
    pois_as_center_each_parent_layer.clear();
    build_geo_tree(geo_tree_node_id, mesh, *root_geo, poi_num, pois_B_tree, pois_as_center_each_parent_layer, pre_distance_poi_to_poi_map, pre_path_poi_to_poi_map, poi_unordered_map);

    std::vector<GeoNode *> partition_tree_to_compressed_partition_tree_to_be_removed_nodes;
    partition_tree_to_compressed_partition_tree_to_be_removed_nodes.clear();
    geo_node_in_partition_tree_unordered_map.clear();
    partition_tree_to_compressed_partition_tree(*root_geo, partition_tree_to_compressed_partition_tree_to_be_removed_nodes, geo_node_in_partition_tree_unordered_map);

    std::unordered_map<int, int> geo_pair_unordered_map;
    geo_pair_unordered_map.clear();
    int pairwise_path_poi_to_poi_size = 0;
    int face_sequence_index_list_size = 0;
    generate_geo_pair(geo_tree_node_id, WSPD_Oracle_edge_num, WSPD_Oracle_weight, mesh, *root_geo, *root_geo, epsilon, geopairs, poi_unordered_map, geo_pair_unordered_map, pre_distance_poi_to_poi_map, pre_path_poi_to_poi_map, pairwise_path_poi_to_poi_size, face_sequence_index_list_size);
    WSPD_Oracle_size = WSPD_Oracle_edge_num * sizeof(double) + pairwise_path_poi_to_poi_size * sizeof(geodesic::SurfacePoint);

    memory_usage += pre_memory_usage + (geo_tree_node_id + 1) * sizeof(GeoNode) + WSPD_Oracle_edge_num * sizeof(double) + pairwise_path_poi_to_poi_size * sizeof(geodesic::SurfacePoint);

    if (run_WSPD_Oracle_Adapt)
    {
        std::vector<geodesic::SurfacePoint> one_source_poi_list2;

        std::vector<geodesic::SurfacePoint> destinations_vertex_list;
        for (int i = 0; i < poi_num; i++)
        {
            one_source_poi_list2.clear();
            one_source_poi_list2.push_back(geodesic::SurfacePoint(&mesh->vertices()[poi_list[i]]));

            for (int j = 0; j < mesh->vertices().size(); j++)
            {
                destinations_vertex_list.push_back(geodesic::SurfacePoint(&mesh->vertices()[j]));
            }
            algorithm.propagate(one_source_poi_list2, geodesic::GEODESIC_INF, &destinations_vertex_list);
            for (int j = 0; j < mesh->vertices().size(); j++)
            {
                double distance;
                algorithm.best_source(destinations_vertex_list[j], distance);
                pairwise_distance_poi_to_vertex[i][j] = distance;
            }
        }
        memory_usage += face_sequence_index_list_size * sizeof(int);
    }

    auto stop_construction_time = std::chrono::high_resolution_clock::now();
    auto duration_construction_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_construction_time - start_construction_time);
    construction_time = duration_construction_time.count();
}

void pre_or_post_WSPD_Oracle_and_pre_WSPD_Oracle_Adapt_query(
    geodesic::Mesh *mesh, int geo_tree_node_id, std::vector<GeoNode *> &all_poi,
    std::unordered_map<int, GeoNode *> &geo_node_in_partition_tree_unordered_map, std::unordered_map<int, GeoPair *> &geopairs, int source_poi_index,
    int destination_poi_index, double &query_time, double &approximate_distance, std::vector<geodesic::SurfacePoint> &approximate_path)
{
    auto start_query_time = std::chrono::high_resolution_clock::now();

    bool one_path;
    std::vector<int> face_sequence_index_list;

    three_paths_query_geo(mesh, geo_tree_node_id, 1, source_poi_index, destination_poi_index, one_path, all_poi, geo_node_in_partition_tree_unordered_map, geopairs, approximate_distance, approximate_path, face_sequence_index_list);
    approximate_distance = round(approximate_distance * 1000000000.0) / 1000000000.0;

    auto stop_query_time = std::chrono::high_resolution_clock::now();
    auto duration_query_time = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_query_time - start_query_time);
    query_time = duration_query_time.count();
    query_time /= 1000000;
}

void post_WSPD_Oracle_Adapt_update(int poi_num, geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
                                   geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list,
                                   int geo_tree_node_id, std::vector<GeoNode *> &all_poi,
                                   std::unordered_map<int, GeoNode *> &geo_node_in_partition_tree_unordered_map,
                                   std::unordered_map<int, GeoPair *> &geopairs,
                                   std::vector<std::vector<double>> &pairwise_distance_poi_to_poi,
                                   std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> &pairwise_path_poi_to_poi,
                                   std::vector<std::vector<double>> &pairwise_distance_poi_to_vertex,
                                   double &update_time, double &memory_usage)
{
    auto start_update_time = std::chrono::high_resolution_clock::now();

    geodesic::GeodesicAlgorithmExact post_algorithm(post_mesh);

    std::vector<geodesic::SurfacePoint> one_source_poi_list;
    std::vector<geodesic::SurfacePoint> destinations_poi_list;
    double const distance_limit = geodesic::GEODESIC_INF;

    int changed_pairwise_path_poi_to_poi_size = 0;
    int changed_pairwise_distance_poi_to_poi_size = 0;

    // pairwise_distance_poi_to_poi: if there is an edge between s and t, store the distance, otherwise, set it to be -1
    // pairwise_distance_poi_to_poi_with_appro: if there is an edge between s and t, store the distance, otherwise, store the approximate distance using three path query
    // this is used for checking whether a poi is close to changed area or not (calculating the radius of the circle)
    // pairwise_distance_poi_to_poi_changed: all false
    // pairwise_path_poi_to_poi: if there is an edge between s and t, store the path, otherwise, set it to be empty
    // pre_face_sequence_index_list: all face sequence using pre WSPD_Adapt query
    std::vector<std::vector<double>> pairwise_distance_poi_to_poi_with_appro;
    std::vector<std::vector<bool>> pairwise_distance_poi_to_poi_changed;
    std::vector<std::vector<std::vector<int>>> pre_face_sequence_index_list;

    pairwise_distance_poi_to_poi.clear();
    pairwise_distance_poi_to_poi_with_appro.clear();
    pairwise_distance_poi_to_poi_changed.clear();
    pairwise_path_poi_to_poi.clear();
    pre_face_sequence_index_list.clear();

    for (int i = 0; i < poi_num; i++)
    {
        std::vector<double> one_distance;
        std::vector<double> one_distance_with_appro;
        std::vector<bool> one_distance_changed;
        std::vector<std::vector<geodesic::SurfacePoint>> one_path;
        std::vector<std::vector<int>> one_poi_to_other_poi_pre_face_sequence_index_list;
        std::vector<int> one_poi_pre_face_sequence_index_list;

        one_distance.clear();
        one_distance_with_appro.clear();
        one_distance_changed.clear();
        one_path.clear();
        one_poi_to_other_poi_pre_face_sequence_index_list.clear();
        one_poi_pre_face_sequence_index_list.clear();

        for (int j = 0; j < i; j++)
        {
            one_poi_to_other_poi_pre_face_sequence_index_list.push_back(one_poi_pre_face_sequence_index_list);
        }

        for (int j = i; j < poi_num; j++)
        {
            double distance;
            std::vector<geodesic::SurfacePoint> path;
            path.clear();
            one_poi_pre_face_sequence_index_list.clear();

            if (i != j)
            {
                bool one_path_bool;
                three_paths_query_geo(pre_mesh, geo_tree_node_id, 1, i, j, one_path_bool, all_poi, geo_node_in_partition_tree_unordered_map, geopairs, distance, path, one_poi_pre_face_sequence_index_list);
                if (one_path_bool)
                {
                    one_distance.push_back(distance);
                    one_distance_with_appro.push_back(distance);
                    if ((path[0].getx() == geodesic::SurfacePoint(&pre_mesh->vertices()[pre_poi_list[j]]).getx()) &&
                        (path[0].gety() == geodesic::SurfacePoint(&pre_mesh->vertices()[pre_poi_list[j]]).gety()) &&
                        (path[0].getz() == geodesic::SurfacePoint(&pre_mesh->vertices()[pre_poi_list[j]]).getz()) &&
                        (path[path.size() - 1].getx() == geodesic::SurfacePoint(&pre_mesh->vertices()[pre_poi_list[i]]).getx()) &&
                        (path[path.size() - 1].gety() == geodesic::SurfacePoint(&pre_mesh->vertices()[pre_poi_list[i]]).gety()) &&
                        (path[path.size() - 1].getz() == geodesic::SurfacePoint(&pre_mesh->vertices()[pre_poi_list[i]]).getz()))
                    {
                        // do nothing
                    }
                    else if ((path[0].getx() == geodesic::SurfacePoint(&pre_mesh->vertices()[pre_poi_list[i]]).getx()) &&
                             (path[0].gety() == geodesic::SurfacePoint(&pre_mesh->vertices()[pre_poi_list[i]]).gety()) &&
                             (path[0].getz() == geodesic::SurfacePoint(&pre_mesh->vertices()[pre_poi_list[i]]).getz()) &&
                             (path[path.size() - 1].getx() == geodesic::SurfacePoint(&pre_mesh->vertices()[pre_poi_list[j]]).getx()) &&
                             (path[path.size() - 1].gety() == geodesic::SurfacePoint(&pre_mesh->vertices()[pre_poi_list[j]]).gety()) &&
                             (path[path.size() - 1].getz() == geodesic::SurfacePoint(&pre_mesh->vertices()[pre_poi_list[j]]).getz()))
                    {
                        std::reverse(path.begin(), path.end());
                    }
                    else
                    {
                        assert(false);
                    }
                }
                else
                {
                    one_distance.push_back(-1);
                    one_distance_with_appro.push_back(distance);
                    path.clear();
                }
            }
            else
            {
                one_distance.push_back(0);
                one_distance_with_appro.push_back(0);
                path.push_back(geodesic::SurfacePoint(&pre_mesh->vertices()[pre_poi_list[i]]));
            }

            one_distance_changed.push_back(false);
            one_path.push_back(path);
            one_poi_to_other_poi_pre_face_sequence_index_list.push_back(one_poi_pre_face_sequence_index_list);
        }
        pairwise_distance_poi_to_poi.push_back(one_distance);
        pairwise_distance_poi_to_poi_with_appro.push_back(one_distance_with_appro);
        pairwise_distance_poi_to_poi_changed.push_back(one_distance_changed);
        pairwise_path_poi_to_poi.push_back(one_path);
        pre_face_sequence_index_list.push_back(one_poi_to_other_poi_pre_face_sequence_index_list);
    }

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
            changed_face_index_list.push_back(i);
            changed_face_index_unordered_map[i] = i;
        }
    }
    assert(changed_face_index_list.size() == changed_face_index_unordered_map.size());

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

    // stores the changed poi index such that (1) the index of poi changed, (2) this poi is in the changed face, (3) the path that connects this poi passes the changed face
    std::vector<int> changed_poi_index_list(poi_num, 0);

    // stores only the poi in the changed face (used in calculation the 2D distance between these pois and other pois)
    std::vector<int> poi_in_the_changed_face_index_list;
    poi_in_the_changed_face_index_list.clear();

    // for the pre and post terrain, the index of poi may change, but the actual vertex on the terrain that these two pois stand for are very close
    // if the poi for the pre and post list changed, then the value in this list becomes 1, which indicates changes
    for (int i = 0; i < poi_num; i++)
    {
        if (pre_poi_list[i] != post_poi_list[i])
        {
            changed_poi_index_list[i] = 1;
        }
    }

    // for the pre and post terrain, even if the index of poi not changes, if the poi is on the changed face, we also indicate it by 3 in the following list
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
                changed_poi_index_list[i] = 3;
                poi_in_the_changed_face_index_list.push_back(i);
                break;
            }
        }
    }

    // if two pois are not in the changed face, but the path passes the changed face, we also indicate the poi by 4
    // note that for a poi, if one of the path that connects this poi passes the changed face, then we need to run SSAD
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
            }
            if (break_loop)
            {
                break;
            }
        }
    }

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
        }

        post_algorithm.propagate(one_source_poi_list, distance_limit, &destinations_poi_list);

        int index = 0;
        for (int j = 0; j < poi_num; j++)
        {
            if ((changed_poi_index_list[j] == 1 || changed_poi_index_list[j] == 3 || changed_poi_index_list[j] == 4) && j < i)
            {
                index++;
                continue;
            }
            std::vector<geodesic::SurfacePoint> path;
            post_algorithm.trace_back(destinations_poi_list[j - index], path);
            changed_pairwise_path_poi_to_poi_size += path.size();
            changed_pairwise_distance_poi_to_poi_size++;

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
    }

    // if two pois are not in the changed face, and the path doesn't pass the changed face, but one of pois is close to the
    // changed face, so we may need to update the new path with these two pois as endpoints on the new terrain, the following
    // is to check whether the original path is too close to the changed face or not, if so, we directly update the path
    // we also indicate this type of poi in changed_poi_index_list, but indicate it as 2 for clarification

    std::vector<double> euclidean_distance_of_poi_to_changed_face(poi_num, 0);
    std::vector<std::pair<double, int>> euclidean_distance_of_poi_to_changed_face_and_original_index;

    for (int i = 0; i < poi_num; i++)
    {
        if (changed_poi_index_list[i] == 1 || changed_poi_index_list[i] == 3 || changed_poi_index_list[i] == 4)
        {
            continue;
        }
        if (poi_in_the_changed_face_index_list.size() > 0)
        {
            for (int j = 0; j < poi_in_the_changed_face_index_list.size(); j++)
            {
                euclidean_distance_of_poi_to_changed_face[i] += euclidean_distance(post_mesh->vertices()[post_poi_list[i]].x(), post_mesh->vertices()[post_poi_list[i]].y(),
                                                                                   post_mesh->vertices()[poi_in_the_changed_face_index_list[j]].x(), post_mesh->vertices()[poi_in_the_changed_face_index_list[j]].y());
            }
        }
        else
        {
            euclidean_distance_of_poi_to_changed_face[i] = euclidean_distance(post_mesh->vertices()[post_poi_list[i]].x(), post_mesh->vertices()[post_poi_list[i]].y(),
                                                                              post_mesh->vertices()[changed_vertex_index_list[0]].x(), post_mesh->vertices()[changed_vertex_index_list[0]].y());
        }
    }

    sort_min_to_max_and_get_original_index(euclidean_distance_of_poi_to_changed_face, euclidean_distance_of_poi_to_changed_face_and_original_index);

    assert(euclidean_distance_of_poi_to_changed_face.size() == euclidean_distance_of_poi_to_changed_face_and_original_index.size());

    for (int i = 0; i < euclidean_distance_of_poi_to_changed_face_and_original_index.size(); i++)
    {
        if (euclidean_distance_of_poi_to_changed_face_and_original_index[i].first == 0)
        {
            continue;
        }
        int current_poi_index = euclidean_distance_of_poi_to_changed_face_and_original_index[i].second;

        assert(pairwise_distance_poi_to_poi[current_poi_index].size() == pairwise_distance_poi_to_poi_changed[current_poi_index].size());
        assert(pairwise_distance_poi_to_poi_with_appro[current_poi_index].size() == pairwise_distance_poi_to_poi_changed[current_poi_index].size());

        double max_distance = 0;
        for (int j = 0; j < euclidean_distance_of_poi_to_changed_face_and_original_index.size(); j++)
        {
            if (euclidean_distance_of_poi_to_changed_face_and_original_index[j].first == 0)
            {
                continue;
            }
            int checking_poi_index = euclidean_distance_of_poi_to_changed_face_and_original_index[j].second;

            if (current_poi_index <= checking_poi_index)
            {
                if (pairwise_distance_poi_to_poi_changed[current_poi_index][checking_poi_index - current_poi_index])
                {
                    continue;
                }
                max_distance = std::max(max_distance, pairwise_distance_poi_to_poi_with_appro[current_poi_index][checking_poi_index - current_poi_index]);
            }
            else
            {
                if (pairwise_distance_poi_to_poi_changed[checking_poi_index][current_poi_index - checking_poi_index])
                {
                    continue;
                }
                max_distance = std::max(max_distance, pairwise_distance_poi_to_poi_with_appro[checking_poi_index][current_poi_index - checking_poi_index]);
            }
        }

        // if the current poi is too close to the changed face, we need to run SSAD for this poi to update it path on the new terrain
        for (int k = 0; k < changed_vertex_index_list.size(); k++)
        {
            if (pairwise_distance_poi_to_vertex[current_poi_index][changed_vertex_index_list[k]] < max_distance)
            {
                changed_poi_index_list[current_poi_index] = 2;

                std::vector<int> destinations_poi_index_list;
                destinations_poi_index_list.clear();
                one_source_poi_list.clear();
                destinations_poi_list.clear();
                one_source_poi_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_poi_list[current_poi_index]]));

                for (int j = 0; j < euclidean_distance_of_poi_to_changed_face_and_original_index.size(); j++)
                {
                    int the_other_poi_index = euclidean_distance_of_poi_to_changed_face_and_original_index[j].second;

                    if ((current_poi_index <= the_other_poi_index && !pairwise_distance_poi_to_poi_changed[current_poi_index][the_other_poi_index - current_poi_index]) ||
                        (current_poi_index > the_other_poi_index && !pairwise_distance_poi_to_poi_changed[the_other_poi_index][current_poi_index - the_other_poi_index]))
                    {
                        destinations_poi_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_poi_list[the_other_poi_index]]));
                        destinations_poi_index_list.push_back(the_other_poi_index);
                    }
                }

                post_algorithm.propagate(one_source_poi_list, distance_limit, &destinations_poi_list);

                assert(destinations_poi_list.size() == destinations_poi_index_list.size());
                for (int j = 0; j < destinations_poi_index_list.size(); j++)
                {
                    int the_other_poi_index = destinations_poi_index_list[j];

                    std::vector<geodesic::SurfacePoint> path;
                    post_algorithm.trace_back(destinations_poi_list[j], path);
                    changed_pairwise_path_poi_to_poi_size += path.size();
                    changed_pairwise_distance_poi_to_poi_size++;

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
    memory_usage += post_algorithm.get_memory() + changed_face_index_list.size() * sizeof(int);

    auto stop_update_time = std::chrono::high_resolution_clock::now();
    auto duration_update_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_update_time - start_update_time);
    update_time = duration_update_time.count();
}

void pre_or_post_EAR_Oracle_and_pre_EAR_Oracle_Adapt_construction(
    int sqrt_num_of_box, int poi_num, geodesic::Mesh *mesh, std::vector<int> &poi_list, double epsilon,
    std::vector<int> &highway_node_list, std::unordered_map<int, std::unordered_map<int, int>> &highway_node_id_with_box_id_map,
    int &geo_tree_node_id, std::vector<GeoNode *> &all_highway_node, std::unordered_map<int, GeoNode *> &geo_node_in_partition_tree_unordered_map,
    std::unordered_map<int, GeoPair *> &geopairs, std::unordered_map<int, int> &highway_node_unordered_map, std::vector<std::vector<std::vector<int>>> &pre_face_sequence_index_list,
    std::vector<std::vector<double>> &pairwise_distance_poi_to_poi, std::vector<std::vector<bool>> &pairwise_distance_poi_to_poi_changed,
    std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> &pairwise_path_poi_to_poi, std::vector<std::vector<double>> &pairwise_distance_poi_to_vertex,
    std::unordered_map<int, double> &distance_poi_to_highway_node_map, std::unordered_map<int, std::vector<geodesic::SurfacePoint>> &path_poi_to_highway_node_map,
    bool run_EAR_Oracle_Adapt, double &construction_time, double &memory_usage, double &EAR_Oracle_size, int &EAR_Oracle_edge_num, double &EAR_Oracle_weight)
{
    auto start_construction_time = std::chrono::high_resolution_clock::now();

    geodesic::GeodesicAlgorithmExact algorithm(mesh);

    std::unordered_map<int, int> highway_node_id_map;
    highway_node_id_map.clear();
    highway_node_id_with_box_id_map.clear();
    divide_mesh_into_box(mesh, sqrt_num_of_box, highway_node_id_map, highway_node_id_with_box_id_map);

    highway_node_list.clear();
    for (auto i : highway_node_id_map)
    {
        highway_node_list.push_back(i.first);
        // std::cout << i.first << std::endl;
    }
    // std::cout << "size:" << highway_node_list.size() << std::endl;

    std::unordered_map<int, double> pre_distance_highway_node_to_highway_node_map;
    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> pre_path_highway_node_to_highway_node_map;
    double pre_memory_usage = 0;

    pre_compute_EAR_Oracle_highway_node(mesh, highway_node_list, pre_distance_highway_node_to_highway_node_map,
                                        pre_path_highway_node_to_highway_node_map, pre_memory_usage);

    geopairs.clear();
    all_highway_node.clear();
    std::vector<std::pair<int, GeoNode *>> highway_nodes;
    highway_nodes.clear();
    highway_node_unordered_map.clear();

    for (int i = 0; i < highway_node_list.size(); i++)
    {
        GeoNode *n = new GeoNode(highway_node_list[i], 0);
        all_highway_node.push_back(n);
        std::pair<int, GeoNode *> m(highway_node_list[i], n);
        highway_nodes.push_back(m);
        highway_node_unordered_map[highway_node_list[i]] = i;
    }

    double radius = 0;
    stx::btree<int, GeoNode *> highway_nodes_B_tree(highway_nodes.begin(), highway_nodes.end());

    for (int i = 0; i < highway_node_list.size(); i++)
    {
        int x_in_highway_node_list = 0;
        int y_in_highway_node_list = i;
        int x_y_in_highway_node_list;
        if (x_in_highway_node_list <= y_in_highway_node_list)
        {
            hash_function_two_keys_to_one_key(highway_node_list.size(), x_in_highway_node_list, y_in_highway_node_list, x_y_in_highway_node_list);
        }
        else
        {
            hash_function_two_keys_to_one_key(highway_node_list.size(), y_in_highway_node_list, x_in_highway_node_list, x_y_in_highway_node_list);
        }
        radius = std::max(pre_distance_highway_node_to_highway_node_map[x_y_in_highway_node_list], radius);
    }
    GeoNode *root_geo = new GeoNode(0, highway_node_list[0], radius);

    stx::btree<int, GeoNode *> highway_nodes_as_center_each_parent_layer;
    highway_nodes_as_center_each_parent_layer.clear();
    build_geo_tree(geo_tree_node_id, mesh, *root_geo, highway_node_list.size(), highway_nodes_B_tree, highway_nodes_as_center_each_parent_layer, pre_distance_highway_node_to_highway_node_map, pre_path_highway_node_to_highway_node_map, highway_node_unordered_map);

    std::vector<GeoNode *> partition_tree_to_compressed_partition_tree_to_be_removed_nodes;
    partition_tree_to_compressed_partition_tree_to_be_removed_nodes.clear();
    geo_node_in_partition_tree_unordered_map.clear();
    partition_tree_to_compressed_partition_tree(*root_geo, partition_tree_to_compressed_partition_tree_to_be_removed_nodes, geo_node_in_partition_tree_unordered_map);

    std::unordered_map<int, int> geo_pair_unordered_map;
    geo_pair_unordered_map.clear();
    int pairwise_path_highway_node_to_highway_node_size = 0;
    int face_sequence_index_list_size = 0;
    generate_geo_pair(geo_tree_node_id, EAR_Oracle_edge_num, EAR_Oracle_weight, mesh, *root_geo, *root_geo, epsilon, geopairs, highway_node_unordered_map, geo_pair_unordered_map, pre_distance_highway_node_to_highway_node_map, pre_path_highway_node_to_highway_node_map, pairwise_path_highway_node_to_highway_node_size, face_sequence_index_list_size);
    EAR_Oracle_size = EAR_Oracle_edge_num * sizeof(double) + pairwise_path_highway_node_to_highway_node_size * sizeof(geodesic::SurfacePoint);

    poi_to_highway_node_path(mesh, sqrt_num_of_box, highway_node_id_with_box_id_map, poi_list, distance_poi_to_highway_node_map, path_poi_to_highway_node_map, EAR_Oracle_size, memory_usage);

    memory_usage += pre_memory_usage + (geo_tree_node_id + 1) * sizeof(GeoNode) + EAR_Oracle_edge_num * sizeof(double) + pairwise_path_highway_node_to_highway_node_size * sizeof(geodesic::SurfacePoint);

    if (run_EAR_Oracle_Adapt)
    {
        std::vector<geodesic::SurfacePoint> one_source_poi_list2;
        std::vector<geodesic::SurfacePoint> destinations_poi_list2;
        std::vector<geodesic::SurfacePoint> destinations_vertex_list2;
        std::vector<std::vector<int>> one_poi_to_other_poi_pre_face_sequence_index_list2;
        std::vector<int> one_poi_pre_face_sequence_index_list2;
        for (int i = 0; i < poi_list.size(); i++)
        {
            std::vector<double> current_poi_to_other_poi_distance2;
            current_poi_to_other_poi_distance2.clear();
            std::vector<bool> current_poi_to_other_poi_distance_changed2;
            current_poi_to_other_poi_distance_changed2.clear();
            std::vector<std::vector<geodesic::SurfacePoint>> current_poi_to_other_poi_path2;
            current_poi_to_other_poi_path2.clear();
            one_source_poi_list2.clear();
            destinations_poi_list2.clear();
            destinations_vertex_list2.clear();
            one_poi_to_other_poi_pre_face_sequence_index_list2.clear();
            one_poi_pre_face_sequence_index_list2.clear();
            one_source_poi_list2.push_back(geodesic::SurfacePoint(&mesh->vertices()[poi_list[i]]));

            for (int j = i; j < poi_num; j++)
            {
                destinations_poi_list2.push_back(geodesic::SurfacePoint(&mesh->vertices()[poi_list[j]]));
            }
            for (int j = 0; j < mesh->vertices().size(); j++)
            {
                destinations_vertex_list2.push_back(geodesic::SurfacePoint(&mesh->vertices()[j]));
            }
            algorithm.propagate(one_source_poi_list2, geodesic::GEODESIC_INF, &destinations_vertex_list2);
            for (int j = 0; j < i; j++)
            {
                one_poi_to_other_poi_pre_face_sequence_index_list2.push_back(one_poi_pre_face_sequence_index_list2);
            }
            for (int j = i; j < poi_num; j++)
            {
                std::vector<geodesic::SurfacePoint> path;
                algorithm.trace_back(destinations_poi_list2[j - i], path);
                current_poi_to_other_poi_distance2.push_back(length(path));
                current_poi_to_other_poi_distance_changed2.push_back(false);
                current_poi_to_other_poi_path2.push_back(path);
                get_face_sequence(path, one_poi_pre_face_sequence_index_list2);
                one_poi_to_other_poi_pre_face_sequence_index_list2.push_back(one_poi_pre_face_sequence_index_list2);
            }
            pairwise_distance_poi_to_poi.push_back(current_poi_to_other_poi_distance2);
            pairwise_distance_poi_to_poi_changed.push_back(current_poi_to_other_poi_distance_changed2);
            pairwise_path_poi_to_poi.push_back(current_poi_to_other_poi_path2);
            pre_face_sequence_index_list.push_back(one_poi_to_other_poi_pre_face_sequence_index_list2);

            for (int j = 0; j < mesh->vertices().size(); j++)
            {
                double distance;
                algorithm.best_source(destinations_vertex_list2[j], distance);
                pairwise_distance_poi_to_vertex[i][j] = distance;
            }
        }
        memory_usage += face_sequence_index_list_size * sizeof(int);
    }

    auto stop_construction_time = std::chrono::high_resolution_clock::now();
    auto duration_construction_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_construction_time - start_construction_time);
    construction_time = duration_construction_time.count();
}

void pre_or_post_EAR_Oracle_and_pre_EAR_Oracle_Adapt_query(
    geodesic::Mesh *mesh, std::vector<int> &poi_list,
    int sqrt_num_of_box, int geo_tree_node_id, std::unordered_map<int, std::unordered_map<int, int>> &highway_node_id_with_box_id_map,
    std::vector<GeoNode *> &all_highway_node, std::unordered_map<int, GeoNode *> &geo_node_in_partition_tree_unordered_map,
    std::unordered_map<int, GeoPair *> &geopairs,
    std::unordered_map<int, double> &distance_poi_to_highway_node_map, std::unordered_map<int, std::vector<geodesic::SurfacePoint>> &path_poi_to_highway_node_map,
    int source_poi_index, int destination_poi_index, double &query_time, double &approximate_distance, std::vector<geodesic::SurfacePoint> &approximate_path)
{
    auto start_query_time = std::chrono::high_resolution_clock::now();

    std::vector<int> face_sequence_index_list;

    EAR_Oracle_query(mesh, poi_list, sqrt_num_of_box, geo_tree_node_id, highway_node_id_with_box_id_map,
                     all_highway_node, geo_node_in_partition_tree_unordered_map,
                     geopairs, distance_poi_to_highway_node_map, path_poi_to_highway_node_map,
                     source_poi_index, destination_poi_index, approximate_distance, approximate_path,
                     face_sequence_index_list);
    approximate_distance = round(approximate_distance * 1000000000.0) / 1000000000.0;

    auto stop_query_time = std::chrono::high_resolution_clock::now();
    auto duration_query_time = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_query_time - start_query_time);
    query_time = duration_query_time.count();
    query_time /= 1000000;
}

void post_EAR_Oracle_Adapt_update(int poi_num, int sqrt_num_of_box, geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
                                  geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list,
                                  int geo_tree_node_id, std::unordered_map<int, std::unordered_map<int, int>> &highway_node_id_with_box_id_map,
                                  std::vector<GeoNode *> &all_highway_node,
                                  std::unordered_map<int, GeoNode *> &geo_node_in_partition_tree_unordered_map,
                                  std::unordered_map<int, GeoPair *> &geopairs,
                                  std::unordered_map<int, double> &distance_poi_to_highway_node_map,
                                  std::unordered_map<int, std::vector<geodesic::SurfacePoint>> &path_poi_to_highway_node_map,
                                  std::vector<std::vector<std::vector<int>>> &pre_face_sequence_index_list,
                                  std::vector<std::vector<double>> &pairwise_distance_poi_to_poi,
                                  std::vector<std::vector<bool>> &pairwise_distance_poi_to_poi_changed,
                                  std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> &pairwise_path_poi_to_poi,
                                  std::vector<std::vector<double>> &pairwise_distance_poi_to_vertex,
                                  double &update_time, double &memory_usage)
{
    auto start_update_time = std::chrono::high_resolution_clock::now();

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
            changed_face_index_list.push_back(i);
            changed_face_index_unordered_map[i] = i;
        }
    }
    assert(changed_face_index_list.size() == changed_face_index_unordered_map.size());

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

    // stores the changed poi index such that (1) the index of poi changed, (2) this poi is in the changed face, (3) the path that connects this poi passes the changed face
    std::vector<int> changed_poi_index_list(poi_num, 0);

    // stores only the poi in the changed face (used in calculation the 2D distance between these pois and other pois)
    std::vector<int> poi_in_the_changed_face_index_list;
    poi_in_the_changed_face_index_list.clear();

    // for the pre and post terrain, the index of poi may change, but the actual vertex on the terrain that these two pois stand for are very close
    // if the poi for the pre and post list changed, then the value in this list becomes 1, which indicates changes
    for (int i = 0; i < poi_num; i++)
    {
        if (pre_poi_list[i] != post_poi_list[i])
        {
            changed_poi_index_list[i] = 1;
        }
    }

    // for the pre and post terrain, even if the index of poi not changes, if the poi is on the changed face, we also indicate it by 3 in the following list
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
                changed_poi_index_list[i] = 3;
                poi_in_the_changed_face_index_list.push_back(i);
                break;
            }
        }
    }

    // if two pois are not in the changed face, but the path passes the changed face, we also indicate the poi by 4
    // note that for a poi, if one of the path that connects this poi passes the changed face, then we need to run SSAD
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
            }
            if (break_loop)
            {
                break;
            }
        }
    }

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
        }

        post_algorithm.propagate(one_source_poi_list, distance_limit, &destinations_poi_list);

        int index = 0;
        for (int j = 0; j < poi_num; j++)
        {
            if ((changed_poi_index_list[j] == 1 || changed_poi_index_list[j] == 3 || changed_poi_index_list[j] == 4) && j < i)
            {
                index++;
                continue;
            }
            std::vector<geodesic::SurfacePoint> path;
            post_algorithm.trace_back(destinations_poi_list[j - index], path);
            changed_pairwise_path_poi_to_poi_size += path.size();
            changed_pairwise_distance_poi_to_poi_size++;

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
    }

    // if two pois are not in the changed face, and the path doesn't pass the changed face, but one of pois is close to the
    // changed face, so we may need to update the new path with these two pois as endpoints on the new terrain, the following
    // is to check whether the original path is too close to the changed face or not, if so, we directly update the path
    // we also indicate this type of poi in changed_poi_index_list, but indicate it as 2 for clarification

    std::vector<double> euclidean_distance_of_poi_to_changed_face(poi_num, 0);
    std::vector<std::pair<double, int>> euclidean_distance_of_poi_to_changed_face_and_original_index;

    for (int i = 0; i < poi_num; i++)
    {
        if (changed_poi_index_list[i] == 1 || changed_poi_index_list[i] == 3 || changed_poi_index_list[i] == 4)
        {
            continue;
        }
        if (poi_in_the_changed_face_index_list.size() > 0)
        {
            for (int j = 0; j < poi_in_the_changed_face_index_list.size(); j++)
            {
                euclidean_distance_of_poi_to_changed_face[i] += euclidean_distance(post_mesh->vertices()[post_poi_list[i]].x(), post_mesh->vertices()[post_poi_list[i]].y(),
                                                                                   post_mesh->vertices()[poi_in_the_changed_face_index_list[j]].x(), post_mesh->vertices()[poi_in_the_changed_face_index_list[j]].y());
            }
        }
        else
        {
            euclidean_distance_of_poi_to_changed_face[i] = euclidean_distance(post_mesh->vertices()[post_poi_list[i]].x(), post_mesh->vertices()[post_poi_list[i]].y(),
                                                                              post_mesh->vertices()[changed_vertex_index_list[0]].x(), post_mesh->vertices()[changed_vertex_index_list[0]].y());
        }
    }

    sort_min_to_max_and_get_original_index(euclidean_distance_of_poi_to_changed_face, euclidean_distance_of_poi_to_changed_face_and_original_index);

    assert(euclidean_distance_of_poi_to_changed_face.size() == euclidean_distance_of_poi_to_changed_face_and_original_index.size());

    for (int i = 0; i < euclidean_distance_of_poi_to_changed_face_and_original_index.size(); i++)
    {
        if (euclidean_distance_of_poi_to_changed_face_and_original_index[i].first == 0)
        {
            continue;
        }
        int current_poi_index = euclidean_distance_of_poi_to_changed_face_and_original_index[i].second;

        assert(pairwise_distance_poi_to_poi[current_poi_index].size() == pairwise_distance_poi_to_poi_changed[current_poi_index].size());

        double max_distance = 0;
        for (int j = 0; j < euclidean_distance_of_poi_to_changed_face_and_original_index.size(); j++)
        {
            if (euclidean_distance_of_poi_to_changed_face_and_original_index[j].first == 0)
            {
                continue;
            }
            int checking_poi_index = euclidean_distance_of_poi_to_changed_face_and_original_index[j].second;

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

        // if the current poi is too close to the changed face, we need to run SSAD for this poi to update it path on the new terrain
        for (int k = 0; k < changed_vertex_index_list.size(); k++)
        {
            if (pairwise_distance_poi_to_vertex[current_poi_index][changed_vertex_index_list[k]] < 2 * max_distance)
            {
                changed_poi_index_list[current_poi_index] = 2;

                std::vector<int> destinations_poi_index_list;
                destinations_poi_index_list.clear();
                one_source_poi_list.clear();
                destinations_poi_list.clear();
                one_source_poi_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_poi_list[current_poi_index]]));

                for (int j = 0; j < euclidean_distance_of_poi_to_changed_face_and_original_index.size(); j++)
                {
                    int the_other_poi_index = euclidean_distance_of_poi_to_changed_face_and_original_index[j].second;

                    if ((current_poi_index <= the_other_poi_index && !pairwise_distance_poi_to_poi_changed[current_poi_index][the_other_poi_index - current_poi_index]) ||
                        (current_poi_index > the_other_poi_index && !pairwise_distance_poi_to_poi_changed[the_other_poi_index][current_poi_index - the_other_poi_index]))
                    {
                        destinations_poi_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_poi_list[the_other_poi_index]]));
                        destinations_poi_index_list.push_back(the_other_poi_index);
                    }
                }

                post_algorithm.propagate(one_source_poi_list, distance_limit, &destinations_poi_list);

                assert(destinations_poi_list.size() == destinations_poi_index_list.size());
                for (int j = 0; j < destinations_poi_index_list.size(); j++)
                {
                    int the_other_poi_index = destinations_poi_index_list[j];

                    std::vector<geodesic::SurfacePoint> path;
                    post_algorithm.trace_back(destinations_poi_list[j], path);
                    changed_pairwise_path_poi_to_poi_size += path.size();
                    changed_pairwise_distance_poi_to_poi_size++;

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
    memory_usage += post_algorithm.get_memory() + changed_face_index_list.size() * sizeof(int);

    auto stop_update_time = std::chrono::high_resolution_clock::now();
    auto duration_update_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_update_time - start_update_time);
    update_time = duration_update_time.count();
}

void pre_or_post_CH_Fly_Algo_query(geodesic::Mesh *mesh, std::vector<int> &poi_list,
                                   int source_poi_index, int destination_poi_index,
                                   double &query_time, double &memory_usage, double &approximate_distance,
                                   std::vector<geodesic::SurfacePoint> &approximate_path)
{
    auto start_query_time = std::chrono::high_resolution_clock::now();

    geodesic::GeodesicAlgorithmExact algorithm(mesh);
    double const distance_limit = geodesic::GEODESIC_INF;
    geodesic::SurfacePoint source(&mesh->vertices()[poi_list[source_poi_index]]);
    geodesic::SurfacePoint destination(&mesh->vertices()[poi_list[destination_poi_index]]);
    std::vector<geodesic::SurfacePoint> one_source_poi_list(1, source);
    std::vector<geodesic::SurfacePoint> one_destination_poi_list(1, destination);
    algorithm.propagate(one_source_poi_list, distance_limit, &one_destination_poi_list);
    algorithm.trace_back(destination, approximate_path);
    approximate_distance = length(approximate_path);
    approximate_distance = round(approximate_distance * 1000000000.0) / 1000000000.0;
    memory_usage += algorithm.get_memory() + approximate_path.size() * sizeof(geodesic::SurfacePoint) + sizeof(double);

    auto stop_query_time = std::chrono::high_resolution_clock::now();
    auto duration_query_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_query_time - start_query_time);
    query_time = duration_query_time.count();
    query_time /= 1000;
}

void pre_or_post_K_Fly_Algo_query(geodesic::Mesh *mesh, std::vector<int> &poi_list,
                                  double epsilon, int source_poi_index, int destination_poi_index,
                                  double &query_time, double &memory_usage, double &approximate_distance,
                                  std::vector<geodesic::SurfacePoint> &approximate_path)
{
    auto start_query_time = std::chrono::high_resolution_clock::now();

    double subdivision_level = epslion_to_subdivision_level(epsilon);
    geodesic::GeodesicAlgorithmSubdivision algorithm(mesh, subdivision_level);
    double const distance_limit = geodesic::GEODESIC_INF;
    geodesic::SurfacePoint source(&mesh->vertices()[poi_list[source_poi_index]]);
    geodesic::SurfacePoint destination(&mesh->vertices()[poi_list[destination_poi_index]]);
    std::vector<geodesic::SurfacePoint> one_source_poi_list(1, source);
    std::vector<geodesic::SurfacePoint> one_destination_poi_list(1, destination);
    algorithm.propagate(one_source_poi_list, distance_limit, &one_destination_poi_list);
    algorithm.trace_back(destination, approximate_path);
    approximate_distance = length(approximate_path);
    approximate_distance = round(approximate_distance * 1000000000.0) / 1000000000.0;
    memory_usage += algorithm.get_memory() + approximate_path.size() * sizeof(geodesic::SurfacePoint) + sizeof(double);

    auto stop_query_time = std::chrono::high_resolution_clock::now();
    auto duration_query_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_query_time - start_query_time);
    query_time = duration_query_time.count();
    query_time /= 1000;
}

void FU_Oracle_RanUpdSeq_NoDistAppr(int poi_num, geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
                                    geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list, double epsilon,
                                    int source_poi_index, int destination_poi_index,
                                    double post_exact_distance, int pre_MST_weight, int post_MST_weight,
                                    std::string write_file_header, bool RanUpdSeq_not_NoDistAppr)
{
    std::vector<std::vector<std::vector<int>>> pre_face_sequence_index_list;
    std::vector<std::vector<double>> pairwise_distance_poi_to_poi;
    std::vector<std::vector<bool>> pairwise_distance_poi_to_poi_changed;
    std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> pairwise_path_poi_to_poi;
    std::vector<std::vector<double>> pairwise_distance_poi_to_vertex(poi_num, std::vector<double>(pre_mesh->vertices().size(), 0));

    double pre_construction_time = 0;
    double pre_hash_mapping_time = 0;
    double pre_memory_usage = 0;
    double pre_complete_graph_size = 0;
    int pre_index_edge_num = 0;
    double pre_index_weight = 0;

    pre_complete_graph_construction(poi_num, pre_mesh, pre_poi_list, pre_face_sequence_index_list,
                                    pairwise_distance_poi_to_poi, pairwise_distance_poi_to_poi_changed,
                                    pairwise_path_poi_to_poi, pairwise_distance_poi_to_vertex,
                                    pre_construction_time, pre_memory_usage);
    std::unordered_map<int, double> pairwise_distance_poi_to_poi_pre_map;
    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> pairwise_path_poi_to_poi_pre_map;
    get_pairwise_distance_and_path_poi_to_poi_map(pairwise_distance_poi_to_poi, pairwise_distance_poi_to_poi_pre_map,
                                                  pairwise_path_poi_to_poi, pairwise_path_poi_to_poi_pre_map,
                                                  pre_complete_graph_size, pre_index_edge_num, pre_index_weight, pre_hash_mapping_time);

    std::cout << "Pre terrain construction time: " << (pre_construction_time + pre_hash_mapping_time) << " ms" << std::endl;
    std::cout << "Pre terrain memory usage: " << (pre_memory_usage + pre_complete_graph_size) / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index size: " << pre_complete_graph_size / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index edge number: " << pre_index_edge_num << " edges" << std::endl;
    std::cout << "Pre terrain index/MST weight: " << pre_index_weight / pre_MST_weight << std::endl;
    std::cout << std::endl;

    std::ofstream ofs1("../output/output.txt", std::ios_base::app);
    if (RanUpdSeq_not_NoDistAppr)
    {
        ofs1 << "== FU_Oracle_RanUpdSeq ==\n";
    }
    else
    {
        ofs1 << "== FU_Oracle_NoDistAppr ==\n";
    }
    ofs1 << write_file_header << "\t"
         << (pre_construction_time + pre_hash_mapping_time) << "\t"
         << (pre_memory_usage + pre_complete_graph_size) / 1e6 << "\t"
         << pre_complete_graph_size / 1e6 << "\t"
         << pre_index_edge_num << "\t"
         << pre_index_weight / pre_MST_weight << "\t";
    ofs1.close();

    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> post_path_poi_to_poi_map;
    post_path_poi_to_poi_map.clear();
    double post_update_time = 0;
    double post_HGS_time = 0;
    double post_query_time = 0;
    double post_memory_usage = 0;
    double post_complete_graph_size = 0;
    double post_hierarchy_greedy_spanner_size = 0;
    int post_index_edge_num = 0;
    double post_index_weight = 0;
    double post_approximate_distance = 0;
    std::vector<geodesic::SurfacePoint> post_approximate_path;
    post_approximate_path.clear();
    Graph graph(poi_num);

    post_complete_graph_update_RanUpdSeq_NoDistAppr(poi_num, pre_mesh, pre_poi_list, post_mesh, post_poi_list,
                                                    pre_face_sequence_index_list, pairwise_distance_poi_to_poi,
                                                    pairwise_distance_poi_to_poi_changed,
                                                    pairwise_path_poi_to_poi, pairwise_distance_poi_to_vertex,
                                                    post_update_time, post_memory_usage);
    hierarchy_greedy_spanner(epsilon, graph, pairwise_distance_poi_to_poi, pairwise_path_poi_to_poi, post_path_poi_to_poi_map,
                             post_hierarchy_greedy_spanner_size, post_index_edge_num, post_index_weight, post_HGS_time, post_complete_graph_size);
    spanner_query(poi_num, graph, post_path_poi_to_poi_map, source_poi_index, destination_poi_index, post_approximate_distance,
                  post_approximate_path, post_query_time);

    std::cout << "Post terrain update time (CG): " << post_update_time << " ms" << std::endl;
    std::cout << "Post terrain update time (HieGreSpan): " << post_HGS_time << " ms" << std::endl;
    std::cout << "Post terrain query time: " << post_query_time << " ms" << std::endl;
    std::cout << "Post terrain memory usage: " << (post_memory_usage + post_complete_graph_size) / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index size: " << post_complete_graph_size / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain output size: " << post_hierarchy_greedy_spanner_size / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index edge number: " << post_index_edge_num << " edges" << std::endl;
    std::cout << "Post terrain index/MST weight: " << post_index_weight / post_MST_weight << std::endl;
    std::cout << "Post terrain approximate distance: " << post_approximate_distance << ", post terrain exact distance: " << post_exact_distance << ", distance error: " << post_approximate_distance / post_exact_distance - 1 << std::endl;

    std::ofstream ofs2("../output/output.txt", std::ios_base::app);
    ofs2 << post_update_time << "\t"
         << post_HGS_time << "\t"
         << post_query_time << "\t"
         << (post_memory_usage + post_complete_graph_size) / 1e6 << "\t"
         << post_complete_graph_size / 1e6 << "\t"
         << post_hierarchy_greedy_spanner_size / 1e6 << "\t"
         << post_index_edge_num << "\t"
         << post_index_weight / post_MST_weight << "\t"
         << post_approximate_distance / post_exact_distance - 1 << "\n\n";
    ofs2.close();
}

void FU_Oracle_FullRad(int poi_num, geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
                       geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list, double epsilon,
                       int source_poi_index, int destination_poi_index,
                       double post_exact_distance, int pre_MST_weight, int post_MST_weight,
                       std::string write_file_header)
{
    std::vector<std::vector<std::vector<int>>> pre_face_sequence_index_list;
    std::vector<std::vector<double>> pairwise_distance_poi_to_poi;
    std::vector<std::vector<bool>> pairwise_distance_poi_to_poi_changed;
    std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> pairwise_path_poi_to_poi;
    std::vector<std::vector<double>> pairwise_distance_poi_to_vertex(poi_num, std::vector<double>(pre_mesh->vertices().size(), 0));

    double pre_construction_time = 0;
    double pre_hash_mapping_time = 0;
    double pre_memory_usage = 0;
    double pre_complete_graph_size = 0;
    int pre_index_edge_num = 0;
    double pre_index_weight = 0;

    pre_complete_graph_construction(poi_num, pre_mesh, pre_poi_list, pre_face_sequence_index_list,
                                    pairwise_distance_poi_to_poi, pairwise_distance_poi_to_poi_changed,
                                    pairwise_path_poi_to_poi, pairwise_distance_poi_to_vertex,
                                    pre_construction_time, pre_memory_usage);
    std::unordered_map<int, double> pairwise_distance_poi_to_poi_pre_map;
    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> pairwise_path_poi_to_poi_pre_map;
    get_pairwise_distance_and_path_poi_to_poi_map(pairwise_distance_poi_to_poi, pairwise_distance_poi_to_poi_pre_map,
                                                  pairwise_path_poi_to_poi, pairwise_path_poi_to_poi_pre_map,
                                                  pre_complete_graph_size, pre_index_edge_num, pre_index_weight, pre_hash_mapping_time);

    std::cout << "Pre terrain construction time: " << (pre_construction_time + pre_hash_mapping_time) << " ms" << std::endl;
    std::cout << "Pre terrain memory usage: " << (pre_memory_usage + pre_complete_graph_size) / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index size: " << pre_complete_graph_size / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index edge number: " << pre_index_edge_num << " edges" << std::endl;
    std::cout << "Pre terrain index/MST weight: " << pre_index_weight / pre_MST_weight << std::endl;
    std::cout << std::endl;

    std::ofstream ofs1("../output/output.txt", std::ios_base::app);
    ofs1 << "== FU_Oracle_FullRad ==\n";
    ofs1 << write_file_header << "\t"
         << (pre_construction_time + pre_hash_mapping_time) << "\t"
         << (pre_memory_usage + pre_complete_graph_size) / 1e6 << "\t"
         << pre_complete_graph_size / 1e6 << "\t"
         << pre_index_edge_num << "\t"
         << pre_index_weight / pre_MST_weight << "\t";
    ofs1.close();

    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> post_path_poi_to_poi_map;
    post_path_poi_to_poi_map.clear();
    double post_update_time = 0;
    double post_HGS_time = 0;
    double post_query_time = 0;
    double post_memory_usage = 0;
    double post_complete_graph_size = 0;
    double post_hierarchy_greedy_spanner_size = 0;
    int post_index_edge_num = 0;
    double post_index_weight = 0;
    double post_approximate_distance = 0;
    std::vector<geodesic::SurfacePoint> post_approximate_path;
    post_approximate_path.clear();
    Graph graph(poi_num);

    post_complete_graph_update_FullRad(poi_num, pre_mesh, pre_poi_list, post_mesh, post_poi_list,
                                       pre_face_sequence_index_list, pairwise_distance_poi_to_poi,
                                       pairwise_distance_poi_to_poi_changed,
                                       pairwise_path_poi_to_poi, pairwise_distance_poi_to_vertex,
                                       post_update_time, post_memory_usage);
    hierarchy_greedy_spanner(epsilon, graph, pairwise_distance_poi_to_poi, pairwise_path_poi_to_poi, post_path_poi_to_poi_map,
                             post_hierarchy_greedy_spanner_size, post_index_edge_num, post_index_weight, post_HGS_time, post_complete_graph_size);
    spanner_query(poi_num, graph, post_path_poi_to_poi_map, source_poi_index, destination_poi_index, post_approximate_distance,
                  post_approximate_path, post_query_time);

    std::cout << "Post terrain update time (CG): " << post_update_time << " ms" << std::endl;
    std::cout << "Post terrain update time (HieGreSpan): " << post_HGS_time << " ms" << std::endl;
    std::cout << "Post terrain query time: " << post_query_time << " ms" << std::endl;
    std::cout << "Post terrain memory usage: " << (post_memory_usage + post_complete_graph_size) / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index size: " << post_complete_graph_size / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain output size: " << post_hierarchy_greedy_spanner_size / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index edge number: " << post_index_edge_num << " edges" << std::endl;
    std::cout << "Post terrain index/MST weight: " << post_index_weight / post_MST_weight << std::endl;
    std::cout << "Post terrain approximate distance: " << post_approximate_distance << ", post terrain exact distance: " << post_exact_distance << ", distance error: " << post_approximate_distance / post_exact_distance - 1 << std::endl;

    std::ofstream ofs2("../output/output.txt", std::ios_base::app);
    ofs2 << post_update_time << "\t"
         << post_HGS_time << "\t"
         << post_query_time << "\t"
         << (post_memory_usage + post_complete_graph_size) / 1e6 << "\t"
         << post_complete_graph_size / 1e6 << "\t"
         << post_hierarchy_greedy_spanner_size / 1e6 << "\t"
         << post_index_edge_num << "\t"
         << post_index_weight / post_MST_weight << "\t"
         << post_approximate_distance / post_exact_distance - 1 << "\n\n";
    ofs2.close();
}

void FU_Oracle_NoEffIntChe(int poi_num, geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
                           geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list, double epsilon,
                           int source_poi_index, int destination_poi_index,
                           double post_exact_distance, int pre_MST_weight, int post_MST_weight,
                           std::string write_file_header)
{
    std::vector<std::vector<std::vector<int>>> pre_face_sequence_index_list;
    std::vector<std::vector<double>> pairwise_distance_poi_to_poi;
    std::vector<std::vector<bool>> pairwise_distance_poi_to_poi_changed;
    std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> pairwise_path_poi_to_poi;
    std::vector<std::vector<double>> pairwise_distance_poi_to_vertex(poi_num, std::vector<double>(pre_mesh->vertices().size(), 0));

    double pre_construction_time = 0;
    double pre_hash_mapping_time = 0;
    double pre_memory_usage = 0;
    double pre_complete_graph_size = 0;
    int pre_index_edge_num = 0;
    double pre_index_weight = 0;

    pre_complete_graph_construction(poi_num, pre_mesh, pre_poi_list, pre_face_sequence_index_list,
                                    pairwise_distance_poi_to_poi, pairwise_distance_poi_to_poi_changed,
                                    pairwise_path_poi_to_poi, pairwise_distance_poi_to_vertex,
                                    pre_construction_time, pre_memory_usage);
    std::unordered_map<int, double> pairwise_distance_poi_to_poi_pre_map;
    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> pairwise_path_poi_to_poi_pre_map;
    get_pairwise_distance_and_path_poi_to_poi_map(pairwise_distance_poi_to_poi, pairwise_distance_poi_to_poi_pre_map,
                                                  pairwise_path_poi_to_poi, pairwise_path_poi_to_poi_pre_map,
                                                  pre_complete_graph_size, pre_index_edge_num, pre_index_weight, pre_hash_mapping_time);

    std::cout << "Pre terrain construction time: " << (pre_construction_time + pre_hash_mapping_time) << " ms" << std::endl;
    std::cout << "Pre terrain memory usage: " << (pre_memory_usage + pre_complete_graph_size) / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index size: " << pre_complete_graph_size / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index edge number: " << pre_index_edge_num << " edges" << std::endl;
    std::cout << "Pre terrain index/MST weight: " << pre_index_weight / pre_MST_weight << std::endl;
    std::cout << std::endl;

    std::ofstream ofs1("../output/output.txt", std::ios_base::app);
    ofs1 << "== FU_Oracle_NoEffIntChe ==\n";
    ofs1 << write_file_header << "\t"
         << (pre_construction_time + pre_hash_mapping_time) << "\t"
         << (pre_memory_usage + pre_complete_graph_size) / 1e6 << "\t"
         << pre_complete_graph_size / 1e6 << "\t"
         << pre_index_edge_num << "\t"
         << pre_index_weight / pre_MST_weight << "\t";
    ofs1.close();

    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> post_path_poi_to_poi_map;
    post_path_poi_to_poi_map.clear();
    double post_update_time = 0;
    double post_HGS_time = 0;
    double post_query_time = 0;
    double post_memory_usage = 0;
    double post_complete_graph_size = 0;
    double post_hierarchy_greedy_spanner_size = 0;
    int post_index_edge_num = 0;
    double post_index_weight = 0;
    double post_approximate_distance = 0;
    std::vector<geodesic::SurfacePoint> post_approximate_path;
    post_approximate_path.clear();
    Graph graph(poi_num);

    post_complete_graph_update_NoEffIntChe(poi_num, pre_mesh, pre_poi_list, post_mesh, post_poi_list,
                                           pre_face_sequence_index_list, pairwise_distance_poi_to_poi,
                                           pairwise_distance_poi_to_poi_changed,
                                           pairwise_path_poi_to_poi, pairwise_distance_poi_to_vertex,
                                           post_update_time, post_memory_usage);
    hierarchy_greedy_spanner(epsilon, graph, pairwise_distance_poi_to_poi, pairwise_path_poi_to_poi, post_path_poi_to_poi_map,
                             post_hierarchy_greedy_spanner_size, post_index_edge_num, post_index_weight, post_HGS_time, post_complete_graph_size);
    spanner_query(poi_num, graph, post_path_poi_to_poi_map, source_poi_index, destination_poi_index, post_approximate_distance,
                  post_approximate_path, post_query_time);

    std::cout << "Post terrain update time (CG): " << post_update_time << " ms" << std::endl;
    std::cout << "Post terrain update time (HieGreSpan): " << post_HGS_time << " ms" << std::endl;
    std::cout << "Post terrain query time: " << post_query_time << " ms" << std::endl;
    std::cout << "Post terrain memory usage: " << (post_memory_usage + post_complete_graph_size) / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index size: " << post_complete_graph_size / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain output size: " << post_hierarchy_greedy_spanner_size / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index edge number: " << post_index_edge_num << " edges" << std::endl;
    std::cout << "Post terrain index/MST weight: " << post_index_weight / post_MST_weight << std::endl;
    std::cout << "Post terrain approximate distance: " << post_approximate_distance << ", post terrain exact distance: " << post_exact_distance << ", distance error: " << post_approximate_distance / post_exact_distance - 1 << std::endl;

    std::ofstream ofs2("../output/output.txt", std::ios_base::app);
    ofs2 << post_update_time << "\t"
         << post_HGS_time << "\t"
         << post_query_time << "\t"
         << (post_memory_usage + post_complete_graph_size) / 1e6 << "\t"
         << post_complete_graph_size / 1e6 << "\t"
         << post_hierarchy_greedy_spanner_size / 1e6 << "\t"
         << post_index_edge_num << "\t"
         << post_index_weight / post_MST_weight << "\t"
         << post_approximate_distance / post_exact_distance - 1 << "\n\n";
    ofs2.close();
}

void FU_Oracle_NoEdgPru(int poi_num, geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
                        geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list,
                        int source_poi_index, int destination_poi_index,
                        double post_exact_distance, int &pre_MST_weight, int &post_MST_weight,
                        std::string write_file_header)
{
    std::vector<std::vector<std::vector<int>>> pre_face_sequence_index_list;
    std::vector<std::vector<double>> pairwise_distance_poi_to_poi;
    std::vector<std::vector<bool>> pairwise_distance_poi_to_poi_changed;
    std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> pairwise_path_poi_to_poi;
    std::vector<std::vector<double>> pairwise_distance_poi_to_vertex(poi_num, std::vector<double>(pre_mesh->vertices().size(), 0));

    double pre_construction_time = 0;
    double pre_hash_mapping_time = 0;
    double pre_memory_usage = 0;
    double pre_complete_graph_size = 0;
    int pre_index_edge_num = 0;
    double pre_index_weight = 0;

    pre_complete_graph_construction(poi_num, pre_mesh, pre_poi_list, pre_face_sequence_index_list,
                                    pairwise_distance_poi_to_poi, pairwise_distance_poi_to_poi_changed,
                                    pairwise_path_poi_to_poi, pairwise_distance_poi_to_vertex,
                                    pre_construction_time, pre_memory_usage);
    std::unordered_map<int, double> pairwise_distance_poi_to_poi_pre_map;
    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> pairwise_path_poi_to_poi_pre_map;
    get_pairwise_distance_and_path_poi_to_poi_map(pairwise_distance_poi_to_poi, pairwise_distance_poi_to_poi_pre_map,
                                                  pairwise_path_poi_to_poi, pairwise_path_poi_to_poi_pre_map,
                                                  pre_complete_graph_size, pre_index_edge_num, pre_index_weight, pre_hash_mapping_time);
    calculate_MST_weight(pairwise_distance_poi_to_poi, pre_MST_weight);

    std::cout << "Pre terrain construction time: " << (pre_construction_time + pre_hash_mapping_time) << " ms" << std::endl;
    std::cout << "Pre terrain memory usage: " << (pre_memory_usage + pre_complete_graph_size) / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index size: " << pre_complete_graph_size / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index edge number: " << pre_index_edge_num << " edges" << std::endl;
    std::cout << "Pre terrain index/MST weight: " << pre_index_weight / pre_MST_weight << std::endl;
    std::cout << std::endl;

    std::ofstream ofs1("../output/output.txt", std::ios_base::app);
    ofs1 << "== FU_Oracle_NoEdgPru ==\n";
    ofs1 << write_file_header << "\t"
         << (pre_construction_time + pre_hash_mapping_time) << "\t"
         << (pre_memory_usage + pre_complete_graph_size) / 1e6 << "\t"
         << pre_complete_graph_size / 1e6 << "\t"
         << pre_index_edge_num << "\t"
         << pre_index_weight / pre_MST_weight << "\t";
    ofs1.close();

    double post_update_time = 0;
    double post_hash_mapping_time = 0;
    double post_query_time = 0;
    double post_memory_usage = 0;
    double post_complete_graph_size = 0;
    int post_index_edge_num = 0;
    double post_index_weight = 0;
    double post_approximate_distance = 0;
    std::vector<geodesic::SurfacePoint> post_approximate_path;
    post_approximate_path.clear();

    post_complete_graph_update(poi_num, pre_mesh, pre_poi_list, post_mesh, post_poi_list,
                               pre_face_sequence_index_list, pairwise_distance_poi_to_poi,
                               pairwise_distance_poi_to_poi_changed,
                               pairwise_path_poi_to_poi, pairwise_distance_poi_to_vertex,
                               post_update_time, post_memory_usage);
    std::unordered_map<int, double> pairwise_distance_poi_to_poi_post_map;
    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> pairwise_path_poi_to_poi_post_map;
    get_pairwise_distance_and_path_poi_to_poi_map(pairwise_distance_poi_to_poi, pairwise_distance_poi_to_poi_post_map,
                                                  pairwise_path_poi_to_poi, pairwise_path_poi_to_poi_post_map,
                                                  post_complete_graph_size, post_index_edge_num, post_index_weight, post_hash_mapping_time);
    complete_graph_query(poi_num, pairwise_distance_poi_to_poi_post_map, pairwise_path_poi_to_poi_post_map, source_poi_index,
                         destination_poi_index, post_approximate_distance, post_approximate_path, post_query_time);
    calculate_MST_weight(pairwise_distance_poi_to_poi, post_MST_weight);

    std::cout << "Post terrain update time (CG): " << post_update_time << " ms" << std::endl;
    std::cout << "Post terrain construction time (hash mapping): " << post_hash_mapping_time << " ms" << std::endl;
    std::cout << "Post terrain query time: " << post_query_time << " ms" << std::endl;
    std::cout << "Post terrain memory usage: " << (post_memory_usage + post_complete_graph_size) / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index size: " << post_complete_graph_size / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain output size: " << post_complete_graph_size / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index edge number: " << post_index_edge_num << " edges" << std::endl;
    std::cout << "Post terrain index/MST weight: " << post_index_weight / post_MST_weight << std::endl;
    std::cout << "Post terrain approximate distance: " << post_approximate_distance << ", post terrain exact distance: " << post_exact_distance << ", distance error: " << post_approximate_distance / post_exact_distance - 1 << std::endl;

    std::ofstream ofs2("../output/output.txt", std::ios_base::app);
    ofs2 << post_update_time << "\t"
         << post_hash_mapping_time << "\t"
         << post_query_time << "\t"
         << (post_memory_usage + post_complete_graph_size) / 1e6 << "\t"
         << post_complete_graph_size / 1e6 << "\t"
         << post_complete_graph_size / 1e6 << "\t"
         << post_index_edge_num << "\t"
         << post_index_weight / post_MST_weight << "\t"
         << post_approximate_distance / post_exact_distance - 1 << "\n\n";
    ofs2.close();
}

void FU_Oracle_NoEffEdgPru(int poi_num, geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
                           geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list, double epsilon,
                           int source_poi_index, int destination_poi_index,
                           double post_exact_distance, int pre_MST_weight, int post_MST_weight,
                           std::string write_file_header)
{
    std::vector<std::vector<std::vector<int>>> pre_face_sequence_index_list;
    std::vector<std::vector<double>> pairwise_distance_poi_to_poi;
    std::vector<std::vector<bool>> pairwise_distance_poi_to_poi_changed;
    std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> pairwise_path_poi_to_poi;
    std::vector<std::vector<double>> pairwise_distance_poi_to_vertex(poi_num, std::vector<double>(pre_mesh->vertices().size(), 0));

    double pre_construction_time = 0;
    double pre_hash_mapping_time = 0;
    double pre_memory_usage = 0;
    double pre_complete_graph_size = 0;
    int pre_index_edge_num = 0;
    double pre_index_weight = 0;

    pre_complete_graph_construction(poi_num, pre_mesh, pre_poi_list, pre_face_sequence_index_list,
                                    pairwise_distance_poi_to_poi, pairwise_distance_poi_to_poi_changed,
                                    pairwise_path_poi_to_poi, pairwise_distance_poi_to_vertex,
                                    pre_construction_time, pre_memory_usage);
    std::unordered_map<int, double> pairwise_distance_poi_to_poi_pre_map;
    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> pairwise_path_poi_to_poi_pre_map;
    get_pairwise_distance_and_path_poi_to_poi_map(pairwise_distance_poi_to_poi, pairwise_distance_poi_to_poi_pre_map,
                                                  pairwise_path_poi_to_poi, pairwise_path_poi_to_poi_pre_map,
                                                  pre_complete_graph_size, pre_index_edge_num, pre_index_weight, pre_hash_mapping_time);

    std::cout << "Pre terrain construction time: " << (pre_construction_time + pre_hash_mapping_time) << " ms" << std::endl;
    std::cout << "Pre terrain memory usage: " << (pre_memory_usage + pre_complete_graph_size) / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index size: " << pre_complete_graph_size / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index edge number: " << pre_index_edge_num << " edges" << std::endl;
    std::cout << "Pre terrain index/MST weight: " << pre_index_weight / pre_MST_weight << std::endl;
    std::cout << std::endl;

    std::ofstream ofs1("../output/output.txt", std::ios_base::app);
    ofs1 << "== FU_Oracle_NoEffEdgPru ==\n";
    ofs1 << write_file_header << "\t"
         << (pre_construction_time + pre_hash_mapping_time) << "\t"
         << (pre_memory_usage + pre_complete_graph_size) / 1e6 << "\t"
         << pre_complete_graph_size / 1e6 << "\t"
         << pre_index_edge_num << "\t"
         << pre_index_weight / pre_MST_weight << "\t";
    ofs1.close();

    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> post_path_poi_to_poi_map;
    post_path_poi_to_poi_map.clear();
    double post_update_time = 0;
    double post_GreSpan_time = 0;
    double post_query_time = 0;
    double post_memory_usage = 0;
    double post_complete_graph_size = 0;
    double post_greedy_spanner_size = 0;
    int post_index_edge_num = 0;
    double post_index_weight = 0;
    double post_approximate_distance = 0;
    std::vector<geodesic::SurfacePoint> post_approximate_path;
    post_approximate_path.clear();
    Graph graph(poi_num);

    post_complete_graph_update(poi_num, pre_mesh, pre_poi_list, post_mesh, post_poi_list,
                               pre_face_sequence_index_list, pairwise_distance_poi_to_poi,
                               pairwise_distance_poi_to_poi_changed,
                               pairwise_path_poi_to_poi, pairwise_distance_poi_to_vertex,
                               post_update_time, post_memory_usage);
    greedy_spanner(epsilon, graph, pairwise_distance_poi_to_poi, pairwise_path_poi_to_poi, post_path_poi_to_poi_map,
                   post_greedy_spanner_size, post_index_edge_num, post_index_weight, post_GreSpan_time, post_complete_graph_size);
    spanner_query(poi_num, graph, post_path_poi_to_poi_map, source_poi_index, destination_poi_index, post_approximate_distance,
                  post_approximate_path, post_query_time);

    std::cout << "Post terrain update time (CG): " << post_update_time << " ms" << std::endl;
    std::cout << "Post terrain update time (GreSpan): " << post_GreSpan_time << " ms" << std::endl;
    std::cout << "Post terrain query time: " << post_query_time << " ms" << std::endl;
    std::cout << "Post terrain memory usage: " << (post_memory_usage + post_complete_graph_size) / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index size: " << post_complete_graph_size / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain output size: " << post_greedy_spanner_size / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index edge number: " << post_index_edge_num << " edges" << std::endl;
    std::cout << "Post terrain index/MST weight: " << post_index_weight / post_MST_weight << std::endl;
    std::cout << "Post terrain approximate distance: " << post_approximate_distance << ", post terrain exact distance: " << post_exact_distance << ", distance error: " << post_approximate_distance / post_exact_distance - 1 << std::endl;

    std::ofstream ofs2("../output/output.txt", std::ios_base::app);
    ofs2 << post_update_time << "\t"
         << post_GreSpan_time << "\t"
         << post_query_time << "\t"
         << (post_memory_usage + post_complete_graph_size) / 1e6 << "\t"
         << post_complete_graph_size / 1e6 << "\t"
         << post_greedy_spanner_size / 1e6 << "\t"
         << post_index_edge_num << "\t"
         << post_index_weight / post_MST_weight << "\t"
         << post_approximate_distance / post_exact_distance - 1 << "\n\n";
    ofs2.close();
}

void FU_Oracle(int poi_num, geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
               geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list, double epsilon,
               int source_poi_index, int destination_poi_index,
               double post_exact_distance, int pre_MST_weight, int post_MST_weight,
               std::string write_file_header)
{
    std::vector<std::vector<std::vector<int>>> pre_face_sequence_index_list;
    std::vector<std::vector<double>> pairwise_distance_poi_to_poi;
    std::vector<std::vector<bool>> pairwise_distance_poi_to_poi_changed;
    std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> pairwise_path_poi_to_poi;
    std::vector<std::vector<double>> pairwise_distance_poi_to_vertex(poi_num, std::vector<double>(pre_mesh->vertices().size(), 0));

    double pre_construction_time = 0;
    double pre_hash_mapping_time = 0;
    double pre_memory_usage = 0;
    double pre_complete_graph_size = 0;
    int pre_index_edge_num = 0;
    double pre_index_weight = 0;

    pre_complete_graph_construction(poi_num, pre_mesh, pre_poi_list, pre_face_sequence_index_list,
                                    pairwise_distance_poi_to_poi, pairwise_distance_poi_to_poi_changed,
                                    pairwise_path_poi_to_poi, pairwise_distance_poi_to_vertex,
                                    pre_construction_time, pre_memory_usage);
    std::unordered_map<int, double> pairwise_distance_poi_to_poi_pre_map;
    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> pairwise_path_poi_to_poi_pre_map;
    get_pairwise_distance_and_path_poi_to_poi_map(pairwise_distance_poi_to_poi, pairwise_distance_poi_to_poi_pre_map,
                                                  pairwise_path_poi_to_poi, pairwise_path_poi_to_poi_pre_map,
                                                  pre_complete_graph_size, pre_index_edge_num, pre_index_weight, pre_hash_mapping_time);

    std::cout << "Pre terrain construction time: " << (pre_construction_time + pre_hash_mapping_time) << " ms" << std::endl;
    std::cout << "Pre terrain memory usage: " << (pre_memory_usage + pre_complete_graph_size) / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index size: " << pre_complete_graph_size / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index edge number: " << pre_index_edge_num << " edges" << std::endl;
    std::cout << "Pre terrain index/MST weight: " << pre_index_weight / pre_MST_weight << std::endl;
    std::cout << std::endl;

    std::ofstream ofs1("../output/output.txt", std::ios_base::app);
    ofs1 << "== FU_Oracle ==\n";
    ofs1 << write_file_header << "\t"
         << (pre_construction_time + pre_hash_mapping_time) << "\t"
         << (pre_memory_usage + pre_complete_graph_size) / 1e6 << "\t"
         << pre_complete_graph_size / 1e6 << "\t"
         << pre_index_edge_num << "\t"
         << pre_index_weight / pre_MST_weight << "\t";
    ofs1.close();

    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> post_path_poi_to_poi_map;
    post_path_poi_to_poi_map.clear();
    double post_update_time = 0;
    double post_HieGreSpan_time = 0;
    double post_query_time = 0;
    double post_memory_usage = 0;
    double post_complete_graph_size = 0;
    double post_hierarchy_greedy_spanner_size = 0;
    int post_index_edge_num = 0;
    double post_index_weight = 0;
    double post_approximate_distance = 0;
    std::vector<geodesic::SurfacePoint> post_approximate_path;
    post_approximate_path.clear();
    Graph graph(poi_num);

    post_complete_graph_update(poi_num, pre_mesh, pre_poi_list, post_mesh, post_poi_list,
                               pre_face_sequence_index_list, pairwise_distance_poi_to_poi,
                               pairwise_distance_poi_to_poi_changed,
                               pairwise_path_poi_to_poi, pairwise_distance_poi_to_vertex,
                               post_update_time, post_memory_usage);
    hierarchy_greedy_spanner(epsilon, graph, pairwise_distance_poi_to_poi, pairwise_path_poi_to_poi, post_path_poi_to_poi_map,
                             post_hierarchy_greedy_spanner_size, post_index_edge_num, post_index_weight, post_HieGreSpan_time, post_complete_graph_size);
    spanner_query(poi_num, graph, post_path_poi_to_poi_map, source_poi_index, destination_poi_index, post_approximate_distance,
                  post_approximate_path, post_query_time);

    std::cout << "Post terrain update time (CG): " << post_update_time << " ms" << std::endl;
    std::cout << "Post terrain update time (HieGreSpan): " << post_HieGreSpan_time << " ms" << std::endl;
    std::cout << "Post terrain query time: " << post_query_time << " ms" << std::endl;
    std::cout << "Post terrain memory usage: " << (post_memory_usage + post_complete_graph_size) / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index size: " << post_complete_graph_size / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain output size: " << post_hierarchy_greedy_spanner_size / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index edge number: " << post_index_edge_num << " edges" << std::endl;
    std::cout << "Post terrain index/MST weight: " << post_index_weight / post_MST_weight << std::endl;
    std::cout << "Post terrain approximate distance: " << post_approximate_distance << ", post terrain exact distance: " << post_exact_distance << ", distance error: " << post_approximate_distance / post_exact_distance - 1 << std::endl;

    std::ofstream ofs2("../output/output.txt", std::ios_base::app);
    ofs2 << post_update_time << "\t"
         << post_HieGreSpan_time << "\t"
         << post_query_time << "\t"
         << (post_memory_usage + post_complete_graph_size) / 1e6 << "\t"
         << post_complete_graph_size / 1e6 << "\t"
         << post_hierarchy_greedy_spanner_size / 1e6 << "\t"
         << post_index_edge_num << "\t"
         << post_index_weight / post_MST_weight << "\t"
         << post_approximate_distance / post_exact_distance - 1 << "\n\n";
    ofs2.close();
}

void WSPD_Oracle(int poi_num, geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
                 geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list, double epsilon,
                 int source_poi_index, int destination_poi_index,
                 double post_exact_distance, int pre_MST_weight, int post_MST_weight,
                 std::string write_file_header)
{
    int pre_geo_tree_node_id = 1;
    std::vector<GeoNode *> pre_all_poi;
    std::unordered_map<int, GeoNode *> pre_geo_node_in_partition_tree_unordered_map;
    std::unordered_map<int, GeoPair *> pre_geopairs;
    std::unordered_map<int, int> pre_poi_unordered_map;
    std::vector<std::vector<double>> pairwise_distance_poi_to_vertex(poi_num, std::vector<double>(pre_mesh->vertices().size(), 0));

    double pre_construction_time = 0;
    double pre_memory_usage = 0;
    double pre_index_size = 0;
    int pre_index_edge_num = 0;
    double pre_index_weight = 0;

    pre_or_post_WSPD_Oracle_and_pre_WSPD_Oracle_Adapt_construction(
        poi_num, pre_mesh, pre_poi_list, epsilon, pre_geo_tree_node_id, pre_all_poi,
        pre_geo_node_in_partition_tree_unordered_map, pre_geopairs, pre_poi_unordered_map,
        pairwise_distance_poi_to_vertex, false, pre_construction_time, pre_memory_usage,
        pre_index_size, pre_index_edge_num, pre_index_weight);

    std::cout << "Pre terrain construction time: " << pre_construction_time << " ms" << std::endl;
    std::cout << "Pre terrain memory usage: " << pre_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index size: " << pre_index_size / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index edge number: " << pre_index_edge_num << " edges" << std::endl;
    std::cout << "Pre terrain index/MST weight: " << pre_index_weight / pre_MST_weight << std::endl;
    std::cout << std::endl;

    std::ofstream ofs1("../output/output.txt", std::ios_base::app);
    ofs1 << "== WSPD_Oracle ==\n";
    ofs1 << write_file_header << "\t"
         << pre_construction_time << "\t"
         << pre_memory_usage / 1e6 << "\t"
         << pre_index_size / 1e6 << "\t"
         << pre_index_edge_num << "\t"
         << pre_index_weight / pre_MST_weight << "\t";
    ofs1.close();

    int post_geo_tree_node_id = 1;
    std::vector<GeoNode *> post_all_poi;
    std::unordered_map<int, GeoNode *> post_geo_node_in_partition_tree_unordered_map;
    std::unordered_map<int, GeoPair *> post_geopairs;
    std::unordered_map<int, int> post_poi_unordered_map;

    double post_construction_time = 0;
    double post_query_time = 0;
    double post_memory_usage = 0;
    double post_index_size = 0;
    int post_index_edge_num = 0;
    double post_index_weight = 0;
    double post_approximate_distance = 0;
    std::vector<geodesic::SurfacePoint> post_approximate_path;
    post_approximate_path.clear();

    pre_or_post_WSPD_Oracle_and_pre_WSPD_Oracle_Adapt_construction(
        poi_num, post_mesh, post_poi_list, epsilon, post_geo_tree_node_id, post_all_poi,
        post_geo_node_in_partition_tree_unordered_map, post_geopairs, post_poi_unordered_map,
        pairwise_distance_poi_to_vertex, false, post_construction_time, post_memory_usage,
        post_index_size, post_index_edge_num, post_index_weight);
    pre_or_post_WSPD_Oracle_and_pre_WSPD_Oracle_Adapt_query(
        post_mesh, post_geo_tree_node_id, post_all_poi, post_geo_node_in_partition_tree_unordered_map,
        post_geopairs, source_poi_index, destination_poi_index,
        post_query_time, post_approximate_distance, post_approximate_path);

    std::cout << "Post terrain construction time: " << post_construction_time << " ms" << std::endl;
    std::cout << "Post terrain query time: " << post_query_time << " ms" << std::endl;
    std::cout << "Post terrain memory usage: " << post_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index size: " << post_index_size / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain output size: " << post_index_size / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index edge number: " << post_index_edge_num << " edges" << std::endl;
    std::cout << "Post terrain index/MST weight: " << post_index_weight / post_MST_weight << std::endl;
    std::cout << "Post terrain approximate distance: " << post_approximate_distance << ", post terrain exact distance: " << post_exact_distance << ", distance error: " << abs(post_approximate_distance / post_exact_distance - 1) << std::endl;

    std::ofstream ofs2("../output/output.txt", std::ios_base::app);
    ofs2 << post_construction_time << "\t"
         << 0 << "\t"
         << post_query_time << "\t"
         << post_memory_usage / 1e6 << "\t"
         << post_index_size / 1e6 << "\t"
         << post_index_size / 1e6 << "\t"
         << post_index_edge_num << "\t"
         << post_index_weight / post_MST_weight << "\t"
         << abs(post_approximate_distance / post_exact_distance - 1) << "\n\n";
    ofs2.close();
}

void WSPD_Oracle_Adapt(int poi_num, geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
                       geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list, double epsilon,
                       int source_poi_index, int destination_poi_index,
                       double post_exact_distance, int pre_MST_weight, int post_MST_weight,
                       std::string write_file_header)
{
    int geo_tree_node_id = 1;
    std::vector<GeoNode *> all_poi;
    std::unordered_map<int, GeoNode *> geo_node_in_partition_tree_unordered_map;
    std::unordered_map<int, GeoPair *> geopairs;
    std::unordered_map<int, int> poi_unordered_map;
    std::vector<std::vector<double>> pairwise_distance_poi_to_vertex(poi_num, std::vector<double>(pre_mesh->vertices().size(), 0));

    double pre_construction_time = 0;
    double pre_memory_usage = 0;
    double pre_index_size = 0;
    int pre_index_edge_num = 0;
    double pre_index_weight = 0;

    pre_or_post_WSPD_Oracle_and_pre_WSPD_Oracle_Adapt_construction(
        poi_num, pre_mesh, pre_poi_list, epsilon, geo_tree_node_id, all_poi,
        geo_node_in_partition_tree_unordered_map, geopairs, poi_unordered_map,
        pairwise_distance_poi_to_vertex, true, pre_construction_time, pre_memory_usage,
        pre_index_size, pre_index_edge_num, pre_index_weight);

    std::cout << "Pre terrain construction time: " << pre_construction_time << " ms" << std::endl;
    std::cout << "Pre terrain memory usage: " << pre_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index size: " << pre_index_size / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index edge number: " << pre_index_edge_num << " edges" << std::endl;
    std::cout << "Pre terrain index/MST weight: " << pre_index_weight / pre_MST_weight << std::endl;
    std::cout << std::endl;

    std::ofstream ofs1("../output/output.txt", std::ios_base::app);
    ofs1 << "== WSPD_Oracle_Adapt ==\n";
    ofs1 << write_file_header << "\t"
         << pre_construction_time << "\t"
         << pre_memory_usage / 1e6 << "\t"
         << pre_index_size / 1e6 << "\t"
         << pre_index_edge_num << "\t"
         << pre_index_weight / pre_MST_weight << "\t";
    ofs1.close();

    std::vector<std::vector<double>> pairwise_distance_poi_to_poi;
    std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> pairwise_path_poi_to_poi;
    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> post_path_poi_to_poi_map;
    post_path_poi_to_poi_map.clear();
    double post_update_time = 0;
    double post_HGS_time = 0;
    double post_query_time = 0;
    double post_memory_usage = 0;
    double post_full_graph_size = 0;
    double post_hierarchy_greedy_spanner_size = 0;
    int post_index_edge_num = 0;
    double post_index_weight = 0;
    double post_approximate_distance = 0;
    std::vector<geodesic::SurfacePoint> post_approximate_path;
    post_approximate_path.clear();
    Graph graph(poi_num);

    post_WSPD_Oracle_Adapt_update(
        poi_num, pre_mesh, pre_poi_list, post_mesh, post_poi_list, geo_tree_node_id, all_poi,
        geo_node_in_partition_tree_unordered_map, geopairs,
        pairwise_distance_poi_to_poi, pairwise_path_poi_to_poi,
        pairwise_distance_poi_to_vertex, post_update_time, post_memory_usage);
    hierarchy_greedy_spanner(epsilon, graph, pairwise_distance_poi_to_poi, pairwise_path_poi_to_poi, post_path_poi_to_poi_map,
                             post_hierarchy_greedy_spanner_size, post_index_edge_num, post_index_weight, post_HGS_time, post_full_graph_size);
    spanner_query(poi_num, graph, post_path_poi_to_poi_map, source_poi_index, destination_poi_index, post_approximate_distance,
                  post_approximate_path, post_query_time);

    std::cout << "Post terrain update time (WSPD_Adapt): " << post_update_time << " ms" << std::endl;
    std::cout << "Post terrain update time (HieGreSpan): " << post_HGS_time << " ms" << std::endl;
    std::cout << "Post terrain query time: " << post_query_time << " ms" << std::endl;
    std::cout << "Post terrain memory usage: " << (post_memory_usage + post_full_graph_size) / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index size: " << post_full_graph_size / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain output size: " << post_hierarchy_greedy_spanner_size / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index edge number: " << post_index_edge_num << " edges" << std::endl;
    std::cout << "Post terrain index/MST weight: " << post_index_weight / post_MST_weight << std::endl;
    std::cout << "Post terrain approximate distance: " << post_approximate_distance << ", post terrain exact distance: " << post_exact_distance << ", distance error: " << post_approximate_distance / post_exact_distance - 1 << std::endl;

    std::ofstream ofs2("../output/output.txt", std::ios_base::app);
    ofs2 << post_update_time << "\t"
         << post_HGS_time << "\t"
         << post_query_time << "\t"
         << (post_memory_usage + post_full_graph_size) / 1e6 << "\t"
         << post_full_graph_size / 1e6 << "\t"
         << post_hierarchy_greedy_spanner_size / 1e6 << "\t"
         << post_index_edge_num << "\t"
         << post_index_weight / post_MST_weight << "\t"
         << post_approximate_distance / post_exact_distance - 1 << "\n\n";
    ofs2.close();
}

void EAR_Oracle(int poi_num, geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
                geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list, double epsilon,
                int source_poi_index, int destination_poi_index,
                double post_exact_distance, int pre_MST_weight, int post_MST_weight,
                std::string write_file_header)
{
    int sqrt_num_of_box = 2;
    std::vector<int> pre_highway_node_list;
    std::unordered_map<int, std::unordered_map<int, int>> pre_highway_node_id_with_box_id_map;
    int pre_geo_tree_node_id = 1;
    std::vector<GeoNode *> pre_all_highway_node;
    std::unordered_map<int, GeoNode *> pre_geo_node_in_partition_tree_unordered_map;
    std::unordered_map<int, GeoPair *> pre_geopairs;
    std::unordered_map<int, int> pre_highway_node_unordered_map;
    std::vector<std::vector<std::vector<int>>> face_sequence_index_list;
    std::vector<std::vector<double>> pairwise_distance_poi_to_poi;
    std::vector<std::vector<bool>> pairwise_distance_poi_to_poi_changed;
    std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> pairwise_path_poi_to_poi;
    std::vector<std::vector<double>> pairwise_distance_poi_to_vertex(poi_num, std::vector<double>(pre_mesh->vertices().size(), 0));
    std::unordered_map<int, double> pre_distance_poi_to_highway_node_map;
    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> pre_path_poi_to_highway_node_map;

    double pre_construction_time = 0;
    double pre_memory_usage = 0;
    double pre_index_size = 0;
    int pre_index_edge_num = 0;
    double pre_index_weight = 0;

    pre_or_post_EAR_Oracle_and_pre_EAR_Oracle_Adapt_construction(
        sqrt_num_of_box, poi_num, pre_mesh, pre_poi_list, epsilon, pre_highway_node_list,
        pre_highway_node_id_with_box_id_map, pre_geo_tree_node_id, pre_all_highway_node,
        pre_geo_node_in_partition_tree_unordered_map, pre_geopairs, pre_highway_node_unordered_map,
        face_sequence_index_list, pairwise_distance_poi_to_poi,
        pairwise_distance_poi_to_poi_changed, pairwise_path_poi_to_poi,
        pairwise_distance_poi_to_vertex, pre_distance_poi_to_highway_node_map,
        pre_path_poi_to_highway_node_map, false, pre_construction_time, pre_memory_usage,
        pre_index_size, pre_index_edge_num, pre_index_weight);

    std::cout << "Pre terrain construction time: " << pre_construction_time << " ms" << std::endl;
    std::cout << "Pre terrain memory usage: " << pre_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index size: " << pre_index_size / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index edge number: " << pre_index_edge_num << " edges" << std::endl;
    std::cout << "Pre terrain index/MST weight: " << pre_index_weight / pre_MST_weight << std::endl;
    std::cout << std::endl;

    std::ofstream ofs1("../output/output.txt", std::ios_base::app);
    ofs1 << "== EAR_Oracle ==\n";
    ofs1 << write_file_header << "\t"
         << pre_construction_time << "\t"
         << pre_memory_usage / 1e6 << "\t"
         << pre_index_size / 1e6 << "\t"
         << pre_index_edge_num << "\t"
         << pre_index_weight / pre_MST_weight << "\t";
    ofs1.close();

    std::vector<int> post_highway_node_list;
    std::unordered_map<int, std::unordered_map<int, int>> post_highway_node_id_with_box_id_map;
    int post_geo_tree_node_id = 1;
    std::vector<GeoNode *> post_all_highway_node;
    std::unordered_map<int, GeoNode *> post_geo_node_in_partition_tree_unordered_map;
    std::unordered_map<int, GeoPair *> post_geopairs;
    std::unordered_map<int, int> post_highway_node_unordered_map;
    std::unordered_map<int, double> post_distance_poi_to_highway_node_map;
    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> post_path_poi_to_highway_node_map;

    double post_construction_time = 0;
    double post_query_time = 0;
    double post_memory_usage = 0;
    double post_index_size = 0;
    int post_index_edge_num = 0;
    double post_index_weight = 0;
    double post_approximate_distance = 0;
    std::vector<geodesic::SurfacePoint> post_approximate_path;
    post_approximate_path.clear();

    pre_or_post_EAR_Oracle_and_pre_EAR_Oracle_Adapt_construction(
        sqrt_num_of_box, poi_num, post_mesh, post_poi_list, epsilon, post_highway_node_list,
        post_highway_node_id_with_box_id_map, post_geo_tree_node_id, post_all_highway_node,
        post_geo_node_in_partition_tree_unordered_map, post_geopairs, post_highway_node_unordered_map,
        face_sequence_index_list, pairwise_distance_poi_to_poi,
        pairwise_distance_poi_to_poi_changed, pairwise_path_poi_to_poi,
        pairwise_distance_poi_to_vertex, post_distance_poi_to_highway_node_map,
        post_path_poi_to_highway_node_map, false, post_construction_time, post_memory_usage,
        post_index_size, post_index_edge_num, post_index_weight);
    pre_or_post_EAR_Oracle_and_pre_EAR_Oracle_Adapt_query(
        post_mesh, post_poi_list, sqrt_num_of_box, post_geo_tree_node_id, post_highway_node_id_with_box_id_map,
        post_all_highway_node, post_geo_node_in_partition_tree_unordered_map,
        post_geopairs, post_distance_poi_to_highway_node_map, post_path_poi_to_highway_node_map,
        source_poi_index, destination_poi_index, post_query_time, post_approximate_distance, post_approximate_path);

    std::cout << "Post terrain construction time: " << post_construction_time << " ms" << std::endl;
    std::cout << "Post terrain query time: " << post_query_time << " ms" << std::endl;
    std::cout << "Post terrain memory usage: " << post_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index size: " << post_index_size / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain output size: " << post_index_size / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index edge number: " << post_index_edge_num << " edges" << std::endl;
    std::cout << "Post terrain index/MST weight: " << post_index_weight / post_MST_weight << std::endl;
    std::cout << "Post terrain approximate distance: " << post_approximate_distance << ", post terrain exact distance: " << post_exact_distance << ", distance error: " << abs(post_approximate_distance / post_exact_distance - 1) << std::endl;

    std::ofstream ofs2("../output/output.txt", std::ios_base::app);
    ofs2 << post_construction_time << "\t"
         << 0 << "\t"
         << post_query_time << "\t"
         << post_memory_usage / 1e6 << "\t"
         << post_index_size / 1e6 << "\t"
         << post_index_size / 1e6 << "\t"
         << post_index_edge_num << "\t"
         << post_index_weight / post_MST_weight << "\t"
         << abs(post_approximate_distance / post_exact_distance - 1) << "\n\n";
    ofs2.close();
}

void EAR_Oracle_Adapt(int poi_num, geodesic::Mesh *pre_mesh, std::vector<int> &pre_poi_list,
                      geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list, double epsilon,
                      int source_poi_index, int destination_poi_index,
                      double post_exact_distance, int pre_MST_weight, int post_MST_weight,
                      std::string write_file_header)
{
    int sqrt_num_of_box = 2;
    std::vector<int> highway_node_list;
    std::unordered_map<int, std::unordered_map<int, int>> highway_node_id_with_box_id_map;
    int geo_tree_node_id = 1;
    std::vector<GeoNode *> all_highway_node;
    std::unordered_map<int, GeoNode *> geo_node_in_partition_tree_unordered_map;
    std::unordered_map<int, GeoPair *> geopairs;
    std::unordered_map<int, int> highway_node_unordered_map;
    std::vector<std::vector<std::vector<int>>> pre_face_sequence_index_list;
    std::vector<std::vector<double>> pairwise_distance_poi_to_poi;
    std::vector<std::vector<bool>> pairwise_distance_poi_to_poi_changed;
    std::vector<std::vector<std::vector<geodesic::SurfacePoint>>> pairwise_path_poi_to_poi;
    std::vector<std::vector<double>> pairwise_distance_poi_to_vertex(poi_num, std::vector<double>(pre_mesh->vertices().size(), 0));
    std::unordered_map<int, double> distance_poi_to_highway_node_map;
    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> path_poi_to_highway_node_map;

    double pre_construction_time = 0;
    double pre_memory_usage = 0;
    double pre_index_size = 0;
    int pre_index_edge_num = 0;
    double pre_index_weight = 0;

    pre_or_post_EAR_Oracle_and_pre_EAR_Oracle_Adapt_construction(
        sqrt_num_of_box, poi_num, pre_mesh, pre_poi_list, epsilon, highway_node_list,
        highway_node_id_with_box_id_map, geo_tree_node_id, all_highway_node,
        geo_node_in_partition_tree_unordered_map, geopairs, highway_node_unordered_map,
        pre_face_sequence_index_list, pairwise_distance_poi_to_poi, pairwise_distance_poi_to_poi_changed,
        pairwise_path_poi_to_poi, pairwise_distance_poi_to_vertex, distance_poi_to_highway_node_map,
        path_poi_to_highway_node_map, true, pre_construction_time, pre_memory_usage,
        pre_index_size, pre_index_edge_num, pre_index_weight);

    std::cout << "Pre terrain construction time: " << pre_construction_time << " ms" << std::endl;
    std::cout << "Pre terrain memory usage: " << pre_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index size: " << pre_index_size / 1e6 << " MB" << std::endl;
    std::cout << "Pre terrain index edge number: " << pre_index_edge_num << " edges" << std::endl;
    std::cout << "Pre terrain index/MST weight: " << pre_index_weight / pre_MST_weight << std::endl;
    std::cout << std::endl;

    std::ofstream ofs1("../output/output.txt", std::ios_base::app);
    ofs1 << "== EAR_Oracle_Adapt ==\n";
    ofs1 << write_file_header << "\t"
         << pre_construction_time << "\t"
         << pre_memory_usage / 1e6 << "\t"
         << pre_index_size / 1e6 << "\t"
         << pre_index_edge_num << "\t"
         << pre_index_weight / pre_MST_weight << "\t";
    ofs1.close();

    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> post_path_poi_to_poi_map;
    post_path_poi_to_poi_map.clear();
    double post_update_time = 0;
    double post_HieGreSpan_time = 0;
    double post_query_time = 0;
    double post_memory_usage = 0;
    double post_full_graph_size = 0;
    double post_hierarchy_greedy_spanner_size = 0;
    int post_index_edge_num = 0;
    double post_index_weight = 0;
    double post_approximate_distance = 0;
    std::vector<geodesic::SurfacePoint> post_approximate_path;
    post_approximate_path.clear();
    Graph graph(poi_num);

    post_EAR_Oracle_Adapt_update(poi_num, sqrt_num_of_box, pre_mesh, pre_poi_list, post_mesh, post_poi_list, geo_tree_node_id, highway_node_id_with_box_id_map, all_highway_node,
                                 geo_node_in_partition_tree_unordered_map, geopairs, distance_poi_to_highway_node_map, path_poi_to_highway_node_map, pre_face_sequence_index_list,
                                 pairwise_distance_poi_to_poi, pairwise_distance_poi_to_poi_changed,
                                 pairwise_path_poi_to_poi, pairwise_distance_poi_to_vertex, post_update_time, post_memory_usage);
    hierarchy_greedy_spanner(epsilon, graph, pairwise_distance_poi_to_poi, pairwise_path_poi_to_poi, post_path_poi_to_poi_map,
                             post_hierarchy_greedy_spanner_size, post_index_edge_num, post_index_weight, post_HieGreSpan_time, post_full_graph_size);
    spanner_query(poi_num, graph, post_path_poi_to_poi_map, source_poi_index, destination_poi_index, post_approximate_distance,
                  post_approximate_path, post_query_time);

    std::cout << "Post terrain update time (EAR_Adapt): " << post_update_time << " ms" << std::endl;
    std::cout << "Post terrain update time (HieGreSpan): " << post_HieGreSpan_time << " ms" << std::endl;
    std::cout << "Post terrain query time: " << post_query_time << " ms" << std::endl;
    std::cout << "Post terrain memory usage: " << (post_memory_usage + post_full_graph_size) / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index size: " << post_full_graph_size / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain output size: " << post_hierarchy_greedy_spanner_size / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain index edge number: " << post_index_edge_num << " edges" << std::endl;
    std::cout << "Post terrain index/MST weight: " << post_index_weight / post_MST_weight << std::endl;
    std::cout << "Post terrain approximate distance: " << post_approximate_distance << ", post terrain exact distance: " << post_exact_distance << ", distance error: " << post_approximate_distance / post_exact_distance - 1 << std::endl;

    std::ofstream ofs2("../output/output.txt", std::ios_base::app);
    ofs2 << post_update_time << "\t"
         << post_HieGreSpan_time << "\t"
         << post_query_time << "\t"
         << (post_memory_usage + post_full_graph_size) / 1e6 << "\t"
         << post_full_graph_size / 1e6 << "\t"
         << post_hierarchy_greedy_spanner_size / 1e6 << "\t"
         << post_index_edge_num << "\t"
         << post_index_weight / post_MST_weight << "\t"
         << post_approximate_distance / post_exact_distance - 1 << "\n\n";
    ofs2.close();
}

void CH_Fly_Algo(geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list,
                 int source_poi_index, int destination_poi_index,
                 double post_exact_distance, std::string write_file_header)
{
    double post_query_time = 0;
    double post_memory_usage = 0;
    double post_approximate_distance = 0;
    std::vector<geodesic::SurfacePoint> post_approximate_path;
    post_approximate_path.clear();

    pre_or_post_CH_Fly_Algo_query(post_mesh, post_poi_list, source_poi_index,
                                  destination_poi_index, post_query_time, post_memory_usage,
                                  post_approximate_distance, post_approximate_path);

    std::cout << "Post terrain query time: " << post_query_time << " ms" << std::endl;
    std::cout << "Post terrain memory usage: " << post_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain approximate distance: " << post_approximate_distance << ", post terrain exact distance: " << post_exact_distance << ", distance error: " << post_approximate_distance / post_exact_distance - 1 << std::endl;

    std::ofstream ofs1("../output/output.txt", std::ios_base::app);
    ofs1 << "== CH_Fly_Algo ==\n";
    ofs1 << write_file_header << "\t"
         << 0 << "\t"
         << 0 << "\t"
         << 0 << "\t"
         << 0 << "\t"
         << 0 << "\t"
         << 0 << "\t"

         << 0 << "\t"
         << 0 << "\t"
         << post_query_time << "\t"
         << post_memory_usage / 1e6 << "\t"
         << 0 << "\t"
         << 0 << "\t"
         << 0 << "\t"
         << 0 << "\t"
         << post_approximate_distance / post_exact_distance - 1 << "\n\n";
    ofs1.close();
}

void K_Fly_Algo(geodesic::Mesh *post_mesh, std::vector<int> &post_poi_list, double epsilon,
                int source_poi_index, int destination_poi_index,
                double post_exact_distance, std::string write_file_header)
{
    double post_query_time = 0;
    double post_memory_usage = 0;
    double post_approximate_distance = 0;
    std::vector<geodesic::SurfacePoint> post_approximate_path;
    post_approximate_path.clear();

    pre_or_post_K_Fly_Algo_query(post_mesh, post_poi_list, epsilon, source_poi_index,
                                 destination_poi_index, post_query_time, post_memory_usage,
                                 post_approximate_distance, post_approximate_path);

    std::cout << "Post terrain query time: " << post_query_time << " ms" << std::endl;
    std::cout << "Post terrain memory usage: " << post_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Post terrain approximate distance: " << post_approximate_distance << ", post terrain exact distance: " << post_exact_distance << ", distance error: " << post_approximate_distance / post_exact_distance - 1 << std::endl;

    std::ofstream ofs1("../output/output.txt", std::ios_base::app);
    ofs1 << "== K_Fly_Algo ==\n";
    ofs1 << write_file_header << "\t"
         << 0 << "\t"
         << 0 << "\t"
         << 0 << "\t"
         << 0 << "\t"
         << 0 << "\t"
         << 0 << "\t"

         << 0 << "\t"
         << 0 << "\t"
         << post_query_time << "\t"
         << post_memory_usage / 1e6 << "\t"
         << 0 << "\t"
         << 0 << "\t"
         << 0 << "\t"
         << 0 << "\t"
         << post_approximate_distance / post_exact_distance - 1 << "\n\n";
    ofs1.close();
}