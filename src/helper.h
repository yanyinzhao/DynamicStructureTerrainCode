#include <stdio.h>
#include <stdlib.h>
#include <list>
#include <cmath>
#include <iostream>
#include <fstream>
#include <bits/stdc++.h>
#include "geodesic_algorithm_exact.h"
#include "geodesic_algorithm_subdivision.h"

class GeoNode
{
public:
    int id;
    double radius;
    int index;
    GeoNode *parent;
    std::list<GeoNode *> children;
    GeoNode *covered_by;
    std::vector<std::pair<int, GeoNode *>> covers;
    GeoNode()
    {
        this->parent = NULL;
        this->children.clear();
        this->covered_by = NULL;
        this->covers.clear();
    }
    GeoNode(const int index, const double r)
    {
        this->index = index;
        this->radius = r;
        this->parent = NULL;
        this->children.clear();
        this->covered_by = NULL;
        this->covers.clear();
    }
    GeoNode(const int id, const int index, const double r)
    {
        this->id = id;
        this->index = index;
        this->radius = r;
        this->parent = NULL;
        this->children.clear();
        this->covered_by = NULL;
        this->covers.clear();
    }
    GeoNode(const GeoNode &n)
    {
        this->id = n.id;
        this->index = n.index;
        this->radius = n.radius;
        this->parent = n.parent;
        this->children.assign(n.children.begin(), n.children.end());
        this->covered_by = n.covered_by;
        this->covers.assign(n.covers.begin(), n.covers.end());
    }
    void set_id(int id)
    {
        this->id = id;
    }
    void set_parent(GeoNode *p)
    {
        this->parent = p;
        p->children.push_back(this);
    }
    void set_covered_by(GeoNode *p)
    {
        this->covered_by = p;
    }
    void set_covers(GeoNode *p)
    {
        std::pair<int, GeoNode *> q(this->index, this);
        p->covers.push_back(q);
    }
    ~GeoNode(){};
};

class GeoPair
{
public:
    GeoNode *node1;
    GeoNode *node2;
    double distance;
    std::vector<geodesic::SurfacePoint> path;
    std::vector<int> face_sequence_index_list;
};

void hash_function_two_keys_to_one_key(int row_or_column, int i, int j, int &i_j)
{
    i_j = i * row_or_column + j;
}

void hash_function_one_key_to_two_keys(int row_or_column, int i_j, int &i, int &j)
{
    i = i_j / row_or_column;
    j = i_j % row_or_column;
}

void get_face_sequence(std::vector<geodesic::SurfacePoint> &path,
                       std::vector<int> &face_sequence_index_list)
{
    face_sequence_index_list.clear();
    for (unsigned i = 0; i < path.size() - 1; ++i)
    {
        geodesic::SurfacePoint &s_1 = path[i];
        geodesic::SurfacePoint &s_2 = path[i + 1];
        for (unsigned j = 0; j < s_2.base_element()->adjacent_faces().size(); ++j)
        {
            for (unsigned k = 0; k < s_1.base_element()->adjacent_faces().size(); ++k)
            {
                if (s_2.base_element()->adjacent_faces()[j]->id() == s_1.base_element()->adjacent_faces()[k]->id())
                {
                    face_sequence_index_list.push_back(s_2.base_element()->adjacent_faces()[j]->id());
                }
            }
        }
    }
}

void build_level(int &geo_tree_node_id, geodesic::Mesh *mesh, int depth,
                 stx::btree<int, GeoNode *> &pois_B_tree, stx::btree<int, GeoNode *> &pois_as_center_each_parent_layer,
                 geodesic::GeodesicAlgorithmExact &algorithm)
{
    std::vector<std::pair<int, GeoNode *>> pois_as_center_each_current_layer;
    pois_as_center_each_current_layer.clear();

    for (stx::btree<int, GeoNode *>::iterator ite = pois_as_center_each_parent_layer.begin(); ite != pois_as_center_each_parent_layer.end(); ite++)
    {
        std::vector<std::pair<int, GeoNode *>> current_parent_covers_but_remained_pois = (*ite).second->covers;

        GeoNode *n = new GeoNode(geo_tree_node_id, (*ite).second->index, ((*ite).second)->radius / 2.0);
        geo_tree_node_id++;
        n->set_parent((*ite).second);

        std::pair<int, GeoNode *> m(n->index, n);
        pois_as_center_each_current_layer.push_back(m);

        geodesic::SurfacePoint source(&mesh->vertices()[n->index]);
        std::vector<geodesic::SurfacePoint> one_source_list(1, source);
        double const distance_limit = n->radius;
        algorithm.propagate(one_source_list, distance_limit * 1.00001);

        auto bite = current_parent_covers_but_remained_pois.begin();

        while (bite != current_parent_covers_but_remained_pois.end())
        {
            double distance;
            geodesic::SurfacePoint p(&mesh->vertices()[(*bite).first]);
            algorithm.best_source(p, distance);
            if (distance <= distance_limit)
            {
                (*bite).second->set_covered_by(n);
                (*bite).second->set_covers(n);
                current_parent_covers_but_remained_pois.erase(bite);
            }
            else
            {
                bite++;
            }
        }

        while (!current_parent_covers_but_remained_pois.empty())
        {
            GeoNode *a = new GeoNode(geo_tree_node_id, (*current_parent_covers_but_remained_pois.begin()).second->index, (*current_parent_covers_but_remained_pois.begin()).second->covered_by->radius / 2.0);
            geo_tree_node_id++;
            a->set_parent((*ite).second);

            current_parent_covers_but_remained_pois.erase(current_parent_covers_but_remained_pois.begin());

            std::pair<int, GeoNode *> b(a->index, a);
            pois_as_center_each_current_layer.push_back(b);

            geodesic::SurfacePoint source(&mesh->vertices()[a->index]);
            std::vector<geodesic::SurfacePoint> one_source_list(1, source);
            double const distance_limit = b.second->radius;
            algorithm.propagate(one_source_list, distance_limit * 1.00001);

            auto bite = current_parent_covers_but_remained_pois.begin();

            while (bite != current_parent_covers_but_remained_pois.end())
            {
                double distance;
                geodesic::SurfacePoint p(&mesh->vertices()[(*bite).first]);
                algorithm.best_source(p, distance);
                if (distance <= distance_limit)
                {
                    (*bite).second->set_covered_by(a);
                    (*bite).second->set_covers(a);
                    current_parent_covers_but_remained_pois.erase(bite);
                }
                else
                {
                    bite++;
                }
            }
        }
    }

    pois_as_center_each_parent_layer.clear();
    pois_as_center_each_parent_layer = stx::btree<int, GeoNode *>(pois_as_center_each_current_layer.begin(), pois_as_center_each_current_layer.end());
}

void build_geo_tree(int &geo_tree_node_id, geodesic::Mesh *mesh, GeoNode &root_geo, int poi_num,
                    stx::btree<int, GeoNode *> &pois_B_tree, stx::btree<int, GeoNode *> &pois_as_center_each_parent_layer,
                    geodesic::GeodesicAlgorithmExact &algorithm)
{
    int depth = 0;
    pois_as_center_each_parent_layer.clear();
    pois_as_center_each_parent_layer.insert(root_geo.index, &root_geo);

    stx::btree<int, GeoNode *> remained_pois = stx::btree<int, GeoNode *>(pois_B_tree);

    remained_pois.erase(root_geo.index);

    geodesic::SurfacePoint source(&mesh->vertices()[root_geo.index]);
    std::vector<geodesic::SurfacePoint> one_source_list(1, source);
    double const distance_limit = root_geo.radius;
    algorithm.propagate(one_source_list, distance_limit * 1.00001);
    for (stx::btree<int, GeoNode *>::iterator bite = remained_pois.begin(); bite != remained_pois.end(); bite++)
    {
        double distance;
        geodesic::SurfacePoint p(&mesh->vertices()[(*bite).first]);
        algorithm.best_source(p, distance);
        if (distance <= distance_limit)
        {
            (*bite).second->set_covered_by(&root_geo);
            (*bite).second->set_covers(&root_geo);
        }
    }

    while (pois_as_center_each_parent_layer.size() != poi_num)
    {
        build_level(geo_tree_node_id, mesh, depth, pois_B_tree, pois_as_center_each_parent_layer, algorithm);
        depth++;
    }
}

void count_single_node(GeoNode &n, std::vector<GeoNode *> &partition_tree_to_compressed_partition_tree_to_be_removed_nodes)
{
    if (!n.children.empty())
    {
        for (std::list<GeoNode *>::iterator ite = n.children.begin(); ite != n.children.end(); ite++)
        {
            if (n.children.size() == 1 && n.parent != NULL)
            {
                partition_tree_to_compressed_partition_tree_to_be_removed_nodes.push_back(&n);
            }
            count_single_node(**ite, partition_tree_to_compressed_partition_tree_to_be_removed_nodes);
        }
    }
}

void remove_single_node(GeoNode &n, std::vector<GeoNode *> &partition_tree_to_compressed_partition_tree_to_be_removed_nodes)
{
    partition_tree_to_compressed_partition_tree_to_be_removed_nodes.clear();
    count_single_node(n, partition_tree_to_compressed_partition_tree_to_be_removed_nodes);

    for (int i = 0; i < partition_tree_to_compressed_partition_tree_to_be_removed_nodes.size(); i++)
    {
        partition_tree_to_compressed_partition_tree_to_be_removed_nodes[i]->parent->children.erase(std::remove(partition_tree_to_compressed_partition_tree_to_be_removed_nodes[i]->parent->children.begin(),
                                                                                                               partition_tree_to_compressed_partition_tree_to_be_removed_nodes[i]->parent->children.end(),
                                                                                                               partition_tree_to_compressed_partition_tree_to_be_removed_nodes[i]),
                                                                                                   partition_tree_to_compressed_partition_tree_to_be_removed_nodes[i]->parent->children.end());
        partition_tree_to_compressed_partition_tree_to_be_removed_nodes[i]->children.front()->set_parent(partition_tree_to_compressed_partition_tree_to_be_removed_nodes[i]->parent);
    }
}

void set_bottom_children_raduis(GeoNode &n, std::unordered_map<int, GeoNode *> &geo_node_in_partition_tree_unordered_map)
{
    if (n.children.empty())
    {
        n.radius = 0;
        geo_node_in_partition_tree_unordered_map[n.index] = &n;
    }
    else
    {
        for (std::list<GeoNode *>::iterator ite = n.children.begin(); ite != n.children.end(); ite++)
        {
            set_bottom_children_raduis(**ite, geo_node_in_partition_tree_unordered_map);
        }
    }
}

void partition_tree_to_compressed_partition_tree(GeoNode &n, std::vector<GeoNode *> &partition_tree_to_compressed_partition_tree_to_be_removed_nodes,
                                                 std::unordered_map<int, GeoNode *> &geo_node_in_partition_tree_unordered_map)
{
    remove_single_node(n, partition_tree_to_compressed_partition_tree_to_be_removed_nodes);
    set_bottom_children_raduis(n, geo_node_in_partition_tree_unordered_map);
}

void print_geo_tree(GeoNode &n)
{
    if (!n.children.empty())
    {
        for (std::list<GeoNode *>::iterator ite = n.children.begin(); ite != n.children.end(); ite++)
        {
            std::cout << "children index: " << (*ite)->index << ", children raduis: " << (*ite)->radius << ", parent index: " << n.index << ", parent raduis: " << n.radius << std::endl;
            print_geo_tree(**ite);
        }
    }
}

void print_geo_pair(std::unordered_map<int, GeoPair *> &geopairs)
{
    FILE *fp = fopen("geopairs.txt", "w");
    for (std::unordered_map<int, GeoPair *>::iterator gite = geopairs.begin(); gite != geopairs.end(); gite++)
    {
        fprintf(fp, "index %d index %d distance %f\n", (*gite).second->node1->index, (*gite).second->node2->index, (*gite).second->distance);
    }
    fclose(fp);
}

void delete_geo_tree(GeoNode &n)
{
    if (n.children.empty())
    {
        delete &n;
        return;
    }
    else
    {
        for (std::list<GeoNode *>::iterator ite = n.children.begin(); ite != n.children.end(); ite++)
            delete_geo_tree(**ite);
        return;
    }
}

int randn(int k)
{
    int a;
    a = rand() % k;
    return a;
}

double max(double x, double y)
{
    if (x > y)
        return x;
    return y;
}

void generate_geo_pair(int geo_tree_node_id, int &WSPD_Oracle_edge_num, double &WSPD_Oracle_weight,
                       geodesic::Mesh *mesh, GeoNode &x, GeoNode &y,
                       geodesic::GeodesicAlgorithmExact &algorithm, double epsilon,
                       std::unordered_map<int, GeoPair *> &geopairs,
                       std::unordered_map<int, int> &poi_unordered_map,
                       std::unordered_map<int, int> &geo_pair_unordered_map,
                       std::unordered_map<int, double> &pairwise_distance_unordered_map,
                       std::unordered_map<int, std::vector<geodesic::SurfacePoint>> &pairwise_path_unordered_map,
                       int &pairwise_path_poi_to_poi_size, int &face_sequence_index_list_size)
{
    int x_in_geo_node_id = x.id;
    int y_in_geo_node_id = y.id;
    int x_y_in_geo_node_id;
    if (x_in_geo_node_id > y_in_geo_node_id)
    {
        int temp1 = y_in_geo_node_id;
        y_in_geo_node_id = x_in_geo_node_id;
        x_in_geo_node_id = temp1;
    }
    hash_function_two_keys_to_one_key(geo_tree_node_id, x_in_geo_node_id, y_in_geo_node_id, x_y_in_geo_node_id);

    if (geo_pair_unordered_map.count(x_y_in_geo_node_id) == 0)
    {
        geo_pair_unordered_map[x_y_in_geo_node_id] = 1; // avoid storing a pair for two times

        assert(poi_unordered_map.count(x.index) != 0 && poi_unordered_map.count(y.index) != 0);
        int x_in_poi_list_for_pairwise_distance_index = poi_unordered_map[x.index];
        int y_in_poi_list_for_pairwise_distance_index = poi_unordered_map[y.index];
        int x_y_in_poi_list_for_pairwise_distance_index;
        if (x_in_poi_list_for_pairwise_distance_index > y_in_poi_list_for_pairwise_distance_index)
        {
            int temp2 = y_in_poi_list_for_pairwise_distance_index;
            y_in_poi_list_for_pairwise_distance_index = x_in_poi_list_for_pairwise_distance_index;
            x_in_poi_list_for_pairwise_distance_index = temp2;
        }
        hash_function_two_keys_to_one_key(poi_unordered_map.size(), x_in_poi_list_for_pairwise_distance_index, y_in_poi_list_for_pairwise_distance_index, x_y_in_poi_list_for_pairwise_distance_index);

        if (pairwise_distance_unordered_map.count(x_y_in_poi_list_for_pairwise_distance_index) == 0)
        {
            pairwise_distance_unordered_map[x_y_in_poi_list_for_pairwise_distance_index] = 0;
            if (&x != &y)
            {
                double const distance_limit = geodesic::GEODESIC_INF;
                geodesic::SurfacePoint source(&mesh->vertices()[x.index]);
                std::vector<geodesic::SurfacePoint> one_source_list(1, source);
                geodesic::SurfacePoint destination(&mesh->vertices()[y.index]);
                std::vector<geodesic::SurfacePoint> one_destination_list(1, destination);
                algorithm.propagate(one_source_list, distance_limit, &one_destination_list);
                std::vector<geodesic::SurfacePoint> path;
                algorithm.trace_back(destination, path);
                pairwise_path_unordered_map[x_y_in_poi_list_for_pairwise_distance_index] = path;
                pairwise_distance_unordered_map[x_y_in_poi_list_for_pairwise_distance_index] = length(path);
            }
        }

        double distancexy = pairwise_distance_unordered_map[x_y_in_poi_list_for_pairwise_distance_index];
        std::vector<geodesic::SurfacePoint> pathxy = pairwise_path_unordered_map[x_y_in_poi_list_for_pairwise_distance_index];

        if (x.radius == 0 && y.radius == 0)
        {
            if (&x == &y)
            {
                return;
            }
            if (x.index == y.index)
            {
                return;
            }
        }
        if (distancexy >= (2.0 / epsilon + 2.0) * max(x.radius, y.radius))
        {
            std::vector<int> face_sequence_index_listxy;
            face_sequence_index_listxy.clear();
            get_face_sequence(pathxy, face_sequence_index_listxy);
            WSPD_Oracle_edge_num++;
            GeoPair *nodepair = new GeoPair();
            nodepair->node1 = &x;
            nodepair->node2 = &y;
            nodepair->distance = distancexy;
            nodepair->path = pathxy;
            nodepair->face_sequence_index_list = face_sequence_index_listxy;
            WSPD_Oracle_weight += distancexy;
            pairwise_path_poi_to_poi_size += pathxy.size();
            face_sequence_index_list_size += face_sequence_index_listxy.size();
            int x_in_geo_node_id_for_geo_pair = x.id;
            int y_in_geo_node_id_for_geo_pair = y.id;
            int x_y_in_geo_node_id_for_geo_pair;
            if (x_in_geo_node_id_for_geo_pair > y_in_geo_node_id_for_geo_pair)
            {
                int temp3 = y_in_geo_node_id_for_geo_pair;
                y_in_geo_node_id_for_geo_pair = x_in_geo_node_id_for_geo_pair;
                x_in_geo_node_id_for_geo_pair = temp3;
            }
            hash_function_two_keys_to_one_key(geo_tree_node_id, x_in_geo_node_id_for_geo_pair, y_in_geo_node_id_for_geo_pair, x_y_in_geo_node_id_for_geo_pair);
            geopairs[x_y_in_geo_node_id_for_geo_pair] = nodepair;
        }
        else
        {
            if (x.radius > y.radius)
            {
                for (std::list<GeoNode *>::iterator ite = x.children.begin(); ite != x.children.end(); ite++)
                {
                    generate_geo_pair(geo_tree_node_id, WSPD_Oracle_edge_num, WSPD_Oracle_weight, mesh, (**ite), y, algorithm, epsilon, geopairs, poi_unordered_map, geo_pair_unordered_map, pairwise_distance_unordered_map, pairwise_path_unordered_map, pairwise_path_poi_to_poi_size, face_sequence_index_list_size);
                }
            }
            else
            {
                for (std::list<GeoNode *>::iterator jte = y.children.begin(); jte != y.children.end(); jte++)
                {
                    generate_geo_pair(geo_tree_node_id, WSPD_Oracle_edge_num, WSPD_Oracle_weight, mesh, x, (**jte), algorithm, epsilon, geopairs, poi_unordered_map, geo_pair_unordered_map, pairwise_distance_unordered_map, pairwise_path_unordered_map, pairwise_path_poi_to_poi_size, face_sequence_index_list_size);
                }
            }
        }
    }
}

double one_path_query_geo(int geo_tree_node_id, GeoNode &x, GeoNode &y,
                          int &returned_source_neighbour_index, int &returned_destination_neighbour_index,
                          std::unordered_map<int, GeoPair *> &geopairs,
                          std::unordered_map<int, int> &poi_unordered_map,
                          std::vector<geodesic::SurfacePoint> &approximate_path,
                          std::vector<int> &face_sequence_index_list)
{
    GeoNode *p;

    if (&x == &y)
    {
        return 0;
    }

    int x_in_geo_node_id_for_geo_pair = x.id;
    int y_in_geo_node_id_for_geo_pair = y.id;
    int x_y_in_geo_node_id_for_geo_pair;
    if (x_in_geo_node_id_for_geo_pair > y_in_geo_node_id_for_geo_pair)
    {
        int temp1 = y_in_geo_node_id_for_geo_pair;
        y_in_geo_node_id_for_geo_pair = x_in_geo_node_id_for_geo_pair;
        x_in_geo_node_id_for_geo_pair = temp1;
    }
    hash_function_two_keys_to_one_key(geo_tree_node_id, x_in_geo_node_id_for_geo_pair, y_in_geo_node_id_for_geo_pair, x_y_in_geo_node_id_for_geo_pair);

    if (geopairs.count(x_y_in_geo_node_id_for_geo_pair) != 0)
    {
        returned_source_neighbour_index = x.index;
        returned_destination_neighbour_index = y.index;
        approximate_path = geopairs[x_y_in_geo_node_id_for_geo_pair]->path;
        face_sequence_index_list = geopairs[x_y_in_geo_node_id_for_geo_pair]->face_sequence_index_list;
        return geopairs[x_y_in_geo_node_id_for_geo_pair]->distance;
    }

    p = &x;
    while (p->parent != NULL)
    {
        p = p->parent;

        x_in_geo_node_id_for_geo_pair = (*p).id;
        y_in_geo_node_id_for_geo_pair = y.id;
        if (x_in_geo_node_id_for_geo_pair > y_in_geo_node_id_for_geo_pair)
        {
            int temp2 = y_in_geo_node_id_for_geo_pair;
            y_in_geo_node_id_for_geo_pair = x_in_geo_node_id_for_geo_pair;
            x_in_geo_node_id_for_geo_pair = temp2;
        }
        hash_function_two_keys_to_one_key(geo_tree_node_id, x_in_geo_node_id_for_geo_pair, y_in_geo_node_id_for_geo_pair, x_y_in_geo_node_id_for_geo_pair);

        if (geopairs.count(x_y_in_geo_node_id_for_geo_pair) != 0)
        {
            returned_source_neighbour_index = (*p).index;
            returned_destination_neighbour_index = y.index;
            approximate_path = geopairs[x_y_in_geo_node_id_for_geo_pair]->path;
            face_sequence_index_list = geopairs[x_y_in_geo_node_id_for_geo_pair]->face_sequence_index_list;
            return geopairs[x_y_in_geo_node_id_for_geo_pair]->distance;
        }
    }
    p = &y;
    while (p->parent != NULL)
    {
        p = p->parent;

        x_in_geo_node_id_for_geo_pair = x.id;
        y_in_geo_node_id_for_geo_pair = (*p).id;
        if (x_in_geo_node_id_for_geo_pair > y_in_geo_node_id_for_geo_pair)
        {
            int temp3 = y_in_geo_node_id_for_geo_pair;
            y_in_geo_node_id_for_geo_pair = x_in_geo_node_id_for_geo_pair;
            x_in_geo_node_id_for_geo_pair = temp3;
        }
        hash_function_two_keys_to_one_key(geo_tree_node_id, x_in_geo_node_id_for_geo_pair, y_in_geo_node_id_for_geo_pair, x_y_in_geo_node_id_for_geo_pair);

        if (geopairs.count(x_y_in_geo_node_id_for_geo_pair) != 0)
        {
            returned_source_neighbour_index = x.index;
            returned_destination_neighbour_index = (*p).index;
            approximate_path = geopairs[x_y_in_geo_node_id_for_geo_pair]->path;
            face_sequence_index_list = geopairs[x_y_in_geo_node_id_for_geo_pair]->face_sequence_index_list;
            return geopairs[x_y_in_geo_node_id_for_geo_pair]->distance;
        }
    }
    return one_path_query_geo(geo_tree_node_id, *x.parent, *y.parent, returned_source_neighbour_index, returned_destination_neighbour_index, geopairs, poi_unordered_map, approximate_path, face_sequence_index_list);
}

void three_paths_query_geo(int geo_tree_node_id, int source_poi_index, int destination_poi_index, bool &one_path,
                           std::vector<GeoNode *> &all_poi,
                           std::unordered_map<int, GeoNode *> &geo_node_in_partition_tree_unordered_map,
                           std::unordered_map<int, GeoPair *> &geopairs,
                           std::unordered_map<int, int> &poi_unordered_map,
                           double &approximate_distance,
                           std::vector<geodesic::SurfacePoint> &approximate_path,
                           std::vector<int> &face_sequence_index_list)
{
    double approximate_distance_source_short;
    double approximate_distance_destination_short;
    std::vector<geodesic::SurfacePoint> approximate_path_source_short;
    std::vector<geodesic::SurfacePoint> approximate_path_destination_short;
    std::vector<int> face_sequence_index_list_source_short;
    std::vector<int> face_sequence_index_list_destination_short;

    int returned_source_neighbour_index;
    int returned_destination_neighbour_index;

    approximate_distance = one_path_query_geo(geo_tree_node_id, *geo_node_in_partition_tree_unordered_map[all_poi[source_poi_index]->index], *geo_node_in_partition_tree_unordered_map[all_poi[destination_poi_index]->index],
                                              returned_source_neighbour_index, returned_destination_neighbour_index, geopairs, poi_unordered_map, approximate_path, face_sequence_index_list);
    one_path = true;

    if (all_poi[destination_poi_index]->index != returned_destination_neighbour_index)
    {
        int destination_returned_source_neighbour_index;
        int destination_returned_destination_neighbour_index;
        approximate_distance_destination_short = one_path_query_geo(geo_tree_node_id, *geo_node_in_partition_tree_unordered_map[returned_destination_neighbour_index], *geo_node_in_partition_tree_unordered_map[all_poi[destination_poi_index]->index],
                                                                    destination_returned_source_neighbour_index, destination_returned_destination_neighbour_index, geopairs, poi_unordered_map, approximate_path_destination_short, face_sequence_index_list_destination_short);
        assert(destination_returned_source_neighbour_index == returned_destination_neighbour_index && destination_returned_destination_neighbour_index == all_poi[destination_poi_index]->index);
        approximate_distance += approximate_distance_destination_short;
        one_path = false;
    }
    if (all_poi[source_poi_index]->index != returned_source_neighbour_index)
    {
        int source_returned_source_neighbour_index;
        int source_returned_destination_neighbour_index;
        approximate_distance_source_short = one_path_query_geo(geo_tree_node_id, *geo_node_in_partition_tree_unordered_map[all_poi[source_poi_index]->index], *geo_node_in_partition_tree_unordered_map[returned_source_neighbour_index],
                                                               source_returned_source_neighbour_index, source_returned_destination_neighbour_index, geopairs, poi_unordered_map, approximate_path_source_short, face_sequence_index_list_source_short);
        assert(source_returned_source_neighbour_index == all_poi[source_poi_index]->index && source_returned_destination_neighbour_index == returned_source_neighbour_index);
        approximate_distance += approximate_distance_source_short;
        one_path = false;
    }

    if (approximate_path_destination_short.size() > 0)
    {
        if (approximate_path_destination_short[0].getx() == approximate_path[0].getx() &&
            approximate_path_destination_short[0].gety() == approximate_path[0].gety() &&
            approximate_path_destination_short[0].getz() == approximate_path[0].getz())
        {
            std::reverse(approximate_path.begin(), approximate_path.end());
            std::reverse(face_sequence_index_list.begin(), face_sequence_index_list.end());
            for (int i = 0; i < approximate_path_destination_short.size(); i++)
            {
                approximate_path.push_back(approximate_path_destination_short[i]);
            }
            for (int i = 0; i < face_sequence_index_list_destination_short.size(); i++)
            {
                face_sequence_index_list.push_back(face_sequence_index_list_destination_short[i]);
            }
        }
        else if (approximate_path_destination_short[0].getx() == approximate_path[approximate_path.size() - 1].getx() &&
                 approximate_path_destination_short[0].gety() == approximate_path[approximate_path.size() - 1].gety() &&
                 approximate_path_destination_short[0].getz() == approximate_path[approximate_path.size() - 1].getz())
        {
            for (int i = 0; i < approximate_path_destination_short.size(); i++)
            {
                approximate_path.push_back(approximate_path_destination_short[i]);
            }
            for (int i = 0; i < face_sequence_index_list_destination_short.size(); i++)
            {
                face_sequence_index_list.push_back(face_sequence_index_list_destination_short[i]);
            }
        }
        else if (approximate_path_destination_short[approximate_path_destination_short.size() - 1].getx() == approximate_path[0].getx() &&
                 approximate_path_destination_short[approximate_path_destination_short.size() - 1].gety() == approximate_path[0].gety() &&
                 approximate_path_destination_short[approximate_path_destination_short.size() - 1].getz() == approximate_path[0].getz())
        {
            std::reverse(approximate_path.begin(), approximate_path.end());
            std::reverse(face_sequence_index_list.begin(), face_sequence_index_list.end());
            for (int i = approximate_path_destination_short.size() - 1; i >= 0; i--)
            {
                approximate_path.push_back(approximate_path_destination_short[i]);
            }
            for (int i = face_sequence_index_list_destination_short.size() - 1; i >= 0; i--)
            {
                face_sequence_index_list.push_back(face_sequence_index_list_destination_short[i]);
            }
        }
        else if (approximate_path_destination_short[approximate_path_destination_short.size() - 1].getx() == approximate_path[approximate_path.size() - 1].getx() &&
                 approximate_path_destination_short[approximate_path_destination_short.size() - 1].gety() == approximate_path[approximate_path.size() - 1].gety() &&
                 approximate_path_destination_short[approximate_path_destination_short.size() - 1].getz() == approximate_path[approximate_path.size() - 1].getz())
        {
            for (int i = approximate_path_destination_short.size() - 1; i >= 0; i--)
            {
                approximate_path.push_back(approximate_path_destination_short[i]);
            }
            for (int i = face_sequence_index_list_destination_short.size() - 1; i >= 0; i--)
            {
                face_sequence_index_list.push_back(face_sequence_index_list_destination_short[i]);
            }
        }
        else
        {
            assert(false);
        }
    }

    if (approximate_path_source_short.size() > 0)
    {
        if (approximate_path_source_short[0].getx() == approximate_path[0].getx() &&
            approximate_path_source_short[0].gety() == approximate_path[0].gety() &&
            approximate_path_source_short[0].getz() == approximate_path[0].getz())
        {
            std::reverse(approximate_path.begin(), approximate_path.end());
            std::reverse(face_sequence_index_list.begin(), face_sequence_index_list.end());
            for (int i = 0; i < approximate_path_source_short.size(); i++)
            {
                approximate_path.push_back(approximate_path_source_short[i]);
            }
            for (int i = 0; i < face_sequence_index_list_source_short.size(); i++)
            {
                face_sequence_index_list.push_back(face_sequence_index_list_source_short[i]);
            }
        }
        else if (approximate_path_source_short[0].getx() == approximate_path[approximate_path.size() - 1].getx() &&
                 approximate_path_source_short[0].gety() == approximate_path[approximate_path.size() - 1].gety() &&
                 approximate_path_source_short[0].getz() == approximate_path[approximate_path.size() - 1].getz())
        {
            for (int i = 0; i < approximate_path_source_short.size(); i++)
            {
                approximate_path.push_back(approximate_path_source_short[i]);
            }
            for (int i = 0; i < face_sequence_index_list_source_short.size(); i++)
            {
                face_sequence_index_list.push_back(face_sequence_index_list_source_short[i]);
            }
        }
        else if (approximate_path_source_short[approximate_path_source_short.size() - 1].getx() == approximate_path[0].getx() &&
                 approximate_path_source_short[approximate_path_source_short.size() - 1].gety() == approximate_path[0].gety() &&
                 approximate_path_source_short[approximate_path_source_short.size() - 1].getz() == approximate_path[0].getz())
        {
            std::reverse(approximate_path.begin(), approximate_path.end());
            std::reverse(face_sequence_index_list.begin(), face_sequence_index_list.end());
            for (int i = approximate_path_source_short.size() - 1; i >= 0; i--)
            {
                approximate_path.push_back(approximate_path_source_short[i]);
            }
            for (int i = face_sequence_index_list_source_short.size() - 1; i >= 0; i--)
            {
                face_sequence_index_list.push_back(face_sequence_index_list_source_short[i]);
            }
        }
        else if (approximate_path_source_short[approximate_path_source_short.size() - 1].getx() == approximate_path[approximate_path.size() - 1].getx() &&
                 approximate_path_source_short[approximate_path_source_short.size() - 1].gety() == approximate_path[approximate_path.size() - 1].gety() &&
                 approximate_path_source_short[approximate_path_source_short.size() - 1].getz() == approximate_path[approximate_path.size() - 1].getz())
        {
            for (int i = approximate_path_source_short.size() - 1; i >= 0; i--)
            {
                approximate_path.push_back(approximate_path_source_short[i]);
            }
            for (int i = face_sequence_index_list_source_short.size() - 1; i >= 0; i--)
            {
                face_sequence_index_list.push_back(face_sequence_index_list_source_short[i]);
            }
        }
        else
        {
            assert(false);
        }
    }
}

double euclidean_distance(double x_1, double y_1, double x_2, double y_2)
{
    return sqrt(pow(x_1 - x_2, 2) + pow(y_1 - y_2, 2));
}

void min_but_greater_than(std::vector<double> &list, double greater_than_value, double &return_min, int &return_index)
{
    double min = 1e100;
    int index = -1;
    for (int i = 0; i < list.size(); i++)
    {
        if (list[i] > greater_than_value)
        {
            if (list[i] < min)
            {
                min = list[i];
                index = i;
            }
        }
    }
    return_min = min;
    return_index = index;
}
void sort_min_to_max_and_get_original_index(std::vector<double> &vec, std::vector<std::pair<double, int>> &vp)
{
    for (int i = 0; i < vec.size(); i++)
    {
        vp.push_back(std::make_pair(vec[i], i));
    }
    std::sort(vp.begin(), vp.end());
}

// for Dijkstra shortest path, https://www.geeksforgeeks.org/dijkstras-shortest-path-algorithm-greedy-algo-7/
// for Kruskal MST, https://www.geeksforgeeks.org/kruskals-minimum-spanning-tree-using-stl-in-c/
#define INF 0x3f3f3f3f

// for Kruskal MST
struct DisjointSets
{
    int *parent, *rnk;
    int n;

    DisjointSets(int n)
    {
        this->n = n;
        parent = new int[n + 1];
        rnk = new int[n + 1];
        for (int i = 0; i <= n; i++)
        {
            rnk[i] = 0;
            parent[i] = i;
        }
    }

    int find(int u)
    {
        if (u != parent[u])
        {
            parent[u] = find(parent[u]);
        }
        return parent[u];
    }

    void merge(int x, int y)
    {
        x = find(x), y = find(y);
        if (rnk[x] > rnk[y])
        {
            parent[y] = x;
        }
        else
        {
            parent[x] = y;
        }
        if (rnk[x] == rnk[y])
        {
            rnk[y]++;
        }
    }
};

class Graph
{
    int V;
    int E;
    double total_weight;                                       // for Dijkstra shortest path
    int total_geo_path_size;                                   // for Dijkstra shortest path
    std::list<std::pair<int, double>> *adj;                    // for Dijkstra shortest path
    std::vector<std::pair<double, std::pair<int, int>>> edges; // for Kruskal MST

public:
    // for Dijkstra shortest path
    Graph(int V);
    void add_edge_Dijkstra(int u, int v, double w);
    void add_edge_and_geo_path_size_Dijkstra(int u, int v, double w, int geo_path_size);
    int get_edge_num_Dijkstra();
    void update_edge_Dijkstra(int u, int v, double w);
    void shortest_distance_Dijkstra(int src, std::vector<double> &dist);
    void shortest_distance_Dijkstra(int src, std::vector<double> &dist, double max_dist);
    void shortest_distance_Dijkstra(int src, int dest, std::vector<double> &dist);
    void shortest_path_Dijkstra(int src, std::vector<double> &dist, std::vector<std::vector<int>> &path);
    void shortest_path_Dijkstra(int src, std::vector<double> &dist, std::vector<std::vector<int>> &path, double max_dist);
    void shortest_path_Dijkstra(int src, int dest, std::vector<double> &dist, std::vector<std::vector<int>> &path);
    double get_total_weight_Dijkstra();
    int get_total_geo_path_size_Dijkstra();

    // for Kruskal MST
    Graph(int V, int E);
    void add_edge_Kruskal(int u, int v, double w);
    int MST_Kruskal();
};

Graph::Graph(int V)
{
    this->V = V;
    E = 0;
    total_weight = 0;
    total_geo_path_size = 0;
    adj = new std::list<std::pair<int, double>>[V];
}

void Graph::add_edge_Dijkstra(int u, int v, double w)
{
    adj[u].push_back(std::make_pair(v, w));
    adj[v].push_back(std::make_pair(u, w));
    E++;
    total_weight += w;
}

void Graph::add_edge_and_geo_path_size_Dijkstra(int u, int v, double w, int geo_path_size)
{
    adj[u].push_back(std::make_pair(v, w));
    adj[v].push_back(std::make_pair(u, w));
    E++;
    total_weight += w;
    total_geo_path_size += geo_path_size;
}

int Graph::get_edge_num_Dijkstra()
{
    return E;
}

void Graph::update_edge_Dijkstra(int u, int v, double w)
{
    for (std::list<std::pair<int, double>>::iterator i = adj[u].begin(); i != adj[u].end(); i++)
    {
        if ((*i).first == v)
        {
            (*i).second = w;
            break;
        }
    }

    for (std::list<std::pair<int, double>>::iterator i = adj[v].begin(); i != adj[v].end(); i++)
    {
        if ((*i).first == u)
        {
            (*i).second = w;
            break;
        }
    }
}

void Graph::shortest_distance_Dijkstra(int src, std::vector<double> &dist)
{
    std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, std::greater<std::pair<int, double>>> pq;
    pq.push(std::make_pair(0, src));
    dist[src] = 0;

    while (!pq.empty())
    {
        int u = pq.top().second;
        pq.pop();
        std::list<std::pair<int, double>>::iterator i;
        for (i = adj[u].begin(); i != adj[u].end(); ++i)
        {
            int v = (*i).first;
            double weight = (*i).second;
            if (dist[v] > dist[u] + weight)
            {
                dist[v] = dist[u] + weight;
                pq.push(std::make_pair(dist[v], v));
            }
        }
    }
}

// given the max searching distance, if reach this distance, break the loop
void Graph::shortest_distance_Dijkstra(int src, std::vector<double> &dist, double max_dist)
{
    std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, std::greater<std::pair<int, double>>> pq;
    pq.push(std::make_pair(0, src));
    dist[src] = 0;

    while (!pq.empty())
    {
        int u = pq.top().second;
        if (pq.top().first > max_dist)
        {
            break;
        }
        pq.pop();
        std::list<std::pair<int, double>>::iterator i;
        for (i = adj[u].begin(); i != adj[u].end(); ++i)
        {
            int v = (*i).first;
            double weight = (*i).second;
            if (dist[v] > dist[u] + weight)
            {
                dist[v] = dist[u] + weight;
                pq.push(std::make_pair(dist[v], v));
            }
        }
    }
}

// given the destination index, if reach this destination, break the loop
void Graph::shortest_distance_Dijkstra(int src, int dest, std::vector<double> &dist)
{
    std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, std::greater<std::pair<int, double>>> pq;
    pq.push(std::make_pair(0, src));
    dist[src] = 0;

    while (!pq.empty())
    {
        int u = pq.top().second;
        if (u == dest)
        {
            break;
        }
        pq.pop();
        std::list<std::pair<int, double>>::iterator i;
        for (i = adj[u].begin(); i != adj[u].end(); ++i)
        {
            int v = (*i).first;
            double weight = (*i).second;
            if (dist[v] > dist[u] + weight)
            {
                dist[v] = dist[u] + weight;
                pq.push(std::make_pair(dist[v], v));
            }
        }
    }
}

void Graph::shortest_path_Dijkstra(int src, std::vector<double> &dist, std::vector<std::vector<int>> &path)
{
    std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, std::greater<std::pair<int, double>>> pq;
    pq.push(std::make_pair(0, src));
    dist[src] = 0;
    std::vector<int> prev(dist.size());

    while (!pq.empty())
    {
        int u = pq.top().second;
        pq.pop();
        std::list<std::pair<int, double>>::iterator i;
        for (i = adj[u].begin(); i != adj[u].end(); ++i)
        {
            int v = (*i).first;
            double weight = (*i).second;
            if (dist[v] > dist[u] + weight)
            {
                dist[v] = dist[u] + weight;
                prev[v] = u;
                pq.push(std::make_pair(dist[v], v));
            }
        }
    }
    for (int i = 0; i < dist.size(); i++)
    {
        // note that the path is in reverse order
        path[i].push_back(i);
        if (i != src)
        {
            int j = i;
            do
            {
                j = prev[j];
                path[i].push_back(j);
            } while (j != src);
        }
    }
}

// given the max searching distance, if reach this distance, break the loop
void Graph::shortest_path_Dijkstra(int src, std::vector<double> &dist, std::vector<std::vector<int>> &path, double max_dist)
{
    std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, std::greater<std::pair<int, double>>> pq;
    pq.push(std::make_pair(0, src));
    dist[src] = 0;
    std::vector<int> prev(dist.size());

    while (!pq.empty())
    {
        int u = pq.top().second;
        if (pq.top().first > max_dist)
        {
            break;
        }
        pq.pop();
        std::list<std::pair<int, double>>::iterator i;
        for (i = adj[u].begin(); i != adj[u].end(); ++i)
        {
            int v = (*i).first;
            double weight = (*i).second;
            if (dist[v] > dist[u] + weight)
            {
                dist[v] = dist[u] + weight;
                prev[v] = u;
                pq.push(std::make_pair(dist[v], v));
            }
        }
    }
    for (int i = 0; i < dist.size(); i++)
    {
        // note that the path is in reverse order
        path[i].push_back(i);
        if (i != src)
        {
            int j = i;
            do
            {
                j = prev[j];
                path[i].push_back(j);
            } while (j != src);
        }
    }
}

// given the destination index, if reach this destination, break the loop
void Graph::shortest_path_Dijkstra(int src, int dest, std::vector<double> &dist, std::vector<std::vector<int>> &path)
{
    std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, std::greater<std::pair<int, double>>> pq;
    pq.push(std::make_pair(0, src));
    dist[src] = 0;
    std::vector<int> prev(dist.size());

    while (!pq.empty())
    {
        int u = pq.top().second;
        if (u == dest)
        {
            break;
        }
        pq.pop();
        std::list<std::pair<int, double>>::iterator i;
        for (i = adj[u].begin(); i != adj[u].end(); ++i)
        {
            int v = (*i).first;
            double weight = (*i).second;
            if (dist[v] > dist[u] + weight)
            {
                dist[v] = dist[u] + weight;
                prev[v] = u;
                pq.push(std::make_pair(dist[v], v));
            }
        }
    }

    int i = dest;
    // note that the path is in reverse order
    path[i].push_back(i);
    if (i != src)
    {
        int j = i;
        do
        {
            j = prev[j];
            path[i].push_back(j);
        } while (j != src);
    }
}

double Graph::get_total_weight_Dijkstra()
{
    return total_weight;
}

int Graph::get_total_geo_path_size_Dijkstra()
{
    return total_geo_path_size;
}

Graph::Graph(int V, int E)
{
    this->V = V;
    this->E = E;
}

void Graph::add_edge_Kruskal(int u, int v, double w)
{
    edges.push_back({w, {u, v}});
}

int Graph::MST_Kruskal()
{
    int mst_wt = 0;
    sort(edges.begin(), edges.end());
    DisjointSets ds(V);
    std::vector<std::pair<double, std::pair<int, int>>>::iterator it;
    for (it = edges.begin(); it != edges.end(); it++)
    {
        int u = it->second.first;
        int v = it->second.second;
        int set_u = ds.find(u);
        int set_v = ds.find(v);

        if (set_u != set_v)
        {
            mst_wt += it->first;
            ds.merge(set_u, set_v);
        }
    }
    return mst_wt;
}

// calculate the mean and standard deviation
void cal_mean_std(std::vector<double> &value, double &mean, double &std)
{
    double sum = 0;
    double variance = 0;
    assert(value.size() > 1);
    for (int i = 0; i < value.size(); i++)
    {
        sum += value[i];
    }
    mean = sum / value.size();
    for (int i = 0; i < value.size(); i++)
    {
        variance += pow((value[i] - mean), 2);
    }
    variance = variance / (value.size() - 1);
    std = sqrt(variance);
}

void remove_outliers(std::vector<double> &with_outliers,
                     std::vector<double> &without_outliers)
{
    double mean;
    double std;
    cal_mean_std(with_outliers, mean, std);
    double upper_limit = mean + 3 * std;
    double lower_limit = mean - 3 * std;

    for (int i = 0; i < with_outliers.size(); i++)
    {
        if (with_outliers[i] <= upper_limit && with_outliers[i] >= lower_limit)
        {
            without_outliers.push_back(with_outliers[i]);
        }
    }
}
