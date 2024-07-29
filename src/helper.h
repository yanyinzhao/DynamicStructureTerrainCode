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
                 std::unordered_map<int, double> &pre_distance_poi_to_poi_map,
                 std::unordered_map<int, std::vector<geodesic::SurfacePoint>> &pre_path_poi_to_poi_map,
                 std::unordered_map<int, int> &poi_unordered_map)
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

        double const distance_limit = n->radius;

        auto bite = current_parent_covers_but_remained_pois.begin();

        while (bite != current_parent_covers_but_remained_pois.end())
        {
            int x_in_poi_list = poi_unordered_map[n->index];
            int y_in_poi_list = poi_unordered_map[(*bite).first];
            int x_y_in_poi_list;
            for (int i = 0; i < poi_unordered_map.size() * 100; i++)
            {
                if (x_in_poi_list <= y_in_poi_list)
                {
                    hash_function_two_keys_to_one_key(poi_unordered_map.size(), x_in_poi_list, y_in_poi_list, x_y_in_poi_list);
                }
                else
                {
                    hash_function_two_keys_to_one_key(poi_unordered_map.size(), y_in_poi_list, x_in_poi_list, x_y_in_poi_list);
                }
            }
            if (pre_distance_poi_to_poi_map[x_y_in_poi_list] <= distance_limit)
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

            double const distance_limit = b.second->radius;

            auto bite = current_parent_covers_but_remained_pois.begin();

            while (bite != current_parent_covers_but_remained_pois.end())
            {
                int x_in_poi_list = poi_unordered_map[a->index];
                int y_in_poi_list = poi_unordered_map[(*bite).first];
                int x_y_in_poi_list;
                for (int i = 0; i < poi_unordered_map.size() * 100; i++)
                {
                    if (x_in_poi_list <= y_in_poi_list)
                    {
                        hash_function_two_keys_to_one_key(poi_unordered_map.size(), x_in_poi_list, y_in_poi_list, x_y_in_poi_list);
                    }
                    else
                    {
                        hash_function_two_keys_to_one_key(poi_unordered_map.size(), y_in_poi_list, x_in_poi_list, x_y_in_poi_list);
                    }
                }
                if (pre_distance_poi_to_poi_map[x_y_in_poi_list] <= distance_limit)
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
                    std::unordered_map<int, double> &pre_distance_poi_to_poi_map,
                    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> &pre_path_poi_to_poi_map,
                    std::unordered_map<int, int> &poi_unordered_map)
{
    int depth = 0;
    pois_as_center_each_parent_layer.clear();
    pois_as_center_each_parent_layer.insert(root_geo.index, &root_geo);

    stx::btree<int, GeoNode *> remained_pois = stx::btree<int, GeoNode *>(pois_B_tree);

    remained_pois.erase(root_geo.index);

    double const distance_limit = root_geo.radius;
    for (stx::btree<int, GeoNode *>::iterator bite = remained_pois.begin(); bite != remained_pois.end(); bite++)
    {
        int x_in_poi_list = poi_unordered_map[root_geo.index];
        int y_in_poi_list = poi_unordered_map[(*bite).first];
        int x_y_in_poi_list;
        for (int i = 0; i < poi_unordered_map.size() * 100; i++)
        {
            if (x_in_poi_list <= y_in_poi_list)
            {
                hash_function_two_keys_to_one_key(poi_num, x_in_poi_list, y_in_poi_list, x_y_in_poi_list);
            }
            else
            {
                hash_function_two_keys_to_one_key(poi_num, y_in_poi_list, x_in_poi_list, x_y_in_poi_list);
            }
        }
        if (pre_distance_poi_to_poi_map[x_y_in_poi_list] <= distance_limit)
        {
            (*bite).second->set_covered_by(&root_geo);
            (*bite).second->set_covers(&root_geo);
        }
    }

    while (pois_as_center_each_parent_layer.size() != poi_num)
    {
        build_level(geo_tree_node_id, mesh, depth, pois_B_tree, pois_as_center_each_parent_layer, pre_distance_poi_to_poi_map, pre_path_poi_to_poi_map, poi_unordered_map);
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

void pre_compute_WSPD_RC_TIN_Oracle(int poi_num, geodesic::Mesh *mesh, std::vector<int> &poi_list,
                                    std::unordered_map<int, double> &distance_poi_to_poi_map,
                                    std::unordered_map<int, std::vector<geodesic::SurfacePoint>> &path_poi_to_poi_map,
                                    double &memory_usage, int notA2A_one_A2A_two)
{
    geodesic::GeodesicAlgorithmExact algorithm(mesh);
    double const distance_limit = geodesic::GEODESIC_INF;
    distance_poi_to_poi_map.clear();
    path_poi_to_poi_map.clear();
    int path_poi_to_poi_size = 0;
    std::vector<geodesic::SurfacePoint> one_source_poi_list;
    std::vector<geodesic::SurfacePoint> destinations_poi_list;

    for (int i = 0; i < poi_num; i++)
    {
        one_source_poi_list.clear();
        destinations_poi_list.clear();
        one_source_poi_list.push_back(geodesic::SurfacePoint(&mesh->vertices()[poi_list[i]]));
        for (int j = i; j < poi_num; j++)
        {
            destinations_poi_list.push_back(geodesic::SurfacePoint(&mesh->vertices()[poi_list[j]]));
        }
        if (notA2A_one_A2A_two == 1)
        {
            algorithm.propagate(one_source_poi_list, distance_limit);
        }
        else if (notA2A_one_A2A_two == 2)
        {
            for (int j = 0; j < mesh->vertices().size() / poi_num; j++)
            {
                algorithm.propagate(one_source_poi_list, distance_limit);
            }
        }
        for (int j = i; j < poi_num; j++)
        {
            std::vector<geodesic::SurfacePoint> path;
            geodesic::SurfacePoint one_destination(&mesh->vertices()[poi_list[j]]);
            algorithm.trace_back(one_destination, path);
            int i_j;
            hash_function_two_keys_to_one_key(poi_num, i, j, i_j);
            distance_poi_to_poi_map[i_j] = length(path);
            path_poi_to_poi_map[i_j] = path;
            path_poi_to_poi_size += path.size();
        }
    }
    memory_usage += algorithm.get_memory() + 0.5 * poi_num * (poi_num - 1) * sizeof(double) + path_poi_to_poi_size * sizeof(geodesic::SurfacePoint);
}

void generate_geo_pair(int geo_tree_node_id, int &WSPD_RC_TIN_oracle_edge_num, double &WSPD_RC_TIN_oracle_weight,
                       geodesic::Mesh *mesh, GeoNode &x, GeoNode &y,
                       double epsilon,
                       std::unordered_map<int, GeoPair *> &geopairs,
                       std::unordered_map<int, int> &poi_unordered_map,
                       std::unordered_map<int, int> &geo_pair_unordered_map,
                       std::unordered_map<int, double> &pre_distance_poi_to_poi_map,
                       std::unordered_map<int, std::vector<geodesic::SurfacePoint>> &pre_path_poi_to_poi_map,
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
        int x_in_poi_list = poi_unordered_map[x.index];
        int y_in_poi_list = poi_unordered_map[y.index];
        int x_y_in_poi_list;
        for (int i = 0; i < poi_unordered_map.size() * 100; i++)
        {
            if (x_in_poi_list <= y_in_poi_list)
            {
                hash_function_two_keys_to_one_key(poi_unordered_map.size(), x_in_poi_list, y_in_poi_list, x_y_in_poi_list);
            }
            else
            {
                hash_function_two_keys_to_one_key(poi_unordered_map.size(), y_in_poi_list, x_in_poi_list, x_y_in_poi_list);
            }
        }

        double distancexy = pre_distance_poi_to_poi_map[x_y_in_poi_list];
        std::vector<geodesic::SurfacePoint> pathxy = pre_path_poi_to_poi_map[x_y_in_poi_list];

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
            WSPD_RC_TIN_oracle_edge_num++;
            GeoPair *nodepair = new GeoPair();
            nodepair->node1 = &x;
            nodepair->node2 = &y;
            nodepair->distance = distancexy;
            nodepair->path = pathxy;
            nodepair->face_sequence_index_list = face_sequence_index_listxy;
            WSPD_RC_TIN_oracle_weight += distancexy;
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
                    generate_geo_pair(geo_tree_node_id, WSPD_RC_TIN_oracle_edge_num, WSPD_RC_TIN_oracle_weight, mesh, (**ite), y, epsilon, geopairs, poi_unordered_map, geo_pair_unordered_map, pre_distance_poi_to_poi_map, pre_path_poi_to_poi_map, pairwise_path_poi_to_poi_size, face_sequence_index_list_size);
                }
            }
            else
            {
                for (std::list<GeoNode *>::iterator jte = y.children.begin(); jte != y.children.end(); jte++)
                {
                    generate_geo_pair(geo_tree_node_id, WSPD_RC_TIN_oracle_edge_num, WSPD_RC_TIN_oracle_weight, mesh, x, (**jte), epsilon, geopairs, poi_unordered_map, geo_pair_unordered_map, pre_distance_poi_to_poi_map, pre_path_poi_to_poi_map, pairwise_path_poi_to_poi_size, face_sequence_index_list_size);
                }
            }
        }
    }
}

double one_path_query_geo(int geo_tree_node_id, GeoNode &x, GeoNode &y,
                          int &returned_source_neighbour_index, int &returned_destination_neighbour_index,
                          std::unordered_map<int, GeoPair *> &geopairs,
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
    return one_path_query_geo(geo_tree_node_id, *x.parent, *y.parent, returned_source_neighbour_index, returned_destination_neighbour_index, geopairs, approximate_path, face_sequence_index_list);
}

void three_paths_query_geo(geodesic::Mesh *mesh, int geo_tree_node_id, int WSPD_RC_TINone_EARtwo,
                           int WSPD_RC_TIN_source_poi_index_EAR_source_vertex_id,
                           int WSPD_RC_TIN_destination_poi_index_EAR_destination_vertex_id, bool &one_path,
                           std::vector<GeoNode *> &all_poi,
                           std::unordered_map<int, GeoNode *> &geo_node_in_partition_tree_unordered_map,
                           std::unordered_map<int, GeoPair *> &geopairs,
                           double &approximate_distance,
                           std::vector<geodesic::SurfacePoint> &approximate_path,
                           std::vector<int> &face_sequence_index_list)
{
    geodesic::GeodesicAlgorithmExact algorithm(mesh);

    double approximate_distance_source_short;
    double approximate_distance_destination_short;
    std::vector<geodesic::SurfacePoint> approximate_path_source_short;
    std::vector<geodesic::SurfacePoint> approximate_path_destination_short;
    std::vector<int> face_sequence_index_list_source_short;
    std::vector<int> face_sequence_index_list_destination_short;

    int returned_source_neighbour_index;
    int returned_destination_neighbour_index;

    int src_vertex_id;
    int dest_vertex_id;

    if (WSPD_RC_TINone_EARtwo == 1)
    {
        src_vertex_id = all_poi[WSPD_RC_TIN_source_poi_index_EAR_source_vertex_id]->index;
        dest_vertex_id = all_poi[WSPD_RC_TIN_destination_poi_index_EAR_destination_vertex_id]->index;
    }
    else if (WSPD_RC_TINone_EARtwo == 2)
    {
        src_vertex_id = WSPD_RC_TIN_source_poi_index_EAR_source_vertex_id;
        dest_vertex_id = WSPD_RC_TIN_destination_poi_index_EAR_destination_vertex_id;
    }

    approximate_distance = one_path_query_geo(geo_tree_node_id, *geo_node_in_partition_tree_unordered_map[src_vertex_id], *geo_node_in_partition_tree_unordered_map[dest_vertex_id],
                                              returned_source_neighbour_index, returned_destination_neighbour_index, geopairs, approximate_path, face_sequence_index_list);
    one_path = true;

    if (dest_vertex_id != returned_destination_neighbour_index)
    {
        int destination_returned_source_neighbour_index;
        int destination_returned_destination_neighbour_index;
        approximate_distance_destination_short = one_path_query_geo(geo_tree_node_id, *geo_node_in_partition_tree_unordered_map[returned_destination_neighbour_index], *geo_node_in_partition_tree_unordered_map[dest_vertex_id],
                                                                    destination_returned_source_neighbour_index, destination_returned_destination_neighbour_index, geopairs, approximate_path_destination_short, face_sequence_index_list_destination_short);
        if (destination_returned_source_neighbour_index != returned_destination_neighbour_index || destination_returned_destination_neighbour_index != dest_vertex_id)
        {
            std::vector<geodesic::SurfacePoint> one_source_list;
            one_source_list.clear();
            one_source_list.push_back(geodesic::SurfacePoint(&mesh->vertices()[dest_vertex_id]));
            std::vector<geodesic::SurfacePoint> one_destination_list;
            one_destination_list.clear();
            one_destination_list.push_back(geodesic::SurfacePoint(&mesh->vertices()[returned_destination_neighbour_index]));
            algorithm.propagate_cover_dest(one_source_list, &one_destination_list);
            geodesic::SurfacePoint one_destination(&mesh->vertices()[returned_destination_neighbour_index]);
            algorithm.trace_back(one_destination, approximate_path_destination_short);
            approximate_distance_destination_short = length(approximate_path_destination_short);
            face_sequence_index_list_destination_short.clear();
            get_face_sequence(approximate_path_destination_short, face_sequence_index_list_destination_short);
        }
        approximate_distance += approximate_distance_destination_short;
        one_path = false;
    }
    if (src_vertex_id != returned_source_neighbour_index)
    {
        int source_returned_source_neighbour_index;
        int source_returned_destination_neighbour_index;
        approximate_distance_source_short = one_path_query_geo(geo_tree_node_id, *geo_node_in_partition_tree_unordered_map[src_vertex_id], *geo_node_in_partition_tree_unordered_map[returned_source_neighbour_index],
                                                               source_returned_source_neighbour_index, source_returned_destination_neighbour_index, geopairs, approximate_path_source_short, face_sequence_index_list_source_short);
        if (source_returned_source_neighbour_index != src_vertex_id || source_returned_destination_neighbour_index != returned_source_neighbour_index)
        {
            std::vector<geodesic::SurfacePoint> one_source_list;
            one_source_list.clear();
            one_source_list.push_back(geodesic::SurfacePoint(&mesh->vertices()[src_vertex_id]));
            std::vector<geodesic::SurfacePoint> one_destination_list;
            one_destination_list.clear();
            one_destination_list.push_back(geodesic::SurfacePoint(&mesh->vertices()[returned_source_neighbour_index]));
            algorithm.propagate_cover_dest(one_source_list, &one_destination_list);
            geodesic::SurfacePoint one_destination(&mesh->vertices()[returned_source_neighbour_index]);
            algorithm.trace_back(one_destination, approximate_path_source_short);
            approximate_distance_source_short = length(approximate_path_source_short);
            face_sequence_index_list_source_short.clear();
            get_face_sequence(approximate_path_source_short, face_sequence_index_list_source_short);
        }
        // assert(source_returned_source_neighbour_index == src_vertex_id && source_returned_destination_neighbour_index == returned_source_neighbour_index);
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

int doubleCmp(const double &x)
{
    if (fabs(x) < 1e-5)
        return 0;
    return x < 0 ? -1 : 1;
}

int cross_product(double p1_x, double p1_y, double p2_x, double p2_y)
{
    return doubleCmp(p1_x * p2_y - p2_x * p1_y);
}

bool segment_intersection(double a_x, double a_y, double b_x, double b_y,
                          double c_x, double c_y, double d_x, double d_y)
{
    if (doubleCmp(std::max(a_x, b_x) - std::min(c_x, d_x)) < 0 || doubleCmp(std::max(a_y, b_y) - std::min(c_y, d_y)) < 0 ||
        doubleCmp(std::max(c_x, d_x) - std::min(a_x, b_x)) < 0 || doubleCmp(std::max(c_y, d_y) - std::min(a_y, b_y)) < 0)
        return false;
    float dir1 = (c_x - a_x) * (b_y - a_y) - (b_x - a_x) * (c_y - a_y);
    float dir2 = (d_x - a_x) * (b_y - a_y) - (b_x - a_x) * (d_y - a_y);
    float dir3 = (a_x - c_x) * (d_y - c_y) - (d_x - c_x) * (a_y - c_y);
    float dir4 = (b_x - c_x) * (d_y - c_y) - (d_x - c_x) * (b_y - c_y);
    return doubleCmp(dir1 * dir2) <= 0 && doubleCmp(dir3 * dir4) <= 0;
}

//  triangle:ABC, segment:PQ
bool triangle_segment_intersection(double a_x, double a_y, double b_x, double b_y,
                                   double c_x, double c_y, double p_x, double p_y,
                                   double q_x, double q_y)
{
    int dir1 = cross_product(b_x - a_x, b_y - a_y, p_x - a_x, p_y - a_y);
    int dir2 = cross_product(c_x - b_x, c_y - b_y, p_x - b_x, p_y - b_y);
    int dir3 = cross_product(a_x - c_x, a_y - c_y, p_x - c_x, p_y - c_y);
    bool flag_p = false, flag_q = false;
    if (dir1 == dir2 && dir2 == dir3)
    {
        flag_p = true;
    }
    dir1 = cross_product(b_x - a_x, b_y - a_y, q_x - a_x, q_y - a_y);
    dir2 = cross_product(c_x - b_x, c_y - b_y, q_x - b_x, q_y - b_y);
    dir3 = cross_product(a_x - c_x, a_y - c_y, q_x - c_x, q_y - c_y);
    if (dir1 == dir2 && dir2 == dir3)
    {
        flag_q = true;
    }
    if (flag_p && flag_q)
    {
        return true;
    }
    return segment_intersection(a_x, a_y, b_x, b_y, p_x, p_y, q_x, q_y) ||
           segment_intersection(b_x, b_y, c_x, c_y, p_x, p_y, q_x, q_y) ||
           segment_intersection(c_x, c_y, a_x, a_y, p_x, p_y, q_x, q_y);
}

void divide_mesh_into_box(geodesic::Mesh *mesh, int sqrt_num_of_box,
                          std::unordered_map<int, int> &highway_node_id_map,
                          std::unordered_map<int, std::unordered_map<int, int>> &highway_node_id_with_box_id_map)
{
    double box_width = mesh->m_width / sqrt_num_of_box;
    double box_height = mesh->m_height / sqrt_num_of_box;
    int box_id = 0;

    for (int i = 0; i < sqrt_num_of_box; i++)
    {
        for (int j = 0; j < sqrt_num_of_box; j++)
        {
            double x_min = mesh->m_xmin + i * box_width;
            double x_max = mesh->m_xmin + (i + 1) * box_width;
            double y_min = mesh->m_ymin + j * box_height;
            double y_max = mesh->m_ymin + (j + 1) * box_height;

            // std::cout << x_min << " " << x_max << " " << y_min << " " << y_max << std::endl;

            std::unordered_map<int, int> one_highway_node_id;
            one_highway_node_id.clear();

            for (int k = 0; k < mesh->faces().size(); k++)
            {
                double a_x = mesh->faces()[k].adjacent_vertices()[0]->getx();
                double a_y = mesh->faces()[k].adjacent_vertices()[0]->gety();
                double b_x = mesh->faces()[k].adjacent_vertices()[1]->getx();
                double b_y = mesh->faces()[k].adjacent_vertices()[1]->gety();
                double c_x = mesh->faces()[k].adjacent_vertices()[2]->getx();
                double c_y = mesh->faces()[k].adjacent_vertices()[2]->gety();

                if (triangle_segment_intersection(a_x, a_y, b_x, b_y, c_x, c_y, x_min, y_min, x_max, y_min) ||
                    triangle_segment_intersection(a_x, a_y, b_x, b_y, c_x, c_y, x_min, y_min, x_min, y_max) ||
                    triangle_segment_intersection(a_x, a_y, b_x, b_y, c_x, c_y, x_min, y_max, x_max, y_max) ||
                    triangle_segment_intersection(a_x, a_y, b_x, b_y, c_x, c_y, x_max, y_min, x_max, y_max))
                {
                    int v1_id = mesh->faces()[k].adjacent_vertices()[0]->id();
                    int v2_id = mesh->faces()[k].adjacent_vertices()[1]->id();
                    int v3_id = mesh->faces()[k].adjacent_vertices()[2]->id();
                    one_highway_node_id[v1_id] = v1_id;
                    one_highway_node_id[v2_id] = v2_id;
                    one_highway_node_id[v3_id] = v3_id;
                    highway_node_id_map[v1_id] = v1_id;
                    highway_node_id_map[v2_id] = v2_id;
                    highway_node_id_map[v3_id] = v3_id;
                    // std::cout << v1_id << std::endl;
                    // std::cout << v2_id << std::endl;
                    // std::cout << v3_id << std::endl;
                }
            }
            highway_node_id_with_box_id_map[box_id] = one_highway_node_id;
            box_id++;
        }
    }
}

void pre_compute_EAR_Oracle_highway_node(geodesic::Mesh *mesh, std::vector<int> &highway_node_list,
                                         std::unordered_map<int, double> &distance_highway_node_to_highway_node_map,
                                         std::unordered_map<int, std::vector<geodesic::SurfacePoint>> &path_highway_node_to_highway_node_map,
                                         double &memory_usage)
{
    geodesic::GeodesicAlgorithmExact algorithm(mesh);
    double const distance_limit = geodesic::GEODESIC_INF;
    distance_highway_node_to_highway_node_map.clear();
    path_highway_node_to_highway_node_map.clear();
    int path_highway_node_to_highway_node_size = 0;
    std::vector<geodesic::SurfacePoint> one_source_highway_node_list;
    std::vector<geodesic::SurfacePoint> destinations_highway_node_list;

    for (int i = 0; i < highway_node_list.size(); i++)
    {
        one_source_highway_node_list.clear();
        destinations_highway_node_list.clear();
        one_source_highway_node_list.push_back(geodesic::SurfacePoint(&mesh->vertices()[highway_node_list[i]]));
        for (int j = i; j < highway_node_list.size(); j++)
        {
            destinations_highway_node_list.push_back(geodesic::SurfacePoint(&mesh->vertices()[highway_node_list[j]]));
        }
        algorithm.propagate(one_source_highway_node_list, distance_limit);
        for (int j = i; j < highway_node_list.size(); j++)
        {
            std::vector<geodesic::SurfacePoint> path;
            geodesic::SurfacePoint one_destination(&mesh->vertices()[highway_node_list[j]]);
            algorithm.trace_back(one_destination, path);
            int i_j;
            hash_function_two_keys_to_one_key(highway_node_list.size(), i, j, i_j);
            distance_highway_node_to_highway_node_map[i_j] = length(path);
            path_highway_node_to_highway_node_map[i_j] = path;
            path_highway_node_to_highway_node_size += path.size();
        }
    }
    memory_usage += algorithm.get_memory() + 0.5 * highway_node_list.size() * (highway_node_list.size() - 1) * sizeof(double) + path_highway_node_to_highway_node_size * sizeof(geodesic::SurfacePoint);
}

void poi_to_highway_node_path(geodesic::Mesh *mesh, int sqrt_num_of_box,
                              std::unordered_map<int, std::unordered_map<int, int>> &highway_node_id_with_box_id_map,
                              std::vector<int> &poi_list, std::unordered_map<int, double> &distance_poi_to_highway_node_map,
                              std::unordered_map<int, std::vector<geodesic::SurfacePoint>> &path_poi_to_highway_node_map,
                              double &EAR_Oracle_size, double &memory_usage)
{
    geodesic::GeodesicAlgorithmExact algorithm(mesh);
    std::vector<geodesic::SurfacePoint> one_source_poi_list;
    std::vector<geodesic::SurfacePoint> destinations_highway_node_list;

    for (int k = 0; k < poi_list.size(); k++)
    {
        one_source_poi_list.clear();
        one_source_poi_list.push_back(geodesic::SurfacePoint(&mesh->vertices()[poi_list[k]]));

        double box_width = mesh->m_width / sqrt_num_of_box;
        double box_height = mesh->m_height / sqrt_num_of_box;
        int box_id = 0;
        double x_min;
        double x_max;
        double y_min;
        double y_max;

        double poi_x = mesh->vertices()[poi_list[k]].getx();
        double poi_y = mesh->vertices()[poi_list[k]].gety();

        for (int i = 0; i < sqrt_num_of_box; i++)
        {
            for (int j = 0; j < sqrt_num_of_box; j++)
            {
                x_min = mesh->m_xmin + i * box_width;
                x_max = mesh->m_xmin + (i + 1) * box_width;
                y_min = mesh->m_ymin + j * box_height;
                y_max = mesh->m_ymin + (j + 1) * box_height;

                if (poi_x >= x_min && poi_x <= x_max & poi_y >= y_min && poi_y <= y_max)
                {
                    i = sqrt_num_of_box;
                    j = sqrt_num_of_box;
                    continue;
                }
                box_id++;
            }
        }

        std::vector<int> same_box_poi_vertex_id;
        same_box_poi_vertex_id.clear();
        for (int i = 0; i < poi_list.size(); i++)
        {
            if (i == k)
            {
                continue;
            }
            double poi2_x = mesh->vertices()[poi_list[i]].getx();
            double poi2_y = mesh->vertices()[poi_list[i]].gety();
            if (poi2_x >= x_min && poi2_x <= x_max & poi2_y >= y_min && poi2_y <= y_max)
            {
                same_box_poi_vertex_id.push_back(poi_list[i]);
            }
        }

        // add highway node in of current box
        for (auto i : highway_node_id_with_box_id_map[box_id])
        {
            destinations_highway_node_list.push_back(geodesic::SurfacePoint(&mesh->vertices()[i.first]));
        }
        // add poi in the current box
        for (int i = 0; i < same_box_poi_vertex_id.size(); i++)
        {
            destinations_highway_node_list.push_back(geodesic::SurfacePoint(&mesh->vertices()[same_box_poi_vertex_id[i]]));
        }

        algorithm.propagate_cover_dest(one_source_poi_list, &destinations_highway_node_list);

        // add highway node in of current box
        for (auto i : highway_node_id_with_box_id_map[box_id])
        {
            geodesic::SurfacePoint one_dest(&mesh->vertices()[i.first]);
            std::vector<geodesic::SurfacePoint> path;
            algorithm.trace_back(one_dest, path);

            int x_vertex_id = poi_list[k];
            int y_vertex_id = i.first;
            int x_y_vertex_id;
            if (x_vertex_id <= y_vertex_id)
            {
                hash_function_two_keys_to_one_key(mesh->vertices().size(), x_vertex_id, y_vertex_id, x_y_vertex_id);
            }
            else
            {
                hash_function_two_keys_to_one_key(mesh->vertices().size(), y_vertex_id, x_vertex_id, x_y_vertex_id);
            }
            if (path_poi_to_highway_node_map.count(x_y_vertex_id) == 0)
            {
                path_poi_to_highway_node_map[x_y_vertex_id] = path;
                distance_poi_to_highway_node_map[x_y_vertex_id] = length(path);
                EAR_Oracle_size += path.size() * sizeof(geodesic::SurfacePoint);
                memory_usage += path.size() * sizeof(geodesic::SurfacePoint);
            }
        }
        // add poi in the current box
        for (int i = 0; i < same_box_poi_vertex_id.size(); i++)
        {
            geodesic::SurfacePoint one_dest(&mesh->vertices()[same_box_poi_vertex_id[i]]);
            std::vector<geodesic::SurfacePoint> path;
            algorithm.trace_back(one_dest, path);

            int x_vertex_id = poi_list[k];
            int y_vertex_id = same_box_poi_vertex_id[i];
            int x_y_vertex_id;
            if (x_vertex_id <= y_vertex_id)
            {
                hash_function_two_keys_to_one_key(mesh->vertices().size(), x_vertex_id, y_vertex_id, x_y_vertex_id);
            }
            else
            {
                hash_function_two_keys_to_one_key(mesh->vertices().size(), y_vertex_id, x_vertex_id, x_y_vertex_id);
            }
            if (path_poi_to_highway_node_map.count(x_y_vertex_id) == 0)
            {
                path_poi_to_highway_node_map[x_y_vertex_id] = path;
                distance_poi_to_highway_node_map[x_y_vertex_id] = length(path);
                EAR_Oracle_size += path.size() * sizeof(geodesic::SurfacePoint);
                memory_usage += path.size() * sizeof(geodesic::SurfacePoint);
            }
        }
    }
}

void EAR_Oracle_query(geodesic::Mesh *mesh, std::vector<int> &poi_list,
                      int sqrt_num_of_box, int geo_tree_node_id, std::unordered_map<int, std::unordered_map<int, int>> &highway_node_id_with_box_id_map,
                      std::vector<GeoNode *> &all_highway_node, std::unordered_map<int, GeoNode *> &geo_node_in_partition_tree_unordered_map,
                      std::unordered_map<int, GeoPair *> &geopairs,
                      std::unordered_map<int, double> &distance_poi_to_highway_node_map, std::unordered_map<int, std::vector<geodesic::SurfacePoint>> &path_poi_to_highway_node_map,
                      int source_poi_index, int destination_poi_index, double &approximate_distance, std::vector<geodesic::SurfacePoint> &approximate_path,
                      std::vector<int> &face_sequence_index_list)
{
    bool one_path;

    double box_width = mesh->m_width / sqrt_num_of_box;
    double box_height = mesh->m_height / sqrt_num_of_box;
    int src_box_id = 0;
    int dest_box_id = 0;
    double src_box_x_min, dest_box_x_min;
    double src_box_x_max, dest_box_x_max;
    double src_box_y_min, dest_box_y_min;
    double src_box_y_max, dest_box_y_max;

    double src_x = mesh->vertices()[poi_list[source_poi_index]].getx();
    double src_y = mesh->vertices()[poi_list[source_poi_index]].gety();
    double dest_x = mesh->vertices()[poi_list[destination_poi_index]].getx();
    double dest_y = mesh->vertices()[poi_list[destination_poi_index]].gety();

    // for source poi
    for (int i = 0; i < sqrt_num_of_box; i++)
    {
        for (int j = 0; j < sqrt_num_of_box; j++)
        {
            src_box_x_min = mesh->m_xmin + i * box_width;
            src_box_x_max = mesh->m_xmin + (i + 1) * box_width;
            src_box_y_min = mesh->m_ymin + j * box_height;
            src_box_y_max = mesh->m_ymin + (j + 1) * box_height;

            if (src_x >= src_box_x_min && src_x <= src_box_x_max & src_y >= src_box_y_min && src_y <= src_box_y_max)
            {
                i = sqrt_num_of_box;
                j = sqrt_num_of_box;
                continue;
            }
            src_box_id++;
        }
    }

    // for destination poi
    for (int i = 0; i < sqrt_num_of_box; i++)
    {
        for (int j = 0; j < sqrt_num_of_box; j++)
        {
            dest_box_x_min = mesh->m_xmin + i * box_width;
            dest_box_x_max = mesh->m_xmin + (i + 1) * box_width;
            dest_box_y_min = mesh->m_ymin + j * box_height;
            dest_box_y_max = mesh->m_ymin + (j + 1) * box_height;

            if (dest_x >= dest_box_x_min && dest_x <= dest_box_x_max & dest_y >= dest_box_y_min && dest_y <= dest_box_y_max)
            {
                i = sqrt_num_of_box;
                j = sqrt_num_of_box;
                continue;
            }
            dest_box_id++;
        }
    }

    if (src_box_id == dest_box_id)
    {
        int x_vertex_id = poi_list[source_poi_index];
        int y_vertex_id = poi_list[destination_poi_index];
        int x_y_vertex_id;
        if (x_vertex_id <= y_vertex_id)
        {
            hash_function_two_keys_to_one_key(mesh->vertices().size(), x_vertex_id, y_vertex_id, x_y_vertex_id);
        }
        else
        {
            hash_function_two_keys_to_one_key(mesh->vertices().size(), y_vertex_id, x_vertex_id, x_y_vertex_id);
        }

        assert(path_poi_to_highway_node_map.count(x_y_vertex_id) != 0);
        approximate_distance = distance_poi_to_highway_node_map[x_y_vertex_id];
        approximate_path = path_poi_to_highway_node_map[x_y_vertex_id];
    }
    else
    {
        double min_distance = 1e10;
        std::vector<geodesic::SurfacePoint> min_src_to_highway_node_path;
        min_src_to_highway_node_path.clear();
        std::vector<geodesic::SurfacePoint> min_dest_to_highway_node_path;
        min_dest_to_highway_node_path.clear();

        for (auto i : highway_node_id_with_box_id_map[src_box_id])
        {
            for (auto j : highway_node_id_with_box_id_map[dest_box_id])
            {
                // highway node to highway node
                double highway_node_to_highway_node_distance = 0;
                std::vector<geodesic::SurfacePoint> highway_node_to_highway_node_path;
                highway_node_to_highway_node_path.clear();
                std::vector<int> highway_node_to_highway_node_face_sequence_index_list;
                if (i.first != j.first)
                {
                    three_paths_query_geo(mesh, geo_tree_node_id, 2, i.first, j.first, one_path, all_highway_node, geo_node_in_partition_tree_unordered_map, geopairs,
                                          highway_node_to_highway_node_distance, highway_node_to_highway_node_path, highway_node_to_highway_node_face_sequence_index_list);
                }
                else
                {
                    highway_node_to_highway_node_path.push_back(geodesic::SurfacePoint(&mesh->vertices()[i.first]));
                }

                // src to highway node
                double src_to_highway_node_distance = 0;
                std::vector<geodesic::SurfacePoint> src_to_highway_node_path;
                src_to_highway_node_path.clear();

                int x_vertex_id1 = poi_list[source_poi_index];
                int y_vertex_id1 = i.first;
                int x_y_vertex_id1;
                if (x_vertex_id1 <= y_vertex_id1)
                {
                    hash_function_two_keys_to_one_key(mesh->vertices().size(), x_vertex_id1, y_vertex_id1, x_y_vertex_id1);
                }
                else
                {
                    hash_function_two_keys_to_one_key(mesh->vertices().size(), y_vertex_id1, x_vertex_id1, x_y_vertex_id1);
                }
                assert(path_poi_to_highway_node_map.count(x_y_vertex_id1) != 0);
                src_to_highway_node_distance = distance_poi_to_highway_node_map[x_y_vertex_id1];
                src_to_highway_node_path = path_poi_to_highway_node_map[x_y_vertex_id1];
                if (src_to_highway_node_path.size() == 0)
                {
                    src_to_highway_node_distance = 0;
                    src_to_highway_node_path.push_back(geodesic::SurfacePoint(&mesh->vertices()[i.first]));
                }

                // dest to highway node
                double dest_to_highway_node_distance = 0;
                std::vector<geodesic::SurfacePoint> dest_to_highway_node_path;
                dest_to_highway_node_path.clear();

                int x_vertex_id2 = poi_list[destination_poi_index];
                int y_vertex_id2 = j.first;
                int x_y_vertex_id2;
                if (x_vertex_id2 <= y_vertex_id2)
                {
                    hash_function_two_keys_to_one_key(mesh->vertices().size(), x_vertex_id2, y_vertex_id2, x_y_vertex_id2);
                }
                else
                {
                    hash_function_two_keys_to_one_key(mesh->vertices().size(), y_vertex_id2, x_vertex_id2, x_y_vertex_id2);
                }
                assert(path_poi_to_highway_node_map.count(x_y_vertex_id2) != 0);
                dest_to_highway_node_distance = distance_poi_to_highway_node_map[x_y_vertex_id2];
                dest_to_highway_node_path = path_poi_to_highway_node_map[x_y_vertex_id2];
                if (dest_to_highway_node_path.size() == 0)
                {
                    dest_to_highway_node_distance = 0;
                    dest_to_highway_node_path.push_back(geodesic::SurfacePoint(&mesh->vertices()[j.first]));
                }

                if (min_distance > highway_node_to_highway_node_distance + src_to_highway_node_distance + dest_to_highway_node_distance)
                {
                    min_distance = highway_node_to_highway_node_distance + src_to_highway_node_distance + dest_to_highway_node_distance;
                    min_src_to_highway_node_path = src_to_highway_node_path;
                    min_dest_to_highway_node_path = dest_to_highway_node_path;
                    approximate_path = highway_node_to_highway_node_path;
                    approximate_distance = min_distance;
                }
            }
        }

        if (min_dest_to_highway_node_path[0].getx() == approximate_path[0].getx() &&
            min_dest_to_highway_node_path[0].gety() == approximate_path[0].gety() &&
            min_dest_to_highway_node_path[0].getz() == approximate_path[0].getz())
        {
            std::reverse(approximate_path.begin(), approximate_path.end());
            for (int i = 0; i < min_dest_to_highway_node_path.size(); i++)
            {
                approximate_path.push_back(min_dest_to_highway_node_path[i]);
            }
        }
        else if (min_dest_to_highway_node_path[0].getx() == approximate_path[approximate_path.size() - 1].getx() &&
                 min_dest_to_highway_node_path[0].gety() == approximate_path[approximate_path.size() - 1].gety() &&
                 min_dest_to_highway_node_path[0].getz() == approximate_path[approximate_path.size() - 1].getz())
        {
            for (int i = 0; i < min_dest_to_highway_node_path.size(); i++)
            {
                approximate_path.push_back(min_dest_to_highway_node_path[i]);
            }
        }
        else if (min_dest_to_highway_node_path[min_dest_to_highway_node_path.size() - 1].getx() == approximate_path[0].getx() &&
                 min_dest_to_highway_node_path[min_dest_to_highway_node_path.size() - 1].gety() == approximate_path[0].gety() &&
                 min_dest_to_highway_node_path[min_dest_to_highway_node_path.size() - 1].getz() == approximate_path[0].getz())
        {
            std::reverse(approximate_path.begin(), approximate_path.end());
            for (int i = min_dest_to_highway_node_path.size() - 1; i >= 0; i--)
            {
                approximate_path.push_back(min_dest_to_highway_node_path[i]);
            }
        }
        else if (min_dest_to_highway_node_path[min_dest_to_highway_node_path.size() - 1].getx() == approximate_path[approximate_path.size() - 1].getx() &&
                 min_dest_to_highway_node_path[min_dest_to_highway_node_path.size() - 1].gety() == approximate_path[approximate_path.size() - 1].gety() &&
                 min_dest_to_highway_node_path[min_dest_to_highway_node_path.size() - 1].getz() == approximate_path[approximate_path.size() - 1].getz())
        {
            for (int i = min_dest_to_highway_node_path.size() - 1; i >= 0; i--)
            {
                approximate_path.push_back(min_dest_to_highway_node_path[i]);
            }
        }
        else
        {
            assert(false);
        }

        if (min_src_to_highway_node_path[0].getx() == approximate_path[0].getx() &&
            min_src_to_highway_node_path[0].gety() == approximate_path[0].gety() &&
            min_src_to_highway_node_path[0].getz() == approximate_path[0].getz())
        {
            std::reverse(approximate_path.begin(), approximate_path.end());
            for (int i = 0; i < min_src_to_highway_node_path.size(); i++)
            {
                approximate_path.push_back(min_src_to_highway_node_path[i]);
            }
        }
        else if (min_src_to_highway_node_path[0].getx() == approximate_path[approximate_path.size() - 1].getx() &&
                 min_src_to_highway_node_path[0].gety() == approximate_path[approximate_path.size() - 1].gety() &&
                 min_src_to_highway_node_path[0].getz() == approximate_path[approximate_path.size() - 1].getz())
        {
            for (int i = 0; i < min_src_to_highway_node_path.size(); i++)
            {
                approximate_path.push_back(min_src_to_highway_node_path[i]);
            }
        }
        else if (min_src_to_highway_node_path[min_src_to_highway_node_path.size() - 1].getx() == approximate_path[0].getx() &&
                 min_src_to_highway_node_path[min_src_to_highway_node_path.size() - 1].gety() == approximate_path[0].gety() &&
                 min_src_to_highway_node_path[min_src_to_highway_node_path.size() - 1].getz() == approximate_path[0].getz())
        {
            std::reverse(approximate_path.begin(), approximate_path.end());
            for (int i = min_src_to_highway_node_path.size() - 1; i >= 0; i--)
            {
                approximate_path.push_back(min_src_to_highway_node_path[i]);
            }
        }
        else if (min_src_to_highway_node_path[min_src_to_highway_node_path.size() - 1].getx() == approximate_path[approximate_path.size() - 1].getx() &&
                 min_src_to_highway_node_path[min_src_to_highway_node_path.size() - 1].gety() == approximate_path[approximate_path.size() - 1].gety() &&
                 min_src_to_highway_node_path[min_src_to_highway_node_path.size() - 1].getz() == approximate_path[approximate_path.size() - 1].getz())
        {
            for (int i = min_src_to_highway_node_path.size() - 1; i >= 0; i--)
            {
                approximate_path.push_back(min_src_to_highway_node_path[i]);
            }
        }
        else
        {
            assert(false);
        }
    }
    get_face_sequence(approximate_path, face_sequence_index_list);
}

double euclidean_distance(double x_1, double y_1, double x_2, double y_2)
{
    return sqrt(pow(x_1 - x_2, 2) + pow(y_1 - y_2, 2));
}

double epslion_to_subdivision_level(double epsilon)
{
    assert(epsilon > 0 && epsilon <= 1);
    double subdivision_level = 0;
    if (epsilon > 0 && epsilon <= 0.05)
    {
        subdivision_level = floor(7 / (10 * epsilon)) - 8;
    }
    else if (epsilon > 0.05 && epsilon <= 0.1)
    {
        subdivision_level = floor(7 / (10 * epsilon)) - 2;
    }
    else if (epsilon > 0.1 && epsilon <= 0.25)
    {
        subdivision_level = floor(1 / epsilon) - 1;
    }
    else if (epsilon > 0.25 && epsilon < 1)
    {
        subdivision_level = floor(1 / epsilon);
    }
    if (subdivision_level < 5)
    {
        subdivision_level++;
    }
    std::cout << "subdivision_level: " << subdivision_level << std::endl;
    return subdivision_level;
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
