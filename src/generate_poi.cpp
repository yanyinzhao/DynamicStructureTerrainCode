#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "algorithms.h"

struct input_struct
{
    int poi_num;
    int face_num;
    int discrete_poi;
    int six_neighbour_poi_num; // five_neighbour_poi_num = discrete_poi - six_neighbour_poi_num
    int poi_index_changed_num;
    int poi_in_changed_area_num;
    int poi_path_pass_changed_area_num;
    int poi_path_close_changed_area_num;

    input_struct() {}

    input_struct(int _poi_num, int _face_num,
                 int _discrete_poi, int _six_neighbour_poi_num,
                 int _poi_index_changed_num, int _poi_in_changed_area_num,
                 int _poi_path_pass_changed_area_num, int _poi_path_close_changed_area_num)
    {
        poi_num = _poi_num;
        face_num = _face_num;
        discrete_poi = _discrete_poi;
        six_neighbour_poi_num = _six_neighbour_poi_num;
        poi_index_changed_num = _poi_index_changed_num;
        poi_in_changed_area_num = _poi_in_changed_area_num;
        poi_path_pass_changed_area_num = _poi_path_pass_changed_area_num;
        poi_path_close_changed_area_num = _poi_path_close_changed_area_num;
    }
};

int main(int argc, char **argv)
{
    std::vector<std::pair<std::string, std::string>> dataset;
    dataset.push_back(std::make_pair("SCpre", "SCpost"));
    dataset.push_back(std::make_pair("AUpre", "AUpost"));
    dataset.push_back(std::make_pair("VSpre", "VSpost"));

    std::map<int, input_struct> input{
        // {0, input_struct(50, 10082, 9, 0, 1, 2, 1, 1)}, // delete this
        {0, input_struct(50, 10082, 9, 4, 1, 0, 0, 0)},
        {1, input_struct(100, 10082, 18, 8, 2, 0, 0, 0)},
        {2, input_struct(150, 10082, 27, 12, 3, 0, 0, 0)},
        {3, input_struct(200, 10082, 36, 16, 4, 0, 0, 0)},
        {4, input_struct(250, 10082, 45, 20, 5, 0, 0, 0)},
        {5, input_struct(50, 20000, 9, 4, 1, 0, 0, 0)},
        {6, input_struct(50, 30258, 9, 4, 1, 0, 0, 0)},
        {7, input_struct(50, 40328, 9, 4, 1, 0, 0, 0)},
        {8, input_struct(50, 50562, 9, 4, 1, 0, 0, 0)},
        {9, input_struct(500, 500000, 90, 45, 1, 3, 1, 0)},
        {10, input_struct(1000, 500000, 180, 90, 1, 8, 1, 0)},
        {11, input_struct(1500, 500000, 270, 135, 1, 13, 1, 0)},
        {12, input_struct(2000, 500000, 360, 180, 1, 18, 1, 0)},
        {13, input_struct(2500, 500000, 450, 225, 1, 23, 1, 0)},
        {14, input_struct(500, 1002528, 90, 45, 1, 3, 1, 0)},
        {15, input_struct(500, 1503378, 90, 45, 1, 3, 1, 0)},
        {16, input_struct(500, 2000000, 90, 45, 1, 3, 1, 0)},
        {17, input_struct(500, 2504322, 90, 45, 1, 3, 1, 0)},
    };

    for (int k = 0; k < input.size(); k++)
    {
        for (int m = 0; m < dataset.size(); m++)
        {
            if (k == 5 || k == 6 || k == 7 || k == 8)
            {
                if (m == 1 || m == 2)
                {
                    continue;
                }
            }

            std::string input_folder = "../exp_input/";
            std::string output_folder = "../exp_input_poi/";
            // std::string output_folder = "../input/";
            std::string input_file_type = ".off";
            std::string output_file_type = ".txt";

            int poi_num = input[k].poi_num;
            int face_num = input[k].face_num;

            double nonchanged_area_over_changed_area_ratio = 0.5;

            int x_vertex_num = sqrt(face_num / 2) + 1;
            int y_vertex_num = sqrt(face_num / 2) + 1;

            std::cout << output_folder + dataset[m].first + "_" + std::to_string(poi_num) + "_poi_on_" + std::to_string(face_num) << std::endl;
            std::cout << output_folder + dataset[m].second + "_" + std::to_string(poi_num) + "_poi_on_" + std::to_string(face_num) << std::endl;

            // std::string read_path = input_folder + dataset + "_" + std::to_string(face_num) + input_file_type;

            std::unordered_map<int, int> pre_poi_list;
            std::unordered_map<int, int> post_poi_list;

            std::vector<int> pre_poi_vector;
            std::vector<int> post_poi_vector;

            int count = 0;
            while (count < input[k].discrete_poi - input[k].six_neighbour_poi_num)
            {
                int b = rand() % (int(y_vertex_num * nonchanged_area_over_changed_area_ratio) - 8) + 4; // 4 to y_vertex_num * nonchanged_area_over_changed_area_ratio - 4
                int a = rand() % (int(x_vertex_num * nonchanged_area_over_changed_area_ratio) - 8) + 4; // 4 to x_vertex_num * nonchanged_area_over_changed_area_ratio - 4
                // std::cout << "a: " << a << ", b: " << b << std::endl;
                if (a + b > 1.5 * y_vertex_num * nonchanged_area_over_changed_area_ratio - 8)
                {
                    continue;
                }
                int center_index = a + b * x_vertex_num;
                int top_index = a + (b - 1) * x_vertex_num;
                int bottom_index = a + (b + 1) * x_vertex_num;
                int left_index = (a - 1) + b * x_vertex_num;
                int right_index = (a + 1) + b * x_vertex_num;
                // std::cout << center_index << " " << top_index << " " << bottom_index << " " << left_index << " " << right_index << std::endl;

                if (pre_poi_list.count(center_index) == 0 &&
                    pre_poi_list.count(top_index) == 0 &&
                    pre_poi_list.count(bottom_index) == 0 &&
                    pre_poi_list.count(left_index) == 0 &&
                    pre_poi_list.count(right_index) == 0 &&
                    post_poi_list.count(center_index) == 0 &&
                    post_poi_list.count(top_index) == 0 &&
                    post_poi_list.count(bottom_index) == 0 &&
                    post_poi_list.count(left_index) == 0 &&
                    post_poi_list.count(right_index) == 0)
                {
                    // std::cout << "not redundent" << std::endl;
                    count++;
                    pre_poi_list[center_index] = center_index;
                    pre_poi_list[top_index] = top_index;
                    pre_poi_list[bottom_index] = bottom_index;
                    pre_poi_list[left_index] = left_index;
                    pre_poi_list[right_index] = right_index;
                    post_poi_list[center_index] = center_index;
                    post_poi_list[top_index] = top_index;
                    post_poi_list[bottom_index] = bottom_index;
                    post_poi_list[left_index] = left_index;
                    post_poi_list[right_index] = right_index;

                    pre_poi_vector.push_back(center_index);
                    pre_poi_vector.push_back(top_index);
                    pre_poi_vector.push_back(bottom_index);
                    pre_poi_vector.push_back(left_index);
                    pre_poi_vector.push_back(right_index);
                    post_poi_vector.push_back(center_index);
                    post_poi_vector.push_back(top_index);
                    post_poi_vector.push_back(bottom_index);
                    post_poi_vector.push_back(left_index);
                    post_poi_vector.push_back(right_index);
                }
            }

            count = 0;
            while (count < input[k].six_neighbour_poi_num)
            {
                int b = rand() % (int(y_vertex_num * nonchanged_area_over_changed_area_ratio) - 8) + 4; // 4 to y_vertex_num * nonchanged_area_over_changed_area_ratio - 4
                int a = rand() % (int(x_vertex_num * nonchanged_area_over_changed_area_ratio) - 8) + 4; // 4 to x_vertex_num * nonchanged_area_over_changed_area_ratio - 4
                // std::cout << "a: " << a << ", b: " << b << std::endl;
                if (a + b > 1.5 * y_vertex_num * nonchanged_area_over_changed_area_ratio - 8)
                {
                    continue;
                }
                int center_index = a + b * x_vertex_num;
                int top_index = a + (b - 1) * x_vertex_num;
                int bottom_index = a + (b + 1) * x_vertex_num;
                int left_index = (a - 1) + b * x_vertex_num;
                int right_index = (a + 1) + b * x_vertex_num;
                int top_right_index = (a + 1) + (b - 1) * x_vertex_num;
                // std::cout << center_index << " " << top_index << " " << bottom_index << " " << left_index << " " << right_index << " " << top_right_index << std::endl;

                if (pre_poi_list.count(center_index) == 0 &&
                    pre_poi_list.count(top_index) == 0 &&
                    pre_poi_list.count(bottom_index) == 0 &&
                    pre_poi_list.count(left_index) == 0 &&
                    pre_poi_list.count(right_index) == 0 &&
                    pre_poi_list.count(top_right_index) == 0 &&
                    post_poi_list.count(center_index) == 0 &&
                    post_poi_list.count(top_index) == 0 &&
                    post_poi_list.count(bottom_index) == 0 &&
                    post_poi_list.count(left_index) == 0 &&
                    post_poi_list.count(right_index) == 0 &&
                    post_poi_list.count(top_right_index) == 0)
                {
                    // std::cout << "not redundent" << std::endl;

                    count++;
                    pre_poi_list[center_index] = center_index;
                    pre_poi_list[top_index] = top_index;
                    pre_poi_list[bottom_index] = bottom_index;
                    pre_poi_list[left_index] = left_index;
                    pre_poi_list[right_index] = right_index;
                    pre_poi_list[top_right_index] = top_right_index;
                    post_poi_list[center_index] = center_index;
                    post_poi_list[top_index] = top_index;
                    post_poi_list[bottom_index] = bottom_index;
                    post_poi_list[left_index] = left_index;
                    post_poi_list[right_index] = right_index;
                    post_poi_list[top_right_index] = top_right_index;

                    pre_poi_vector.push_back(center_index);
                    pre_poi_vector.push_back(top_index);
                    pre_poi_vector.push_back(bottom_index);
                    pre_poi_vector.push_back(left_index);
                    pre_poi_vector.push_back(right_index);
                    pre_poi_vector.push_back(top_right_index);
                    post_poi_vector.push_back(center_index);
                    post_poi_vector.push_back(top_index);
                    post_poi_vector.push_back(bottom_index);
                    post_poi_vector.push_back(left_index);
                    post_poi_vector.push_back(right_index);
                    post_poi_vector.push_back(top_right_index);
                }
            }

            std::cout << "1: " << pre_poi_list.size() << " " << post_poi_list.size() << std::endl;

            std::string pre_input_file_path = input_folder + dataset[m].first + "_" + std::to_string(face_num) + input_file_type;
            std::string post_input_file_path = input_folder + dataset[m].second + "_" + std::to_string(face_num) + input_file_type;

            std::vector<double> pre_points;
            std::vector<unsigned> pre_faces;
            geodesic::read_mesh_from_file(&pre_input_file_path[0], pre_points, pre_faces);
            geodesic::Mesh pre_mesh;
            pre_mesh.initialize_mesh_data(pre_points, pre_faces);
            std::vector<double> post_points;
            std::vector<unsigned> post_faces;
            geodesic::read_mesh_from_file(&post_input_file_path[0], post_points, post_faces);
            geodesic::Mesh post_mesh;
            post_mesh.initialize_mesh_data(post_points, post_faces);

            std::unordered_map<int, int> changed_vertex_index_unordered_map;
            changed_vertex_index_unordered_map.clear();

            assert(pre_mesh.vertices().size() == post_mesh.vertices().size());
            for (int i = 0; i < pre_mesh.vertices().size(); i++)
            {
                if (pre_mesh.vertices()[i].x() != post_mesh.vertices()[i].x() ||
                    pre_mesh.vertices()[i].y() != post_mesh.vertices()[i].y() ||
                    pre_mesh.vertices()[i].z() != post_mesh.vertices()[i].z())
                {
                    changed_vertex_index_unordered_map[pre_mesh.vertices()[i].id()] = pre_mesh.vertices()[i].id();
                }
            }

            assert(input[k].poi_in_changed_area_num <= changed_vertex_index_unordered_map.size());

            // POI in the changed area
            count = 0;
            for (auto it = changed_vertex_index_unordered_map.begin(); it != changed_vertex_index_unordered_map.end(); ++it)
            {
                if (count >= input[k].poi_in_changed_area_num)
                {
                    break;
                }
                if (pre_poi_list.count(it->second) == 0 &&
                    post_poi_list.count(it->second) == 0)
                {
                    count++;
                    // std::cout << "changed POI: " << it->second << std::endl;
                    pre_poi_list[it->second] = it->second;
                    post_poi_list[it->second] = it->second;

                    pre_poi_vector.push_back(it->second);
                    post_poi_vector.push_back(it->second);
                }
            }

            std::cout << "2: " << pre_poi_list.size() << " " << post_poi_list.size() << std::endl;

            // POI path pass changed area or close to changed area
            count = 0;
            while (count < input[k].poi_path_close_changed_area_num + input[k].poi_path_pass_changed_area_num)
            {
                int b = rand() % (int(y_vertex_num * (1 - nonchanged_area_over_changed_area_ratio))) + y_vertex_num * nonchanged_area_over_changed_area_ratio; // y_vertex_num * nonchanged_area_over_changed_area_ratio to y_vertex_num
                int a = rand() % (int(x_vertex_num * (1 - nonchanged_area_over_changed_area_ratio))) + x_vertex_num * nonchanged_area_over_changed_area_ratio; // x_vertex_num * nonchanged_area_over_changed_area_ratio to x_vertex_num
                int center_index = a + b * x_vertex_num;
                if (changed_vertex_index_unordered_map.count(center_index) != 0)
                {
                    continue;
                }
                if (pre_poi_list.count(center_index) == 0 &&
                    post_poi_list.count(center_index) == 0)
                {
                    count++;
                    pre_poi_list[center_index] = center_index;
                    post_poi_list[center_index] = center_index;

                    pre_poi_vector.push_back(center_index);
                    post_poi_vector.push_back(center_index);
                }
            }

            std::cout << "3: " << pre_poi_list.size() << " " << post_poi_list.size() << std::endl;

            // POI index changed
            count = 0;
            while (count < input[k].poi_index_changed_num)
            {
                int b = rand() % (int(y_vertex_num * (1 - nonchanged_area_over_changed_area_ratio) - 4)) + y_vertex_num * nonchanged_area_over_changed_area_ratio + 2; // y_vertex_num * nonchanged_area_over_changed_area_ratio + 2 to y_vertex_num - 2
                int a = rand() % (int(x_vertex_num * (1 - nonchanged_area_over_changed_area_ratio) - 4)) + x_vertex_num * nonchanged_area_over_changed_area_ratio + 2; // x_vertex_num * nonchanged_area_over_changed_area_ratio + 2 to x_vertex_num - 2
                int center_index = a + b * x_vertex_num;
                int left_index = (a - 1) + b * x_vertex_num;
                if (changed_vertex_index_unordered_map.count(center_index) != 0)
                {
                    continue;
                }
                if (pre_poi_list.count(center_index) == 0 &&
                    post_poi_list.count(left_index) == 0)
                {
                    count++;
                    pre_poi_list[center_index] = center_index;
                    post_poi_list[left_index] = left_index;

                    pre_poi_vector.push_back(center_index);
                    post_poi_vector.push_back(left_index);
                }
            }

            std::cout << "4: " << pre_poi_list.size() << " " << post_poi_list.size() << std::endl;
            std::cout << std::endl;

            assert(pre_poi_list.size() == poi_num && post_poi_list.size() == poi_num);

            std::string pre_poi_write_path = output_folder + dataset[m].first + "_" + std::to_string(poi_num) + "_poi_on_" + std::to_string(face_num) + output_file_type;
            std::string post_poi_write_path = output_folder + dataset[m].second + "_" + std::to_string(poi_num) + "_poi_on_" + std::to_string(face_num) + output_file_type;

            std::ofstream pre_ofs(&pre_poi_write_path[0], std::ofstream::trunc);
            pre_ofs << std::to_string(poi_num) << "\n";
            // for (auto it = pre_poi_list.begin(); it != pre_poi_list.end(); ++it)
            // {
            //     pre_ofs << it->second << " ";
            // }
            for (int i = 0; i < pre_poi_vector.size(); i++)
            {
                pre_ofs << pre_poi_vector[i] << " ";
            }
            pre_ofs.close();

            std::ofstream post_ofs(&post_poi_write_path[0], std::ofstream::trunc);
            post_ofs << std::to_string(poi_num) << "\n";
            // for (auto it = post_poi_list.begin(); it != post_poi_list.end(); ++it)
            // {
            //     post_ofs << it->second << " ";
            // }
            for (int i = 0; i < post_poi_vector.size(); i++)
            {
                post_ofs << post_poi_vector[i] << " ";
            }
            post_ofs.close();
        }
    }
    return 0;
}
