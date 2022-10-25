#include <iostream>
#include <fstream>
#include <string>

#include "geodesic_algorithm_subdivision.h"

std::vector<std::string> split(const std::string &str, char delim)
{
    std::vector<std::string> strings;
    size_t start;
    size_t end = 0;
    while ((start = str.find_first_not_of(delim, end)) != std::string::npos)
    {
        end = str.find(delim, start);
        strings.push_back(str.substr(start, end - start));
    }
    return strings;
}

double area(double x1, double y1, double x2, double y2, double x3, double y3)
{
    return abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
}

bool point_inside_or_on_edge_tri(double x, double y, double x1, double y1, double x2, double y2, double x3, double y3)
{
    double A = area(x1, y1, x2, y2, x3, y3);
    double A1 = area(x, y, x2, y2, x3, y3);
    double A2 = area(x1, y1, x, y, x3, y3);
    double A3 = area(x1, y1, x2, y2, x, y);
    return (A == A1 + A2 + A3);
}

// double print_area(double x1, double y1, double x2, double y2, double x3, double y3)
// {
//     std::cout << "abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)): " << abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) << std::endl;
//     return abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
// }

// void print_point_inside_or_on_edge_tri(double x, double y, double x1, double y1, double x2, double y2, double x3, double y3)
// {
//     double A = print_area(x1, y1, x2, y2, x3, y3);
//     double A1 = print_area(x, y, x2, y2, x3, y3);
//     double A2 = print_area(x1, y1, x, y, x3, y3);
//     double A3 = print_area(x1, y1, x2, y2, x, y);
//     std::cout << "A1 + A2 + A3: " << A1 + A2 + A3 << std::endl;
//     std::cout << "A == A1 + A2 + A3: " << (A == A1 + A2 + A3) << std::endl;
// }

int main(int argc, char **argv)
{

    // **** Change here ****
    std::string dataset = "AUpre";
    // **** Change here ****

    std::vector<int> intend_face_num_list;

    if (dataset.compare("SCpre") == 0)
    {
        intend_face_num_list.push_back(10000);
        intend_face_num_list.push_back(20000);
        intend_face_num_list.push_back(30000);
        intend_face_num_list.push_back(40000);
        intend_face_num_list.push_back(50000);
        intend_face_num_list.push_back(500000);
        intend_face_num_list.push_back(1000000);
        intend_face_num_list.push_back(1500000);
        intend_face_num_list.push_back(2000000);
        intend_face_num_list.push_back(2500000);
    }
    else if (dataset.compare("SCpost") == 0)
    {
        intend_face_num_list.push_back(10000);
        intend_face_num_list.push_back(20000);
        intend_face_num_list.push_back(30000);
        intend_face_num_list.push_back(40000);
        intend_face_num_list.push_back(50000);
        intend_face_num_list.push_back(500000);
        intend_face_num_list.push_back(1000000);
        intend_face_num_list.push_back(1500000);
        intend_face_num_list.push_back(2000000);
        intend_face_num_list.push_back(2500000);
    }
    else if (dataset.compare("AUpre") == 0)
    {
        intend_face_num_list.push_back(10000);
        intend_face_num_list.push_back(500000);
        intend_face_num_list.push_back(1000000);
        intend_face_num_list.push_back(1500000);
        intend_face_num_list.push_back(2000000);
        intend_face_num_list.push_back(2500000);
    }
    else if (dataset.compare("AUpost") == 0)
    {
        intend_face_num_list.push_back(10000);
        intend_face_num_list.push_back(500000);
        intend_face_num_list.push_back(1000000);
        intend_face_num_list.push_back(1500000);
        intend_face_num_list.push_back(2000000);
        intend_face_num_list.push_back(2500000);
    }
    else if (dataset.compare("VSpre") == 0)
    {
        intend_face_num_list.push_back(10000);
        intend_face_num_list.push_back(500000);
        intend_face_num_list.push_back(1000000);
        intend_face_num_list.push_back(1500000);
        intend_face_num_list.push_back(2000000);
        intend_face_num_list.push_back(2500000);
    }
    else if (dataset.compare("VSpost") == 0)
    {
        intend_face_num_list.push_back(10000);
        intend_face_num_list.push_back(500000);
        intend_face_num_list.push_back(1000000);
        intend_face_num_list.push_back(1500000);
        intend_face_num_list.push_back(2000000);
        intend_face_num_list.push_back(2500000);
    }

    for (int k = 0; k < intend_face_num_list.size(); ++k)
    {
        int intend_face_num = intend_face_num_list[k];

        std::string input_folder = "../input/";
        std::string output_folder = "../exp_input/";
        std::string file_type = ".off";

        std::string read_path = input_folder + dataset + file_type;

        std::vector<double> points;
        std::vector<unsigned> faces;

        geodesic::read_mesh_from_file(&read_path[0], points, faces);
        geodesic::Mesh mesh;
        mesh.initialize_mesh_data(points, faces);

        std::ifstream ifs;
        ifs.open(&read_path[0]);
        std::string oneline;
        std::vector<std::string> oneword;

        int total_original_vertex_num = 0;
        int count = 0;

        double x_min = 1e10;
        double x_max = 0;
        double y_min = 1e10;
        double y_max = 0;

        if (ifs.is_open())
        {
            std::getline(ifs, oneline);
            std::getline(ifs, oneline);
            oneword = split(oneline, ' ');
            total_original_vertex_num = std::stoi(oneword[0]);
            std::getline(ifs, oneline);
            while (ifs)
            {
                if (count == total_original_vertex_num)
                {
                    break;
                }

                oneword = split(oneline, ' ');
                double x = std::stod(oneword[0]);
                double y = std::stod(oneword[1]);

                x_min = std::min(x, x_min);
                x_max = std::max(x, x_max);
                y_min = std::min(y, y_min);
                y_max = std::max(y, y_max);

                std::getline(ifs, oneline);
                count++;
            }
        }
        double x_length = x_max - x_min;
        double y_length = y_max - y_min;

        std::cout << "make sure (x_new_vertex_num - 1)/(y_new_vertex_num - 1) is x/y ratio: " << x_length / y_length << std::endl;

        int x_new_vertex_num = (int)ceil(sqrt(intend_face_num / 2 * x_length / y_length)) + 1;
        int y_new_vertex_num = (int)ceil((y_length * sqrt(intend_face_num / 2 * x_length / y_length) / x_length)) + 1;

        std::cout << "x_new_vertex_num: " << x_new_vertex_num << std::endl;
        std::cout << "y_new_vertex_num: " << y_new_vertex_num << std::endl;

        int total_new_vertex_num = x_new_vertex_num * y_new_vertex_num;
        int total_new_face_num = (x_new_vertex_num - 1) * (y_new_vertex_num - 1) * 2;
        int total_new_edge_num = (x_new_vertex_num - 1) * y_new_vertex_num + x_new_vertex_num * (y_new_vertex_num - 1) + (x_new_vertex_num - 1) * (y_new_vertex_num - 1);

        assert(total_new_face_num >= intend_face_num);

        std::string write_path = output_folder + dataset + "_" + std::to_string(total_new_face_num) + file_type;

        std::ofstream ofs(&write_path[0], std::ofstream::trunc);
        // ofs << std::setprecision(1);
        ofs << "OFF\n";
        ofs << total_new_vertex_num << " " << total_new_face_num << " " << total_new_edge_num << " \n";

        double x_increment = x_length / (x_new_vertex_num - 1);
        double y_increment = y_length / (y_new_vertex_num - 1);

        for (int j = 0; j < y_new_vertex_num; ++j)
        {
            for (int i = 0; i < x_new_vertex_num; ++i)
            {
                double x_new = x_min + x_increment * i;
                double y_new = y_min + y_increment * j;

                double x_new_round = round(x_new);
                double y_new_round = round(y_new);
                // double x_new_round_print = round(x_new * 10) / 10;
                // double y_new_round_print = round(y_new * 10) / 10;
                double z_new_round = 0;

                for (int l = 0; l < mesh.faces().size(); l++)
                {
                    double a_x = mesh.faces()[l].adjacent_vertices()[0]->x();
                    double a_y = mesh.faces()[l].adjacent_vertices()[0]->y();
                    double a_z = mesh.faces()[l].adjacent_vertices()[0]->z();
                    double b_x = mesh.faces()[l].adjacent_vertices()[1]->x();
                    double b_y = mesh.faces()[l].adjacent_vertices()[1]->y();
                    double b_z = mesh.faces()[l].adjacent_vertices()[1]->z();
                    double c_x = mesh.faces()[l].adjacent_vertices()[2]->x();
                    double c_y = mesh.faces()[l].adjacent_vertices()[2]->y();
                    double c_z = mesh.faces()[l].adjacent_vertices()[2]->z();

                    if (point_inside_or_on_edge_tri(x_new_round, y_new_round, a_x, a_y, b_x, b_y, c_x, c_y))
                    {
                        // Barycentric Coordinates https://codeplea.com/triangular-interpolation
                        double w_1 = ((b_y - c_y) * (x_new_round - c_x) + (c_x - b_x) * (y_new_round - c_y)) /
                                     ((b_y - c_y) * (a_x - c_x) + (c_x - b_x) * (a_y - c_y));
                        double w_2 = ((c_y - a_y) * (x_new_round - c_x) + (a_x - c_x) * (y_new_round - c_y)) /
                                     ((b_y - c_y) * (a_x - c_x) + (c_x - b_x) * (a_y - c_y));
                        double w_3 = 1 - w_1 - w_2;

                        if (w_1 < -1e-5 || w_1 > 1 + 1e-5 || w_2 < -1e-5 || w_2 > 1 + 1e-5 || w_3 < -1e-5 || w_3 > 1 + 1e-5)
                        {
                            // print_point_inside_or_on_edge_tri(x_new_round, y_new_round, a_x, a_y, b_x, b_y, c_x, c_y);
                            // std::cout << "a_x: " << a_x << "\ta_y: " << a_y << "\tb_x: " << b_x << "\tb_y: " << b_y << "\tc_x: " << c_x << "\tc_y: " << c_y << "\tx_new_round: " << x_new_round << "\ty_new_round: " << y_new_round << std::endl;
                            // std::cout << "w_1: " << w_1 << "\tw_2: " << w_2 << "\tw_3: " << w_3 << std::endl;
                            // std::cout << "(c_y - a_y): " << (c_y - a_y) << std::endl;
                            // std::cout << "(x_new_round - c_x): " << (x_new_round - c_x) << std::endl;
                            // std::cout << "(a_x - c_x): " << (a_x - c_x) << std::endl;
                            // std::cout << "(y_new_round - c_y): " << (y_new_round - c_y) << std::endl;
                            // std::cout << "(c_y - a_y) * (x_new_round - c_x): " << (c_y - a_y) * (x_new_round - c_x) << std::endl;
                            // std::cout << "(a_x - c_x) * (y_new_round - c_y): " << (a_x - c_x) * (y_new_round - c_y) << std::endl;
                            // std::cout << "((c_y - a_y) * (x_new_round - c_x) + (a_x - c_x) * (y_new_round - c_y)): " << ((c_y - a_y) * (x_new_round - c_x) + (a_x - c_x) * (y_new_round - c_y)) << std::endl;
                            // std::cout << "((b_y - c_y) * (a_x - c_x) + (c_x - b_x) * (a_y - c_y)): " << ((b_y - c_y) * (a_x - c_x) + (c_x - b_x) * (a_y - c_y)) << std::endl;
                            exit(0);
                        }

                        z_new_round = round(w_1 * a_z + w_2 * b_z + w_3 * c_z);
                        // std::cout << z_new_round << std::endl;
                        break;
                    }
                }
                ofs << std::fixed << (int)round(x_new_round) << " " << (int)round(y_new_round) << " " << (int)round(z_new_round) << " \n";
            }
        }

        for (int j = 0; j < y_new_vertex_num - 1; ++j)
        {
            for (int i = 0; i < x_new_vertex_num - 1; ++i)
            {
                int top_left_index = i + j * x_new_vertex_num;
                int top_right_index = (i + 1) + j * x_new_vertex_num;
                int bottom_left_index = i + (j + 1) * x_new_vertex_num;
                int bottom_right_index = (i + 1) + (j + 1) * x_new_vertex_num;

                ofs << "3 " << top_left_index << " " << top_right_index << " " << bottom_left_index << " \n";
                ofs << "3 " << top_right_index << " " << bottom_right_index << " " << bottom_left_index << " \n";
            }
        }

        ifs.close();
        ofs.close();
    }
    return 0;
}
