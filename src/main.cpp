#include "algorithms.h"

struct input_struct
{
       std::string pre_input;
       std::string post_input;
       std::string pre_poi;
       std::string post_poi;
       int source_poi_index;
       int destination_poi_index;

       input_struct() {}

       input_struct(std::string _pre_input, std::string _post_input,
                    std::string _pre_poi, std::string _post_poi,
                    int _source_poi_index, int _destination_poi_index)
       {
              pre_input = _pre_input;
              post_input = _post_input;
              pre_poi = _pre_poi;
              post_poi = _post_poi;
              source_poi_index = _source_poi_index;
              destination_poi_index = _destination_poi_index;
       }
};

int main(int argc, char **argv)
{
       int input_file_index = std::stoi(argv[1]);
       double epsilon = std::stod(argv[2]);

       std::string input_folder = "../input/";

       std::map<int, input_struct> input_file{
           {0, input_struct("SCpre_10082.off", "SCpost_10082.off", "SCpre_50_poi_on_10082.txt", "SCpost_50_poi_on_10082.txt", 0, 1)},
           {1, input_struct("SCpre_10082.off", "SCpost_10082.off", "SCpre_100_poi_on_10082.txt", "SCpost_100_poi_on_10082.txt", 0, 1)},
           {2, input_struct("SCpre_10082.off", "SCpost_10082.off", "SCpre_150_poi_on_10082.txt", "SCpost_150_poi_on_10082.txt", 0, 1)},
           {3, input_struct("SCpre_10082.off", "SCpost_10082.off", "SCpre_200_poi_on_10082.txt", "SCpost_200_poi_on_10082.txt", 0, 1)},
           {4, input_struct("SCpre_10082.off", "SCpost_10082.off", "SCpre_250_poi_on_10082.txt", "SCpost_250_poi_on_10082.txt", 0, 1)},
           {5, input_struct("SCpre_20000.off", "SCpost_20000.off", "SCpre_50_poi_on_20000.txt", "SCpost_50_poi_on_20000.txt", 0, 1)},
           {6, input_struct("SCpre_30258.off", "SCpost_30258.off", "SCpre_50_poi_on_30258.txt", "SCpost_50_poi_on_30258.txt", 0, 1)},
           {7, input_struct("SCpre_40328.off", "SCpost_40328.off", "SCpre_50_poi_on_40328.txt", "SCpost_50_poi_on_40328.txt", 0, 1)},
           {8, input_struct("SCpre_50562.off", "SCpost_50562.off", "SCpre_50_poi_on_50562.txt", "SCpost_50_poi_on_50562.txt", 0, 1)},
           {9, input_struct("AUpre_10082.off", "AUpost_10082.off", "AUpre_50_poi_on_10082.txt", "AUpost_50_poi_on_10082.txt", 0, 1)},
           {10, input_struct("AUpre_10082.off", "AUpost_10082.off", "AUpre_100_poi_on_10082.txt", "AUpost_100_poi_on_10082.txt", 0, 1)},
           {11, input_struct("AUpre_10082.off", "AUpost_10082.off", "AUpre_150_poi_on_10082.txt", "AUpost_150_poi_on_10082.txt", 0, 1)},
           {12, input_struct("AUpre_10082.off", "AUpost_10082.off", "AUpre_200_poi_on_10082.txt", "AUpost_200_poi_on_10082.txt", 0, 1)},
           {13, input_struct("AUpre_10082.off", "AUpost_10082.off", "AUpre_250_poi_on_10082.txt", "AUpost_250_poi_on_10082.txt", 0, 1)},
           {14, input_struct("VSpre_10082.off", "VSpost_10082.off", "VSpre_50_poi_on_10082.txt", "VSpost_50_poi_on_10082.txt", 0, 1)},
           {15, input_struct("VSpre_10082.off", "VSpost_10082.off", "VSpre_100_poi_on_10082.txt", "VSpost_100_poi_on_10082.txt", 0, 1)},
           {16, input_struct("VSpre_10082.off", "VSpost_10082.off", "VSpre_150_poi_on_10082.txt", "VSpost_150_poi_on_10082.txt", 0, 1)},
           {17, input_struct("VSpre_10082.off", "VSpost_10082.off", "VSpre_200_poi_on_10082.txt", "VSpost_200_poi_on_10082.txt", 0, 1)},
           {18, input_struct("VSpre_10082.off", "VSpost_10082.off", "VSpre_250_poi_on_10082.txt", "VSpost_250_poi_on_10082.txt", 0, 1)},

           {19, input_struct("SCpre_500000.off", "SCpost_500000.off", "SCpre_500_poi_on_500000.txt", "SCpost_500_poi_on_500000.txt", 0, 1)},
           {20, input_struct("SCpre_500000.off", "SCpost_500000.off", "SCpre_1000_poi_on_500000.txt", "SCpost_1000_poi_on_500000.txt", 0, 1)},
           {21, input_struct("SCpre_500000.off", "SCpost_500000.off", "SCpre_1500_poi_on_500000.txt", "SCpost_1500_poi_on_500000.txt", 0, 1)},
           {22, input_struct("SCpre_500000.off", "SCpost_500000.off", "SCpre_2000_poi_on_500000.txt", "SCpost_2000_poi_on_500000.txt", 0, 1)},
           {23, input_struct("SCpre_500000.off", "SCpost_500000.off", "SCpre_2500_poi_on_500000.txt", "SCpost_2500_poi_on_500000.txt", 0, 1)},
           {24, input_struct("SCpre_1002528.off", "SCpost_1002528.off", "SCpre_500_poi_on_1002528.txt", "SCpost_500_poi_on_1002528.txt", 0, 1)},
           {25, input_struct("SCpre_1503378.off", "SCpost_1503378.off", "SCpre_500_poi_on_1503378.txt", "SCpost_500_poi_on_1503378.txt", 0, 1)},
           {26, input_struct("SCpre_2000000.off", "SCpost_2000000.off", "SCpre_500_poi_on_2000000.txt", "SCpost_500_poi_on_2000000.txt", 0, 1)},
           {27, input_struct("SCpre_2504322.off", "SCpost_2504322.off", "SCpre_500_poi_on_2504322.txt", "SCpost_500_poi_on_2504322.txt", 0, 1)},
           {28, input_struct("AUpre_500000.off", "AUpost_500000.off", "AUpre_500_poi_on_500000.txt", "AUpost_500_poi_on_500000.txt", 0, 1)},
           {29, input_struct("AUpre_500000.off", "AUpost_500000.off", "AUpre_1000_poi_on_500000.txt", "AUpost_1000_poi_on_500000.txt", 0, 1)},
           {30, input_struct("AUpre_500000.off", "AUpost_500000.off", "AUpre_1500_poi_on_500000.txt", "AUpost_1500_poi_on_500000.txt", 0, 1)},
           {31, input_struct("AUpre_500000.off", "AUpost_500000.off", "AUpre_2000_poi_on_500000.txt", "AUpost_2000_poi_on_500000.txt", 0, 1)},
           {32, input_struct("AUpre_500000.off", "AUpost_500000.off", "AUpre_2500_poi_on_500000.txt", "AUpost_2500_poi_on_500000.txt", 0, 1)},
           {33, input_struct("AUpre_1002528.off", "AUpost_1002528.off", "AUpre_500_poi_on_1002528.txt", "AUpost_500_poi_on_1002528.txt", 0, 1)},
           {34, input_struct("AUpre_1503378.off", "AUpost_1503378.off", "AUpre_500_poi_on_1503378.txt", "AUpost_500_poi_on_1503378.txt", 0, 1)},
           {35, input_struct("AUpre_2000000.off", "AUpost_2000000.off", "AUpre_500_poi_on_2000000.txt", "AUpost_500_poi_on_2000000.txt", 0, 1)},
           {36, input_struct("AUpre_2504322.off", "AUpost_2504322.off", "AUpre_500_poi_on_2504322.txt", "AUpost_500_poi_on_2504322.txt", 0, 1)},
           {37, input_struct("VSpre_500000.off", "VSpost_500000.off", "VSpre_500_poi_on_500000.txt", "VSpost_500_poi_on_500000.txt", 0, 1)},
           {38, input_struct("VSpre_500000.off", "VSpost_500000.off", "VSpre_1000_poi_on_500000.txt", "VSpost_1000_poi_on_500000.txt", 0, 1)},
           {39, input_struct("VSpre_500000.off", "VSpost_500000.off", "VSpre_1500_poi_on_500000.txt", "VSpost_1500_poi_on_500000.txt", 0, 1)},
           {40, input_struct("VSpre_500000.off", "VSpost_500000.off", "VSpre_2000_poi_on_500000.txt", "VSpost_2000_poi_on_500000.txt", 0, 1)},
           {41, input_struct("VSpre_500000.off", "VSpost_500000.off", "VSpre_2500_poi_on_500000.txt", "VSpost_2500_poi_on_500000.txt", 0, 1)},
           {42, input_struct("VSpre_1002528.off", "VSpost_1002528.off", "VSpre_500_poi_on_1002528.txt", "VSpost_500_poi_on_1002528.txt", 0, 1)},
           {43, input_struct("VSpre_1503378.off", "VSpost_1503378.off", "VSpre_500_poi_on_1503378.txt", "VSpost_500_poi_on_1503378.txt", 0, 1)},
           {44, input_struct("VSpre_2000000.off", "VSpost_2000000.off", "VSpre_500_poi_on_2000000.txt", "VSpost_500_poi_on_2000000.txt", 0, 1)},
           {45, input_struct("VSpre_2504322.off", "VSpost_2504322.off", "VSpre_500_poi_on_2504322.txt", "VSpost_500_poi_on_2504322.txt", 0, 1)},
       };

       std::string pre_input_file = input_folder + input_file[input_file_index].pre_input;
       std::string post_input_file = input_folder + input_file[input_file_index].post_input;
       std::string pre_input_poi = input_folder + input_file[input_file_index].pre_poi;
       std::string post_input_poi = input_folder + input_file[input_file_index].post_poi;
       int source_poi_index = input_file[input_file_index].source_poi_index;
       int destination_poi_index = input_file[input_file_index].destination_poi_index;
       assert(source_poi_index != destination_poi_index);

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
       assert(epsilon > 0 && epsilon <= 1);

       std::string write_file_header = input_file[input_file_index].pre_input + "\t" +
                                       input_file[input_file_index].post_input + "\t" +
                                       std::to_string(pre_mesh.faces().size()) + "\t" +
                                       std::to_string(poi_num) + "\t" +
                                       std::to_string(epsilon);

       std::cout << "pre_dataset: " << input_file[input_file_index].pre_input << "\tpost_dataset: " << input_file[input_file_index].post_input << "\tdatasize: " << pre_mesh.faces().size() << "\tpoi_num: " << poi_num << "\tepsilon: " << epsilon << std::endl;
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

       std::cout << "== UE_N1 ==" << std::endl;
       UE_N1(poi_num, &pre_mesh, pre_poi_list, &post_mesh, post_poi_list, source_poi_index,
             destination_poi_index, pre_exact_distance, post_exact_distance, pre_MST_weight,
             post_MST_weight, write_file_header);
       std::cout << std::endl;

       std::cout << "== UE_N2 ==" << std::endl;
       UE_N2(poi_num, &pre_mesh, pre_poi_list, &post_mesh, post_poi_list, epsilon,
             source_poi_index, destination_poi_index, pre_exact_distance, post_exact_distance,
             pre_MST_weight, post_MST_weight, write_file_header);
       std::cout << std::endl;

       std::cout << "== UE ==" << std::endl;
       UE(poi_num, &pre_mesh, pre_poi_list, &post_mesh, post_poi_list, epsilon,
          source_poi_index, destination_poi_index, pre_exact_distance, post_exact_distance,
          pre_MST_weight, post_MST_weight, write_file_header);
       std::cout << std::endl;

       if (input_file_index >= 0 && input_file_index <= 18)
       {
              std::cout << "== WSPD_oracle ==" << std::endl;
              WSPD_oracle(poi_num, &pre_mesh, pre_poi_list, &post_mesh, post_poi_list, epsilon,
                          source_poi_index, destination_poi_index, pre_exact_distance, post_exact_distance,
                          pre_MST_weight, post_MST_weight, write_file_header);
              std::cout << std::endl;

              std::cout << "== SP_oracle ==" << std::endl;
              SP_oracle(&pre_mesh, pre_poi_list, &post_mesh, post_poi_list, epsilon,
                        source_poi_index, destination_poi_index, pre_exact_distance, post_exact_distance,
                        write_file_header);
              std::cout << std::endl;
       }

       std::cout << "== KF ==" << std::endl;
       KF(&pre_mesh, pre_poi_list, &post_mesh, post_poi_list, epsilon,
          source_poi_index, destination_poi_index, pre_exact_distance, post_exact_distance,
          write_file_header);
       std::cout << std::endl;
}
