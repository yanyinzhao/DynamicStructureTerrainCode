# An Efficiently Updatable Path Oracle for Terrain Surfaces
## Overview

This project provides the implementation of the algorithm for calculating an efficiently updatable path oracle on an updated terrain surface. We refer the readers to our paper for more details.

We compared 13 algorithms as follows:

- WSPD-Oracle (oracle based baseline)
- WSPD-Oracle-Adapt (adapted oracle based baseline)
- EAR-Oracle (oracle based baseline)
- EAR-Oracle-Adapt (adapted oracle based baseline)
- UP-Oracle-RanUpdSeq (variation)
- UP-Oracle-FullRad (variation)
- UP-Oracle-NoDistAppr (variation)
- UP-Oracle-NoEffIntChe (variation)
- UP-Oracle-NoEdgPru (variation)
- UP-Oracle-NoEffEdgPru (variation)
- UP-Oracle (our oracle)
- CH-Fly-Algo (on-the-fly baseline)
- K-Fly-Algo (on-the-fly baseline)

Make sure there is a folder called "input/" and a folder called "output/" under the working directory. They will be used for storing the input/output files.

The source code are stored in "src/" folder.

## Dataset

The dataset are stored in "input/" folder.

The datasets are as follows, where XX could be {TJ, SC, GI, AU, LH, VS}:

- "TJpre_1058.off" (multiresolution resolution TJ pre earthquake terrain dataset with dataset size of 1058)
- "TJpost_1058.off" (multiresolution resolution TJ post earthquake terrain dataset with dataset size of 1058)
- "XXpre_500000.off" (default resolution XX pre earthquake terrain dataset with dataset size of 500000)
- "XXpost_500000.off" (default resolution XX post earthquake terrain dataset with dataset size of 500000)
- "XXpre_1002528.off" (multiresolution XX pre earthquake terrain dataset with dataset size of 1002528)
- "XXpost_1002528.off" (multiresolution XX post earthquake terrain dataset with dataset size of 1002528)
- "XXpre_1503378.off" (multiresolution XX pre earthquake terrain dataset with dataset size of 1503378)
- "XXpost_1503378.off" (multiresolution XX post earthquake terrain dataset with dataset size of 1503378)
- "XXpre_2000000.off" (multiresolution XX pre earthquake terrain dataset with dataset size of 2000000)
- "XXpost_2000000.off" (multiresolution XX post earthquake terrain dataset with dataset size of 2000000)
- "XXpre_2504322.off" (multiresolution XX pre earthquake terrain dataset with dataset size of 2504322)
- "XXpost_2504322.off" (multiresolution XX post earthquake terrain dataset with dataset size of 2504322)
- "TJpre_50_poi_on_1058.txt" (POI list with POI number of 50 on "TJpre_1058.off")
- "TJpost_50_poi_on_1058.txt" (POI list with POI number of 50 on "TJpost_1058.off")
- "XXpre_500_poi_on_500000.txt" (POI list with POI number of 500 on "XXpre_500000.off")
- "XXpost_500_poi_on_500000.txt" (POI list with POI number of 500 on "XXpost_500000.off")
- "XXpre_500_poi_on_1002528.txt" (POI list with POI number of 500 on "XXpre_1002528.off")
- "XXpost_500_poi_on_1002528.txt" (POI list with POI number of 500 on "XXpost_1002528.off")
- "XXpre_500_poi_on_1503378.txt" (POI list with POI number of 500 on "XXpre_1503378.off")
- "XXpost_500_poi_on_1503378.txt" (POI list with POI number of 500 on "XXpost_1503378.off")
- "XXpre_500_poi_on_2000000.txt" (POI list with POI number of 500 on "XXpre_2000000.off")
- "XXpost_500_poi_on_2000000.txt" (POI list with POI number of 500 on "XXpost_2000000.off")
- "XXpre_500_poi_on_2504322.txt" (POI list with POI number of 500 on "XXpre_2504322.off")
- "XXpost_500_poi_on_2504322.txt" (POI list with POI number of 500 on "XXpost_2504322.off")
- "XXpre_1000_poi_on_500000.txt" (POI list with POI number of 1000 on "XXpre_500000.off")
- "XXpost_1000_poi_on_500000.txt" (POI list with POI number of 1000 on "XXpost_500000.off")
- "XXpre_1500_poi_on_500000.txt" (POI list with POI number of 1500 on "XXpre_500000.off")
- "XXpost_1500_poi_on_500000.txt" (POI list with POI number of 1500 on "XXpost_500000.off")
- "XXpre_2000_poi_on_500000.txt" (POI list with POI number of 2000 on "XXpre_500000.off")
- "XXpost_2000_poi_on_500000.txt" (POI list with POI number of 2000 on "XXpost_500000.off")
- "XXpre_2500_poi_on_500000.txt" (POI list with POI number of 2500 on "XXpre_500000.off")
- "XXpost_2500_poi_on_500000.txt" (POI list with POI number of 2500 on "XXpost_500000.off")
- "XXpre_50_poi_on_500000.txt" (POI list with POI number of 50 on "XXpre_500000.off")
- "XXpost_50_poi_on_500000.txt" (POI list with POI number of 50 on "XXpost_500000.off")
- "XXpre_50_poi_on_1002528.txt" (POI list with POI number of 50 on "XXpre_1002528.off")
- "XXpost_50_poi_on_1002528.txt" (POI list with POI number of 50 on "XXpost_1002528.off")
- "XXpre_50_poi_on_1503378.txt" (POI list with POI number of 50 on "XXpre_1503378.off")
- "XXpost_50_poi_on_1503378.txt" (POI list with POI number of 50 on "XXpost_1503378.off")
- "XXpre_50_poi_on_2000000.txt" (POI list with POI number of 50 on "XXpre_2000000.off")
- "XXpost_50_poi_on_2000000.txt" (POI list with POI number of 50 on "XXpost_2000000.off")
- "XXpre_50_poi_on_2504322.txt" (POI list with POI number of 50 on "XXpre_2504322.off")
- "XXpost_50_poi_on_2504322.txt" (POI list with POI number of 50 on "XXpost_2504322.off")
- "XXpre_100_poi_on_500000.txt" (POI list with POI number of 100 on "XXpre_500000.off")
- "XXpost_100_poi_on_500000.txt" (POI list with POI number of 100 on "XXpost_500000.off")
- "XXpre_150_poi_on_500000.txt" (POI list with POI number of 150 on "XXpre_500000.off")
- "XXpost_150_poi_on_500000.txt" (POI list with POI number of 150 on "XXpost_500000.off")
- "XXpre_200_poi_on_500000.txt" (POI list with POI number of 200 on "XXpre_500000.off")
- "XXpost_200_poi_on_500000.txt" (POI list with POI number of 200 on "XXpost_500000.off")
- "XXpre_250_poi_on_500000.txt" (POI list with POI number of 250 on "XXpre_500000.off")
- "XXpost_250_poi_on_500000.txt" (POI list with POI number of 250 on "XXpost_500000.off")

Data Format:

For the terrain dataset, we used the .off format in the experiment. The content of the .off file is as follows:

```
OFF

vertices_num faces_num edges_num

1st_vertex_x_coord 1st_vertex_y_coord 1st_vertex_z_coord

2nd_vertex_x_coord 2nd_vertex_y_coord 2nd_vertex_z_coord

......

last_vertex_x_coord last_vertex_y_coord last_vertex_z_coord

1st_face_1st_vertex_ID 1st_face_2nd_vertex_ID 1st_face_3td_vertex_ID

2nd_face_1st_vertex_ID 2nd_face_2nd_vertex_ID 2nd_face_3td_vertex_ID

......

last_face_1st_vertex_ID last_face_2nd_vertex_ID last_face_3td_vertex_ID
```

For the POI list, we used the .txt format in the experiment. The content of the .txt file is as follows:

```
POI_num
1st_POI_index 2nd_POI_index ......
```

## Compile command

```
cd src
g++ -o main main.cpp -std=c++11
```

## Run command

```
./main [terrain_data_and_dataset_size_and_poi_number_map_index] [epsilon]
```

The meaning for each parameter is as follows:

- [terrain_data_and_dataset_size_and_poi_number_map_index]: an index for the map of terrain data and dataset size and poi number (a integer from 0 to 108)
- [epsilon]: the epsilon value (0 < epsilon <= 1)

For the [terrain_data_and_dataset_size_and_poi_number_map_index], each index value corresponding to a terrain data, the dataset size of the terrain and the poi number on the terrain, their relationships are as follows:

| Index | Terrain data | Dataset size | POI number |
| ----------- | ----------- | ----------- | ----------- |
| 0 | TJ | 1058 | 50 |
| 1 | TJ | 500000 | 50 |
| 2 | TJ | 500000 | 100 |
| 3 | TJ | 500000 | 150 |
| 4 | TJ | 500000 | 200 |
| 5 | TJ | 500000 | 250 |
| 6 | TJ | 1002528 | 50 |
| 7 | TJ | 1503378 | 50 |
| 8 | TJ | 2000000 | 50 |
| 9 | TJ | 2504322 | 50 |
| 10 | SC | 500000 | 50 |
| 11 | SC | 500000 | 100 |
| 12 | SC | 500000 | 150 |
| 13 | SC | 500000 | 200 |
| 14 | SC | 500000 | 250 |
| 15 | SC | 1002528 | 50 |
| 16 | SC | 1503378 | 50 |
| 17 | SC | 2000000 | 50 |
| 18 | SC | 2504322 | 50 |
| 19 | GI | 500000 | 50 |
| 20 | GI | 500000 | 100 |
| 21 | GI | 500000 | 150 |
| 22 | GI | 500000 | 200 |
| 23 | GI | 500000 | 250 |
| 24 | GI | 1002528 | 50 |
| 25 | GI | 1503378 | 50 |
| 26 | GI | 2000000 | 50 |
| 27 | GI | 2504322 | 50 |
| 28 | AU | 500000 | 50 |
| 29 | AU | 500000 | 100 |
| 30 | AU | 500000 | 150 |
| 31 | AU | 500000 | 200 |
| 32 | AU | 500000 | 250 |
| 33 | AU | 1002528 | 50 |
| 34 | AU | 1503378 | 50 |
| 35 | AU | 2000000 | 50 |
| 36 | AU | 2504322 | 50 |
| 37 | LH | 500000 | 50 |
| 38 | LH | 500000 | 100 |
| 39 | LH | 500000 | 150 |
| 40 | LH | 500000 | 200 |
| 41 | LH | 500000 | 250 |
| 42 | LH | 1002528 | 50 |
| 43 | LH | 1503378 | 50 |
| 44 | LH | 2000000 | 50 |
| 45 | LH | 2504322 | 50 |
| 46 | VS | 500000 | 50 |
| 47 | VS | 500000 | 100 |
| 48 | VS | 500000 | 150 |
| 49 | VS | 500000 | 200 |
| 50 | VS | 500000 | 250 |
| 51 | VS | 1002528 | 50 |
| 52 | VS | 1503378 | 50 |
| 53 | VS | 2000000 | 50 |
| 54 | VS | 2504322 | 50 |
| 55 | SC | 500000 | 500 |
| 56 | SC | 500000 | 1000 |
| 57 | SC | 500000 | 1500 |
| 58 | SC | 500000 | 2000 |
| 59 | SC | 500000 | 2500 |
| 60 | SC | 1002528 | 500 |
| 61 | SC | 1503378 | 500 |
| 62 | SC | 2000000 | 500 |
| 63 | SC | 2504322 | 500 |
| 64 | TJ | 500000 | 500 |
| 65 | TJ | 500000 | 1000 |
| 66 | TJ | 500000 | 1500 |
| 67 | TJ | 500000 | 2000 |
| 68 | TJ | 500000 | 2500 |
| 69 | TJ | 1002528 | 500 |
| 70 | TJ | 1503378 | 500 |
| 71 | TJ | 2000000 | 500 |
| 72 | TJ | 2504322 | 500 |
| 73 | GI | 500000 | 500 |
| 74 | GI | 500000 | 1000 |
| 75 | GI | 500000 | 1500 |
| 76 | GI | 500000 | 2000 |
| 77 | GI | 500000 | 2500 |
| 78 | GI | 1002528 | 500 |
| 79 | GI | 1503378 | 500 |
| 80 | GI | 2000000 | 500 |
| 81 | GI | 2504322 | 500 |
| 82 | AU | 500000 | 500 |
| 83 | AU | 500000 | 1000 |
| 84 | AU | 500000 | 1500 |
| 85 | AU | 500000 | 2000 |
| 86 | AU | 500000 | 2500 |
| 87 | AU | 1002528 | 500 |
| 88 | AU | 1503378 | 500 |
| 89 | AU | 2000000 | 500 |
| 90 | AU | 2504322 | 500 |
| 91 | LH | 500000 | 500 |
| 92 | LH | 500000 | 1000 |
| 93 | LH | 500000 | 1500 |
| 94 | LH | 500000 | 2000 |
| 95 | LH | 500000 | 2500 |
| 96 | LH | 1002528 | 500 |
| 97 | LH | 1503378 | 500 |
| 98 | LH | 2000000 | 500 |
| 99 | LH | 2504322 | 500 |
| 100 | VS | 500000 | 500 |
| 101 | VS | 500000 | 1000 |
| 102 | VS | 500000 | 1500 |
| 103 | VS | 500000 | 2000 |
| 104 | VS | 500000 | 2500 |
| 105 | VS | 1002528 | 500 |
| 106 | VS | 1503378 | 500 |
| 107 | VS | 2000000 | 500 |
| 108 | VS | 2504322 | 500 |

By default, the project will run WSPD-Oracle, WSPD-Oracle-Adapt, EAR-Oracle, EAR-Oracle-Adapt, UP-Oracle-RanUpdSeq, UP-Oracle-FullRad, UP-Oracle-NoDistAppr, UP-Oracle-NoEffIntChe, UP-Oracle-NoEdgPru, UP-Oracle-NoEffEdgPru, UP-Oracle, CH-Fly-Algo, and K-Fly-Algo. But as mentioned in our paper, WSPD-Oracle, WSPD-Oracle-Adapt, EAR-Oracle, EAR-Oracle-Adapt, UP-Oracle-RanUpdSeq, UP-Oracle-FullRad, and UP-Oracle-NoDistAppr are very time consuming. So when the POI number is large, i.e., [terrain_data_and_dataset_size_and_poi_number_map_index] > 54, the project will only run UP-Oracle-NoEffIntChe,UP-Oracle-NoEdgPru, UP-Oracle-NoEffEdgPru, UP-Oracle, CH-Fly-Algo, and K-Fly-Algo.

An example:

```
./main 0 0.5
```

In this example, [terrain_data_and_dataset_size_and_poi_number_map_index] is 0, [epsilon] is 0.5. So, it will run TJ pre earthquake terrain dataset and TJ post earthquake terrain dataset, with dataset size equal to 1058 and poi number equal to 50, and epsilon is 0.5. It will run 13 algorithms, i.e., WSPD-Oracle, WSPD-Oracle-Adapt, FEAR-Oracle, EAR-Oracle-Adapt, U-Oracle-RanUpdSeq, UP-Oracle-FullRad, UP-Oracle-NoDistAppr, UP-Oracle-NoEffIntChe, UP-Oracle-NoEdgPru, UP-Oracle-NoEffEdgPru, UP-Oracle, CH-Fly-Algo, and K-Fly-Algo.

## Output

The output will be stored in "output/output.txt" file. The format will be as follows:

```

[pre_dataset] [post_dataset] [datasize] [poi_num] [epsilon] [pre_construction_time (ms)] [pre_memory_usage (MB)]  [pre_index_size (MB)] [pre_index_edge_num] [pre_index/MST_weight] [post_update_time1 (ms)] [post_update_time2 (ms)] [post_query_time (ms)] [post_memory_usage (MB)] [post_index_size (MB)] [post_output_size (MB)] [post_index_edge_num] [post_index/MST_weight] [post_distance_error]

```

These information will also be shown in the terminal. 

[post_update_time1] means the update time for terrain surface and POIs change detection step and pairwise P2P exact shortest path updating time, and [post_update_time2] means the update time for sub-graph generating step. 