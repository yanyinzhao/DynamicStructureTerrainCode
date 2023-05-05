# Fast Update Path Oracle on Updated Terrain Surface

## Overview

This project provides the implementation of the algorithm for calculating a fast update path oracle on updated terrain surface. We refer the readers to our paper for more details.

We compared 10 algorithms as follows:

- WSPD-Oracle (oracle based baseline)
- WSPD-Oracle-Adapt (adapted oracle based baseline)
- FU-Oracle-RanUpdSeq (variation)
- FU-Oracle-FullRad (variation)
- FU-Oracle-NoDistAppr (variation)
- FU-Oracle-NoEffIntChe (variation)
- FU-Oracle-NoEdgPru (variation)
- FU-Oracle-NoEffEdgPru (variation)
- FU-Oracle (our oracle)
- K-Fly-Algo (on-the-fly baseline)

Make sure there is a folder called "input/" and a folder called "output/" under the working directory. They will be used for storing the input/output files.

The source code are stored in "src/" folder.

## Dataset

The dataset are stored in "input/" folder.

The datasets are as follows, where XX could be {TJ, SC, GI, AU, LH, VS}:

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

- [terrain_data_and_dataset_size_and_poi_number_map_index]: an index for the map of terrain data and dataset size and poi number (a integer from 0 to 107)
- [epsilon]: the epsilon value (0 < epsilon <= 1)

For the [terrain_data_and_dataset_size_and_poi_number_map_index], each index value corresponding to a terrain data, the dataset size of the terrain and the poi number on the terrain, their relationships are as follows:

| Index | Terrain data | Dataset size | POI number |
| ----------- | ----------- | ----------- | ----------- |
| 0 | TJ | 500000 | 50 |
| 1 | TJ | 500000 | 100 |
| 2 | TJ | 500000 | 150 |
| 3 | TJ | 500000 | 200 |
| 4 | TJ | 500000 | 250 |
| 5 | TJ | 1002528 | 50 |
| 6 | TJ | 1503378 | 50 |
| 7 | TJ | 2000000 | 50 |
| 8 | TJ | 2504322 | 50 |
| 9 | SC | 500000 | 50 |
| 10 | SC | 500000 | 100 |
| 11 | SC | 500000 | 150 |
| 12 | SC | 500000 | 200 |
| 13 | SC | 500000 | 250 |
| 14 | SC | 1002528 | 50 |
| 15 | SC | 1503378 | 50 |
| 16 | SC | 2000000 | 50 |
| 17 | SC | 2504322 | 50 |
| 18 | GI | 500000 | 50 |
| 19 | GI | 500000 | 100 |
| 20 | GI | 500000 | 150 |
| 21 | GI | 500000 | 200 |
| 22 | GI | 500000 | 250 |
| 23 | GI | 1002528 | 50 |
| 24 | GI | 1503378 | 50 |
| 25 | GI | 2000000 | 50 |
| 26 | GI | 2504322 | 50 |
| 27 | AU | 500000 | 50 |
| 28 | AU | 500000 | 100 |
| 29 | AU | 500000 | 150 |
| 30 | AU | 500000 | 200 |
| 31 | AU | 500000 | 250 |
| 32 | AU | 1002528 | 50 |
| 33 | AU | 1503378 | 50 |
| 34 | AU | 2000000 | 50 |
| 35 | AU | 2504322 | 50 |
| 36 | LH | 500000 | 50 |
| 37 | LH | 500000 | 100 |
| 38 | LH | 500000 | 150 |
| 39 | LH | 500000 | 200 |
| 40 | LH | 500000 | 250 |
| 41 | LH | 1002528 | 50 |
| 42 | LH | 1503378 | 50 |
| 43 | LH | 2000000 | 50 |
| 44 | LH | 2504322 | 50 |
| 45 | VS | 500000 | 50 |
| 46 | VS | 500000 | 100 |
| 47 | VS | 500000 | 150 |
| 48 | VS | 500000 | 200 |
| 49 | VS | 500000 | 250 |
| 50 | VS | 1002528 | 50 |
| 51 | VS | 1503378 | 50 |
| 52 | VS | 2000000 | 50 |
| 53 | VS | 2504322 | 50 |
| 54 | SC | 500000 | 500 |
| 55 | SC | 500000 | 1000 |
| 56 | SC | 500000 | 1500 |
| 57 | SC | 500000 | 2000 |
| 58 | SC | 500000 | 2500 |
| 59 | SC | 1002528 | 500 |
| 60 | SC | 1503378 | 500 |
| 61 | SC | 2000000 | 500 |
| 62 | SC | 2504322 | 500 |
| 63 | TJ | 500000 | 500 |
| 64 | TJ | 500000 | 1000 |
| 65 | TJ | 500000 | 1500 |
| 66 | TJ | 500000 | 2000 |
| 67 | TJ | 500000 | 2500 |
| 68 | TJ | 1002528 | 500 |
| 69 | TJ | 1503378 | 500 |
| 70 | TJ | 2000000 | 500 |
| 71 | TJ | 2504322 | 500 |
| 72 | GI | 500000 | 500 |
| 73 | GI | 500000 | 1000 |
| 74 | GI | 500000 | 1500 |
| 75 | GI | 500000 | 2000 |
| 76 | GI | 500000 | 2500 |
| 77 | GI | 1002528 | 500 |
| 78 | GI | 1503378 | 500 |
| 79 | GI | 2000000 | 500 |
| 80 | GI | 2504322 | 500 |
| 81 | AU | 500000 | 500 |
| 82 | AU | 500000 | 1000 |
| 83 | AU | 500000 | 1500 |
| 84 | AU | 500000 | 2000 |
| 85 | AU | 500000 | 2500 |
| 86 | AU | 1002528 | 500 |
| 87 | AU | 1503378 | 500 |
| 88 | AU | 2000000 | 500 |
| 89 | AU | 2504322 | 500 |
| 90 | LH | 500000 | 500 |
| 91 | LH | 500000 | 1000 |
| 92 | LH | 500000 | 1500 |
| 93 | LH | 500000 | 2000 |
| 94 | LH | 500000 | 2500 |
| 95 | LH | 1002528 | 500 |
| 96 | LH | 1503378 | 500 |
| 97 | LH | 2000000 | 500 |
| 98 | LH | 2504322 | 500 |
| 99 | VS | 500000 | 500 |
| 100 | VS | 500000 | 1000 |
| 101 | VS | 500000 | 1500 |
| 102 | VS | 500000 | 2000 |
| 103 | VS | 500000 | 2500 |
| 104 | VS | 1002528 | 500 |
| 105 | VS | 1503378 | 500 |
| 106 | VS | 2000000 | 500 |
| 107 | VS | 2504322 | 500 |

By default, the project will run WSPD-Oracle, WSPD-Oracle-Adapt, FU-Oracle-RanUpdSeq, FU-Oracle-FullRad, FU-Oracle-NoDistAppr, FU-Oracle-NoEffIntChe, FU-Oracle-NoEdgPru, FU-Oracle-NoEffEdgPru, FU-Oracle, and K-Fly-Algo. But as mentioned in our paper, WSPD-Oracle and WSPD-Oracle-Adapt, FU-Oracle-RanUpdSeq, FU-Oracle-FullRad, FU-Oracle-NoDistAppr are very time consuming. So when the POI number is large, i.e., [terrain_data_and_dataset_size_and_poi_number_map_index] > 53, the project will only run FU-Oracle-NoEffIntChe,FU-Oracle-NoEdgPru, FU-Oracle-NoEffEdgPru, FU-Oracle, and K-Fly-Algo.

An example:

```
./main 0 0.5
```

In this example, [terrain_data_and_dataset_size_and_poi_number_map_index] is 0, [epsilon] is 0.5. So, it will run TJ pre earthquake terrain dataset and TJ post earthquake terrain dataset, with dataset size equal to 500000 and poi number equal to 50, and epsilon is 0.5. It will run 10 algorithms, i.e., WSPD-Oracle, WSPD-Oracle-Adapt, FU-Oracle-RanUpdSeq, FU-Oracle-FullRad, FU-Oracle-NoDistAppr, FU-Oracle-NoEffIntChe, FU-Oracle-NoEdgPru, FU-Oracle-NoEffEdgPru, FU-Oracle, and K-Fly-Algo.

## Output

The output will be stored in "output/output.txt" file. The format will be as follows:

```

[pre_dataset] [post_dataset] [datasize] [poi_num] [epsilon] [pre_construction_time (ms)] [pre_memory_usage (MB)]  [pre_index_size (MB)] [pre_index_edge_num] [pre_index/MST_weight] [post_update_time1 (ms)] [post_update_time2 (ms)] [post_query_time (ms)] [post_memory_usage (MB)] [post_index_size (MB)] [post_output_size (MB)] [post_index_edge_num] [post_index/MST_weight] [post_distance_error]

```

These information will also be shown in the terminal. 

[post_update_time1] means the update time for terrain surface and POIs change detection step and pairwise P2P exact shortest path updating time, and [post_update_time2] means the update time for sub-graph generating step. 