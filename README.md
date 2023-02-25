# Path Oracle on Updated Terrain Surface

## Overview

This project provides the implementation of the algorithm for calculating a shortest path oracle on updated terrain surface.

Our oralce FU-Oracle, and the baselines, i.e., WSPD-Oracle, WSPD-Oracle-Adapt, and K-Algo are studied in the experiments. In order to conduct the ablation study, i.e., show that algorithm HGS could significantly reduce the running time compared with algorithm GS, and could significantly reduce the oracle size and oracle weight compared with original complete graph, we also studied FU-Oracle-Naive1 and FU-Oracle-Naive2 in the experiments. In total, we compared six algorithms, namely, WSPD-Oracle, WSPD-Oracle-Adapt, FU-Oracle-Naive1, FU-Oracle-Naive2, FU-Oracle, and K-Algo. Since WSPD-Oracle and WSPD-Oracle-Adapt are not feasible on large-version POI due to their expensive oracle construction time, so we (1) compared these six algorithms on our 15 datasets with small-version POI (default 50 POIs), and (2) compared FU-Oracle-Naive1, FU-Oracle-Naive2, FU-Oracle, and K-Algo on our 15 datasets with large-version POI (default 500 POIs). We refer the readers to our paper for more details.

In total, we compared six algorithms as follows:

- WSPD-Oracle (oracle based baseline)
- WSPD-Oracle-Adapt (adapted oracle based baseline)
- FU-Oracle-Naive1 (variation)
- FU-Oracle-Naive2 (variation)
- FU-Oracle (our oracle)
- K-Algo (on-the-fly baseline)

Make sure there is a folder called "input/" and a folder called "output/" under the working directory. They will be used for storing the input/output files.

The source code are stored in "src/" folder.

## Dataset

The dataset are stored in "input/" folder.

The datasets are as follows:

- "SCpre_500000.off" (default resolution SC pre earthquake terrain dataset with dataset size of 500000)
- "SCpost_500000.off" (default resolution SC post earthquake terrain dataset with dataset size of 500000)
- "SCpre_1002528.off" (multiresolution SC pre earthquake terrain dataset with dataset size of 1002528)
- "SCpost_1002528.off" (multiresolution SC post earthquake terrain dataset with dataset size of 1002528)
- "SCpre_1503378.off" (multiresolution SC pre earthquake terrain dataset with dataset size of 1503378)
- "SCpost_1503378.off" (multiresolution SC post earthquake terrain dataset with dataset size of 1503378)
- "SCpre_2000000.off" (multiresolution SC pre earthquake terrain dataset with dataset size of 2000000)
- "SCpost_2000000.off" (multiresolution SC post earthquake terrain dataset with dataset size of 2000000)
- "SCpre_2504322.off" (multiresolution SC pre earthquake terrain dataset with dataset size of 2504322)
- "SCpost_2504322.off" (multiresolution SC post earthquake terrain dataset with dataset size of 2504322)
- "AUpre_500000.off" (default resolution AU pre earthquake terrain dataset with dataset size of 500000)
- "AUpost_500000.off" (default resolution AU post earthquake terrain dataset with dataset size of 500000)
- "AUpre_1002528.off" (multiresolution AU pre earthquake terrain dataset with dataset size of 1002528)
- "AUpost_1002528.off" (multiresolution AU post earthquake terrain dataset with dataset size of 1002528)
- "AUpre_1503378.off" (multiresolution AU pre earthquake terrain dataset with dataset size of 1503378)
- "AUpost_1503378.off" (multiresolution AU post earthquake terrain dataset with dataset size of 1503378)
- "AUpre_2000000.off" (multiresolution AU pre earthquake terrain dataset with dataset size of 2000000)
- "AUpost_2000000.off" (multiresolution AU post earthquake terrain dataset with dataset size of 2000000)
- "AUpre_2504322.off" (multiresolution AU pre earthquake terrain dataset with dataset size of 2504322)
- "AUpost_2504322.off" (multiresolution AU post earthquake terrain dataset with dataset size of 2504322)
- "VSpre_500000.off" (default resolution VS pre earthquake terrain dataset with dataset size of 500000)
- "VSpost_500000.off" (default resolution VS post earthquake terrain dataset with dataset size of 500000)
- "VSpre_1002528.off" (multiresolution VS pre earthquake terrain dataset with dataset size of 1002528)
- "VSpost_1002528.off" (multiresolution VS post earthquake terrain dataset with dataset size of 1002528)
- "VSpre_1503378.off" (multiresolution VS pre earthquake terrain dataset with dataset size of 1503378)
- "VSpost_1503378.off" (multiresolution VS post earthquake terrain dataset with dataset size of 1503378)
- "VSpre_2000000.off" (multiresolution VS pre earthquake terrain dataset with dataset size of 2000000)
- "VSpost_2000000.off" (multiresolution VS post earthquake terrain dataset with dataset size of 2000000)
- "VSpre_2504322.off" (multiresolution VS pre earthquake terrain dataset with dataset size of 2504322)
- "VSpost_2504322.off" (multiresolution VS post earthquake terrain dataset with dataset size of 2504322)
- "SCpre_500_poi_on_500000.txt" (POI list with POI number of 500 on "SCpre_500000.off")
- "SCpost_500_poi_on_500000.txt" (POI list with POI number of 500 on "SCpost_500000.off")
- "SCpre_500_poi_on_1002528.txt" (POI list with POI number of 500 on "SCpre_1002528.off")
- "SCpost_500_poi_on_1002528.txt" (POI list with POI number of 500 on "SCpost_1002528.off")
- "SCpre_500_poi_on_1503378.txt" (POI list with POI number of 500 on "SCpre_1503378.off")
- "SCpost_500_poi_on_1503378.txt" (POI list with POI number of 500 on "SCpost_1503378.off")
- "SCpre_500_poi_on_2000000.txt" (POI list with POI number of 500 on "SCpre_2000000.off")
- "SCpost_500_poi_on_2000000.txt" (POI list with POI number of 500 on "SCpost_2000000.off")
- "SCpre_500_poi_on_2504322.txt" (POI list with POI number of 500 on "SCpre_2504322.off")
- "SCpost_500_poi_on_2504322.txt" (POI list with POI number of 500 on "SCpost_2504322.off")
- "SCpre_1000_poi_on_500000.txt" (POI list with POI number of 1000 on "SCpre_500000.off")
- "SCpost_1000_poi_on_500000.txt" (POI list with POI number of 1000 on "SCpost_500000.off")
- "SCpre_1500_poi_on_500000.txt" (POI list with POI number of 1500 on "SCpre_500000.off")
- "SCpost_1500_poi_on_500000.txt" (POI list with POI number of 1500 on "SCpost_500000.off")
- "SCpre_2000_poi_on_500000.txt" (POI list with POI number of 2000 on "SCpre_500000.off")
- "SCpost_2000_poi_on_500000.txt" (POI list with POI number of 2000 on "SCpost_500000.off")
- "SCpre_2500_poi_on_500000.txt" (POI list with POI number of 2500 on "SCpre_500000.off")
- "SCpost_2500_poi_on_500000.txt" (POI list with POI number of 2500 on "SCpost_500000.off")
- "SCpre_50_poi_on_500000.txt" (POI list with POI number of 50 on "SCpre_500000.off")
- "SCpost_50_poi_on_500000.txt" (POI list with POI number of 50 on "SCpost_500000.off")
- "SCpre_50_poi_on_1002528.txt" (POI list with POI number of 50 on "SCpre_1002528.off")
- "SCpost_50_poi_on_1002528.txt" (POI list with POI number of 50 on "SCpost_1002528.off")
- "SCpre_50_poi_on_1503378.txt" (POI list with POI number of 50 on "SCpre_1503378.off")
- "SCpost_50_poi_on_1503378.txt" (POI list with POI number of 50 on "SCpost_1503378.off")
- "SCpre_50_poi_on_2000000.txt" (POI list with POI number of 50 on "SCpre_2000000.off")
- "SCpost_50_poi_on_2000000.txt" (POI list with POI number of 50 on "SCpost_2000000.off")
- "SCpre_50_poi_on_2504322.txt" (POI list with POI number of 50 on "SCpre_2504322.off")
- "SCpost_50_poi_on_2504322.txt" (POI list with POI number of 50 on "SCpost_2504322.off")
- "SCpre_100_poi_on_500000.txt" (POI list with POI number of 100 on "SCpre_500000.off")
- "SCpost_100_poi_on_500000.txt" (POI list with POI number of 100 on "SCpost_500000.off")
- "SCpre_150_poi_on_500000.txt" (POI list with POI number of 150 on "SCpre_500000.off")
- "SCpost_150_poi_on_500000.txt" (POI list with POI number of 150 on "SCpost_500000.off")
- "SCpre_200_poi_on_500000.txt" (POI list with POI number of 200 on "SCpre_500000.off")
- "SCpost_200_poi_on_500000.txt" (POI list with POI number of 200 on "SCpost_500000.off")
- "SCpre_250_poi_on_500000.txt" (POI list with POI number of 250 on "SCpre_500000.off")
- "SCpost_250_poi_on_500000.txt" (POI list with POI number of 250 on "SCpost_500000.off")
- "AUpre_500_poi_on_500000.txt" (POI list with POI number of 500 on "AUpre_500000.off")
- "AUpost_500_poi_on_500000.txt" (POI list with POI number of 500 on "AUpost_500000.off")
- "AUpre_500_poi_on_1002528.txt" (POI list with POI number of 500 on "AUpre_1002528.off")
- "AUpost_500_poi_on_1002528.txt" (POI list with POI number of 500 on "AUpost_1002528.off")
- "AUpre_500_poi_on_1503378.txt" (POI list with POI number of 500 on "AUpre_1503378.off")
- "AUpost_500_poi_on_1503378.txt" (POI list with POI number of 500 on "AUpost_1503378.off")
- "AUpre_500_poi_on_2000000.txt" (POI list with POI number of 500 on "AUpre_2000000.off")
- "AUpost_500_poi_on_2000000.txt" (POI list with POI number of 500 on "AUpost_2000000.off")
- "AUpre_500_poi_on_2504322.txt" (POI list with POI number of 500 on "AUpre_2504322.off")
- "AUpost_500_poi_on_2504322.txt" (POI list with POI number of 500 on "AUpost_2504322.off")
- "AUpre_1000_poi_on_500000.txt" (POI list with POI number of 1000 on "AUpre_500000.off")
- "AUpost_1000_poi_on_500000.txt" (POI list with POI number of 1000 on "AUpost_500000.off")
- "AUpre_1500_poi_on_500000.txt" (POI list with POI number of 1500 on "AUpre_500000.off")
- "AUpost_1500_poi_on_500000.txt" (POI list with POI number of 1500 on "AUpost_500000.off")
- "AUpre_2000_poi_on_500000.txt" (POI list with POI number of 2000 on "AUpre_500000.off")
- "AUpost_2000_poi_on_500000.txt" (POI list with POI number of 2000 on "AUpost_500000.off")
- "AUpre_2500_poi_on_500000.txt" (POI list with POI number of 2500 on "AUpre_500000.off")
- "AUpost_2500_poi_on_500000.txt" (POI list with POI number of 2500 on "AUpost_500000.off")
- "AUpre_50_poi_on_500000.txt" (POI list with POI number of 50 on "AUpre_500000.off")
- "AUpost_50_poi_on_500000.txt" (POI list with POI number of 50 on "AUpost_500000.off")
- "AUpre_50_poi_on_1002528.txt" (POI list with POI number of 50 on "AUpre_1002528.off")
- "AUpost_50_poi_on_1002528.txt" (POI list with POI number of 50 on "AUpost_1002528.off")
- "AUpre_50_poi_on_1503378.txt" (POI list with POI number of 50 on "AUpre_1503378.off")
- "AUpost_50_poi_on_1503378.txt" (POI list with POI number of 50 on "AUpost_1503378.off")
- "AUpre_50_poi_on_2000000.txt" (POI list with POI number of 50 on "AUpre_2000000.off")
- "AUpost_50_poi_on_2000000.txt" (POI list with POI number of 50 on "AUpost_2000000.off")
- "AUpre_50_poi_on_2504322.txt" (POI list with POI number of 50 on "AUpre_2504322.off")
- "AUpost_50_poi_on_2504322.txt" (POI list with POI number of 50 on "AUpost_2504322.off")
- "AUpre_100_poi_on_500000.txt" (POI list with POI number of 100 on "AUpre_500000.off")
- "AUpost_100_poi_on_500000.txt" (POI list with POI number of 100 on "AUpost_500000.off")
- "AUpre_150_poi_on_500000.txt" (POI list with POI number of 150 on "AUpre_500000.off")
- "AUpost_150_poi_on_500000.txt" (POI list with POI number of 150 on "AUpost_500000.off")
- "AUpre_200_poi_on_500000.txt" (POI list with POI number of 200 on "AUpre_500000.off")
- "AUpost_200_poi_on_500000.txt" (POI list with POI number of 200 on "AUpost_500000.off")
- "AUpre_250_poi_on_500000.txt" (POI list with POI number of 250 on "AUpre_500000.off")
- "AUpost_250_poi_on_500000.txt" (POI list with POI number of 250 on "AUpost_500000.off")
- "VSpre_500_poi_on_500000.txt" (POI list with POI number of 500 on "VSpre_500000.off")
- "VSpost_500_poi_on_500000.txt" (POI list with POI number of 500 on "VSpost_500000.off")
- "VSpre_500_poi_on_1002528.txt" (POI list with POI number of 500 on "VSpre_1002528.off")
- "VSpost_500_poi_on_1002528.txt" (POI list with POI number of 500 on "VSpost_1002528.off")
- "VSpre_500_poi_on_1503378.txt" (POI list with POI number of 500 on "VSpre_1503378.off")
- "VSpost_500_poi_on_1503378.txt" (POI list with POI number of 500 on "VSpost_1503378.off")
- "VSpre_500_poi_on_2000000.txt" (POI list with POI number of 500 on "VSpre_2000000.off")
- "VSpost_500_poi_on_2000000.txt" (POI list with POI number of 500 on "VSpost_2000000.off")
- "VSpre_500_poi_on_2504322.txt" (POI list with POI number of 500 on "VSpre_2504322.off")
- "VSpost_500_poi_on_2504322.txt" (POI list with POI number of 500 on "VSpost_2504322.off")
- "VSpre_1000_poi_on_500000.txt" (POI list with POI number of 1000 on "VSpre_500000.off")
- "VSpost_1000_poi_on_500000.txt" (POI list with POI number of 1000 on "VSpost_500000.off")
- "VSpre_1500_poi_on_500000.txt" (POI list with POI number of 1500 on "VSpre_500000.off")
- "VSpost_1500_poi_on_500000.txt" (POI list with POI number of 1500 on "VSpost_500000.off")
- "VSpre_2000_poi_on_500000.txt" (POI list with POI number of 2000 on "VSpre_500000.off")
- "VSpost_2000_poi_on_500000.txt" (POI list with POI number of 2000 on "VSpost_500000.off")
- "VSpre_2500_poi_on_500000.txt" (POI list with POI number of 2500 on "VSpre_500000.off")
- "VSpost_2500_poi_on_500000.txt" (POI list with POI number of 2500 on "VSpost_500000.off")
- "VSpre_50_poi_on_500000.txt" (POI list with POI number of 50 on "VSpre_500000.off")
- "VSpost_50_poi_on_500000.txt" (POI list with POI number of 50 on "VSpost_500000.off")
- "VSpre_50_poi_on_1002528.txt" (POI list with POI number of 50 on "VSpre_1002528.off")
- "VSpost_50_poi_on_1002528.txt" (POI list with POI number of 50 on "VSpost_1002528.off")
- "VSpre_50_poi_on_1503378.txt" (POI list with POI number of 50 on "VSpre_1503378.off")
- "VSpost_50_poi_on_1503378.txt" (POI list with POI number of 50 on "VSpost_1503378.off")
- "VSpre_50_poi_on_2000000.txt" (POI list with POI number of 50 on "VSpre_2000000.off")
- "VSpost_50_poi_on_2000000.txt" (POI list with POI number of 50 on "VSpost_2000000.off")
- "VSpre_50_poi_on_2504322.txt" (POI list with POI number of 50 on "VSpre_2504322.off")
- "VSpost_50_poi_on_2504322.txt" (POI list with POI number of 50 on "VSpost_2504322.off")
- "VSpre_100_poi_on_500000.txt" (POI list with POI number of 100 on "VSpre_500000.off")
- "VSpost_100_poi_on_500000.txt" (POI list with POI number of 100 on "VSpost_500000.off")
- "VSpre_100_poi_on_500000.txt" (POI list with POI number of 150 on "VSpre_500000.off")
- "VSpost_100_poi_on_500000.txt" (POI list with POI number of 150 on "VSpost_500000.off")
- "VSpre_200_poi_on_500000.txt" (POI list with POI number of 200 on "VSpre_500000.off")
- "VSpost_200_poi_on_500000.txt" (POI list with POI number of 200 on "VSpost_500000.off")
- "VSpre_250_poi_on_500000.txt" (POI list with POI number of 250 on "VSpre_500000.off")
- "VSpost_250_poi_on_500000.txt" (POI list with POI number of 250 on "VSpost_500000.off")

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

- [terrain_data_and_dataset_size_and_poi_number_map_index]: an index for the map of terrain data and dataset size and poi number (a integer from 0 to 45)
- [epsilon]: the epsilon value (0 < epsilon <= 1)

For the [terrain_data_and_dataset_size_and_poi_number_map_index], each index value corresponding to a terrain data, the dataset size of the terrain and the poi number on the terrain, their relationships are as follows:

| Index | Terrain data | Dataset size | POI number |
| ----------- | ----------- | ----------- | ----------- |
| 0 | SC | 500000 | 50 |
| 1 | SC | 500000 | 100 |
| 2 | SC | 500000 | 150 |
| 3 | SC | 500000 | 200 |
| 4 | SC | 500000 | 250 |
| 5 | SC | 1002528 | 50 |
| 6 | SC | 1503378 | 50 |
| 7 | SC | 2000000 | 50 |
| 8 | SC | 2504322 | 50 |
| 9 | AU | 500000 | 50 |
| 10 | AU | 500000 | 100 |
| 11 | AU | 500000 | 150 |
| 12 | AU | 500000 | 200 |
| 13 | AU | 500000 | 250 |
| 14 | AU | 1002528 | 50 |
| 15 | AU | 1503378 | 50 |
| 16 | AU | 2000000 | 50 |
| 17 | AU | 2504322 | 50 |
| 18 | VS | 500000 | 50 |
| 19 | VS | 500000 | 100 |
| 20 | VS | 500000 | 150 |
| 21 | VS | 500000 | 200 |
| 22 | VS | 500000 | 250 |
| 23 | VS | 1002528 | 50 |
| 24 | VS | 1503378 | 50 |
| 25 | VS | 2000000 | 50 |
| 26 | VS | 2504322 | 50 |
| 27 | SC | 500000 | 50 |
| 28 | SC | 500000 | 1000 |
| 29 | SC | 500000 | 1500 |
| 30 | SC | 500000 | 2000 |
| 31 | SC | 500000 | 2500 |
| 32 | SC | 1002528 | 500 |
| 33 | SC | 1503378 | 500 |
| 34 | SC | 2000000 | 500 |
| 35 | SC | 2504322 | 500 |
| 36 | AU | 500000 | 500 |
| 37 | AU | 500000 | 1000 |
| 38 | AU | 500000 | 1500 |
| 39 | AU | 500000 | 2000 |
| 40 | AU | 500000 | 2500 |
| 41 | AU | 1002528 | 500 |
| 42 | AU | 1503378 | 500 |
| 43 | AU | 2000000 | 500 |
| 44 | AU | 2504322 | 500 |
| 45 | VS | 500000 | 500 |
| 46 | VS | 500000 | 1000 |
| 47 | VS | 500000 | 1500 |
| 48 | VS | 500000 | 2000 |
| 49 | VS | 500000 | 2500 |
| 50 | VS | 1002528 | 500 |
| 51 | VS | 1503378 | 500 |
| 52 | VS | 2000000 | 500 |
| 53 | VS | 2504322 | 500 |

By default, the project will run FU-Oracle-Naive1, FU-Oracle-Naive2, FU-Oracle, WSPD-Oracle, WSPD-Oracle-Adapt, and K-Algo. But as mentioned in our paper, WSPD-Oracle and WSPD-Oracle-Adapt are very time consuming. So when the dataset size is large, i.e., [terrain_data_and_dataset_size_and_poi_number_map_index] > 27, the project will only run FU-Oracle-Naive1, FU-Oracle-Naive2, FU-Oracle, and K-Algo.

An example:

```
./main 0 0.5
```

In this example, [terrain_data_and_dataset_size_and_poi_number_map_index] is 0, [epsilon] is 0.5. So, it will run SC pre earthquake terrain dataset and SC post earthquake terrain dataset, with dataset size equal to 500000 and poi number equal to 50, and epsilon is 0.5. It will run six algorithms, i.e., FU-Oracle-Naive1, FU-Oracle-Naive2, FU-Oracle, WSPD-Oracle, WSPD-Oracle-Adapt, and K-Algo.

## Output

The output will be stored in "output/output.txt" file. The format will be as follows:

```

[pre_dataset] [post_dataset] [datasize] [poi_num] [epsilon] [pre_construction_time (ms)] [pre_memory_usage (MB)]  [pre_index_size (MB)] [pre_index_edge_num] [pre_index/MST_weight] [post_update_time1 (ms)] [post_update_time2 (ms)] [post_query_time (ms)] [post_memory_usage (MB)] [post_index_size (MB)] [post_output_size (MB)] [post_index_edge_num] [post_index/MST_weight] [post_distance_error]

```

These information will also be shown in the terminal. 