# Update Efficient Shortest Path Oracle on Dynamic Terrain Surface

## Overview

This project provides the implementation of the algorithm for calculating an update efficient shortest path on dynamic terrain surface.

Our oralce UE, and the baselines, i.e., WSPD-oracle, SP-oracle, and KF are studied in the experiments. In order to conduct the ablation study, i.e., show that algorithm HGS could significantly reduce the running time compared with algorithm GS, and could significantly reduce the oracle size and oracle weight compared with original complete graph, we also studied UE-N1 and UE-N2 in the experiments. In total, we compared six algorithms, namely, WSPD-oracle, SP-oracle, UE-N1, UE-N2, UE, and KF. Since WSPD-oracle, SP-oracle are not feasible on large datasets due to their expensive running time, so we (1) compared these six algorithms on SC-small, AU-small and VS-small datasets, and the set of small-version datasets, and (2) compared UE-N1, UE-N2, UE, and KF on SC, AU and VS datasets, and the set of large-version datasets. We refer the readers to our paper for more details.

In total, we compared six algorithms as follows:

- WSPD-oracle (oracle based baseline)
- SP-oracle (oracle based baseline)
- UE-N1 (variation)
- UE-N2 (variation)
- UE (our oracle)
- KF (on-the-fly baseline)

Make sure there is a folder called "input/" and a folder called "output/" under the working directory. They will be used for storing the input/output files.

The source code are stored in "src/" folder.

## Dataset

The dataset are stored in "input/" folder.

The datasets are as follows:

- "SCpre_500000.off" (large version default resolution SC pre earthquake terrain dataset with dataset size of 500000)
- "SCpost_500000.off" (large version default resolution SC post earthquake terrain dataset with dataset size of 500000)
- "SCpre_1002528.off" (large version multiresolution SC pre earthquake terrain dataset with dataset size of 1002528)
- "SCpost_1002528.off" (large version multiresolution SC post earthquake terrain dataset with dataset size of 1002528)
- "SCpre_1503378.off" (large version multiresolution SC pre earthquake terrain dataset with dataset size of 1503378)
- "SCpost_1503378.off" (large version multiresolution SC post earthquake terrain dataset with dataset size of 1503378)
- "SCpre_2000000.off" (large version multiresolution SC pre earthquake terrain dataset with dataset size of 2000000)
- "SCpost_2000000.off" (large version multiresolution SC post earthquake terrain dataset with dataset size of 2000000)
- "SCpre_2504322.off" (large version multiresolution SC pre earthquake terrain dataset with dataset size of 2504322)
- "SCpost_2504322.off" (large version multiresolution SC post earthquake terrain dataset with dataset size of 2504322)
- "SCpre_10082.off" (small version default resolution SC pre earthquake terrain dataset with dataset size of 10082)
- "SCpost_10082.off" (small version default resolution SC post earthquake terrain dataset with dataset size of 10082)
- "SCpre_20000.off" (small version multiresolution SC pre earthquake terrain dataset with dataset size of 20000)
- "SCpost_20000.off" (small version multiresolution SC post earthquake terrain dataset with dataset size of 20000)
- "SCpre_30258.off" (small version multiresolution SC pre earthquake terrain dataset with dataset size of 30258)
- "SCpost_30258.off" (small version multiresolution SC post earthquake terrain dataset with dataset size of 30258)
- "SCpre_40328.off" (small version multiresolution SC pre earthquake terrain dataset with dataset size of 40328)
- "SCpost_40328.off" (small version multiresolution SC post earthquake terrain dataset with dataset size of 40328)
- "SCpre_50562.off" (small version multiresolution SC pre earthquake terrain dataset with dataset size of 50562)
- "SCpost_50562.off" (small version multiresolution SC post earthquake terrain dataset with dataset size of 50562)
- "AUpre_500000.off" (large version default resolution AU pre earthquake terrain dataset with dataset size of 500000)
- "AUpost_500000.off" (large version default resolution AU post earthquake terrain dataset with dataset size of 500000)
- "AUpre_1002528.off" (large version multiresolution AU pre earthquake terrain dataset with dataset size of 1002528)
- "AUpost_1002528.off" (large version multiresolution AU post earthquake terrain dataset with dataset size of 1002528)
- "AUpre_1503378.off" (large version multiresolution AU pre earthquake terrain dataset with dataset size of 1503378)
- "AUpost_1503378.off" (large version multiresolution AU post earthquake terrain dataset with dataset size of 1503378)
- "AUpre_2000000.off" (large version multiresolution AU pre earthquake terrain dataset with dataset size of 2000000)
- "AUpost_2000000.off" (large version multiresolution AU post earthquake terrain dataset with dataset size of 2000000)
- "AUpre_2504322.off" (large version multiresolution AU pre earthquake terrain dataset with dataset size of 2504322)
- "AUpost_2504322.off" (large version multiresolution AU post earthquake terrain dataset with dataset size of 2504322)
- "AUpre_10082.off" (small version default resolution AU pre earthquake terrain dataset with dataset size of 10082)
- "AUpost_10082.off" (small version default resolution AU post earthquake terrain dataset with dataset size of 10082)
- "VSpre_500000.off" (large version default resolution VS pre earthquake terrain dataset with dataset size of 500000)
- "VSpost_500000.off" (large version default resolution VS post earthquake terrain dataset with dataset size of 500000)
- "VSpre_1002528.off" (large version multiresolution VS pre earthquake terrain dataset with dataset size of 1002528)
- "VSpost_1002528.off" (large version multiresolution VS post earthquake terrain dataset with dataset size of 1002528)
- "VSpre_1503378.off" (large version multiresolution VS pre earthquake terrain dataset with dataset size of 1503378)
- "VSpost_1503378.off" (large version multiresolution VS post earthquake terrain dataset with dataset size of 1503378)
- "VSpre_2000000.off" (large version multiresolution VS pre earthquake terrain dataset with dataset size of 2000000)
- "VSpost_2000000.off" (large version multiresolution VS post earthquake terrain dataset with dataset size of 2000000)
- "VSpre_2504322.off" (large version multiresolution VS pre earthquake terrain dataset with dataset size of 2504322)
- "VSpost_2504322.off" (large version multiresolution VS post earthquake terrain dataset with dataset size of 2504322)
- "VSpre_10082.off" (small version default resolution VS pre earthquake terrain dataset with dataset size of 10082)
- "VSpost_10082.off" (small version default resolution VS post earthquake terrain dataset with dataset size of 10082)
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
- "SCpre_50_poi_on_10082.txt" (POI list with POI number of 50 on "SCpre_10082.off")
- "SCpost_50_poi_on_10082.txt" (POI list with POI number of 50 on "SCpost_10082.off")
- "SCpre_50_poi_on_20000.txt" (POI list with POI number of 50 on "SCpre_20000.off")
- "SCpost_50_poi_on_20000.txt" (POI list with POI number of 50 on "SCpost_20000.off")
- "SCpre_50_poi_on_30258.txt" (POI list with POI number of 50 on "SCpre_30258.off")
- "SCpost_50_poi_on_30258.txt" (POI list with POI number of 50 on "SCpost_30258.off")
- "SCpre_50_poi_on_40328.txt" (POI list with POI number of 50 on "SCpre_40328.off")
- "SCpost_50_poi_on_40328.txt" (POI list with POI number of 50 on "SCpost_40328.off")
- "SCpre_50_poi_on_50562.txt" (POI list with POI number of 50 on "SCpre_50562.off")
- "SCpost_50_poi_on_50562.txt" (POI list with POI number of 50 on "SCpost_50562.off")
- "SCpre_100_poi_on_10082.txt" (POI list with POI number of 100 on "SCpre_10082.off")
- "SCpost_100_poi_on_10082.txt" (POI list with POI number of 100 on "SCpost_10082.off")
- "SCpre_150_poi_on_10082.txt" (POI list with POI number of 150 on "SCpre_10082.off")
- "SCpost_150_poi_on_10082.txt" (POI list with POI number of 150 on "SCpost_10082.off")
- "SCpre_200_poi_on_10082.txt" (POI list with POI number of 200 on "SCpre_10082.off")
- "SCpost_200_poi_on_10082.txt" (POI list with POI number of 200 on "SCpost_10082.off")
- "SCpre_250_poi_on_10082.txt" (POI list with POI number of 250 on "SCpre_10082.off")
- "SCpost_250_poi_on_10082.txt" (POI list with POI number of 250 on "SCpost_10082.off")
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
- "AUpre_50_poi_on_10082.txt" (POI list with POI number of 50 on "AUpre_10082.off")
- "AUpost_50_poi_on_10082.txt" (POI list with POI number of 50 on "AUpost_10082.off")
- "AUpre_100_poi_on_10082.txt" (POI list with POI number of 100 on "AUpre_10082.off")
- "AUpost_100_poi_on_10082.txt" (POI list with POI number of 100 on "AUpost_10082.off")
- "AUpre_150_poi_on_10082.txt" (POI list with POI number of 150 on "AUpre_10082.off")
- "AUpost_150_poi_on_10082.txt" (POI list with POI number of 150 on "AUpost_10082.off")
- "AUpre_200_poi_on_10082.txt" (POI list with POI number of 200 on "AUpre_10082.off")
- "AUpost_200_poi_on_10082.txt" (POI list with POI number of 200 on "AUpost_10082.off")
- "AUpre_250_poi_on_10082.txt" (POI list with POI number of 250 on "AUpre_10082.off")
- "AUpost_250_poi_on_10082.txt" (POI list with POI number of 250 on "AUpost_10082.off")
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
- "VSpre_50_poi_on_10082.txt" (POI list with POI number of 50 on "VSpre_10082.off")
- "VSpost_50_poi_on_10082.txt" (POI list with POI number of 50 on "VSpost_10082.off")
- "VSpre_100_poi_on_10082.txt" (POI list with POI number of 100 on "VSpre_10082.off")
- "VSpost_100_poi_on_10082.txt" (POI list with POI number of 100 on "VSpost_10082.off")
- "VSpre_150_poi_on_10082.txt" (POI list with POI number of 150 on "VSpre_10082.off")
- "VSpost_150_poi_on_10082.txt" (POI list with POI number of 150 on "VSpost_10082.off")
- "VSpre_200_poi_on_10082.txt" (POI list with POI number of 200 on "VSpre_10082.off")
- "VSpost_200_poi_on_10082.txt" (POI list with POI number of 200 on "VSpost_10082.off")
- "VSpre_250_poi_on_10082.txt" (POI list with POI number of 250 on "VSpre_10082.off")
- "VSpost_250_poi_on_10082.txt" (POI list with POI number of 250 on "VSpost_10082.off")

Since the file size for the dataset "SCpre_500000.off", "SCpost_500000.off", "SCpre_1002528.off", "SCpost_1002528.off", "SCpre_1503378.off", "SCpost_1503378.off", "SCpre_2000000.off", "SCpost_2000000.off", "SCpre_2504322.off", "SCpost_2504322.off","AUpre_500000.off", "AUpost_500000.off", "AUpre_1002528.off", "AUpost_1002528.off", "AUpre_1503378.off", "AUpost_1503378.off", "AUpre_2000000.off", "AUpost_2000000.off", "AUpre_2504322.off", "AUpost_2504322.off", "VSpre_500000.off", "VSpost_500000.off", "VSpre_1002528.off", "VSpost_1002528.off", "VSpre_1503378.off", "VSpost_1503378.off", "VSpre_2000000.off", "VSpost_2000000.off", "VSpre_2504322.off", and "VSpost_2504322.off" are too large (i.e., they exceed the maximum file size for Github), please download these five files at https://, and put them back in "input/" folder.

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
| 0 | SC | 10082 | 50 |
| 1 | SC | 10082 | 100 |
| 2 | SC | 10082 | 150 |
| 3 | SC | 10082 | 200 |
| 4 | SC | 10082 | 250 |
| 5 | SC | 20000 | 50 |
| 6 | SC | 30258 | 50 |
| 7 | SC | 40328 | 50 |
| 8 | SC | 50562 | 50 |
| 9 | AU | 10082 | 50 |
| 10 | AU | 10082 | 100 |
| 11 | AU | 10082 | 150 |
| 12 | AU | 10082 | 200 |
| 13 | AU | 10082 | 250 |
| 14 | VS | 10082 | 50 |
| 15 | VS | 10082 | 100 |
| 16 | VS | 10082 | 150 |
| 17 | VS | 10082 | 200 |
| 18 | VS | 10082 | 250 |
| 19 | SC | 500000 | 500 |
| 20 | SC | 500000 | 1000 |
| 21 | SC | 500000 | 1500 |
| 22 | SC | 500000 | 2000 |
| 23 | SC | 500000 | 2500 |
| 24 | SC | 1002528 | 500 |
| 25 | SC | 1503378 | 500 |
| 26 | SC | 2000000 | 500 |
| 27 | SC | 2504322 | 500 |
| 28 | AU | 500000 | 500 |
| 29 | AU | 500000 | 1000 |
| 30 | AU | 500000 | 1500 |
| 31 | AU | 500000 | 2000 |
| 32 | AU | 500000 | 2500 |
| 33 | AU | 1002528 | 500 |
| 34 | AU | 1503378 | 500 |
| 35 | AU | 2000000 | 500 |
| 36 | AU | 2504322 | 500 |
| 37 | VS | 500000 | 500 |
| 38 | VS | 500000 | 1000 |
| 39 | VS | 500000 | 1500 |
| 40 | VS | 500000 | 2000 |
| 41 | VS | 500000 | 2500 |
| 42 | VS | 1002528 | 500 |
| 43 | VS | 1503378 | 500 |
| 44 | VS | 2000000 | 500 |
| 45 | VS | 2504322 | 500 |

By default, the project will run UE-N1, UE-N2, UE, WSPD-oracle, SP-oracle, and KF. But as mentioned in our paper, WSPD_oracle and SP_oracle are very time consuming. So when the dataset size is large, i.e., [terrain_data_and_dataset_size_and_poi_number_map_index] > 18, the project will only run UE-N1, UE-N2, UE, and KF.

An example:

```
./main 0 0.5
```

In this example, [terrain_data_and_dataset_size_and_poi_number_map_index] is 0, [epsilon] is 0.5. So, it will run SC pre earthquake terrain dataset and SC post earthquake terrain dataset, with dataset size equal to 10082 and poi number equal to 50, and epsilon is 0.5. It will run siz algorithms, i.e., UE-N1, UE-N2, UE, WSPD-oracle, SP-oracle, and KF.

## Output

The output will be stored in "output/output.txt" file. The format will be as follows:

```

[pre_dataset] [post_dataset] [datasize] [poi_num] [epsilon] [pre_preprocessing_time1 (ms)] [pre_preprocessing_time2 (ms)] [pre_query_time (ms)] [pre_memory_usage1 (MB)] [pre_memory_usage2 (MB)] [pre_index_size] [pre_index_edge_num] [pre_index/MST_weight] [pre_distance_error] [post_updating_time1 (ms)] [post_updating_time2 (ms)] [post_query_time (ms)] [post_memory_usage1 (MB)] [post_memory_usage2 (MB)] [post_index_size] [post_index_edge_num] [post_index/MST_weight] [post_distance_error]

```

These information will also be shown in the terminal. 