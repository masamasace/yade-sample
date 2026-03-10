#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <filesystem>
#include <list>
#include <algorithm>
#include <cmath>
#include <sys/stat.h>
#include <cuda_runtime.h>


namespace fs = std::filesystem;

__global__ void checkVoxels(float *particles, int num_particles, float *voxel_centers, int num_voxels, int *occupied) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_voxels) return;

    float voxel_x = voxel_centers[3 * idx];
    float voxel_y = voxel_centers[3 * idx + 1];
    float voxel_z = voxel_centers[3 * idx + 2];

    for (int i = 0; i < num_particles; i++) {
        float radius = particles[4 * i];
        float px = particles[4 * i + 1];
        float py = particles[4 * i + 2];
        float pz = particles[4 * i + 3];

        float dist2 = (voxel_x - px) * (voxel_x - px) + (voxel_y - py) * (voxel_y - py) + (voxel_z - pz) * (voxel_z - pz);
        if (dist2 <= radius * radius) {
            occupied[idx] = 1;
            return;
        }
    }
    occupied[idx] = 0;
}

std::vector<float> readCSV(const std::string &filename) {
    std::vector<float> particles;
    std::ifstream file(filename);
    std::string line, value, header_label;
    std::vector<std::string> target_header_labels;
    std::vector<int> target_header_label_indexes;
    
    target_header_labels = {"radii", "Points:0", "Points:1", "Points:2"};
    
    std::getline(file, line);
    std::istringstream header(line);
    int cur_index = 0;

    while (std::getline(header, header_label, ',')) {

        for (std::string target_header_label : target_header_labels) {

            if (header_label.find(target_header_label) != std::string::npos) {
                target_header_label_indexes.push_back(cur_index);
            }
        }

        cur_index++;
    }

    while (std::getline(file, line)) {
        std::istringstream iss(line);

        std::vector<std::string> values;
        while (std::getline(iss, value, ',')) {
            values.push_back(value);
        }

        particles.push_back(std::stof(values[target_header_label_indexes[0]]));
        particles.push_back(std::stof(values[target_header_label_indexes[1]]));
        particles.push_back(std::stof(values[target_header_label_indexes[2]]));
        particles.push_back(std::stof(values[target_header_label_indexes[3]])); 
    }
    return particles;
}

int main(int argc, char* argv[]) {

    const std::list<std::list<std::list<float>>> regions = {
        {{0.005, 0, 0.015}, {0.015, 0.01, 0.025}},
        {{0.005, 0.09, 0.015}, {0.015, 0.1, 0.025}},
    };

    const int voxel_resolution = 200;
    const std::filesystem::path target_dir = argv[1];
    const std::filesystem::path dir_output_path = target_dir.parent_path() / "VOXEL";
    std::filesystem::create_directory(dir_output_path);

    std::vector<std::filesystem::path> sorted_file_paths;

    for (const auto &entry : fs::directory_iterator(target_dir)) {
        sorted_file_paths.push_back(entry.path());
    }

    std::sort(sorted_file_paths.begin(), sorted_file_paths.end());
    
    for (const auto &path : sorted_file_paths) {
        if (path.extension() == ".csv") {
            std::vector<float> particles = readCSV(path);
            int counter_region = 0;

            for (const auto &region : regions) {  // regionsは部分区画のリスト

                std::vector<float> region_min_vec(region.front().begin(), region.front().end());
                std::vector<float> region_max_vec(region.back().begin(), region.back().end());
            
                float offset_x = region_max_vec[0] - region_min_vec[0];
                float offset_y = region_max_vec[1] - region_min_vec[1];
                float offset_z = region_max_vec[2] - region_min_vec[2];

                std::vector<float> voxel_centers(std::pow(voxel_resolution, 3) * 3);

                for (int i = 0; i < voxel_resolution; i++) {
                    for (int j = 0; j < voxel_resolution; j++) {
                        for (int k = 0; k < voxel_resolution; k++) {
                            voxel_centers[(i * std::pow(voxel_resolution, 2) + j * voxel_resolution + k) * 3] = 
                                region_min_vec[0] + float(i + 0.5) * offset_x / voxel_resolution;
                            voxel_centers[(i * std::pow(voxel_resolution, 2) + j * voxel_resolution + k) * 3 + 1] = 
                                region_min_vec[1] + float(j + 0.5) * offset_y / voxel_resolution;
                            voxel_centers[(i * std::pow(voxel_resolution, 2) + j * voxel_resolution + k) * 3 + 2] = 
                                region_min_vec[2] + float(k + 0.5) * offset_z / voxel_resolution;
                        }
                    }
                }

                // CUDA処理...
                int *d_occupied;
                cudaMalloc(&d_occupied, voxel_centers.size() * sizeof(int));

                float *d_particles, *d_voxel_centers;
                cudaMalloc(&d_particles, particles.size() * sizeof(float));
                cudaMalloc(&d_voxel_centers, voxel_centers.size() * sizeof(float));

                cudaMemcpy(d_particles, particles.data(), particles.size() * sizeof(float), cudaMemcpyHostToDevice);
                cudaMemcpy(d_voxel_centers, voxel_centers.data(), voxel_centers.size() * sizeof(float), cudaMemcpyHostToDevice);

                // スレッドの数とブロックの数を設定
                int threads_per_block = 1024;
                int blocks = (voxel_centers.size() / 3 + threads_per_block - 1) / threads_per_block;

                // CUDAカーネル関数を呼び出し
                checkVoxels<<<blocks, threads_per_block>>>(d_particles, particles.size() / 4, d_voxel_centers, voxel_centers.size() / 3 , d_occupied);

                // GPUメモリからホストメモリへ結果をコピー
                std::vector<int> occupied(voxel_centers.size());
                cudaMemcpy(occupied.data(), d_occupied, voxel_centers.size() * sizeof(int), cudaMemcpyDeviceToHost);
                
                // 空隙率を計算...
                int sum = 0;
                for (int o : occupied) {
                    sum += o;
                }
                float porosity = 1.0 - float(sum) / voxel_centers.size() ;
                std::cout << path.filename() << ", " << counter_region << " : " << porosity << std::endl;

                // GPUメモリを解放
                cudaFree(d_particles);
                cudaFree(d_voxel_centers);
                cudaFree(d_occupied);

                counter_region++;

            }
        }
    }
    return 0;

}
