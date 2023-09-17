#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <filesystem>
#include <list>
#include <cuda_runtime.h>


namespace fs = std::filesystem;

__global__ void checkVoxels(float *particles, int num_particles, float *voxel_centers, int num_voxels, int *occupied) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_voxels) return;

    float voxel_x = voxel_centers[3 * idx];
    float voxel_y = voxel_centers[3 * idx + 1];
    float voxel_z = voxel_centers[3 * idx + 2];

    for (int i = 0; i < num_particles; i++) {
        float px = particles[4 * i];
        float py = particles[4 * i + 1];
        float pz = particles[4 * i + 2];
        float radius = particles[4 * i + 3];

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
    std::string line, value;
    std::getline(file, line);  // headerを読み飛ばす

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<std::string> values;
        while (std::getline(iss, value, ',')) {
            values.push_back(value);
        }
        particles.push_back(std::stof(values[0]));  // Points:0
        particles.push_back(std::stof(values[1]));  // Points:1
        particles.push_back(std::stof(values[2]));  // Points:2
        particles.push_back(std::stof(values[3]));  // radii

        std::cout << values[3] << std::endl;
    }
    return particles;
}

int main() {
    const std::list<std::list<std::list<float>>> regions = {
        {{0.005, 0, 0.015}, {0.015, 0.01, 0.025}},
        {{0.005, 0.01, 0.015}, {0.015, 0.02, 0.025}},
    };

    const int voxel_resolution = 100;
    
    for (const auto &entry : fs::directory_iterator(".")) {
        if (entry.path().extension() == ".csv") {
            std::vector<float> particles = readCSV(entry.path());

            for (const auto &region : regions) {  // regionsは部分区画のリスト

                auto [region_min_vec, region_max_vec] = region

                
                // 部分区画の座標とボクセルの総数を計算...

                // CUDA処理...

                // 空隙率を計算...
                int sum = 0;
                for (int o : occupied) {
                    sum += o;
                }
                float porosity = 1.0 - float(sum) / num_voxels;
                std::cout << entry.path() << ": " << porosity << std::endl;
            }
        }
    }
    return 0;
}
