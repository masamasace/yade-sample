#include <iostream>
#include <vector>
#include <cuda_runtime.h>

// 定数としてボクセルのサイズを設定
const float VOXEL_SIZE = 1.0;

// CUDAカーネル関数
// 各ボクセルが粒状体によって占められているかどうかをチェックする
__global__ void checkVoxels(float *particles, int num_particles, float *voxel_centers, int num_voxels, int *occupied) {
    // 現在のスレッドのインデックスを計算
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    // インデックスがボクセルの総数以上の場合、このスレッドは何もしない
    if (idx >= num_voxels) return;

    // 現在のボクセルの中心座標を取得
    float voxel_x = voxel_centers[3 * idx];
    float voxel_y = voxel_centers[3 * idx + 1];
    float voxel_z = voxel_centers[3 * idx + 2];

    // 全ての粒状体に対して、ボクセルが占められているかどうかをチェック
    for (int i = 0; i < num_particles; i++) {
        // 粒状体の中心座標と半径を取得
        float px = particles[5 * i + 1];
        float py = particles[5 * i + 2];
        float pz = particles[5 * i + 3];
        float radius = particles[5 * i + 4];

        // ボクセルの中心と粒状体の中心との距離の二乗を計算
        float dist2 = (voxel_x - px) * (voxel_x - px) + (voxel_y - py) * (voxel_y - py) + (voxel_z - pz) * (voxel_z - pz);

        // 距離が粒状体の半径以下ならボクセルは占められていると判断
        if (dist2 <= radius * radius) {
            occupied[idx] = 1;
            return;
        }
    }

    // どの粒状体にも占められていない場合、occupiedを0に設定
    occupied[idx] = 0;
}

int main() {
    std::vector<float> particles; // [id, x, y, z, radius]
    // データの読み込み...

    // 部分区画の座標を設定
    float region_start_x = ...;
    float region_end_x = ...;
    // y, zも同様に

    // 部分区画内のボクセルの総数を計算
    int num_voxels = ...;
    // 部分区画内の各ボクセルの中心座標を保存するための配列を準備
    std::vector<float> voxel_centers(3 * num_voxels);
    // 中心座標の計算...

    // GPUメモリ上にボクセルが占められているかどうかを保存するための配列を確保
    int *d_occupied;
    cudaMalloc(&d_occupied, num_voxels * sizeof(int));

    // GPUメモリ上に粒状体とボクセルのデータをコピーするためのポインタを確保
    float *d_particles, *d_voxel_centers;
    cudaMalloc(&d_particles, particles.size() * sizeof(float));
    cudaMalloc(&d_voxel_centers, voxel_centers.size() * sizeof(float));

    // ホストメモリからGPUメモリへデータをコピー
    cudaMemcpy(d_particles, particles.data(), particles.size() * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_voxel_centers, voxel_centers.data(), voxel_centers.size() * sizeof(float), cudaMemcpyHostToDevice);

    // スレッドの数とブロックの数を設定
    int threads_per_block = 256;
    int blocks = (num_voxels + threads_per_block - 1) / threads_per_block;

    // CUDAカーネル関数を呼び出し
    checkVoxels<<<blocks, threads_per_block>>>(d_particles, particles.size() / 5, d_voxel_centers, num_voxels, d_occupied);

    // GPUメモリからホストメモリへ結果をコピー
    std::vector<int> occupied(num_voxels);
    cudaMemcpy(occupied.data(), d_occupied, num_voxels * sizeof(int), cudaMemcpyDeviceToHost);

    // 空隙率を計算
    int sum = 0;
    for (int o : occupied) {
        sum += o;
    }
    float porosity = 1.0 - float(sum) / num_voxels;
    std::cout << "Porosity: " << porosity << std::endl;

    // GPUメモリを解放
    cudaFree(d_particles);
    cudaFree(d_voxel_centers);
    cudaFree(d_occupied);

    return 0;
}
