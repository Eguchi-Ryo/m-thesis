#include <iostream>
#include <fstream>
#include <vector>
#include <array>

std::array<double, 3> transformToRightHanded(const std::array<double, 3>& displacement) {
    // 左手系から右手系への変換
	return { displacement[0], displacement[1], displacement[2] };
    //return { displacement[0], displacement[1], displacement[2] };  //patient, (A,S,C)=(96, 200, 180)の時完全一致
}

int main() {
    // データのサイズ情報
    int width = 576;
    int height = 576;
    int depth = 233;

    // 3次元配列の作成（各ボクセルに3つの値: x, y, zの変位）
    std::vector<std::vector<std::vector<std::array<double, 3>>>> displacement(
        depth, std::vector<std::vector<std::array<double, 3>>>(
            height, std::vector<std::array<double, 3>>(width)));

    // RAWバイナリファイルを読み込む
    std::ifstream file("/Users/lelab/Downloads/exam/patient/debug/registration/deformationField_3D.raw", std::ios::binary);
    if (!file) {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }

    // ファイルからデータを3次元配列に読み込む
    for (int z = 0; z < depth; ++z) {
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
				std::array<double, 3> df;
				file.read(reinterpret_cast<char*>(&df), sizeof(double) *3);
                displacement[z][y][x] = transformToRightHanded(df); 
            }
        }
    }
    file.close();

    // 特定の位置[0,0,0]の値を取得
    auto displacement1 = displacement[66][235][282];
    std::cout << "Displacement at [279,244,66], (dx,dy,dz)= (" 
              << displacement1[0] << ", " 
              << displacement1[1] << ", " 
              << displacement1[2] << ")" << std::endl;

    auto displacement2 = displacement[66][244][283];
    std::cout << "Displacement at [283,244,66], (dx,dy,dz)= (" 
              << displacement2[0] << ", " 
              << displacement2[1] << ", " 
              << displacement2[2] << ")" << std::endl;

    auto displacement3 = displacement[232][0][0];
    std::cout << "Displacement at [233, 0, 0], (dx,dy,dz)= ("
            << displacement3[0] << ", "
            << displacement3[1] << ", "
            << displacement3[2] << ",) " << std::endl;

    auto displacement4 = displacement[116][287][287];
    std::cout << "Displacement at [116, 287, 287], (dx,dy,dz)= ("
            << displacement4[0] << ", "
            << displacement4[1] << ", "
            << displacement4[2] << ",) " << std::endl;

    return 0;
}
