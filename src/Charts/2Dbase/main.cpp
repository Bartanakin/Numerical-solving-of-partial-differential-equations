#include <matplot/matplot.h>
#include <filesystem>
#include <fstream>
#include <iostream>

using namespace matplot;

int main(int argc, char* argv[]) {
    std::filesystem::path file_path = "ChartData";
    file_path /= argv[1];
    std::ifstream file(file_path.string());
    if (file.fail()) {
        std::cout << "File corrupt ot doesn't exist." << std::endl;

        return 1;
    }

    int x_size;
    int y_size;
    file >> x_size;
    file >> y_size;

    vector_1d X_axis(x_size);
    vector_1d Y_axis(y_size);
    vector_2d Z(x_size, vector_1d(y_size, 0.));

    for (auto& x : X_axis) {
        file >> x;
    }
    
    for (auto& y : Y_axis) {
        file >> y;
    }

    for (auto& z_col : Z) {
        for (auto& z : z_col) {
            file >> z;
        }
    }

    auto [X, Y] = meshgrid(Y_axis, X_axis);
    // std::cout << X.size() << " " << Y.size() << " " << Z.size() << std::endl;
    auto s3 = mesh(X, Y, Z);

    show();


    return 0;
}