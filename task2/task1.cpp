#include <cmath>
#include <vector>
#include <iostream>
#include <functional>
#include <fstream>
#include <algorithm>

class Matrix {
private:
    std::vector<std::vector<double>> matrix;
    std::vector<double> rightPart;
    std::vector<double> x;
    std::vector<double> f;
    std::vector<double> h;
    int n;
    
public:
    Matrix(std::vector<double> x, std::vector<double> f);
    void printMatrix();
    std::vector<double> prognat();
    void printFunction(std::ofstream& output, std::vector<double>& gamma);
};

Matrix::Matrix(std::vector<double> x, std::vector<double> f) {
    this->x = x;
    this->f = f;
    this->n = x.size();
    
    h.resize(n - 1);
    for (int i = 0; i < n - 1; i++) {
        h[i] = x[i + 1] - x[i];
    }

    int systemSize = n - 2;
    matrix.resize(systemSize);
    for (int i = 0; i < systemSize; i++) {
        matrix[i].resize(systemSize);
    }
    rightPart.resize(systemSize);
    
    for (int i = 0; i < systemSize; i++) {
        int idx = i + 1;
        
        matrix[i][i] = 2.0 * (h[idx - 1] + h[idx]) / 3.0;

        if (i > 0) {
            matrix[i][i - 1] = h[idx - 1] / 3.0;
        }
        if (i < systemSize - 1) {
            matrix[i][i + 1] = h[idx] / 3.0;
        }
        
        rightPart[i] = (f[idx + 1] - f[idx]) / h[idx] - (f[idx] - f[idx - 1]) / h[idx - 1];
    }
}

void Matrix::printMatrix() {
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << "| " << rightPart[i] << "\n";
    }
}

std::vector<double> Matrix::prognat() {
    int size = matrix.size();
    
    if (size == 0) {
        return std::vector<double>();
    }
    
    std::vector<double> ksi(size);
    std::vector<double> etta(size);
    std::vector<double> x(size);

    ksi[0] = -matrix[0][1] / matrix[0][0];
    etta[0] = rightPart[0] / matrix[0][0];
    
    for (int i = 1; i < size - 1; i++) {
        double denominator = matrix[i][i] + matrix[i][i-1] * ksi[i-1];
        ksi[i] = -matrix[i][i+1] / denominator;
        etta[i] = (rightPart[i] - matrix[i][i-1] * etta[i-1]) / denominator;
    }

    double last_denominator = matrix[size-1][size-1] + matrix[size-1][size-2] * ksi[size-2];
    etta[size-1] = (rightPart[size-1] - matrix[size-1][size-2] * etta[size-2]) / last_denominator;

    x[size-1] = etta[size-1];
    for (int i = size - 2; i >= 0; i--) {
        x[i] = ksi[i] * x[i+1] + etta[i];
    }
    
    return x;
}

void Matrix::printFunction(std::ofstream& output, std::vector<double>& gamma) {
    const int numPoints = 1000;
    
    std::vector<double> fullGamma(n);
    fullGamma[0] = 0.0;
    for (int i = 0; i < gamma.size(); i++) {
        fullGamma[i + 1] = gamma[i];
    }
    fullGamma[n - 1] = 0.0;
    
    for (int i = 0; i < n - 1; i++) {
        double a = f[i];
        double b = (f[i + 1] - f[i]) / h[i] - h[i] * (fullGamma[i + 1] + 2 * fullGamma[i]) / 3.0;
        double c = fullGamma[i];
        double d = (fullGamma[i + 1] - fullGamma[i]) / (3.0 * h[i]);

        int pointsPerInterval = numPoints / (n - 1);
        for (int j = 0; j <= pointsPerInterval; j++) {
            double t = static_cast<double>(j) / pointsPerInterval;
            double x_val = x[i] + t * h[i];
            double dx = x_val - x[i];
            double s_x = a + b * dx + c * dx * dx + d * dx * dx * dx;
            
            output << x_val << " " << s_x << "\n";
        }
    }
}

int main() {
    std::vector<double> x = {-1, 0, 2, 3, 5}; 
    std::vector<double> f = {1, 2, 4, 1, -3}; 
    
    Matrix m(x, f);
    std::vector<double> gamma = m.prognat();
    
    std::ofstream out("out.txt");
    m.printFunction(out, gamma);
    out.close();
    
    std::cout << "Кубический сплайн построен и сохранен в out.txt\n";
    return 0;
}
