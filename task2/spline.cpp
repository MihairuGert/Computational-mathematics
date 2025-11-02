#include <cmath>
#include <vector>
#include <iostream>
#include <functional>
#include <fstream>

class Matrix {
private:
    std::vector<std::vector<double>> matrix;
    std::vector<double> rightPart;
    std::function<double(double x)> f;
    
    double h;
    double a;
    double b;
    int n;
public:
    Matrix(int n, double a, double b);
    void printMatrix();
    std::vector<double> prognat();
    void setFunction(std::function<double(double x)> f);
    void printFunction(std::ofstream& output, std::vector<double>& gamma);
};

void Matrix::setFunction(std::function<double(double x)> f) {
    this->f = f;
    rightPart.resize(n - 2); 
    
    for (int i = 0; i < n - 2; i++) {
        double x_im1 = a + (i) * h;     
        double x_i = a + (i + 1) * h;   
        double x_ip1 = a + (i + 2) * h; 
        
        rightPart[i] = (f(x_ip1) - 2*f(x_i) + f(x_im1)) / h;
    }
}

Matrix::Matrix(int n, double a, double b) {
    this->a = a;
    this->b = b;
    this->n = n;
    
    if (b <= a || n <= 0) {
        throw 1;
    }

    h = (b - a) / (n - 1);
    
    int systemSize = n - 2;
    matrix.resize(systemSize);
    for (int i = 0; i < systemSize; i++) {
        matrix[i].resize(systemSize);
    }

    for (int i = 0; i < systemSize; i++) {
        for (int j = 0; j < systemSize; j++) {
            if (i == j) {
                matrix[i][j] = (4.0 * h) / 6;  
            } else if (abs(i - j) == 1) {
                matrix[i][j] = h / 6;  
            } else {
                matrix[i][j] = 0.0;
            }
        }
    }
}

void Matrix::printMatrix() {
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << "\n";
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

double modulus(double x) {
    return std::abs(x);
}

void Matrix::printFunction(std::ofstream& output, std::vector<double>& gamma) {
    const int numPoints = 100000;
    double step = (b - a) / numPoints;
    
    std::vector<double> fullGamma(n);
    fullGamma[0] = 0.0; 
    for (int i = 0; i < gamma.size(); i++) {
        fullGamma[i + 1] = gamma[i];
    }
    fullGamma[n - 1] = 0.0; 
    
    for (int i = 0; i <= numPoints; i++) {
        double x = a + i * step;
        
        int j = std::min((int)((x - a) / h), n - 2);
        if (j < 0) j = 0;
        
        double x_j = a + j * h;
        double x_j1 = a + (j + 1) * h;
        
        double term1 = f(x_j) * (x_j1 - x) / h;
        double term2 = f(x_j1) * (x - x_j) / h;
        double term3 = fullGamma[j] * (std::pow(x_j1 - x, 3) - h * h * (x_j1 - x)) / (6.0 * h);
        double term4 = fullGamma[j + 1] * (std::pow(x - x_j, 3) - h * h * (x - x_j)) / (6.0 * h);
        
        double s_x = term1 + term2 + term3 + term4;
        
        output << x << " " << s_x << "\n";
    }
}

int main() {
    Matrix m = Matrix(100, -1, 1);
    m.setFunction(modulus);
    // leave here for debug purpose.
    // m.printMatrix();
    std::vector<double> gamma = m.prognat();
    
    std::ofstream out("out.txt");
    m.printFunction(out, gamma);
    std::cout << "\n";
    return 0;
}
