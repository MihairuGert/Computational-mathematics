#include <vector>
#include <iostream>
#include <cmath>

class MethodSettings {
public:
    double epsilon;
    double delta;
    MethodSettings(double epsilon, double delta) : epsilon(epsilon), delta(delta) {};
};

class Equation {
public:
    double a;
    double b;
    double c;
    Equation(double a, double b, double c) : a(a), b(b), c(c) {};
    double getCubicValue(double x) const;
    std::vector<double> solve_quadratic(MethodSettings& methodSettings) const;
};

double Equation::getCubicValue(double x) const {
    return x*x*x + a*x*x + b * x + c;
}

std::vector<double> Equation::solve_quadratic(MethodSettings &methodSettings) const {
    double discriminant = b*b - 4 * a * c;
    std::vector<double> roots;

    if (discriminant > methodSettings.epsilon) {
        double x1 = (-b - sqrt(discriminant)) / (2*a);
        double x2 = (-b + sqrt(discriminant)) / (2*a);
        roots.push_back(x1);
        roots.push_back(x2);
        return roots;
    }

    if (discriminant <= methodSettings.epsilon && discriminant >= -methodSettings.epsilon) {
        double x = (-b + sqrt(discriminant)) / (2*a);
        roots.push_back(x);
        return roots;
    }

    return roots;
}

double dichotomy(Equation& equation, MethodSettings& methodSettings, double a, double b, bool isAInf, bool isBinf) {
    for (;;) {
        double c = (a + b) / 2;
        double value = equation.getCubicValue(c);
        if (value <= methodSettings.epsilon && value >= -methodSettings.epsilon) {
            return c;
        }

        if (value < -methodSettings.epsilon) {
            a = c;
            if (isBinf)
                b = b + methodSettings.delta;
        }

        if (value > methodSettings.epsilon) {
            if (isAInf)
                a = a - methodSettings.delta;
            b = c;
        }
    }
}

void find_roots(Equation& cubicEquation, MethodSettings& methodSettings) {
    Equation quad_equation = {3, 2*cubicEquation.a, cubicEquation.b};
    std::vector<double> extremum_points = quad_equation.solve_quadratic(methodSettings);

    if (extremum_points.empty()) {
        std::cout << dichotomy(cubicEquation, methodSettings, 0, 0, true, true);
        return;
    }

    if (extremum_points.size() == 1) {
        std::cout << dichotomy(cubicEquation, methodSettings, extremum_points[0], extremum_points[0], true, true);
        return;
    }

    if (extremum_points.size() == 2) {
        double value0 = cubicEquation.getCubicValue(extremum_points[0]);
        double value1 = cubicEquation.getCubicValue(extremum_points[1]);
        if (value0 > 0 && value1 > 0) {
            std::cout << dichotomy(cubicEquation, methodSettings, extremum_points[0], extremum_points[0], true, false);
            return;
        }

        if (value0 < 0 && value1 < 0) {
            std::cout << dichotomy(cubicEquation, methodSettings, extremum_points[1], extremum_points[1], false, true);
            return;
        }

        if (value0 < methodSettings.epsilon &&
                value0 > -methodSettings.epsilon &&
                value1 < -methodSettings.epsilon) {
            std::cout << extremum_points[0] << " " << dichotomy(cubicEquation, methodSettings, extremum_points[1], extremum_points[1], false, true);
            return;
        }

        if (value1 < methodSettings.epsilon &&
                value1 > -methodSettings.epsilon &&
                value0 > methodSettings.epsilon) {
            std::cout << extremum_points[1] << " " << dichotomy(cubicEquation, methodSettings, extremum_points[0], extremum_points[0], true, false);
            return;
        }

        if (value0 > methodSettings.epsilon &&
                value1 < -methodSettings.epsilon) {
            std::cout << dichotomy(cubicEquation, methodSettings, extremum_points[0], extremum_points[0], true, false)
                    << " " << dichotomy(cubicEquation, methodSettings, extremum_points[0], extremum_points[1], false, false)
                    << " " << dichotomy(cubicEquation, methodSettings, extremum_points[1], extremum_points[1], false, true);
            return;
        }
    }
}

int main() {

    double a = 1, b = 1, c = 1, epsilon = 1, delta = 1;
    std::cin >> a >> b >> c >> epsilon >> delta;

    Equation cubicEquation = Equation(a, b, c);

    if (epsilon <= 0) {
        return 1;
    }
    MethodSettings methodSettings = MethodSettings(epsilon, delta);

    find_roots(cubicEquation, methodSettings);

    return 0;
}
