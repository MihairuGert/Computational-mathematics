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

double dichotomy(Equation& equation, MethodSettings& methodSettings, double a, double b, bool isAInf, bool isBInf) {
    if (isAInf && isBInf) {
        while (equation.getCubicValue(b)*equation.getCubicValue(a) > 0) {
            b += methodSettings.delta;
            a -= methodSettings.delta;
        }
    }
    if (isAInf)
        while (equation.getCubicValue(b)*equation.getCubicValue(a) > 0) {
            b = a;
            a -= methodSettings.delta;
        }
    if (isBInf) {
        while (equation.getCubicValue(b)*equation.getCubicValue(a) > 0) {
            a = b;
            b += methodSettings.delta;
        }
    }
    for (;;) {
        double c = (a + b) / 2;
        double value = equation.getCubicValue(c);
        if (value <= methodSettings.epsilon && value >= -methodSettings.epsilon) {
            return c;
        }
        if(equation.getCubicValue(b) * equation.getCubicValue(c) < 0) {
            a = c;
        } else {
            b = c;
        }
    }
}

int multiplicity_test(Equation& equation, MethodSettings& methodSettings, double root) {
    double first_der = 3*root*root + 2*equation.a*root + equation.b;
    double secnd_der = 6*root + 2*equation.a;
    if (first_der > methodSettings.epsilon || first_der < -methodSettings.epsilon) {
        return 1;
    }
    if (secnd_der > methodSettings.epsilon || secnd_der < -methodSettings.epsilon) {
        return 2;
    }
    return 3;
}

std::vector<double> find_roots(Equation& cubicEquation, MethodSettings& methodSettings) {
    Equation quad_equation = {3, 2*cubicEquation.a, cubicEquation.b};
    std::vector<double> extremum_points = quad_equation.solve_quadratic(methodSettings);
    std::vector<double> roots;

    std::cout << "Roots are: \n";
    if (extremum_points.empty()) {
        double x = dichotomy(cubicEquation, methodSettings, 0, 0, true, true);
        roots.push_back(x);
        return roots;
    }

    if (extremum_points.size() == 1) {
        double x = dichotomy(cubicEquation, methodSettings, extremum_points[0], extremum_points[0], true, true);
        roots.push_back(x);
        return roots;
    }

    if (extremum_points.size() == 2) {
        double value0 = cubicEquation.getCubicValue(extremum_points[0]);
        double value1 = cubicEquation.getCubicValue(extremum_points[1]);
        if (value0 > 0 && value1 > 0) {
            double x = dichotomy(cubicEquation, methodSettings, extremum_points[0], extremum_points[0], true, false);
            roots.push_back(x);
            return roots;
        }

        if (value0 < 0 && value1 < 0) {
            double x = dichotomy(cubicEquation, methodSettings, extremum_points[1], extremum_points[1], false, true);
            roots.push_back(x);
            return roots;
        }

        if (value0 < methodSettings.epsilon &&
                value0 > -methodSettings.epsilon &&
                value1 < -methodSettings.epsilon) {
            double x1 = extremum_points[0];
            double x2 = dichotomy(cubicEquation, methodSettings, extremum_points[1], extremum_points[1], false, true);
            roots.push_back(x1);
            roots.push_back(x2);
            return roots;
        }

        if (value1 < methodSettings.epsilon &&
                value1 > -methodSettings.epsilon &&
                value0 > methodSettings.epsilon) {
            double x1 = extremum_points[1];
            double x2 = dichotomy(cubicEquation, methodSettings, extremum_points[0], extremum_points[0], true, false);
            roots.push_back(x1);
            roots.push_back(x2);
            return roots;
        }

        if (value0 > methodSettings.epsilon &&
                value1 < -methodSettings.epsilon) {
            double x1 = dichotomy(cubicEquation, methodSettings, extremum_points[0], extremum_points[0], true, false);
            double x2 = dichotomy(cubicEquation, methodSettings, extremum_points[0], extremum_points[1], false, false);
            double x3 = dichotomy(cubicEquation, methodSettings, extremum_points[1], extremum_points[1], false, true);
            roots.push_back(x1);
            roots.push_back(x2);
            roots.push_back(x3);
            return roots;
        }

        if (value0 > -methodSettings.epsilon && value0 < methodSettings.epsilon &&
            value1 > -methodSettings.epsilon && value1 < methodSettings.epsilon) {
            double x1 = (extremum_points[0] + extremum_points[1]) / 2;
            roots.push_back(x1);
            return roots;
        }
    }
    return {};
}

int main() {

    double a = 1, b = 1, c = 1, epsilon = 1, delta = 1;
    std::cin >> a >> b >> c >> epsilon >> delta;

    Equation cubicEquation = Equation(a, b, c);

    if (epsilon <= 0) {
        return 1;
    }
    MethodSettings methodSettings = MethodSettings(epsilon, delta);

    std::vector<double> roots = find_roots(cubicEquation, methodSettings);

    for (double root : roots) {
        std::cout << "Root=" << root << " Multiplicity=" << multiplicity_test(cubicEquation, methodSettings, root) << '\n';
    }

    return 0;
}
