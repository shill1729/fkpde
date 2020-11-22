#ifndef TRIDIAG_H
#define TRIDIAG_H
#include <vector>

class Tridiag{
public:

  Tridiag();
  virtual ~Tridiag();

  // Solve tridiagonal linear systems
  static std::vector<double> thomasAlgorithm(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> d);

  // For constant tridiagonals
  static std::vector<double> thomasAlgorithm(double a, double b, double c, std::vector<double> d);
};

#endif
