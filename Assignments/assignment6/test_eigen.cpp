#include <eigen3/Eigen/Eigenvalues> // header file
#include <iostream>

int main(){
    Eigen::Matrix<double, 2, 2> A; // declare a real (double) 2x2 matrix
    A << 0, 2, 1, 0; // defined the matrix A

    Eigen::EigenSolver<Eigen::Matrix<double, 2,2> > s(A); // the instance s(A) includes the eigensystem
    std::cout << A << std::endl;
    std::cout << "eigenvalues:" << std::endl;
    std::cout << s.eigenvalues() << std::endl;
    std::cout << "eigenvectors=" << std::endl;
    std::cout << s.eigenvectors() << std::endl;

    return(0);
}