#include <vector> 
#include <iostream>

namespace zich
{
    class Matrix
    {
        public:

            int rows;
            int cols;
            // 2D Vector
            // Reference https://www.cplusplus.com/reference/vector/vector/vector/
            std::vector<std::vector<double>> myMat;
            Matrix(const std::vector<double>& v, int n, int m);
            Matrix();
            Matrix(const Matrix& other);

            // Operators
            // Reference https://en.cppreference.com/w/cpp/language/operators

            // Output
            friend std::ostream& operator << (std::ostream& out, const Matrix& mat);
            friend std::istream& operator >> (std::istream& out, Matrix& mat);

            friend Matrix operator*(double num, const Matrix& mat);

            // Modifications, A and B are matrices, c is a scalar.
            // used https://www.tutorialspoint.com/cplusplus/unary_operators_overloading.htm


            Matrix& operator+(const Matrix& other);
            Matrix operator+();

            Matrix& operator-(const Matrix& other);
            Matrix& operator-();
            Matrix& operator*(const double num);
            Matrix& operator*(const Matrix& other);
            Matrix& operator=(const Matrix& other);
            Matrix& operator*=(const double num);

            Matrix& operator--();
            Matrix operator--(int);

            Matrix& operator++();
            Matrix operator++(int);

            Matrix& operator+=(const Matrix& other);
            Matrix& operator-=(const Matrix& other);
            Matrix& operator*=(const Matrix& other);


            // Comparison
            friend bool operator==(const Matrix& m1, const Matrix& m2);
            friend bool operator!=(const Matrix& m1, const Matrix& m2);
            friend bool operator>(const Matrix& m1, const Matrix& m2);
            friend bool operator>=(const Matrix& m1, const Matrix& m2);
            friend bool operator<=(const Matrix& m1, const Matrix& m2);
            friend bool operator<(const Matrix& m1, const Matrix& m2);
            
        
    };
}