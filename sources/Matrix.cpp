#include "Matrix.hpp"
#include <stdexcept>
#include <string>
#include <iostream>


namespace zich
{
    bool notCompatible(const Matrix& m1, const Matrix& m2)
    {
        return m1.rows != m2.rows || m1.cols != m2.cols;
    }

    bool notCompatibleMult(const Matrix& m1, const Matrix& m2)
    {
        return m1.cols != m2.rows;
    }

    double multOfVectors(const std::vector<double>v1, const std::vector<double>v2)
    {
        if (v1.size() != v2.size())
            {
                throw std::invalid_argument("Error while trying to sum two vectors");
            }
        double ans = 0;
        for (size_t i = 0; i < v1.size(); i++)
        {
            ans += v1.at(i) * v2.at(i);
        }
        return ans;
    }
    
    Matrix::Matrix(const std::vector<double>& v, int n, int m)
    {
        if (v.size() != n * m)
            {
                throw std::invalid_argument("Different row \\ col number than given vector.");
            }
            
        if (n < 0 || m < 0) 
            throw std::invalid_argument("Negative matrix size");

        this->rows = n;
        this->cols = m;

        int i = 0;
        std::vector<double> tempCol;
        for (double curr : v)
        {
            tempCol.push_back(curr);
            i++;

            if (i == m)
            {
                i = 0;
                myMat.push_back(tempCol);
                tempCol.clear();
            }
        }
    }

    Matrix::Matrix()
    {
        
    }
    
    Matrix::Matrix(const Matrix& other)
    {
        this->cols = other.cols;
        this->rows = other.rows;
        for (std::vector<double> r : other.myMat)
        {
            std::vector<double> v;
            for (double d : r)
            {
                v.push_back(d);
            }
            myMat.push_back(v);
        }
    }

    std::ostream& operator << (std::ostream& out, const Matrix& mat)
    {
        std::string ans;
        for (size_t i = 0; i < mat.rows; i++)
        {
            ans += "[";
            std::vector<double> my_row = mat.myMat.at(i);
            for (size_t j = 0; j < mat.cols; j++)
            {
                if (j != 0)
                    ans += " ";
                ans += std::to_string(my_row.at(j));
            }

            ans += "]\n";

        }
        return out << ans ;
    }


    std::istream& operator >> (std::istream& out, Matrix& mat)
    {
        return out;
    }

    Matrix operator * (double num, const Matrix& mat)
    {
        Matrix m;

        for (size_t i = 0; i < mat.rows; i++)
        {
            std::vector<double> m_row;
            std::vector<double> other_row = mat.myMat.at(i);

            for (size_t j = 0; j < mat.cols; j++)
            {
                m_row.push_back(other_row.at(j) * num);
            }
            m.myMat.push_back(m_row);
        }
        return m;
    }

    // Operators
    // Reference https://en.cppreference.com/w/cpp/language/operators

    // Output

    // Matrix& zich::Matrix::operator*(double num, /*const*/ Matrix& mat)
    // {
    // return *this;
    // }

    Matrix& Matrix::operator+(const Matrix& other)
    {
        if (notCompatible(*this, other))
            {
                throw std::invalid_argument("Different row \\ col number.");
            }
        Matrix& m(*this);
        for (size_t i = 0; i < other.rows; i++)
        {
            for (size_t j = 0; j < other.cols; j++)
            {
                m.myMat.at(i).at(j) += other.myMat.at(i).at(j);
            }
        }
        return m;
    }
    Matrix Matrix::operator+()
    {
        Matrix old = *this;
        operator++();
        return old;
    }
    Matrix& Matrix::operator-( const Matrix& other)
    {
        if (notCompatible(*this, other))
            throw std::invalid_argument("Different row \\ col number.");

        Matrix& m(*this);
        for (size_t i = 0; i < other.rows; i++)
        {
            for (size_t j = 0; j < other.cols; j++)
            {
                m.myMat.at(i).at(j) -= other.myMat.at(i).at(j);
            }
        }
    return m;
    }
    Matrix& Matrix::operator-()
    {
    for (size_t i = 0; i < this->rows; i++)
        {
            for (size_t j = 0; j < this->cols; j++)
            {
                this->myMat.at(i).at(j) *= -1;
            }
        }
    return *this;
    }
    
    Matrix& Matrix::operator++()
    {
    for (size_t i = 0; i < this->rows; i++)
        {
            for (size_t j = 0; j < this->cols; j++)
            {
                this->myMat.at(i).at(j)++; // At this i column, at this j row, increment value.
            }
        }
    return *this;
    }
    
    Matrix Matrix::operator++(int)
    {
        Matrix old = *this;
        operator++();
        return old;
    }

    Matrix& Matrix::operator--()
    {
    for (size_t i = 0; i < this->rows; i++)
        {
            for (size_t j = 0; j < this->cols; j++)
            {
                this->myMat.at(i).at(j)--; // At this i column, at this j row, decrement value.
            }
        }
    return *this;
    }

    Matrix Matrix::operator--(int)
    {
        Matrix old = *this;
        operator--();
        return old;
    }

    Matrix& Matrix::operator*(double num)
    {
        Matrix& m(*this);
        for (size_t i = 0; i < this->rows; i++)
            {
                for (size_t j = 0; j < this->cols; j++)
                {
                    this->myMat.at(i).at(j) *= num;
                }
            }
    return m;
    }

    Matrix& Matrix::operator*(const Matrix& other)
    {
        if (notCompatibleMult(*this, other))
            {
                throw std::invalid_argument("Different row \\ col number.");
            }

        // Mult of m x n by n x k is a new matrix of m x k 
        size_t r = (size_t) this->rows;
        size_t c = (size_t) other.cols;        
        std::vector<double> v;

        for (size_t i = 0; i < c * r; i++)
        {
            std::vector<double> v1;
            std::vector<double> v2;
            for (size_t k = 0; k < other.rows; k++)
            {
                v1.push_back(this->myMat.at(size_t (i / c)).at(k));     // Line vector of this matrix
                v2.push_back(other.myMat.at(k).at(size_t (i % c)));     // Row vector of other matrix
            }
            double ans = multOfVectors (v1, v2);
            v.push_back(ans);
        }
        
        Matrix* m = new Matrix(v, c, r);
        return *m;
    }

    Matrix& Matrix::operator=(const Matrix& other)
    {
        if (notCompatible(*this, other))
            {
                throw std::invalid_argument("Different row \\ col number.");
            }

        for (size_t i = 0; i < this->rows; i++)
            {
                for (size_t j = 0; j < this->cols; j++)
                {
                    this->myMat.at(i).at(j) = other.myMat.at(i).at(j);
                }
            }
        return *this;
    }
    Matrix& Matrix::operator*=(const double num)
    {
    for (size_t i = 0; i < this->rows; i++)
        {
            for (size_t j = 0; j < this->cols; j++)
            {
                this->myMat.at(i).at(j) *= num;
            }
        }
        return *this;
    }

    Matrix& Matrix::operator*=(const Matrix& other)
    {
        if (notCompatibleMult(*this, other))
            throw std::invalid_argument("Different row \\ col number.");
        return *this;
    }

    Matrix& Matrix::operator+=(const Matrix& other)
    {
        if (notCompatible(*this, other))
            throw std::invalid_argument("Different row \\ col number.");
        for (size_t i = 0; i < other.rows; i++)
        {
            std::vector<double> my_row = this->myMat.at(i);
            std::vector<double> other_row = other.myMat.at(i);

            for (size_t j = 0; j < other.cols; j++)
            {
                my_row.at(j) += other_row.at(j);
            }
        }
        return *this;
    }

    Matrix& Matrix::operator-=(const Matrix& other)
    {
        if (notCompatible(*this, other))
            {
                throw std::invalid_argument("Different row \\ col number.");
            }
        for (size_t i = 0; i < other.rows; i++)
        {
            std::vector<double> my_row = this->myMat.at(i);
            std::vector<double> other_row = other.myMat.at(i);

            for (size_t j = 0; j < other.cols; j++)
            {
                my_row.at(j) -= other_row.at(j);
            }
        }
        return *this;
    }

    // Matrix& operator*(const double num, const Matrix& other)
    // {

    // }


    // // Comparison
    bool operator==(const Matrix& m1, const Matrix& m2)
    {
        if (notCompatible(m1, m2))
            {
                throw std::invalid_argument("Different row \\ col number.");
            }
        for (size_t i = 0; i < m1.rows; i++)
        {
            for (size_t j = 0; j < m1.cols; j++)
            {
                if ( m1.myMat.at(i).at(j) != m2.myMat.at(i).at(j) )
                    {
                        return false;
                    }
            }
        }
        return true;
    }
    bool operator!=(const Matrix& m1, const Matrix& m2)
    {
        if (notCompatible(m1, m2))
            throw std::invalid_argument("Different row \\ col number.");

        return !(operator==(m1,m2));
    }
    bool operator>(const Matrix& m1, const Matrix& m2)
    {
        if (notCompatible(m1, m2))
            throw std::invalid_argument("Different row \\ col number.");

        double sum1 = 0;
        double sum2 = 0;
        for (size_t i = 0; i < m1.rows; i++)
        {
            std::vector<double> my_row = m1.myMat.at(i);
            std::vector<double> other_row = m2.myMat.at(i);

            for (size_t j = 0; j < m1.cols; j++)
            {
                sum1 += my_row.at(j); 
                sum2 += other_row.at(j);
            }
        }
        return sum1 > sum2;
    }
    bool operator>=(const Matrix& m1, const Matrix& m2)
    {
        if (notCompatible(m1, m2))
            throw std::invalid_argument("Different row \\ col number.");

        return operator==(m1,m2) || operator>(m1,m2);
    }

    bool operator<(const Matrix& m1, const Matrix& m2)
    {
        if (notCompatible(m1, m2))
            throw std::invalid_argument("Different row \\ col number.");
        double sum1 = 0;
        double sum2 = 0;
        for (size_t i = 0; i < m1.rows; i++)
        {
            std::vector<double> my_row = m1.myMat.at(i);
            std::vector<double> other_row = m2.myMat.at(i);

            for (size_t j = 0; j < m1.cols; j++)
            {
                sum1 += my_row.at(j); 
                sum2 += other_row.at(j);
            }
        }
        return sum1 < sum2;
    }
    bool operator<=(const Matrix& m1, const Matrix& m2)
    {
        if (notCompatible(m1, m2))
            throw std::invalid_argument("Different row \\ col number.");
        return operator==(m1,m2) || operator<(m1,m2);

    }
}