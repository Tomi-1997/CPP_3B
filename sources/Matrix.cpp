#include "Matrix.hpp"
#include <stdexcept>
#include <string>
#include <iostream>
#include <ctype.h> // for isDigit() for input


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

    double multOfVectors(const std::vector<double>& v1, const std::vector<double>& v2)
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
            {
                throw std::invalid_argument("Negative matrix size");
            }

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
    
    Matrix::Matrix(const Matrix& other)
    {
        this->cols = other.cols;
        this->rows = other.rows;
        for (const std::vector<double>& r : other.myMat)
        {
            std::vector<double> v;
            v.reserve(size_t(other.cols));
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
                    {
                        ans += " ";
                    }
                ans += std::to_string( (int) my_row.at(j));
            }

            ans += "]";
            if (i != mat.rows - 1)
            {
                ans += "\n";
            }

        }
        return out << ans ;
    }

    /**
     * @brief Given index, return from string a double up to next ']' or whitespace.
     * 
     * @param in 
     * @param ind 
     * @return string 
     */
    std::string getNum(const std::string& in, int ind)
    {
        std::string ans_str;
        for (size_t i = (size_t)ind; i < in.length(); i++)
        {
            if (in[i] == ' ' || in[i] == ']')
            {
                break;
            }
            if(in[i] == ',')
            {
                throw std::invalid_argument("Invalid input of matrix");
            }
            ans_str += in[i];
        }

        return ans_str;
    }
    std::istream& operator >> (std::istream& in, Matrix& mat)
    {
        int col = 0;
        int row = 0;
        int new_row = 0;
        std::string input(std::istreambuf_iterator<char>(in), {}); // Load input to string
        std::vector<double> v;

        bool input_start_mat = true;
        bool input_nums = false;
        bool input_single_num = false;
        bool input_end_mat = false;

        for (size_t i = 0; i < input.length() - 1; i++)
        {
            char c = input[i];
            if (c == '[' && input_start_mat)
            {
                input_start_mat = false;
                input_nums = true;

                std::string num = getNum(input, (int)i + 1);
                i = i + num.length();

                v.push_back(std::stod(num));
                new_row++;
                continue;
            }

            if (c == ' ' && input_nums)
            {
                input_single_num = true;
                continue;
            }

            if (input_single_num && isdigit(c) != 0)
            {
                std::string num = getNum(input, i);
                i = i + num.length() - 1;
                v.push_back(std::stod(num));
                new_row++;
                continue;                
            }

            if (c == ']' && input_nums)
            {
                input_nums = false;
                input_end_mat = true;
                input_single_num = false;
                col++;
                continue;
            }

            if (input_end_mat && c == ',' && input[i+1] == ' ')
            {
                input_end_mat = false;
                input_start_mat = true;
                new_row++;
                i++;
                continue;
            }

            if (c == '\\' && input[i+1] == 'n')
            {
                break;
            }

            // Different amount of numbers in a different row
            if (row == 0 && new_row != 0 && new_row == row )
            {
                row = new_row;
                new_row = 0;
                continue;
            }

            throw std::invalid_argument("Invalid input of matrix");          
        }

        return in;
    }

    Matrix operator * (double num, const Matrix& mat)
    {
        Matrix m(mat);
        for (size_t i = 0; i < m.rows; i++)
        {
            for (size_t j = 0; j < m.cols; j++)
            {
                m.myMat.at(i).at(j) *= num;
            }
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

    Matrix Matrix::operator+(const Matrix& other)
    {
        if (notCompatible(*this, other))
            {
                throw std::invalid_argument("Different row \\ col number.");
            }
        Matrix m(*this);
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
        // operator++();
        return old;
    }
    Matrix Matrix::operator-( const Matrix& other)
    {
        if (notCompatible(*this, other))
{            throw std::invalid_argument("Different row \\ col number.");
}
        Matrix m(*this);
        for (size_t i = 0; i < other.rows; i++)
        {
            for (size_t j = 0; j < other.cols; j++)
            {
                m.myMat.at(i).at(j) -= other.myMat.at(i).at(j);
            }
        }
    return m;
    }
    Matrix Matrix::operator-()
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

    Matrix Matrix::operator*(double num)
    {
        Matrix m(*this);
        for (size_t i = 0; i < this->rows; i++)
            {
                for (size_t j = 0; j < this->cols; j++)
                {
                    m.myMat.at(i).at(j) *= num;
                }
            }
        return m;
    }

    Matrix Matrix::operator*(const Matrix& other)
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
        
        Matrix m = Matrix(v, r, c);
        return m;
    }

    // Matrix& Matrix::operator=(const Matrix& other)
    // {
    //     this->rows = other.rows;
    //     this->cols = other.cols;

    //     std::vector<std::vector<double>> newv;
    //     for (size_t i = 0; i < other.rows; i++)
    //         {
    //             std::vector<double> v;
    //             for (size_t j = 0; j < other.cols; j++)
    //             {
    //                 v.push_back(other.myMat.at(i).at(j));
    //             }
    //             newv.push_back(v);
    //         }

    //     this->myMat = newv;
    //     return *this;
    // }
    Matrix Matrix::operator*=(const double num)
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

    Matrix Matrix::operator*=(const Matrix& other)
    {
        if (notCompatibleMult(*this, other))
            {
                throw std::invalid_argument("Different row \\ col number at A*=B.");
            }

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
        
        Matrix m = Matrix(v, r, c);
        this->myMat = m.myMat;
        this->rows = r;
        this->cols = c;
        return m;
    }

    Matrix Matrix::operator+=(const Matrix& other)
    {
        if (notCompatible(*this, other))
            {
                throw std::invalid_argument("Different row \\ col number.");
            }
        for (size_t i = 0; i < other.rows; i++)
        {
            for (size_t j = 0; j < other.cols; j++)
            {
                this->myMat.at(i).at(j) += other.myMat.at(i).at(j);
            }
        }
        return *this;
    }

    Matrix Matrix::operator-=(const Matrix& other)
    {
        if (notCompatible(*this, other))
            {
                throw std::invalid_argument("Different row \\ col number.");
            }
        for (size_t i = 0; i < other.rows; i++)
        {
            for (size_t j = 0; j < other.cols; j++)
            {
                this->myMat.at(i).at(j) -= other.myMat.at(i).at(j);
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
                        //std::cout << "FALSE MAN";
                        return false;
                    }
            }
        }
        return true;
    }
    bool operator!=(const Matrix& m1, const Matrix& m2)
    {
        if (notCompatible(m1, m2))
{            throw std::invalid_argument("Different row \\ col number.");
}
        
        // std::cout << m1 << std::endl;
        // std::cout << m2 << std::endl;
        return !(m1 == m2);
    }
    bool operator>(const Matrix& m1, const Matrix& m2)
    {
        if (notCompatible(m1, m2))
{            throw std::invalid_argument("Different row \\ col number.");
}
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
{            throw std::invalid_argument("Different row \\ col number.");
}
        return operator==(m1,m2) || operator>(m1,m2);
    }

    bool operator<(const Matrix& m1, const Matrix& m2)
    {
        if (notCompatible(m1, m2))
{            throw std::invalid_argument("Different row \\ col number.");
}        double sum1 = 0;
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
            {
                throw std::invalid_argument("Different row \\ col number.");
            }
        return operator==(m1,m2) || operator<(m1,m2);

    }
}