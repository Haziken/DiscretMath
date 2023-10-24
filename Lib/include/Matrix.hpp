#pragma once

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <ostream>
#include <span>
#include <vector>

namespace dml
{
template<typename Tt = bool,
         size_t Tx = 4,
         size_t Ty = 4,
         size_t Tmin = std::min(Tx, Ty),
         size_t Tmax = std::max(Tx, Ty)>
class Matrix {
public:
    using Point = std::pair<size_t, size_t>;
    // standart constructors;
    Matrix() : m_matrix(std::vector<std::vector<Tt>>(Tx, std::vector<Tt>(Ty, 0))) {
    }
    Matrix(const Matrix& matrix) = default;
    Matrix(Matrix&& matrix) = default;

    Matrix& operator=(const Matrix& matrix) = default;
    Matrix& operator=(Matrix&& matrix) = default;

    // Load data
    Matrix(std::vector<std::tuple<size_t, size_t, Tt>> data)
      : m_matrix(std::vector<std::vector<Tt>>(Tx, std::vector<Tt>(Ty, 0))) {
        for (const auto& [_x, _y, _data] : data) {
            m_matrix[_x][_y] = _data;
        }
    }
    Matrix(std::vector<std::vector<Tt>> data) : m_matrix(data) {
    }

    void set(Point p, Tt value) {
        m_matrix[p.first][p.second] = value;
    }

    void set(const std::vector<Tt>& data) {
        for (size_t i = 0; i < Tx; ++i) {
            for (size_t j = 0; j < Ty; ++j) {
                m_matrix.at(i).at(j) = data.at(j + i * Ty);
            }
        }
    }

    Tt& get(Point p) {
        return m_matrix.at(p.first).at(p.second);
    }

    // Math functions
    bool isReflectivity() {
        bool result = true;
        for (size_t i = 0; i < Tmin; ++i)
            result = result && m_matrix.at(i).at(i);
        return result;
    }
    bool isAntiReflexivity() {
        bool result = true;
        for (size_t i = 0; i < Tmin; ++i)
            result = result && !(m_matrix.at(i).at(i));
        return result;
    }
    bool isSymmetry() {
        bool result = true;
        auto all_points = getAllPoints();
        for (auto& point : all_points)
            result = result && (m_matrix.at(point.second).at(point.first) && m_matrix.at(point.first).at(point.second));
        return result;
    }
    bool isAntiSymmetry() {
        return !isSymmetry();
    }
    bool isAsymmetry() {
        bool result = true;
        auto all_points = getAllPoints();
        for (auto& point : all_points)
            result = result && (m_matrix.at(point.first).at(point.second) != m_matrix.at(point.second).at(point.first));
        return result;
    }
    bool isTransitivity() {
        bool result = true;
        auto all_points = getAllPoints();
        for (size_t a = 0; a < all_points.size(); ++a)
            for (size_t b = 0; b < all_points.size(); ++b) {
                if (a == b)
                    continue;
                Point& p_a = all_points.at(a);
                Point& p_b = all_points.at(b);
                if (p_a.second == p_b.first)
                    result = result && m_matrix.at(p_a.first).at(p_b.second);
            }
        return result;
    }

    Matrix getReflexiveClosureMatrix() {
        Matrix result = *this;
        for (size_t i = 0; i < Tmin; ++i)
            result.set({i, i}, true);
        return result;
    }
    Matrix getSymmetricClosureMatrix() {
        Matrix result = *this;
        for (size_t i = 0; i < Tmin; ++i)
            for (size_t j = i; j < Tmin; ++j)

                if (m_matrix.at(i).at(j) || m_matrix.at(j).at(i)) {
                    result.set({i, j}, true);
                    result.set({j, i}, true);
                }
        return result;
    }
    Matrix getTransitiveClosureMatrix() {
        Matrix last_iter = *this;
        Matrix result;
        do {
            result = last_iter;
            last_iter = result * result;
            last_iter = last_iter + result;
        } while (result != last_iter);
        return result;
    }

    Matrix getTransitiveClosureMatrixWarshall() {
        Matrix result = *this;
        for (size_t i = 0; i < Tx; ++i)
            for (size_t j = 0; j < Ty; ++j) {
                if (i == j)
                    continue;
                if (result.m_matrix.at(i).at(j))
                    for (size_t k = 0; k < Ty; ++k)
                        result.set({i, k}, result.m_matrix.at(i).at(k) || result.m_matrix.at(j).at(k));
            }

        return result;
    }

    size_t getPower() {
        size_t sum = 0;
        for (auto& i : m_matrix)
            for (auto j : i)
                sum += j;
        return sum;
    }

    std::vector<Point> getBinaryRelation();

    std::vector<Point> getAllPoints() {
        std::vector<Point> result;
        for (size_t x = 0; x < Tx; ++x)
            for (size_t y = 0; y < Ty; ++y)
                if (m_matrix.at(x).at(y))
                    result.push_back({x, y});
        return result;
    }

    template<typename _Tt, size_t _Tx, size_t _Ty>
    Matrix<Tt, Ty, _Tx> operator*(Matrix<_Tt, _Tx, _Ty> matrix) {
        Matrix<Tt, Ty, _Tx> result;
        for (size_t i = 0; i < Ty; ++i) {
            for (size_t j = 0; j < _Tx; ++j) {
                Tt sum = 0;
                for (size_t k = 0; k < Tx; ++k) {
                    sum += m_matrix[i][k] * matrix.m_matrix[k][j];
                }
                result.set({i, j}, sum);
            }
        }
        return result;
    }

    Matrix operator+(const Matrix& matrix) {
        Matrix result = *this;
        for (size_t i = 0; i < Tx; ++i)
            for (size_t j = 0; j < Ty; ++j)
                result.set({i, j}, m_matrix.at(i).at(j) + matrix.m_matrix.at(i).at(j));
        return result;
    }

    bool operator==(const Matrix& matrix) {
        bool result = true;
        for (size_t i = 0; i < Tx; ++i)
            for (size_t j = 0; j < Ty; ++j)
                result = result && (m_matrix.at(i).at(j) == matrix.m_matrix.at(i).at(j));
        return result;
    }

    const std::vector<std::vector<Tt>>& getMatrix() const {
        return m_matrix;
    }

private:
    std::vector<std::vector<Tt>> m_matrix;
};

} // namespace dml

template<typename Tt, size_t Tx, size_t Ty>
std::ostream& operator<<(std::ostream& stream, const dml::Matrix<Tt, Tx, Ty>& matrix) {
    for (auto& i : matrix.getMatrix()) {
        for (auto j : i) {
            stream << j << " ";
        }
        stream << std::endl;
    }
    return stream;
}