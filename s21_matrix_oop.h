#ifndef SRC_S21_MATRIX_OOP_H
#define SRC_S21_MATRIX_OOP_H

#include <algorithm>  // std::copy | std::min
#include <cmath>      // std::abs
#include <cstring>
#include <iostream>
#include <stdexcept>  // length_error | out_of_range | logic_error
#include <utility>    // std::move | std::swap

using namespace std;

class S21Matrix final {
 private:
  int rows_, cols_;  // Rows and columns
  double **matrix_;  // Pointer to the memory where the matrix is allocated
  const double kEpsilon = 1e-7;  // Constant to compare with zero
 public:
  S21Matrix() noexcept;                    // Default constructor
  explicit S21Matrix(int rows, int cols);  // Parameterized Constructor
  S21Matrix(const S21Matrix &other);       // Copy constructor
  S21Matrix(S21Matrix &&other) noexcept;   // Move constructor
  S21Matrix &operator=(const S21Matrix &other);  // Присваивание копированием
  S21Matrix &operator=(
      S21Matrix &&other) noexcept;  // Присваивание перемещением
  ~S21Matrix() noexcept;            // Деструктор

  int get_rows() const;  // A read-only function
  int get_cols() const;  // A read-only function
  void set_rows(int new_rows);  // Изменение количества строк матрицы
  void set_cols(int new_cols);  // Изменение количества столбцов матрицы

  // Basic task functions we should to create
  bool EqMatrix(const S21Matrix &other) const;
  void SumMatrix(const S21Matrix &other);
  void SubMatrix(const S21Matrix &other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix &other);
  S21Matrix Transpose() const;
  S21Matrix CalcComplements() const;
  S21Matrix GetMinorMatrix(int rows, int cols) const;
  double Determinant() const;
  S21Matrix InverseMatrix() const;

  // Operators
  double &operator()(int row, int col);
  const double &operator()(int row, int col) const;
  bool operator==(const S21Matrix &other) const;
  bool operator!=(const S21Matrix &other) const;
  S21Matrix operator+(const S21Matrix &other) const;
  S21Matrix &operator+=(const S21Matrix &other);
  S21Matrix operator-(const S21Matrix &other) const;
  S21Matrix &operator-=(const S21Matrix &other);
  S21Matrix operator*(double num) const;
  friend S21Matrix operator*(double number, const S21Matrix &matrix) noexcept;
  S21Matrix operator*(const S21Matrix &other) const;
  S21Matrix &operator*=(const S21Matrix &other);
  S21Matrix &operator*=(double num);

  // Other methods..
  bool IsEqualSize(const S21Matrix &other) const;
  bool IsSquare() const;
  void CreateMatrix();
  void RemoveMatrix();
  void SwapRows(int row1, int row2) noexcept;
};

#endif  // SRC_S21_MATRIX_OOP_H
