#include "s21_matrix_oop.h"

/*******************************************************************************
 * must-have functions for any class
 ******************************************************************************/

S21Matrix::S21Matrix() noexcept : rows_(0), cols_(0), matrix_(nullptr) {}

S21Matrix::S21Matrix(int rows, int cols) {
  rows_ = rows;
  cols_ = cols;
  if (rows_ < 0 || cols_ < 0)
    throw std::out_of_range("Error: Wrong matrix size");
  matrix_ = nullptr;
  CreateMatrix();
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = 0;
    }
  }
}

S21Matrix::S21Matrix(const S21Matrix &other) {
  // S21Matrix b(a) -> создать объект b по образу и подобию объекта а
  matrix_ = nullptr;
  this->rows_ = other.rows_;
  this->cols_ = other.cols_;
  CreateMatrix();
  for (int row = 0; row < rows_; row++) {
    std::memcpy(matrix_[row], other.matrix_[row], cols_ * sizeof(double));
  }
}

S21Matrix::S21Matrix(S21Matrix &&other) noexcept
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.matrix_ = nullptr;
  other.rows_ = 0;
  other.cols_ = 0;
}

S21Matrix &S21Matrix::operator=(const S21Matrix &other) {
  S21Matrix copy{other};
  // Ниже вызовется конструктор присваивания перемещением
  // std::move() преобразовывает переданное значение в R-value
  *this = std::move(copy);
  return *this;
}

S21Matrix &S21Matrix::operator=(S21Matrix &&other) noexcept {
  if (this != &other) {
    RemoveMatrix();
    std::swap(rows_, other.rows_);
    std::swap(cols_, other.cols_);
    std::swap(matrix_, other.matrix_);
  }
  return *this;
}

S21Matrix::~S21Matrix() noexcept { RemoveMatrix(); }

/*******************************************************************************
 * Task functions
 ******************************************************************************/

bool S21Matrix::IsEqualSize(const S21Matrix &other) const {
  return (rows_ == other.rows_ && cols_ == other.cols_);
}

bool S21Matrix::EqMatrix(const S21Matrix &other) const {
  bool res = true;
  if (IsEqualSize(other)) {
    for (int i = 0; i < this->rows_; i++) {
      for (int j = 0; j < this->cols_; j++) {
        if (abs(this->matrix_[i][j] - other.matrix_[i][j]) > kEpsilon) {
          res = false;
        }
      }
    }
  } else {
    res = false;
  }
  return res;
}

void S21Matrix::SumMatrix(const S21Matrix &other) {
  if (!IsEqualSize(other)) {
    throw std::range_error("Error: Matrices must have the same size");
  }
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < this->cols_; j++) {
      this->matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  if (!IsEqualSize(other)) {
    throw std::range_error("Error: Matrices must have the same size");
  }
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < this->cols_; j++) {
      this->matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < this->cols_; j++) {
      this->matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix &other) {
  if (this->cols_ != other.rows_)
    throw std::range_error(
        "Error: the number of columns of the first "
        "matrix is not equal to the number of rows of the second matrix");
  S21Matrix result(rows_, other.cols_);
  for (int row = 0; row < result.rows_; row++) {
    for (int col = 0; col < result.cols_; col++) {
      for (int row_t = 0; row_t < cols_; row_t++) {
        result.matrix_[row][col] +=
            matrix_[row][row_t] * other.matrix_[row_t][col];
      }
    }
  }
  // Используем перемещение вместо копирования (тк матрица result нам больше не
  // нужна)
  *this = std::move(result);
}

S21Matrix S21Matrix::Transpose() const {
  S21Matrix result(cols_, rows_);
  for (int row = 0; row < rows_; row++) {
    for (int col = 0; col < cols_; col++) {
      result.matrix_[col][row] = matrix_[row][col];
    }
  }
  return result;
}

bool S21Matrix::IsSquare() const {
  if (rows_ == cols_) {
    return true;
  } else {
    return false;
  }
}

S21Matrix S21Matrix::GetMinorMatrix(int rows, int cols) const {
  S21Matrix minor(rows_ - 1, cols_ - 1);
  for (int row = 0, row_n = 0; row < rows_; row++) {
    if (row == rows) continue;
    for (int col = 0, col_n = 0; col < cols_; col++) {
      if (col == cols) continue;
      minor.matrix_[row_n][col_n] = matrix_[row][col];
      col_n++;
    }
    row_n++;
  }
  return minor;
}

double S21Matrix::Determinant() const {
  if (rows_ != cols_) {
    throw std::range_error("Error: matrix is not square");
  }
  double det = 1.0;
  S21Matrix tmp{*this};

  for (int i = 0; i < rows_; i++) {
    // Выбираем наибольший элемент в столбце. max = индекс макс эл-та.
    int max = i;
    for (int j = i + 1; j < rows_; j++) {
      if (abs(tmp.matrix_[j][i]) > abs(tmp.matrix_[max][i])) {
        max = j;
      }
    }
    // Если получили в качестве макс значения ноль, то det сразу = 0
    if (abs(tmp.matrix_[max][i]) < kEpsilon) {
      return 0.0;
    } else {
      tmp.SwapRows(i, max);
      // Текущая строка уже не будет меняться => можем сразу умножить det на
      // элемент (i, i) матрицы
      det *= tmp.matrix_[i][i];
      // Если был обмен строк, то меняем знак определителя
      if (i != max) {
        det *= -1;
      }
      // Осуществляем вычитание текущей строки из всех следующих
      for (int j = i + 1; j < rows_; j++) {
        // Считаем коэффициент, чтобы занулить элементы строки
        double koef = tmp(j, i) / tmp(i, i);
        for (int k = i; k < rows_; k++) {
          // Зануляем, начиная со столбца i, т.к. предыдущие столбцы уже
          // занулены на предыдущих шагах и не изменятся
          tmp.matrix_[j][k] -= tmp.matrix_[i][k] * koef;
        }
      }
    }
  }
  return det;
}

void S21Matrix::SwapRows(int row1, int row2) noexcept {
  if (row1 != row2) {
    for (int i = 0; i < cols_; i++) {
      std::swap((*this)(row1, i), (*this)(row2, i));
    }
  }
}

S21Matrix S21Matrix::CalcComplements() const {
  if (!IsSquare() || rows_ < 2)
    throw std::range_error("Error: matrix is not square");
  S21Matrix result(rows_, cols_);
  for (int row = 0; row < rows_; row++) {
    for (int col = 0; col < cols_; col++) {
      double det = GetMinorMatrix(row, col).Determinant();
      result.matrix_[row][col] = det * pow(-1, row + col);
    }
  }
  return result;
}

S21Matrix S21Matrix::InverseMatrix() const {
  if (rows_ != cols_) {
    throw std::logic_error("Incorrect matrix size for Inverse");
  }
  double det = Determinant();
  if (fabs(det) < kEpsilon) {
    throw std::range_error("Error: matrix Determinant is zero");
  }
  return Transpose().CalcComplements() * (1.0 / det);
}

/*******************************************************************************
 * getters, setters, operators
 ******************************************************************************/
int S21Matrix::get_rows() const { return rows_; }

int S21Matrix::get_cols() const { return cols_; }

void S21Matrix::set_rows(int new_rows) {
  if (new_rows < 0) {
    throw std::length_error("Matrix rows count must be non-negative");
  }
  if (new_rows != rows_) {
    // Работа с копией - гарантия безопасности исключений
    S21Matrix tmp{new_rows, cols_};
    int min = std::min(rows_, new_rows);
    for (int i = 0; i < min; i++) {
      for (int j = 0; j < cols_; j++) {
        tmp.matrix_[i][j] = this->matrix_[i][j];
      }
    }
    // Используем перемещение вместо копирования (тк матрица tmp нам больше не
    // нужна)
    *this = std::move(tmp);
  }
}

void S21Matrix::set_cols(int new_cols) {
  if (new_cols < 0) {
    throw std::length_error("Matrix cols count must be non-negative");
  }
  if (new_cols != cols_) {
    // Работа с копией - гарантия безопасности исключений
    S21Matrix tmp{rows_, new_cols};
    int min = std::min(cols_, new_cols);
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < min; j++) {
        tmp.matrix_[i][j] = this->matrix_[i][j];
      }
    }
    // Используем перемещение вместо копирования (тк матрица tmp нам больше не
    // нужна)
    *this = std::move(tmp);
  }
}

// Индексация по элементам матрицы (строка, колонка)
double &S21Matrix::operator()(int row, int col) {
  if (row < 0 || col < 0 || row >= rows_ || col >= cols_)
    throw std::out_of_range("Error: index is out of range");
  return matrix_[row][col];
}

// Индексация по элементам матрицы (строка, колонка) для const объектов
const double &S21Matrix::operator()(int row, int col) const {
  if (row < 0 || col < 0 || row >= rows_ || col >= cols_)
    throw std::out_of_range("Error: index is out of range");
  return matrix_[row][col];
}

bool S21Matrix::operator==(const S21Matrix &other) const {
  return EqMatrix(other);
}

bool S21Matrix::operator!=(const S21Matrix &other) const {
  return !EqMatrix(other);
}

S21Matrix S21Matrix::operator+(const S21Matrix &other) const {
  S21Matrix result(*this);
  result.SumMatrix(other);
  return result;
}

S21Matrix &S21Matrix::operator+=(const S21Matrix &other) {
  SumMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator-(const S21Matrix &other) const {
  S21Matrix result(*this);
  result.SumMatrix(other);
  return result;
}

S21Matrix &S21Matrix::operator-=(const S21Matrix &other) {
  SubMatrix(other);
  return *this;
}

// Умножение матрицы на число
S21Matrix S21Matrix::operator*(double num) const {
  S21Matrix result(*this);
  result.MulNumber(num);
  return result;
}

// Дружественная функция (делаем так для переменной number)
S21Matrix operator*(double number, const S21Matrix &matrix) noexcept {
  S21Matrix tmp = matrix * number;
  return tmp;
}

S21Matrix S21Matrix::operator*(const S21Matrix &other) const {
  S21Matrix result(*this);
  result.MulMatrix(other);
  return result;
}

S21Matrix &S21Matrix::operator*=(const S21Matrix &other) {
  MulMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator*=(double num) {
  MulNumber(num);
  return *this;
}

void S21Matrix::CreateMatrix() {
  matrix_ = new double *[rows_];
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_];
  }
}

void S21Matrix::RemoveMatrix() {
  if (matrix_ != nullptr) {
    for (int i = 0; i < rows_; i++) {
      if (matrix_[i] != nullptr) {
        delete[] matrix_[i];
        matrix_[i] = nullptr;
      }
    }
    delete[] matrix_;
    matrix_ = nullptr;
  }
}