//Copyright 2020 Robert

#ifndef INCLUDE_MATRIX_HPP_
#define INCLUDE_MATRIX_HPP_

#include <stdio.h>
#include <math.h>

template<class T>
class Matrix{
 T **p;
 int rows, cols;
public:
 ~Matrix();
 Matrix(int rows, int cols);
 Matrix(const Matrix& copy);
 Matrix Inverse();
 int Rows() const;
 int Cols() const;
 Matrix& operator =(Matrix &m2);
 Matrix operator +(Matrix &m2);
 Matrix operator -(Matrix &m2);
 Matrix operator *(Matrix &m2);
 template<class V>
 friend bool operator ==(const Matrix<V> &m1, const Matrix<V> &m2);
 template<class V>
 friend bool operator !=(const Matrix<V> &m1, const Matrix<V> &m2);
 T* operator [](size_t i) const;
 double determinant(Matrix mat);
 Matrix deleteRowsAndCols(Matrix mat, int nRow, int nCol);
};

template<class T>
Matrix<T>::~Matrix(){
    for (int i = 0; i < rows; i++){
        free(p[i]);
    }
    free(p);
}

template<class T>
int Matrix<T>::Rows() const {
 return rows;
}

template<class T>
int Matrix<T>::Cols() const {
 return cols;
}

template<class T>
Matrix<T>::Matrix(int rows, int cols){
 this -> p = reinterpret_cast<T**>(malloc(rows * sizeof(T *)));
 this -> rows = rows;
 this -> cols = cols;
 for (int i = 0 ; i < rows ; i++){
  p[i] = reinterpret_cast<T*>(malloc(cols * sizeof(T)));
  for (int j = 0 ; j < cols ; j++){
   p[i][j] = 0;
  }
 }
}

template<class T>
Matrix<T>::Matrix(const Matrix& copy){
 cols = copy.Cols();
 rows = copy.Rows();
 this -> p = reinterpret_cast<T**>(malloc(rows * sizeof(T*)));
 for (int i = 0 ; i < copy.Rows() ; i++){
  p[i] = reinterpret_cast<T*>(malloc(cols * sizeof(T)));
  for (int j = 0 ; j < copy.Cols() ; j++){
   p[i][j] = copy[i][j];
  }
 }
}

template<class T>
T* Matrix<T>::operator [](size_t i) const{
 return p[i];
}

template<class T>
Matrix<T>& Matrix<T>::operator =(Matrix &m2){
 this->rows = m2.Rows();
 this->cols = m2.Cols();
 for (int i = 0 ; i < this->rows ; i++){
  for (int j = 0 ; j < this->Cols() ; j++){
   this->p[i][j] = m2[i][j];
  }
 }
 return *this;
}

template <class T>
Matrix<T> Matrix<T>::operator +(Matrix& m2){
 if (this->Cols() != m2.Cols() || this->Rows() != m2.Rows()){
  Matrix<T> res(0, 0);
  return res;
 }
 Matrix<T> res(rows, cols);
 for (int i = 0; i < rows ; i++){
  for (int j = 0 ; j < cols; j++){
   res[i][j] = res[i][j] + m2[i][j];
  }
 }
 return res;
}

template <class T>
Matrix<T> Matrix<T>::operator -(Matrix& m2){
 if (this->Cols() != m2.Cols() || this->Rows() != m2.Rows()){
  Matrix<T> res(0, 0);
  return res;
 }
 Matrix<T> res(rows, cols);
 for (int i = 0; i < rows ; i++){
  for (int j = 0 ; j < cols; j++){
   res[i][j] = (*this)[i][j] - m2[i][j];
  }
 }
 return res;
}

template <class T>
Matrix<T> Matrix<T>::operator *(Matrix& m2){
 if (this->Cols() != m2.Rows()){
  Matrix<T> res(0, 0);
  return res;
 }
 Matrix<T> res(rows, m2.Cols());
 for (int i = 0 ; i < res.Rows(); i++){
  for (int j = 0 ; j < res.Cols(); j++){
   for (int k = 0 ; k < cols ; k++){
    res[i][j] += (*this)[i][k] * m2[k][j];
   }
  }
 }
 return res;
}

template<class T>
Matrix<T> Matrix<T>::deleteRowsAndCols(Matrix<T> mat, int nRow, int nCol){
 Matrix<T> res(mat.Rows() - 1, mat.Cols() - 1);
 int numInRow = 0;
 int numInCol = 0;
 for (int i = 0 ; i < mat.Rows() ; i ++){
  if (i != nRow){
   for (int j = 0 ; j < mat.Cols() ; j++){
    if (j != nCol){
     res[numInRow][numInCol] = mat[i][j];
     numInCol += 1;
    }else{
     continue;
    }
   }
   numInRow += 1;
   numInCol = 0;
  }else{
   continue;
  }
 }
 return res;
}

template<class T>
double Matrix<T>::determinant(Matrix<T> mat){
 double det = 0;
 if (mat.Rows() > 2){
  for (int i = 0 ; i < mat.Rows() ; i++){
   det += pow(-1, i) * mat[0][i] * determinant(deleteRowsAndCols(mat, 0, i));
  }
 }else{
  det = mat[0][0] * mat[1][1] - mat[0][1]*mat[1][0];
 }
 return det;
}


template<class T>
Matrix<T> Matrix<T>::Inverse(){
 Matrix<T> inv(this->rows, this->cols);
 double det = determinant(*this);
 for (int i = 0; i < inv.Rows() ; i++){
  for (int j = 0; j < inv.Cols() ; j++){
   inv[i][j] = pow(-1, i+j)*determinant(deleteRowsAndCols(*this, i, j));
  }
 }
 Matrix<T> invT(inv.Rows(), inv.Cols());
 for (int i = 0 ; i < invT.Rows(); i++){
  for (int j = 0 ; j < invT.Cols(); j++){
   T temp = inv[j][i];
   T dividedTemp = temp / det;
   invT[i][j] = dividedTemp;
  }
 }
 return invT;
}

template<class T>
bool operator==(const Matrix<T> &m1, const Matrix<T> &m2){
 for (int i = 0; i < m1.Rows(); i++){
  for (int j = 0; j < m1.Cols(); j++){
   if (m1[i][j] != m2[i][j]){
    return false;
   }
  }
 }
 return true;
}

template<class T>
bool operator!=(const Matrix<T> &m1, const Matrix<T> &m2){
 return !(m1 == m2);
}

#endif // INCLUDE_MATRIX_HPP_
