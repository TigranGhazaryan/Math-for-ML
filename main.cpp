#include <iostream>
#include "Matrix.hpp"
#include "Vector.hpp"
#include <vector>

int main()
{


  //////////////////////////////// INVERSE MATRIX, MATRIX ITERATOR ////////////////////////
  /*
  As an Example Matrix for inverse 
  1 2 3        -24  18  5
  0 1 4         20 -15 -4
  5 6 0         -5  4  1
  */

 
  Matrix<double> Matrix1(2,2);
  Matrix1.Print();
  for(auto& i : Matrix1)
    std::cin >> i;
 
 /* for(Matrix<int>::Iterator it = Matrix1.begin(); it != Matrix1.end(); ++it)
  {
    std::cout << "Matrix1.begin() = " << Matrix1.begin().Get_m_ptr() << '\n';
    std::cout << "it = " << it.Get_m_ptr() << '\n';
    std::cout << "Matrix1.end() = " << Matrix1.end().Get_m_ptr() << '\n';
   
    std::cin >> *it;
  }
*/


  std::cout << "Range-based for loop:\n";
  for(auto& i : Matrix1)
    std::cout << i << " ";

  std::cout << '\n';
  std::cout << "Iterator:\n";
  for(Matrix<double>::Iterator it = Matrix1.begin(); it != Matrix1.end(); it++)
    std::cout << *it << " ";
  std::cout << '\n';
 
  Matrix1.Print();
  Matrix1.Determinant();
  std::cout << "Matrix1 Determinant = " << Matrix1.Get_Determinant() << '\n';

  Matrix<double>Matrix2;
  Matrix2 = Matrix1.Inverse_Matrix();

  std::cout << "Matrix1 * Matrix2 = Identity Matrix\n";
  Matrix<double> Matrix3 = Matrix1 * Matrix2;
  Matrix3.Print();

  std::vector<std::pair<double, int>> Pivots = Matrix1.Get_Matrix_Pivots();
  std::cout << "Matrix1 Pivots are: \n";
  for(auto& i : Pivots)
    std::cout << i.first << " " << i.second << " ";
  std::cout << "\n";

  

  

//////////////////////// LINEAR INDEPENDENCE ///////////////////////  

  
  
  Matrix<double> Matrix4;
  int size = 2;
  Vector<double> V1(size);
  std::cout << "Input V1\n";   // 1 1 1
  for(int i = 0; i < V1.Get_size(); ++i)
    std::cin >> V1[i];
  Matrix4.Add_Vector(V1);

  std::cout << "Input V2\n";   // 1 2 3
  Vector<double> V2(size);
  for(int i = 0; i < V2.Get_size(); ++i)
    std::cin >> V2[i];
  Matrix4.Add_Vector(V2);

  std::cout << "Input V3\n";   // 1 3 5
  Vector<double> V3(size);
  for(int i = 0; i < V3.Get_size(); ++i)
    std::cin >> V3[i];
  //Matrix4.Add_Vector(V3);

  // V1.Print();
  // V2.Print();
  // V3.Print();  
  
  std::cout << "Linear Dependence Matrix:\n";
  Matrix4.Print();
  Matrix<double>Matrix5 = Matrix4.Gauss_Jordan_Elemination();
  std::cout << "Gauss-Jordan Elemination\n";
  Matrix5.Print();
  Matrix4.Linear_Dependence_Check();
  

  Matrix<double> Matrix6 = Matrix4.Basis();
  std::cout << "Basis for the Matrix4 is: \n";
  Matrix6.Print();

 std::cout <<"Matrix Dimension = Matrix rank = " << Matrix4.Get_Dimension() << '\n';
 
  Matrix<double> Inverse = Matrix4.Inverse_Matrix();

  
  // Ax = b; x = A^-1 * b
  Vector<double> Solution(size);
  std::cout << "Input Solution vector\n";
  for (int i = 0; i < Solution.Get_size(); ++i)
    std::cin >> Solution[i];

  Vector<double> Fin_Solution = Inverse * Solution;

  std::cout << "Solution: \n";
  Fin_Solution.Print(); 

}

