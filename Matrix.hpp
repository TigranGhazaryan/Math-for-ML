#ifndef _MATRIX
#define _MATRIX

#include <iostream>
#include "Vector.hpp"
#include <iterator> 	// For std::forward_iterator_tag
#include <cstddef> 	// For std::ptrdiff_t
#include <vector>  	// For Pivots and Linear Independence check
#include <utility>  	// For std::pair



template<typename T>
class Matrix
{
private:
	int row = 0;
	int column = 0;
	double determinant = 1;
	int Dimension = 0;
	int Rank = 0;
	T** matrix;
	std::vector<int> Kernel_Columns;
	int nullity = 0;
	std::vector<std::pair<T, int>> Matrix_Pivot_List;
	
public:
	typedef T Value_Type;
public:
	Matrix();
	Matrix(int r, int c);
	Matrix(const Matrix<T>& M);
    ~Matrix();
	void Set_row(int r) {this->row = r; };
	void Set_column(int c){ this->column = c; };
	void Set_Matrix(int r, int c);
	T** Get_Matrix() const { return this->matrix; };
	int Get_row() const { return this->row; };
	int Get_column() const { return this->column; };
	void Set_Determinant(double d) { this->determinant = d;}
	double Get_Determinant() { return this->determinant; }
	void Set_Matrix_Pivots(const std::vector<std::pair<T, int>>& P) { this->Matrix_Pivot_List = P; }
	std::vector<std::pair<T, int>> Get_Matrix_Pivots() { return this-> Matrix_Pivot_List;}
	void Print() const;

	T* operator[](int index) const;
	Matrix<T>& operator=(const Matrix<T>& M);
    Matrix<T> operator+ (const Matrix<T>& M);
    Matrix<T> operator-(const Matrix<T>& M);
	bool operator==(const Matrix<T>& M);
	Matrix<T> operator*(T s);
	Matrix<T> operator* (const Matrix<T>& M);
	Vector<T> operator* (Vector<T>& V);
	Matrix<T> Identity_Matrix();
	Matrix<T> Transpose_Matrix();
	bool Symmetric();
	Matrix<T> Matrix_i_j(Matrix<T>& M, int i, int j);
	Matrix<T> Gauss_Jordan_Form();
	void Determinant();
	void Swap_Rows(Matrix<T>& M, int row1, int row2);
	void Check_Pivot(Matrix<T>& M, int i, int j);
	Matrix<T> Gauss_Jordan_Augmented_Elemination(Matrix<T>& M);
	Matrix<T> Gauss_Jordan_Elemination ();
	Matrix<T> Inverse_Matrix();
	Matrix<T> Add_Vector(const Vector<T>& v);
	void Linear_Dependence_Check();
	Vector<T> Matrix_to_Vector(Vector<T>& v, int col);
	Matrix<T> Basis();
	void Set_Dimension(int dim) { this->Dimension = dim; }
	void Set_Rank(int rank) { this->Rank = rank; }
	int Get_Dimension() { return this->Dimension; }
	int Get_Rank() { return this->Rank; }
	Matrix<T> Kernel_Basis();
	void Kernel_Print();
	void Set_Kernel_Columns(const std::vector<int>& K) { this->Kernel_Columns = K; }
	std::vector<int> Get_Kernel_Columns() { return this->Kernel_Columns; }
	void Set_nullity(int _nullity) { this->nullity = _nullity; };
	int Get_nullity() { return this->nullity; }


	// Inspiration / Sources for the Iterator: https://internalpointers.com/post/writing-custom-iterators-modern-cpp
	//					   					   https://www.youtube.com/watch?v=F9eDv-YIOQ0 by The Cherno

	class Iterator
	{
	public:
		typedef std::forward_iterator_tag iterator_category;
		typedef std::ptrdiff_t difference_type;
		typedef T value_type;
		typedef value_type* pointer;
		typedef value_type& reference;

	public:
		Iterator(pointer ptr, const Matrix<value_type>& Mat): m_ptr(ptr),  m_ptr2(ptr)
		{
			this->M = Mat;
			this->row = Mat.Get_row();
			this->column = Mat.Get_column();
			this->Pointer_Array = new pointer*[this->row];
			for(int i = 0; i < this->row; ++i)
				this->Pointer_Array[i] = new pointer[this->column];
			for(int i = 0; i < this->row; ++i)
				for(int j = 0; j < this->column; ++j)
					this->Pointer_Array[i][j] = &Mat[i][j];	

			//	M.Print();
		}; 
		Iterator (const Iterator& a) 
		{ 
			this->m_ptr = a.Get_m_ptr(); 
			this->m_ptr2 = a.Get_m_ptr2();
			this->row = a.Get_M().Get_row();
			this->column = a.Get_M().Get_column();
			this->Pointer_Array = new pointer*[this->row];
			for(int i = 0; i < this->row; ++i)
				this->Pointer_Array[i] = new pointer[this->column];
			for(int i = 0; i < this->row; ++i)
				for(int j = 0; j < this->column; ++j)
					this->Pointer_Array[i][j] = &a.Get_M()[i][j];
		};
		pointer Get_m_ptr() const 
		{ 
			return this->m_ptr; 
		}
		pointer Get_m_ptr2() const
		{
			return this->m_ptr2;
		}
		Matrix<value_type> Get_M() const
		{
			return this->M;
		}
		Iterator operator++() 
		{ 
			if(index == (this->row - 1))
			{
				m_ptr++;
				return *this;
			}
			if(index != (this->row - 1) && m_ptr - m_ptr2 == 
			this->Pointer_Array[this->index][this->column - 1] - 
			this->Pointer_Array[this->index][0])
			{
				this->index++;
				m_ptr = this->Pointer_Array[this->index][0];
				m_ptr2 = this->Pointer_Array[this->index][0];
				//std::cout << m_ptr << " =? " << m_ptr2 << '\n';
				//exit(0);
				return *this;
			} 
			
			m_ptr++;
			return *this; 
		};
		Iterator operator++(int) 
		{ 
			Iterator tmp = *this; 
			++(*this); 
			return tmp; 
		};
		Iterator operator--() 
		{ 
			m_ptr--;
			return *this; 
		};
		Iterator operator--(int) 
		{ 
			Iterator tmp = *this; 
			--(*this); 
			return tmp; 
		};
		reference operator[](int index) 
		{ 
			return *(m_ptr + index);
		}
		pointer operator->() 
		{ 
			return m_ptr; 
		};
		reference operator*() 
		{ 
			return *m_ptr; 
		};
		friend bool operator==(const Iterator& a, const Iterator& b)
		{ 
			return a.Get_m_ptr() == b.Get_m_ptr(); 
		}
		friend bool operator!=(const Iterator& a, const Iterator& b)
		{ 
			return a.Get_m_ptr() != b.Get_m_ptr();
		}
	private:
		pointer m_ptr;
		pointer m_ptr2;
		int index = 0;
		Matrix<value_type> M;
		int row;
		int column;
		pointer** Pointer_Array;
	};

	Iterator begin() 
	{ 
		return Iterator(&this->matrix[0][0], *this);
	};
	Iterator end() 
	{ 
		return Iterator(&(this->matrix[this->row - 1][this->column]), *this);
	};
	
};

//Default Constructor
template<typename T>
Matrix<T>::Matrix()
{
//	std::cout << "Default Constructor Called\n";
	Set_Matrix(0,0);
}
//Constructor
template<typename T>
Matrix<T>::Matrix(int r, int c)
{
//	std::cout << "Main Constructor Called\n";
	Set_Matrix(r,c);	
}

// Copy Constructor
template <typename T>
Matrix<T>::Matrix(const Matrix<T>& M)
{
//	std::cout << "Copy Constructor Called\n";
	*this = M;
}

//Destructor
template<typename T>
Matrix<T>::~Matrix()
{
	//std::cout << "Destructor Called\n";
	for(int i = 0; i < this->row; ++i)
		delete[] this->matrix[i];
	delete[] this->matrix;		
}

//Access Operator [] const
template<typename T>
T* Matrix<T>::operator[](int index) const
{
    if(index >= this->row){
        std::cout << "Index Out of bounds\n";
        exit(0);
    }
    return this->matrix[index];
}

//Set Matrix
template<typename T>
void Matrix<T>::Set_Matrix(int r, int c)
{
	Set_row(r);
	Set_column(c);
	this->matrix = new T*[this->row];
	for(int i = 0; i < this->row; ++i)
		this->matrix[i] = new T[this->column];

	for(int i = 0; i < this->row; ++i)
		for(int j = 0; j < this->column; ++j)
			this->matrix[i][j] = 0;	
}

//Print
template<typename T>
void Matrix<T>::Print() const
{
	for(int i = 0; i < this->row; ++i)
	{
		for(int j = 0; j < this->column; ++j)
			std::cout << this->matrix[i][j] << " ";
		std::cout << '\n';
	}
	std::cout << '\n';
}

//Assignment Operator =
template <typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& M)
{
	if(this->row < M.Get_row() || this->column < M.Get_column())
		Set_Matrix(M.Get_row(), M.Get_column());

	for(size_t i = 0; i < M.Get_row(); ++i)
		for(size_t j = 0; j < M.Get_column(); ++j)
			this->matrix[i][j] = M[i][j];
	return *this;
}

//Operator +
template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& M)
{
	if(this->row != M.Get_row() || this->column != M.Get_column())
	{
		std::cout << "Rows and Columns must be equal\n";
		exit(0);
	}
	Matrix<T> M2(this->row, this->column);
	for(size_t i = 0; i < this->row; ++i)
		for(size_t j = 0; j < this->column; ++j)
			M2[i][j] = this->matrix[i][j] + M[i][j];

	return M2;
}

//Operator -
template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& M)
{
	if(this->row != M.Get_row() || this->column != M.Get_column())
	{
		std::cout << "Rows and Columns must be equal\n";
		exit(0);
	}
	Matrix<T> M2(this->row, this->column);
	for(size_t i = 0; i < this->row; ++i)
		for(size_t j = 0; j < this->column; ++j)
			M2[i][j] = this->matrix[i][j] - M[i][j];
	return M2;
}

// Operator == (equality)
template<typename T>
bool Matrix<T>::operator==(const Matrix<T>& M)
{
	bool Are_Equal = true;
	if(this->row != M.Get_row() || this->column != M.Get_column())
	{
		Are_Equal = false;
		return Are_Equal;
	}

	for(int i = 0; i < this->row; ++i)
	{
		for(int j = 0; j < this->column; ++j)
		{
			if(matrix[i][j] != M[i][j])
			{
				Are_Equal = false;
				return Are_Equal;
			}
		}
	}

	return Are_Equal;
}

//Multiplication by Scalar *
template <typename T>
Matrix<T> Matrix<T>::operator*(T s)
{
	Matrix<T> temp(this->row, this->column);
	for(int i = 0; i < this->row; ++i)
		for(int j = 0; j < this->column; ++j)
			temp[i][j] = this->matrix[i][j] * s;
	return temp;
}

//Matrix Multiplication
template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& M2)
{
	if(Get_column()!= M2.Get_row())
	{
		std::cout << "Row and Column must be equal\n";
		exit(0);
	}
	Matrix<T> M3(this->row, M2.Get_column());
	for(size_t i = 0; i < this->row; ++i)
	{
		for(size_t j = 0; j < M2.Get_column(); ++j)
		{	
			int m = 0; 
			int n = 0;
			while(m != M2.Get_row() && n != this->column)
			{
			// Commented pieces are for Printing the Process if Needed	
			//	std::cout << "M3[ " << i << " ][ " << j << " ] += ";
			//	std::cout << matrix[i][m] << " * " << M2[n][j] << " = ";
				M3[i][j] += this->matrix[i][m] * M2[n][j]; 
			//	std::cout << M3[i][j] << std::endl;
				++m;
				++n;
			}
		}
	}
	return M3;
}

// Multiplication (Transformation) of a Vector
template<typename T>
Vector<T> Matrix<T>::operator* (Vector<T>& V)
{
	Matrix<T> temp (V.Get_size(), 1);
	if(V.Get_size() != this->column)
	{
		std::cout << "Matrix Column must be equal to the Rows of the Vector\n";
		exit(0);
	}
	for(int i = 0; i < V.Get_size(); ++i)
		temp[i][0] = V[i];

	Matrix<T> temp2 = *this;
	Matrix<T> temp3 = temp2 * temp;
	V.Set_Vec(temp3.Get_row());
	for(int i = 0; i < V.Get_size(); ++i)
		V[i] = temp3[i][0];

	return V;
}


// Identity Matrix
template<typename T>
Matrix<T> Matrix<T>::Identity_Matrix()
{
	Matrix<T> Identity = *this;
	for(int i = 0; i < Identity.Get_row(); ++i)
	{
		for(int j = 0; j < Identity.Get_column(); ++j)
		{
			if(i != j)
				Identity[i][j] = 0;
			else
				Identity[i][j] = 1;
		}
	}
	return Identity;
}

// Matrix Transpose
template<typename T>
Matrix<T> Matrix<T>::Transpose_Matrix()
{
	Matrix<T> Transpose (this->column, this->row);
	for(int i = 0; i < this->row; ++i)
		for(int j = 0; j < this->column; ++j)
			Transpose[j][i] = this->matrix[i][j];
		
	return Transpose;
}

// Symmetricity 
template<typename T>
bool Matrix<T>::Symmetric()
{
	Matrix<T> M1 = *this;
	Matrix<T> M2 = Transpose_Matrix();
	return M1 == M2;

}
//
template <typename T>
Matrix<T> Matrix<T>::Matrix_i_j(Matrix<T>& M, int i, int j)
{
	if(i > M.Get_row() || j > M.Get_column())
	{
		std::cout << "Indexes out of bounds\n";
		exit(0);
	}

	Matrix<T> Matrix_i_j(M.Get_row() - 1, M.Get_column() - 1);
	for(int m = 0; m < Matrix_i_j.Get_row(); ++m)
	{
		if(m == i)
		{
			for(int m = i; m < Matrix_i_j.Get_row(); ++m)
				for(int n = 0; n < Matrix_i_j.Get_column(); ++n)
					Matrix_i_j[m][n] = M[m + 1][n + 1];
			break;
		}
		else 
			for(int n = 0; n < Matrix_i_j.Get_column(); ++n)
				Matrix_i_j[m][n] = M[m][n + 1];			
	}
			 
//	Matrix_i_j.Print();
//	std::cout << "\n";
	return Matrix_i_j;
}


// Gauss Jordan Form
template <typename T>
Matrix<T> Matrix<T>::Gauss_Jordan_Form() 
{
	Matrix<T> Gauss_Jordan_F(this->row, 2 * this->column);
	Gauss_Jordan_F = *this;
	for(int i = 0; i < Gauss_Jordan_F.Get_row(); ++i)
		for(int j = 0; j < Gauss_Jordan_F.Get_column(); ++j)
			if(i == j)
				Gauss_Jordan_F[i][(this->column) + i] = 1;
	
//	std::cout << "Gauss-Jordan Form:\n";
//	Gauss_Jordan_F.Print();
	return Gauss_Jordan_F;
}


// Swap Rows of a Matrix
template<typename T>
void Matrix<T>::Swap_Rows(Matrix<T>& M, int row1, int row2)
{
	if(row1 >= M.Get_row() || row2 >= M.Get_row())
	{
		std::cout << "Input rows are out of bounds\n";
		exit(0);
	}

	for(int j = 0; j < M.Get_column(); ++j)
	{
			T temp = M[row1][j];
			M[row1][j] = M[row2][j];
			M[row2][j] = temp;
	}
//	std::cout << "Swap rows " << ++row1 << " and " << ++row2 << "\n";
//	M.Print();
}


// Check Pivots of a Matrix
template <typename T>
void Matrix<T>::Check_Pivot(Matrix<T>& M, int i, int j)
{
	if(M[i][j] == 0)
	{
		for(int n = i + 1; n < M.Get_row(); ++n)
		{
			if(M[n][j] != 0)
			{
				Swap_Rows(M, i, n);
				break;
			}
		}
	}
}



//Gauss Jordan Augmented Elemination 
template <typename T>
Matrix<T> Matrix<T>::Gauss_Jordan_Augmented_Elemination(Matrix<T>& M)
{
	Matrix<T> Gauss_Jordan_Aug_Elem;
	Gauss_Jordan_Aug_Elem = Gauss_Jordan_Form();	
	
	for(int i = 0; i < Gauss_Jordan_Aug_Elem.Get_row(); ++i)
	{
		for(int j = 0; j < M.Get_column(); ++j)
		{
			// If that row has a pivot, dont check anymore 
			if((j == 0 && i == 0) || i == M.Get_Matrix_Pivots().size()) 
			{
				// Commented pieces are for Printing the Process if Needed
			
				if(Gauss_Jordan_Aug_Elem[i][j] == 0)
				{
					while(Gauss_Jordan_Aug_Elem[i][j] == 0)
					{
						Check_Pivot(Gauss_Jordan_Aug_Elem, i, j);
						
						//	Set Determinant
						if(i == j)
						{
							double Determ = M.Get_Determinant();
							Determ *= -1;
							if(Gauss_Jordan_Aug_Elem[i][j] == 0)
								Determ *= 0;
							M.Set_Determinant(Determ);
						}
						if(Gauss_Jordan_Aug_Elem[i][j] == 0)
						{
							++j;
							if(j == M.Get_column())
							{
							//	std::cout << "Gauss-Jordan Augmented Partial Elemination\n";
							//	Gauss_Jordan_Aug_Elem.Print();
								return Gauss_Jordan_Aug_Elem;
							}
						}
					}
				}
				if(Gauss_Jordan_Aug_Elem[i][j] != 0)
				{
					if(i == j)
						{
							double Determ = M.Get_Determinant();
							Determ *= Gauss_Jordan_Aug_Elem[i][j];
							M.Set_Determinant(Determ);
						}
				}

				T Pivot = Gauss_Jordan_Aug_Elem[i][j];
				if(Pivot != 0)
				{
				//	Set Pivots
					
					std::vector<std::pair<T, int>> Piv = M.Get_Matrix_Pivots();
					std::cout << std::endl;
					Piv.emplace_back(Pivot, j);
					M.Set_Matrix_Pivots(Piv);
				}

				// Gauss-Jordan method
				// Step 1. Make the Pivot = 1
				for(int n = 0; n < Gauss_Jordan_Aug_Elem.Get_column(); ++n)
				{
					if(Pivot != 0)
						Gauss_Jordan_Aug_Elem[i][n] /= Pivot;
				}

				// Step 2. Change all members under the Pivot to 0
				for(int m = 0; m < Gauss_Jordan_Aug_Elem.Get_row(); ++m)
				{
					for(int n = 0; n < Gauss_Jordan_Aug_Elem.Get_column(); ++n)
					{
						if(m != i && Gauss_Jordan_Aug_Elem[m][j] != 0)
						{
							T Change = Gauss_Jordan_Aug_Elem[i][n] * Gauss_Jordan_Aug_Elem[m][j];
							if(n != j)
								Gauss_Jordan_Aug_Elem[m][n] = Gauss_Jordan_Aug_Elem[m][n] - Change;
						}
					}
					if(m != i)
						Gauss_Jordan_Aug_Elem[m][j] = 0;
				}				
			}
		}
	}
//	std::cout << "Gauss-Jordan Augmented Elemination\n";
//	Gauss_Jordan_Aug_Elem.Print();
	return Gauss_Jordan_Aug_Elem;
}



// Matrix Determinant
template <typename T>
void Matrix<T>::Determinant()
{
//	std::cout << "\n\n\n";
//	std::cout << "///////////// MATRIX DETERMINANT CALCULATION /////////\n\n";
	Matrix<T> ThisMatrix;
	ThisMatrix = *this;
	Gauss_Jordan_Augmented_Elemination(ThisMatrix);	
	if(this->row == this->column)		
		Set_Determinant(ThisMatrix.Get_Determinant());
	else
		this->determinant = 0;
	Set_Matrix_Pivots(ThisMatrix.Get_Matrix_Pivots());
}



// Gauss Jordan Elemination Form
template <typename T>
Matrix<T> Matrix<T>::Gauss_Jordan_Elemination()
{
	Matrix<T> Gauss_Jordan_Elem(this->row, this->column);
	Matrix<T> Gauss_Jordan_Aug;
	Matrix<T> temp = *this;
	Gauss_Jordan_Aug = Gauss_Jordan_Augmented_Elemination(temp);
	for(int i = 0; i < Gauss_Jordan_Elem.Get_row(); ++i)
		for(int j = 0; j < Gauss_Jordan_Elem.Get_column(); ++j)
			Gauss_Jordan_Elem[i][j] = Gauss_Jordan_Aug[i][j];

	//std::cout << "Gauss-Jordan Elemination\n";
	//Gauss_Jordan_Elem.Print();
	return Gauss_Jordan_Elem;
}


// Matrix Inverse
template <typename T>
Matrix<T> Matrix<T>::Inverse_Matrix()
{
	Determinant();
	if(this->row != this->column || this->determinant == 0)
	{
		std::cout << "This Matrix Doesnt have an Inverse\n";
		std::cout << "Determinant is " << this->determinant << '\n';
		
		Matrix<T> Inverse = Identity_Matrix();
		return Inverse;
	}

//	std::cout << "\n\n\n";
//	std::cout << "///////////// MATRIX INVERSE CALCULATION /////////\n\n";
	Matrix<T> ThisMatrix;
	ThisMatrix = *this;

	Matrix<T> Gauss_Jordan_Aug;
	Gauss_Jordan_Aug = Gauss_Jordan_Augmented_Elemination(ThisMatrix);
	
	Matrix<T> Inverse(this->row, this->column);
	for(int i = 0; i < this->row; ++i)
		for(int j = 0; j < this->column; ++j)
			Inverse[i][j] = Gauss_Jordan_Aug[i][j + (this->column)];
	
//	std::cout << "Inverse Matrix\n";
//	Inverse.Print();
	return Inverse;
}


// Add Vector to a Matrix
template<typename T>
Matrix<T> Matrix<T>::Add_Vector(const Vector<T>& v)
{
	Matrix<T> Add_Vec(v.Get_size(), this->column + 1);
	Add_Vec = *this;
	for(int i = 0; i < Add_Vec.Get_row(); ++i)
		Add_Vec[i][Add_Vec.Get_column() - 1] = v[i];
	
	*this = Add_Vec;
	return *this;
}


// Linear (In)Dependence Check for a Matrix
template<typename T>
void Matrix<T>::Linear_Dependence_Check()
{
	Determinant();
	for(int i = 0; i < Get_Matrix_Pivots().size(); ++i)
	{
		std::cout << "Column" << Get_Matrix_Pivots()[i].second << " Pivot = "
					  << Get_Matrix_Pivots()[i].first << "\n"
					  << "V" << Get_Matrix_Pivots()[i].second  << " is Linearly Independent\n\n";
	}
		
}


// Matrix to Vector(s) from Columns
template<typename T>
Vector<T> Matrix<T>::Matrix_to_Vector(Vector<T>& v, int col)
{
	Matrix<T> temp = *this;
	if(v.Get_size() < temp.Get_row())
		v.Set_Vec(temp.Get_row());
	for(int i = 0; i < v.Get_size(); ++i)
		v[i] = temp[i][col];

	return v;
}


// Matrix Basis
template<typename T>
Matrix<T> Matrix<T>::Basis()
{
	Determinant();
	Matrix<T> Basis_Vectors;
	for(int i = 0; i < Get_Matrix_Pivots().size(); ++i)
	{
		if(Get_Matrix_Pivots()[i].first != 0)
		{
			Vector<T> temp(Get_row());
			Vector<T> v = Matrix_to_Vector(temp, Get_Matrix_Pivots()[i].second);
			Basis_Vectors.Add_Vector(v);
		}
	}

	// Set Dimension and Rank
	Set_Dimension(Basis_Vectors.Get_column());
	Set_Rank(Basis_Vectors.Get_column());
	Set_nullity(this->column - Get_Rank());
	return Basis_Vectors;
}


// Kernel | Null Space of a Matrix
template<typename T>
Matrix<T> Matrix<T>::Kernel_Basis()
{
	Matrix<T> Bas = Basis();
	Matrix<T> temp = Gauss_Jordan_Elemination();
	Matrix<T> _Kernel;
	if(Get_Rank() == this->column)
	{
		Vector<T> v(this->column);
		_Kernel.Add_Vector(v);
		return _Kernel;
	}

	// Keep the non-Independent columns
	std::vector<int> Kernel_columns;
	
	for(int i = 0; i < temp.Get_column(); ++i)
	{
		bool Kernel_col = true;
		for(int j = 0; j < Matrix_Pivot_List.size(); ++j)
		{
			if(i == this->Matrix_Pivot_List[j].second )
			{
				Kernel_col = false;
				break;
			}
		}
		if(Kernel_col == true)
		{	
			Kernel_columns.push_back(i);	
		}
	}
	
	Set_Kernel_Columns(Kernel_columns);
	for(int m = 0; m < Kernel_columns.size(); ++m)
	{
		
		Vector<T> v(temp.Get_column());
		v[Kernel_columns[m]] = 1;
		for(int i = 0; i < temp.Get_row(); ++i)
		{
			for(int j = 0; j < temp.Get_column(); ++j)
			{
				if(temp[i][j] == 1)
				{	
					v[j] += temp[i][j] * (-1) * temp[i][Kernel_columns[m]];
					break;
				}
			}
		}
		_Kernel.Add_Vector(v);
	}


	return _Kernel;
}


// Print Null Space of a Matrix
template<typename T>
void Matrix<T>::Kernel_Print() 
{
	Matrix<T> Kern = Kernel_Basis();
	std::cout << "Matrix Kernel / Null Space\n";
	if(Get_Rank() == this->column)
	{
		std::cout << "The system has a unique solution:\n";
		Kern.Print();
		std::cout << "\nnullity = " << Get_nullity() << '\n';
	}
	else
	{
		std::cout << "Vectors {";
		for(int i = 0; i < Kern.Get_column(); ++i)
		{
			int col = Get_Kernel_Columns()[i] + 1;
			std::cout << " v(" << col <<")";
			if(i != Kern.Get_column() - 1)
				std::cout << " , ";
			else
				std::cout << " } :\n";	
		}
		Kern.Print();

		std::cout << "The system has Infinitely many Solutions.\n";
		std::cout << "For: \n";

		for(int j = 0; j < this->column; ++j)
		{
			int col = j + 1;
			if(j == this->column / 2)
			{
				std::cout << "c(" << col << ") = ";
				for(int i = 0; i < Kern.Get_column(); ++i)
				{
					int col2 = Get_Kernel_Columns()[i] + 1;
					std::cout << "c(" << col2 << ") * v(" << col2 <<")";
					if(i != Kern.Get_column() - 1)
						std::cout << " + ";
					else
						std::cout << '\n';		
				}
			}
			else	
				std::cout << "c(" << col << ")\n";
		
		}
		std::cout << "\nnullity = " << Get_nullity() << '\n';
	}
	
}

#endif
