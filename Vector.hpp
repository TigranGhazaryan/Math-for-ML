#ifndef _VECTOR
#define _VECTOR

#include <iostream>
#include <iterator>
#include "Matrix.hpp"
#include <cmath> 	// For Vector Norm

template<typename T>
class Matrix;

template<typename T>
class Vector
{
private:
    T* vector;
    int size = 0;
    Matrix<T>* Bases;

public:
    Vector();
    Vector(int s);
    Vector(const Vector<T>& v);
    ~Vector();
    void Set_size(int size) { this->size = size; }
    void Set_Vec(int size);
    int Get_size() const { return this->size; }
    T* Get_Vec() const { return this->vector; }
    void Set_Default_Bases();
    void Set_Bases(const Matrix<T>& M);
    Matrix<T> Get_Bases() const { return *this->Bases; }
    void Print_Bases() const;
    void Print () const;
    T& operator[](int index) const;
    Vector<T>& operator=(const Vector<T>& v);
    Vector<T> operator+(const Vector<T>& v);
    Vector<T> operator-(const Vector<T>& v);
    Vector<T>& operator*(const T& s); 
    Matrix<T> Transpose();
    double Manhathan_Norm();
    double Euclidean_Norm();
    double Scalar_Product(Vector<T>& v);

    // Inspiration / Sources for the Iterator: https://internalpointers.com/post/writing-custom-iterators-modern-cpp
	//					                       https://www.youtube.com/watch?v=F9eDv-YIOQ0 by The Cherno
    class Iterator
    {
    public:
        typedef std::forward_iterator_tag iterator_category;
		typedef std::ptrdiff_t difference_type;
		typedef T value_type;
		typedef value_type* pointer;
		typedef value_type& reference;
    public:
        Iterator (pointer ptr)
        { 
            this->m_ptr = ptr;
        }
        reference operator*()const 
        { 
            return *this->m_ptr; 
        }
        pointer operator->() 
        { 
            return this->m_ptr; 
        }
        Iterator& operator++()
        { 
            this->m_ptr++; 
            return *this;
        }
        Iterator operator++(int) 
        { 
            Iterator temp = *this; 
            ++(*this); 
            return temp;
        }
        friend bool operator==(const Iterator& a, const Iterator& b) 
        { 
            return a.m_ptr == b.m_ptr; 
        }
        friend bool operator!=(const Iterator& a, const Iterator& b)
        { 
            return a.m_ptr != b.m_ptr; 
        }
    private:
        pointer m_ptr;

    };

    Iterator begin() { return Iterator(&(this->vector[0])); }
    Iterator end()   { return Iterator(&(this->vector[size])); }

};

//Default Constructor
template<typename T>
Vector<T>::Vector()
{
    Set_Vec(1);
    Set_Default_Bases();
}
//Constructor
template<typename T>
Vector<T>::Vector(int s)
{
    Set_Vec(s);
    Set_Default_Bases();
}

//Copy Constructor
template<typename T>
Vector<T>::Vector(const Vector<T>& v)
{
    *this = v;
}

//Destructor
template<typename T>
Vector<T>::~Vector()
{
    delete[] this->vector;
}

//Access Operator[]
template<typename T>
T& Vector<T>::operator[](int index) const
{
    if(index >= this->size){
        std::cout << "Index Out of bounds\n";
        exit(0);
    }
    return this->vector[index];
}

//Setter
template<typename T>
void Vector<T>::Set_Vec(int size)
{
    Set_size(size);
    this->vector = new T[this->size];
    for(int i = 0; i < this->size; ++i)
        this->vector[i] = 0;
}

// Set Bases
template<typename T>
void Vector<T>::Set_Bases(const Matrix<T>& M)
{
    this->Bases = new Matrix<T>(M.Get_row(), M.Get_column());
    *this->Bases = M;
    
    Vector<T> temp = Get_Bases() * (*this);
    Set_Vec(temp.Get_size());
    *this = temp;
}

template<typename T>
void Vector<T>::Set_Default_Bases()
{
    Matrix<T> temp(this->size, this->size);
    temp = temp.Identity_Matrix();
    Set_Bases(temp);

    *this = Get_Bases() * (*this);
}

//Print Bases
template<typename T>
void Vector<T>::Print_Bases() const
{
    std::cout << "Vector Bases are:\n";
    this->Bases->Print();
}

//Print
template<typename T>
void Vector<T>::Print() const
{
   for(int i = 0; i < this->size; ++i)
        std::cout << this->vector[i] << std::endl;
    std::cout << std::endl;
}

//Operator =
template<typename T>
Vector<T>& Vector<T>::operator=(const Vector<T>& v)
{
    if(this->size < v.Get_size())
        Set_Vec(v.Get_size());

    for(size_t i = 0; i < this->size; ++i)
       this->vector[i] = v[i];
    return *this;
}

//Operator +
template<typename T>
Vector<T> Vector<T>::operator+(const Vector<T>& v)
{
    Vector<T> Vec3(this->size);
    for(int i = 0; i < this->size; ++i)
        Vec3[i] = this->vector[i] + v[i];
    return Vec3;
}

//Operator -
template<typename T>
Vector<T> Vector<T>::operator-(const Vector<T>& v)
{
    Vector<T> Vec3(this->size);
    for(int i = 0; i < this->size; ++i)
        Vec3[i] = this->vector[i] - v[i];
    return Vec3;
}

//Multiplication by a Scalar
template<typename T>
Vector<T>& Vector<T>::operator*(const T& s)
{
    for(size_t i = 0; i < this->size; ++i)
        this->vector[i]*= s;
    return *this;
}

// Vector Transpose
template<typename T>
Matrix<T> Vector<T>::Transpose()
{
    Vector<T> temp = *this;
    Matrix<T> temp2;
    temp2.Add_Vector(temp);
    Matrix<T> Transpose = temp2.Transpose_Matrix();
    //Transpose.Print();
    return Transpose;
}

// Manhathan Norm
template<typename T>
double Vector<T>::Manhathan_Norm()
{
    Vector<T> temp = *this;
    T Norm = 0;
    for(auto i : temp)
        Norm += std::abs(i);
    return Norm;
}

//Euclidean Norm
template<typename T>
double Vector<T>::Euclidean_Norm()
{
    Matrix<T> temp = Transpose();
    Vector<T> temp2 = *this;
    Vector<T> outcome = temp * temp2;
    double Norm = std::sqrt(outcome[0]);
    return Norm;
}

//Scalar Product
template<typename T>
double Vector<T>::Scalar_Product(Vector<T>& v)
{
    if(v.Get_size() != this->size)
    {
        std::cout << "Vector size must be equal\n";
        exit(0);
    }
    Matrix<T> temp = v.Transpose();
    Vector<T> temp2 = *this;

    Vector<T> Product = temp * temp2;
    double _Scalar = Product[0];
    return _Scalar;
}



#endif
