#pragma once
#include "Matrix.hpp"

template<typename T>
class Vector
{
private:
    T* vector;
    int size = 0;
    
public:
    Vector();
    Vector(int s);
    Vector(const Vector<T>& v);
    ~Vector();
    void Set_size(int size) { this->size = size; }
    void Set_Vec(int size);
    int Get_size() const { return this->size; }
    T* Get_Vec() const { return this->vector; }
    void Print () const;
    T& operator[](int index) const;
    Vector<T>& operator=(const Vector<T>& v);
    Vector<T> operator+(const Vector<T>& v);
    Vector<T> operator-(const Vector<T>& v);
    Vector<T>& operator*(const T& s);
};

//Default Constructor
template<typename T>
Vector<T>::Vector()
{
    Set_Vec(1);
}
//Constructor
template<typename T>
Vector<T>::Vector(int s)
{
    Set_Vec(s);
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

