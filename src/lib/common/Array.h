#ifndef ARRAY_H
#define ARRAY_H

#include <iostream>
#include <vector>

template <class T> 
class Array1D
{
    public:
        Array1D(int width) :
        width(width), data(width){}

        Array1D() :
        width(0){}

        inline T& operator()(int x)
        { return data[x]; }

        inline int size()
        { return data.size(); }

        inline void resize(int n)
        { width = n; data.resize(n); }

        inline void setZero()
        { for(T& n : data) n = 0; }

    private:
    	int width; 
	    std::vector<T> data;
};


template <class T> 
class Array2D
{
    public:
        Array2D(int height, int width) :
        width(width), height(height), data(width * height){}

        Array2D() :
        width(0), height(0){}

        inline T& operator()(int n)
        { return data[n]; }

        inline T& operator()(int y, int x)
        { return data[y * width + x]; }

        inline int size()
        { return data.size(); }

        inline void resize(int n, int m)
        { data.resize(n * m); }

    private:
    	int width, height; 
	    std::vector<T> data;

};


template <class T> 
class Array3D
{
    public:
        Array3D(int depth, int height, int width) :
        width(width), height(height), depth(depth), data(depth * width * height){}

        Array3D() :
        width(0), height(0), depth(0){}

        inline T& operator()(int n)
        { return data[n]; }

        inline T& operator()(int z, int y, int x)
        { return data[z * width * height + y * width + x]; }

        inline int size()
        { return data.size(); }

    private:
    	int width, height, depth; 
	    std::vector<T> data;

};

#endif

