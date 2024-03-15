#ifdef ARRAY_H
#define ARRAY_H

#include <iostream>
#include <vector>

template <class T> 
class Array1D
{
    public:
        Array1D(size_t width) :
        width(width), data(width){}

        Array1D() :
        width(0){}

        inline T& operator()(size_t x)
        { return data[x]; }

        inline size_t size()
        { return data.size(); }

    private:
    	size_t width; 
	    vector<T> data;

};


template <class T> 
class Array2D
{
    public:
        Array2D(size_t width, size_t height) :
        width(width), height(height), data(width * height){}

        Array2D() :
        width(0), height(0){}

        inline T& operator()(size_t n)
        { return data[y * width + x]; }

        inline T& operator()(size_t x, size_t y)
        { return data[y * width + x]; }

        inline size_t size()
        { return data.size(); }

    private:
    	size_t width, height; 
	    vector<T> data;

};


template <class T> 
class Array3D
{
    public:
        Array3D(size_t width, size_t height, size_t depth) :
        width(width), height(height), depth(depth), data(length * width * height){}

        Array3D() :
        width(0), height(0), depth(0){}

        inline T& operator()(size_t x, size_t y, size_t z)
        { return data[z * width * height + y * width + x]; }

        inline size_t size()
        { return data.size(); }

    private:
    	size_t width, height, depth; 
	    vector<T> data;

};

#endif