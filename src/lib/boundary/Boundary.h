#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <iostream>
#include <vector>

class Boundary
{
    public:
        Boundary(){};
        virtual ~Boundary(){};

        void applyBCs();
        void assignBCs();
};

class InletBoundary : public Boundary
{
    public:
};

class OutletBoundary : public Boundary
{
    public:
};

class WallBoundary : public Boundary
{
    public:
};


class EdgeBoundary
{
    public:
        size_t size;
        std::vector<size_t> node;
};

class LeftEdgeBoundary : public EdgeBoundary
{
    public:
        size_t size;
        std::vector<size_t> node;

    private:
        void setStructuredEdgeBoundary();
};

class RightEdgeBoundary : public EdgeBoundary
{
    public:
        size_t size;
        std::vector<size_t> node;

    private:
        void setStructuredEdgeBoundary();
};

class TopEdgeBoundary : public EdgeBoundary
{
    public:
        size_t size;
        std::vector<size_t> node;

    private:
        void setStructuredEdgeBoundary();
};

class BottomEdgeBoundary : public EdgeBoundary
{
    public:
        size_t size;
        std::vector<size_t> node;
    
    private:
        void setStructuredEdgeBoundary();
};

class FrontEdgeBoundary : public EdgeBoundary
{
    public:
        size_t size;
        std::vector<size_t> node;

    private:
        void setStructuredEdgeBoundary();
};

class BackEdgeBoundary : public EdgeBoundary
{
    public:
        size_t size;
        std::vector<size_t> node;
    
    private:
        void setStructuredEdgeBoundary();
};

#endif



