
#ifndef BOX_HEADER
#define BOX_HEADER


#include <sstream>
#include <string>

#include <misc.h>


/*
 * Structure to implement a bounding box and operations with it. The box coordinates (xMin,xMax,yMin,yMax,...) are stored in the 'coords' array.
 */

template<typename T, size_t N>
struct Box
{
    T coords[2*N];
    
    //! Class functions
    Box();                                              //initializes the elements to 0
    Box(T *values, size_t const noElements);            //initialize the elements to a set of boundaries (xMin,xMax,yMin,yMax,...)
    Box(Box &other);
    template <typename T2>
    void assign(T2 *values, size_t const noElements);    //assigns values to the elements
    bool isPointInBox(T *pos) const;                    //checks if a particle is inside the box
    bool isBoxOverlaping(Box const &other) const;       //checks if two boxes are overlapping
    void translate(T *offset);                          // translates the box according to the offset
    void addPadding(T *padding);                        // adds padding to the box
    void addPadding(T padding);                         // adds constant padding to the box (same on each side)
    T volume() const;                                   // computes the volume of the box
    void length(T *result);                             // returns the length of each side of the box
    std::string print(int const type=1) const;          // prints the box coordinates
    Box& operator =(Box const &other);                  // access to the box elements
    T& operator [](size_t const i) { return coords[i];} // access to the box elements
    size_t size() const { return 2*N;}                  // returns the number of coordinates (2*N)
    bool validBox();                                    // checks if the box volume >0 = valid box, oterwise returns false
};


template<typename T, size_t N> Box<T,N>::Box()
{
    for (size_t i=0; i<2*N; ++i)
        coords[i] = T(0.);
}
template<typename T, size_t N> Box<T,N>::Box(T *values, size_t const noElements)
{
    this->assign( values, noElements );
}
template<typename T, size_t N> Box<T,N>::Box(Box<T,N> &other)
{
    this->assign( other.coords, 2*N );
}

// assign values to the box coordinates
template<typename T, size_t N>
template <typename T2>
void Box<T,N>::assign(T2 *values, size_t const noElements)
{
    if (noElements!=2*N) throwError( "In function 'Box::assign'. The number of elements for the values argument should be 2*N, with N the dimension of the box. The program found '", noElements, "' input values when the box dimension is N='", N, "'." );
    for (size_t i=0; i<2*N; ++i)
        coords[i] = T(values[i]);
}

// Define the equality operator
template<typename T, size_t N> Box<T,N>& Box<T,N>::operator =(Box<T,N> const &other)
{
    for (size_t i=0; i<2*N; ++i)
        coords[i] = other.coords[i];
    return (*this);
}

// checks if a point ( an array with N values) is inside the box
template<typename T, size_t N> bool Box<T,N>::isPointInBox(T *pos) const
{
    if (N==3)
        return (pos[0]>=coords[0] and pos[0]<=coords[1] and pos[1]>=coords[2] and pos[1]<=coords[3] and pos[2]>=coords[4] and pos[2]<=coords[5]);
    else if (N==2)
        return (pos[0]>=coords[0] and pos[0]<=coords[1] and pos[1]>=coords[2] and pos[1]<=coords[3]);
    else if (N==1)
        return (pos[0]>=coords[0] and pos[0]<=coords[1]);
    else
    {
        for (int i=0; i<N; ++i)
            if (pos[i]<coords[2*i] or pos[i]>coords[2*i+1])
                return false;
        return true;
    }
}

// checks if second box is overlaping with this one
template<typename T, size_t N> bool Box<T,N>::isBoxOverlaping(Box const &other) const
{
    for (size_t i=0; i<N; ++i)
    {
        if ( this->coords[2*i]<=other.coords[2*i] and this->coords[2*i+1]>other.coords[2*i] )
            return true;
        if ( this->coords[2*i]<=other.coords[2*i+1] and this->coords[2*i+1]>other.coords[2*i+1] )
            return true;
    }
    return false;
}

// translates the box according to the offset
template<typename T, size_t N> void Box<T,N>::translate(T *offset)
{
    for (size_t i=0; i<N; ++i)
    {
        coords[2*i] += offset[i];
        coords[2*i+1] += offset[i];
    }
}

// increase box size by padding[i] along each edge
template<typename T, size_t N> void Box<T,N>::addPadding(T *padding)
{
    for (size_t i=0; i<N; ++i)
    {
        coords[2*i] -= padding[2*i];
        coords[2*i+1] += padding[2*i+1];
    }
}
template<typename T, size_t N> void Box<T,N>::addPadding(T padding)
{
    for (size_t i=0; i<N; ++i)
    {
        coords[2*i] -= padding;
        coords[2*i+1] += padding;
    }
}

// returns box volume
template<typename T, size_t N> T Box<T,N>::volume() const
{
    T temp = T(1.);
    for (size_t i=0; i<N; ++i)
        temp *= coords[2*i+1] - coords[2*i];
    return temp;
}
// returns box length along each edge
template<typename T, size_t N> void Box<T,N>::length(T *result)
{
    for (size_t i=0; i<N; ++i)
        result[i] = coords[2*i+1] - coords[2*i];
}

// outputs the box boundaries
template<typename T, size_t N> std::string Box<T,N>::print(int const type) const
{
    std::ostringstream buffer;
    if (type==1)
    {
        buffer << "{";
        for (size_t i=0; i<NO_DIM; ++i)
            buffer << "[" << coords[2*i] << "," << coords[2*i+1] << "],";
        buffer << char(8) << "}";
    }
    else if (type==2)
    {
        for (size_t i=0; i<NO_DIM; ++i)
            buffer << i << "-th box coordinate: minimum = " << coords[2*i] << "   to    maximum = " << coords[2*i+1] << "\n";
    }
    return buffer.str();
}

// check if there are specified coordinates for the box (i.e. if box !=[0,0] )
template<typename T, size_t N> bool Box<T,N>::validBox()
{
    for (size_t i=0; i<N; ++i)
        if ( coords[2*i]>=coords[2*i+1] )
            return false;
    if (this->volume()<=Real(0.) )
        return false;
    return true;
}
    



#endif
