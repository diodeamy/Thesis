
/* Returns the minimum of two numbers.*/
template <typename T1> inline T1 Max(T1 x, T1 y)
{
    return (y>x) ? y : x;
}


/* Returns the maximum of an array */
template <typename T, typename T_INT> inline T Max(T *x, T_INT const size )
{
    T temp = x[0];
    for (T_INT i=1; i<size; ++i)
        if ( temp<x[i] )
            temp = x[i];
    return temp;
}




/* Returns the maximum of an array and the index of the maximum. */
template <typename T, typename T_INT> inline void Max(T *x, T_INT const size, T *maxValue, T_INT *indexMaximum )
{
    *maxValue = x[0];
    *indexMaximum = 0;
    for (T_INT i=1; i<size; ++i)
        if ( *maxValue < x[i] )
    {
        *maxValue = x[i];
        *indexMaximum = i;
    }
}

