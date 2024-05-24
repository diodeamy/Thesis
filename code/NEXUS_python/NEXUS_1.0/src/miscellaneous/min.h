
/* Returns the minimum of two numbers.*/
template <typename T1> inline T1 Min(T1 x, T1 y)
{
    return (x>y) ? y : x;
}


/* Returns the minimum of an array */
template <typename T, typename T_INT> inline T Min(T *x, T_INT const size )
{
    T temp = x[0];
    for (T_INT i=1; i<size; ++i)
        if ( temp>x[i] )
            temp = x[i];
    return temp;
}




/* Returns the minimum of an array and the index of the minimum. */
template <typename T, typename T_INT> inline void Min(T *x, T_INT const size, T *minValue, T_INT *indexMinimum )
{
    *minValue = x[0];
    *indexMinimum = 0;
    for (T_INT i=1; i<size; ++i)
        if ( *minValue > x[i] )
        {
            *minValue = x[i];
            *indexMinimum = i;
        }
}

