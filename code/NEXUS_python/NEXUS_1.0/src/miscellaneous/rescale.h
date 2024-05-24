
/* Rescales an array to have values in the given interval. */
template <typename T1>
void rescale(T1 *array,
             int const arraySize,
             T1 const minValue, T1 const maxValue)
{
    T1 min = array[0];
    T1 max = array[0];
    for (int i=1; i<arraySize; ++i) //find min and max values of the array
    {
        if ( min>array[i] ) min = array[i];
        if ( max<array[i] ) max = array[i];
    }
    
//     cout << "\nIn function 'rescale' found: minimum = " << min << " and maximum = " << max << " values inside the array.\n";
    T1 a = (maxValue-minValue) / (max-min);
    T1 b = minValue - a*min;
    
    for (int i=0; i<arraySize; ++i)
        array[i] = a * array[i] + b;
}

template <typename T1>
void rescale(T1 *array,
                int const arraySize,
                T1 const minValue, T1 const maxValue,
                T1 *constantValue, T1 *slope)
{
    T1 min = array[0];
    T1 max = array[0];
    for (int i=1; i<arraySize; ++i) //find min and max values of the array
    {
        if ( min>array[i] ) min = array[i];
        if ( max<array[i] ) max = array[i];
    }
    
//     cout << "\nIn function 'rescale' found: minimum = " << min << " and maximum = " << max << " values inside the array.\n";
    T1 a = (maxValue-minValue) / (max-min);
    T1 b = minValue - a*min;
    
    for (int i=0; i<arraySize; ++i)
        array[i] = a * array[i] + b;
    
    *constantValue = b;
    *slope = a;
}
