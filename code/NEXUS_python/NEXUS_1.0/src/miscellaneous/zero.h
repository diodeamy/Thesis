

/* Initializes all the values of an array to 0.
*/
template <typename T1, typename T2>
void zero(T1 array[],
          T2 const size)
{
    for (T2 i=0; i<size; ++i)
        array[i] = 0;
}


