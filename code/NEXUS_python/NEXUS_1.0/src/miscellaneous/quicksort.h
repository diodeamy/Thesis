


/* Quicksort algorithm for sorting an array in increasing order according to the elements. It sorts the array only between elements iMin to iMax.
For further details see: http://en.wikipedia.org/wiki/Quicksort
*/
template <typename T1, typename T2>
void quicksort( T1 values[], T2 iMin, T2 iMax )
{
    if ( iMin>=iMax )   //condition to stop the iterative computation
        return;
    
    T1 temp = values[iMax]; //the value against which we compare = the pivot
    T2 i = iMin, i1=iMin, i2=iMax;
    do
    {
        if ( values[i]>temp )   //if value larger than pivot
        {
            if ( i==i1 )    //move value to the right of the pivot
                values[i2] = values[i];
            //now the value is already to the right of the pivot
            --i2;
            i = i2;
        }
        else
        {
            if ( i==i2 )    //move value to the left of the pivot
                values[i1] = values[i];
            ++i1;
            i = i1;
        }
    } while ( i1!=i2 );
    
    //now copy the pivot to i=i1=i2 since that array element is free
    values[i] = temp;
    
    //we still have to order again the elements:
    quicksort( values, iMin, i-1);    //at the left of the pivot
    quicksort( values, i+1, iMax);    //at the right of the pivot
}


