#ifndef YARRAY_H
#define YARRAY_H
// Possible primative data types of `yarr.data`
typedef enum {INT, FLOAT, LONG, DOUBLE} dataType;
typedef struct yarr yarr;

// Complex array type implemented in C
// Modeled after `ndarray` from Numpy
struct yarr {
  int dims;
  // Widths hold the number of elements in the next dimension
  int *widths;
  // Strides are organized from outer dimensions to inner dimensions
  // eg: arr[5][4][3] -> strides{12, 3, 1}
  // where each stride denotes how many elements are in the nested dimensions
  // To access arr[2][1][0] ->
  //   *arr + ((2*strides[0] + 1*strides[1] + 0*strides[2]) * sizeof(int))
  int *strides;
  dataType tag;
  
  union data {
    int    *idata;
    float  *fdata;
    long   *ldata;
    double *ddata;
  } data;
};
#endif /* YARRAY_H */
