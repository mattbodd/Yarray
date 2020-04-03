#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "yarray.h"

void allocate_contiguous(yarr *y, int total_elems) {
  // Allocate space for the contiguous array
  switch (y->tag) {
  case INT:
    y->data.idata = malloc(total_elems * sizeof(int));
    // Handle case where idata cannot be allocated
    if (y->data.idata == 0) {
      printf("%sCould not allocate space for %d ints%s\n",
	     RED, total_elems, NC);
      free(y);
    }
    
    break;
  case FLOAT:
    y->data.fdata = malloc(total_elems * sizeof(float));
    // Handle case where idata cannot be allocated
    if (y->data.fdata == 0) {
      printf("%sCould not allocate space for %d floats%s\n",
	     RED, total_elems, NC);
      free(y);
    }

    break;
  case LONG:
    y->data.ldata = malloc(total_elems * sizeof(long));
    // Handle case where idata cannot be allocated
    if (y->data.fdata == 0) {
      printf("%sCould not allocate space for %d longs%s\n",
	     RED, total_elems, NC);      
      free(y);
    }

    break;
  case DOUBLE:
    y->data.ddata = malloc(total_elems * sizeof(double));
    // Handle case where ddata cannot be allocated
    if (y->data.ddata == 0) {
      printf("%sCould not allocate space for %d doubles%s\n",
             RED, total_elems, NC);
      free(y);
    }

    break;
  default:
    #ifdef DEBUG
    printf("Invalid dataType received\n");
    #endif
    return;
  }
}

void realloc_contiguous(yarr *y, int total_elems) {
  // Allocate space for the contiguous array
  switch (y->tag) {
  case INT:
    {
      int *new_idata = realloc(y->data.idata, total_elems*sizeof(int));
      if (new_idata == NULL) {
        printf("%sCould not reallocate space for %d ints%s\n",
               RED, total_elems, NC);
        break;
      }
      y->data.idata = new_idata;
    
      break;
    }
  case FLOAT:
    {
      float *new_fdata = realloc(y->data.fdata, total_elems*sizeof(float));
      if (new_fdata == NULL) {
        printf("%sCould not reallocate space for %d floats%s\n",
               RED, total_elems, NC);
        break;
      }
      y->data.fdata = new_fdata;
    
      break;
    }
  case LONG:
    {
      long *new_ldata = realloc(y->data.ldata, total_elems*sizeof(long));
      if (new_ldata == NULL) {
        printf("%sCould not reallocate space for %d longs%s\n",
               RED, total_elems, NC);
        break;
      }
      y->data.ldata = new_ldata;

      break;
    }
  case DOUBLE:
    {
      double *new_ddata = realloc(y->data.ddata, total_elems*sizeof(double));
      if (new_ddata == NULL) {
        printf("%sCould not reallocate space for %d doubles%s\n",
               RED, total_elems, NC);
        break;
      }
      y->data.ddata = new_ddata;
    
      break;
    }
  default:
    #ifdef DEBUG
    printf("Invalid dataType received\n");
    #endif
    return;
  }
}

void fill_array(yarr *y, int total_elems, double fill_val) {
  if (y->tag == INT) {
    // Initialize int array
    #ifdef DEBUG
    printf("%sInitializing int array with %d elements\n%s",
           YELLOW, total_elems, NC);
    #endif

    //#pragma scop
    for (int i = 0; i < total_elems; i++) {
      for (int j = 0; j < 10; j++) {
             y->data.idata[i] = (int) fill_val;
      }
    }
    //#pragma endscop
    
    #ifdef DEBUG
    printf("%sDone initializing int array with %d elements%s\n",
           YELLOW, total_elems, NC);
    #endif
    
  } else if (y->tag == FLOAT) {
    // Initialize float array
    #ifdef DEBUG
    printf("%sInitializing float array with %d elements\n%s",
           YELLOW, total_elems, NC);
    #endif

    //#pragma scop
    for (int i = 0; i < total_elems; i++) {
      y->data.fdata[i] = (float) fill_val;
    }
    //#pragma endscop
    
    #ifdef DEBUG
    printf("%sDone initializing float array with %d elements%s\n",
           YELLOW, total_elems, NC);
    #endif
  } else if (y->tag == LONG) {
    // Initialize long array
    #ifdef DEBUG
    printf("%sInitializing float array with %d elements\n%s",
           YELLOW, total_elems, NC);
    #endif

    //#pragma scop
    for (int i = 0; i < total_elems; i++) {
      y->data.ldata[i] = (long) fill_val;
    }
    //#pragma endscop
    
    #ifdef DEBUG
    printf("%sDone initializing long array with %d elements%s\n",
           YELLOW, total_elems, NC);
    #endif
  } else if (y->tag == DOUBLE) {
    // Initialize double array
    #ifdef DEBUG
    printf("%sInitializing int array with %d elements\n%s",
           YELLOW, total_elems, NC);
    #endif

    //#pragma scop
    for (int i = 0; i < total_elems; i++) {
      y->data.ddata[i] = fill_val;
    }
    //#pragma endscop
    
    #ifdef DEBUG
    printf("%sDone initializing double array with %d elements%s\n",
           YELLOW, total_elems, NC);
    #endif
  } else {
    #ifdef DEBUG
    printf("Invalid dataType received\n");
    #endif
    return;
  }
}

yarr *C_array(double fill_val, dataType tag, int *widths, int dims) {
  // Create pointer to yarr strut where `yarr` strcut is simply an array of
  // primatives and meta information (eg: tag, dims, *strides)
  yarr *y = malloc(1 * sizeof(yarr));

  // Set the tag to be the value passed to function
  y->tag = tag;
  // Set dims to be the value passed to function
  y->dims = dims;
  // Widths represents the number of elements of successive dimension
  // eg: arr[3][2] -> width(arr[][]) = 2
  y->widths = widths;
  // Each dimension has a stride to access the next element in that dimension
  y->strides = malloc(dims * sizeof(int));
  // Strides are the product of widths of nested dimensions
  // Building up from the inner most dimension utilizes a growing product
  // Eg: widths of [ 6,  5, 4, 3] yield
  //               [60, 12, 3, 1]
  // Inner most dimension always has a stride of 1
  y->strides[dims-1] = 1;
  for (int i = dims-2; i >= 0; i--) {
    y->strides[i] = widths[i+1] * y->strides[i+1];
  }

  // Total number of primative elements held in `yarr.data`
  // The outer-most stride represents the total number of primative elements
  // nested in that slice of the outer-most dimension; the sum off all the slices
  // represents the total number of slices (in 2D thinking: L*W)
  int total_elems = y->strides[0]*widths[0];
  
  // Allocate contiguous array
  allocate_contiguous(y, total_elems);
  // Handle case where allocation fails
  if (y == 0) {
    return NULL;
  }

  // Initialize array
  fill_array(y, total_elems, fill_val);

  return y;
}

// Alternate call to `C_array` which expects `in_widths` as `{#, #, ..., #}`
yarr *_C_array(dataType tag, double fill_val, int dims, int *in_widths) {
  int *widths = malloc(dims * sizeof(int));

  memcpy(widths, in_widths, (dims * sizeof(int)));

  return C_array(fill_val, tag, widths, dims);
}

void dealloc_yarr(yarr *y) {
  // Free strides field
  free(y->strides);
  // Free correct data type
  switch (y->tag) {
  case INT:
    free(y->data.idata);
    break;
  case FLOAT:
    free(y->data.fdata);
    break;
  case LONG:
    free(y->data.ldata);
    break;
  case DOUBLE:
    free(y->data.ddata);
    break;
  default:
    printf(RED"Unsupported data type while deallocating\n"NC);
    break;
  }
  // Free yarr itself
  free(y);
}
