#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "yarray.h"
#include "yutilities.h"

// Macro definitions
#define max(a,b) ({__typeof__ (a) _a = (a);     \
      __typeof__ (b) _b = (b);                  \
      _a > _b ? _a : _b; })

#define min(a,b) ({__typeof__ (a) _a = (a);     \
      __typeof__ (b) _b = (b);                  \
      _a < _b ? _a : _b; })

// Possible pairwise operators
typedef enum {ADD, SUB, MULT, DIV, MIN, MAX} pairwise_op;

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
      printf("Test\n");
      int *new_idata = realloc(y->data.idata, total_elems);
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
      float *new_fdata = realloc(y->data.fdata, total_elems);
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
      long *new_ldata = realloc(y->data.ldata, total_elems);
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
      double *new_ddata = realloc(y->data.ddata, total_elems);
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

// Reshaping an array maintains the total number of elements and simply requires
// readjusting strides
// This is equivalent to Yorick `reform`
void reform_C_array(yarr *y, int *new_widths, int new_dims) {
    // Total number of elements in existing array
  int total_existing_elems = y->strides[0] * y->widths[0];
  
  // Total number of elements in reshaped array is the product of each width
  int total_proposed_elems = 1;
  for (int i = 0; i < new_dims; i++) {
    total_proposed_elems *= new_widths[i];
  }

  // If the number of elements remains unchanged, simply change strides
  if (total_existing_elems == total_proposed_elems) {
    // If the number of dimensions are changing, ensure that widths, strides and
    // dims are updated appropriately
    if (y->dims != new_dims) {
      y->dims = new_dims;
      // Do not need to free or realloc `y->widths` as `y->widths` was set equal
      // to a pointer that was allocated elsewhere
      
      // Try to make use of previously allocated memory
      int *new_strides = realloc(y->strides, new_dims);
      if (new_strides == NULL) {
        printf("%sCould not reallocate for new strides of size %d%s\n",
               RED, new_dims, NC);
        return;
      }
      y->strides = new_strides;
    }
    // Overwrite value pointed to by `y->widths`
    y->widths = new_widths;
    // Recompute strides
    // Inner most dimension always has a stride of 1
    y->strides[new_dims-1] = 1;
    for (int i = new_dims-2; i >= 0; i--) {
      y->strides[i] = new_widths[i+1] * y->strides[i+1];
    }
  } else {
    printf("%sCannot reform to shape with different total number of elements: %d != %d%s\n",
           RED, total_existing_elems, total_proposed_elems, NC);
  }
}
// Resizing such that the total number of elements changes decreases is a more
// involved adjustment
// Eventually, it may be more time efficient to shrink by simply reinterpreting
// the existing array (eg: keep track of start and end and shift these values
// accordingly)
// TODO: This function contains old code meant to shrink and/or grow an array
void shrink_C_array(yarr *y, int *new_widths, int new_dims) {
  // Total number of elements in existing array
  int total_existing_elems = y->strides[0] * y->widths[0];
  
  // Total number of elements in resized array is the product of each width
  int total_proposed_elems = 1;
  for (int i = 0; i < new_dims; i++) {
    total_proposed_elems *= new_widths[i];
  }

  if (total_existing_elems > total_proposed_elems) {
    // If the number of elements is changed, it is necessary to realloc memory
    // If the array is shrinking, preserve values
    // If the array is growing, overwrite values
    realloc_contiguous(y, total_proposed_elems);
  } else {
    printf("%sCannot grow array with shrink operation: %d != %d%s\n",
           RED, total_existing_elems, total_proposed_elems, NC);
  }
}

/* Operations
 * DONE: array()
 * TODO: indexing of various kinds
 * TODO: matrixMult()
 * DONE: max() & min()
 * TODO: transpose()
 * DONE: primative pairwise operations (+, -, *, /, max, min)
 */

// TODO: Determine if casting should happen after operation is applied or
// before?
// Generate a new `yarr` with same shape as inputs and internal values computed
// using pairwise operators defined as `op`
yarr *apply_pairwise_op(yarr *y1, yarr *y2, pairwise_op op) {
  // Ensure that shapes are consistent
  // Start by ensuring number of dimensions is consistent
  if (y1->dims != y2->dims) {
    #ifdef DEBUG
    printf("%sInconsistent number of dimensions for pairwise operation%s\n",
           RED, NC);
    #endif
    return NULL;
  }
  
  int total_elems = y1->strides[0] * y1->widths[0];
  // Ensure the total number of elements is consistent
  if ((y2->strides[0] * y2->widths[0]) != total_elems) {
    #ifdef DEBUG
    printf("%sInconsistent number of total elements for pairwise operation%s\n",
           RED, NC);
    #endif
    return NULL;
  }

  // Array to return
  yarr *res =  malloc(1 * sizeof(yarr));
  // List of widths will have same number of dimensions as y1 and/or y2
  int *res_widths = malloc(y1->dims * sizeof(int));
  // List of strides will have same number of dimensions as y1 and/or y2
  int *res_strides = malloc(y1->dims * sizeof(int));
  
  // Perform operation unless inconsistent width is found
  for (int depth = 0; depth < y1->dims; depth++) {
    if (y1->widths[depth] != y2->widths[depth]) {
      #ifdef DEBUG
      printf("%sInconsistent numer of widths at depth %d for pairwise operation%s\n",
             RED, depth, NC);
      #endif
      // Free partially allocated `res`
      free(res);
      // Free widths list
      free(res_widths);
      return NULL;
    }
    // If width at current depth matches, update res_widths entry
    res_widths[depth] = y1->widths[depth];
    res_strides[depth] = y1->strides[depth];
  }

  // Set res->dims from y1->dims or y2->dims
  res->dims = y1->dims;
  // Set res->widths
  res->widths = res_widths;
  // Set res->strides
  res->strides = res_strides;
  // Set tag to be widest tag of (y1,y2)
  res->tag = widest_dt(y1->tag, y2->tag);
  // Initialize data
  allocate_contiguous(res, total_elems);
  
  // Perform pairwise operation
  for (int i = 0; i < total_elems; i++) {
    switch (res->tag) {
    case DOUBLE:
      if (op == ADD) {
        res->data.ddata[i] = get_double(y1,i) + get_double(y2,i);
      }
      else if (op == SUB) {
        res->data.ddata[i] = get_double(y1,i) - get_double(y2,i);
      }
      else if (op == MULT) {
        res->data.ddata[i] = get_double(y1,i) * get_double(y2,i);
      }
      else if (op == DIV) {
        res->data.ddata[i] = get_double(y1, i), get_double(y2, i);
      }
      else if (op == MAX) {
        res->data.ddata[i] = max(get_double(y1, i), get_double(y2, i));
      }
      else { res->data.ddata[i] = min(get_double(y1, i), get_double(y2, i)); }
      break;
    case LONG:
      if (op == ADD) {
        res->data.ldata[i] = get_long(y1,i) + get_long(y2,i);
      }
      else if (op == SUB) {
        res->data.ldata[i] = get_long(y1,i) - get_long(y2,i);
      }
      else if (op == MULT) {
        res->data.ldata[i] = get_long(y1,i) * get_long(y2,i);
      }
      else if (op == DIV) {
        res->data.ldata[i] = get_long(y1, i) / get_long(y2, i);
      }
      else if (op == MAX) {
        res->data.ldata[i] = max(get_long(y1, i), get_long(y2, i));
      }
      else { res->data.ldata[i] = min(get_long(y1, i), get_long(y2, i)); }
      break;
    case FLOAT:
      if (op == ADD) {
        res->data.fdata[i] = get_float(y1,i) + get_float(y2,i);
      }
      else if (op == SUB) {
        res->data.fdata[i] = get_float(y1,i) - get_float(y2,i);
      }
      else if (op == MULT) {
        res->data.fdata[i] = get_float(y1,i) * get_float(y2,i);
      }
      else if (op == DIV) {
        res->data.fdata[i] = get_float(y1, i) / get_float(y2, i);
      }
      else if (op == MAX) {
        res->data.fdata[i] = max(get_float(y1, i), get_float(y2, i));
      }
      else { res->data.fdata[i] = min(get_float(y1, i), get_float(y2, i)); }
      break;
    case INT:
      if (op == ADD) {
        res->data.idata[i] = get_int(y1,i) + get_int(y2,i);
      }
      else if (op == SUB) {
        res->data.idata[i] = get_int(y1,i) - get_int(y2,i);
      }
      else if (op == MULT) {
        res->data.idata[i] = get_int(y1,i) * get_int(y2,i);
      }
      else if (op == DIV) {
        res->data.idata[i] = get_int(y1, i) / get_int(y2, i);
      }
      else if (op == MAX) {
        res->data.idata[i] = max(get_int(y1, i), get_int(y2, i));
      }
      else { res->data.idata[i] = min(get_int(y1, i), get_int(y2, i)); }
      break;
    }
  }

  return res;
}

// Given a list of dimensions, return reference to widest set of values
// Arrays in Yorick are indexed from the inner-most dimension outwards
// eg: `a = array(1.0, 4,3,2)` creates a 2x3x4 (4 is the inner-most dimension)
// array where if we want to update the 2nd element in the 1st column of the
// second row, we would use:
// `a(2,1,2) = 7.0`
// visually this would look like:
// `[ [1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0] ]
//  [ [1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0] ]`
// becomes:
// `[ [1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0] ]
//  [ [1.0, 7.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0] ]`
// TODO: Implement for types beyond `int`
void *grab_point(yarr *y, int *dim_list, int dim_list_size) {
  int offset = 0;
  for (int i = 0; i < dim_list_size; i++) {    
    offset += (dim_list[i]-1) * y->strides[y->dims - (i+1)];
  }

  // Yorick allows indexing to cross dimensional boundries
  // eg: a = array(1, 4,3,2) is a 2x3x4 matrix
  // a(5,2,1) is interpreted to a(1,3,1) as the inner-most dimension is
  // only 4 elements wide
  // Ensure that offset is within bounds of contiguous array
  int total_elems = y->strides[0] * y->widths[0];
  if (offset > total_elems) {
    #ifdef DEBUG
    printf("%sAccess to element %d exceeds boundary of %d\n%s",
           RED, offset, total_elems, NC);
    #endif
    return NULL;
  }
  
  #ifdef DEBUG
  char access[dim_list_size * sizeof(int)];
  int d_list_pos = 0;
  for (int i = 0; i < dim_list_size; i++) {
    d_list_pos += sprintf(&access[d_list_pos], "%d ", dim_list[i]);
  }
  printf("%sOffset is measured as %d for access to ( %s)%s\n",
         YELLOW, offset, access, NC);
  #endif
  
  return (int *)(y->data.idata + offset);
}

int main() {
  int dims = 3;
  int *widths = malloc(dims * sizeof(int));
  widths[0] = 2;
  widths[1] = 3;
  widths[2] = 4;
  yarr *y = C_array(1.0, INT, widths, dims);
  print_C_array(y);

  #ifdef DEBUG
  printf("%supdating %d element%s\n",
         YELLOW,
         ((2*y->widths[0])+(1*y->widths[1])+(2*y->widths[2])),
         NC);
  #endif
  
  y->data.idata[2] = 7;
  // Simulate a(2,1,2) = 8
  int *make_8 = ((int *)y->data.idata +
    (2-1)*y->strides[0] +
    (1-1)*y->strides[1] +
    (2-1)*y->strides[2]);
  *make_8 = 8;

  // Simulate a(3,3,1) = 4
  int *make_4 = ((int *)y->data.idata +
    (1-1)*y->strides[0] +
    (3-1)*y->strides[1] +
    (3-1)*y->strides[2]);
  *make_4 = 4;

  int dim_list[3] = {3,3,1};
  int *make_7 = grab_point(y, dim_list, 3);
  *make_7 = 2;

  int dim_list2[3] = {5,2,1};
  int *make_0 = grab_point(y, dim_list2, 3);
  *make_0 = 0;
  
  #ifdef DEBUG
  printf("\n");
  print_yarr(y);
  printf("\n");
  #endif
  
  print_C_array(y);

  // Test out addition
  yarr *augend = C_array(12.0, INT, widths, dims);
  yarr *addend = C_array(3.0, INT, widths, dims);
  yarr *sum = apply_pairwise_op(augend, addend, ADD);

  print_C_array(sum);
  
  // Test out min and max
  #ifdef DEBUG
  printf("%sPrinting out min of two arrays:%s\n", YELLOW, NC);
  #endif
  yarr *mins = apply_pairwise_op(y, addend, MIN);
  print_C_array(mins);

  #ifdef DEBUG
  printf("%sPrinting out max of two arrays:%s\n", YELLOW, NC);
  #endif
  yarr *maxs = apply_pairwise_op(y, addend, MAX);
  print_C_array(maxs);

  #ifdef DEBUG
  printf("%sReform an array%s\n", YELLOW, NC);
  printf("%sOriginal:%s\n", YELLOW, NC);
  #endif
  // Total number of elements in y is 24 (2*3*4)
  print_C_array(y);
  #ifdef DEBUG
  printf("%s Reformed array:%s\n", YELLOW, NC);
  #endif
  int reformed_dims = 3;
  int reformed_widths[3] = {2, 2, 6};
  reform_C_array(y, reformed_widths, reformed_dims);
  print_C_array(y);

  #ifdef DEBUG
  printf("%sReform array again:%s\n", YELLOW, NC);
  #endif
  reformed_dims = 2;
  int another_reformed_widths[2] = {2, 12};
  reform_C_array(y, another_reformed_widths, reformed_dims);
  print_C_array(y);
  
  return 0;
}
