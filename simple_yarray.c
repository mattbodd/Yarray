#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// DEBUG information
#define YELLOW "\033[0;33m"
#define RED    "\033[0;31m"
#define NC     "\033[0m"

// Possible primative data types of `yarr.data`
typedef enum {INT, FLOAT, LONG, DOUBLE} dataType;
// Possible scalar operators
typedef enum {ADD, SUB, MULT, DIV} scalar_op;

// Forward declaration of yarr type
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

// Function to compute the widest dataType between two yarrs
dataType widest_dt(dataType y1_tag, dataType y2_tag) {
  // Double is the widest supported dataType
  if      (y1_tag == DOUBLE || y2_tag == DOUBLE) { return DOUBLE; }
  else if (y1_tag == LONG   || y2_tag == LONG)   { return LONG; }
  else if (y1_tag == FLOAT  || y2_tag == FLOAT)  { return FLOAT; }
  else                                           { return INT; }
}

double get_double(yarr *y, int index) {
  switch (y->tag) {
  case DOUBLE:
    return y->data.ddata[index];
    break;
  case LONG:
    return (double)y->data.ldata[index];
    break;
  case FLOAT:
    return (double)y->data.fdata[index];    
    break;
  case INT:
    return (double)y->data.idata[index];    
    break;
  }
}

double get_long(yarr *y, int index) {
  switch (y->tag) {
  case DOUBLE:
    return (long)y->data.ddata[index];
    break;
  case LONG:
    return y->data.ldata[index];
    break;
  case FLOAT:
    return (long)y->data.fdata[index];    
    break;
  case INT:
    return (long)y->data.idata[index];    
    break;
  }
}

double get_float(yarr *y, int index) {
  switch (y->tag) {
  case DOUBLE:
    return (float)y->data.ddata[index];
    break;
  case LONG:
    return (float)y->data.ldata[index];
    break;
  case FLOAT:
    return y->data.fdata[index];    
    break;
  case INT:
    return (float)y->data.idata[index];    
    break;
  }
}

double get_int(yarr *y, int index) {
  switch (y->tag) {
  case DOUBLE:
    return (int)y->data.ddata[index];
    break;
  case LONG:
    return (int)y->data.ldata[index];
    break;
  case FLOAT:
    return (int)y->data.fdata[index];    
    break;
  case INT:
    return y->data.idata[index];    
    break;
  }
}

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

void fill_array(yarr *y, int total_elems, double fill_val) {
  switch (y->tag) {
  case INT:
    // Initialize int array
    #ifdef DEBUG
    printf("%sInitializing int array with %d elements\n%s",
	   YELLOW, total_elems, NC);
    #endif
    for (int i = 0; i < total_elems; i++) {
      y->data.idata[i] = (int) fill_val;
    }
    #ifdef DEBUG
    printf("%sDone initializing int array with %d elements%s\n",
	   YELLOW, total_elems, NC);
    #endif
    
    break;
  case FLOAT:
    // Initialize float array
    for (int i = 0; i < total_elems; i++) {
      y->data.fdata[i] = (float) fill_val;
    }
    #ifdef DEBUG
    printf("%sDone initializing float array with %d elements%s\n",
	   YELLOW, total_elems, NC);
    #endif
    
    break;
  case LONG:
    // Initialize long array
    for (int i = 0; i < total_elems; i++) {
      y->data.ldata[i] = (long) fill_val;
    }
    #ifdef DEBUG
    printf("%sDone initializing long array with %d elements%s\n",
	   YELLOW, total_elems, NC);
    #endif
    
    break;
  case DOUBLE:
    // Initialize double array
    for (int i = 0; i < total_elems; i++) {
      y->data.ddata[i] = fill_val;
    }
    #ifdef DEBUG
    printf("%sDone initializing double array with %d elements%s\n",
	   YELLOW, total_elems, NC);
    #endif
    
    break;
  default:
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

// TODO: implement resize operation
// Resizing an array in such a way that it maintains the total number of elements requires
// readjusting strides
// Resizing such that the total number of elements changes requires a more involved adjustment
void resize_C_array(yarr *y) {}

/* Operations
 * DONE: array()
 * TODO: indexing of various kinds
 * TODO: matrixMult()
 * DONE: max() & min()
 * TODO: transpose()
 * DONE: primative scalar operations (+, -, *, /)
 */

// TODO: Determine if casting should happen after operation is applied or
// before?
yarr *apply_scalar_op(yarr *y1, yarr *y2, scalar_op op) {
  // Ensure that shapes are consistent
  // Start by ensuring number of dimensions is consistent
  if (y1->dims != y2->dims) {
    #ifdef DEBUG
    printf("%sInconsistent number of dimensions for scalar operation%s\n",
           RED, NC);
    #endif
    return NULL;
  }
  
  int total_elems = y1->strides[0] * y1->widths[0];
  // Ensure the total number of elements is consistent
  if ((y2->strides[0] * y2->widths[0]) != total_elems) {
    #ifdef DEBUG
    printf("%sInconsistent number of total elements for scalar operation%s\n",
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
      printf("%sInconsistent numer of widths at depth %d for scalar operation%s\n",
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
  
  // Perform primative operation
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
      else { res->data.ddata[i] = get_double(y1, i), get_double(y2, i); }
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
      else { res->data.ldata[i] = get_long(y1, i) / get_long(y2, i); }
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
      else { res->data.fdata[i] = get_float(y1, i) / get_float(y2, i); }
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
      else { res->data.idata[i] = get_int(y1, i) / get_int(y2, i); }
      break;
    }
  }

  return res;
}

void print_indent(int level) {
  // Indentation is a double space
  for (int i = 0; i < level*2; i++) { printf(" "); }
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
    return NULL;
    #endif
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

void print_C_array_helper(yarr *y, int depth, int offset) {
  for (int i = 0; i < y->widths[depth]; i++) {
    if (depth < (y->dims-1)) {
      // Indent line
      print_indent(depth);
      // Add array indicator
      printf("{");
      // If not the second-most inner dimension, print next dimension on newline
      if (depth+1 < y->dims-1) { printf("\n"); }

      // Recursively print next dimension
      print_C_array_helper(y, depth+1, offset + ((i)*y->strides[depth]));

      // If not the second-most inner dimension, add proper newline indentation
      if (depth+1 < y->dims-1) { print_indent(depth); }
      printf("}\n");
    } else {
      printf("%d", y->data.idata[offset+i]);
      // Append ", " for each inner-most element that is not the final element
      if (i < y->widths[depth]-1) { printf(", "); }
    }
  }
}


void print_C_array(yarr *y) {
  print_C_array_helper(y, 0, 0);
}

void print_yarr(yarr *y) {
  char widths[sizeof(int) * y->dims + 1];
  char strides[sizeof(int) * y->dims + 1];
  // Longest tag name is 6 chars ("double")
  char tag[sizeof(char) * 6 + 1];
  char data[sizeof(int) * (y->strides[0]*y->widths[0])];
  int w_pos = 0;
  int s_pos = 0;
  int d_pos = 0;

  // Populate widths and strides strings
  for (int i = 0; i < y->dims; i++) {
    w_pos += sprintf(&widths[w_pos], "%d ", y->widths[i]);
    s_pos += sprintf(&strides[s_pos], "%d ", y->strides[i]);
  }

  // Populate tag string
  switch (y->tag) {
  case INT:
    sprintf(tag, "int");
    break;
  case FLOAT:
    sprintf(tag, "float");
    break;
  case LONG:
    sprintf(tag, "long");
    break;
  case DOUBLE:
    sprintf(tag, "double");
    break;
  }

  // Populate data string
  for (int i = 0; i < y->widths[0]*y->strides[0]; i++) {
    d_pos += sprintf(&data[d_pos], "%d ", y->data.idata[i]);
  }
  
  printf("%sdims:      %d\n"
         "widths:    { %s}\n"
         "strides:   { %s}\n"
         "tag:       %s\n"
         "data:      %s%s\n",
         YELLOW,
         y->dims, widths, strides, tag, data,
         NC
         );
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
  yarr *sum = apply_scalar_op(augend, addend, ADD);

  print_C_array(sum);
  
  return 0;
}
