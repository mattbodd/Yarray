#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "yoperations.h"

// Macro definitions
#define max(a,b) ({__typeof__ (a) _a = (a);     \
      __typeof__ (b) _b = (b);                  \
      _a > _b ? _a : _b; })

#define min(a,b) ({__typeof__ (a) _a = (a);     \
      __typeof__ (b) _b = (b);                  \
      _a < _b ? _a : _b; })

// TODO: decide if `source_index` is necessary as `fill_val` may always be a single
// element `yarray`
void update_point(yarr *y, int dest_index, yarr *fill_val, int source_index) {
  if (y->tag == INT) {
    y->data.idata[dest_index] = fill_val->data.idata[source_index];
  } else if (y->tag == LONG) {
    y->data.ldata[dest_index] = fill_val->data.ldata[source_index];
  } else if (y->tag == FLOAT) {
    y->data.fdata[dest_index] = fill_val->data.fdata[source_index];
  } else if (y->tag == DOUBLE) {
    y->data.ddata[dest_index] = fill_val->data.ddata[source_index];
  } else {
    printf(RED"Invalid dataType recieved in `update_point\n`"NC);
  }
}

// Update elements in a `yarray` given bounds for each specified dimension
// Bounds are taken as a two tuple
// When a dimension is meant to be an observer (eg: a(:,...)), its bounds will
// be passed as (0,-1)
// TODO: determine if `bounds_size` will always be equivalent to `y->dims`
// TODO: expand `bounds` to be a three tuple: (low, high, step)
void update_C_array(yarr *y, yarr *fill_vals, int bounds_size, int *bounds) {  
  // Keep track of index within each dimension
  int dims_index[bounds_size/2];
  memset(dims_index, 0, (bounds_size/2) * sizeof(int));
  // Initialize index into each dim as the lower bound for a given dim
  for (int dim = 0; dim < bounds_size/2; dim++) {
    dims_index[dim] = bounds[dim*2];
  }
  
  // Determine the number of elements to update past the lower bound
  int widths[bounds_size/2];
  memset(widths, 0, (bounds_size/2) * sizeof(int));
  int upper = 0;
  for (int dim = 0; dim < (bounds_size/2); dim++) {
    upper = bounds[(dim*2)+1];
    // If dimension is meant to spectated over
    if (upper == -1) {
      upper = y->widths[dim] - 1;
      // Overwrite -1 to be correct upper boundry
      bounds[(dim*2)+1] = upper;
    }
    // Ensure start <= stop
    if (dims_index[dim] > upper) {
      printf(YELLOW"Invalid range in `update_C_array`: %d..%d"NC,
             dims_index[dim], upper);
      return;
    }
    // Calculate width
    widths[dim] = (upper - dims_index[dim]) + 1;
  }

  int total_updates = 1;
  for (int dim = 0; dim < (bounds_size/2); dim++) {
    total_updates *= widths[dim];
  }

  // Determine if `fill_vals` is a single value or some number of elements
  // which is equal to the number of elements being updated
  int total_fill_vals = fill_vals->strides[0] * fill_vals->widths[0];
  if (total_fill_vals != 1) {
    // Ensure that the number of values to update is equal to the number of
    // fill values
    if (total_fill_vals != total_updates) {
      printf(RED"Cannot update %d elements with %d values\n"NC,
             total_updates, total_fill_vals);
      return NULL;
    }
  }

  #ifdef DEBUG
  printf(YELLOW"Updating %d elements, with widths:"NC, total_updates);
  for (int i = 0; i < (bounds_size/2); i++) {
    printf(YELLOW"%d "NC, widths[i]);
  }
  printf("\n");
  printf(YELLOW"With dims_index:"NC);
  for (int i = 0; i < (bounds_size/2); i++) {
    printf(YELLOW"%d "NC, dims_index[i]);
  }
  printf("\n");
  #endif
  
  // Index into contiguous representation
  int index = 0;
  // Index to track which fill value to use
  int fill_index = 0;
  // Start with inner most dimension
  int curr_dim;
  
  for (int i = 0; i < total_updates; i++) {
    // Reset curr_dim to point to the inner most dimension
    curr_dim = (bounds_size/2)-1;

    // Reset index
    index = 0;
    for (int dim = 0; dim < (bounds_size/2); dim++) {
      index += y->strides[dim]*dims_index[dim];
    }
    // Perform actual update
    // `get_element` will always return a single elemenet `yarray`
    update_point(y, index, get_element(fill_vals, fill_index), 0);
    // Conditionally update `fill_index`
    if (total_fill_vals > 1) {
      fill_index++;
    }

    // Update the `dim_index` in the inner most dimension
    ++dims_index[curr_dim];
    // Increment higher dimensions if necessary
    while (curr_dim >= 0) {
      // If the current index exceeds the upper bound
      if (dims_index[curr_dim] > bounds[(curr_dim*2) + 1]) {
        dims_index[curr_dim] = bounds[curr_dim*2];
        // Increment the index of the next highest dimension
        dims_index[--curr_dim]++;
      } else {
        // Stop incrementing dimension indices
        break;
      }
    }
  }
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
      int *new_strides = realloc(y->strides, new_dims*sizeof(int));
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
// TODO: Handle negative indexing
void shrink_C_array(yarr *y, int *new_widths, int new_dims, int start, int stop) {
  // Total number of elements in existing array
  int total_existing_elems = y->strides[0] * y->widths[0];
  
  // Total number of elements in resized array is the product of each width
  int total_proposed_elems = 1;
  for (int i = 0; i < new_dims; i++) {
    total_proposed_elems *= new_widths[i];
  }

  if ((stop-start) != total_proposed_elems) {
    printf("%sValues to preserve do not match total number of proposed values%s\n",
           RED, NC);
  }

  // Ensure that values start..stop occupy the first (stop-start) values in
  // the contiguous representation so that reallocation preserves necessary
  // values
  for (int i = 0; i < total_proposed_elems; i++) {
    move_elem(y, start+i, i);
  }

  #ifdef DEBUG
  printf("%s%p[ ", YELLOW, &y->data.idata);
  for (int i = 0; i < total_proposed_elems; i++) {
    printf("%d(%d) ", y->data.idata[i], i);
  }
  printf("]%s\n", NC);
  #endif
  
  if (total_existing_elems > total_proposed_elems) {
    // If the number of elements is changed, it is necessary to realloc memory
    // If the array is shrinking, preserve values
    // If the array is growing, overwrite values
    realloc_contiguous(y, total_proposed_elems);
    // Update dims and widths
    if (new_dims != y->dims) {
      y->dims = new_dims;

      // Try to make use of previously allocated memory
      int *new_strides = realloc(y->strides, new_dims*sizeof(int));
      if (new_strides == NULL) {
        printf("%sCould not reallocate for new strides of size %d%s\n",
               RED, new_dims, NC);
        return;
      }
      y->strides = new_strides;
    }

    #ifdef DEBUG
    printf("%s%p[ ", YELLOW, &y->data.idata);
    for (int i = 0; i < total_proposed_elems; i++) {
      printf("%d(%d) ", y->data.idata[i], i);
    }
    printf("]%s\n", NC);
    #endif

    y->widths = new_widths;
    
    // Recompute strides
    // Inner most dimension always has a stride of 1
    y->strides[new_dims-1] = 1;
    for (int i = new_dims-2; i >= 0; i--) {
      y->strides[i] = new_widths[i+1] * y->strides[i+1];
    }
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

// Operation on two-dimensional array where relationship between sizes is:
// (m*n)(n*k)
yarr *matrix_mul(yarr *multiplicand, yarr *multiplier) {
  if (multiplicand->widths[1] != multiplier->widths[0]) {
    printf(RED"Matrix dimensions are incompatible:"NC);
    return NULL;
  }
}


yarr *broadcast(yarr *y, yarr *val, Op op) {
  yarr *copy = copy_yarray(y);
  
  int total_elems = y->strides[0] * y->widths[0];
  
  for (int i = 0; i < total_elems; i++) {
    apply_bin_op(copy, i, op, get_element(y, i), get_element(val, 0));
  }

  return copy;
}

// DONE: Determine if casting should happen after operation is applied or
// before? - After
// Generate a new `yarr` with same shape as inputs and internal values computed
// using pairwise operators defined as `op`
yarr *apply_Op(yarr *y1, yarr *y2, Op op) {
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
    apply_bin_op(res, i, op, get_element(y1, i), get_element(y2, i));
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

/*
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
  yarr *sum = apply_Op(augend, addend, ADD);

  print_C_array(sum);
  
  // Test out min and max
  #ifdef DEBUG
  printf("%sPrinting out min of two arrays:%s\n", YELLOW, NC);
  #endif
  yarr *mins = apply_Op(y, addend, MIN);
  print_C_array(mins);
  // Free unused memory
  //free(mins);

  #ifdef DEBUG
  printf("%sPrinting out max of two arrays:%s\n", YELLOW, NC);
  #endif
  yarr *maxs = apply_Op(y, addend, MAX);
  print_C_array(maxs);
  // Free unused memory
  //free(maxs);

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

  #ifdef DEBUG
  printf("%sShrink an array%s\n", YELLOW, NC);
  printf(YELLOW"Original:\n"NC);
  #endif
  print_C_array(y);
  #ifdef DEBUG
  printf(YELLOW"Shrunken array:\n"NC);
  #endif
  int shrunken_dims = 2;
  int shrunken_widths[2] = {2, 5};
  shrink_C_array(y, shrunken_widths, shrunken_dims, 2, 12);
  print_C_array(y);
  print_yarr(y);

  printf("\n");

  // Test alternate invocation of `_C_array`
  yarr *alt_y = _C_array(INT, 0, 2, (int []){4, 3});
  print_yarr(alt_y);
  print_C_array(alt_y);

  // Perform update on array
  printf("About to update array\n");
  update_C_array(alt_y, 1.0, 4, (int []){2,-1 , 0,-1});
  print_C_array(alt_y);

  // Perform additional update on array
  printf("About to update array\n");
  update_C_array(alt_y, 2.0, 4, (int []){2,2 , 1,2});
  print_C_array(alt_y);

  // Create yarray from slice
  printf("Creating slice\n");
  yarr *slice = get_slice(alt_y, 4, (int[]){2,3 , 1,2});
  print_C_array(slice);

  // Three dimensional yarray examples of `update_C_array` and `get_slice`
  // Initialize three dimensional array
  printf(YELLOW"Three dimensional array initialization\n"NC);  
  yarr *three_dim_ex = _C_array(INT, 0, 3, (int []){5,4,3});
  print_C_array(three_dim_ex);
  // Update values in three dimensional array
  printf(YELLOW"Update three dimensional array\n"NC);
  update_C_array(three_dim_ex, 1.0, 6, (int []){2,2 , 1,2 , 0,-1});
  print_C_array(three_dim_ex);
  // Create slice of array and store into new array
  printf(YELLOW"Create three dimensional slice\n"NC);
  yarr *three_dim_slice = get_slice(three_dim_ex, 6, (int []){2,2 , 1,2 , 0,-1});
  print_C_array(three_dim_slice);
  print_yarr(three_dim_slice);
  
  // Free unused memory
  dealloc_yarr(augend);
  dealloc_yarr(addend);
  dealloc_yarr(sum);
  free(y);
  free(widths);
  
  return 0;
}
*/
