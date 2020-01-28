#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// internals of yarr
// dataType is used to signal the type of values in the inner most dimension of yarr
typedef enum {INT, FLOAT, LONG, YARR} dataType;
// data represents a single dimension of a yarr object

// Forward declaration of yarr type
typedef struct yarr yarr;

// yarr is the overarching array type used to hold all values
// even single values will be stored using yarr
// eg: x = 1; -> info,x = array(long,1)
struct yarr {
  int length;
  int dim;
  dataType tag;
  union data {
    int *idata;
    float *fdata;
    long *ldata;
    yarr *ydata;
  } data;
};

/* Arugments:
 * `fill_val`: value to set inner-most elements to (default as float to allow for truncation)
 * `lengths*:  list of length for each dimension starting with innermost dimension
 * `dims`:     number of dimensions (size of lengths)
 */
/* Example usage:
 * Invocation: C_array(0, [3, 2, 1], 3)
 * Result:
 * [ 
 *   [ 
 *     [0], [0], [0]
 *   ],
 *   [
 *     [0], [0], [0]
 *   ]
 * ]
 * Recursive calls: C_array(0, [3, 2, 1], 3)
 *                  C_array(0, [3, 2, 1], 2)
 *                  C_array(0, [3, 2, 1], 1)
 */
yarr *C_array(float fill_val, dataType tag, int *lengths, int dims) {
  // Create pointer to yarr struct
  // Yarr type is always used to encapsulate either another set of Yarrs or a set
  // of primative values
  yarr *y;

  // Allocate space for the dim-th dimension
  y = malloc(sizeof(yarr) * lengths[dims-1]);
  // handle case where new Yarr pointer cannot be allocated
  if (y == 0) {
    return NULL;
  }

  // DEBUG
  printf("Dim: %d\n", dims);
  // GUBED

  // Set type field
  // Outermost Yarrs may encapsulate other Yarr's which ultimately hold primative data
  // types, however, the innermost type is what is displayed pervasively
  y->tag = tag;
  // Set current dimension
  y->dim = dims;
  // Set length field
  y->length = lengths[dims-1];
  
  // Outer dimensions will be of type YARR
  if (dims > 1) {
    // For each inner dimmension, allocate a Yarr to hold data
    for (int i = 0; i < y->length; i++) {
      y->data.ydata = malloc(y->length * sizeof(yarr));
      // Handle case where new Yarr pointer cannot be allocated
      if (y->data.ydata == 0) {
	free(y);
	return NULL;
      }
      // Recursively 'fill' inner dimensions
      y->data.ydata = C_array(fill_val, tag, lengths, dims-1);
    }
  } else {
    // Allocate space for inner-most dimension
    switch (tag) {
    case INT:
      // allocate memory for data
      y->data.idata = malloc(y->length * sizeof(int));
      // handle case where data pointer cannot be allocated
      if (y->data.idata == 0) {
	free(y);
	return NULL;
      }
      // initialize data
      for (int i = 0; i < y->length; i++) {
	y->data.idata[i] = fill_val;
      }
      break;
    case FLOAT:
      // allocate memory for data
      y->data.fdata = malloc(y->length * sizeof(float));
      // handle case where data pointer cannot be allocated
      if (y->data.idata == 0) {
	free(y);
	return NULL;
      }
      // initialize data
      for (int i = 0; i < y->length; i++) {
	y->data.fdata[i] = fill_val;
      }
      break;
    case LONG:
      // allocate memory for data
      y->data.ldata = malloc(y->length * sizeof(long));
      // handle case where data pointer cannot be allocated
      if (y->data.idata == 0) {
	free(y);
	return NULL;
      }
      // initialize data
      for (int i = 0; i < y->length; i++) {
	y->data.ldata[i] = fill_val;
      }
      break;
    default:
      printf("Invalid dataType received\n");
      return NULL;
    }
  }
  
  return y;
}

void free_C_array(yarr *y) {
  if (y->dim > 1) {
    free_C_array(y->data.ydata);
    printf("Freeing %d items from dim %d\n", y->length, y->dim);
    free(y);
  } else {
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
    }
    printf("Freeing %d items from dim %d\n", y->length, y->dim);
    free(y);
  }
}

void print_yarr_helper(yarr *y, int indentation) {
  // Allocate buffer to store indentation string
  char indent_str[2*indentation];
  // Initialize string buffer
  sprintf(indent_str, "");
  for (int i = indentation; i > 0; i--) {
    // Indentations are printed as two spaces
    sprintf(indent_str + strlen(indent_str), "  ");
  }
  
  // Observing outer dimension
  if (y->dim > 1) {
    // DEBUG
    //printf("In outer dimension with length: %d\n", y->length);
    // GUBED
    
    printf("%s[\n", indent_str);
    
    // Recurisvely display inner dimensions
    for (int inner_y = 0; inner_y < y->length; inner_y++) {
      // DEBUG
      //printf("Recusively calling print_yarr_helper\n");
      // GUBED
      print_yarr_helper(y->data.ydata, indentation+1); 
    }
    
    printf("%s]\n", indent_str);
  } else {
    // DEBUG
    //printf("In inner dimension\n");
    // GUBED
    
    // Pad with proper indentation
    printf("%s[ ", indent_str);
    // Innermost dimension has been reached
    for (int inner_data = 0; inner_data < y->length; inner_data++) {
      switch (y->tag) {
      case INT:
	printf("%d ", y->data.idata[inner_data]);
	break;
      case FLOAT:
	printf("%f ", y->data.fdata[inner_data]);
	break;
      case LONG:
	printf("%lu ", y->data.ldata[inner_data]);
	break;
      }
    }
    printf("%s", "]\n");
  }
}

void print_yarr(yarr *y) {
  print_yarr_helper(y, 0);
}

int main(int argc, char **argv) {
  int lengths[3] = {2, 2, 3};
  yarr *y = C_array(5, INT, lengths, 3);
  
  print_yarr(y);

  free_C_array(y);

  return 0;
}
