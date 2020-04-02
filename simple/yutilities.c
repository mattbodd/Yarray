#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "yarray.h"
#include "yutilities.h"

// Macro definitions
#define max(a,b) ({__typeof__ (a) _a = (a);     \
      __typeof__ (b) _b = (b);                  \
      _a > _b ? _a : _b; })

#define min(a,b) ({__typeof__ (a) _a = (a);     \
      __typeof__ (b) _b = (b);                  \
      _a < _b ? _a : _b; })

// Function to compute the widest dataType between two yarrs
dataType widest_dt(dataType y1_tag, dataType y2_tag) {
  // Double is the widest supported dataType
  if      (y1_tag == DOUBLE || y2_tag == DOUBLE) { return DOUBLE; }
  else if (y1_tag == LONG   || y2_tag == LONG)   { return LONG; }
  else if (y1_tag == FLOAT  || y2_tag == FLOAT)  { return FLOAT; }
  else                                           { return INT; }
}

// Safely assume `val` can be cast to the type of `y` as this function does not
// imply any changes to underlying dataType of `y`
void update_yarr(yarr *y, int index, yarr *val) {
  if (y->tag == DOUBLE) {
    y->data.ddata[index] = val->data.ddata[0];
  } else if (y->tag == LONG) {
    y->data.ldata[index] = val->data.ldata[0];
  } else if (y->tag == FLOAT) {
    y->data.fdata[index] = val->data.fdata[0];
  } else if (y->tag == INT) {
    y->data.idata[index] = val->data.idata[0];
  } else {
    printf("ERROR: unexpected dataType in `update_yarr`\n");
  }
}

// Safely assume the widest dataType between `arg_1` and `arg_2` is the `tag`
// of `y`
// In the case where `y` is a `yarr` containing `arg_1`, and `arg_2` is a
// non-associated value (is not contained in a `yarr`), the dataType of `y`
// cannot be changed using this function
// In the cahse where `y` is a new `yarr` meant to hold the results of a pairwise
// operation on `arg_1` and `arg_2` which are both `yarr`, `y` will have already
// been initialized with the widest dataType as its tag form `arg_1` and `arg_2`
void apply_bin_op(yarr *y, int index, Op op, yarr *arg_1, yarr *arg_2) {
  if (y->tag == DOUBLE) {
    double primitive_res = 0.0;
    if (op == ADD) {
      primitive_res = (arg_1->data.ddata[0] + arg_2->data.ddata[0]);
    } else if (op == SUB) {
      primitive_res = (arg_1->data.ddata[0] - arg_2->data.ddata[0]);
    } else if (op == MUL) {
      primitive_res = (arg_1->data.ddata[0] * arg_2->data.ddata[0]);
    } else if (op == DIV) {
      primitive_res = (arg_1->data.ddata[0] / arg_2->data.ddata[0]);
    } else if (op == MAX) {
      primitive_res = max(arg_1->data.ddata[0], arg_2->data.ddata[0]);
    } else if (op == MIN) {
      primitive_res = min(arg_1->data.ddata[0], arg_2->data.ddata[0]);
    } else if (op == POW) {
      primitive_res = (arg_1->data.ddata[0] * arg_2->data.ddata[0]);
    } else {
      printf("ERROR: unexpected operation in `apply_bin_op`\n");
      return;
    }
    yarr *yarr_res = _C_array(DOUBLE, primitive_res, 1, (int []){1});
    // Update value
    update_yarr(y, index, yarr_res);
  } else if (y->tag == LONG) {
    long primitive_res = 0.0;
    if (op == ADD) {
      primitive_res = (arg_1->data.ldata[0] + arg_2->data.ldata[0]);
    } else if (op == SUB) {
      primitive_res = (arg_1->data.ldata[0] - arg_2->data.ldata[0]);
    } else if (op == MUL) {
      primitive_res = (arg_1->data.ldata[0] * arg_2->data.ldata[0]);
    } else if (op == DIV) {
      primitive_res = (arg_1->data.ldata[0] / arg_2->data.ldata[0]);
    } else if (op == MAX) {
      primitive_res = max(arg_1->data.ldata[0], arg_2->data.ldata[0]);
    } else if (op == MIN) {
      primitive_res = min(arg_1->data.ldata[0], arg_2->data.ldata[0]);
    } else if (op == POW) {
      primitive_res = (arg_1->data.ldata[0] * arg_2->data.ldata[0]);
    } else {
      printf("ERROR: unexpected operation in `apply_bin_op`\n");
      return;
    }
    yarr *yarr_res = _C_array(LONG, primitive_res, 1, (int []){1});
    // Update value
    update_yarr(y, index, yarr_res);
  } else if (y->tag == FLOAT) {
    float primitive_res = 0.0;
    if (op == ADD) {
      primitive_res = (arg_1->data.fdata[0] + arg_2->data.fdata[0]);
    } else if (op == SUB) {
      primitive_res = (arg_1->data.fdata[0] - arg_2->data.fdata[0]);
    } else if (op == MUL) {
      primitive_res = (arg_1->data.fdata[0] * arg_2->data.fdata[0]);
    } else if (op == DIV) {
      primitive_res = (arg_1->data.fdata[0] / arg_2->data.fdata[0]);
    } else if (op == MAX) {
      primitive_res = max(arg_1->data.fdata[0], arg_2->data.fdata[0]);
    } else if (op == MIN) {
      primitive_res = min(arg_1->data.fdata[0], arg_2->data.fdata[0]);
    } else if (op == POW) {
      primitive_res = (arg_1->data.fdata[0] * arg_2->data.fdata[0]);
    } else {
      printf("ERROR: unexpected operation in `apply_bin_op`\n");
      return;
    }
    yarr *yarr_res = _C_array(FLOAT, primitive_res, 1, (int []){1});
    // Update value
    update_yarr(y, index, yarr_res);
  } else if (y->tag == INT) {
    int primitive_res = 0;
    if (op == ADD) {
      primitive_res = (arg_1->data.idata[0] + arg_2->data.idata[0]);
    } else if (op == SUB) {
      primitive_res = (arg_1->data.idata[0] - arg_2->data.idata[0]);
    } else if (op == MUL) {
      primitive_res = (arg_1->data.idata[0] * arg_2->data.idata[0]);
    } else if (op == DIV) {
      primitive_res = (arg_1->data.idata[0] / arg_2->data.idata[0]);
    } else if (op == MAX) {
      primitive_res = max(arg_1->data.idata[0], arg_2->data.idata[0]);
    } else if (op == MIN) {
      primitive_res = min(arg_1->data.idata[0], arg_2->data.idata[0]);
    } else if (op == POW) {
      primitive_res = (arg_1->data.idata[0] * arg_2->data.idata[0]);
    } else {
      printf("ERROR: unexpected operation in `apply_bin_op`\n");
      return;
    }
    yarr *yarr_res = _C_array(INT, primitive_res, 1, (int []){1});
    // Update value
    update_yarr(y, index, yarr_res);
  } else {
    printf("ERROR: unexpected dataType in `apply_bin_op`\n");
    return;
  }
}

// Create and return a single member `yarr` containing the relevant value
yarr *get_element(yarr *y, int index) {
  // Create the single value `yarr`
  yarr *res;
  if (y->tag == DOUBLE) {
    res = _C_array(DOUBLE, y->data.ddata[index], 1, (int []){1});
  } else if (y->tag == LONG) {
    res = _C_array(LONG, y->data.ldata[index], 1, (int []){1});
  } else if (y->tag == FLOAT) {
    res = _C_array(FLOAT, y->data.fdata[index], 1, (int []){1});
  } else if (y->tag == INT) {
    res = _C_array(INT, y->data.idata[index], 1, (int []){1});
  } else {
    printf("ERROR: unexpected data type in yarray\n");
    return NULL;
  }

  return res;
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
  default:
    printf(RED"Invalid dataType detected in `get_double`"NC);
    return -1.0;
    break;
  }
}

long get_long(yarr *y, int index) {
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
  default:
    printf(RED"Invalid dataType detected in `get_long`"NC);
    break;
  }
}

float get_float(yarr *y, int index) {
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
  default:
    printf(RED"Invalid dataType detected in `get_float`"NC);
  }
}

int get_int(yarr *y, int index) {
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
  default:
    printf(RED"Invalid dataType detected in `get_int`"NC);
    break;
  }
}

void move_elem(yarr *y, int source, int dest) {
  switch (y->tag) {
  case INT:
    #ifdef DEBUG
    printf("%sSwapping %d(%d) with %d(%d)%s | ", YELLOW,
           y->data.idata[dest], dest,
           y->data.idata[source], source, NC);
    #endif
    y->data.idata[dest] = y->data.idata[source];
    #ifdef DEBUG
    printf("%sCheck after: %d%s\n", YELLOW, y->data.idata[dest], NC);
    #endif
    break;
  case FLOAT:
    y->data.fdata[dest] = y->data.fdata[source];
    break;
  case LONG:
    y->data.ldata[dest] = y->data.ldata[source];
    break;
  case DOUBLE:
    y->data.ddata[dest] = y->data.ddata[source];
    break;
  default:
    printf("%sInvalid type associated with y->tag%s\n", RED, NC);
    break;
  }
}

void print_indent(int level) {
  // Indentation is a double space
  for (int i = 0; i < level*2; i++) { printf(" "); }
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
