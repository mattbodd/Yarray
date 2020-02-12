#include <stdlib.h>
#include <stdio.h>
#include "yarray.h"
#include "yutilities.h"

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