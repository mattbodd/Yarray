#include <stdlib.h>
#include <stdio.h>
#include "yarray.h"

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
