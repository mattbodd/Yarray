#ifndef YUTILITIES_H
#define YUTILITIES_H

#include "yarray.h"

dataType widest_dt(dataType, dataType);
double get_double(yarr *, int);
long get_long(yarr *, int);
float get_float(yarr *, int);
int get_int(yarr *, int);
void move_elem(yarr *, int, int);
void print_C_array(yarr *);
void print_yarr(yarr *);
#endif /* YUTILITIES_H */
