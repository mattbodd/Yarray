#ifndef YOPERATIONS_H
#define YOPERATIONS_H

#include "yarray.h"
#include "yutilities.h"

void update_point(yarr *, int, yarr *, int);
void update_C_array(yarr *, yarr *, int, int *);
void reform_C_array(yarr *, int *, int);
void shrink_C_array(yarr *, int *, int, int, int);
double matrix_mult(yarr *, yarr *);
yarr *broadcast(yarr *, yarr*, Op);
yarr *apply_op(yarr *, yarr*, Op);
void *grab_point(yarr *, int *, int);

#endif /* YOPERATIONS_H */
