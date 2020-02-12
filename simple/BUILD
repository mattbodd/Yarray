cc_library(
  name="object",
  hdrs=["yarray.h"],
  srcs=["yarray.c"],
)

cc_library(
  name="utilities",
  hdrs=["yutilities.h"],
  srcs=["yutilities.c"],
  deps=[":object"],
)

cc_binary(
  name="simple-yarray",
  srcs=["yoperations.c"],
  deps=[":utilities",":object"],
)
