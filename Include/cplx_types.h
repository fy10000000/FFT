#pragma once

#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>


typedef struct { float r, i; } c32;

c32 get_conj(const c32 x);

c32 mult(const c32 a, const c32 b);

c32 add(const c32 a, const c32 b);

float mag(const c32 in);

