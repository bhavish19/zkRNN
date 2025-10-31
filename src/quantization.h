#pragma once
#include "config_pc.hpp"

#define Q 32
#define M_exp 8
#ifndef QUANT_Q_BITS
#define QUANT_Q_BITS 30  // q，比特位宽（含符号），按需改
#endif

F quantize(float num);
F exp(F num);
vector<F> shift(F num1, int depth);
F divide(F num1,F num2);
float dequantize(F num,int depth);
