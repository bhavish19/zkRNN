#include <math.h>
#include "quantization.h"
#include "config_pc.hpp"



F quantize(float num){
	long long int s = (long long int)(num * (1ULL << Q));
#if QUANT_Q_BITS > 0
    const long long lo = -(1LL << (QUANT_Q_BITS - 1));
    const long long hi =  (1LL << (QUANT_Q_BITS - 1)) - 1LL;
    if (s < lo) s = lo;
    if (s > hi) s = hi;
#endif
    return F(s);
}

float dequantize(F num,int depth){
    unsigned long long int div = 1;
    long long int r_num;

    for(int i = 0; i < depth-1; i++){
        div = div*(1ULL<<Q);
    }

    if(num.getBit(60)){
        num = F(0) - num;
        r_num = num.toint128();
        r_num = r_num/div;
        r_num = -r_num;
    }
    else{
        r_num = num.toint128();
        r_num = r_num/div;
    }
    return (float)(r_num)/(float)(1ULL<<Q);
}

vector<F> shift(F num1, int depth){
    char buff[257];
    long long int dividend,divisor,quotient,remainder;
    vector<F> r;

    F num2 = F(1);
    for(int i = 0; i < depth; i++){
        num2 = num2*F(1ULL<<Q);
    }

    if(num1.getBit(60)){
        dividend = (F(0)-num1).toint128();
        divisor = num2.toint128();
        quotient = dividend/divisor;

        // long int tmp = mpz_get_ui(quotient);
        // r.push_back(F(tmp));
        r.push_back(F(0) - F(quotient) - F(1));
        r.push_back(F(divisor));
        r.push_back(divisor - (dividend - quotient*divisor));
        return r;
    }
    else{
        dividend = num1.toint128();
        divisor = num2.toint128();
        quotient = dividend/divisor;

        r.push_back(F(quotient));
        r.push_back(F(divisor));
        r.push_back(dividend - divisor*quotient);
        return r;
    }
}

F divide(F num1,F num2){
    char buff[257];

    long long int dividend,divisor,quotient;

    if(num1.getBit(60) == num2.getBit(60)){

        if(num1.getBit(60) == 1){
            num1 = F(-1)*num1;
            num2 = F(-1)*num2;
        }
        dividend = num1.toint128();
        divisor = num2.toint128();
        quotient = dividend/divisor;

        return F(quotient);
        // return F(tmp);
    }
    else{
        if(num1.getBit(60) == 1){
            num1 = F(-1)*num1;
        }
        else{
            num2 = F(-1)*num2;
        }

        dividend = num1.toint128();
        divisor = num2.toint128();
        quotient = dividend/divisor;

        // long int tmp = mpz_get_ui(quotient);
        return F(0)-F(quotient);
        // return F(-tmp);
    }
    // return F(0);
}

F exp(F x) {
    const F SCALE = F(1ULL << Q);
    const int M = M_exp;               // 建议 M_exp 是 2 的幂

    F step = divide(x, F(M));          // 仍是 Q 定点：表示 (x/M)
    F base = SCALE + step;             // 表示 (1 + x/M) * SCALE

    int rounds = (int)log2(M);
    for (int i = 0; i < rounds; ++i) {
        base = divide(base * base, SCALE);   // ((base/SCALE)^2) * SCALE
    }

    return base;                       // Q 定点输出，表示 e^x
}