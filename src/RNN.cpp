#include "quantization.h"
#include "utils.hpp"
#include "RNN.h"

const F SCALE = F(1ULL<<Q);

static const F INV_SCALE = F(1)/SCALE;

const F ONE = SCALE;

static vector<F> init_row(int dim)
{
    vector<F> row(dim);
    for (int i = 0; i < dim; ++i)
    {
        row[i] = quantize((float)rand() / (float)RAND_MAX);
    }
    return row;
}

static vector<vector<F>> init_input(int T, int n)
{
    vector<vector<F>> in(T);
    for (int t = 0; t < T; ++t)
    {
        in[t] = init_row(n);
    }
    return in;
}

softmax_layer _softmax_layer(vector<F> z){
    softmax_layer sl;
    int n = (int)z.size();
    if (n == 0) return sl;
    sl.M_in = z;
    sl.M_out.assign(n, F_ZERO);
    sl.z_shift.assign(n, F_ZERO);
    sl.exp_shift.assign(n, F_ZERO);
    sl.remainder.assign(n, F_ZERO);

    F z_max = z[0];
    for(int i = 1; i <n; i++){
        if (!(z[i] - z_max).isNegative()) z_max = z[i];
    }
    sl.z_max = z_max;


    for (int i = 0; i < n; i++){
        sl.z_shift[i] = z[i] -sl.z_max;
        sl.exp_shift[i] = exp(sl.z_shift[i]);
    }

    sl.den = F_ZERO;
    for (int i = 0; i < n; i++){
        sl.den +=sl.exp_shift[i];
    }

    for (int i = 0; i < n; i++){
        F y = divide(sl.exp_shift[i], sl.den);
        sl.M_out[i] = y;
        sl.remainder[i] = sl.exp_shift[i] * SCALE - y * sl.den;
    }


    return sl;
}

tanh_layer _tanh_layer(vector<F> input){
    tanh_layer tl;
    tl.input = input;
    tl.output.resize(input.size(), F_ZERO);
    tl.two_x.resize(input.size(), F_ZERO);
    tl.e2x.resize(input.size(),   F_ZERO);
    tl.num.resize(input.size(),   F_ZERO);
    tl.den.resize(input.size(),   F_ZERO);
    int R = (int)input.size();
    for (int i = 0; i<R; i++){
        tl.two_x[i] = F(2)*input[i];
        tl.e2x[i] = exp(tl.two_x[i]);
        tl.num[i] = tl.e2x[i]-ONE;
        tl.den[i] = tl.e2x[i]+ONE;
        tl.output[i] = tl.num[i]* SCALE/tl.den[i];
    }
    return tl;
}

tanh_layer_backprop tanh_backprop(vector<F> dx, tanh_layer tl){
    tanh_layer_backprop der;
    int R = (int)dx.size();
    der.dx_prev.resize(R);
    der.dx_prev = dx;

    der.dx.resize(R, F_ZERO);
    der.one_minus_y2.resize(R);
    der.y2.resize(R);
    der.one_minus_y2.resize(R); 

    for(int i = 0; i < R; i++ ){
        F y = tl.output[i];
        F yy = y*y*INV_SCALE;
        F dy = ONE - yy;
        if (dy < F_ZERO) dy = F_ZERO;
        der.y2[i] = yy;
        der.one_minus_y2[i] = dy;

        der.dx[i] = dx[i]*dy*INV_SCALE;
    }
    return der;
}



rnn_layer rnn_init(
    int seq_len,
    int input_size,
    int hidden_size,
    int output_size)
{
    rnn_layer net;

    net.seq_len = seq_len;
    net.input_size = input_size;
    net.hidden_size = hidden_size;
    net.output_size = output_size;
    printf("RNN init, sequence length : %d \n", net.seq_len);

    net.W_x.resize(hidden_size);
    for (int i = 0; i < hidden_size; ++i)
    {
        net.W_x[i] = init_row(input_size);
    }

    net.W_h.resize(hidden_size);
    for (int i = 0; i < hidden_size; ++i)
    {
        net.W_h[i] = init_row(hidden_size);
    }
    net.W_y.resize(output_size);
    for (int i = 0; i < output_size; ++i)
    {
        net.W_y[i] = init_row(hidden_size);
    }

    net.b1.assign(hidden_size, F(0));
    net.b2.assign(output_size, F(0));

    
    net.fwd.x.clear();
    net.fwd.a.clear();
    net.fwd.h.clear();
    net.fwd.h0.clear();
    net.fwd.z.clear();
    net.fwd.yHat.clear();
    net.fwd.Wh.clear();
    net.fwd.Wx.clear();
    net.fwd.WhWx.clear();
    net.fwd.Zy.clear();
    net.fwd.softmaxes.clear();
    net.fwd.tanhs.clear();
    net.fwd.exp_tanh.clear();
    net.fwd.exp_softmax.clear();

    net.bwd.delta_y.clear();
    net.bwd.dh.clear();
    net.bwd.da.clear();
    net.bwd.WhT_da.clear();
    net.bwd.WyT_delta_y.clear();
    net.bwd.tanh_backprops.clear();

    net.bwd.dW_x.clear();
    net.bwd.dW_h.clear();
    net.bwd.dW_y.clear();
    net.bwd.db1.clear();
    net.bwd.db2.clear();

    return net;
}

rnn_layer rnn_update(rnn_layer net, F lr){
    printf("Begin Rnn Update.\n");
    int T = net.seq_len;
    int n = net.input_size;
    int m = net.hidden_size;
    int k = net.output_size;
    net.learning_rate = lr;

    net.W_h_next.assign(m, vector<F>(m, F(0)));
    net.W_x_next.assign(m, vector<F>(n, F(0)));
    net.W_y_next.assign(k, vector<F>(m, F(0)));
    net.b1_next.assign(m, F(0));
    net.b2_next.assign(k, F(0));

    for (int i=0;i<m;i++){
        for (int j=0;j<n;j++)
            net.W_x_next[i][j] = net.W_x[i][j] - lr * net.bwd.dW_x[i][j] *INV_SCALE;
        for (int j=0;j<m;j++)
            net.W_h_next[i][j] = net.W_h[i][j] - lr * net.bwd.dW_h[i][j] *INV_SCALE;
        net.b1_next[i] = net.b1[i] - lr * net.bwd.db1[i] *INV_SCALE;
    }
    for (int i=0;i<k;i++){
        for (int j=0;j<m;j++)
            net.W_y_next[i][j] = net.W_y[i][j] - lr * net.bwd.dW_y[i][j] *INV_SCALE;
        net.b2_next[i] = net.b2[i] - lr * net.bwd.db2[i] *INV_SCALE;
    }

    return net;
}

rnn_layer rnn_forward(rnn_layer net, vector<vector<F>> x0, const vector<F> h0)
{
    printf("Begin Rnn Forward.\n");
    // x0: [T x n], h0: [m]
    int T = net.seq_len;
    int n = net.input_size;
    int m = net.hidden_size;
    int k = net.output_size;

    // allocate: a[T x m], h[T x m], z[T x k], yHat[T x k]
    net.fwd.x.assign(T, vector<F>(n, F(0)));
    net.fwd.a.assign(T, vector<F>(m, F(0)));
    net.fwd.h.assign(T, vector<F>(m, F(0)));
    net.fwd.h0.assign(m, F(0));
    net.fwd.z.assign(T, vector<F>(k, F(0)));
    net.fwd.yHat.assign(T, vector<F>(k, F(0)));
    net.fwd.Wh.assign(T, vector<F>(m, F(0)));
    net.fwd.Wx.assign(T, vector<F>(m, F(0)));
    net.fwd.WhWx.assign(T, vector<F>(m, F(0)));
    net.fwd.Zy.assign(T, vector<F>(k, F(0)));
    net.fwd.exp_tanh.reserve(static_cast<size_t>(T) * static_cast<size_t>(m));
    net.fwd.exp_softmax.reserve(static_cast<size_t>(T) * static_cast<size_t>(k));
        
    net.fwd.x = x0;
    vector<F> prev_h = h0; // [m]
    net.fwd.h0 = h0;

    for (int t = 0; t < T; ++t)
    {
        // 1) Wh = W_h * prev_h
        auto Wh = _mat_vec_mul(net.W_h, prev_h);
        net.fwd.Wh[t] = Wh;

        // 2) Wx = W_x * x[t]
        auto Wx = _mat_vec_mul(net.W_x, net.fwd.x[t]);
        net.fwd.Wx[t] = Wx;

        // 3) WhWx = Wh + Wx
        for (int i = 0; i < m; ++i)
        {
            net.fwd.WhWx[t][i] = net.fwd.Wh[t][i] + net.fwd.Wx[t][i];
        }
        
        // 4) a = WhWx + b1
        for (int i = 0; i < m; ++i)
        {
            net.fwd.a[t][i] = net.fwd.WhWx[t][i] + net.b1[i];
        }

        // 4) h = tanh(a)
        auto tl = _tanh_layer(net.fwd.a[t]);
        net.fwd.tanhs.push_back(tl);
        net.fwd.h[t] = tl.output;
        for (int i = 0; i < m; ++i) {
            net.fwd.exp_tanh.emit(tl.two_x[i], tl.e2x[i]);  // x=2a_i, y=exp(x)
        }

        // 5) Zy = W_y * h[t]
        auto Zy = _mat_vec_mul(net.W_y, net.fwd.h[t]);
        net.fwd.Zy[t] = Zy;

        // 6) z = Zy + b2
        for (int i = 0; i < k; ++i)
        {
            net.fwd.z[t][i] = Zy[i] + net.b2[i];
        }

        auto soft = _softmax_layer(net.fwd.z[t]);
        net.fwd.softmaxes.push_back(soft);
        for (int i = 0; i < k; ++i) {
            net.fwd.exp_softmax.emit(soft.z_shift[i], soft.exp_shift[i]); // x=z_i-zmax, y=exp(x)
        }

        net.fwd.yHat[t] = soft.M_out;

        prev_h = net.fwd.h[t];
    }
    return net;
}

rnn_layer rnn_bptt(rnn_layer net, const vector<vector<F>> &y_true)
{
    printf("Begin Rnn Backward.\n");

    int T = net.seq_len;
    int n = net.input_size;
    int m = net.hidden_size;
    int k = net.output_size;

    net.bwd.delta_y.assign(T, vector<F>(k, F(0)));
    net.bwd.dh.assign(T, vector<F>(m, F(0)));
    net.bwd.da.assign(T, vector<F>(m, F(0)));
    net.bwd.dW_x.assign(m, vector<F>(n, F(0)));
    net.bwd.dW_h.assign(m, vector<F>(m, F(0)));
    net.bwd.dW_y.assign(k, vector<F>(m, F(0)));
    net.bwd.WyT_delta_y.assign(T, vector<F>(m, F(0)));
    net.bwd.WhT_da.assign(T, vector<F>(m, F(0)));
    net.bwd.db1.assign(m, F(0));
    net.bwd.db2.assign(k, F(0));
    net.bwd.tanh_backprops.clear();
    net.bwd.tanh_backprops.resize(T);


    vector<F> delta_h_next(m, F(0));

    for (int t = T - 1; t >= 0; --t)
    {
        // delta_y
        for (int i = 0; i < k; ++i)
            net.bwd.delta_y[t][i] = net.fwd.yHat[t][i] - y_true[t][i];
        // W_y, b2 梯度
        for (int i = 0; i < k; ++i)
        {
            for (int j = 0; j < m; ++j)
                net.bwd.dW_y[i][j] += net.bwd.delta_y[t][i] * net.fwd.h[t][j]*INV_SCALE;
            net.bwd.db2[i] += net.bwd.delta_y[t][i];
        }
        // dh = W_y^T * delta_y + delta_h_next
        auto dh_cols = _mat_vec_mul(transpose(net.W_y),net.bwd.delta_y[t]);
        assert((int)net.bwd.WyT_delta_y[t].size() == m);
        net.bwd.WyT_delta_y[t] = dh_cols;

        for(int i = 0; i < m; i++){
            net.bwd.dh[t][i] = net.bwd.WyT_delta_y[t][i] + delta_h_next[i];
        }

        auto tlb = tanh_backprop(net.bwd.dh[t], net.fwd.tanhs[t]);
        net.bwd.tanh_backprops[t]= tlb;
        net.bwd.da[t] = tlb.dx;

        // e) W_h, W_x, b1 梯度
        vector<F> prev_h = (t > 0 ? net.fwd.h[t - 1]
                                  : net.fwd.h0);
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < m; ++j)
                net.bwd.dW_h[i][j] += net.bwd.da[t][i] * prev_h[j]*INV_SCALE;
            for (int j = 0; j < n; ++j)
                net.bwd.dW_x[i][j] += net.bwd.da[t][i] * net.fwd.x[t][j]*INV_SCALE;
            net.bwd.db1[i] += net.bwd.da[t][i];
        }

        // f) delta_h_next = W_h^T * da
        delta_h_next = _mat_vec_mul(transpose(net.W_h), net.bwd.da[t]);
        net.bwd.WhT_da[t] = delta_h_next;
    }
    return net;
}