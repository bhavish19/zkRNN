#include "config_pc.hpp"

struct ExpBatch{
    vector<F> Xq;
    vector<F> Yq;

    void clear(){
        Xq.clear();
        Yq.clear();
    }

    void reserve(size_t n){
        Xq.reserve(n);
        Yq.reserve(n);
    }
    
    void emit(const F& x, const F& y){
        Xq.push_back(x);
        Yq.push_back(y);
    }
}; 

struct softmax_layer{
    // 输入矩阵 M（同你 _softmax 里传入的 M）
    vector<F> M_in;

    // 输出 softmax 后的值 M_softmax
    vector<F> M_out;

    vector<F> z_shift;
    vector<F> exp_shift;
    vector<F> remainder;

    F z_max;
    F den;
};

struct tanh_layer{
    vector<F> input;
    vector<F> output;
    vector<F> two_x;
    vector<F> e2x;
    vector<F> num, den;
};

struct tanh_layer_backprop{
    vector<F> dx;
    vector<F> dx_prev;
    vector<F> y2;
    vector<F> one_minus_y2;
};


struct rnn_fwd {
    // 原始输入与激活
    vector<vector<F>> x;      // [T][n]
    vector<vector<F>> a;      // [T][m]
    vector<vector<F>> h;      // [T][m]
    vector<vector<F>> z;      // [T][k]
    vector<vector<F>> yHat;   // [T][k]
    vector<F> h0;

    // 中间矩阵乘积
    vector<vector<F>> Wh;     // [T][m] = W_h * prev_h
    vector<vector<F>> Wx;     // [T][m] = W_x * x[t]
    vector<vector<F>> Zy;     // [T][k] = W_y * h[t]
    vector<vector<F>> WhWx;   // [T][m] = Wh + Wx

    vector<softmax_layer> softmaxes;
    vector<tanh_layer> tanhs;

    ExpBatch exp_tanh;
    ExpBatch exp_softmax;
};

//—— 反向缓存结构 ——  
struct rnn_bwd {
    // 输出误差与隐藏误差
    vector<vector<F>> delta_y;       // [T][k]
    vector<vector<F>> dh;            // [T][m]
    vector<vector<F>> da;            // [T][m]

    vector<vector<F>> WyT_delta_y;  // [T][m]
    vector<vector<F>> WhT_da;


    vector<tanh_layer_backprop> tanh_backprops;


    // 累积梯度
    vector<vector<F>> dW_x;  // [m][n]
    vector<vector<F>> dW_h;  // [m][m]
    vector<vector<F>> dW_y;  // [k][m]
    vector<F>         db1;   // [m]
    vector<F>         db2;   // [k]
};

//—— 完整 RNN 层定义 ——  
struct rnn_layer {
    //—— 超参数 & 可训练参数 ——  
    int seq_len;       // T
    int input_size;    // n
    int hidden_size;   // m
    int output_size;   // k

    vector<vector<F>> W_x;  // [m][n]
    vector<vector<F>> W_h;  // [m][m]
    vector<vector<F>> W_y;  // [k][m]
    vector<F>         b1;   // [m]
    vector<F>         b2;   // [k]

    vector<vector<F>> W_x_next;  // [m][n]
    vector<vector<F>> W_h_next;  // [m][m]
    vector<vector<F>> W_y_next;  // [k][m]
    vector<F>         b1_next;   // [m]
    vector<F>         b2_next;   // [k]

    F learning_rate;

    //—— 前向与反向缓存 ——  
    rnn_fwd fwd;
    rnn_bwd bwd;
};



// Initialize an RNN layer with random weights & zero biases
rnn_layer rnn_init(
    int seq_len,
    int input_size,
    int hidden_size,
    int output_size
);

/**
 * 单层前向计算: x^(t) -> a^(t) -> h^(t) -> z^(t) -> yHat^(t)
 * @param net 已初始化参数的 rnn_layer
 * @param x0  输入序列, 大小 seq_len x input_size
 * @param h0  初始隐藏状态, 大小 hidden_size
 * @return yHat 预测序列, 大小 seq_len x output_size
 */
rnn_layer rnn_forward(
    rnn_layer net,
    vector<vector<F>> x0,
    const vector<F> h0
);

/**
 * 单层 BPTT 反向传播: 计算梯度并更新参数
 * @param net    包含前向结果的 rnn_layer
 * @param y_true 真实标签, 大小 seq_len x output_size
 */
rnn_layer rnn_bptt(rnn_layer net,const vector<vector<F>> &y_true);

rnn_layer rnn_update(rnn_layer net, F lr);

softmax_layer _softmax_layer(vector<F> M);

tanh_layer _tanh_layer(vector<F> input);

tanh_layer_backprop tanh_backprop(vector<F> dx, tanh_layer tl);


