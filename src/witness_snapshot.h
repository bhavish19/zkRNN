#pragma once

#include "types.h"
#include "RNN.h"
#include <string>

struct WitnessSnapshot {
  int seq_len = 0;
  int input_size = 0;
  int hidden_size = 0;
  int output_size = 0;

  FieldVector input_t0;
  FieldVector h_prev_t0;
  FieldVector a_t0;
  FieldVector h_t0;
  FieldVector z_t0;
  FieldVector y_hat_t0;

  ExpBatch exp_tanh;
  ExpBatch exp_softmax;

  std::vector<FieldVector> W_x;
  std::vector<FieldVector> W_h;
  std::vector<FieldVector> W_y;
  FieldVector b1;
  FieldVector b2;
};

WitnessSnapshot LoadWitnessSnapshot(const std::string &path);


