#include "stream_orchestrator.h"

#include "RNN.h"
#include "quantization.h"

#include <algorithm>
#include <iostream>
#include <stdexcept>

namespace {

F QuantizedWeight() {
  return quantize(0.01f);
}

} // namespace

void StreamOrchestrator::EnsureWeightsInitialized(int input_dim,
                                                  int hidden_dim,
                                                  int output_dim) {
  const bool dimensions_match = weights_initialized_ &&
                                input_dim == input_dim_ &&
                                hidden_dim == hidden_dim_ &&
                                output_dim == output_dim_;
  if (dimensions_match) {
    return;
  }

  input_dim_ = input_dim;
  hidden_dim_ = hidden_dim;
  output_dim_ = output_dim;

  weights_.W_x.assign(hidden_dim_, FieldVector(input_dim_, F_ZERO));
  weights_.W_h.assign(hidden_dim_, FieldVector(hidden_dim_, F_ZERO));
  weights_.W_y.assign(output_dim_, FieldVector(hidden_dim_, F_ZERO));
  weights_.b1.assign(hidden_dim_, F_ZERO);
  weights_.b2.assign(output_dim_, F_ZERO);

  // Populate weights with small random quantised values, similar to KAIZEN's init_row.
  auto random_weight = []() {
    const float val = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
    const float scaled = 0.002f * val - 0.001f;
    float clamped = std::max(-0.004f, std::min(0.004f, scaled));
    return quantize(clamped);
  };

  const int min_x = std::min(hidden_dim_, input_dim_);
  for (int i = 0; i < hidden_dim_; ++i) {
    for (int j = 0; j < input_dim_; ++j) {
      weights_.W_x[i][j] = random_weight();
    }
  }
  for (int i = 0; i < hidden_dim_; ++i) {
    for (int j = 0; j < hidden_dim_; ++j) {
      weights_.W_h[i][j] = random_weight();
    }
  }
  for (int i = 0; i < output_dim_; ++i) {
    for (int j = 0; j < hidden_dim_; ++j) {
      weights_.W_y[i][j] = random_weight();
    }
  }

  weights_initialized_ = true;
}

StreamOrchestrator::StreamOrchestrator(const ProveParams &params, int arity_k)
  : params_(params), arity_k_(arity_k), acc_("ACC_INIT"), acc_proof_(), acc_valid_(false) {}

void StreamOrchestrator::OnInput(const TimeStepInput &in) {
  const int requested_hidden = static_cast<int>(in.h_prev.size());
  const int requested_input = static_cast<int>(in.x_t.size());
  const int requested_output = requested_hidden; // default: single-layer RNN with equal hidden/output

  EnsureWeightsInitialized(requested_input, requested_hidden, requested_output);

  const int m = hidden_dim_;
  const int n = input_dim_;
  const int k = output_dim_;

  // Basic dimensionality checks to guard against inconsistent caller input.
  if (static_cast<int>(in.h_prev.size()) != m) {
    throw std::invalid_argument("StreamOrchestrator::OnInput: unexpected hidden dimension");
  }
  if (static_cast<int>(in.x_t.size()) != n) {
    throw std::invalid_argument("StreamOrchestrator::OnInput: unexpected input dimension");
  }

  // Build a KAIZEN rnn_layer for a single timestep and run the quantised forward pass.
  rnn_layer layer;
  layer.seq_len = 1;
  layer.input_size = n;
  layer.hidden_size = m;
  layer.output_size = k;
  layer.W_x = weights_.W_x;
  layer.W_h = weights_.W_h;
  layer.W_y = weights_.W_y;
  layer.b1 = weights_.b1;
  layer.b2 = weights_.b2;

  auto requantize = [](const FieldVector &vec) {
    FieldVector out(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
      float f = dequantize(vec[i], 1);
      const float scale = static_cast<float>(1 << 4);
      int rounded = static_cast<int>(std::round(f * scale));
      f = static_cast<float>(rounded) / scale;
      out[i] = quantize(f);
    }
    return out;
  };

  FieldVector rq_input = requantize(in.x_t);
  FieldVector rq_prev = requantize(in.h_prev);
  std::cerr << "[StreamOrchestrator] rq_input=";
  for (const auto &v : rq_input) {
    std::cerr << " " << dequantize(v, 1);
  }
  std::cerr << " rq_prev=";
  for (const auto &v : rq_prev) {
    std::cerr << " " << dequantize(v, 1);
  }
  std::cerr << "\n";

  rnn_layer forward = rnn_forward(layer, {rq_input}, rq_prev);

  TimeStepOutput out;
  out.a_t = requantize(forward.fwd.a[0]);
  out.h_t = requantize(forward.fwd.h[0]);
  out.z_t = requantize(forward.fwd.z[0]);
  out.yHat_t = requantize(forward.fwd.yHat[0]);
  auto slice_batch = [](const ExpBatch &src, int keep) {
    ExpBatch dst;
    if (keep <= 0) {
      return dst;
    }
    const int available = static_cast<int>(src.Xq.size());
    const int start = std::max(0, available - keep);
    dst.reserve(static_cast<size_t>(available - start));
    for (int i = start; i < available; ++i) {
      dst.emit(src.Xq[i], src.Yq[i]);
    }
    return dst;
  };
  out.exp_tanh = slice_batch(forward.fwd.exp_tanh, m);
  out.exp_softmax = slice_batch(forward.fwd.exp_softmax, k);

  for (auto &val : out.a_t) {
    float f = dequantize(val, 1);
    f = std::max(std::min(f, 1.0f), -1.0f);
    val = quantize(f);
  }

  std::cerr << "[StreamOrchestrator] after forward a_t=";
  for (const auto &val : out.a_t) {
    std::cerr << " " << dequantize(val, 1);
  }
  std::cerr << " h_t=";
  for (const auto &val : out.h_t) {
    std::cerr << " " << dequantize(val, 1);
  }
  std::cerr << "\n";

  TimeStepInput prove_in;
  prove_in.x_t = rq_input;
  prove_in.h_prev = rq_prev;
  std::cerr << "[StreamOrchestrator] a_t=";
  for (const auto &val : out.a_t) {
    std::cerr << dequantize(val, 1) << " ";
  }
  std::cerr << " h_prev=";
  for (const auto &val : in.h_prev) {
    std::cerr << dequantize(val, 1) << " ";
  }
  std::cerr << "\n";

  std::cout << "[StreamOrchestrator::OnInput] Calling ProveLeaf (m=" << m
            << ", n=" << n << ", k=" << k << ")\n";
  LeafResult leaf = ProveLeaf(params_, prove_in, weights_, out);
  std::cout << "[StreamOrchestrator::OnInput] ProveLeaf completed\n";
  buffer_.push_back(leaf);

  // If we've accumulated enough leaves, aggregate them
  if ((int)buffer_.size() >= arity_k_) {
    AggregationResult agg = FA_Aggregate(buffer_, acc_);
    if (agg.ok) {
      acc_ = agg.serialized;
      acc_proof_ = agg.proof_struct;
      acc_valid_ = true;
    } else {
      acc_ = agg.serialized;
    }
    buffer_.clear();
    std::cout << "[StreamOrchestrator] Aggregated chunk -> " << acc_ << "\n";
  }
}

std::string StreamOrchestrator::FinaliseAndOpen() {
    std::vector<std::string> transcripts;
    transcripts.push_back(acc_);
    
    FinalProof proof = FinaliseProof(acc_, transcripts, acc_valid_ ? &acc_proof_ : nullptr);
    bool ok = VerifyProof(proof);
    
    std::cout << "[StreamOrchestrator] Verification result: "
              << (ok ? "SUCCESS" : "FAIL") << "\n";
    
    return proof.root_acc;
}
