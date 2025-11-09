#include "activation_circuit_loader.h"

#include <cstdlib>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif

namespace {

constexpr const char *kCircuitDirName = "circuits";

// Supported circuit fan-outs extracted from upstream KAIZEN artefacts.
constexpr std::size_t kTanhDims[] = {2, 6, 12};
constexpr std::size_t kSoftmaxDims[] = {6, 12};

std::string ActivationName(ActivationType type) {
  switch (type) {
  case ACTIVATION_TANH:
    return "tanh";
  case ACTIVATION_SOFTMAX:
    return "softmax";
  default:
    return "unknown";
  }
}

std::string BuildCircuitFileName(ActivationType type, std::size_t dim) {
  std::ostringstream oss;
  oss << "activation_" << ActivationName(type) << "_d" << dim << ".pws";
  return oss.str();
}

void CopyFile(const std::string &src, const std::string &dst) {
  std::ifstream in(src, std::ios::binary);
  if (!in) {
    throw std::runtime_error("Unable to open source circuit asset: " + src);
  }
  std::ofstream out(dst, std::ios::binary);
  if (!out) {
    throw std::runtime_error("Unable to create circuit file: " + dst);
  }
  out << in.rdbuf();
  if (!out.good()) {
    throw std::runtime_error("Failed to copy circuit artefact from " + src + " to " + dst);
  }
}

void EnsureDirectoryExists(const std::string &dir) {
#ifdef _WIN32
  _mkdir(dir.c_str());
#else
  mkdir(dir.c_str(), 0755);
#endif
}

bool FileExists(const std::string &path) {
  std::ifstream in(path);
  return in.good();
}

std::vector<std::string> CandidateCircuitRoots() {
  std::vector<std::string> roots;
  if (const char *env = std::getenv("KAIZEN_CIRCUIT_DIR")) {
    roots.emplace_back(env);
  }
#ifdef PROJECT_SOURCE_DIR
  roots.emplace_back(std::string(PROJECT_SOURCE_DIR) + "/" + kCircuitDirName);
#endif
  roots.emplace_back(kCircuitDirName);
  roots.emplace_back(std::string("..") + "/" + kCircuitDirName);
  roots.emplace_back(std::string("..") + "/.." + "/" + kCircuitDirName);
  return roots;
}

std::string LocateCircuitAsset(const std::string &file_name) {
  for (const auto &root : CandidateCircuitRoots()) {
    if (root.empty()) continue;
    std::string candidate = root + "/" + file_name;
    if (FileExists(candidate)) {
      return candidate;
    }
  }
  return {};
}

std::size_t SelectCircuitDimension(ActivationType type, std::size_t requested) {
  const std::size_t *begin = nullptr;
  const std::size_t *end = nullptr;
  switch (type) {
  case ACTIVATION_TANH:
    begin = std::begin(kTanhDims);
    end = std::end(kTanhDims);
    break;
  case ACTIVATION_SOFTMAX:
    begin = std::begin(kSoftmaxDims);
    end = std::end(kSoftmaxDims);
    break;
  default:
    break;
  }

  if (begin == nullptr || begin == end) {
    throw std::runtime_error("No activation circuits registered for type: " + ActivationName(type));
  }

  for (const std::size_t *it = begin; it != end; ++it) {
    if (*it >= requested) {
      return *it;
    }
  }

  throw std::runtime_error(
      "Requested " + ActivationName(type) + " circuit of dimension " + std::to_string(requested) +
      " exceeds available KAIZEN assets");
}

} // namespace

ActivationCircuitInfo EnsureActivationCircuit(ActivationType type, std::size_t dimension) {
  if (dimension == 0) {
    throw std::invalid_argument("Activation circuit dimension must be positive");
  }

  std::string circuit_dir = kCircuitDirName;
  EnsureDirectoryExists(circuit_dir);

  std::size_t selected_dim = SelectCircuitDimension(type, dimension);
  std::string circuit_path = circuit_dir + "/" + BuildCircuitFileName(type, selected_dim);
  if (!FileExists(circuit_path)) {
    const std::string file_name = BuildCircuitFileName(type, selected_dim);
    std::string asset_source = LocateCircuitAsset(file_name);
    if (asset_source.empty()) {
      throw std::runtime_error("Missing activation circuit file: " + circuit_path +
                               ". Please import the KAIZEN artefact.");
    }
    CopyFile(asset_source, circuit_path);
  }

  ActivationCircuitInfo info;
  info.circuit_path = circuit_path;
  info.input_label = ActivationName(type);
  info.lanes = selected_dim;
  return info;
}


