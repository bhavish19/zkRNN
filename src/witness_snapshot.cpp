#include "witness_snapshot.h"
#include "proof_serialization.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cctype>

namespace {

struct Cursor {
  const std::string &s;
  size_t pos{0};
  explicit Cursor(const std::string &str) : s(str) {}

  void SkipWhitespace() {
    while (pos < s.size() && std::isspace(static_cast<unsigned char>(s[pos])))
      ++pos;
  }

  void Expect(char c) {
    SkipWhitespace();
    if (pos >= s.size() || s[pos] != c) {
      throw std::runtime_error("JSON parse error: expected character missing");
    }
    ++pos;
  }

  bool TryConsume(char c) {
    SkipWhitespace();
    if (pos < s.size() && s[pos] == c) {
      ++pos;
      return true;
    }
    return false;
  }

  std::string ParseString() {
    SkipWhitespace();
    if (pos >= s.size() || s[pos] != '"') {
      throw std::runtime_error("JSON parse error: expected string");
    }
    ++pos;
    std::string out;
    while (pos < s.size()) {
      char c = s[pos++];
      if (c == '"')
        break;
      if (c == '\\') {
        if (pos >= s.size())
          throw std::runtime_error("JSON parse error: bad escape");
        char esc = s[pos++];
        out.push_back(esc);
      } else {
        out.push_back(c);
      }
    }
    return out;
  }

  long long ParseInt() {
    SkipWhitespace();
    size_t start = pos;
    if (pos < s.size() && (s[pos] == '-' || s[pos] == '+'))
      ++pos;
    while (pos < s.size() && std::isdigit(static_cast<unsigned char>(s[pos])))
      ++pos;
    if (start == pos)
      throw std::runtime_error("JSON parse error: expected integer");
    return std::stoll(s.substr(start, pos - start));
  }
};

std::vector<std::string> ParseStringArray(Cursor &cur);
std::vector<std::vector<std::string>> ParseStringMatrix(Cursor &cur);
void ParseObject(Cursor &cur, WitnessSnapshot &snapshot);

FieldVector ConvertStringsToFields(const std::vector<std::string> &raw) {
  FieldVector out;
  out.reserve(raw.size());
  for (const auto &token : raw) {
    out.push_back(StringToFieldElement(token));
  }
  return out;
}

std::vector<FieldVector> ConvertMatrix(const std::vector<std::vector<std::string>> &raw) {
  std::vector<FieldVector> out;
  out.reserve(raw.size());
  for (const auto &row : raw) {
    out.push_back(ConvertStringsToFields(row));
  }
  return out;
}

std::vector<std::string> ParseStringArray(Cursor &cur) {
  std::vector<std::string> arr;
  cur.Expect('[');
  while (true) {
    cur.SkipWhitespace();
    if (cur.TryConsume(']'))
      break;
    arr.push_back(cur.ParseString());
    cur.SkipWhitespace();
    if (cur.TryConsume(']'))
      break;
    cur.Expect(',');
  }
  return arr;
}

std::vector<long long> ParseIntArray(Cursor &cur) {
  std::vector<long long> arr;
  cur.Expect('[');
  while (true) {
    cur.SkipWhitespace();
    if (cur.TryConsume(']'))
      break;
    arr.push_back(cur.ParseInt());
    cur.SkipWhitespace();
    if (cur.TryConsume(']'))
      break;
    cur.Expect(',');
  }
  return arr;
}

std::vector<std::vector<std::string>> ParseStringMatrix(Cursor &cur) {
  std::vector<std::vector<std::string>> matrix;
  cur.Expect('[');
  while (true) {
    cur.SkipWhitespace();
    if (cur.TryConsume(']'))
      break;
    auto row = ParseStringArray(cur);
    matrix.push_back(std::move(row));
    cur.SkipWhitespace();
    if (cur.TryConsume(']'))
      break;
    cur.Expect(',');
  }
  return matrix;
}

void ParseExpBatch(Cursor &cur, ExpBatch &batch) {
  cur.Expect('{');
  while (true) {
    cur.SkipWhitespace();
    if (cur.TryConsume('}'))
      break;
    std::string key = cur.ParseString();
    cur.Expect(':');
    if (key == "Xq") {
      auto raw = ParseStringArray(cur);
      batch.Xq = ConvertStringsToFields(raw);
    } else if (key == "Yq") {
      auto raw = ParseStringArray(cur);
      batch.Yq = ConvertStringsToFields(raw);
    } else {
      // Skip unhandled field
      int brace_depth = 0;
      cur.SkipWhitespace();
      if (cur.TryConsume('{')) {
        ++brace_depth;
        while (brace_depth > 0 && cur.pos < cur.s.size()) {
          if (cur.s[cur.pos] == '{') ++brace_depth;
          else if (cur.s[cur.pos] == '}') --brace_depth;
          ++cur.pos;
        }
      } else if (cur.TryConsume('[')) {
        int depth = 1;
        while (depth > 0 && cur.pos < cur.s.size()) {
          if (cur.s[cur.pos] == '[') ++depth;
          else if (cur.s[cur.pos] == ']') --depth;
          ++cur.pos;
        }
      } else {
        cur.ParseString();
      }
    }
    cur.SkipWhitespace();
    if (cur.TryConsume('}'))
      break;
    cur.Expect(',');
  }
}

void ParseWeights(Cursor &cur, WitnessSnapshot &snapshot) {
  cur.Expect('{');
  while (true) {
    cur.SkipWhitespace();
    if (cur.TryConsume('}'))
      break;
    std::string key = cur.ParseString();
    cur.Expect(':');
    if (key == "W_x") {
      auto mat = ParseStringMatrix(cur);
      snapshot.W_x = ConvertMatrix(mat);
    } else if (key == "W_h") {
      auto mat = ParseStringMatrix(cur);
      snapshot.W_h = ConvertMatrix(mat);
    } else if (key == "W_y") {
      auto mat = ParseStringMatrix(cur);
      snapshot.W_y = ConvertMatrix(mat);
    } else if (key == "b1") {
      auto raw = ParseStringArray(cur);
      snapshot.b1 = ConvertStringsToFields(raw);
    } else if (key == "b2") {
      auto raw = ParseStringArray(cur);
      snapshot.b2 = ConvertStringsToFields(raw);
    } else {
      // Skip unknown field
      int brace_depth = 0;
      cur.SkipWhitespace();
      if (cur.TryConsume('{')) {
        ++brace_depth;
        while (brace_depth > 0 && cur.pos < cur.s.size()) {
          if (cur.s[cur.pos] == '{') ++brace_depth;
          else if (cur.s[cur.pos] == '}') --brace_depth;
          ++cur.pos;
        }
      } else {
        (void)ParseStringArray(cur);
      }
    }
    cur.SkipWhitespace();
    if (cur.TryConsume('}'))
      break;
    cur.Expect(',');
  }
}

void ParseObject(Cursor &cur, WitnessSnapshot &snapshot) {
  cur.Expect('{');
  while (true) {
    cur.SkipWhitespace();
    if (cur.TryConsume('}'))
      break;
    std::string key = cur.ParseString();
    cur.Expect(':');
    if (key == "seq_len") {
      snapshot.seq_len = static_cast<int>(cur.ParseInt());
    } else if (key == "input_size") {
      snapshot.input_size = static_cast<int>(cur.ParseInt());
    } else if (key == "hidden_size") {
      snapshot.hidden_size = static_cast<int>(cur.ParseInt());
    } else if (key == "output_size") {
      snapshot.output_size = static_cast<int>(cur.ParseInt());
    } else if (key == "input_t0") {
      snapshot.input_t0 = ConvertStringsToFields(ParseStringArray(cur));
    } else if (key == "h_prev_t0") {
      snapshot.h_prev_t0 = ConvertStringsToFields(ParseStringArray(cur));
    } else if (key == "a_t0") {
      snapshot.a_t0 = ConvertStringsToFields(ParseStringArray(cur));
    } else if (key == "h_t0") {
      snapshot.h_t0 = ConvertStringsToFields(ParseStringArray(cur));
    } else if (key == "z_t0") {
      snapshot.z_t0 = ConvertStringsToFields(ParseStringArray(cur));
    } else if (key == "y_hat_t0") {
      snapshot.y_hat_t0 = ConvertStringsToFields(ParseStringArray(cur));
    } else if (key == "exp_tanh") {
      ParseExpBatch(cur, snapshot.exp_tanh);
    } else if (key == "exp_softmax") {
      ParseExpBatch(cur, snapshot.exp_softmax);
    } else if (key == "weights") {
      ParseWeights(cur, snapshot);
    } else {
      // Skip unknown entries
      cur.SkipWhitespace();
      if (cur.TryConsume('{')) {
        int depth = 1;
        while (depth > 0 && cur.pos < cur.s.size()) {
          if (cur.s[cur.pos] == '{') ++depth;
          else if (cur.s[cur.pos] == '}') --depth;
          ++cur.pos;
        }
      } else if (cur.TryConsume('[')) {
        int depth = 1;
        while (depth > 0 && cur.pos < cur.s.size()) {
          if (cur.s[cur.pos] == '[') ++depth;
          else if (cur.s[cur.pos] == ']') --depth;
          ++cur.pos;
        }
      } else {
        cur.ParseString();
      }
    }
    cur.SkipWhitespace();
    if (cur.TryConsume('}'))
      break;
    cur.Expect(',');
  }
}

} // namespace

WitnessSnapshot LoadWitnessSnapshot(const std::string &path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("Cannot open snapshot file: " + path);
  }
  std::string json((std::istreambuf_iterator<char>(in)),
                   std::istreambuf_iterator<char>());
  Cursor cur(json);
  WitnessSnapshot snapshot;
  ParseObject(cur, snapshot);
  return snapshot;
}


