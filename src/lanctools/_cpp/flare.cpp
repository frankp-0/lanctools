#include "flare.hpp"

#include <cstdint>
#include <iostream>
#include <pybind11/stl.h>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>
#include <zlib.h>

namespace py = pybind11;

struct AncestryTract {
  std::string chrom;
  uint32_t spos;
  uint32_t epos;
  uint8_t anc0;
  uint8_t anc1;
};

void split_views(std::string_view s, char delim,
                 std::vector<std::string_view> &fields) {
  fields.clear();

  size_t start = 0;

  while (start <= s.size()) {
    size_t end = s.find(delim, start);

    if (end == std::string_view::npos) {
      fields.emplace_back(s.substr(start));
      break;
    }

    fields.emplace_back(s.substr(start, end - start));
    start = end + 1;
  }
}

uint8_t parse_ancestry_value(std::string_view value) {
  if (value.empty() || value == ".")
    return 255;

  unsigned int result = 0;

  for (char c : value) {
    if (c < '0' || c > '9')
      return 255;

    result = result * 10 + static_cast<unsigned int>(c - '0');

    if (result > 255)
      return 255;
  }

  return static_cast<uint8_t>(result);
}

uint8_t extract_format_value(std::string_view sample_field, int target_idx) {
  if (target_idx < 0)
    return 255;

  int field_idx = 0;
  size_t start = 0;

  while (start <= sample_field.size()) {
    size_t end = sample_field.find(':', start);

    if (end == std::string_view::npos)
      end = sample_field.size();

    if (field_idx == target_idx) {
      return parse_ancestry_value(sample_field.substr(start, end - start));
    }

    if (end == sample_field.size())
      break;

    ++field_idx;
    start = end + 1;
  }

  return 255;
}

int find_format_index(std::string_view format, std::string_view target) {
  int field_idx = 0;
  size_t start = 0;

  while (start <= format.size()) {
    size_t end = format.find(':', start);

    if (end == std::string_view::npos)
      end = format.size();

    if (format.substr(start, end - start) == target)
      return field_idx;

    if (end == format.size())
      break;

    ++field_idx;
    start = end + 1;
  }

  return -1;
}

std::string gz_readline(gzFile file) {
  const size_t chunk_size = 1024 * 1024;
  char buffer[chunk_size];

  std::string line;

  while (gzgets(file, buffer, chunk_size)) {
    line += buffer;

    if (!line.empty() && line.back() == '\n')
      break;
  }

  while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) {
    line.pop_back();
  }

  return line;
}

void finalize_open_tracts(
    const std::vector<std::string> &sample_ids,
    const std::vector<uint8_t> &prev_anc,
    const std::vector<uint32_t> &prev_spos,
    std::vector<std::vector<AncestryTract>> &sample_tracts,
    const std::string &chrom, uint32_t final_pos) {

  for (size_t i = 0; i < sample_ids.size(); ++i) {
    uint8_t anc0 = prev_anc[i * 2];
    uint8_t anc1 = prev_anc[i * 2 + 1];

    if (anc0 != 255 || anc1 != 255) {
      sample_tracts[i].push_back(
          {chrom, prev_spos[i * 2], final_pos, anc0, anc1});
    }
  }
}

py::dict read_flare(const std::string &flare_file) {
  gzFile file = gzopen(flare_file.c_str(), "rb");

  if (!file)
    throw std::runtime_error("Failed to open input VCF file");

  std::string line;
  bool found_header = false;

  while (!(line = gz_readline(file)).empty()) {
    if (line.size() >= 6 && line.compare(0, 6, "#CHROM") == 0) {
      found_header = true;
      break;
    }
  }

  if (!found_header) {
    gzclose(file);
    throw std::runtime_error("Missing #CHROM header line");
  }

  std::vector<std::string_view> header_fields;
  header_fields.reserve(64);

  split_views(line, '\t', header_fields);

  int chrom_idx = -1;
  int pos_idx = -1;
  int format_idx = -1;

  for (size_t i = 0; i < header_fields.size(); ++i) {
    if (header_fields[i] == "#CHROM")
      chrom_idx = static_cast<int>(i);
    else if (header_fields[i] == "POS")
      pos_idx = static_cast<int>(i);
    else if (header_fields[i] == "FORMAT")
      format_idx = static_cast<int>(i);
  }

  if (chrom_idx == -1 || pos_idx == -1 || format_idx == -1) {
    gzclose(file);
    throw std::runtime_error("Missing essential VCF columns");
  }

  std::vector<std::string> sample_ids;

  if (format_idx + 1 < static_cast<int>(header_fields.size())) {
    sample_ids.reserve(header_fields.size() - format_idx - 1);

    for (size_t i = format_idx + 1; i < header_fields.size(); ++i) {
      sample_ids.emplace_back(header_fields[i]);
    }
  }

  const size_t n_samples = sample_ids.size();

  std::vector<std::vector<AncestryTract>> sample_tracts(n_samples);

  std::vector<uint8_t> prev_anc(n_samples * 2, 255);

  std::vector<uint32_t> prev_spos(n_samples * 2, 0);

  uint32_t prev_pos = 0;
  std::string cur_chrom = "chr0";

  bool is_first_record = true;

  size_t variant_count = 0;

  std::vector<std::string_view> fields;
  fields.reserve(format_idx + n_samples + 8);

  while (!(line = gz_readline(file)).empty()) {
    if (line.empty() || line[0] == '#')
      continue;

    split_views(line, '\t', fields);

    if (fields.size() < static_cast<size_t>(format_idx + 1 + n_samples)) {
      continue;
    }

    ++variant_count;

    if (variant_count % 10000 == 0) {
      std::cerr << "Processed " << variant_count << " variants\n";
    }

    std::string_view chrom_view = fields[chrom_idx];
    std::string_view pos_view = fields[pos_idx];
    std::string_view format_view = fields[format_idx];

    uint32_t pos = 0;

    for (char c : pos_view) {
      if (c < '0' || c > '9')
        break;

      pos = pos * 10 + static_cast<uint32_t>(c - '0');
    }

    int an1_idx = find_format_index(format_view, "AN1");
    int an2_idx = find_format_index(format_view, "AN2");

    if (is_first_record) {
      cur_chrom.assign(chrom_view);

      for (size_t i = 0; i < n_samples; ++i) {
        std::string_view sample_field = fields[format_idx + 1 + i];

        uint8_t anc0 = extract_format_value(sample_field, an1_idx);

        uint8_t anc1 = extract_format_value(sample_field, an2_idx);

        prev_spos[i * 2] = pos;
        prev_spos[i * 2 + 1] = pos;

        prev_anc[i * 2] = anc0;
        prev_anc[i * 2 + 1] = anc1;
      }

      prev_pos = pos;
      is_first_record = false;
      continue;
    }

    if (chrom_view != cur_chrom) {
      finalize_open_tracts(sample_ids, prev_anc, prev_spos, sample_tracts,
                           cur_chrom, pos - 1);

      cur_chrom.assign(chrom_view);

      for (size_t i = 0; i < n_samples; ++i) {
        std::string_view sample_field = fields[format_idx + 1 + i];

        uint8_t anc0 = extract_format_value(sample_field, an1_idx);

        uint8_t anc1 = extract_format_value(sample_field, an2_idx);

        prev_spos[i * 2] = pos;
        prev_spos[i * 2 + 1] = pos;

        prev_anc[i * 2] = anc0;
        prev_anc[i * 2 + 1] = anc1;
      }

      prev_pos = pos;
      continue;
    }

    for (size_t i = 0; i < n_samples; ++i) {
      size_t idx0 = i * 2;
      size_t idx1 = i * 2 + 1;

      std::string_view sample_field = fields[format_idx + 1 + i];

      uint8_t new_anc0 = extract_format_value(sample_field, an1_idx);

      uint8_t new_anc1 = extract_format_value(sample_field, an2_idx);

      if (new_anc0 != prev_anc[idx0] || new_anc1 != prev_anc[idx1]) {

        uint32_t midpoint = prev_pos + (pos - prev_pos) / 2;

        sample_tracts[i].push_back({cur_chrom, prev_spos[idx0], midpoint,
                                    prev_anc[idx0], prev_anc[idx1]});

        prev_spos[idx0] = midpoint + 1;
        prev_spos[idx1] = midpoint + 1;

        prev_anc[idx0] = new_anc0;
        prev_anc[idx1] = new_anc1;
      }
    }

    prev_pos = pos;
  }

  gzclose(file);

  if (!is_first_record) {
    finalize_open_tracts(sample_ids, prev_anc, prev_spos, sample_tracts,
                         cur_chrom, prev_pos);
  }

  py::dict result;

  std::vector<std::string> samples;
  std::vector<std::string> chroms;
  std::vector<uint32_t> spos_vec;
  std::vector<uint32_t> epos_vec;
  std::vector<int> anc0_vec;
  std::vector<int> anc1_vec;

  size_t total_tracts = 0;

  for (const auto &tracts : sample_tracts)
    total_tracts += tracts.size();

  samples.reserve(total_tracts);
  chroms.reserve(total_tracts);
  spos_vec.reserve(total_tracts);
  epos_vec.reserve(total_tracts);
  anc0_vec.reserve(total_tracts);
  anc1_vec.reserve(total_tracts);

  for (size_t i = 0; i < n_samples; ++i) {
    for (const auto &tract : sample_tracts[i]) {
      samples.push_back(sample_ids[i]);
      chroms.push_back(tract.chrom);
      spos_vec.push_back(tract.spos);
      epos_vec.push_back(tract.epos);
      anc0_vec.push_back(static_cast<int>(tract.anc0));
      anc1_vec.push_back(static_cast<int>(tract.anc1));
    }
  }

  result["sample"] = samples;
  result["chrom"] = chroms;
  result["spos"] = spos_vec;
  result["epos"] = epos_vec;
  result["anc0"] = anc0_vec;
  result["anc1"] = anc1_vec;

  return result;
}

void bind_flare(py::module_ &m) {
  m.def("read_flare", &read_flare, "Read FLARE VCF and return ancestry tracts");
}
