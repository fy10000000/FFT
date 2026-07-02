#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <complex>
#include <unordered_map>
#include <cstdint>
#include <algorithm>

constexpr size_t LS4_HEADER_SIZE = 4096;

//////////////////////////////////////////////////////////////
// Simple INI Parser
//////////////////////////////////////////////////////////////

class IniParser
{
public:
  bool load(const std::string& filename)
  {
    std::ifstream file(filename);
    if (!file)
      return false;

    std::string line;
    std::string currentSection;

    while (std::getline(file, line))
    {
      trim(line);
      if (line.empty() || line[0] == '#')
        continue;

      if (line.front() == '[' && line.back() == ']')
      {
        currentSection = line.substr(1, line.size() - 2);
        continue;
      }

      auto pos = line.find('=');
      if (pos == std::string::npos)
        continue;

      std::string key = line.substr(0, pos);
      std::string value = line.substr(pos + 1);

      trim(key);
      trim(value);

      data[currentSection + "." + key] = value;
    }

    return true;
  }

  std::string get(const std::string& section,
    const std::string& key,
    const std::string& defaultVal = "") const
  {
    auto it = data.find(section + "." + key);
    if (it != data.end())
      return it->second;
    return defaultVal;
  }

private:
  std::unordered_map<std::string, std::string> data;

  static void trim(std::string& s)
  {
    s.erase(s.begin(),
      std::find_if(s.begin(), s.end(),
        [](int ch) { return !std::isspace(ch); }));
    s.erase(std::find_if(s.rbegin(), s.rend(),
      [](int ch) { return !std::isspace(ch); }).base(),
      s.end());
  }
};

//////////////////////////////////////////////////////////////
// LabSat Configuration Struct
//////////////////////////////////////////////////////////////

struct LabSatConfig
{
  int channels = 0;
  int quantization = 0;
  int quantA = 0;
  int quantB = 0;
  int sampleRate = 0;
  int shiftFormat = 0;

  void print() const
  {
    std::cout << "Channels      : " << channels << "\n";
    std::cout << "Quantization  : " << quantization << " bit\n";
    std::cout << "Sample Rate   : " << sampleRate << " Hz\n";
    std::cout << "Shift Format  : " << shiftFormat << "\n";
  }
};

//////////////////////////////////////////////////////////////
// 2-bit Decoder
//////////////////////////////////////////////////////////////

inline float decode_2bit(uint8_t bits)
{
  switch (bits & 0x03)
  {
  case 0: return -3.0f / 3.0f;
  case 1: return -1.0f / 3.0f;
  case 2: return  1.0f / 3.0f;
  case 3: return  3.0f / 3.0f;
  }
  return 0.0f;
}

//////////////////////////////////////////////////////////////
// LabSat4 Parser
//////////////////////////////////////////////////////////////

class LabSat4Parser
{
public:
  LabSat4Parser(const std::string& ls4file,
    const LabSatConfig& cfg)
    : config(cfg)
  {
    file.open(ls4file, std::ios::binary);
    if (!file)
      throw std::runtime_error("Cannot open LS4 file");

    file.seekg(LS4_HEADER_SIZE, std::ios::beg);

    validateConfiguration();
  }

  bool readSamples(std::vector<std::complex<float>>& chA,
    std::vector<std::complex<float>>& chB,
    size_t count)
  {
    if (config.quantization == 2 && config.channels == 2)
      return read2bitDual(chA, chB, count);

    printf("Unsupported format"); exit(0);
    return false;
  }

private:
  std::ifstream file;
  LabSatConfig config;

  void validateConfiguration()
  {
    if (config.channels != 2)
      throw std::runtime_error("This parser expects 2 channels");

    if (config.quantization != 2)
      throw std::runtime_error("This version supports 2-bit only");
  }

  bool read2bitDual(std::vector<std::complex<float>>& chA,
    std::vector<std::complex<float>>& chB,
    size_t count)
  {
    chA.resize(count);
    chB.resize(count);

    std::vector<uint8_t> buffer(count);

    if (!file.read(reinterpret_cast<char*>(buffer.data()), count))
      return false;

    for (size_t i = 0; i < count; ++i)
    {
      uint8_t byte = buffer[i];

      uint8_t A_I = (byte >> 6) & 0x03;
      uint8_t A_Q = (byte >> 4) & 0x03;
      uint8_t B_I = (byte >> 2) & 0x03;
      uint8_t B_Q = byte & 0x03;

      chA[i] = {
          decode_2bit(A_I),
          decode_2bit(A_Q)
      };

      chB[i] = {
          decode_2bit(B_I),
          decode_2bit(B_Q)
      };
    }

    return true;
  }
};

//////////////////////////////////////////////////////////////
// Main
//////////////////////////////////////////////////////////////
/* 
int main()
{
  try
  {
    IniParser ini;
    if (!ini.load("X1_L1_X5_20260618T050502Z.ini"))
      throw std::runtime_error("Failed to load INI");

    LabSatConfig cfg;
    cfg.channels = std::stoi(ini.get("config", "CHN"));
    cfg.quantization = std::stoi(ini.get("config", "QUA"));
    cfg.quantA = std::stoi(ini.get("config", "QUAN_A"));
    cfg.quantB = std::stoi(ini.get("config", "QUAN_B"));
    cfg.sampleRate = std::stoi(ini.get("config", "SMP"));
    cfg.shiftFormat = std::stoi(ini.get("config", "SFT"));

    cfg.print();

    LabSat4Parser parser("recording.LS4", cfg);

    std::vector<std::complex<float>> chA, chB;

    const size_t blockSize = 16384;

    while (parser.readSamples(chA, chB, blockSize))
    {
      for (size_t i = 0; i < chA.size(); ++i)
      {
        std::cout
          << chA[i].real() << ","
          << chA[i].imag() << ","
          << chB[i].real() << ","
          << chB[i].imag() << "\n";
      }
    }
  }
  catch (const std::exception& e)
  {
    std::cerr << "Error: " << e.what() << "\n";
  }

  return 0;
}
*/