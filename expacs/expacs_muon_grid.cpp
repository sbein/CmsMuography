#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits.h>
#include <stdexcept>
#include <string>
#include <unistd.h>

double getHPcpp(int, int, int);
double getrcpp(double, double);
double getdcpp(double, double);
double getSpecCpp(int, double, double, double, double, double);
double getSpecAngFinalCpp(int, double, double, double, double, double, double);

namespace {

struct Options {
  int year = 2026;
  int month = 8;
  int day = 21;
  int ne = 120;
  int ncos = 100;
  double lat = 46.310819;
  double lon = 6.078311;
  double altitude_km = 0.510;
  double ground_water_fraction = 0.15;
  double emin_mev = 10.0;
  double emax_mev = 1.0e7;
  double w_index = -1.0;
  std::string expacs_dir;
  std::string output = "expacs_muon_flux_grid.csv";
};

double parse_double(const char* text) {
  char* end = nullptr;
  const double value = std::strtod(text, &end);
  if (!end || *end != '\0') {
    throw std::runtime_error(std::string("bad floating point value: ") + text);
  }
  return value;
}

int parse_int(const char* text) {
  char* end = nullptr;
  const long value = std::strtol(text, &end, 10);
  if (!end || *end != '\0') {
    throw std::runtime_error(std::string("bad integer value: ") + text);
  }
  return static_cast<int>(value);
}

Options parse_options(int argc, char** argv) {
  Options options;
  for (int iarg = 1; iarg < argc; ++iarg) {
    const std::string key = argv[iarg];
    auto need_value = [&]() -> const char* {
      if (iarg + 1 >= argc) {
        throw std::runtime_error("missing value for " + key);
      }
      return argv[++iarg];
    };
    if (key == "--year") options.year = parse_int(need_value());
    else if (key == "--month") options.month = parse_int(need_value());
    else if (key == "--day") options.day = parse_int(need_value());
    else if (key == "--lat") options.lat = parse_double(need_value());
    else if (key == "--lon") options.lon = parse_double(need_value());
    else if (key == "--altitude-km") options.altitude_km = parse_double(need_value());
    else if (key == "--ground-water-fraction") options.ground_water_fraction = parse_double(need_value());
    else if (key == "--w-index") options.w_index = parse_double(need_value());
    else if (key == "--expacs-dir") options.expacs_dir = need_value();
    else if (key == "--emin-mev") options.emin_mev = parse_double(need_value());
    else if (key == "--emax-mev") options.emax_mev = parse_double(need_value());
    else if (key == "--ne") options.ne = parse_int(need_value());
    else if (key == "--ncos") options.ncos = parse_int(need_value());
    else if (key == "--output") options.output = need_value();
    else if (key == "--help" || key == "-h") {
      std::cout
        << "Usage: expacs_muon_grid [options]\n"
        << "  --output PATH\n"
        << "  --lat DEG --lon DEG --altitude-km KM\n"
        << "  --year YYYY --month M --day D\n"
        << "  --w-index W  optional override for PARMA solar W-index\n"
        << "  --expacs-dir DIR  chdir to the unpacked EXPACS/PARMA directory before loading input tables\n"
        << "  --emin-mev E --emax-mev E --ne N --ncos N\n";
      std::exit(0);
    } else {
      throw std::runtime_error("unknown option: " + key);
    }
  }
  if (options.ne <= 0 || options.ncos <= 0) {
    throw std::runtime_error("--ne and --ncos must be positive");
  }
  if (options.emin_mev <= 0.0 || options.emax_mev <= options.emin_mev) {
    throw std::runtime_error("bad energy range");
  }
  return options;
}

std::string absolute_path(const std::string& path) {
  if (!path.empty() && path[0] == '/') {
    return path;
  }
  char cwd[PATH_MAX];
  if (!getcwd(cwd, sizeof(cwd))) {
    throw std::runtime_error("could not determine current working directory");
  }
  return std::string(cwd) + "/" + path;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Options options = parse_options(argc, argv);
    const std::string output_path = absolute_path(options.output);
    if (!options.expacs_dir.empty() && chdir(options.expacs_dir.c_str()) != 0) {
      throw std::runtime_error("could not chdir to EXPACS directory: " + options.expacs_dir);
    }
    const int muon_plus = 29;
    const int muon_minus = 30;
    const int angular_family = 4;

    const double solar = options.w_index >= 0.0
      ? options.w_index
      : getHPcpp(options.year, options.month, options.day);
    const double rigidity = getrcpp(options.lat, options.lon);
    const double depth = getdcpp(options.altitude_km, options.lat);
    const double log_emin = std::log10(options.emin_mev);
    const double log_step = (std::log10(options.emax_mev) - log_emin) / options.ne;
    const double cos_step = 1.0 / options.ncos;

    std::ofstream out(output_path);
    if (!out) {
      throw std::runtime_error("could not create " + output_path);
    }
    out << "# EXPACS/PARMA muon flux grid from parma_cpp.zip\n";
    out << "# EXPACS URL: https://phits.jaea.go.jp/expacs/\n";
    out << "# Cite: T. Sato, PLOS ONE 10(12): e0144679 (2015), doi:10.1371/journal.pone.0144679\n";
    out << "# Cite: T. Sato, PLOS ONE 11(8): e0160390 (2016), doi:10.1371/journal.pone.0160390\n";
    out << "# lat_deg=" << options.lat << " lon_deg=" << options.lon
        << " altitude_km=" << options.altitude_km << "\n";
    out << "# year=" << options.year << " month=" << options.month
        << " day=" << options.day << " g=" << options.ground_water_fraction << "\n";
    out << "# W_index=" << solar << " cutoff_rigidity_GV=" << rigidity
        << " atmospheric_depth_g_cm2=" << depth << "\n";
    out << "# flux units: cm^-2 s^-1 MeV^-1 sr^-1, summed over mu+ and mu-\n";
    out << "energy_low_mev,energy_high_mev,energy_mid_mev,cos_low,cos_high,cos_mid,flux\n";
    out << std::scientific << std::setprecision(8);

    for (int ie = 0; ie < options.ne; ++ie) {
      const double e_low = std::pow(10.0, log_emin + log_step * ie);
      const double e_high = std::pow(10.0, log_emin + log_step * (ie + 1));
      const double e_mid = std::sqrt(e_low * e_high);
      const double spectral =
        getSpecCpp(muon_plus, solar, rigidity, depth, e_mid, options.ground_water_fraction) +
        getSpecCpp(muon_minus, solar, rigidity, depth, e_mid, options.ground_water_fraction);
      for (int ic = 0; ic < options.ncos; ++ic) {
        const double cos_low = cos_step * ic;
        const double cos_high = cos_step * (ic + 1);
        const double cos_mid = 0.5 * (cos_low + cos_high);
        const double angular =
          getSpecAngFinalCpp(angular_family, solar, rigidity, depth, e_mid,
                             options.ground_water_fraction, cos_mid);
        out << e_low << "," << e_high << "," << e_mid << ","
            << cos_low << "," << cos_high << "," << cos_mid << ","
            << spectral * angular << "\n";
      }
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << "\n";
    return 1;
  }
}
