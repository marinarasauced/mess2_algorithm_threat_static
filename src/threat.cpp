
#include "mess2_algorithm_threat_static/threat.hpp"

namespace mess2_algorithms
{
    std::vector<double> generate_threat_static(const arma::mat& x_mesh, const arma::mat& y_mesh)
    {
        if (x_mesh.n_rows != y_mesh.n_rows || x_mesh.n_cols != y_mesh.n_cols) {
            throw std::invalid_argument("x_mesh and y_mesh must have the same dimensions.");
        }

        const int64_t n_rows = x_mesh.n_rows;
        const int64_t n_cols = x_mesh.n_cols;
        const int64_t n_elem = x_mesh.n_elem;

        std::vector<double> threat;
        threat.reserve(n_elem);

        arma::mat matrix = arma::mat(n_rows, n_cols, arma::fill::zeros);

        const int16_t n_peaks = 11;
        arma::mat coeff_peaks = {
            { 0.4218, 1.4253, 0.8271,  1.5330,  1.7610,  0.4533,  0.2392,  0.8364,  0.5060, 1.5190,  2},
            {-4.7996, 6.8744, 1.6560,  4.3881, -3.1295, -9.8145, -6.1511,  0.1478, -9.5152, 1.0008,  9},
            {-4.4763, 3.9248, 7.6271, -9.5064, -3.1759, -1.5719, -8.3998, -8.4129, -8.5525, 8.0068, -1},
            { 3.2579, 1.5239, 1.2908,  2.0099,  2.7261,  2.7449,  2.9398,  2.7439,  0.3691, 3.1097,  3},
            { 0.4039, 0.4382, 2.4844,  1.9652,  1.9238,  1.8567,  0.5470,  1.0401,  0.7011, 3.3193,  3}
        };

        for (int16_t iter = 0; iter < n_peaks; ++iter) {
            arma::mat c_xym = exp(
                -pow((x_mesh - coeff_peaks(1, iter)), 2) / (2 * pow(coeff_peaks(3, iter), 2))
                -pow((y_mesh - coeff_peaks(2, iter)), 2) / (2 * pow(coeff_peaks(4, iter), 2))
            );
            matrix += coeff_peaks(0, iter) * c_xym;
        }
        const double c1 = 5;
        matrix = c1 * matrix + 1;

        for (int64_t iter = 0; iter < n_rows; ++iter) {
            for (int64_t jter = 0; jter < n_cols; ++jter) {
                threat.push_back(matrix(iter, jter));
            }
        }

        return threat;
    }

    void save_threat_static(const std::vector<double>& threat, const arma::mat& x_mesh, const arma::mat& y_mesh, const std::string& path_file)
    {
        if (x_mesh.n_elem != y_mesh.n_elem || x_mesh.n_elem != threat.size()) {
            throw std::invalid_argument("x_mesh, y_mesh, and threat must have the same number of elements");
        }

        std::string path_new;
        if (!path_file.empty() && path_file[0] == '~') {
            const char* home = getenv("HOME");
            if (home) {
                path_new = std::string(home) + path_file.substr(1);
            } else {
                throw std::runtime_error("could not determine the home directory");
            }
        }

        std::ofstream file(path_new);
        if (!file.is_open()) {
            throw std::runtime_error("could not open file for writing");
        }

        const int64_t n_rows = x_mesh.n_rows;
        const int64_t n_cols = x_mesh.n_cols;

        for (int64_t iter = 0; iter < n_rows; ++iter) {
            for (int64_t jter = 0; jter < n_cols; ++jter) {
                int64_t index = iter * n_cols + jter;
                file << x_mesh(iter, jter) << "," << y_mesh(iter, jter) << "," << threat[index] << "\n";
            }
        }

        file.close();
    }

} // namespace mess2_algorithms
