#ifndef MESS2_ALGORITHM_THREAT_STATIC_HPP
#define MESS2_ALGORITHM_THREAT_STATIC_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <queue>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

#include </usr/include/armadillo>

namespace mess2_algorithms
{
    /**
     * 
     */
    std::vector<double> generate_threat_static(const arma::mat& x_mesh, const arma::mat& y_mesh);

    /**
     * 
     */
    void save_threat_static(const std::vector<double>& threat, const arma::mat& x_mesh, const arma::mat& y_mesh, const std::string& path_file);

} // namespace mess2_algorithms

#endif // MESS2_ALGORITHM_THREAT_STATIC_HPP
