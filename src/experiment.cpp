#include <iostream>
#include <vector>
#include <map>
#include <omp.h>
#include <string>
#include "isd.hpp"
#include "binop.hpp"
#include "lcode.hpp"
#include "misc.hpp"

bool run_algorithm(ISD& isd, const std::vector<size_t>& pi, const std::string& output_dir) {
    
    std::string name = isd.get_name();
    std::map<size_t, std::pair<double, long double>> stats;

    // Run the algorithm to find errors
    binop::binvec e;
    isd.find_min_error(e, pi, stats);

    // Print the smallest error weight
    std::cout << "Smallest error weight: " << e.count() << ". ";

    // Obtain predictions
    std::vector<size_t> profile;
    if (isd.get_name() == "LBB" || isd.get_name() == "SternBabai" || isd.get_name() == "WagnerBabai") {
        profile = isd.get_C().get_profile();
    }

    std::vector<long double> predictions = isd.errors_dist_predict(profile);
    for (size_t i = 0; i < predictions.size(); ++i) {
        stats[i].second = predictions[i];
    }

    // Write the statistics to a file
    const std::string output_file = "experiment_" + name + "_" + std::to_string(binop::n) + ".txt";
    try {
        write_statistics_to_file(name, output_dir, output_file, stats);
        std::cout << "Statistics written to " << output_dir + "/" + output_file << "\n";
    }
    catch (const std::ios_base::failure& e) {
        std::cerr << "Error writing statistics to file: " << e.what() << "\n";
    }

    return true;
}

int main() {
    
    // Parameters for the experiments
    size_t k = binop::n / 2;
    size_t w2 = 3;
    size_t a = 2;

    // Obtain optimal c for Stern
    size_t c_stern = static_cast<size_t>(std::log2(static_cast<long double>(binop::util::count_combinations(k / 2, w2))));
    std::cout << "Optimal c for Stern: " << c_stern << "\n\n";

    // Obtain optimal c for Wagner
    size_t num_lists = static_cast<size_t>(pow(2, a));
    size_t c_wagner = static_cast<size_t>(std::log2(static_cast<long double>(binop::util::count_combinations(k / num_lists, w2))));
    c_wagner *= a;
    std::cout << "Optimal c for Wagner: " << c_wagner << "\n\n";

    // Create the data directory
    const std::string data_dir = "./data";
    #ifdef _WIN32
        if (_mkdir(data_dir.c_str()) != 0 && errno != EEXIST) {
            throw std::runtime_error("Failed to create directory: " + data_dir);
        }
    #else
        if (mkdir(dir.c_str(), 0755) != 0 && errno != EEXIST) {
            throw std::runtime_error("Failed to create directory: " + data_dir);
        }
    #endif

    // Directory for output files
    const std::string output_dir = "./data/n" + std::to_string(binop::n) 
        + "_k" + std::to_string(k)
        + "_c" + std::to_string(c_stern);

    // Generate a shared instance of the problem
    binop::binmat G;
    binop::util::gen_rnd(k, G);
    LCode code(G);

    binop::binvec t;
    binop::util::gen_rnd(t);

    // Generate permutation
    std::vector<size_t> pi;
    binop::util::gen_rnd_perm(pi);

    // Profile
    std::vector<size_t> profile;

    /************************* Wagner ***************************/

    // Run Wagner-Babai experiment on the shared problem instance
    code.set_G(G);
    Wagner wagnerbabai(code, t, w2, c_wagner, a, true);
    run_algorithm(wagnerbabai, pi, output_dir);

    // Run Wagner experiment on the shared problem instance
    code.set_G(G);
    Wagner wagner(code, t, w2, c_wagner, a);
    run_algorithm(wagner, pi, output_dir);

    /************************* Stern ***************************/

    // Run Stern-Babai experiment on the shared problem instance
    code.set_G(G);
    Stern sternbabai(code, t, w2, c_stern, true);
    run_algorithm(sternbabai, pi, output_dir);

    // Run Stern experiment on the shared problem instance
    code.set_G(G);
    Stern stern(code, t, w2, c_stern);
    run_algorithm(stern, pi, output_dir);

    /*********************** Lee-Brickell ***********************/

    // Run Lee-Brickell-Babai experiment on the shared problem instance
    code.set_G(G);
    LeeBrickell lbb(code, t, w2, true);
    run_algorithm(lbb, pi, output_dir);

    // Run Lee-Brickell experiment on the shared problem instance
    code.set_G(G);
    LeeBrickell lb(code, t, w2);
    run_algorithm(lb, pi, output_dir);

    return 0;
}

