#include "Birch.h"
#include "abstract_birch.h"
#include "_CFSubcluster.h"
#include "_CFNode.h"
#include <xtensor/xnpy.hpp> // for testing

int main() {
    // Simple example running random sets with 1000 to 50000 molecules
    std::string z;
    for (int n = 1000; n < 2001; n += 1000) {
        std::cout << n << std::endl;
        // xt::random::seed(0); // for testing
        auto dat = xt::random::randint<int>({n, 100}, 0, 2);
        // std::string filename = "../bitbirch/development-code/bitbirchX/random_data_" + std::to_string(n) + ".npy";
        // xt::xarray<int64_t> dat = xt::load_npy<int64_t>(filename);
        // std::cout << dat << std::endl;
        auto brc = Birch(0.50, 50);
        auto v = std::chrono::high_resolution_clock::now();
        brc.fit(dat);
        auto leaves = brc._get_leaves();
        // for (_CFNode* leave : leaves) {
        //     for (_CFSubcluster* subcluster : leave->subclusters_) {
        //         std::cout << subcluster->mol_indices.size() << std::endl;
        //     }
        // }
        std::ostringstream oss;
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = end - v;
        oss << std::setw(10) << n << "    " << std::fixed << std::setprecision(6) << std::setw(10) << std::setfill('0') << std::setfill(' ') << time.count() << std::endl;
        z += oss.str();
    }

    std::ofstream outfile("../bitbirch/development-code/bitbirchX/jt_fit_label.txt");
    if(outfile.is_open()) {
        outfile << z;
        outfile.close();
    }
    else {
        std::cout << "Unable to open file" << std::endl;
    }

    return 0;
}