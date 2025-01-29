/* Authors: Manoj Kumar <manojkumarsivaraj334@gmail.com>
          Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
          Joel Nothman <joel.nothman@gmail.com>
 License: BSD 3 clause */

#pragma once
#include <iostream>
#include <string>
#include <cmath>
#include <ctime>
#include <vector>
#include <cstdlib>
#include <tuple>
#include <utility>
#include <fstream>
#include <xtensor/xarray.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xsort.hpp> //for argmin
#include <xtensor/xview.hpp> //for view
#include <xtensor-blas/xlinalg.hpp> //for dot
#include <xtensor/xrandom.hpp> //for randint
#include "_CFNode.h"
#include "_CFSubcluster.h"
#include "Birch.h"
class _CFNode;
class _CFSubcluster;

xt::xarray<float> find_centroid(xt::xarray<float> X);
float isim(xt::xarray<float> X);
xt::xarray<float> sim_search(xt::xarray<float> point, xt::xarray<float> X);
float jt_isim(xt::xarray<float> c_total, int n_objects);
std::tuple<std::pair<int, int>, xt::xarray<float>, xt::xarray<float>> max_separation(xt::xarray<float> X);
xt::xarray<float> calc_centroid(xt::xarray<float> linear_sum, int n_samples, std::string sim_index = "JT");
std::pair<_CFSubcluster*, _CFSubcluster*> _split_node(_CFNode* node, double threshold, int branching_factor);
float cluster_radius(xt::xarray<float> cluster);