cmake_minimum_required(VERSION 3.5)
project(bitbirchX VERSION 0.1.0 LANGUAGES C CXX)

set(BLAS_LIBRARIES "C:\\msys64\\mingw64\\lib\\libopenblas.a") # need openblas for xtensor-blas
find_package(BLAS REQUIRED)
set(xtensor_DIR "C:\\Users\\Kelly\\mambaforge\\Library\\share\\cmake\\xtensor")
find_package(xtensor REQUIRED)
set(xtensor-blas_DIR "C:\\Users\\Kelly\\mambaforge\\Library\\lib\\cmake\\xtensor-blas")
find_package(xtensor-blas REQUIRED)

add_executable(bitbirchX src/main.cpp src/Birch.cpp src/abstract_birch.cpp src/_CFNode.cpp src/_CFSubcluster.cpp)

target_include_directories(bitbirchX PUBLIC "C:\\Users\\Kelly\\mambaforge\\Library\\include\\xtensor")
target_include_directories(bitbirchX PUBLIC "C:\\Users\\Kelly\\mambaforge\\Library\\include\\xtensor-blas")
target_link_libraries(bitbirchX PUBLIC xtensor)
target_link_libraries(bitbirchX PUBLIC xtensor-blas)