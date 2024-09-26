cmake_minimum_required(VERSION 3.0.0)
project(bitbirchX VERSION 0.1.0 LANGUAGES C CXX)

include(CTest)
enable_testing()

add_executable(bitbirchX main.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

set(xtensor_DIR "C:\\Users\\Kelly\\mambaforge\\Library\\share\\cmake\\xtensor")
find_package(xtensor REQUIRED)
set(xtensor-blas_DIR "C:\\Users\\Kelly\\mambaforge\\Library\\lib\\cmake\\xtensor-blas")
find_package(xtensor-blas REQUIRED)
target_include_directories(bitbirchX PUBLIC ${xtensor_INCLUDE_DIRS})
target_link_libraries(bitbirchX PUBLIC xtensor)



# cmake_minimum_required(VERSION 3.1)
# project(bitbirchX)

# set(xtl_DIR "C:\\Users\\Kelly\\mambaforge\\Library\\share\\cmake\\xtl")
# find_package(xtl REQUIRED)
# set(xtensor_DIR "C:\\Users\\Kelly\\mambaforge\\Library\\share\\cmake\\xtensor")
# find_package(xtensor REQUIRED)

# add_executable(bitbirchX main.cpp)

# if(MSVC)
#     set(CMAKE_EXE_LINKER_FLAGS /MANIFEST:NO)
# endif()

# target_link_libraries(bitbirchX xtensor::optimize)
# target_link_libraries(bitbirchX xtensor)