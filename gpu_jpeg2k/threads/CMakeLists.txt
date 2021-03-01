set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS})

#SET(CMAKE_CXX_FLAGS "-std=c++0x -lpthread")
cuda_add_library(threads
ThreadPool.cpp
thpool.c
)
TARGET_LINK_LIBRARIES(threads 
	-pthread
)
