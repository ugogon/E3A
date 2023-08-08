project(e3a)

file(GLOB SRC_FILES *.cpp)

add_library(core STATIC ${SRC_FILES})
target_include_directories(core PUBLIC ${CMAKE_CURRENT_LIST_DIR}/..)

target_link_libraries(core PUBLIC
    Eigen3::Eigen igl::core igl_restricted::triangle gmp
)

