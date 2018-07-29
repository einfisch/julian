cmake_minimum_required (VERSION 2.7)
project (ModSim2)



link_directories("/mnt/c/Users/julian/Documents/info/bachelor/opencmiss_sdk/lib/")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
## targets

#add_executable(fiber_angles_exe fiber_angles.cpp)
#add_executable(exfile_reader_exe exfile_reader.cpp)
#add_executable(vtu_write_exe vtu_write.cpp)
#add_executable(ptoe_exe prolate_spheroidal_to_euclidian.cpp)
#add_executable(testexe vtu_converter.cpp)
add_executable(xd mesh_generator.cpp)


#target_include_directories(fiber_angles_exe PRIVATE "/mnt/c/Users/julian/Documents/info/bachelor/heart/cpp")
#target_include_directories(fiber_angles_exe PRIVATE "/mnt/c/Users/julian/Documents/info/bachelor/opencmiss_sdk/include")	
#target_include_directories(fiber_angles_exe PRIVATE "/mnt/c/Users/julian/Documents/info/bachelor/opencmiss_sdk/include/opencmiss")
#target_include_directories(exfile_reader_exe PRIVATE "/mnt/c/Users/julian/Documents/info/bachelor/opencmiss_sdk/include")	
#target_include_directories(exfile_reader_exe PRIVATE "/mnt/c/Users/julian/Documents/info/bachelor/opencmiss_sdk/include/opencmiss")
target_include_directories(xd PRIVATE "/mnt/c/Users/julian/Documents/info/bachelor/opencmiss_sdk/include")	
target_include_directories(xd PRIVATE "/mnt/c/Users/julian/Documents/info/bachelor/opencmiss_sdk/include/opencmiss")

#target_link_libraries(fiber_angles_exe zinc )
#target_link_libraries(exfile_reader_exe zinc )
target_link_libraries(xd zinc )
target_link_libraries(xd zinc )
