add_executable(hadamard hadamard.cc)
target_link_libraries(hadamard snark legobasic)

add_executable(cplink cplink.cc)
target_link_libraries(cplink snark legobasic)

add_executable(legogrothmatrix legogrothmatrix.cc)
target_link_libraries(legogrothmatrix snark legobasic)

add_executable(matrixsc matrixsc.cc)
target_link_libraries(matrixsc snark legobasic)


#add_executable(matrixAC matrixAC.cc)
#target_link_libraries(matrixAC snark legobasic)