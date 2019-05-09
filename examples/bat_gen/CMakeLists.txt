include_directories(${YAP_SOURCE_DIR}/include ${BAT_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS})

set(SOURCES
    bat_fit.cxx
    bat_gen.cxx
    bat_yap_base.cxx
    fit_fitFraction.cxx)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)
add_custom_target(link-macros
    COMMAND ${CMAKE_COMMAND} -E create_symlink ../../../examples/bat_gen/macros macros)

set(PROGRAMS
  runBtoDpipi0
    runBatGen
    runBatFit
    runBatFitFitFractions
    runFitFractions)

# Add the program to the `make example` dependencies
set(examples_depends ${examples_depends} ${PROGRAMS} PARENT_SCOPE)

file(MAKE_DIRECTORY data)
configure_file(${YAP_SOURCE_DIR}/data/evt.pdl ${CMAKE_CURRENT_BINARY_DIR}/data/evt.pdl COPYONLY)

foreach(program ${PROGRAMS})
    add_executable(${program} EXCLUDE_FROM_ALL ${program}.cxx ${SOURCES})
    target_link_libraries(${program} ${ROOT_LIBRARIES} ${BAT_LIBRARIES} YAP)
endforeach(program)
