include_directories(${YAP_SOURCE_DIR}/include ${BAT_INCLUDE_DIRS} ${ROOT_INCLUDE_DIR})

set(SOURCES
    bat_fit.cxx
    bat_gen.cxx
    bat_yap_base.cxx
    fit_base.cxx
    fit_fitFraction.cxx
    mi_fit.cxx
    )

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)
add_custom_target(link-macros
    COMMAND ${CMAKE_COMMAND} -E create_symlink ../../../examples/bat_gen/macros macros)

set(PROGRAMS
    runBatGen
    runBatFit
    runBatFitFitFractions
    runFitFractions)

# Add the program to the `make example` dependencies
set(examples_depends ${examples_depends} ${PROGRAMS} PARENT_SCOPE)

foreach(program ${PROGRAMS})
    add_executable(${program} EXCLUDE_FROM_ALL ${program}.cxx ${SOURCES})
    target_link_libraries(${program} ${ROOT_LIBS} ${BAT_LIBRARIES} YAP)
endforeach(program)
