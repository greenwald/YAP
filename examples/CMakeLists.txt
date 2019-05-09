# Empty variable to keep track of the `make example` dependencies:
# it will be updated in each of the add_subdirectory.
set(examples_depends)

add_subdirectory(programs)

# Check if $ROOTSYS is set
if(DEFINED ENV{ROOTSYS})

  list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
  # Use the ROOT-provided FindROOT.cmake
   # set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{ROOTSYS}/etc/cmake)

    find_package(ROOT)
    if(ROOT_FOUND)

        # Set the path where FindBAT.cmake is found
        set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmakeModules)

        find_package(BAT)
        if(BAT_FOUND)

            # If ROOT and BAT have been found, add this directory
            add_subdirectory(bat_gen)

        else(BAT_FOUND)
            message(WARNING "BAT not found: bat_gen examples won't be compilable.")
        endif(BAT_FOUND)

    else(ROOT_FOUND)
        message(WARNING "ROOT not found: bat_gen examples won't be compilable.")
    endif(ROOT_FOUND)
else(DEFINED ENV{ROOTSYS})
    message(WARNING "$ROOTSYS is not set: bat_gen examples won't be compilable because I can't find ROOT.")
endif(DEFINED ENV{ROOTSYS})

file(MAKE_DIRECTORY data)
configure_file(${YAP_SOURCE_DIR}/data/evt.pdl ${CMAKE_BINARY_DIR}/data/evt.pdl COPYONLY)

# Create custom 'examples' target
add_custom_target(examples DEPENDS ${examples_depends})
