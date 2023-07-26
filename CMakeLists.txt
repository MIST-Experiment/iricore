cmake_minimum_required(VERSION 3.27)

project(irilib)
enable_language(Fortran)


file(GLOB IRI2016_SRC CONFIGURE_DEPENDS "src/iricore/python_interface.f90" "src/iricore/python_stec.f90" "src/iricore/iri2016/*.for")
file(GLOB IRI2020_SRC CONFIGURE_DEPENDS "src/iricore/python_interface.f90" "src/iricore/python_stec.f90" "src/iricore/iri2020/*.for")

add_compile_options(-fpic -O0 -w -g)

add_library(iri2016 SHARED ${IRI2016_SRC})
add_library(iri2020 SHARED ${IRI2020_SRC})

install(TARGETS iri2016 iri2020 LIBRARY DESTINATION irilib)

set(CMAKE_VERBOSE_MAKEFILE OFF)