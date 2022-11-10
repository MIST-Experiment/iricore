cmake_minimum_required(VERSION 3.16)

project(IRI)
enable_language(Fortran)

file(GLOB IRI2016_SRC CONFIGURE_DEPENDS "iri2016/python_interface.f90" "iri2016/*.for")
file(GLOB IRI2020_SRC CONFIGURE_DEPENDS "iri2020/python_interface.f90" "iri2020/*.for")

add_compile_options(-fpic -O0 -w -g)

add_library(iri2016 SHARED ${IRI2016_SRC})
add_library(iri2020 SHARED ${IRI2020_SRC})


set(CMAKE_VERBOSE_MAKEFILE OFF)