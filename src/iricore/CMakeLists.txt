cmake_minimum_required(VERSION 3.16)

project(IRI)
enable_language(Fortran)

file(GLOB IRI2016_SRC CONFIGURE_DEPENDS "iri2016/*.f90" "src_iri2016/*.for")

add_compile_options(-fpic -O3 -w)
add_library(iri2016 SHARED ${IRI2016_SRC})

#[[
file(GLOB IRI2020_SRC CONFIGURE_DEPENDS "iri2020/*.f90" "src_iri2020/*.for")
add_library(iri2020 SHARED ${IRI2020_SRC})
]]

set(CMAKE_VERBOSE_MAKEFILE OFF)