FUNCTION(PREPEND var prefix)
   SET(listVar "")
   FOREACH(f ${ARGN})
      LIST(APPEND listVar "${prefix}/${f}")
   ENDFOREACH(f)
   SET(${var} "${listVar}" PARENT_SCOPE)
ENDFUNCTION(PREPEND)


include_directories(.)
include_directories("../depends/eigen-git-mirror/")
include_directories("gadgets" "prototools" "examples" "utils")


# TODO: Include Hana
include_directories(
  ${DEPENDS_DIR}/libsnark
  ${DEPENDS_DIR}/libsnark/depends/libfqfft
  ${DEPENDS_DIR}/fmt/include
  )

link_libraries(fmt-header-only)


## == LEGOFW library ==
#add_subdirectory("legoframework")

set(LEGOFW_SRC
        nizk.h nizk.cc
        relation.h relation.cc
        vartypes.h vartypes.cc
        fw_util.h fw_util.cc
        commit.h commit.cc
    )

PREPEND(FULL_LEGOFW_SRC "legoframework" ${LEGOFW_SRC})

#add_library(
 #       legofw STATIC
 #       ${FULL_LEGOFW_SRC}
#)


## == LEGOBASIC library ==

add_subdirectory("gadgets")
add_subdirectory("prototools")
add_subdirectory("utils")

set(GADGETS_SRC
  lipmaa.h lipmaa.cc
  snark.h snark.cc
  subspace.h subspace.cc
  sumcheck.h sumcheck.cc
  hadamardsc.h hadamardsc.cc
  sigma.h sigma.cc
  matrixsc.h matrixsc.cc
  #arithcirc.h arithcirc.cc
)


set(PROTOTOOLS_SRC
  mle.h mle.cc
  commit.h commit.cc
  interp.h interp.cc
  polytools.h polytools.cc
)


set(UTILS_SRC
  matrix.h matrix.cc
  sparsemexp.h sparsemexp.cc
  benchmark.h benchmark.cc
  dbgutil.h dbgutil.cc
  util.h util.cc
  bp_circuits.h bp_circuits.cc
)

PREPEND(FULL_GADGETS_SRC "gadgets" ${GADGETS_SRC})
PREPEND(FULL_PROTOTOOLS_SRC "prototools" ${PROTOTOOLS_SRC})
PREPEND(FULL_UTILS_SRC "utils" ${UTILS_SRC})

add_library(
  legobasic STATIC
  ${FULL_GADGETS_SRC}
  ${FULL_PROTOTOOLS_SRC}
  ${FULL_UTILS_SRC}
 )

add_subdirectory("examples")

