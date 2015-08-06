set(INCLUDE_DIRECTORIES
	${ROOT_INCLUDE_DIR}
	${YAP_SOURCE_DIR}/include)
include_directories(${INCLUDE_DIRECTORIES})

set(YAP_SOURCES
  BlattWeisskopf.cxx
  DataAccessor.cxx
  DataPoint.cxx
  DecayChannel.cxx
  Particle.cxx
  Resonance.cxx
)

make_static_library(YAP "${YAP_SOURCES}")# "${ROOT_LIBS}")
#make_shared_library(YAP "${YAP_SOURCES}")# "${ROOT_LIBS}")