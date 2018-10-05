include_directories(
	${YAP_SOURCE_DIR}/include 
	${BAT_INCLUDE_DIRS} 
	${ROOT_INCLUDE_DIRS}
	"/nfs/freenas/tuph/e18/project/belle/users/jrauch/SmokeyJoesBelle/analysis/include"
	)

set(SOURCES
    bat_fit.cxx
    bat_gen.cxx
    bat_yap_base.cxx
    fit_fitFraction.cxx)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)
add_custom_target(link-macros
    COMMAND ${CMAKE_COMMAND} -E create_symlink ../../../examples/bat_gen/macros macros)

set(PROGRAMS
	a1_integrals
	drawMassShapes
    runBatGen
    runBatFit
    runBatFit_4pi
    runBatFit_4pi_FF
    runBatFitFitFractions
    runFitFractions
    runRootGen_4pi
    runWeighEvents)

# Add the program to the `make example` dependencies
set(examples_depends ${examples_depends} ${PROGRAMS} PARENT_SCOPE)

file(MAKE_DIRECTORY data)
configure_file(${YAP_SOURCE_DIR}/data/evt.pdl ${CMAKE_CURRENT_BINARY_DIR}/data/evt.pdl COPYONLY)

foreach(program ${PROGRAMS})
    add_executable(${program} EXCLUDE_FROM_ALL ${program}.cxx ${SOURCES})
    target_link_libraries(${program} ${ROOT_LIBS} ${BAT_LIBRARIES} YAP
    	/nfs/freenas/tuph/e18/project/belle/users/jrauch/SmokeyJoesBelle/analysis/lib/libanalysis.so
    )
endforeach(program)
