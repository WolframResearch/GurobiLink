cmake_minimum_required(VERSION 3.13)
project(GUROBILink)

string(TIMESTAMP BUILD_TIME "%d %b %Y at %H:%M")

set(GUROBI_VERSION "9.0.2" CACHE STRING "GUROBI version")
set(GUROBI_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/Libraries/Gurobi/${GUROBI_VERSION}" CACHE PATH "GUROBI installation")

set(INSTALL_PDB OFF CACHE BOOL "Include PDB file in paclet layout")

set(BUILD_STUB OFF CACHE BOOL "Build stub paclet")
set(PACLET_VERSION 1.0.1 CACHE STRING "Paclet version")
set(TARGET_WOLFRAM_VERSION 12.2+ CACHE STRING "Target Mathematica version")

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/install" CACHE PATH "Install prefix" FORCE)
endif()

if(WIN32)
	if(CMAKE_SIZEOF_VOID_P EQUAL 8)
		set(SYSTEMID Windows-x86-64)
		set(PLATFORM win64)
		set(QUALIFIER Win64)
		set(GUROBI_LIB gurobi90.lib)
		set(GUROBI_DLL gurobi90.dll)
	else()
		set(SYSTEMID Windows)
		set(PLATFORM win32)
		set(QUALIFIER Win32)
		set(GUROBI_LIB gurobi90.lib)
		set(GUROBI_DLL gurobi90.dll)
		add_definitions(-DMINT_32)
	endif()
	set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd /Zi")
	set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /MTd /Zi")
	set(MATHEMATICA_LAYOUT "C:/Program Files/Wolfram Research/Mathematica/12.2" CACHE PATH "Mathematica installation")
elseif(APPLE)
	set(SYSTEMID MacOSX-x86-64)
	set(PLATFORM mac64)
	set(QUALIFIER Mac64)
	set(GUROBI_LIB libgurobi90.dylib)
	set(GUROBI_LIB_INSTALL_NAME_PATH /Library/gurobi902/mac64/lib)
	set(GUROBI_DLL ${GUROBI_LIB})
	set(MATHEMATICA_LAYOUT "/Applications/Mathematica.app/Contents" CACHE PATH "Mathematica installation")
elseif(UNIX)
	set(SYSTEMID Linux-x86-64)
	set(PLATFORM linux64)
	set(QUALIFIER Lin64)
	set(GUROBI_LIB libgurobi90.so)
	set(GUROBI_DLL ${GUROBI_LIB})
	set(CMAKE_INSTALL_RPATH "\$ORIGIN")
	set(MATHEMATICA_LAYOUT "/usr/local/Wolfram/Mathematica/12.2" CACHE PATH "Mathematica installation")
endif()

if(WIN32)
	set(GUROBI_DLL_PATH "${GUROBI_ROOT}/${PLATFORM}/bin")
else()
	set(GUROBI_DLL_PATH "${GUROBI_ROOT}/${PLATFORM}/lib")
endif()
set(GUROBI_LIBRARY_PATH "${GUROBI_ROOT}/${PLATFORM}/lib")
set(GUROBI_INCLUDE_PATH "${GUROBI_ROOT}/${PLATFORM}/include")

if(NOT WOLFRAM_INCLUDE_PATH)
	set(WOLFRAM_INCLUDE_PATH "${MATHEMATICA_LAYOUT}/SystemFiles/IncludeFiles/C")
endif()

set(GUROBILINK_SOURCES
	CSource/GUROBILink.cpp
	CSource/GUROBISolution.cpp
)

function(copy_wl_sources SRC DEST)
	string(RANDOM temp)
	configure_file(
		${SRC}/PacletInfo.m.in
		${CMAKE_BINARY_DIR}/generated/${temp}/PacletInfo.m
	)
	set(WL_SOURCES
		${SRC}/${PROJECT_NAME}.m
		${CMAKE_BINARY_DIR}/generated/${temp}/PacletInfo.m
		)
	install(FILES
		${WL_SOURCES}
		DESTINATION ${DEST}
		)
endfunction()

if(BUILD_STUB)
	set(STUB_DESTINATION ${CMAKE_INSTALL_PREFIX}/stub CACHE PATH "Stub build destination")
	copy_wl_sources(${CMAKE_CURRENT_SOURCE_DIR}/PacletStub ${STUB_DESTINATION}/${PROJECT_NAME} )
	set(FULL full)
endif()

set(PACLET_DESTINATION ${CMAKE_INSTALL_PREFIX}/${FULL} CACHE PATH "Full build destination")
set(PACLET_LIBDIR ${PACLET_DESTINATION}/${PROJECT_NAME}/LibraryResources/${SYSTEMID})



add_library(${PROJECT_NAME} SHARED ${GUROBILINK_SOURCES})
set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX "")
target_include_directories(${PROJECT_NAME} PRIVATE
	${CMAKE_SOURCE_DIR}/CSource
	${GUROBI_INCLUDE_PATH}
	${WOLFRAM_INCLUDE_PATH}
	)
target_link_directories(${PROJECT_NAME} PRIVATE ${GUROBI_LIBRARY_PATH})
target_link_libraries(${PROJECT_NAME} PRIVATE ${GUROBI_LIB})
target_compile_definitions(${PROJECT_NAME} PRIVATE "$<$<CONFIG:DEBUG>:DEBUG>")

if(APPLE)
	add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
		COMMAND install_name_tool -change ${GUROBI_LIB_INSTALL_NAME_PATH}/${GUROBI_LIB} @loader_path/${GUROBI_LIB} $<TARGET_FILE_NAME:${PROJECT_NAME}>
		COMMENT "Setting rpath for ${GUROBI_LIB}"
	)
endif()

install(TARGETS ${PROJECT_NAME}
		LIBRARY DESTINATION ${PACLET_LIBDIR}
		RUNTIME DESTINATION ${PACLET_LIBDIR}
)

install(FILES
	${GUROBI_DLL_PATH}/${GUROBI_DLL}
	DESTINATION ${PACLET_LIBDIR}
)

if(INSTALL_PDB)
	install(FILES $<TARGET_PDB_FILE:${PROJECT_NAME}> DESTINATION ${PACLET_LIBDIR} OPTIONAL)
endif()

copy_wl_sources(${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME} ${PACLET_DESTINATION}/${PROJECT_NAME})
