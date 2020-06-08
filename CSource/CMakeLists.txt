cmake_minimum_required(VERSION 3.0)
project(GUROBILink)

set(GUROBI_VERSION "9.0" CACHE STRING "GUROBI version")
#set(GUROBI_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../Libraries/GUROBI/${GUROBI_VERSION}/tools/platform" CACHE PATH "GUROBI installation")
set(GUROBI_ROOT "C:/gurobi902" CACHE PATH "GUROBI installation")

set(MATHEMATICA_LAYOUT "C:/Program Files/Wolfram Research/Mathematica/12.2" CACHE PATH "Mathematica installation")
set(INSTALL_PDB OFF CACHE BOOL "Include PDB file in paclet layout")

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/install" CACHE PATH "Install prefix" FORCE)
endif()

if(WIN32)
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(SYSTEMID Windows-x86-64)
        set(PLATFORM win64)
        set(GUROBI_LIB gurobi90.lib)
        set(GUROBI_DLL gurobi90.dll)
    else()
        set(SYSTEMID Windows)
        set(PLATFORM win32)
        set(GUROBI_LIB gurobi90.dll)
        add_definitions(-DMINT_32)
    endif()
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd /Zi")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /MTd /Zi")
elseif(APPLE)
    set(SYSTEMID MacOSX-x86-64)
    set(PLATFORM NOT_FOUND)
    set(GUROBI_LIB NOT_FOUND)
elseif(UNIX)
    set(SYSTEMID Linux-x86-64)
    set(PLATFORM NOT_FOUND)
    set(GUROBI_LIB NOT_FOUND)
    set(CMAKE_INSTALL_RPATH "\$ORIGIN")
endif()

set(GUROBI_BIN_PATH "${GUROBI_ROOT}/${PLATFORM}/bin")
set(GUROBI_LIBRARY_PATH "${GUROBI_ROOT}/${PLATFORM}/lib")
set(GUROBI_INCLUDE_PATH "${GUROBI_ROOT}/${PLATFORM}/include")

if(NOT WOLFRAM_INCLUDE_PATH)
    set(WOLFRAM_INCLUDE_PATH "${MATHEMATICA_LAYOUT}/SystemFiles/IncludeFiles/C")
endif()

include_directories(${CMAKE_SOURCE_DIR} ${GUROBI_INCLUDE_PATH} ${WOLFRAM_INCLUDE_PATH})
link_directories(${GUROBI_LIBRARY_PATH})
link_libraries(${GUROBI_LIB})

add_compile_options("$<$<CONFIG:DEBUG>:-DDEBUG>")

set(GUROBILINK_SOURCES
    GUROBILink.cpp
    GUROBISolution.cpp
)

set(WL_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../${PROJECT_NAME})
set(WL_SOURCES
    ${WL_ROOT}/GUROBILink.m
    ${WL_ROOT}/PacletInfo.m
)
    

add_library(${PROJECT_NAME} SHARED ${GUROBILINK_SOURCES})
set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX "")

if(APPLE)
    add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND install_name_tool -change ${GUROBI_LIB} @loader_path/${GUROBI_LIB} $<TARGET_FILE_NAME:${PROJECT_NAME}>
        COMMENT "Setting rpath for ${GUROBI_LIB}"
    )
endif()

set(PACLET_ROOT ${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME})
set(PACLET_LIBDIR ${PACLET_ROOT}/LibraryResources/${SYSTEMID})

install(TARGETS ${PROJECT_NAME}
        LIBRARY DESTINATION ${PACLET_LIBDIR}
        RUNTIME DESTINATION ${PACLET_LIBDIR}
)
install(FILES 
    ${GUROBI_BIN_PATH}/${GUROBI_DLL}
    DESTINATION ${PACLET_LIBDIR}
)
install(FILES ${WL_SOURCES} DESTINATION ${PACLET_ROOT})
if(INSTALL_PDB)
    install(FILES $<TARGET_PDB_FILE:${PROJECT_NAME}> DESTINATION ${PACLET_LIBDIR} OPTIONAL)
endif()
