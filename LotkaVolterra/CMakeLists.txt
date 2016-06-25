

cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)

project(Lotkavolterra)


if ( win32 )
	include_directories( "C:/Program Files (x86)/GnuWin32/include" 
						 "C:/Users/Takahiro/bin/freeglut-3.0.0/include" )
	link_directories( "C:/Program Files (x86)/GnuWin32/lib"
		              "C:/Users/Takahiro/bin/freeglut-3.0.0/buildGnu/lib" )
endif ( win32 )

if ( UNIX )
	find_library(GLUT_LIBRARY glut)
	find_library(OPENGL_LIBRARY opengl)

	include_directories("/usr/local/include")
	link_directories("/usr/local/lib")
endif ( UNIX )

add_executable(lotkavolterra src/main.cpp)

if ( win32 )
	target_link_libraries(lotkavolterra gsl glu32 freeglut opengl32 )
endif ( win32 )

if ( UNIX )
	target_link_libraries(lotkavolterra gsl ${GLUT_LIBRARY} ${OPENGL_LIBRARY})
endif ( UNIX )

