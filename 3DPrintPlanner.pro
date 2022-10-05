QT += core gui opengl xml

TARGET = 3DPrintPlanner

TEMPLATE = app

HEADERS += \
    $$PWD\\*.h \
    $$PWD\\GLKLib\\*.h \
    $$PWD\\PQPLib\\*.h \
    $$files($$PWD/QMeshLib/*.h, true) \
    $$files($$PWD/QMeshLib/*.hh, true) \
    $$files($$PWD/Utils/*.h, true) \
    $$files($$PWD/Eigen/*.h, true) \

SOURCES +=  \
    $$PWD\\*.cpp \
    $$PWD\\GLKLib\\*.cpp \
    $$PWD\\PQPLib\\*.cpp \
    $$files($$PWD/QMeshLib/*.cpp, true) \
    $$files($$PWD/QMeshLib/*.cc, true) \
    $$files($$PWD/Utils/*.cpp, true) \
    $$files($$PWD/Eigen/*cpph, true) \

INCLUDEPATH += $$PWD \  # for root category
               $$PWD/QMeshLib/ \
               $$PWD/Eigen/ \
               $$PWD/GLLib/incldue \

LIBS += $$PWD/GLLib/lib/glut32.lib \
        $$PWD/GLLib/lib/glew32.lib \

#LIBS += -fopenmp

CONFIG -= app_bundle

CONFIG += console c++11 c++14
#QMAKE_CXXFLAGS+=/openmp

CONFIG+=debug_and_release
CONFIG(Debug, Debug|Release){DESTDIR = $$PWD/build/debug}
CONFIG(Release, Debug|Release){DESTDIR = $$PWD/build/release}
message(DESTDIR: ($$DESTDIR))

OBJECTS += build\\debug\\SupportStructure.obj \

DISTFILES +=  $$PWD\\SupportStructure.cu \
    
#-------------------------------------------------
# CUDA settings

CUDA_DIR = "C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v10.2"                # Path to cuda toolkit install
SYSTEM_NAME = x64                 # Depending on your system either 'Win32', 'x64', or 'Win64'
SYSTEM_TYPE = 64                    # '32' or '64', depending on your system
CUDA_ARCH = compute_86                 # Type of CUDA architecture
CUDA_CODE = sm_86
NVCC_OPTIONS = --use_fast_math
# include paths
INCLUDEPATH += "$$CUDA_DIR\\include" \
"C:\\ProgramData\\NVIDIA Corporation\\CUDA Samples\\v10.2\\common\\inc"
# library directories
QMAKE_LIBDIR += "$$CUDA_DIR\\lib\\x64"
# The following makes sure all path names (which often include spaces) are put between quotation marks
CUDA_INC = $$join(INCLUDEPATH,'" -I"','-I"','"')
# Add the necessary libraries
CUDA_LIB_NAMES += \cublas \cuda \cudadevrt \cudart \cudart_static \cufft \cufftw \curand \cusolver \cusparse \nppc \nppial \nppicc \nppidei \nppif \nppig \nppim \nppist \nppisu \nppitc \npps \nvblas \nvml \nvrtc \OpenCL \kernel32 \user32 \gdi32 \winspool \comdlg32 \advapi32 \shell32 \ole32 \oleaut32 \uuid \odbc32 \odbccp32 

for(lib, CUDA_LIB_NAMES) {
    CUDA_LIBS += $$lib.lib
}
for(lib, CUDA_LIB_NAMES) {
    NVCC_LIBS += -l$$lib
}
LIBS += $$NVCC_LIBS
# The following library conflicts with something in Cuda
QMAKE_LFLAGS_RELEASE = /NODEFAULTLIB:msvcrt.lib
QMAKE_LFLAGS_DEBUG   = /NODEFAULTLIB:msvcrtd.lib
# MSVCRT link option (static or dynamic, it must be the same with your Qt SDK link option)
MSVCRT_LINK_FLAG_DEBUG   = "/MTd"
MSVCRT_LINK_FLAG_RELEASE = "/MD"
# Configuration of the Cuda compiler
CONFIG(debug, debug|release) {
    # Debug mode
    DESTDIR =  $$PWD\\build\\debug
    OBJECTS_DIR =  $$PWD\\build\\debug\\obj
    CUDA_OBJECTS_DIR =  $$PWD\\build\\debug\\cuda
    cuda_d.input = CUDA_SOURCES
    cuda_d.output = $$CUDA_OBJECTS_DIR\\${QMAKE_FILE_BASE}_cuda.o
    cuda_d.commands = $$CUDA_DIR\\bin\\nvcc.exe -D_DEBUG $$NVCC_OPTIONS $$CUDA_INC $$LIBS \
                      --machine $$SYSTEM_TYPE -arch=$$CUDA_ARCH -code=$$CUDA_CODE \
                      --compile -cudart static -g -DWIN32 -D_MBCS \
                      -Xcompiler "/wd4819,/EHsc,/W3,/nologo,/Od,/Zi,/RTC1" \
                      -Xcompiler $$MSVCRT_LINK_FLAG_DEBUG \
                      -c -o ${QMAKE_FILE_OUT} ${QMAKE_FILE_NAME}
    cuda_d.dependency_type = TYPE_C
    QMAKE_EXTRA_COMPILERS += cuda_d
}
else {
    # Release mode
    DESTDIR =  $$PWD\\build\\release
    OBJECTS_DIR =  $$PWD\\build\\release\\obj
    CUDA_OBJECTS_DIR =  $$PWD\\build\\release\\cuda
    cuda.input = CUDA_SOURCES
    cuda.output = $$CUDA_OBJECTS_DIR\\${QMAKE_FILE_BASE}_cuda.o
    cuda.commands = $$CUDA_DIR\\bin\\nvcc.exe $$NVCC_OPTIONS $$CUDA_INC $$LIBS \
                    --machine $$SYSTEM_TYPE -arch=$$CUDA_ARCH -code=$$CUDA_CODE \
                    --compile -cudart static -D_MBCS \
                    -Xcompiler "/wd4819,/EHsc,/W3,/nologo,/O2,/Zi" \
                    -Xcompiler $$MSVCRT_LINK_FLAG_RELEASE \
                    -c -o ${QMAKE_FILE_OUT} ${QMAKE_FILE_NAME}
    cuda.dependency_type = TYPE_C
    QMAKE_EXTRA_COMPILERS += cuda
}

