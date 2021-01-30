QT += core gui opengl xml

TARGET = 3DPrintPlanner

TEMPLATE = app

HEADERS += \
    $$PWD/*.h \
    $$PWD/GLKLib/*.h \
    $$PWD/PQPLib/*.h \
    $$files($$PWD/QMeshLib/*.h, true) \
    $$files($$PWD/QMeshLib/*.hh, true) \
    $$files($$PWD/Utils/*.h, true) \
    $$files($$PWD/Eigen/*.h, true) \

SOURCES +=  \
    $$PWD/*.cpp \
    $$PWD/GLKLib/*.cpp \
    $$PWD/PQPLib/*.cpp \
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
