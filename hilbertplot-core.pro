#-------------------------------------------------
#
# Project created by QtCreator 2020-06-08T15:05:19
#
#-------------------------------------------------


QT -= gui
QT += core
CONFIG += c++11 debug
TARGET = hilbertplot-core
TEMPLATE = lib

DEFINES += HILBERTPLOT_LIBRARY
DEFINES += QT_CORE


# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

INCLUDEPATH += ./headers \
               ../common

SOURCES += \
        src/hilbertcurve.cpp \
    src/hpoint.cpp \
    src/datasequence.cpp \
    src/hilbertplot.cpp

HEADERS += \
        headers/hilbertcurve.h \
        headers/hpoint.h \
        headers/hilbertdefines.h \
        headers/hilbertdefines.h \
        headers/datasequence.h \
        headers/hilbertplot.h
        ../common/threads_utility.h

LIBS += -lfftw3

unix {
    target.path = /usr/lib
    INSTALLS += target
}
