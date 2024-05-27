TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += "/media/sf_Tasks/eigen-3.4.0"

SOURCES += \
        eulerpecesolver.cpp \
        main.cpp

HEADERS += \
    eulerpecesolver.h
