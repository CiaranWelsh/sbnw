//
// Created by cwelsh on 2/22/2020.
//

#ifndef SBNW_EXPORT_H
#define SBNW_EXPORT_H

#define BUILD_SHARED

#ifdef  BUILD_SHARED
/*Enabled as "export" while compiling the dll project*/
    #define EXPORT __declspec(dllexport)
#else
/*Enabled as "import" in the Client side for using already created dll file*/
#define EXPORT __declspec(dllimport)
#endif

#endif //SBNW_EXPORT_H
