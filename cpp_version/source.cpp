/* ************************************************************************ */
/* Coded by Takuro Tokunaga                                                 */
/* Two-dimensional heat conduction equation solved by Finite Element Method */
/* About this code:                                                         */
/* Liner interpolation                                                      */
/* Required files:                                                          */
/* 1. rectangle1.msh: number of coord and nord                              */
/* 2. rectangle2.msh: cord information                                      */
/* 3. rectangle3.msh: nord information                                      */
/* 4. parameters.txt                                                        */
/* 5. rectangle-t.bc                                                        */
/* 6. rectangle-num.bc                                                      */
/* 7. rectangle-n.bc                                                        */
/* Updated: September 20, 2026                                              */
/* ************************************************************************ */

/* g++-13 -std=c++20 -Iinclude source.cpp Input.cpp -o bmt */
#include "Input.hpp"
int main(){
    Input calculation("conduction");
    calculation.displayParameters();
    return 0;
}