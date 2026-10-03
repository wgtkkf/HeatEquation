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
/* Updated: September 27, 2026                                              */
/* ************************************************************************ */
/* For functions */
#include "../include/Parameter.hpp"
#include "../include/Mesh.hpp"

/*  */
#include <format>   /* for format*/
#include <iostream> /* for std cout*/

int main(){
    ParameterInput param("Parameter information"); /* instance of a class */
    MeshInput mesh("Mesh information");            /* instance of a class */

    param.displayParameters();
    mesh.displayParameters();
    param.readParameters("../inputs/parameters.txt"); 
    mesh.readMeshParameters("../inputs/rectangle1.msh"); 
    mesh.readMeshCoordinates("../inputs/rectangle2.msh"); 
    
    // display parameters    
    std::cout << std::format("--- Loaded Parameters ---\n");
    std::cout << std::format("Thermal conductivity: {}\n", param.thermal_conductivity);
    std::cout << std::format("Density: {}\n", param.density);
    std::cout << std::format("Specific Heat: {}\n", param.specific_heat);

    // display parameters    
    std::cout << std::format("--- Loaded Mesh Parameters ---\n");    
    std::cout << std::format("Number of Nodes: {}\n", mesh.np);    
    std::cout << std::format("Number of Elements: {}\n", mesh.ne);    
    std::cout << std::format("Time Step: {}\n", mesh.dt);    

    // display coordinates
    mesh.displayCoordinates(mesh.np);
    
    return 0;
}