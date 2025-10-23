#version 3.7;
//EXAMPLE OF A "CUBE BUILDING"

#include "colors.inc"
#include "woods.inc"
#include "stones.inc"
#include "metals.inc"
#include "golds.inc"
#include "glass.inc"

global_settings {
    assumed_gamma 1.0
}


//Place the camera
camera {
  sky <0,0,1>          //Don't change this
  direction <-1,0,0>   //Don't change this  
  location  <10,10,5>  //Change this to move the camera to a different point
  look_at   <0,0,0>    //Change this to aim the camera at a different point
  right <-4/3,0,0>     //Don't change this
  angle 30
}

//Place a light
light_source {
  <10,-10,10>  //Change this if you want to put the light at a different point
  color White*3
}

//Set a background color
background { color White*2 }

//Create a "floor"
plane {
  <0,0,-1>,0
  texture { T_Silver_1A }
}

//Create a box that extends between the 2 specified points
#declare mycube = box {
  <0,0,0>  // one corner position <X1 Y1 Z1>
  <1,1,1>  // other corner position <X2 Y2 Z2>
}

//Change cubes and their locations below this point.

object { mycube texture {T_Stone32} }
object { mycube translate <1,0,0> texture {T_Stone32} }
object { mycube translate <0,1,0> texture {T_Stone32} }
object { mycube translate <0,2,0> texture {T_Stone32} }
object { mycube translate <0,1,1> texture {T_Stone32} }