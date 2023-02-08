// Main POV-Ray file for displaying the output of "povwfc"
// Last edited on 2005-08-27 02:27:57 by stolfi

background { color rgb < 0.850,0.870,0.900 > }

global_settings { max_trace_level 8 }

#include "pov-params.inc"

#include "pov-textures.inc"

#declare wfc = union{
  #include modelfile
}

object { wfc }

#include "pov-lights-camera.inc" 
