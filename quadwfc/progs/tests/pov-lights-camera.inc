// Last edited on 2005-08-27 00:53:06 by stolfi

// Lights

light_source {
  100000 * < -3, -3, 5 > 
  color rgb 0.7 * < 1.000, 0.980, 0.960 > 
}

light_source {
  100000 * < 5, -2, 5 > 
  color rgb 0.7 * < 1.000, 0.980, 0.960 > 
}

light_source {
  100000 * < -2, 5, 5 > 
  color rgb 1.0 * < 1.000, 0.980, 0.960 > 
}

// Camera

#local asp = image_width/image_height;
#local dst = 5.0*scene_radius;
#local rad = scene_radius/dst;

camera { 
  location camera_ctr + dst*camera_dir/vlength(camera_dir)
  right   -2.0*rad*sqrt(asp)*x
  up       2.0*rad/sqrt(asp)*y 
  sky z
  look_at camera_ctr
} 
