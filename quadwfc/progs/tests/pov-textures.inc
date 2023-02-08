// Textures for povwfc
// Last edited on 2005-08-26 16:42:08 by stolfi

// Creates a texture of color {clr} and transmittance {trs}
#macro txc(clr, trs)
  texture {
    pigment { color rgb clr transmit trs }
    finish { diffuse 0.7 ambient 0.3 }
  }
#end

// Textures for reflectors

#declare nrf = 5; // Max number of reflectors

// Reflector colors:
#declare clr_rf = array[nrf]
#declare clr_rf[0] = < 1.000, 0.300, 0.000 >;
#declare clr_rf[1] = < 0.700, 0.900, 0.000 >;
#declare clr_rf[2] = < 0.000, 0.800, 0.800 >;
#declare clr_rf[3] = < 0.300, 0.000, 1.000 >;
#declare clr_rf[4] = < 0.700, 0.000, 0.500 >;

// Reflector textures:
#declare tx_rf_v = array[nrf]
#declare tx_rf_e = array[nrf]
#declare tx_rf_f = array[nrf]
#local i = 0;
#while (i < nrf)
  #declare tx_rf_v[i] = txc( clr_rf[i], 0.00 )
  #declare tx_rf_e[i] = txc( clr_rf[i], 0.00 )
  #declare tx_rf_f[i] = txc( clr_rf[i], 0.80 )
  #local i = i + 1;
#end

// Mesh textures:
#declare clr_mesh = < 0.300, 0.300, 1.000 >;
#declare tx_mesh_v = txc( clr_mesh, 0.00 )
#declare tx_mesh_e = txc( clr_mesh, 0.00 )
#declare tx_mesh_f = txc( clr_mesh, 0.90 )
