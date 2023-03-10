// ======================================================================
// C?MERA 

camera {
  location  <  34.00, 8.00, 5.00 >  // Posi??o do observador.
  right     -1.0*x                 // Largura RELATIVA da imagem.
  up        0.75*y                 // Altura RELATIVA da imagem.      
  sky       z                      // Qual dire??o ? "para cima"?
  look_at   <  0.00, 0.00, 0.00 >  // Para onde a c?mera est? apontando.
} 
// Nota: os par?metros "right" e "up" devem ter a mesma propor??o
// que os par?metros ${WIDTH} e ${HEIGHT} no Makefile.

// ======================================================================
// DESCRI??O DA CENA 

background{ color rgb < 1.0, 1.0, 1.0 > }

background{ color rgb < 1.0, 1.0, 1.0 > }

#declare vidro = 
  texture {
    pigment { color rgb 0.9 * < 1.0, 1.0, 1.0 > filter 0.3}
    finish { diffuse 0.0 specular 0.0 reflection 0.0 roughness 0.0 ambient 1.0 }
  }


#declare preto = 
  texture {
    pigment { color rgb < 0.10, 0.10, 0.10 > }
    finish { diffuse 0.0 specular 0.0 roughness 0.0 ambient 0.1 }
  }


#declare branco = 
  texture {
    pigment { color rgb < 0.9, 0.9, 0.9 > }
    finish { diffuse 0.5 specular 0.5 roughness 0.05 ambient 0.1 }
  }


#macro Make_Vertex(v1)
 sphere {v1 0.08 texture{preto}}
#end

//#macro Make_Face(v1, v2, v3)
//union{
//   triangle { v1 v2 v3 texture{vidro}}
//   cylinder { v1 v2 0.1 texture {branco}}
//   cylinder { v2 v3 0.1 texture {branco}}
//   cylinder { v1 v3 0.1 texture {branco}}
// }
//#end
  
#macro Make_Face(numPoints, points)
 union{
  polygon { numPoints+1
    #declare i=0;
    #while (i < numPoints)
      points[i]
      #declare i=i+1;
    #end
    points[0]
    texture{vidro}
  }
  #declare i=0;
  #while (i < numPoints-1)
    cylinder { points[i] points[i+1] 0.04 texture {preto}}
    #declare i=i+1;
  #end
  cylinder { points[numPoints-1] points[0] 0.04 texture {preto}}
 }
#end

                                                                
//object{Make_Face(4, array[4] {<0.0 0.0 0.0>,<10.0 0.0 0.0>,<10.0 10.0 0.0>,<0.0 10.0 0.0>})}  
