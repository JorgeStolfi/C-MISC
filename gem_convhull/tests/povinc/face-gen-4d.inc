// Last edited on 2014-07-22 12:38:50 by stolfilocal
// ======================================================================
// CÂMERA 

camera {
  location  <  40.00, 5.00, 5.00 >  // Posição do observador.
  right     -1.0*x                 // Largura RELATIVA da imagem.
  up        0.75*y                 // Altura RELATIVA da imagem.      
  sky       z                      // Qual direção é "para cima"?
  look_at   <  0.00, 0.00, 0.00 >  // Para onde a câmera está apontando.
} 
// Nota: os parâmetros "right" e "up" devem ter a mesma proporção
// que os parâmetros ${WIDTH} e ${HEIGHT} no Makefile.

// ======================================================================
// FONTES DE LUZ

light_source {
  < 40.0, 0.0, 0.0 >              // Posição da lâmpada.
  color rgb 1 * < 1.00, 1.00, 1.00 >   // Intensidade e corda luz.
} 

light_source {
  < 00.0, 40.0, 0.0 >              // Posição da lâmpada.
  color rgb 1 * < 1.00, 1.00, 1.00 >   // Intensidade e corda luz.
} 

light_source {
  < 00.0, -40.0, 0.0 >              // Posição da lâmpada.
  color rgb 1 * < 1.00, 1.00, 1.00 >   // Intensidade e corda luz.
} 

light_source {
  < 0.0, 0.0, 40.0 >              // Posição da lâmpada.
  color rgb 2 * < 1.00, 1.00, 1.00 >   // Intensidade e corda luz.
} 


// ======================================================================
// DESCRIÇÃO DA CENA 

//background{ color rgb < 0.8, 0.95, 0.95 > }
background{ color rgb < 1.0, 1.0, 1.0 > }

#declare vidro = 
  texture {
    pigment { color rgb < 0.9, 0.9, 0.9 > filter 0.80}
    finish { diffuse 0.3 specular 0.1 reflection 0.0 roughness 0.01 ambient 0.3 }
  }


#declare preto = 
  texture {
    pigment { color rgb < 0.10, 0.10, 0.50 > }
    finish { diffuse 0.5 specular 0.5 roughness 0.05 ambient 0.1 }
  }


#declare branco = 
  texture {
    pigment { color rgb < 0.7, 0.5, 0.5 > }
    finish { diffuse 0.5 specular 0.1 roughness 0.05 ambient 0.1 }
  }


#macro Make_Vertex(v1)
 sphere {v1 0.3 texture{preto}}
#end

#macro Make_Face(numPoints, points)
 union{
  polygon { numPoints
//  prism { 0, 1, numPoints

    #declare i=0;
    #while (i < numPoints)
      points[i]
      #declare i=i+1;
    #end
    texture{vidro}
  }
  #declare i=0;
  #while (i < numPoints-1)
    cylinder { points[i] points[i+1] 0.1 texture {branco}}
    #declare i=i+1;
  #end
  cylinder { points[numPoints-1] points[0] 0.1 texture {branco}}
 }
#end

                                                                
//object{Make_Face(4, array[4] {<0.0 0.0 0.0>,<10.0 0.0 0.0>,<10.0 10.0 0.0>,<0.0 10.0 0.0>})}  
