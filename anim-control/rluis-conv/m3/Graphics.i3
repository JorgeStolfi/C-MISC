  Texture = RECORD
      aR, aG, aB,  (* ambient color *)
      dR, dG, dB,  (* diffuse color *)
      sR, sG, sB,  (* specular color *)
      tR, tG, tB,  (* transparent color *)
      ir,          (* index of refraction *)
      n: REAL;     (* Phong specular exponent *)
    END;
