function res = unplot_wave(wave)

l = length(wave);
for i=1:l
   if(wave(i)~=-1)
      delete(wave(i));
   end;
end;
