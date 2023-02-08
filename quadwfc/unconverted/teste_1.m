inivel = 50;

global all_edges;
a = all_edges;

figure(gcf);
hold on;

for i=1:length(a)
   for j=1:2:3
      all_edges(i).data{j}.pvel = inivel;
      all_edges(i).data{j}.vel = inivel;

   end;
end;


reflectors(1).velbelow = 55;
reflectors(1).velabove = 50;

for i=1:1000000
a = wave_prop(a,1000,100,reflectors,0);
wave = plot_wave(a);
if(mod(i,5)==0)
     print('-depsc2',['onda_' num2str(i)]);
end;
pause(0.00001);
%pause;

unplot_wave(wave);
%clf;
%my_mesh;

end;
