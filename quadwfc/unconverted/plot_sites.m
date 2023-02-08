function plot_sites(sites)

x=[];
y=[];
z=[];

for i=1:length(sites)
   x = [x, sites{i}.x];
   y = [y,sites{i}.y];
   z = [z,sites{i}.z];

end;

plot3(x,y,z,'gp');
%plot(x,y,'gp');

hold on;