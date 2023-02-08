function plot_edge(e,cor)

o = Org(e);
d = Dest(e);

plot3([o.x,d.x],[o.y,d.y],[o.z,d.z],cor);
%plot([o.x,d.x],[o.y,d.y],cor);

