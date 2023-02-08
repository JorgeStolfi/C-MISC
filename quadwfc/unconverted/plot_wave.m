function res = plot_wave(all_edges)
%global all_edges;

l = length(all_edges);
res = -1*ones(l);
for i=1:l
   edge_a = all_edges(i);
   if(edge_a.destroyed==0)
      res(i) = a_plot_edge([i 1],all_edges,'r');
   end;
end;
