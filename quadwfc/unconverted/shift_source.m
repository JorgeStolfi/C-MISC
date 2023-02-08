source = [5000,5000,4000];
global all_edges;


for i=1:length(all_edges)
   v = all_edges(i).data{1};
   v.px = source(1);
   v.py = source(2);
   v.pz = source(3);
   all_edges(i).data{1} = v;
   v = all_edges(i).data{3};
   v.px = source(1);
   v.py = source(2);
   v.pz = source(3);
   all_edges(i).data{3} = v;
   
end;

