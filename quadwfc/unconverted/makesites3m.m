function [sites,n] = makesites3(op)

wf = init_wf(19);
sites = cell(1,length(wf.dest));

for i=1:length(wf.dest);
   sites{i}.x = wf.dest{i}(1);
   sites{i}.y = wf.dest{i}(2);
   sites{i}.z = wf.dest{i}(3);
end;
