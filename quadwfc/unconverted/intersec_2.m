function res = intersec_2(seg,triangule)

err=10;

%if(seg.z<triangule(1,3))
%   res =1;
%   return;
%end;


w=1;%puttin in homogeneous coord system

A = [w,triangule(1,:)];
B = [w,triangule(2,:)];
C = [w,triangule(3,:)];

D = [w,seg.px,seg.py,seg.pz];
E = [w,seg.x,seg.y,seg.z];

d1 = det([A;B;C;D]);
d2 = det([A;B;C;E]);
cond1 = (((d1<=err)&(d2>=-1*err))|...
   ((d2<=err)&(d1>=-1*err)));
if(~cond1)
   res = 0;
   return;
end;

d3 = det([D;E;A;B]);
d4 = det([D;E;A;C]);
cond2 = (((d3<=err)&(d4>=-1*err))|...
   ((d4<=err)&(d3>=-1*err)));
if(~cond2)
   res = 0;
   return;
end;


d5 = det([D;E;B;C]);
d6 = det([D;E;B;A]);
cond3 = (((d5<=err)&(d6>=-1*err))|...
   ((d6<=err)&(d5>=-1*err)));
if(~cond3)
   res = 0;
   return;
end;

d7 = det([D;E;C;A]);
d8 = det([D;E;C;B]);
cond4 = (((d7<=err)&(d8>=-1*err))|...
   ((d8<=err)&(d7>=-1*err)));
if(~cond4)
   res = 0;
   return;
end;

res = 1;
