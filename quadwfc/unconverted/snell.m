function [r,t] = snell(seg,triangule)

res = intersec_3(seg,triangule);

%settin the reflected ray
r.pvel = seg.pvel;
r.vel = seg.pvel;

r.px = res.px;
r.py = res.py;
r.pz = res.pz;

r.x = res.rx;
r.y = res.ry;
r.z = res.rz;


%settin the transmitted ray
t.pvel = seg.pvel;
t.vel = seg.vel;

t.px = res.px;
t.py = res.py;
t.pz = res.pz;

t.x = res.tx;
t.y = res.ty;
t.z = res.tz;