function res = intersec_3(seg,triangule)

vel1 = seg.pvel;
vel2 = seg.vel;

o = triangule(1,:);

v1 = triangule(2,:) - o;
v2 = triangule(3,:) - o;

vn = cross(v2,v1); %normal vector
if(vn(3)<0)
   disp(vn(3));
end;

v = [seg.x-o(1),seg.y-o(2),seg.z-o(3)];
c1 = dot(v,vn);
c2 = dot(vn,vn);
c3 = c1/c2;
projvn = (c3)*vn;

% eliminatin numerical trash
projvn = 10*projvn;
projvn = round(projvn);
projvn = projvn/10;

vref = v -2*projvn;%q'

res.rx = vref(1) + o(1);
res.ry = vref(2) + o(2);
res.rz = vref(3) + o(3);


pq = [seg.x - seg.px,...
      seg.y - seg.py,...
      seg.z - seg.pz];

p = [seg.px,seg.py,seg.pz];


d = -1*(dot(vn,o));

lambda = (-1*d - dot(p,vn))/(dot(pq,vn));

u = p + lambda*pq; %u

% verifying obtention of the intersection
qqq = u + (1-lambda)*pq;
qqqqq = [seg.x,seg.y,seg.z];


qprime = vref + o;
uq = qprime - u; %uq'
uq = 10*uq;
uq = round(uq);
uq = uq/10;

%incidence angle
cosiaux = dot(uq,vn)/(norm(uq)*norm(vn));
cosi = sqrt(1- cosiaux^2);

cost = (vel2/vel1)*cosi;%cosine of transmission


uqprojvn = (dot(uq,vn)/dot(vn,vn))*vn;
%eliminatin numerical trash
uqprojvn = 10*(uqprojvn);
uqprojvn = round(uqprojvn);
uqprojvn = 0.1*(uqprojvn);

uq2primes = uq-uqprojvn; %uq''

% Total reflection not allowed!!
if(cost<-1)
   cost=-1;
end;
if(cost>1)
   cost=1;
end;

if(cost<0)
   disp('obtuso');
end;

if(cost~=0)
   tan2 =(1/cost)^2 - 1;
   if(tan2<0)
      disp('tan2 negativo');
   end;
   
   normx = norm(uq2primes)*sqrt(tan2);%|x| = |uq''|*tan(t)
   x = -1*(normx/norm(vn))*vn;% vector x
   
   x = 10*x;
   x = round(x);
   x = 0.1*x;

   uqtv = (uq2primes + x)/(norm(uq2primes + x));%unit vector in uqt direction
   if(tan2 == 0)
       disp(tan2);
   end;
   
else
   uqtv = -1*vn/norm(vn);
end;

uqt = (norm(uq)*vel2/vel1)*uqtv; %uqt
uqt = 10*uqt;
uqt = round(uqt);
uqt = 0.1*uqt;

t = u + uqt;

res.tx = t(1);
res.ty = t(2);
res.tz = t(3);

res.px = u(1);
res.py = u(2);
res.pz = u(3);
