function [wf] = init_wf(ap)
ap = 180;

ap = ap*pi/180;  %converting to rad

% input parameters
d_ap = 0.3;
d = 0.3;
l = 500;



% First point
h=0;
rM = [];
j=0;
i=0;
r = l*[ cos(j)*sin(i) , sin(j)*sin(i) , -cos(i) ];
h = h+1;
wf.orig{h} = [5000, 5000, 4000];
wf.dest{h} = wf.orig{h}+r;
rM = [rM;r];


%other points
disp(ap/2);
for i=d_ap:d_ap:ap/2
	for j=0:d:2*pi
		r = l*[ cos(j)*sin(i) , sin(j)*sin(i) , -cos(i) ];
		
		h = h+1;
		wf.orig{h} = [5000, 5000, 4000];
		wf.dest{h} = wf.orig{h}+r;
		
		rM = [rM;r];
	end;
end;

wf.rMo = zeros(size(rM));
wf.rMd = rM;
wf.raynum = h;


%plot3(rM(:,1) , rM(:,2) , rM(:,3),'b.');
%axis equal;
