function res = intersec_1(seg,reflec,op)
%inc =0

err = 0; %sobreposicao 

res=[];

xmin = min([seg.px seg.x]);
xmax = max([seg.px seg.x]);

ymin = min([seg.py seg.y]);
ymax = max([seg.py seg.y]);

zmin = min([seg.pz seg.z]);
zmax = max([seg.pz seg.z]);

if((zmin>max(reflec.zz))|(zmax<min(reflec.zz))...
      |(xmin>max(reflec.xx))|(xmax<min(reflec.xx))...
      |(ymin>max(reflec.yy))|(ymax<min(reflec.yy)))
   return;
end;


a_lim_inf_x = find(reflec.xx>=xmin);
a_lim_sup_x = find(reflec.xx>xmax);
lim_inf_x = min(a_lim_inf_x)-1-1;
lim_sup_x = min(a_lim_inf_x)+1;

a_lim_inf_y = find(reflec.yy>=ymin);
a_lim_sup_y = find(reflec.yy>ymax);
lim_inf_y = min(a_lim_inf_y)-1-1;
lim_sup_y = min(a_lim_sup_y)+1;

dx = lim_sup_x- lim_inf_x;
dy = lim_sup_y- lim_inf_y;


zz = reflec.zz(lim_inf_y+1:lim_sup_y-1,lim_inf_x+1:lim_sup_x-1);
xx = reflec.xx(lim_inf_x+1:lim_sup_x-1);
yy = reflec.yy(lim_inf_y+1:lim_sup_y-1);

if((zmin>max(zz))|(zmax<min(zz)))
     return;
end;



% -------Collectin possible triangules-----------
triangules = cell(1,2*abs(dx*dy));
trs = 1;

for i=lim_inf_x:lim_sup_x-1
   for j=lim_inf_y:lim_sup_y-1
      s = [reflec.xx(i),reflec.yy(j),reflec.zz(j,i)];
      s = [s;reflec.xx(i),reflec.yy(j+1),reflec.zz(j+1,i)];
      s = [s;reflec.xx(i+1),reflec.yy(j),reflec.zz(j,i+1)];
      triangules{trs}=s;
      trs = trs+1;
      s = [reflec.xx(i),reflec.yy(j+1),reflec.zz(j+1,i)];
      s = [s;reflec.xx(i+1),reflec.yy(j+1),reflec.zz(j+1,i+1)];
      s = [s;reflec.xx(i+1),reflec.yy(j),reflec.zz(j,i+1)];
      triangules{trs}=s;
      trs = trs+1;
   end;
end;

triang_num = length(triangules);
aux = cell(1,triang_num);
k=1;

%triangulos = [];
%for g=1:triang_num
%    s2 = triangules{g};
%    s2 = [s2;s2(1,:)];
%    triangulos = [triangulos,plot3(s2(:,1),s2(:,2),s2(:,3),'r*')];
%end;

for i=1:length(triangules)
   aux_2 = triangules{i};
   if(intersec_2(seg,aux_2))
      if(seg.z<seg.pz)
         seg.vel = reflec.velbelow;
      else
         seg.vel = reflec.velabove;
      end;
      
      aux{k} = triangules{i};
      [r,t] = snell(seg,aux_2);
      if(op==0)
         res = r; %reflect
      else
         res = t; %transmit
      end;
      
      return;
 %     inc = inc+1
      k=k+1;
    end;
end;

if(k<=1)
   return;
end;

res = cell(1,k-1);
for j=1:k-1;
   res{j}=aux{j};
end;
