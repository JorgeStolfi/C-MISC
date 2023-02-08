function new_wf = wave_prop(all_edges,vel,time,reflectors,op)
% Propagates the wavefront represented by "all_edges" in a medium 
% with velocity "vel" m/s durin "time" seconds

%global reflectors;
tri=[];

k = 50;

l = length(all_edges);

data_out = cell(1,3);

interc=0;

for i=1:l
%   if((all_edges(i).destroyed)==0)
   
   stop=0;   
   for j=1:2:3
      vectors = all_edges(i).data{j};
      
      vel = vectors.vel;
      
   d_vector = [vectors.x - vectors.px ,... %difference vector
         vectors.y - vectors.py ,...
         vectors.z - vectors.pz];
   
   d_norm_vector = d_vector / norm(d_vector); %difference vector normalized

   %pause;
   inc_vector = d_norm_vector * vel; % increment vector
   
   data_out{j}.px = vectors.x;
   data_out{j}.py = vectors.y;
   data_out{j}.pz = vectors.z;
   data_out{j}.pvel = vel;
   
   
   data_out{j}.x = vectors.x + inc_vector(1);
   data_out{j}.y = vectors.y + inc_vector(2);
   data_out{j}.z = vectors.z + inc_vector(3);
   data_out{j}.vel = vel;

   
   tri = intersec_1(data_out{j},reflectors(1),op);
   if(length(tri)>0)
      stop=1;
      data_out{j}= tri;
      interc=interc+1;
   end;
   
	end;
   new_wf(i) = all_edges(i);
   %if(~stop)
      new_wf(i).data = data_out;
   %end;
   
%   end;
end;
%interc
   
