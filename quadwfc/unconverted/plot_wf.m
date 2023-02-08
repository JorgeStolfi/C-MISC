function [res] = plot_wf(wf)

plot3(wf.rMd(:,1), wf.rMd(:,2), wf.rMd(:,3),'b.');

hold on;
%axis equal;
