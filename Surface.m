clear;close;clc;
figure;
hold;

[X,Y] = meshgrid(0:1:100*pi,0:1:100*pi);
Z = 400*sin(X./100) + 400*(cos(Y./100-pi)+1);
%contour(X,Y,Z)
surface(X,Y,Z,'EdgeColor','none','FaceAlpha',0.5);

[A,B] = meshgrid(0:20*pi:100*pi,0:20*pi:100*pi);
xypoints = [A(:), B(:)];
A = xypoints(:,1);
B = xypoints(:,2);
C = 400*sin(A./100) + 400*(cos(B./100-pi)+1)+50;
surfPoints = [A,B,C];
scatter3(A, B, C);

k=1:length(A); 
text(A,B,C,num2str(k'))



view(3)