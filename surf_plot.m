## Copyright (C) 2017 Sunder
## Author: Sunder <sunder@sai-Vostro-3546>
## Created: 2017-11-09

function  surf_plot (plot)
 M = load("plot.txt");
 x = M(:,1);
 y = M(:,2);
 phi = M(:,3);

 
xlin=linspace(min(x),max(x),200);
ylin=linspace(min(y),max(y),200);
[X,Y]=meshgrid(xlin,ylin);
Z=griddata(x,y,phi,X,Y,'linear');
mesh(X,Y,Z);
axis tight; 
grid off
hold on

surf(X,Y,Z)
colormap jet
shading interp
colorbar
endfunction
 
