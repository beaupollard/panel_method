function [phi] = Potential2(m,gamma,midx2,midy2,vortx2,vorty2,del,n,x2,y2...
    ,gamma_w,x_wp2,y_wp2,delta,Vt,theta,k,pos,Vb)

%% Update on 6/24/18------ Rotated the phi into body fixed frame

%Calculate the potential function
x3 = x2;
y3 = y2;
% global point_upstream
midx = midx2 - pos(1,1);
midy = midy2 - pos(1,2);
vortx = vortx2 - pos(1,1);
vorty = vorty2 - pos(1,2);
x = x2 - pos(1,1);
y = y2 - pos(1,2);
x_wp = x_wp2 - pos(1,1);
y_wp = y_wp2 - pos(1,2);
midpoint = round(length(midx)/2);

LEx = midx(midpoint) + n(midpoint,1)*0.00001;
LEy = midy(midpoint) + n(midpoint,2)*0.00001;
%% New phi
L1 = abs(200);
L2 = abs(7);

y2 = LEy-L1*sin(theta);
y1 = y2+L2*cos(theta);
x2 = LEx-L1*cos(theta);
x1 = x2-L2*sin(theta);

Nd = 40;
Na = 200;

% down
yd = linspace(y1,y2,Nd+1)';
xd = linspace(x1,x2,Nd+1)';
for i=1:Nd
    del1(i,1) = (yd(i+1)-yd(i))/2;
end
    
xpts1 = (xd(1:end-1)+xd(2:end))/2;
ypts1 = (yd(1:end-1)+yd(2:end))/2;


% across
r = abs(L1);
ya = linspace(y2,LEy,Na+1)';
xa = linspace(x2,LEx,Na+1)';

xpts2 = (xa(1:end-1)+xa(2:end))/2;
ypts2 = (ya(1:end-1)+ya(2:end))/2;
if k == 44
end
%% Convert to BFF
for i=1:length(xd)
    xd2(i,1) = [cos(theta) sin(theta)]*[xd(i);yd(i)];
    yd2(i,1) = [-sin(theta) cos(theta)]*[xd(i);yd(i)];
   
end
for i=1:length(x)
    x2(i,1) = [cos(theta) sin(theta)]*[x(i);y(i)];
    y2(i,1) = [-sin(theta) cos(theta)]*[x(i);y(i)];    
end
for i=1:length(xa)
    xa2(i,1) = [cos(theta) sin(theta)]*[xa(i);ya(i)];
    ya2(i,1) = [-sin(theta) cos(theta)]*[xa(i);ya(i)];    
end
for i=1:length(vortx)
    vortx2(i,1) = [cos(theta) sin(theta)]*[vortx(i);vorty(i)];
    vorty2(i,1) = [-sin(theta) cos(theta)]*[vortx(i);vorty(i)];    
end

xpt = [xd2;xa2];
ypt = [yd2;ya2];
x_wp2 = [cos(theta) sin(theta)]*[x_wp;y_wp];
y_wp2 = [-sin(theta) cos(theta)]*[x_wp;y_wp];

xpt = xpt+pos(1,1);
ypt = ypt+pos(1,2);
x2 = x2 + pos(1,1);
y2 = y2 + pos(1,2);
x_wp2 = x_wp2 + pos(1,1);
y_wp2 = y_wp2 + pos(1,2);
vortx = vortx2 + pos(1,1);
vorty = vorty2 + pos(1,2);

[Vxp,Vyp] = point_velo(m,gamma,x2,y2,xpt,ypt,...
    vortx,vorty,gamma_w,x_wp2,y_wp2,delta,Vb);

phi_inf = 0;

for i=1:Nd
    phi_inf = phi_inf + Vyp(i)*del1(i);
end
count1 = 1;
for i=Nd+1:Nd+Na
    phi_inf = phi_inf + Vxp(i)*del2(count1);
    count1 = count1 + 1;
end

phin(midpoint,1) = phi_inf;

for i=midpoint-1:-1:1
    phin(i,1) = phin(i+1) - Vt(i)*del(i);
end
for i=midpoint+1:length(x)
    phin(i,1) = phin(i-1,1) + Vt(i-1)*del(i-1);
end

phi = (phin(1:i-1)+phin(2:i))/2;

    



end

