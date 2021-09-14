function [Vxp,Vyp] = point_velo(m,gamma,x,y,pointsx,pointsy,...
    vortx,vorty,gamma_w,x_wp,y_wp,delta,Vb)
%Velocity of the points of integration
global k
global L
global V_inf

%% From body panels
[XX,XX,Ax,Ay,XX,XX,Bx2,By2] = Panel_Influence(pointsx,pointsy,x,y);
for i=1:length(pointsx)-1
    Bx(i,1) = sum(Bx2(i,:));
    By(i,1) = sum(By2(i,:));
end


%% From wake panel
[XX,XX,XX,XX,XX,XX,Bxwp,Bywp] = Panel_Influence(pointsx,pointsy,[x(end) x_wp]...
    ,[y(end) y_wp]);

%% From free vortices
for i=1:1:length(pointsx)-1
    mpx(i,1) = (pointsx(i+1)+pointsx(i))/2;
    mpy(i,1) = (pointsy(i+1)+pointsy(i))/2;
end
    
if k > 3
    n1 = 0;
    t1 = 0;
    pan = 1;
    [Cx,Cy,XX,XX] = Free_Vort_Affect(vortx,vorty,[mpx mpy],gamma_w,n1,t1,pan);
else
    Cx = 0;
    Cy = Cx;
end

%% Sum Velocities
Vxp = Ax*m + Bx.*gamma(k) + Bxwp.*(gamma(k-1)-gamma(k))*L/delta(k) + Cx - Vb(1) + V_inf(1);
Vyp = Ay*m + By.*gamma(k) + Bywp.*(gamma(k-1)-gamma(k))*L/delta(k) + Cy - Vb(2) + V_inf(2);



end

