function [m,gamma,alpha,delta,D1,D2,x_wp,y_wp,Vxw,Vyw,conv] = fluid_solver(x,y,...
    midx,midy,delta,alpha,...
    gamma,gamma_w,vortx,vorty,n1,t1,vbn,vbt)
%Solve for the fluid conditions
global k
global L
global V_inf
global interval

gammanew = gamma(k);
deltanew = delta(k);
alphanew = alpha(k);
conv = 0;

gamma_tol = 1*10^-3;    gamma_err = 1;
alpha_tol = 1*10^-3;    alpha_err = 1;
delta_tol = 1*10^-3;    delta_err = 1;

%% Affect of panels on each other
[An,At,XX,XX,Bnv,Btv,XX,XX] = Panel_Influence(x,y,x,y);

Bn = Bnv*ones(length(midx),1);
Bt = Btv*ones(length(midx),1);

%% Affect from free vorticies
if k > 2

    pan = 0;
    [XX,XX,Ct,Cn] = Free_Vort_Affect(vortx,vorty,[midx midy],gamma_w,n1,t1,pan);
else
    Cn = zeros(length(midx),1);
    Ct = Cn;
end

count = 0;
while abs(gamma_err)>gamma_tol || abs(delta_err)>delta_tol...
        || abs(alpha_err)>alpha_tol  
    gamma(k) = gammanew;
    alpha(k) = alphanew;
    delta(k) = deltanew;
    
    
    x_wp = x(end)+delta(k)*cos(alpha(k));
    y_wp = y(end)-delta(k)*sin(alpha(k));
    
    %% Affect from wake panel
    [XX,XX,XX,XX,Wn,Wt,XX,XX] = Panel_Influence(x,y,[x(end),x_wp],[y(end),y_wp]);
    
    %% Normal Velocity
    var = An\((Wn.*L)/delta(k) - Bn);
    const = An\(-Wn.*L*(gamma(k-1))/delta(k) - Cn - n1*V_inf' + vbn);
    
    %% Tangent Velocity
    D1 = At*var + Bt - Wt*L/delta(k);
    D2 = At*const + (Wt*gamma(k-1)*L)/delta(k) + Ct + t1*V_inf';
    
    a = D1(1)^2 - D1(end)^2;
    b = 2*(D1(1)*D2(1) - D1(end)*D2(end) - L/interval);
    c = D2(1)^2-D2(end)^2 + 2*gamma(k-1)*L/interval + vbn(1)^2 - vbn(end)^2;
    
    %% Solve for gamma
    [gammanew,N_Conv,gold] = Newton_It(a,b,c,D1,D2);
    
    %% Source Strength
    m = var*gammanew + const;
    
    %% Recompute velo at trailing edge panels
    if k > 2
        pan = 1;
        [Cxw,Cyw,XX,XX] = Free_Vort_Affect(vortx,vorty,...
            [(x_wp+x(end))/2 (y_wp+y(end))/2],gamma_w,n1,t1,pan);
    else 
        Cxw = 0;
        Cyw = 0;
    end
    
    [XX,XX,Axw,Ayw,XX,XX,Bxw,Byw] = Panel_Influence([x(end),x_wp],[y(end),y_wp],x,y);
    
    Vxw = Axw*m + sum(Bxw)*gammanew + Cxw + V_inf(1);
    Vyw = Ayw*m + sum(Byw)*gammanew + Cyw + V_inf(2);
    
    alphanew = -atan2(Vyw,Vxw);
    deltanew = sqrt(Vyw^2 + Vxw^2)*interval;

    gamma_err = (gamma(k) - gammanew)/gamma(k);
    delta_err = (delta(k) - deltanew)/delta(k);
    alpha_err = (alpha(k) - alphanew)/alpha(k);
   
    count = count + 1;
    cerr(count,:) = [gamma_err,delta_err,alpha_err];

    if count > 100
        dif_err=zeros(count-1,3);
        for h=1:3
            dif_err(:,h) = cerr(2:count,h)-cerr(1:count-1,h);
        end
    end    
    
    prev_alpha(count,1) = alphanew;
    prev_delta(count,1) = deltanew;
    prev_gamma(count,1) = gammanew;
    prev_gammaold(count,1) = gold;
   

    if count > 200
        disp('Wake parameters not converging')
        alpha(k) = mean(prev_alpha(100:end));
        delta(k) = mean(prev_delta(100:end));
        gamma(k) = mean(prev_gamma(100:end));
        conv = 1;
        break
    end
end
    
    
    






end

