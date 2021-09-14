%
clear
clc

global k
global interval
global mass

Initialize;

[mass,Cx,Cy] = Centroid(xi,yi); % Determine the centroid location of the body
[I] = Inertia(xi,yi,Cx,Cy);     % Determine the moment of inertia of the body

t_length = 20000;

load('Mb.mat')                  % Load in the torque that will be applied
flagkk=1;

count_the=1;
counterp = 1;

for k=2:t_length
    Fx = [];
    Fy = [];
    time(k) = time(k-1)+interval;
    text = sprintf('Time = %g',time(k)); 
    disp(text);
    count = 1;
    
    %% Get the initial position and parameters
    [pos,vxp,vyp,vwp,alpha,delta,gamma,midx1,midy1] = Initial_Guess(pos,Vb,alpha, ...
        delta,gamma,x(:,k-1),y(:,k-1),xi,yi,time);
    
    %% Initialize the force/moment errors and vortex position 
    Fx_err = 0;
    Mf_err = 0;
    conv2 = 0;
    if k > 2 
        vortx = vort_posx(1:(k-2),k-1);
        vorty = vort_posy(1:(k-2),k-1);
    else
        vortx = 0;
        vorty = 0;
    end    

    %% Iterates between fluid solver and force solver until convergence
    while Fx_err == 0 || count < 4  || Fy_err == 0 || Mf_err == 0
        
        % Update body position
        [n1,t1,del,midx,midy,zx,zy] = body_orientation(pos,xi,yi);
        x(:,k) = zx;
        y(:,k) = zy;
        
        % Determine body velocities for fluid solver
        Vbx(:,k) = (midx-midx1)/interval;
        Vby(:,k) = (midy-midy1)/interval;

        % Convert velocities to normal and tangent velocities
        Vbn(:,k) = Vbx(:,k).*n1(:,1) + Vby(:,k).*n1(:,2);
        Vbt(:,k) = Vbx(:,k).*t1(:,1) + Vby(:,k).*t1(:,2);
        
        % Satisfy the no flow through BC to determine vortex panel
        % strengths
        [m,gamma,alpha,delta,D1,D2,x_wp,y_wp,Vxw,Vyw,conv2] = fluid_solver(zx,zy,midx,midy...
            ,delta,alpha,gamma,gamma_w,vortx,vorty,n1,t1,Vbn(:,k),Vbt(:,k));
        if conv2 == 1
            conv(k,1) = conv(k,1)+1;
        end
        Vt(:,k) = D1*gamma(k) + D2;
        VW(k,:) = [Vxw,Vyw];
        
        % Use updated vortex strengths to determine body forces
        [phi,cp,Pres] = Determine_Forces(Vt(:,k),Vbn(:,k),n1,t1,pos,m,gamma,midx...
                ,midy,vortx,vorty,del,x_wp,y_wp,delta,phi,Vbx(:,k),Vby(:,k),...
                x(:,k),y(:,k),gamma_w,Vbt(:,k),Vb(k,:));
        weight = 0.75;
        
        % Calculate forces from pressure coefficients
        Fx(count,1) = sum(-cp.*del.*n1(:,1));
        Fy(count,1) = sum(-cp.*del.*n1(:,2));
        
        % Calculate the moment acting on the body
        [M_f(count,1)] = Moment(cp,del,n1,midx1,midy1,pos(k,:));
        
        % Calculate the error between the current and previous iteration
        if count > 1
            Mf_err = abs(M_f(count-1)-M_f(count))<M_f_tol;
            Fx_err = abs(Fx(count-1)-Fx(count))<Fx_tol;
            Fy_err = abs(Fy(count-1)-Fy(count))<Fy_tol;
            Fxw = Fx(count)*weight + (1-weight)*Fx(count-1);
            Fyw = Fy(count)*weight + (1-weight)*Fy(count-1);
            Mfw = M_f(count)*weight + (1-weight)*M_f(count-1);
        else
            Fxw = Fx(count)*weight + (1-weight)*cd(k-1,1);
            Fyw = Fy(count)*weight + (1-weight)*cl(k-1,1);
            Mfw = M_f(count);
        end
        
        % Average the forces/moments if it is taking to long to converge
        % (this really only helps get past the first few iterations)
        if count > 30
            Mfw = mean(M_f(20:30));
            Fxw = mean(Fx(20:30));
            Fyw = mean(Fy(20:30));
            Fx_err = 0; Fy_err = 0; Mf_err = 0;
            dragx = -sign(-V_inf(1)+vxp)*5/2*0.0146*2*L*(-V_inf(1)+vxp)^2;
            dragy = -sign(-V_inf(2)+vyp)*5/2*0.0146*2*L*(-V_inf(2)+vyp)^2;
            vxp = Vb(k-1,1) + (Fxw+dragx)/1*interval;
            vyp = Vb(k-1,2) + (Fyw+dragy)/1*interval;            
            omegap = weight*(Vb(k-1,3)) + (1-weight)*(Mfw+Mb(k,1))/(2*I(1,1))*interval;
            pos(k,1) = pos(k-1,1) + vxp*interval;
            pos(k,2) = pos(k-1,2) + vyp*interval;
            pos(k,3) = pos(k-1,3) + omegap*interval;
            break
        end
        
        % Update the velocities/new body positions
        dragx = -sign(-V_inf(1)+vxp)*5/2*0.0146*2*L*(-V_inf(1)+vxp)^2;
        dragy = -sign(-V_inf(2)+vyp)*5/2*0.0146*2*L*(-V_inf(2)+vyp)^2;
        accx(k,1)=(Fxw+dragx);
        accy(k,1)=(Fyw+dragy);
        omegad(k,1)=(Mfw+Mb(k,1))/(2*I(1,1));
        vxp = Vb(k-1,1) + (Fxw+dragx)/1*interval;
        vyp = Vb(k-1,2) + (Fyw+dragy)/1*interval;
        omegap = weight*(Vb(k-1,3)) + (1-weight)*(Mfw+Mb(k,1))/(2*I(1,1))*interval;
        pos(k,1) = pos(k-1,1) + vxp*interval;
        pos(k,2) = pos(k-1,2) + vyp*interval;
        pos(k,3) = pos(k-1,3) + omegap*interval;
        
        % Print the body velocities and angular orientation
        text2 = sprintf('Vx = %g \t Vy = %g \t Pos = %g',vxp,vyp,pos(k,3).*180/pi);
        disp(text2);  
        count = count + 1;
    end

    %% Store the useful information from the current iteration
    Fxp(:,k) = sum(-cp.*del.*n1(:,1));
    Fyp(:,k) = sum(-cp.*del.*n1(:,2));
    Mf_p(k,1) = Mfw;
    Vb(k,1) = vxp;
    Vb(k,2) = vyp;
    Vb(k,3) = omegap;
    cp2(:,k) = cp;
    Pressure(:,k) = Pres;
    gamma_w(k,1) = L*(gamma(k-1)-gamma(k));
    
    %% Shed Vorticies
    [vortx2,vorty2] = Vortex_Shed2(vortx,vorty,gamma_w ...
    ,gamma,m,x_wp,y_wp,Vxw,Vyw,x(:,k),y(:,k),delta,pos(k,3),n1,alpha(k)); 

    % This makes it easier to call the vortices in the next iteration but
    % it is a pretty bad way of doing it
    vort_posx(1:length(vortx2),k) = vortx2;
    vort_posy(1:length(vortx2),k) = vorty2;
    
    % Store the CL, CD, and wake panel location 
    cl(k,1) = sum(-cp.*del.*n1(:,2));
    cd(k,1) = sum(-cp.*del.*n1(:,1));
    WP_p(k,:) = [x_wp,y_wp];
    
    %% Determine the Velocity at the tail
    leng = sqrt((pos(k,1)-zx(1))^2+(pos(k,2)-zy(1))^2);
    Vb_tip(k,:) = [-Vb(k,1)*cos(pos(k,3))-Vb(k,2)*sin(pos(k,3)),...
        -Vb(k,1)*sin(pos(k,3))+Vb(k,2)*cos(pos(k,3))+leng*Vb(k,3)];
    delta_gamma(k,1) = (gamma(k)-gamma(k-1))/interval;
    if k > 100
        aveux = mean(Vb_tip(k-100:k,1));
    else
        aveux = Vb_tip(k,1);
    end
     text3 = sprintf('Vbx = %g ',aveux);
     disp(text3);
     
    %% Rotate Coordinates %% 
    if k > 100
        theta = mean(pos(k-100:k,3));
        if k > 100 && theta > 1*pi/180
            theta2(count_the) = theta;
            pos2(k,3) = pos(k,3) + theta;
            theta = mean(pos(k-100:k,3));

            vort3x = vortx2 - pos(k,1);
            vort3y = vorty2 - pos(k,2);  
            for j=k:k
                x2 = x(:,j) - pos(j,1);
                y2 = y(:,j) - pos(j,2);

                pos(j,3) = pos(j,3)-theta;
                for i=1:length(x2)
                    x(i,j) = [cos(theta) sin(theta)]*[x2(i);y2(i)] + pos(j,1);
                    y(i,j) = [-sin(theta) cos(theta)]*[x2(i);y2(i)] + pos(j,2);

                end
            end
            for i=1:length(vortx2)
                vort2x(i,1) = [cos(theta) sin(theta)]*[vort3x(i);vort3y(i)] + pos(k,1);
                vort2y(i,1) = [-sin(theta) cos(theta)]*[vort3x(i);vort3y(i)] + pos(k,2); 
            end
            vort_posx(1:length(vortx2),k) = vort2x;
            vort_posy(1:length(vortx2),k) = vort2y; 
            count_the = count_the + 1;
        end
    end

    %% Save run information after a number of iterations
    if k == 20000
        dat = sprintf('Amp_%.2f.mat',16);
        Velo=Vb(1:k,:);
        tau=Mb(1:k);
        position=pos(1:k,:);
        
        save(dat,'Velo','tau','position','accx','accy','omegad');
    end
end
    
    

