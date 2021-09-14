function [gamma,N_Conv,gold] = Newton_It(a,b,c,D1,D2)

global k

%Newton iteration method
N_Conv = 0;
tol = 0.0001;
%% Find the zero
gz = -b/(2*a);
fgz = a*gz^2+b*gz+c;
gold = 0;
rfgz = 2*a*(gz-1)+b;
if k == 148
end
if fgz*rfgz < 0
    N_Conv = 1;
    disp('No real solution')
end

if N_Conv == 0
    err_x = 1;
    x1(1) = gz+1;
    i = 2;
    counter = 1;
    while(err_x > tol)
        x1(i) = x1(i-1)-(a*x1(i-1)^2+b*x1(i-1)+c)/(a*2*x1(i-1)+b);
        err_x = abs(x1(i)-x1(i-1));
        err_x2(i) = err_x;
        i = i+1;
        if i > 10000 && N_Conv ~= 1 && err_x/x1(end) < tol
            err_x = 0;
        end
    end
    if k == 13
    end
    dg(1,1) = x1(end);
    err_x = 1;
    x1 = [];
    x1(1) = gz-1;
    i = 2;
    while(err_x > tol)
        x1(i) = x1(i-1)-(a*x1(i-1)^2+b*x1(i-1)+c)/(a*2*x1(i-1)+b);
        err_x = abs(x1(i)-x1(i-1));
        i = i+1;
        if i > 10000 && N_Conv ~= 1 && err_x/x1(end) < tol
            err_x = 0;
        end
    end

    dg(1,2) = x1(end);

    Vt1 = dg(1,1).*D1+D2;
    Vt2 = dg(1,2).*D1+D2;
    Vt1_s = Vt1(1)*Vt1(end);
    Vt2_s = Vt2(1)*Vt2(end);

    check1 = Vt1(1)*Vt1(end);
    check2 = Vt2(1)*Vt2(end);
    if check1 < 0
        gamma = dg(1,1);
        gold = dg(1,2);
    elseif check2 < 0
        gamma = dg(1,2);
        gold = dg(1,1);
    else
        if check2 > check1
            gamma = dg(1,1);
            gold = dg(1,2);
        else
            gamma = dg(1,2);
            gold = dg(1,1);
        end

    end
else
    gamma = 0;
end

end

