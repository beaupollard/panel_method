
global panels
global L
global V_inf 
global point_upstream


point_upstream = [-200 7];
xi = load('x.txt');
yi = load('y.txt');

Fx_tol = 1*10^(-5);
Fy_tol = 1*10^(-5);
M_f_tol = 1*10^(-5);

L = 0;
for i=1:length(xi)-1
    L = L + sqrt((xi(i+1)-xi(i))^2+(yi(i+1)-yi(i))^2);
end

interval = 0.01;
tf = 200;
t_length = round(tf/interval);

panels = length(xi)-1;

x = zeros(panels+1,t_length);
y = x;
x(:,1) = xi;
y(:,1) = yi;

cp = zeros(panels,t_length);
cl = zeros(t_length,1);
cd = cl;
phi = cp;

vort_posx = zeros(t_length,t_length);
vort_posy = zeros(t_length,t_length);

gamma = zeros(t_length,1);
gamma_w = gamma;
alpha = gamma;
delta = gamma;
time = gamma;

Vb = zeros(t_length,3);
pos = Vb;
V_inf = [0 0];
Vb(1,1) = 0;
Vb(2,1) = 0;
w = 2*pi;
A = 12;
Mb = zeros(t_length,1);
Mf_p = Mb;
U_Tailx = zeros(t_length,30);
U_Taily = U_Tailx;
Vb_tip = zeros(t_length,2);
WP_p = zeros(t_length,2);
conv = zeros(t_length,1);

