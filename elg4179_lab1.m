clear all;
clc;

%Part 1
R = 1:1:10^4;
figure(1)
semilogx(R,-10*log10(compute_path_loss(3*10^8, 2, 2, 16, 'v', 1:1:10^4,'2-ray')))
grid on
hold on
semilogx(R,-10*log10(compute_path_loss(3*10^8, 2, 2, 16, 'h', 1:1:10^4,'2-ray')))
hold on
semilogx(R,-10*log10(compute_path_loss(3*10^8, 2, 2, 16, 'v', 1:1:10^4,'freespace')))

%Part 2

Pr_2_ray = [];
Pr_LOS = [];
Path_Loss_2_ray = [];
Path_Loss_LOS = [];
Pt = 100;
Gt = 1;
Gr = 1;
N = 100; 
rng(1,'twister');
a = 100;
b = 1000;
r = (b-a).*rand(N,1) + a;

for i = 1:size(r)
    Path_Loss_2_ray(i) = compute_path_loss(1*10^9, 10, 1, 15, 'v', r(i), '2-ray');
    Path_Loss_LOS(i) = compute_path_loss(1*10^9, 10, 1, 15, 'v', r(i), 'freespace');

end
    
Pr_2_ray = Pt*Gt*Gr*(1./Path_Loss_2_ray);
Pr_LOS = Pt*Gt*Gr*(1./Path_Loss_LOS);

figure(2);
histogram(Pr_2_ray);
hold on;
histogram(Pr_LOS);

figure(3)
cdfplot(Pr_2_ray);
hold on;
cdfplot(Pr_LOS);


%Part 7 
N = 100000
a = 0;
b = 2000;
x_array = [];
y_array = [];
x = 1000.*ones(1,N);
y = 1000.*ones(1,N);


for i = 1:N
    
    while((x(i)>950 && x(i)<1050) && (y(i)>950 && y(i)<1050))
        
        x(i) = (b-a).*rand(1,1) + a;
        y(i) = (b-a).*rand(1,1) + a;
        
    end 
       
    x_array(i) = x(i);
    y_array(i) = y(i);

end 

figure(4)
histogram2(x_array, y_array, [40 40]);



%hist(10*log10(100*compute_path_loss(1*10^9, 2, 2, 15, 'v', (1000-100).*rand(1000,1) + 0.1,'2-ray')))
%hold on
%hist(10*log10(compute_path_loss(1*10^9, 2, 2, 15, 'v', (1000-0.1).*rand(1000,1) + 0.1,'freespace')))


    
function [path_loss] = compute_path_loss(freq, tx_height, rx_height, rel_permittivity, field_polarization, R, link_type)

if strcmp(link_type, 'freespace')
    %path_loss = (((3*10^8)/freq)^2)./(16*pi^2*R.^2);
    path_loss = ((4*pi.*R)/((3*10^8)/freq)).^2;
    return
    
else 
    
    dD = sqrt(R.^2+(tx_height-rx_height)^2); %R1
    incidence_angle = atan((tx_height+rx_height)./R);
    d1_d2 = sqrt((tx_height+rx_height)^2 +R.^2); %R2

    if strcmp(field_polarization,'v')
        z = (1/rel_permittivity)*sqrt(rel_permittivity-(cos(incidence_angle)).^2);

    else
        z = sqrt(rel_permittivity-(cos(incidence_angle)).^2);
    end

    reflection_coeff = (sin(incidence_angle)-z)./(sin(incidence_angle)+z);
    phase_diff = 2*pi*(1/(3*10^8/freq)).*(d1_d2-dD);
    path_loss = 1./(((abs(1+reflection_coeff.*(dD./(d1_d2)).*exp(1i*phase_diff))).^2).*(1./((4*pi.*R)./(3*10^8/freq))).^2);
    
  
end 

end 