clear all;
clc;

%Part 1
R = 1:1:10^4;
figure(1)
semilogx(R,10*log10(compute_path_loss(1*10^9, 2, 2, 16, 'v', 1:1:10^4,'2-ray')))
grid on
hold on
semilogx(R,10*log10(compute_path_loss(1*10^9, 2, 2, 16, 'h', 1:1:10^4,'2-ray')))
hold on
semilogx(R,-10*log10(compute_path_loss(1*10^9, 2, 2, 16, 'v', 1:1:10^4,'freespace')))

%Part 2

Path_Loss_2_ray = []
Path_Loss_free_space = []

rng(1,'twister');
figure(2)
a = 100;
b = 1000;
r = (b-a).*rand(100,1) + a;

%for i = 1:size(r)
    %Path_Loss_2_ray[i] = compute_path_loss(1*10^9, 10, 1, 15, 'v', r(i)
    
    

%edges = [100 1000];
%plot(histogram(r),20);
%histogram(r,36)
histogram(r)

Pt = 100
Gt = 1
Gr = 1


%figure(2)
%x = 1:100;
%plot(unifpdf(x,100,1000));
%xlim([50 1100])



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
    path_loss = ((abs(1+reflection_coeff.*(dD./(d1_d2)).*exp(1i*phase_diff))).^2).*(1./((4*pi.*R)./(3*10^8/freq))).^2;
    
  
end 

end 