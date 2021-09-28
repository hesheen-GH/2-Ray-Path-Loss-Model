%Part 1 -------------------------------------------------------------------
clear all;
clc;

Lp_2_ray_v_polarization = [];
Lp_2_ray_h_polarization = [];
Lp_2_ray_approx = [];
Lp_LOS = [];
function_obj = Functions_Class;
f = 3*10^8;

for R=1:1:10^4
    Lp_2_ray_v_polarization(R)=function_obj.main_menu(f, 2, 2, 16, 'v', R,'2-ray','n');
    Lp_2_ray_h_polarization(R)=function_obj.main_menu(f, 2, 2, 16, 'h', R,'2-ray','n');
    Lp_LOS(R) = function_obj.main_menu(f, 2, 2, 16, 'h', R,'freespace','n');
    Lp_2_ray_approx(R) = function_obj.main_menu(f, 2, 2, 16, 'h', R,'2-ray_approx','n');
end 

figure(1)
semilogx(-10*log10(Lp_2_ray_v_polarization));
grid on;
hold on
semilogx(-10*log10(Lp_2_ray_h_polarization));
hold on;
semilogx(-10*log10(Lp_LOS));
hold on;
semilogx(-10*log10(Lp_2_ray_approx));
xlim([1 10000]);
ylim([-160 1]);
legend('2-ray vertical polarization', '2-ray horizontal polarization', 'LOS', '2-ray approximation');
xlabel('distance (m)');
ylabel('Path gain (dB)');
title('2-ray and LOS for Part 1');
hold off; 
