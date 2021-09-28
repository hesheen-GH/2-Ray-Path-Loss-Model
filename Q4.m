%Part 4 -------------------------------------------------------------------
clear all;
clc; 

Pt = 100;
Gt = 1;
Gr = 1; 
N = 100;
f = 2*10^9;
a = 100;
b = 1000;
ht = 10;
hr = 1;
rel_permittivity = 15;
field_polarization = 'v';


function_obj = Functions_Class;
function_obj.generate_pdf_and_cdf(Pt,Gt,Gr,N,f,a,b,ht,hr,rel_permittivity,field_polarization);

f = 10*10^9;
function_obj.generate_pdf_and_cdf(Pt,Gt,Gr,N,f,a,b,ht,hr,rel_permittivity,field_polarization);



