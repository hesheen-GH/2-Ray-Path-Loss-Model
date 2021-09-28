classdef Functions_Class
    properties (Constant)
        c = 3*10^8;
    end
    
    methods
        
        function [] = q7_pdf_cdf(obj, N)
        Pt = 100;
        Gt = 1;
        Gr = 1;    
        a = 0;
        b = 2000;
        x_array = [];
        y_array = [];
        Pr_2_ray = [];
        Pr_LOS = [];
        Path_Loss_2_ray = [];
        Path_Loss_LOS = [];
        x = 1000.*ones(1,N);
        y = 1000.*ones(1,N);

        for i = 1:N

            while((x(i)>950 && x(i)<1050) && (y(i)>950 && y(i)<1050))

                x(i) = (b-a).*rand(1,1) + a;
                y(i) = (b-a).*rand(1,1) + a;

            end 

            x_array(i) = x(i);
            y_array(i) = y(i);

            Path_Loss_2_ray(i) = obj.main_menu(1*10^9, 10, 1, 15, 'v', sqrt((1000-x_array(i)).^2+(1000-y_array(i)).^2), '2-ray', 'n');
            Path_Loss_LOS(i) = obj.main_menu(1*10^9, 10, 1, 15, 'v', sqrt((1000-x_array(i)).^2+(1000-y_array(i)).^2), 'freespace','n');


        end 

        Pr_2_ray = 10*log10(1000*(Pt*Gt*Gr*(1./Path_Loss_2_ray)));
        Pr_LOS = 10*log10(1000*(Pt*Gt*Gr*(1./Path_Loss_LOS)));

        figure;
        histogram2(x_array, y_array, [40 40]);

        figure;
        histogram(Pr_2_ray, 'Normalization','probability');
        hold on;
        histogram(Pr_LOS, 'Normalization','probability');
        xlabel('Recieved Power (dBm)');
        ylabel('Probability');
        title('Q7 Empirical PDF for N = ' + string(N));

        figure;
        cdfplot(Pr_2_ray);
        hold on;
        cdfplot(Pr_LOS);
        xlabel('Recieved Power (dBm)');
        ylabel('Probability');
        title('Q7 Empirical CDF for N = ' + string(N))
            
        end 
        
        
        function [] = generate_pdf_and_cdf(obj, Pt, Gt, Gr, N, f, a, b, tx_height, rx_height, rel_permittivity, field_polarization)
            
        Pr_2_ray = [];
        Pr_LOS = [];
        Path_Loss_2_ray = [];
        Path_Loss_LOS = [];
        rng(1,'twister');
        r = (b-a).*rand(N,1) + a;
        
        for i = 1:size(r)
            Path_Loss_2_ray(i) = obj.main_menu(f, tx_height, rx_height, rel_permittivity, field_polarization, r(i), '2-ray','n');
            Path_Loss_LOS(i) = obj.main_menu(f, tx_height, rx_height, rel_permittivity, field_polarization, r(i), 'freespace','n');
        end
        
        figure1 = figure;
        Pr_2_ray = 10*log10(1000*(Pt*Gt*Gr*(1./Path_Loss_2_ray)));
        Pr_LOS = 10*log10(1000*(Pt*Gt*Gr*(1./Path_Loss_LOS)));
        
        histogram(Pr_2_ray, 'Normalization','probability');
        hold on;
        histogram(Pr_LOS, 'Normalization','probability');
        xlabel('Recieved Power (dBm)');
        ylabel('Probability');
        str = 'Empirical PDF for N = ' + string(N)+ ' f = ' + (string(f/(10^9))) + ' GHz  Field polarization = ' + string(field_polarization) + ' Rmin = ' + string(a)+ 'm';
        title(str); 
        legend('2-ray vertical polar.', 'LOS')
        hold off;
        %saveas(figure1,str+'.jpg');

        figure2 = figure;
        cdfplot(Pr_2_ray);
        hold on;
        cdfplot(Pr_LOS);
        xlabel('Recieved Power (dBm)');
        ylabel('Probability');
        str = 'Empirical CDF for N = ' + string(N)+ ' f = ' + (string(f/(10^9))) + ' GHz  Field polarization = ' + string(field_polarization) + ' Rmin = ' + string(a)+ 'm';
        title(str);
        legend('2-ray vertical polar.', 'LOS');
        hold off;
        %saveas(figure2,str+'.jpg');
        
        end 

        function [path_loss] = main_menu(obj, freq, tx_height, rx_height, rel_permittivity, field_polarization, R, link_type, far_field_approx)

        if strcmp(link_type, 'freespace')
            path_loss = obj.LOS_path_loss(freq,R);
        end 

        if strcmp(link_type, '2-ray_approx')
            path_loss = obj.two_ray_path_loss_approx(tx_height,rx_height,R);
        end 

        if strcmp(link_type, '2-ray')
                if strcmp(far_field_approx, 'y')
                    if R<((20*tx_height*rx_height)/((obj.c)/freq)) %far-field approximation
                        path_loss = 0;
                    else
                        path_loss = obj.two_ray_path_loss(freq,tx_height,rx_height, rel_permittivity, field_polarization, R);
                    end

                else
                    path_loss = obj.two_ray_path_loss(freq,tx_height,rx_height, rel_permittivity, field_polarization, R);
                end  
        end 
        end
        
        function [path_loss] = two_ray_path_loss(obj, freq, tx_height, rx_height, rel_permittivity, field_polarization, R)

        dD = sqrt(R.^2+(tx_height-rx_height)^2); %R1
        incidence_angle = atan((tx_height+rx_height)./R);
        d1_d2 = sqrt((tx_height+rx_height)^2 +R.^2); %R2

        if strcmp(field_polarization,'v')
            z = (1/rel_permittivity)*sqrt(rel_permittivity-(cos(incidence_angle)).^2);

        else
            z = sqrt(rel_permittivity-(cos(incidence_angle)).^2);
        end

        reflection_coeff = (sin(incidence_angle)-z)./(sin(incidence_angle)+z);
        phase_diff = 2*pi*(1/(obj.c/freq)).*(d1_d2-dD);
        path_loss = 1./(((abs(1+reflection_coeff.*(dD./(d1_d2)).*exp(1i*phase_diff))).^2).*(1./((4*pi.*R)./(obj.c/freq))).^2);

        return
        end 


        function [path_loss] = two_ray_path_loss_approx(obj, tx_height, rx_height, R)

        path_loss = (R.^4)/((tx_height^2)*(rx_height^2)); 

        end


        function [path_loss] = LOS_path_loss(obj, freq, R)

        path_loss = ((4*pi.*R)/((obj.c)/freq)).^2;
        return

        end
        
    end
end