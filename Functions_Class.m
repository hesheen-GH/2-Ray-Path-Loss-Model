classdef Functions_Class
    properties (Constant)
        c = 3*10^8; %speed of light
    end
    
    methods
        
        function [] = q7_pdf_cdf(obj, N)
        Pt = 100; %transmitter power, antenna gains
        Gt = 1;
        Gr = 1;    
        a = 0; 
        b = 2000; % used to generate numbers between b and a
        x_array = []; %initalize random variables array on the x-axis
        y_array = []; % initalize random variables array on the y-axis
        Pr_2_ray = []; %initalize recieved power for 2 ray model
        Pr_LOS = []; %initalize recieved power for line of sight model
        Path_Loss_2_ray = []; %initalize path loss arrays
        Path_Loss_LOS = [];
        x = 1000.*ones(1,N); %initalize x,y to 1000 so it can be randomized right away
        y = 1000.*ones(1,N);

        for i = 1:N

            while((x(i)>950 && x(i)<1050) && (y(i)>950 && y(i)<1050)) %check if x,y coordinates inside 0.1x0.1 central square
                
                x(i) = (b-a).*rand(1,1) + a; %if inside central square , randomize again between 0 and 2000meters 
                y(i) = (b-a).*rand(1,1) + a;

            end 

            x_array(i) = x(i); %add randomly generated x,y to x,y array
            y_array(i) = y(i);
            
            %compute path loss, R is calculated using pythagorean theorem 
            % R = sqrt((x2-x1)^2+(y2-y1)^2) where x2 and y2 is the base
            % transmitter located at 1000,1000
            Path_Loss_2_ray(i) = obj.main_menu(1*10^9, 10, 1, 15, 'v', sqrt((1000-x_array(i)).^2+(1000-y_array(i)).^2), '2-ray', 'n');  
            Path_Loss_LOS(i) = obj.main_menu(1*10^9, 10, 1, 15, 'v', sqrt((1000-x_array(i)).^2+(1000-y_array(i)).^2), 'freespace','n');


        end 
        
        %compute recieved power using link budget equations
        Pr_2_ray = 10*log10(1000*(Pt*Gt*Gr*(1./Path_Loss_2_ray)));
        Pr_LOS = 10*log10(1000*(Pt*Gt*Gr*(1./Path_Loss_LOS)));

        figure;
        histogram2(x_array, y_array, [40 40]);

        figure;
        h1 = histogram(Pr_2_ray, 'Normalization','probability');
        hold on;
        h2 = histogram(Pr_LOS, 'Normalization','probability');
        morebins(h1);
        morebins(h1);
        morebins(h1);
        morebins(h2);
        morebins(h2);
        morebins(h2);
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
        
        Path_Loss_2_ray = rmmissing(Path_Loss_2_ray); %used to remove NaN values in the array in the case of invalidating the far field condition
        
        Pr_2_ray = 10*log10(1000*(Pt*Gt*Gr*(1./Path_Loss_2_ray)));
        Pr_LOS = 10*log10(1000*(Pt*Gt*Gr*(1./Path_Loss_LOS)));
        
        
        figure1 = figure;
        h1 = histogram(Pr_2_ray, 'Normalization','probability');
        hold on;
        h2 = histogram(Pr_LOS, 'Normalization','probability');
        morebins(h1);
        morebins(h1);
        morebins(h1);
        morebins(h2);
        morebins(h2);
        morebins(h2);
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
        %user selection, depending on propagation path model chosen

        if strcmp(link_type, 'freespace')
            path_loss = obj.LOS_path_loss(freq,R);
        end 

        if strcmp(link_type, '2-ray_approx')
            if strcmp(far_field_approx, 'y')
                if R<((20*tx_height*rx_height)/((obj.c)/freq)) %far field condition
                    path_loss = NaN;
                else
                    path_loss = obj.two_ray_path_loss_approx(tx_height,rx_height,R);
                end 
                
            else
                path_loss = obj.two_ray_path_loss_approx(tx_height,rx_height,R);
                
            end
        end 
        
        if strcmp(link_type, '2-ray')
                
            path_loss = obj.two_ray_path_loss(freq,tx_height,rx_height, rel_permittivity, field_polarization, R);

        end 
        end
        
        function [path_loss] = two_ray_path_loss(obj, freq, tx_height, rx_height, rel_permittivity, field_polarization, R)
        %compute 2-ray path loss using percise formula

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
        %compute 2-ray path loss approximation

        path_loss = (R.^4)/((tx_height^2)*(rx_height^2)); 

        end


        function [path_loss] = LOS_path_loss(obj, freq, R)
        %compute line of sight path loss

        path_loss = ((4*pi.*R)/((obj.c)/freq)).^2;
        return

        end
        
    end
end