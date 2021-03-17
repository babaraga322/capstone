clc
clear all
close all

M = 1;
steps = 2*10^3;
 



x_min = 0;
x_max = 150;

z_min = 10;
z_max = 120;

x_BS = 0;
z_BS = 25;

x_D = 150;
z_D = 25;




sigma2 = 1;
L = 100;
values_K = 0:10:100;
values_r_p = 0:2:20;

SNR_Tx_dB = 30;

Sec_cap_3d = [];
for K = values_K
    K
    
    Avg_cap_UAV = [];
    Avg_cap_E = [];

    for  r_p = values_r_p 

        SNR_Tx = 10^(SNR_Tx_dB/10); 

        cap_UAV = [];
        cap_E = [];
        for step = 1:steps

            active = 1:K;

            x_K = unifrnd(x_min,x_max,1,K);
            z_K = unifrnd(z_min,z_max,1,K);

            x_L = unifrnd(x_min,x_max,1,L);
            z_L = unifrnd(z_min,z_max,1,L);

            x_E = unifrnd(x_min,x_max,1,1);
            z_E = 25;

            for i=1:K
               for j=1:L
                   dist_K_L = sqrt((x_K(i) - x_L(j))^2 + (z_K(i) - z_L(j))^2);
                   if dist_K_L < r_p
                      active(active == i) = [];
                   end
               end
            end

            for i=active
                dist_K_E = sqrt((x_K(i) - x_E)^2 + (z_K(i) - z_E)^2);
                if dist_K_L < r_p
                    active(active == i) = [];
                end
            end

            H_BS_K = [];
            H_K_D = [];
            H_K_E = [];

                for i = 1:K
                    h = sqrt(1/2)*(randn(M,1) + 1i*randn(M,1));
                    dist_BS_K = sqrt((x_K(i) - x_BS)^2 + (z_K(i) - z_BS)^2);
                    P_tx = SNR_Tx * sigma2;
                    S_tx = P_tx * abs(h)^2/ (sigma2 * dist_BS_K);
                    H_BS_K = [H_BS_K S_tx];
                end

                for i = 1:K
                    h = sqrt(1/2)*(randn(M,1) + 1i*randn(M,1));
                    dist_K_D = sqrt((x_K(i) - x_D)^2 + (z_K(i) - z_D)^2);
                    P_tx = SNR_Tx * sigma2;
                    S_tx = P_tx * abs(h)^2/ (sigma2 * dist_K_D);
                    H_K_D = [H_K_D S_tx];
                end
            if ~isempty(active)
                H = [H_BS_K; H_K_D];
                active_H = min(H(:,active));
                selected_UAV = max(active_H);
                index_of_selected = active((active_H==selected_UAV));


                g = sqrt(1/2)*(randn(M,1) + 1i*randn(M,1));
                dist_BS_E = sqrt((x_BS - x_E)^2 + (z_BS - z_E)^2);
                P_tx = SNR_Tx * sigma2;
                S_tx = (P_tx * abs(g)^2)/ (P_tx * abs(g)^2+(sigma2 * dist_BS_E));

                H_K_E = [H_K_E S_tx];

                cap_UAV = [cap_UAV log2(1+selected_UAV)];
                cap_E = [cap_E log2(1+H_K_E)];
            else
                cap_UAV = [cap_UAV 0];
                cap_E = [cap_E 0];

            end
        end
        Avg_cap_UAV = [Avg_cap_UAV sum(cap_UAV)/steps];
        Avg_cap_E = [Avg_cap_E sum(cap_E)/steps];

    end
    Sec_capac = max(Avg_cap_UAV - Avg_cap_E,0);
    Sec_cap_3d = [Sec_cap_3d; Sec_capac];
end



surf(values_r_p,values_K,Sec_cap_3d)
xlabel('r_p')
ylabel('number of authorized UAV, K')
zlabel('Secrecy capacity')
