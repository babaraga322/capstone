clc
clear all
close all

%% INITIAL DATA
M = 64; % antennas at BS to serve Ku single-antenna UAVs
Ku = 10; % number of served UAVs
K = 10:10:150; % number of authorized UAVs, Ku < K
Ke = 10:10:150; % number of unathorized UAVs

r = [0 5 10 15];
Markers = {'-+','-o','-*','-x'};
Iterations = 1000;

d_min = 0;
d_max = 100;

z_min = 10;
z_max = 120;

x_BS = 0;
z_BS = 25; % the height of BS

fc = 4.2*10^6;
c_light = 3*10^8;

k_Rician_dB = 10; % Rician K-factor in dB
k_Rician = 10^(k_Rician_dB/10); % Rician K-factor in linear scale

SNR_Tx_dB = 10;
SNR_Tx = 10^(SNR_Tx_dB/10);
sigma2 = 1;

env_a = 9.61;
env_b = 0.28;

eta_LOS_dB = 1;
eta_LOS = 10^(eta_LOS_dB/10);
eta_NLOS_dB = 20;
eta_NLOS = 10^(eta_NLOS_dB/10);

%% SIMULATION
for rr=1:4
rp = r(rr);
mk = Markers{rr};
for ii=1:length(Ke)
    P_Tx = SNR_Tx * sigma2;
%     number_K = K(ii);
    number_K = K(ii);
    number_Ke = Ke(ii);
    
    R_SUM = [];
    ave_K = [];
    %% Authorized UAVs
    for ee=1:Iterations
    d_K = unifrnd(d_min,d_max,1,number_Ke);
    z_K = unifrnd(z_min,z_max,1,number_Ke);
    theta_K = atan((z_K - z_BS)./d_K);
    P_LOS_K = 1./(1 + env_a*exp(-env_b*(theta_K - env_a)));
    P_NLOS_K = 1 - P_LOS_K;
    
    beta_K_NLOS = 10*log10(d_K.^2 + (z_K - z_BS).^2) + 20*log10(4*pi*fc/c_light) + eta_NLOS_dB;
    beta_K_LOS = 10*log10(d_K.^2 + (z_K - z_BS).^2) + 20*log10(4*pi*fc/c_light) + eta_LOS_dB;
    
    %% LOS part
    H_K = [];
    H_Ke = [];
    for pp=1:number_K
        theta_K_pp = theta_K(pp);
        
        beta_K_NLOS_pp = beta_K_NLOS(pp);
        beta_K_LOS_pp = beta_K_LOS(pp);
        
        P_LOS_K_pp = P_LOS_K(pp);
        
        for jj=1:M
            h_bar_K_LOS(jj,1) = exp(-1j*pi*(jj - 1)*cos(theta_K_pp));
        end
        hh_bar_K_LOS = sqrt(k_Rician/(1 + k_Rician))*sqrt(beta_K_LOS_pp)*h_bar_K_LOS;
        
        %% NLOS part
        h_tilde_K_NLOS = sqrt(1/2)*(randn(M,1) + 1i*randn(M,1));
%         R_K_NLOS = h_tilde_K_NLOS * h_tilde_K_NLOS'; % actual spatial correlation
        R_K_NLOS = eye(M);
        hh_tilde_K_NLOS = sqrt(1/(1 + k_Rician))*sqrt(R_K_NLOS)*h_tilde_K_NLOS*sqrt(beta_K_LOS_pp);
        
        h_overall_K = P_LOS_K_pp * hh_bar_K_LOS + (1 - P_LOS_K_pp) * hh_tilde_K_NLOS;
        
        H_K = [H_K h_overall_K]; 
    end
    %% Unauthorized UAVs
%     number_Ke = Ke(ii);
    d_Ke = unifrnd(d_min,d_max,1,Ke(ii));
    z_Ke = unifrnd(z_min,z_max,1,Ke(ii));
    theta_Ke = atan((z_Ke - z_BS)./d_Ke);
    P_LOS_Ke = 1./(1 + env_a * exp(-env_b*(theta_Ke - env_a)));
    P_NLOS_Ke = 1 - P_LOS_Ke;
    
    beta_Ke_NLOS = 10*log10(d_Ke.^2 + (z_Ke - z_BS).^2) + 20*log10(4*pi*fc/c_light) + eta_NLOS_dB;
    beta_Ke_LOS = 10*log10(d_Ke.^2 + (z_Ke - z_BS).^2) + 20*log10(4*pi*fc/c_light) + eta_LOS_dB;
    %% LOS part of Unauthorized UAVs
    H_Ke = [];
    for pp=1:number_Ke
        theta_Ke_pp = theta_Ke(pp);
        
        beta_Ke_NLOS_pp = beta_Ke_NLOS(pp);
        beta_Ke_LOS_pp = beta_Ke_LOS(pp);
        
        P_LOS_Ke_pp = P_LOS_Ke(pp);
        
        for jj=1:M
            h_bar_Ke_LOS(jj,1) = sqrt(k_Rician/(1 + k_Rician))*exp(-1j*pi*(jj - 1)*cos(theta_Ke_pp));
        end
        hh_bar_Ke_LOS = sqrt(k_Rician/(1 + k_Rician))*sqrt(beta_K_LOS_pp)*h_bar_Ke_LOS;
        
        %% NLOS part
        h_tilde_Ke_NLOS = sqrt(1/2)*(randn(M,1) + 1i*randn(M,1));
        R_Ke_NLOS = h_tilde_Ke_NLOS * h_tilde_Ke_NLOS';
        hh_tilde_Ke_NLOS = sqrt(1/(1 + k_Rician))*sqrt(R_Ke_NLOS)*h_tilde_Ke_NLOS*sqrt(beta_K_LOS_pp);
        
        h_overall_Ke = P_LOS_Ke_pp * hh_bar_Ke_LOS + (1 - P_LOS_Ke_pp) * hh_tilde_Ke_NLOS;
        
        H_Ke = [H_Ke h_overall_Ke]; 
    end
    %% Protected zone
    number_Kp = number_Ke;
    index_valid_Kp = 1:number_Ke;
    for oo=1:number_K
        oo_d = d_K(oo);
        oo_z = z_K(oo);
        for ll=1:number_Ke
            ll_d = d_Ke(ll);
            ll_z = z_Ke(ll);
            res_ll = (oo_d - ll_d)^2 + (oo_z - ll_z)^2;
            if res_ll < rp^2
                number_Kp = number_Kp - 1;
                index_valid_Kp(index_valid_Kp==oo) = []; 
                break              
            end
        end
    end
    if number_Kp>Ku
        number_Kp = Ku;
        index_valid_Kp = index_valid_Kp(1:number_Kp);
    end
    
    H_Kp = H_K(:,index_valid_Kp);
    H_K_final = H_Kp.';
    W = H_K_final'/(H_K_final*H_K_final' + (number_Kp/(SNR_Tx))*eye(number_Kp));
    %% SNR Calculation of served UAVs
    for tt=1:number_Kp
        p_k(tt) = sqrt(P_Tx/number_Kp) * (norm(W(:,tt)))^(-1);
    end
    for ll=1:number_Kp
        numerator_ll = (p_k(ll))^2 * abs(H_K_final(ll,:)*W(:,ll))^2;
        denominator_ll = 0;
        for qq=1:number_Kp
            denominator_qq = (p_k(qq))^2 * abs(H_K_final(ll,:)*W(:,qq))^2;
            denominator_ll = denominator_ll + denominator_qq;
        end
        denominator_ll = denominator_ll - numerator_ll;
        SINR_K_UAV(ll) = numerator_ll/(denominator_ll + sigma2);
    end
    Rate_SINR_K_UAV = log2(1 + SINR_K_UAV);
    %% SNR Calculation of unauthorized UAVs
    no_of_gains_Ke = sum(abs(H_Ke));
    [sort_val_Ke, sort_ind_Ke] = sort(no_of_gains_Ke,'descend');
    new_ind_H_Ke = sort_ind_Ke(1:number_Kp);
    H_Ke_new = H_Ke(:,new_ind_H_Ke);
    H_Ke_new = H_Ke_new.';
    
    for ww=1:number_Kp
        numerator_ww = (p_k(ww))^2 * abs(H_Ke_new(ww,:)*W(:,ww))^2;
        denominator_ww = 0;
        for qqq=1:number_Kp
            denominator_qqq = (p_k(qqq))^2 * abs(H_Ke_new(ww,:)*W(:,qqq))^2;
            denominator_ww = denominator_ww + denominator_qqq;
        end
        denominator_ww = denominator_ww - numerator_ww;
        SINR_Ke_UAV(ww) = numerator_ww/(denominator_ww + sigma2);
    end
    Rate_SINR_Ke_UAV = log2(1 + SINR_Ke_UAV);
    R_sum = Rate_SINR_K_UAV - Rate_SINR_Ke_UAV;
    for www=1:number_Kp
        if R_sum(www) < 0
            R_sum(www) = 0;
        end
    end
    R_sum = sum(R_sum);
    R_SUM = [R_SUM R_sum];
    ave_K = [ave_K number_Kp];
    end
    Sum_secrecy(ii) = sum(R_SUM)/Iterations;
    Average_K(ii) = sum(ave_K)/Iterations;
    ii
end 

figure(1);
plot(Ke,Sum_secrecy,mk,'LineWidth',1);
grid on;
hold on
xlim([10 150])
xlabel('Number of Unauthorized UAVs (K_e)')
ylabel('Average sum of secrecy rates')
legend('r_p = 0 m','r_p = 5 m','r_p = 10 m','r_p = 15 m')

figure(2);
plot(Ke,Average_K,mk,'LineWidth',1);
hold on
xlim([10 150])
xlabel('Number of Unauthorized UAVs (K_e)')
ylabel('Average number of serving authorized UAVs (K_u)')
legend('r_p = 0 m','r_p = 5 m','r_p = 10 m','r_p = 15 m')
end
