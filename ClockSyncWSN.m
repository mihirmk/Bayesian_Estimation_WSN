% Code for Clock Synchronization for Wireless Sensor Networks
clc
clear all
close all
load("ClockSync.mat");

tx = tx_timestamps;
rx = rx_timestamps;

% Estimate of phi2 based on the Clock Synchronisation Equation
x = rx - tx - tau;
% (x has a Gaussian noise distribution (contributions from rx & tx measurement errors)

% Phi_2 Scaled as per Data set
phi_2s = phi_2 .* ones(K,trial,6);
phi2 = phi_2 .* ones(10,6);
p = phi_2.* ones(1,6);

% Estimated phi - Mean of timestamps of K messages with trials
phi_mean = squeeze(mean(x,2));
% Standard Deviation for K messages with trials
phi_std  = squeeze(std(x,1,2));
% Bias in Measurements for K messages with trials
phi_bias = (phi_mean - phi2).^2;
% Variance in Measurements for K messages
phi_var = (1/trial) * phi_std.^2;

% Calculated MSE & CRLB for the Dataset (average over K messages)
se = (x - phi_2s).^2 ;% Squared Error
mse_c = mean(squeeze(mean(se,2)));
crlb_c = mean(phi_var);

% Theoretical MSE & CRLB for the Dataset (average over K messages)
mse_t = mean(phi_var + phi_bias);
crlb_t = (1/trial) * noise_var;

% Phi : Estimates and Standard deviations
phi_mm = mean(phi_mean);
phi_sm = mean(phi_std);
 
% Plot of Theoretical & Estimated CRLB vs. Noise Variance 
figure
loglog(noise_var, crlb_t,'Marker', 'o')
%hold on
%loglog(noise_var, crlb_c,'Marker', 'x')
title('CRLB of \phi with noise variance','FontName', 'Times')
xlabel('Noise \sigma^2','FontName', 'Times')
ylabel('CRLB','FontName', 'Times')
legend('Theoretical CRLB') %, 'Calculated CRLB')
saveas(gca,['CRLB_Variance.png']);

% Plot of MSE vs. Noise Variance 
figure
loglog(noise_var, mse_t,'Marker', 'o')
hold on
loglog(noise_var, mse_c,'Marker', 'x')
title('MSE of \phi with noise variance','FontName', 'Times')
xlabel('Noise \sigma^2','FontName', 'Times')
ylabel('MSE','FontName', 'Times')
legend('Theoretical MSE', 'Calculated MSE')
saveas(gca,['MSE_Variance.png']);

% Plot of Phi_Estimate vs. Noise Variance 
figure
semilogx(noise_var, p,'Marker','o')
hold on
semilogx(noise_var, phi_mm,'Marker', 'x')
title('Estimate of \phi with noise variance','FontName', 'Times')
xlabel('Noise \sigma^2','FontName', 'Times')
ylabel('Estimate E[\phi]','FontName', 'Times')
legend('Reference', 'Numerical Estimate')
saveas(gca,['Estimate_Variance.png']);

% Plot of Data Sets with different noise variances
min_x = min(x,[],'all');
max_x = max(x,[],'all');
figure('Renderer', 'painters', 'Position', [10 10 1200 600])
t = tiledlayout(2,3,'TileSpacing','none');
for i = 1 : 6
    nexttile
    x1 = x(:,:,i);
    grid minor
    % pdf distribution for plotting Observed Data
    x_pdf(i,:) = linspace(min_x,max_x,K*trial);
    y_pdf(i,:) = normpdf(x_pdf(i,:), phi_mm(i), phi_sm(i));
    xline(phi_mm(i),'-','Estimate');
    hold on
    histogram(x1(:), 'Normalization', 'pdf', 'EdgeColor','k' ,'FaceColor' , 'none');
    hold on
    plot(x_pdf(i,:), y_pdf(i,:), 'color', 'b' , 'LineWidth',1);
    title(['E[\phi] for 10k MC run; Noise \sigma^2  = ', num2str(phi_sm(i))],'FontName', 'Times');
    xlabel('Distribution','FontName', 'Times');
    ylabel(['PDF'],'FontName', 'Times');
end
saveas(gca,['Var_Comb.png']);

% Analysis of Variance for Different Standard Deviations of Noise 
figure
boxplot(phi_mean,'Labels', [noise_var]);
grid on
title('Variance of E[\phi] against \sigma^2')
ylabel('pdf Distribution of estimate of \phi')
xlabel('Noise \sigma^2')
saveas(gca,['VarBox.png']);

clear all