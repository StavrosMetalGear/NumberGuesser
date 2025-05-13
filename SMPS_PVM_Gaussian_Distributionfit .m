clear all;
close all;
filename='SMPS_data_20241120';
filename2='SMPS_data_20241120_inside';
SMPS_data = xlsread(filename);
x_diam=SMPS_data(2,1:40);
SMPS_data_incloud=xlsread(filename2, 2);
%%3th BEFORE CLOUD EVENT
D3a=SMPS_data(19000:19193,1:40);
D3b=SMPS_data_incloud(26:1345,1:40);
S3a=size(D3a);
rws=S3a(1,1);

mm3a=[];
for i=1:12:rws
    mm3a=[mm3a;mean(D3a(i:i+11,1:end))]
    if i+24>rws
        disp('break')
        break;
    end
    disp(i+11)
end
S3b=size(D3b);
rws3b=S3b(1,1);

mm3b_m=[];
for i=1:12:rws3b
mm3b_m=[mm3b_m;mean(D3b(i:i+11,1:end))]
end
D4b=SMPS_data_incloud(1345:1718,1:40);
D4a=SMPS_data(20968:21000,1:40);
D4c=SMPS_data(21374:21406,1:40);

%%%BEFORE CLOUD

S4a=size(D4a);
rws4a=S4a(1,1);
mm4a_m=[];
for i=1:12:rws4a
    mm4a_m=[mm4a_m;mean(D4a(i,1:end))]
end
S4b=size(D4b);
rws4b=S4b(1,1);
mm4b_m=[];
for i=1:12:rws4b
    mm4b_m=[mm4b_m;mean(D4b(i:i+11,1:end))]
    if i+24>rws4b
        disp('break')
        break;
    end
    disp(i+11)
end
S4c=size(D4c);
rws4c=S4c(1,1);
mm4c_m=[];
for i=1:12:rws4c
    mm4c_m=[mm4c_m;mean(D4c(i:i+11,1:end))]
    if i+24>rws4c
        disp('break')
        break;
    end
    disp(i+11)
end
D5a=SMPS_data(21453:21582,1:40);
D5b=SMPS_data_incloud(1718:2004,1:40);
S5a=size(D5a);
rws5a=S5a(1,1);
mm5a_m=[];
for i=1:12:rws5a
    mm5a_m=[mm5a_m;mean(D5a(i:i+11,1:end))]
    if i+24>rws5a
        disp('break')
        break;
    end
    disp(i+11)
end
S5b=size(D5b);
rws5b=S5b(1,1);
mm5b_m=[];
for i=1:12:rws5b
    mm5b_m=[mm5b_m;mean(D5b(i:i+11,1:end))]
    if i+24>rws5b
        disp('break')
        break;
    end
    disp(i+11)
end
all_D=SMPS_data(1:end,1:40);

all_S5a=size(all_D);
all_rws=all_S5a(1,1);
all_mma_m=[];
for i=1:12:all_rws
    all_mma_m=[all_mma_m;mean(all_D(i:i+11,1:end))]
    if i+24>all_rws
        disp('break')
        break;
    end
    disp(i+11)
end
All_mma_m_round=round([x_diam;all_mma_m]);

cloud_mma_m_round = [x_diam; mm3a; mm3b_m; mm4b_m; mm4c_m; mm5a_m; mm5b_m];
f_cloud = [];
coefficientall_cloud = [];
szall_cloud = size(cloud_mma_m_round);
srall_cloud = szall_cloud(1, 1);
coefficientall_nuc_cloud = [];
coefficientall_aitk_cloud = [];
coefficientall_accum_cloud = [];
mode1_cloud = [];
mode2_cloud = [];
mode3_cloud = [];
total_hour_number = [];
fit_parameters = [];
goodnessoffit = [];
cloud_total_hour_number= [];
coeffvaluess_cloud=[];
title = ["1/sigmaswrt(2pi)", "geom.mean1", "2sigma^2", "1/sigma1swrt(2pi)", "geom.mean2", "2sigma2^2", "1/sigma3swrt(2pi)", "geom.mean3", "2sigma3^2", "Total Number", "Mode 1%", "Mode 2%", "Mode 3%", "goodnessoffit"];

for i = 1:srall_cloud - 1
    if any(isnan(cloud_mma_m_round(1 + i, 1:end)))
        continue;
    elseif sum(cloud_mma_m_round(1 + i, 1:end),'all')>40000
        continue;
        
        
    end
    
    f_all_i_cloud = fit(x_diam.', cloud_mma_m_round(1 + i, 1:end)', 'gauss3');
    yfit = f_all_i_cloud(x_diam);
    ydata = cloud_mma_m_round(1 + i, 1:end)';
    gof = goodnessOfFit(yfit, ydata, 'NRMSE');
    goodnessoffit = [goodnessoffit; gof];

    coeffvaluess_cloud = coeffvalues(f_all_i_cloud);
    fit_parameters = [fit_parameters; coeffvaluess_cloud];
    x = 0;
    y = 0;
    z = 0;
    
    % Process coefficients < 40 and >= 0
    for j = [2, 5, 8]
        if coeffvaluess_cloud(1, j) < 40 && coeffvaluess_cloud(1, j) >= 0
            x = x + 1;
            if x == 1
                coefficientall_nuc_cloud = [coefficientall_nuc_cloud; coeffvaluess_cloud(1, j)];
                disp("Adding to coefficientall_nuc_cloud: " + coeffvaluess_cloud(1, j))
            elseif x == 2
                coefficientall_nuc_cloud(end) = (coefficientall_nuc_cloud(end) + coeffvaluess_cloud(1, j)) / 2;
                disp("Updating coefficientall_nuc_cloud: " + coefficientall_nuc_cloud(end))
            elseif x == 3
                coefficientall_nuc_cloud(end) = (coefficientall_nuc_cloud(end) * 2 + coeffvaluess_cloud(1, j)) / 3;
                disp("Updating coefficientall_nuc_cloud: " + coefficientall_nuc_cloud(end))
            end
        elseif coeffvaluess_cloud(1, j) < 100 && coeffvaluess_cloud(1, j) >= 0
            y = y + 1;
            if y == 1
                coefficientall_aitk_cloud = [coefficientall_aitk_cloud; coeffvaluess_cloud(1, j)];
                disp("Adding to coefficientall_aitk_cloud: " + coeffvaluess_cloud(1, j))
            elseif y == 2
                coefficientall_aitk_cloud(end) = (coefficientall_aitk_cloud(end) + coeffvaluess_cloud(1, j)) / 2;
                disp("Updating coefficientall_aitk_cloud: " + coefficientall_aitk_cloud(end))
            elseif y == 3
                coefficientall_aitk_cloud(end) = (coefficientall_aitk_cloud(end) * 2 + coeffvaluess_cloud(1, j)) / 3;
                disp("Updating coefficientall_aitk_cloud: " + coefficientall_aitk_cloud(end))
            end
        elseif coeffvaluess_cloud(1, j) > 100 && coeffvaluess_cloud(1, j) <= 300
            z = z + 1;
            if z == 1
                coefficientall_accum_cloud = [coefficientall_accum_cloud; coeffvaluess_cloud(1, j)];
                disp("Adding to coefficientall_accum_cloud: " + coeffvaluess_cloud(1, j))
            elseif z == 2
                coefficientall_accum_cloud(end) = (coefficientall_accum_cloud(end) + coeffvaluess_cloud(1, j)) / 2;
                disp("Updating coefficientall_accum_cloud: " + coefficientall_accum_cloud(end))
            elseif z == 3
                coefficientall_accum_cloud(end) = (coefficientall_accum_cloud(end) * 2 + coeffvaluess_cloud(1, j)) / 3;
                disp("Updating coefficientall_accum_cloud: " + coefficientall_accum_cloud(end))
            end
        end
    end
    
    if x == 0
        coefficientall_nuc_cloud = [coefficientall_nuc_cloud; 0];
        disp("Adding 0 to coefficientall_nuc_cloud")
    end
    if y == 0
        coefficientall_aitk_cloud = [coefficientall_aitk_cloud; 0];
        disp("Adding 0 to coefficientall_aitk_cloud")
    end
    if z == 0
        coefficientall_accum_cloud = [coefficientall_accum_cloud; 0];
        disp("Adding 0 to coefficientall_accum_cloud")
    end
    c3bi_cloud = sort([coeffvaluess_cloud(1,2), coeffvaluess_cloud(1,5), coeffvaluess_cloud(1,8)]);
    st3bi_cloud = sort([coeffvaluess_cloud(1,3), coeffvaluess_cloud(1,6), coeffvaluess_cloud(1,9)]);
    
    % Ensure disjoint intervals for modes
    mode1i_cloud = cloud_mma_m_round(i, cloud_mma_m_round(1,1:end) < c3bi_cloud(1,1) - abs(sqrt(st3bi_cloud(1,1)/2)));
    mode2i_cloud = cloud_mma_m_round(i, ...
        (cloud_mma_m_round(1,1:end) > c3bi_cloud(1,2) - abs(sqrt(st3bi_cloud(1,2)/2))) & ...
        (cloud_mma_m_round(1,1:end) < c3bi_cloud(1,2) + abs(sqrt(st3bi_cloud(1,2)/2))));
    mode3i_cloud = cloud_mma_m_round(i, ...
        (cloud_mma_m_round(1,1:end) > c3bi_cloud(1,3) - abs(sqrt(st3bi_cloud(1,3)/2))) & ...
        (cloud_mma_m_round(1,1:end) < c3bi_cloud(1,3) + abs(sqrt(st3bi_cloud(1,3)/2))));
    
    Tot_numi_cloud = sum(cloud_mma_m_round(i,1:end), 'all'); % Total number for the row
    cloud_total_hour_number = [cloud_total_hour_number; Tot_numi_cloud];
    
    % Mode 1
    if ~isempty(mode1i_cloud)
        sum_mode1i_cloud = mean(mode1i_cloud); % Take the average of values
        number_mode1i_cloud = sum_mode1i_cloud / Tot_numi_cloud;
        mode1_cloud = [mode1_cloud; number_mode1i_cloud];
    else
        mode1_cloud = [mode1_cloud; 0];
    end
    
    % Mode 2
    if ~isempty(mode2i_cloud)
        sum_mode2i_cloud = mean(mode2i_cloud); % Take the average of values
        number_mode2i_cloud = sum_mode2i_cloud / Tot_numi_cloud;
        mode2_cloud = [mode2_cloud; number_mode2i_cloud];
    else
        mode2_cloud = [mode2_cloud; 0];
    end
    
    % Mode 3
    if ~isempty(mode3i_cloud)
        sum_mode3i_cloud = mean(mode3i_cloud); % Take the average of values
        number_mode3i_cloud = sum_mode3i_cloud / Tot_numi_cloud;
        mode3_cloud = [mode3_cloud; number_mode3i_cloud];
    else
        mode3_cloud = [mode3_cloud; 0];
    end
end

% Post-processing
Nuclei_geom_mean_cloud = round(coefficientall_nuc_cloud);
Aitken_geom_mean_cloud = round(coefficientall_aitk_cloud);
Accumula_geom_mean_cloud = round(coefficientall_accum_cloud);

Nuclei_persentage_cloud = round(mode1_cloud, 2);
Aitken_persentage_cloud = round(mode2_cloud, 2);
Accumula_persentage_cloud = round(mode3_cloud, 2);
filename7 = 'fit_parameters.txt';
%writetable(fit_parameters,filename7);

%%%%Gaussian fit for all the data
% Nuclei_persentage_cloud = round(mode1_cloud, 2);
% Aitken_persentage_cloud = round(mode2_cloud, 2);
% Accumula_persentage_cloud = round(mode3_cloud, 2);
f_all = [];
coefficientall = [];
szall = size(All_mma_m_round);
srall = szall(1,1);
coefficientall_nuc = [];
coefficientall_aitk = [];
coefficientall_accum = [];
mode1 = [];
mode2 = [];
mode3 = [];
allgoodnessoffit = [];
All_fit_parameters = [];
All_total_hour_number = [];
coeffvaluess=[];
title = ["1/sigmaswrt(2pi)", "geom.mean1", "2sigma^2", "1/sigma1swrt(2pi)", "geom.mean2", "2sigma2^2", "1/sigma3swrt(2pi)", "geom.mean3", "2sigma3^2", "Total Number", "Mode 1%", "Mode 2%", "Mode 3%", "goodnessoffit"];

for i = 1:srall-1
    if any(isnan(All_mma_m_round(1+i, 1:end)))
        continue;
    elseif sum(All_mma_m_round(1 + i, 1:end),'all')>40000
        continue;
    end
    
    f_all_i = fit(x_diam.', All_mma_m_round(1+i,1:end)', 'gauss3');
    coeffvaluess = coeffvalues(f_all_i);
    All_fit_parameters = [All_fit_parameters; coeffvaluess];
    yfit = f_all_i(x_diam); 
    ydata = All_mma_m_round(1+i,1:end)'; 
    gof = goodnessOfFit(yfit, ydata, 'NRMSE'); 
    allgoodnessoffit = [allgoodnessoffit; gof];
    
    % Initialize counters
    x = 0;
    y = 0;
    z = 0;
    
    % Process coefficients < 40 and >= 0
    for j = [2, 5, 8]
        if coeffvaluess(1, j) < 40 && coeffvaluess(1, j) >= 0
            x = x + 1;
            if x == 1
                coefficientall_nuc = [coefficientall_nuc; coeffvaluess(1, j)];
                disp("Adding to coefficientall_nuc: " + coeffvaluess(1, j))
            elseif x == 2
                coefficientall_nuc(end) = (coefficientall_nuc(end) + coeffvaluess(1, j)) / 2;
                disp("Updating coefficientall_nuc: " + coefficientall_nuc(end))
            elseif x == 3
                coefficientall_nuc(end) = (coefficientall_nuc(end) * 2 + coeffvaluess(1, j)) / 3;
                disp("Updating coefficientall_nuc: " + coefficientall_nuc(end))
            end
        elseif coeffvaluess(1, j) < 100 && coeffvaluess(1, j) >= 0
            y = y + 1;
            if y == 1
                coefficientall_aitk = [coefficientall_aitk; coeffvaluess(1, j)];
                disp("Adding to coefficientall_aitk: " + coeffvaluess(1, j))
            elseif y == 2
                coefficientall_aitk(end) = (coefficientall_aitk(end) + coeffvaluess(1, j)) / 2;
                disp("Updating coefficientall_aitk: " + coefficientall_aitk(end))
            elseif y == 3
                coefficientall_aitk(end) = (coefficientall_aitk(end) * 2 + coeffvaluess(1, j)) / 3;
                disp("Updating coefficientall_aitk: " + coefficientall_aitk(end))
            end
        elseif coeffvaluess(1, j) > 100 && coeffvaluess(1, j) <= 300
            z = z + 1;
            if z == 1
                coefficientall_accum = [coefficientall_accum; coeffvaluess(1, j)];
                disp("Adding to coefficientall_accum: " + coeffvaluess(1, j))
            elseif z == 2
                coefficientall_accum(end) = (coefficientall_accum(end) + coeffvaluess(1, j)) / 2;
                disp("Updating coefficientall_accum: " + coefficientall_accum(end))
            elseif z == 3
                coefficientall_accum(end) = (coefficientall_accum(end) * 2 + coeffvaluess(1, j)) / 3;
                disp("Updating coefficientall_accum: " + coefficientall_accum(end))
            end
        end
    end

    if x == 0
        coefficientall_nuc = [coefficientall_nuc; 0];
        disp("Adding 0 to coefficientall_nuc")
    end
    if y == 0
        coefficientall_aitk = [coefficientall_aitk; 0];
        disp("Adding 0 to coefficientall_aitk")
    end
    if z == 0
        coefficientall_accum = [coefficientall_accum; 0];
        disp("Adding 0 to coefficientall_accum")
    end
    
    c3bi = sort([coeffvaluess(1,2), coeffvaluess(1,5), coeffvaluess(1,8)]);
    st3bi = sort([coeffvaluess(1,3), coeffvaluess(1,6), coeffvaluess(1,9)]);
    
    % Ensure disjoint intervals for modes
    mode1i = All_mma_m_round(i, All_mma_m_round(1,1:end) < c3bi(1,1) - abs(sqrt(st3bi(1,1)/2)));
    mode2i = All_mma_m_round(i, ...
        (All_mma_m_round(1,1:end) > c3bi(1,2) - abs(sqrt(st3bi(1,2)/2))) & ...
        (All_mma_m_round(1,1:end) < c3bi(1,2) + abs(sqrt(st3bi(1,2)/2))));
    mode3i = All_mma_m_round(i, ...
        (All_mma_m_round(1,1:end) > c3bi(1,3) - abs(sqrt(st3bi(1,3)/2))) & ...
        (All_mma_m_round(1,1:end) < c3bi(1,3) + abs(sqrt(st3bi(1,3)/2))));
    
    Tot_numi = sum(All_mma_m_round(i,1:end), 'all'); % Total number for the row
    All_total_hour_number = [All_total_hour_number; Tot_numi];
    
    % Mode 1
    if ~isempty(mode1i)
        sum_mode1i = mean(mode1i); % Take the average of values
        number_mode1i = sum_mode1i / Tot_numi;
        mode1 = [mode1; number_mode1i];
    else
        mode1 = [mode1; 0];
    end
    
    % Mode 2
    if ~isempty(mode2i)
        sum_mode2i = mean(mode2i); % Take the average of values
        number_mode2i = sum_mode2i / Tot_numi;
        mode2 = [mode2; number_mode2i];
    else
        mode2 = [mode2; 0];
    end
    
    % Mode 3
    if ~isempty(mode3i)
        sum_mode3i = mean(mode3i); % Take the average of values
        number_mode3i = sum_mode3i / Tot_numi;
        mode3 = [mode3; number_mode3i];
    else
        mode3 = [mode3; 0];
    end
end

% Final rounded results
Nuclei_geom_mean = round(coefficientall_nuc);
Aitken_geom_mean = round(coefficientall_aitk);
Accumula_geom_mean = round(coefficientall_accum);

Nuclei_persentage = round(mode1, 2);
Aitken_persentage = round(mode2, 2);
Accumula_persentage = round(mode3, 2);
All_fit_parameters=[All_fit_parameters,All_total_hour_number,mode1,mode2,mode3,allgoodnessoffit];
All_fit_parameters=[title;All_fit_parameters];
%%%Geometric mean for all the hourly averaged data 
M_A=mean(mode1(mode1 ~= 0),'omitnan');
M_C=mean(mode2(mode2 ~= 0),'omitnan');
M_E=mean(mode3(mode3 ~= 0),'omitnan');
M_B=mean(mode1_cloud(mode1_cloud~= 0) ,"omitnan");
M_D=mean(mode2_cloud(mode2_cloud~= 0) ,"omitnan");
M_F=mean(mode3_cloud(mode3_cloud~= 0) ,"omitnan");
G=mean(mode1,'omitnan');
H=mean(mode1_cloud,'omitnan');
I=mode1.*All_total_hour_number;
J=mode1_cloud.*cloud_total_hour_number;
cloud_time=linspace(0,189,1);
all_time=linspace(0,1001,1);
%figure()
%surf(coefficientall_nuc_cloud,coefficientall_aitk_cloud,coefficientall_accum_cloud);

figure()

plot(mode1_cloud);
figure()
plot(mode1);

%%%Statistical Significance

%%1. Έλεγχος t-Student



% Δεδομένα
data1 = All_mma_m_round(~isnan(All_mma_m_round) & All_mma_m_round>0 );
data2 = cloud_mma_m_round(~isnan(cloud_mma_m_round) & cloud_mma_m_round>0 );

%%% t-Student
[h, p] = ttest2(data1, data2);

% Αποτελέσματα
if h == 0
    disp('No statistical differnce between the data.');
else
    disp('There is statistical difference between the data.');
end
disp(['p-value: ', num2str(p)]);







%%% Χ^2 (Chi-Square)



observed = All_mma_m_round(~isnan(All_mma_m_round) & All_mma_m_round>0 );
expected = cloud_mma_m_round(~isnan(cloud_mma_m_round) & cloud_mma_m_round>0 );

if ~isvector(observed) || ~isreal(observed)
    error('observed must be a vector of real values.');
end
if ~isvector(expected) || ~isreal(expected)
    error('expected must be a vector of real values.');
end
% if length(observed) != length(expected)
%     error('observed and expected must have the same length.');
% end
if any(expected < 0)
    error('The ''expected'' values must be non-negative.');
end

% Χ^2 Goodness-of-fit test
[h, p] = chi2gof(observed, 'Expected', expected);


% Χ^2 Goodness-of-fit test
[h, p] = chi2gof(observed, 'Expected', expected);

if h == 0
    disp('No statistical difference.');
else
    disp('There is statistical difference.');
end
disp(['p-value: ', num2str(p)]);

%%%Anoval

% group1 = All_mma_m_round(~isnan(All_mma_m_round) & All_mma_m_round>0 );
% group2 = cloud_mma_m_round(~isnan(cloud_mma_m_round) & cloud_mma_m_round>0 );
% %group3 = [1.5, 2.7, 3.8, 4.9, 5.2];
% 
% 
% data = [group1, group2];
% groups = [ones(size(group1)), 2*ones(size(group2))];
% 
% %  ANOVA
% [p, tbl, stats] = anova1(data, groups);
% 
% 
% if p > 0.05
%     disp('No statistical difference .');
% else
%     disp('There is statistical difference.');
% end
% disp(['p-value: ', num2str(p)]);




