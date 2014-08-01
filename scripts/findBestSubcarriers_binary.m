% What this will do is look at the VARIABILITY and VARIANCE, and determine
% which subcarriers we will look at for our calculations of the GLRT.
% Therefore, the OUTPUT of this function will be nothing but an ARRAY OF
% ZEROS AND ONES.
%       example: return = [1 1 1 0 0 1 1 0 1 1]
%               ones indicate that that subcarrier is to be used
%               zeros indicate that that subcarrier is to be ignored



function good_ind = findBestSubcarriers_binary(scen1dat,scen2dat,min_carriers,Nsub)




good_ind = zeros(1,Nsub);
variability_pts = zeros(1,Nsub);
variance1_pts = zeros(1,Nsub);
variance2_pts = zeros(1,Nsub);


% Exploration of VARIABILITY: the further that each of the scenarios is
% from the other, the better the subcarrier will perform.

variability = zeros(1,Nsub);
variance1 = zeros(1,Nsub);
variance2 = zeros(1,Nsub);
for subnum = 1:Nsub
    scen1 = scen1dat(subnum,:);
    scen2 = scen2dat(subnum,:);
    
    variance1(subnum) = var(scen1);
    variance2(subnum) = var(scen2);
    
    means = [mean(scen1) mean(scen2)];
    
    % Calculate the difference between the maximum scenario and the minimum
    variability(subnum) = max(means)-min(means);
    
end; clear subnum means scen1 scen2 scen3;

sorted_variability = sort(variability);
sorted_variance1 = sort(variance1);
sorted_variance2 = sort(variance2);


% Go through each of the variability valies and assign an ordered number to each
% index in the array. The most-variable subcarrier will be assigned a value
% of 16, then 15, etc, etc for the M highest variabilities
for i = Nsub-min_carriers+1:Nsub
    variability_index = sorted_variability(i) == variability;
    variability_pts(variability_index) = i;
    
end; clear i variability_index

% Similar process as above, but points assigned to the least variance, and
% the max points = 16/3
for i = 1:min_carriers
    variance1_index = sorted_variance1(i) == variance1;
    variance1_pts(variance1_index) = (1/3)*(Nsub-i+1);
    
    variance2_index = sorted_variance2(i) == variance2;
    variance2_pts(variance2_index) = (1/3)*(Nsub-i+1);
    
end; clear i variance1_index variance2_index

%sum_pts = variance1_pts + variance2_pts + variance3_pts + variability_pts;
sum_pts = variance1_pts + variance2_pts + variability_pts;
sorted_pts = sort(sum_pts);

% Strict selection
for i = Nsub-min_carriers+1:Nsub
    good_index = sorted_pts(i) == sum_pts;
    good_ind(good_index) = 1;
end

% Lenient Selection
% mean_pts = mean(sum_pts);
% good_ind(sum_pts >= mean_pts) = 1;

end