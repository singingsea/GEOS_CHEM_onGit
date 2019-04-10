function find_GEOS_CHEM_mean_profile
% this function can find and plot all GEOS-Chem profiles
load('C:\Projects\GEOS_CHEM\output\summary_files_test2\VCDs_Profile_Downsview.mat');

figure;hold all;
time_unique = unique(profiles.UTC);
profiles_P = profiles.P;
profile_P = reshape(profiles_P,[47, 411720/47]);
profiles_vmr = profiles.NO2_vmr;
profile_vmr = reshape(profiles_vmr,[47, 411720/47]);


% for i =1:length(time_unique)    
%     TF_single_profile = profiles.UTC == time_unique(i,:);
%     profile = profiles(TF_single_profile,:);
%     %plot(profile.NO2_vmr,profile.P,'.-');
%     mean_profile_vmr(:,i) = profile.NO2_vmr;    
%     mean_profile_P(:,i) = profile.P;  
% end

mean_vmr = mean(profile_vmr,2);
mean_P = mean(profile_P,2);
figure;
plot(mean_vmr,mean_P,'-');
ylim([600 1010]);
ylabel(['hPa']);
xlabel(['ppbv']);
    
    

