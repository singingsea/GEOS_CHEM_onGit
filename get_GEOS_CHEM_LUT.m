function get_GEOS_CHEM_LUT()
% calculate LUT from GEOS-CHEM summary file

%site = 'Downsview';timezone_offset = -5;
%site = 'Egbert';timezone_offset = -5;
%site = 'FortMcKay';timezone_offset = -7;
%site = 'Beijing';timezone_offset = 8;
site = 'LosAngeles';timezone_offset = -7;

summary_file_output_file_path = 'E:\Projects\GEOS-Chem\output\summary_files\';
matfile_nm = [summary_file_output_file_path 'VCDs_Profile_' site '.mat'];
load(matfile_nm);

% test total VCD CV
VCDs.no2_CV2 = VCDs.no2_surf./VCDs.no2;

VCDs.LTC = VCDs.UTC + hours(timezone_offset);
VCDs(isnan(VCDs.no2_CV),:)=[];
VCDs(isinf(VCDs.no2_CV),:)=[];
TF_0 = VCDs.no2_CV == 0;
VCDs(TF_0,:)=[];
LUT = table;
for i = 0:23
    TF = VCDs.LTC.Hour == i;
    LUT.LST_hour(i+1,:) = i;
    LUT.UTC_hour(i+1,:) = i-timezone_offset;
    LUT.no2_CV(i+1,:) = mean(VCDs.no2_CV(TF,:));
    LUT.no2_CV_err(i+1,:) = std(VCDs.no2_CV(TF,:))/(sum(TF))^0.5;
    LUT.no2_CV2(i+1,:) = mean(VCDs.no2_CV2(TF,:));
    LUT.no2_CV2_err(i+1,:) = std(VCDs.no2_CV2(TF,:))/(sum(TF))^0.5;
end

figure;hold all;
plot(LUT.LST_hour,LUT.no2_CV);
plot(LUT.LST_hour,LUT.no2_CV2);
xlabel('LST [hour]');
ylabel('Converstion ratio [ppbv/DU]');
xlim([0 23]);
legend({'no_2 surface/VCD_P_B_L','no_2 surface/VCD'});
title(site);
