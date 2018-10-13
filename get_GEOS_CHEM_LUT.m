function get_GEOS_CHEM_LUT()
% calculate LUT from GEOS-CHEM summary file

site = 'Downsview';timezone_offset = -5;
%site = 'Egbert'
%site = 'FortMcKay'

summary_file_output_file_path = 'E:\Projects\GEOS-Chem\output\summary_files\';
matfile_nm = [summary_file_output_file_path 'VCDs_Profile_' site '.mat'];
load(matfile_nm);

VCDs.LTC = VCDs.UTC + hours(timezone_offset);
VCDs(isnan(VCDs.no2_CV),:)=[];
VCDs(isinf(VCDs.no2_CV),:)=[];
LUT = table;
for i = 0:23
    TF = VCDs.LTC.Hour == i;
    LUT.LST_hour(i+1,:) = i;
    LUT.UTC_hour(i+1,:) = i-timezone_offset;
    LUT.no2_CV(i+1,:) = mean(VCDs.no2_CV(TF,:));
    LUT.no2_CV_err(i+1,:) = std(VCDs.no2_CV(TF,:))/(sum(TF))^0.5;
end

figure;hold all;
plot(LUT.LST_hour,LUT.no2_CV);
xlabel('LST [hour]');
ylabel('Converstion ratio [ppbv/DU]');
xlim([0 23]);
