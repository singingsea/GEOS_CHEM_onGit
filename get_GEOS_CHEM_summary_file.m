function [VCDs,profiles] = get_GEOS_CHEM_summary_file()
% after running "get_GEOS_CHEM_all", one can use this function to combine
% outputs into single file
if ispc
    addpath('C:\Users\ZhaoX\Documents\MATLAB\matlab');
    output_file_path = 'C:\Projects\GEOS_CHEM\output\temp\';
    summary_file_output_file_path = 'C:\Projects\GEOS_CHEM\output\summary_files\';
else
    addpath('/export/data/home/xizhao/matlab/');
    output_file_path = '/export/data/home/xizhao/GEOS_CHEM/output/temp/';
    summary_file_output_file_path = '/export/data/home/xizhao/GEOS_CHEM/output/summary_files/';
end
mkdir(summary_file_output_file_path);

%site = 'Downsview';
%site = 'Egbert'
site = 'FortMcKay'
%site = 'Beijing'
%site = 'LosAngeles'
plot_path = ['C:\Projects\GEOS_CHEM\plots\' site '\'];
save_fig = 1;


%list_VCD = ls([output_file_path 'VCD*Downsview*']);
list_VCD = ls([output_file_path 'VCD*' site '*']);
N = size(list_VCD);
for i = 1:N(1)
    VCD = readtable([output_file_path list_VCD(i,:)]);
    if i == 1
        VCDs = VCD;
    else
        VCDs = [VCDs;VCD];
    end
    p_finished = (i/N(1))*100;
    disp(['Extracting VCD files:  @' site '  --' num2str(p_finished) '%-- ']);
end
VCDs = sortrows(VCDs,'UTC');
writetable(VCDs,[summary_file_output_file_path 'VCD_' site '.csv']);

list_profile = ls([output_file_path 'profile*' site '*']);
N = size(list_profile);
for i = 1:N(1)
    profile = readtable([output_file_path list_profile(i,:)]);
    if i == 1
        profiles = profile;
    else
        profiles = [profiles;profile];
    end
    p_finished = (i/N(1))*100;
    disp(['Extracting profile files:  @' site '  --' num2str(p_finished) '%-- ']);
end
profiles = sortrows(profiles,'UTC');
writetable(profiles,[summary_file_output_file_path 'profiles_' site '.csv']);

matfile_nm = [summary_file_output_file_path 'VCDs_Profile_' site '.mat'];
save(matfile_nm,'profiles','VCDs');

fig_name = 'ozone_VCD';
figure;hold all;
plot(VCDs.UTC,VCDs.o3,'.');
xlabel('UTC');
ylabel('Ozone VCD [DU]');
print_setting(1/4,save_fig,[plot_path fig_name]);

fig_name = 'ozone_surface';
figure;hold all;
plot(VCDs.UTC,VCDs.o3_surf,'.');
xlabel('UTC');
ylabel('O_3 surface [ppbv]');
print_setting(1/4,save_fig,[plot_path fig_name]);

fig_name = 'ozone_CVratio';
figure;hold all;
plot(VCDs.UTC,VCDs.o3_CV,'-.');
xlabel('UTC');
ylabel('O_3 surface/VCD [ppbv/DU]');
print_setting(1/4,save_fig,[plot_path fig_name]);

fig_name = 'NO2_VCD';
figure;hold all;
plot(VCDs.UTC,VCDs.no2,'.');
xlabel('UTC');
ylabel('NO_2 VCD [DU]');
print_setting(1/4,save_fig,[plot_path fig_name]);

fig_name = 'NO2_surface';
figure;hold all;
plot(VCDs.UTC,VCDs.no2_surf,'.');
xlabel('UTC');
ylabel('NO_2 surface [ppbv]');
print_setting(1/4,save_fig,[plot_path fig_name]);

fig_name = 'NO2_CVratio';
figure;hold all;
plot(VCDs.UTC,VCDs.no2_CV,'-.');
xlabel('UTC');
ylabel('NO_2 surface/VCD [ppbv/DU]');
print_setting(1/4,save_fig,[plot_path fig_name]);