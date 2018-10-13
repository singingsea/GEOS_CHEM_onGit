function [VCDs,profiles] = get_GEOS_CHEM_summary_file()
% after running "get_GEOS_CHEM_all", one can use this function to combine
% outputs into single file

%output_file_path = 'C:\Projects\GEOS_CHEM\data\reformat\';
output_file_path = 'E:\Projects\GEOS-Chem\output\';
summary_file_output_file_path = 'E:\Projects\GEOS-Chem\output\summary_files\';

site = 'Downsview';
%site = 'Egbert'
%site = 'FortMcKay'


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
end
profiles = sortrows(profiles,'UTC');
writetable(profiles,[summary_file_output_file_path 'profiles_' site '.csv']);

matfile_nm = [summary_file_output_file_path 'VCDs_Profile_' site '.mat'];
save(matfile_nm,'profiles','VCDs');

figure;hold all;
plot(VCDs.UTC,VCDs.o3,'.');
xlabel('UTC');
ylabel('Ozone VCD [DU]');

figure;hold all;
plot(VCDs.UTC,VCDs.o3_surf,'.');
xlabel('UTC');
ylabel('O_3 surface [ppbv]');

figure;hold all;
plot(VCDs.UTC,VCDs.o3_CV,'-.');
xlabel('UTC');
ylabel('O_3 surface/VCD [ppbv/DU]');

figure;hold all;
plot(VCDs.UTC,VCDs.no2,'.');
xlabel('UTC');
ylabel('NO_2 VCD [DU]');

figure;hold all;
plot(VCDs.UTC,VCDs.no2_surf,'.');
xlabel('UTC');
ylabel('NO_2 surface [ppbv]');

figure;hold all;
plot(VCDs.UTC,VCDs.no2_CV,'-.');
xlabel('UTC');
ylabel('NO_2 surface/VCD [ppbv/DU]');