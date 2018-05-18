function [VCDs,profiles] = get_GEOS_CHEM_summary_file()
% after running "get_GEOS_CHEM_all", one can use this function to combine
% outputs into single file
output_file_path = 'C:\Projects\GEOS_CHEM\data\reformat\';
summary_file_output_file_path = 'C:\Projects\GEOS_CHEM\data\summary_files\';

list_VCD = ls([output_file_path 'VCD*Downsview*']);
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
writetable(VCDs,[summary_file_output_file_path 'VCD_Downsview.csv']);

list_profile = ls([output_file_path 'profile*Downsview*']);
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
writetable(profiles,[summary_file_output_file_path 'profiles_Downsview.csv']);