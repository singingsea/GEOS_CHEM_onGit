function get_GEOS_CHEM_all_v2()
% this is the "main" function, make sure fill out the folder paths
% location of interest
% this 2nd version read in Tailong's latest data, which has air density and
% box-height info integarted.
% after run this function, use get_GEOS_CHEM_summary_file to generate final
% output file and plots

% Xiaoyi --- 2018/10/12

site = 'Downsview';
%site = 'Egbert'
%site = 'FortMcKay'

% default lon/lat information for Pandora site
if strcmp(site,'Downsview')
    user_lat=43.7810; % Downsview
    user_lon=-79.4680;
elseif strcmp(site,'Egbert')
    user_lat=44.2300; 
    user_lon=-79.7800;
elseif strcmp(site,'FortMcKay')
    user_lat=57.1836;
    user_lon=-111.6400;
end

% %trace_gas_file_path = 'C:\Projects\GEOS_CHEM\data\nc_ts\201501\';% this version has almost no spin up!!!
% trace_gas_file_path = 'C:\Projects\GEOS_CHEM\data\nc_ts\01_new_new\';% this version has spin up for 6 month
% support_file_path = 'C:\Projects\GEOS_CHEM\data\support\';
data_file_path = 'E:\Projects\GEOS-Chem\data\GEOS-CHEM\';
%output_file_path = 'C:\Projects\GEOS_CHEM\data\reformat\';
output_file_path = 'E:\Projects\GEOS-Chem\output\';
list = ls(data_file_path);
N = size(list);
for i_file = 3:N(1)
    filename = list(i_file,:);
    p_finished = ((i_file-2)/N(1))*100;
    disp(['Extracting file: ' filename ' @' site '  --' num2str(p_finished) '%-- ']);
    f_nm = [data_file_path filename];
    %support_file_date = filename(3:10); % get the date stamp from supporting file name
    %f_nm_support = [support_file_path 'T_PS_' support_file_date '.nc'];
    %f_nm_support_BLH = [support_file_path 'TROP_PBLH_' support_file_date '.nc'];
    %f_nm_airdensity = ['C:\Projects\GEOS_CHEM\data\air_density\' 'AD_' support_file_date '.nc'];
    
%     % if these tables are not exist, then we need read them from files for sure
%     support_info_missing = ~exist('PS_a_24hr') || ~exist('T_a47_24hr') || ~exist('PBLH_a') || ~exist('TROPPT_a');
%     if  support_info_missing | strcmp(filename(12:13),'00') % only read support data if these are missing or it's the first file of the day
%         try
%             disp(['Read supporting data from : ' 'T_PS_' support_file_date '.nc']);
%             % read in surface pressure and T profiles (all 24 hours)
%             [PS_a_24hr,T_a47_24hr] = read_GEOS_CHEM_meteo(f_nm_support,user_lat,user_lon);
%             
%             disp(['Read supporting data from : ' 'TROP_PBLH_' support_file_date '.nc']);
%             % read in boundary layer height and Tropause pressure (all 24 hours)
%             [PBLH_a,TROPPT_a] = read_GEOS_CHEM_TROP_PBLH(f_nm_support_BLH,user_lat,user_lon);
%         catch
%             disp(['Warnning: failed to extracting suportting file for ' filename]);
%         end
%     end
    
    try
        %[profile_1hr,VCD_1hr] = read_GEOS_CHEM(f_nm,f_nm_support,f_nm_support_BLH,site);
        [profile_1hr,VCD_1hr] = read_GEOS_CHEM_v2(f_nm,user_lat,user_lon);
        
%         %[BXHGHTS_BXHEIGHT_a,BXHGHTS_AIRNUMDE_a] = read_GEOS_CHEM_AirDensity(f_nm_airdensity,user_lat,user_lon);
%         profile_1hr.BXHGHTS_BXHEIGHT_a = BXHGHTS_BXHEIGHT_a';
%         profile_1hr.BXHGHTS_AIRNUMDE_a = BXHGHTS_AIRNUMDE_a';
%         DU = 2.6870e+20;%[molec/m^2]
%         VCD_1hr.o3_test = sum(profile_1hr.BXHGHTS_BXHEIGHT_a.*profile_1hr.BXHGHTS_AIRNUMDE_a.*profile_1hr.O3_vmr*1e-9)/DU;
%         VCD_1hr.no2_test = sum(profile_1hr.BXHGHTS_BXHEIGHT_a.*profile_1hr.BXHGHTS_AIRNUMDE_a.*profile_1hr.NO2_vmr*1e-9)/DU;
        % export data to csv files
        writetable(VCD_1hr,[output_file_path 'VCD_' filename(3:17) '_' site '.csv']);
        writetable(profile_1hr,[output_file_path 'profile_' filename(3:17) '_' site '.csv']);
    catch
        disp(['Warnning: failed to extracting file from ' filename]);
    end
    
end

