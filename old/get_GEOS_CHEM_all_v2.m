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
%site = 'Beijing'
%site = 'LosAngeles'


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
elseif strcmp(site,'Beijing')
    user_lat=39.9;
    user_lon=116.4;    
elseif strcmp(site,'LosAngeles')
    user_lat=34.03;
    user_lon=-118.23;    
end

general_data_file_path = 'C:\Projects\GEOS_CHEM\data\2014\';
month_list = ls(general_data_file_path);
output_file_path = 'C:\Projects\GEOS_CHEM\output\';

for i_month = 1:12
    data_file_path = [general_data_file_path month_list(i_month+3,:) '\'];
    


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
        if i_file == 3
            %LON = ncread(f_nm,'LON');
            LON = ncread(f_nm,'lon');
            %LAT = ncread(f_nm,'LAT');
            LAT = ncread(f_nm,'lat');
            % find the index of the profiles at given location
            [user_lat_index,user_lon_index,d_min] = find_profiles_at_location(LON,LAT,user_lat,user_lon);
        end
        disp(['Find closest grid at: Lon = ' num2str(LON(user_lat_index)) ' ; Lat = ' num2str(LAT(user_lon_index))]);
        disp(['Distance from site = ' num2str(d_min/1000) 'km']);
        
        %[profile_1hr,VCD_1hr] = read_GEOS_CHEM(f_nm,f_nm_support,f_nm_support_BLH,site);
        [profile_1hr,VCD_1hr] = read_GEOS_CHEM_v2(f_nm,user_lat_index,user_lon_index);
        
%         %[BXHGHTS_BXHEIGHT_a,BXHGHTS_AIRNUMDE_a] = read_GEOS_CHEM_AirDensity(f_nm_airdensity,user_lat,user_lon);
%         profile_1hr.BXHGHTS_BXHEIGHT_a = BXHGHTS_BXHEIGHT_a';
%         profile_1hr.BXHGHTS_AIRNUMDE_a = BXHGHTS_AIRNUMDE_a';
%         DU = 2.6870e+20;%[molec/m^2]
%         VCD_1hr.o3_test = sum(profile_1hr.BXHGHTS_BXHEIGHT_a.*profile_1hr.BXHGHTS_AIRNUMDE_a.*profile_1hr.O3_vmr*1e-9)/DU;
%         VCD_1hr.no2_test = sum(profile_1hr.BXHGHTS_BXHEIGHT_a.*profile_1hr.BXHGHTS_AIRNUMDE_a.*profile_1hr.NO2_vmr*1e-9)/DU;
        % export data to csv files
        %writetable(VCD_1hr,[output_file_path 'VCD_' filename(3:17) '_' site '.csv']);
        %writetable(profile_1hr,[output_file_path 'profile_' filename(3:17) '_' site '.csv']);
        k = strfind(f_nm,'\ts');
        writetable(VCD_1hr,[output_file_path 'VCD_' f_nm(k+3:end-3) '_' site '.csv']);
        writetable(profile_1hr,[output_file_path 'profile_' f_nm(k+3:end-3) '_' site '.csv']);
    catch
        disp(['Warnning: failed to extracting file from ' filename]);
    end
    
end

function [i_min,j_min,d_min] = find_profiles_at_location(LON,LAT,user_lat,user_lon)
% find the cloest profile for the site, output the index for lon and lat
LON_dim = size(LON);
LAT_dim = size(LAT);
start_searching = true;
for i=1:LON_dim(1)% dim of lon
    for j=1:LAT_dim(1)% dim of lat
        d = get_distance(user_lat,user_lon,LAT(j),LON(i));
        if start_searching == true
            d_min = d;
            i_min = i;
            j_min = j;
            start_searching = false;
        else
            if (d_min > d)
                d_min = d;
                i_min = i;
                j_min = j;
            else
            end
        end                
    end
end

%%
function d = get_distance(user_lat,user_lon,lat,lon)
% sub function to calculate distance
R=6371000;%radius of the earth in meters
lat1=degtorad(user_lat);
lat2=degtorad(lat);
delta_lat=degtorad(lat-user_lat);
delta_lon=degtorad(lon-user_lon);
a=(sin(delta_lat/2))*(sin(delta_lat/2))+(cos(lat1))*(cos(lat2))*(sin(delta_lon/2))*(sin(delta_lon/2));
c=2.*asin(sqrt(a));
d=R*c;

