function get_GEOS_CHEM_all_v5()
% this is the "main" function, make sure fill out the folder paths
% location of interest
% this 3rd version read in Tailong's latest data, which has air density and
% box-height info integarted. 
% To create convertion LUT, please follow the steps:
% 1.run this function: get_GEOS_CHEM_all_v5()
% 2. run get_GEOS_CHEM_summary_file ---> to generate final output file and plots
% 3. run monthly_table = get_GEOS_CHEM_LUT() ---> create LUT

% Xiaoyi --- 2019/03/31

start_date = '20150102';% this is the start date of files to be processed

site = 'Downsview';
%site = 'Egbert'
%site = 'FortMcKay';
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

if ispc
    general_data_file_path = ['C:\Projects\GEOS_CHEM\data\MERRA2_05x0625_47L_NA_V12\'];%MERRA2_05x0625_47L_NA_V12 %GEOSFP_2x25_47L_V12
    output_file_path = 'C:\Projects\GEOS_CHEM\output\temp\MERRA2_05x0625_47L_NA_V12\';
    mkdir(output_file_path);
else    
    %general_data_file_path = ['/net/aurora/model_datasets/geoschem/2x25/' year '/'];% data path on Aurora
    general_data_file_path = ['/net/aurora/model_datasets/geoschem/MERRA2_05x0625_47L_NA_V12/'];%GEOSFP_2x25_47L_V12
    %month_list = ls(general_data_file_path);% this only works for Windows! 
    %month_list = dir(general_data_file_path);% this works for both Windows and UNIX! 
    output_file_path = '/export/data/home/xizhao/GEOS_CHEM/output/temp/MERRA2_05x0625_47L_NA_V12/';
    
    mkdir(output_file_path);
end



if ispc
    data_file_path = [general_data_file_path 'chem\'];
    met_file_path = [general_data_file_path 'met\'];
else
    %data_file_path = [general_data_file_path month_list(i_month+2).name '/'];% this works for both Windows and UNIX! 
    data_file_path = [general_data_file_path 'chem/'];
    met_file_path = [general_data_file_path 'met/'];
end
% get list of files in that month
%list = ls(data_file_path);% this only works for Windows! 
list = dir(data_file_path);% this works for both Windows and UNIX! 
list(1:2) = [];% remove the first two empty entries
filenames = cell2mat({list.name}');% get all file names into mat char
filedates = filenames(:,numel('GEOSChem.SpeciesConc.')+1:end - numel('_0000z.nc4'));% convert file names to strings
tf_process = datetime(filedates,'InputFormat','yyyyMMdd') >= datetime(start_date,'InputFormat','yyyyMMdd');% get tf for processing
list = list(tf_process);% update file list that will be processed

N = size(list);
for i_file = 1:N(1)
    %filename = list(i_file,:);% this only works for Windows! 
    filename = list(i_file).name;% this works for both Windows and UNIX! 
    p_finished = ((i_file)/N(1))*100;
    disp(['Extracting file: ' filename ' @' site '  --' num2str(p_finished) '%-- ']);
    f_nm = [data_file_path filename];

    if i_file == 1
        LON = ncread(f_nm,'lon');
        LAT = ncread(f_nm,'lat');
        % find the index of the profiles at given location
        [user_lon_index,user_lat_index,d_min] = find_profiles_at_location(LON,LAT,user_lat,user_lon);
    end
    disp(['Find closest grid at: Lon = ' num2str(LON(user_lon_index)) ' ; Lat = ' num2str(LAT(user_lat_index))]);
    disp(['Distance from site = ' num2str(d_min/1000) 'km']);

    % read profiles and VCDs for a whole day
    [profiles_1day,VCDs_1day] = read_GEOS_CHEM_v5(f_nm,met_file_path,user_lat_index,user_lon_index,output_file_path,site);


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

