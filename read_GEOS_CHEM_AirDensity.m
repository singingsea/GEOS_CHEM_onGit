function [BXHGHTS_BXHEIGHT_a,BXHGHTS_AIRNUMDE_a] = read_GEOS_CHEM_AirDensity(f_nm,user_lat,user_lon)
%function [PBLH_a,TROPPT_a] = read_GEOS_CHEM_TROP_PBLH()
%f_nm = 'C:\Projects\GEOS_CHEM\data\support\TROP_PBLH_20150101.nc';

% location of interest
%site = 'Downsview';
%site = 'Egbert'
%site = 'FortMcKay'

% default lon/lat information for Pandora site
% if strcmp(site,'Downsview')
%     user_lat=43.7810; % Downsview
%     user_lon=-79.4680;
% elseif strcmp(site,'Egbert')
%     user_lat=44.2300; 
%     user_lon=-79.7800;
% elseif strcmp(site,'FortMcKay')
%     user_lat=57.1836;
%     user_lon=-111.6400;
% end



%info = ncinfo(f_nm);

%PBLH = ncread(f_nm,'PBLH');% 'Planetary boundary layer height above surface' [m]
%TROPPT = ncread(f_nm,'TROPPT');%'Temperature-based tropopause pressure' [hPa]
lat = ncread(f_nm,'LAT');
lon = ncread(f_nm,'LON');
BXHGHTS_BXHEIGHT = ncread(f_nm,'BXHGHT-S__BXHEIGHT');
BXHGHTS_AIRNUMDE = ncread(f_nm,'BXHGHT-S__AIRNUMDE');

% find the index of the profiles at given location
[i,j] = find_profiles_at_location(lon,lat,user_lat,user_lon);
% get the profile over the site
BXHGHTS_BXHEIGHT_a = reshape(BXHGHTS_BXHEIGHT(i,j,:),1,[]);
BXHGHTS_AIRNUMDE_a = reshape(BXHGHTS_AIRNUMDE(i,j,:),1,[]);

%%
function [i,j] = find_profiles_at_location(LON,LAT,user_lat,user_lon)
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
        else
            start_searching = false;
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
