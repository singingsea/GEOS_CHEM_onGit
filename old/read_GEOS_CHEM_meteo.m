function [PS_a_24hr,T_a47_24hr] = read_GEOS_CHEM_meteo(f_nm,user_lat,user_lon)
%f_nm = 'C:\Projects\GEOS_CHEM\data\support\T_PS_20150101.nc';
%f_nm2 = 'C:\Projects\GEOS_CHEM\data\support\TROP_PBLH_20150101.nc';

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
%info2 = ncinfo(f_nm2);
PS = ncread(f_nm,'PS');% the surface pressure for "current" day (see time stamp of the file name)
T = ncread(f_nm,'T');% the T profile for "current" day
current_timestamp = f_nm(end-10:end-3);% get time stamp "yyyymmdd" from file name
previous_timestamp = datestr((datetime(datevec(current_timestamp,'yyyymmdd'))-days(1)),'yyyymmdd');% convert it to previous day's time stamp
f_nm_previousday = [f_nm(1:end-11) previous_timestamp '.nc'];
next_timestamp = datestr((datetime(datevec(current_timestamp,'yyyymmdd'))+days(1)),'yyyymmdd');% convert it to next day's time stamp
f_nm_next = [f_nm(1:end-11) next_timestamp '.nc'];
try 
    PS_previousday = ncread(f_nm_previousday,'PS');% the surface pressure for "previous" day (see time stamp of the file name)
    T_previousday = ncread(f_nm_previousday,'T');% the T profile for "previous" day
    previousday_found = true;
catch
    disp('Warnning: failed to find surface P and T profile from previous day!');
    previousday_found = false;
end
try 
    PS_nextday = ncread(f_nm_next,'PS');% the surface pressure for "previous" day (see time stamp of the file name)
    T_nextday = ncread(f_nm_next,'T');% the T profile for "previous" day
    next_found = true;
catch
    disp('Warnning: failed to find surface P and T profile from next day!');
    next_found = false;
end

lev = ncread(f_nm,'lev');
lat = ncread(f_nm,'lat');
lon = ncread(f_nm,'lon');
time = ncread(f_nm,'time');

% find the index of the profiles at given location
[i,j] = find_profiles_at_location(lon,lat,user_lat,user_lon);

if (previousday_found == false) | (next_found == false)
    % get the profile over the site
    PS_a = reshape(PS(i,j,:),1,[]);
    T_a = reshape(T(i,j,:,:),72,[]);

    i_72 = 1;
    for i_47 = 1:47
        if i_47 <= 36 % for level below 36, GEOS-chem use same grid as reanalysis
            T_a47(i_47,:) = T_a(i_72,:);
            i_72 = i_72 + 1;
        elseif (i_47 >= 37) & (i_47 <= 40) % for level between 37 and 40, GEOS-chem average every two level into one
            T_a47(i_47,:) = (T_a(i_72,:) + T_a(i_72 + 1,:))/2;
            i_72 = i_72 + 2;
        elseif (i_47 >= 41) % for level from 41 to 47, GEOS-chem average every four level into one
            T_a47(i_47,:) = (T_a(i_72,:) + T_a(i_72 + 1,:) + T_a(i_72 + 2,:) + T_a(i_72 + 3,:))/4;
            i_72 = i_72 + 4;    
        end
    end

    % the surface pressure and T profile has temeprol res as 3-hr, so need to
    % be interp to 1-hr
    PS_a_24hr = interp1(1:3:24,PS_a,1:1:24);
    for i =1:47
        T_a47_24hr(i,:) = interp1(1:3:24,T_a47(i,:),1:1:24);
    end
    PS_a_24hr(:,23) = PS_a_24hr(:,22);
    PS_a_24hr(:,24) = PS_a_24hr(:,22);
    T_a47_24hr(:,23) = T_a47_24hr(:,22);
    T_a47_24hr(:,24) = T_a47_24hr(:,22);
    
elseif (previousday_found == true) & (next_found == true) % if we have information for both previous day and next day
    % get the profile over the site
    PS_a = reshape(PS(i,j,:),1,[]);
    T_a = reshape(T(i,j,:,:),72,[]);
    PS_previousday_a = reshape(PS_previousday(i,j,:),1,[]);
    T_previousday_a = reshape(T_previousday(i,j,:,:),72,[]);
    PS_nextday_a = reshape(PS_nextday(i,j,:),1,[]);
    T_nextday_a = reshape(T_nextday(i,j,:,:),72,[]);
    
    i_72 = 1;
    for i_47 = 1:47
        if i_47 <= 36 % for level below 36, GEOS-chem use same grid as reanalysis
            T_a47(i_47,:) = T_a(i_72,:);
            T_previousday_a47(i_47,:) = T_previousday_a(i_72,:);
            T_nextday_a47(i_47,:) = T_nextday_a(i_72,:);
            i_72 = i_72 + 1;
        elseif (i_47 >= 37) & (i_47 <= 40) % for level between 37 and 40, GEOS-chem average every two level into one
            T_a47(i_47,:) = (T_a(i_72,:) + T_a(i_72 + 1,:))/2;
            T_previousday_a47(i_47,:) = (T_previousday_a(i_72,:) + T_previousday_a(i_72 + 1,:))/2;
            T_nextday_a47(i_47,:) = (T_nextday_a(i_72,:) + T_nextday_a(i_72 + 1,:))/2;
            i_72 = i_72 + 2;
        elseif (i_47 >= 41) % for level from 41 to 47, GEOS-chem average every four level into one
            T_a47(i_47,:) = (T_a(i_72,:) + T_a(i_72 + 1,:) + T_a(i_72 + 2,:) + T_a(i_72 + 3,:))/4;
            T_previousday_a47(i_47,:) = (T_previousday_a(i_72,:) + T_previousday_a(i_72 + 1,:) + T_previousday_a(i_72 + 2,:) + T_previousday_a(i_72 + 3,:))/4;
            T_nextday_a47(i_47,:) = (T_nextday_a(i_72,:) + T_nextday_a(i_72 + 1,:) + T_nextday_a(i_72 + 2,:) + T_nextday_a(i_72 + 3,:))/4;
            i_72 = i_72 + 4;    
        end
    end

    % the surface pressure and T profile has temeprol res as 3-hr, so need to
    % be interp to 1-hr
    PS_a_4interp = [PS_previousday_a(:,8), PS_a,PS_nextday_a(:,1)];
    PS_a_30min_res = interp1(1:1:10,PS_a_4interp,1:1/6:10);
    PS_a_24hr = PS_a_30min_res(4:2:50);% from 00:00, 01:00, ..., 22:00, 23:00
    T_a47_4interp = [T_previousday_a47(:,8), T_a47,T_nextday_a47(:,1)];
    for i =1:47
        T_a47_30min_res(i,:) = interp1(1:1:10,T_a47_4interp(i,:),1:1/6:10);
        T_a47_24hr(i,:) = T_a47_30min_res(i,4:2:50);
    end

end

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
