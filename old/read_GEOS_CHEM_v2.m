function [profile,VCD] = read_GEOS_CHEM_v2(f_nm,user_lat_index,user_lon_index)
% this function can read GEOS CHEM O3/NO2/SO2/HCHO data
% this function need inputs from "read_GEOS_CHEM_meteo" and "read_GEOS_CHEM_TROP_PBLH", which will provide
% supporting data (surface pressure, T profile, BLH, and tropopause
% height); these information are used to perform column integration and
% trop/strat columns seperation

% inputs:
    %f_nm --> file name (include path) that has trace gas vmr [ppbv] profile
    %PS_a_24hr,T_a47_24hr --> surface pressure [hPa] and temepature [K] profiles (interp
    %to 24 hr), these two data are generated from "read_GEOS_CHEM_meteo"
    %PBLH_a,TROPPT_a --> boundary layer height [m] and tropopause pressure
    %[hPa], these two data are generated from "read_GEOS_CHEM_TROP_PBLH"
    %user_lat,user_lon --> desired location to readout from GEOS-chem output
% outputs:
    % profile --> single day (1hr) trace gas profiles at the given location
    % VCD --> single day (1hr) integareted VCDs (total columns [DU], trop/start columns [DU], trop height [m])

% by Xiaoyi 2018/05/18

%f_nm = 'C:\Projects\GEOS_CHEM\data\nc_ts\ts20150120-170000.nc';
%f_nm_support = 'C:\Projects\GEOS_CHEM\data\support\T_PS_20150120.nc';



% % read in surface pressure and T profiles (all 24 hours)
% [PS_a_24hr,T_a47_24hr] = read_GEOS_CHEM_meteo(f_nm_support,user_lat,user_lon);
% % read in boundary layer height and Tropause pressure (all 24 hours)
% [PBLH_a,TROPPT_a] = read_GEOS_CHEM_TROP_PBLH(f_nm_support_BLH,user_lat,user_lon);
% hr = str2num(f_nm(end-8:end-7));
% P_surf = PS_a_24hr(:,hr+1); % surface pressure in this hour
% PBLH = PBLH_a(:,hr+1); % boundary layer height in this hour
% TROPPT = TROPPT_a(:,hr+1); % Tropause pressure in this hour

% read in file info and variables
info = ncinfo(f_nm);
LON = ncread(f_nm,'lon');
LAT = ncread(f_nm,'lat');
%Ap = ncread(f_nm,'Ap');
%Bp = ncread(f_nm,'Bp');
P = ncread(f_nm,'PEDGE_S_PSURF');
AIRDEN = ncread(f_nm,'TIME_SER_AIRDEN');%molec/cm^3
BXHEIGHT = ncread(f_nm,'BXHGHT_S_BXHEIGHT');% box height m
PBL_M = ncread(f_nm,'PBLDEPTH_PBL_M');% PBL m
PBL_L = ncread(f_nm,'PBLDEPTH_PBL_L');%PBL level

% read in trace gas 
NO2 = ncread(f_nm,'IJ_AVG_S_NO2');
O3 = ncread(f_nm,'IJ_AVG_S_O3');
HCHO = ncread(f_nm,'IJ_AVG_S_CH2O');
SO2 = ncread(f_nm,'IJ_AVG_S_SO2');
NO = ncread(f_nm,'IJ_AVG_S_NO');



% check the date format, and then generate a UTC timestamp
% if strcmp(info.Attributes(6).Name,'Start_Date')
%     date = info.Attributes(6).Value;
%     time = info.Attributes(7).Value;
%     if time == 0 
%         UTC = datevec([num2str(date) '/000000'],'yyyymmdd/HHMMSS');
%     else
%         UTC = datevec([num2str(date) '/' num2str(time)],'yyyymmdd/HHMMSS');
%     end
% else
%     disp("Warnning, the format of the nc file info might be different, please check Attributes to locate start time");
% end

k = strfind(f_nm,'\ts');
datetime_str = f_nm(k+3:end-3);
UTC = datevec([datetime_str '/000000'],'yyyymmdd/HHMMSS');

i = user_lat_index;
j = user_lon_index;
% get the profile over the site
O3_a = O3(i,j,:);% O3 VMR in [ppbv]
NO2_a = NO2(i,j,:);% NO2 VMR in [ppbv]
NO_a = NO(i,j,:);% NO VMR in [ppbv]
SO2_a = SO2(i,j,:);% SO2 VMR in [ppbv]
HCHO_a = HCHO(i,j,:);% HCHO VMR in [ppbv]
AIRDEN_a = AIRDEN(i,j,:);% air density 
BXHEIGHT_a = BXHEIGHT(i,j,:);% box height m
PBL_M_a = PBL_M(i,j); % PBL height [m]
PBL_L_a = PBL_L(i,j); % PBL height level
%Ap_a = Ap(:);% pressure info. that used to calculate Pressure at edge of box
%Bp_a = Bp(:);% pressure info. that used to calculate Pressure at edge of box
P_a = P(i,j,:);

O3_a = reshape(O3_a,[],1);% reshape to simple matrix
NO2_a = reshape(NO2_a,[],1);
NO_a = reshape(NO_a,[],1);
SO2_a = reshape(SO2_a,[],1);
HCHO_a = reshape(HCHO_a,[],1);
AIRDEN_a = reshape(AIRDEN_a,[],1);
BXHEIGHT_a = reshape(BXHEIGHT_a,[],1);
P_a = reshape(P_a,[],1);

% Ap_a = reshape(Ap_a,[],1);
% Bp_a = reshape(Bp_a,[],1);
% P_edge = Ap_a + Bp_a.*P_surf;% this is the pressure at the edge of each box
% T = T_a47_24hr(:,hr+1);% T profile in this hour
Rd = 287.0;
g0 = 9.8;
% for j = 1:numel(P_edge)-1
%     P(j) = (P_edge(j) + P_edge(j+1))/2;% this is the pressure at the centre of each box
%     box_height(j) = Rd/g0*T(j)*log(P_edge(j)/P_edge(j+1));% calculate the box height (note GEOS-CHEM has 47 boxes in vertical direction)
%     if j == 1
%         Z(j) = box_height(j)./2;% centre of the first box
%     else
%         Z(j) = Z(j-1)+box_height(j-1)./2 +box_height(j)./2;% centre of the rest ones
%     end
% end

% prepare output profile table
profile = table;
%profile.Z =  Z';
profile.P =  P_a;
profile.boxheight =  BXHEIGHT_a;
profile.airdensity =  AIRDEN_a;
%profile.pressure = P';
%profile.T = T;
profile.O3_vmr = O3_a;
profile.NO2_vmr = NO2_a;
profile.NO_vmr = NO_a;
profile.SO2_vmr = SO2_a;
profile.HCHO_vmr = HCHO_a;
%profile.box_height =  box_height';
profile.UTC = repmat(datetime(UTC),height(profile),1);

% calculate VCD
Aav = 6.02e23;
R = 8.314;
DU = 2.69e20;% molec/m^3 note here we use al SI units in the following equs

o3_vcd = sum(profile.airdensity.*1e6.*profile.O3_vmr.*1e-9.*profile.boxheight)/DU;
no2_vcd = sum(profile.airdensity.*1e6.*profile.NO2_vmr.*1e-9.*profile.boxheight)/DU;
no_vcd = sum(profile.airdensity.*1e6.*profile.NO_vmr.*1e-9.*profile.boxheight)/DU;
so2_vcd = sum(profile.airdensity.*1e6.*profile.SO2_vmr.*1e-9.*profile.boxheight)/DU;
hcho_vcd = sum(profile.airdensity.*1e6.*profile.HCHO_vmr.*1e-9.*profile.boxheight)/DU;


% TF = (profile.Z<=50e3); % just integrate VMR below 60 km.
% o3_vcd = sum(profile.pressure(TF,:).*100.*profile.O3_vmr(TF,:).*1e-9.*profile.box_height(TF,:).*Aav./(R.*profile.T(TF,:))./DU);
% no2_vcd = sum(profile.pressure(TF,:).*100.*profile.NO2_vmr(TF,:).*1e-9.*profile.box_height(TF,:).*Aav./(R.*profile.T(TF,:))./DU);
% so2_vcd = sum(profile.pressure(TF,:).*100.*profile.SO2_vmr(TF,:).*1e-9.*profile.box_height(TF,:).*Aav./(R.*profile.T(TF,:))./DU);
% hcho_vcd = sum(profile.pressure(TF,:).*100.*profile.HCHO_vmr(TF,:).*1e-9.*profile.box_height(TF,:).*Aav./(R.*profile.T(TF,:))./DU);
% airdensity = profile.pressure.*100.*Aav./(R.*profile.T);
% profile.airdensity = airdensity;

% o3_vcd = sum(profile.pressure.*100.*profile.O3_vmr.*1e-9.*profile.box_height.*Aav./(R.*profile.T)./DU);
% no2_vcd = sum(profile.pressure.*100.*profile.NO2_vmr.*1e-9.*profile.box_height.*Aav./(R.*profile.T)./DU);
% so2_vcd = sum(profile.pressure.*100.*profile.SO2_vmr.*1e-9.*profile.box_height.*Aav./(R.*profile.T)./DU);
% hcho_vcd = sum(profile.pressure.*100.*profile.HCHO_vmr.*1e-9.*profile.box_height.*Aav./(R.*profile.T)./DU);

% calculate VCDtrop
% TF = profile.pressure >= TROPPT;
% TROPZ = max(profile.Z(TF,:));% get tropopause height
% o3_vcdtrop = sum(profile.pressure(TF,:).*100.*profile.O3_vmr(TF,:).*1e-9.*profile.box_height(TF,:).*Aav./(R.*profile.T(TF,:))./DU);
% no2_vcdtrop = sum(profile.pressure(TF,:).*100.*profile.NO2_vmr(TF,:).*1e-9.*profile.box_height(TF,:).*Aav./(R.*profile.T(TF,:))./DU);
% so2_vcdtrop = sum(profile.pressure(TF,:).*100.*profile.SO2_vmr(TF,:).*1e-9.*profile.box_height(TF,:).*Aav./(R.*profile.T(TF,:))./DU);
% hcho_vcdtrop = sum(profile.pressure(TF,:).*100.*profile.HCHO_vmr(TF,:).*1e-9.*profile.box_height(TF,:).*Aav./(R.*profile.T(TF,:))./DU);

% calculate VCD_PBL 
%PBL_level = round(PBL_L_a);
PBL_level = fix(PBL_L_a);% this is the integer part of PBL level
PBL_level_decimal = PBL_L_a - fix(PBL_L_a);% this is the decimal part of PBL level
o3_vcd_pbl = sum(profile.airdensity(1:PBL_level).*1e6.*profile.O3_vmr(1:PBL_level).*1e-9.*profile.boxheight(1:PBL_level))/DU ...
    + PBL_level_decimal.*profile.airdensity(PBL_level+1).*1e6.*profile.O3_vmr(PBL_level+1).*1e-9.*profile.boxheight(PBL_level+1)/DU;
no2_vcd_pbl = sum(profile.airdensity(1:PBL_level).*1e6.*profile.NO2_vmr(1:PBL_level).*1e-9.*profile.boxheight(1:PBL_level))/DU ...
    + PBL_level_decimal.*profile.airdensity(PBL_level+1).*1e6.*profile.NO2_vmr(PBL_level+1).*1e-9.*profile.boxheight(PBL_level+1)/DU;
no_vcd_pbl = sum(profile.airdensity(1:PBL_level).*1e6.*profile.NO_vmr(1:PBL_level).*1e-9.*profile.boxheight(1:PBL_level))/DU ...
    + PBL_level_decimal.*profile.airdensity(PBL_level+1).*1e6.*profile.NO_vmr(PBL_level+1).*1e-9.*profile.boxheight(PBL_level+1)/DU;
so2_vcd_pbl = sum(profile.airdensity(1:PBL_level).*1e6.*profile.SO2_vmr(1:PBL_level).*1e-9.*profile.boxheight(1:PBL_level))/DU ...
    + PBL_level_decimal.*profile.airdensity(PBL_level+1).*1e6.*profile.SO2_vmr(PBL_level+1).*1e-9.*profile.boxheight(PBL_level+1)/DU;
hcho_vcd_pbl = sum(profile.airdensity(1:PBL_level).*1e6.*profile.HCHO_vmr(1:PBL_level).*1e-9.*profile.boxheight(1:PBL_level))/DU ...
    + PBL_level_decimal.*profile.airdensity(PBL_level+1).*1e6.*profile.HCHO_vmr(PBL_level+1).*1e-9.*profile.boxheight(PBL_level+1)/DU;

% % calculate VCDstart
% o3_vcdstart = o3_vcd - o3_vcdtrop;
% no2_vcdstart = no2_vcd - no2_vcdtrop;
% so2_vcdstart = so2_vcd - so2_vcdtrop;
% hcho_vcdstart = hcho_vcd - hcho_vcdtrop;

VCD = table;
VCD.UTC = datetime(UTC);
VCD.o3 = o3_vcd;% all VCDs are in DU
VCD.no2 = no2_vcd;
VCD.no = no_vcd;
VCD.so2 = so2_vcd;
VCD.hcho = hcho_vcd;
VCD.PBLH = PBL_M_a;% boundary layer height [m]
%VCD.TROPPT = TROPPT;% tropause pressure [hPa]
%VCD.TROPZ = TROPZ;% trop height [m]

VCD.o3_pbl = o3_vcd_pbl;% all tropospheric columns are in DU
VCD.no2_pbl = no2_vcd_pbl;
VCD.no_pbl = no_vcd_pbl;
VCD.so2_pbl = so2_vcd_pbl;
VCD.hcho_pbl = hcho_vcd_pbl;

VCD.o3_surf = profile.O3_vmr(1);% all tropospheric columns are in DU
VCD.no2_surf = profile.NO2_vmr(1);
VCD.no_surf = profile.NO_vmr(1);
VCD.so2_surf = profile.SO2_vmr(1);
VCD.hcho_surf = profile.HCHO_vmr(1);

VCD.o3_CV = VCD.o3_surf./VCD.o3_pbl;
VCD.no2_CV = VCD.no2_surf./VCD.no2_pbl;
VCD.no_CV = VCD.no_surf./VCD.no_pbl;
VCD.so2_CV = VCD.so2_surf./VCD.so2_pbl;
VCD.hcho_CV = VCD.hcho_surf./VCD.hcho_pbl;


% 
% VCD.o3_start = o3_vcdstart;% all startospheric columns are in DU
% VCD.no2_start = no2_vcdstart;
% VCD.so2_start = so2_vcdstart;
% VCD.hcho_start = hcho_vcdstart;






%%





