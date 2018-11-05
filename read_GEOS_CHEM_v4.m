function [profiles_1day,VCDs_1day] = read_GEOS_CHEM_v4(f_nm,user_lat_index,user_lon_index,output_file_path,site)
% this is a sub-function for get_GEOS_CHEM_all_v3 to read GEOS CHEM O3/NO2/SO2/HCHO data

% inputs:
    %f_nm --> file name (include path) that has trace gas vmr [ppbv] profile
    %user_lat_index,user_lon_index --> desired location index to readout from GEOS-chem output
% outputs:
    % profile --> single day (1hr) trace gas profiles at the given location
    % VCD --> single day (1hr) integareted VCDs (total columns [DU], trop/start columns [DU], trop height [m])

% by Xiaoyi 2018/11/05


% read in file info and variables
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


% each file has 24 hours of data, so need loop for each hour
for i_hour = 1:24
    if ispc
        k = strfind(f_nm,'\ts');
    else
        k = strfind(f_nm,'/ts');
    end
    
    datetime_str = f_nm(k+3:end-3);
    UTC = datevec([num2str(datetime_str) '/000000'],'yyyymmdd/HHMMSS');% this version of GEOS-Chem data do not have "START_TIME" in the attributes, thus we only can use the file name for time stamp!
    UTC(4) = i_hour-1;

    i = user_lat_index;% the index of Lat will be used
    j = user_lon_index;% the index of Lon will be used
    % get the profile over the site
    O3_a = O3(:,i,j,i_hour);% O3 VMR in [ppbv] ---> the dim of each trace gas is [lev,lat,lon,hour]
    NO2_a = NO2(:,i,j,i_hour);% NO2 VMR in [ppbv]
    NO_a = NO(:,i,j,i_hour);% NO VMR in [ppbv]
    SO2_a = SO2(:,i,j,i_hour);% SO2 VMR in [ppbv]
    HCHO_a = HCHO(:,i,j,i_hour);% HCHO VMR in [ppbv]
    P_a = P(:,i,j,i_hour);% surface pressure
    AIRDEN_a = AIRDEN(:,i,j,i_hour);% air density 
    BXHEIGHT_a = BXHEIGHT(:,i,j,i_hour);% box height m
    PBL_M_a = PBL_M(i,j,i_hour); % PBL height [m]
    PBL_L_a = PBL_L(i,j,i_hour); % PBL height level


    O3_a = reshape(O3_a,[],1);% reshape to simple matrix
    NO2_a = reshape(NO2_a,[],1);
    NO_a = reshape(NO_a,[],1);
    SO2_a = reshape(SO2_a,[],1);
    HCHO_a = reshape(HCHO_a,[],1);
    AIRDEN_a = reshape(AIRDEN_a,[],1);
    BXHEIGHT_a = reshape(BXHEIGHT_a,[],1);
    P_a = reshape(P_a,[],1);


    % prepare output profile table
    profile = table;
    profile.P =  P_a;
    profile.boxheight =  BXHEIGHT_a;
    profile.airdensity =  AIRDEN_a;
    profile.O3_vmr = O3_a;
    profile.NO2_vmr = NO2_a;
    profile.NO_vmr = NO_a;
    profile.SO2_vmr = SO2_a;
    profile.HCHO_vmr = HCHO_a;
    profile.UTC = repmat(datetime(UTC),height(profile),1);

    % calculate VCD
    DU = 2.69e20;% molec/m^3 note here we use al SI units in the following equs
    o3_vcd = sum(profile.airdensity.*1e6.*profile.O3_vmr.*1e-9.*profile.boxheight)/DU;
    no2_vcd = sum(profile.airdensity.*1e6.*profile.NO2_vmr.*1e-9.*profile.boxheight)/DU;
    no_vcd = sum(profile.airdensity.*1e6.*profile.NO_vmr.*1e-9.*profile.boxheight)/DU;
    so2_vcd = sum(profile.airdensity.*1e6.*profile.SO2_vmr.*1e-9.*profile.boxheight)/DU;
    hcho_vcd = sum(profile.airdensity.*1e6.*profile.HCHO_vmr.*1e-9.*profile.boxheight)/DU;


    % calculate VCD_PBL 
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


    % prepare output table
    VCD = table;
    VCD.UTC = datetime(UTC);
    VCD.o3 = o3_vcd;% all VCDs are in DU
    VCD.no2 = no2_vcd;
    VCD.no = no_vcd;
    VCD.so2 = so2_vcd;
    VCD.hcho = hcho_vcd;
    VCD.PBLH = PBL_M_a;% boundary layer height [m]


    VCD.o3_pbl = o3_vcd_pbl;% all tropospheric columns are in DU
    VCD.no2_pbl = no2_vcd_pbl;
    VCD.no_pbl = no_vcd_pbl;
    VCD.so2_pbl = so2_vcd_pbl;
    VCD.hcho_pbl = hcho_vcd_pbl;

    VCD.o3_surf = profile.O3_vmr(1);% all surface mixing ratios are in ppbv
    VCD.no2_surf = profile.NO2_vmr(1);
    VCD.no_surf = profile.NO_vmr(1);
    VCD.so2_surf = profile.SO2_vmr(1);
    VCD.hcho_surf = profile.HCHO_vmr(1);

    VCD.o3_CV = VCD.o3_surf./VCD.o3_pbl;% all convertion ratios are in ppbv/DU
    VCD.no2_CV = VCD.no2_surf./VCD.no2_pbl;
    VCD.no_CV = VCD.no_surf./VCD.no_pbl;
    VCD.so2_CV = VCD.so2_surf./VCD.so2_pbl;
    VCD.hcho_CV = VCD.hcho_surf./VCD.hcho_pbl;
    

    if i_hour == 1
        VCDs_1day = VCD;
        profiles_1day = profile;
    else
        VCDs_1day = [VCDs_1day;VCD];
        profiles_1day = [profiles_1day;profile];
    end
    VCD = []; profile =[];
    
end

% print 1-day output files
writetable(VCDs_1day,[output_file_path 'VCD_' f_nm(k+3:end-3) '_' site '.csv']);
writetable(profiles_1day,[output_file_path 'profile_' f_nm(k+3:end-3) '_' site '.csv']);






