function [profiles_1day,VCDs_1day] = read_GEOS_CHEM_v5(f_nm,met_file_path,user_lat_index,user_lon_index,output_file_path,site)
% this is a sub-function for get_GEOS_CHEM_all_v3 to read GEOS CHEM O3/NO2/SO2/HCHO data
% v5 update: the new GEOS-Chem data has trace gas concenration and
% meterological data saved in two different files

% inputs:
    %f_nm --> file name (include path) that has trace gas vmr [ppbv] profile
    % Note: in v5, the f_nm is the file that has trace gas concentration, the
    % corresponding meterological file will be refered based on name of of
    % trace gas file
    
    %user_lat_index,user_lon_index --> desired location index to readout from GEOS-chem output
% outputs:
    % profile --> single day (1hr) trace gas profiles at the given location
    % VCD --> single day (1hr) integareted VCDs (total columns [DU], trop/start columns [DU], trop height [m])

% by Xiaoyi 2018/11/05

%% check if the meterological file exist
metfile = [met_file_path 'GEOSChem.StateMet.' f_nm(strfind(f_nm,'SpeciesConc.')+length('SpeciesConc.'):end)];
if isfile(metfile)
    disp(['Find corresponding met file: ' metfile]);
else
    disp(['Cannot find the corresponding met file: ' metfile ' !']);
end

%%
% read in file info and variables

AIRDEN = ncread(metfile,'Met_AIRDEN');%'Dry air density' [kg m-3]
Av = 6.022140857e26;%[kmol-1]
M_Air = 28.97; % [kg/kmol]
AIRDEN = AIRDEN./M_Air.*Av; %convert air density to [molec/m3] 
BXHEIGHT = ncread(metfile,'Met_BXHEIGHT');% box height m
PBL_M = ncread(metfile,'Met_PBLH');% PBL m
PBL_L = ncread(metfile,'Met_PBLTOPL');%PBL level

% read in trace gas 
NO2 = ncread(f_nm,'SpeciesConc_NO2');% vmr original ratio mol/mol, not ppb!
O3 = ncread(f_nm,'SpeciesConc_O3');
HCHO = ncread(f_nm,'SpeciesConc_CH2O');
SO2 = ncread(f_nm,'SpeciesConc_SO2');
NO = ncread(f_nm,'SpeciesConc_NO');

% each file has 24 hours of data, so need loop for each hour
for i_hour = 1:24
    if ispc
        %k = strfind(f_nm,'\ts');
        k = strfind(f_nm,'\GEOSChem.SpeciesConc.');
    else
        %k = strfind(f_nm,'/ts');
        k = strfind(f_nm,'/GEOSChem.SpeciesConc.');
    end
    
    datetime_str = f_nm(k+length('\GEOSChem.SpeciesConc.'):end-4);
    %UTC = datevec([num2str(datetime_str) '/000000'],'yyyymmdd/HHMMSS');% this version of GEOS-Chem data do not have "START_TIME" in the attributes, thus we only can use the file name for time stamp!
    UTC = datevec([num2str(datetime_str) ],'yyyymmdd_HHMMz');% this version of GEOS-Chem data do not have "START_TIME" in the attributes, thus we only can use the file name for time stamp!
    UTC(4) = i_hour-1;

    i = user_lat_index;% the index of Lat will be used
    j = user_lon_index;% the index of Lon will be used
    % get the profile over the site
    O3_a = O3(j,i,:,i_hour);% O3 VMR in [mol/mol] ---> the dim of each trace gas is [lon,lat,lev,hour]
    NO2_a = NO2(j,i,:,i_hour);% NO2 VMR in [mol/mol]
    NO_a = NO(j,i,:,i_hour);% NO VMR in [mol/mol]
    SO2_a = SO2(j,i,:,i_hour);% SO2 VMR in [mol/mol]
    HCHO_a = HCHO(j,i,:,i_hour);% HCHO VMR in [mol/mol]

    AIRDEN_a = AIRDEN(j,i,:,i_hour);% air density  [molec/m3] 
    BXHEIGHT_a = BXHEIGHT(j,i,:,i_hour);% box height m
    PBL_M_a = PBL_M(j,i,i_hour); % PBL height [m] --> data size  [144,94,24]
    PBL_M_a = squeeze(PBL_M_a);
    PBL_L_a = PBL_L(j,i,i_hour); % PBL height level
    PBL_L_a = squeeze(PBL_L_a);

    O3_a = reshape(O3_a,[],1);% reshape to simple matrix
    NO2_a = reshape(NO2_a,[],1);
    NO_a = reshape(NO_a,[],1);
    SO2_a = reshape(SO2_a,[],1);
    HCHO_a = reshape(HCHO_a,[],1);
    AIRDEN_a = reshape(AIRDEN_a,[],1);
    BXHEIGHT_a = reshape(BXHEIGHT_a,[],1);

    % prepare output profile table
    profile = table;
    profile.boxheight =  BXHEIGHT_a;
    profile.airdensity =  AIRDEN_a;
    profile.O3_vmr = O3_a;
    profile.NO2_vmr = NO2_a;
    profile.NO_vmr = NO_a;
    profile.SO2_vmr = SO2_a;
    profile.HCHO_vmr = HCHO_a;
    profile.UTC = repmat(datetime(UTC),height(profile),1);

    % calculate VCD
    DU = 2.69e20;% molec/m^2 note here we use al SI units in the following equs
    o3_vcd = sum(profile.airdensity.*profile.O3_vmr.*profile.boxheight)/DU; %[DU]
    no2_vcd = sum(profile.airdensity.*profile.NO2_vmr.*profile.boxheight)/DU;
    no_vcd = sum(profile.airdensity.*profile.NO_vmr.*profile.boxheight)/DU;
    so2_vcd = sum(profile.airdensity.*profile.SO2_vmr.*profile.boxheight)/DU;
    hcho_vcd = sum(profile.airdensity.*profile.HCHO_vmr.*profile.boxheight)/DU;


    % calculate VCD_PBL 
    PBL_level = fix(PBL_L_a);% this is the integer part of PBL level
    PBL_level_decimal = PBL_L_a - fix(PBL_L_a);% this is the decimal part of PBL level
    o3_vcd_pbl = sum(profile.airdensity(1:PBL_level).*profile.O3_vmr(1:PBL_level).*profile.boxheight(1:PBL_level))/DU ...
        + sum(PBL_level_decimal.*profile.airdensity(PBL_level+1).*profile.O3_vmr(PBL_level+1).*profile.boxheight(PBL_level+1))/DU;
    no2_vcd_pbl = sum(profile.airdensity(1:PBL_level).*profile.NO2_vmr(1:PBL_level).*profile.boxheight(1:PBL_level))/DU ...
        + PBL_level_decimal.*profile.airdensity(PBL_level+1).*profile.NO2_vmr(PBL_level+1).*profile.boxheight(PBL_level+1)/DU;
    no_vcd_pbl = sum(profile.airdensity(1:PBL_level).*profile.NO_vmr(1:PBL_level).*profile.boxheight(1:PBL_level))/DU ...
        + PBL_level_decimal.*profile.airdensity(PBL_level+1).*profile.NO_vmr(PBL_level+1).*profile.boxheight(PBL_level+1)/DU;
    so2_vcd_pbl = sum(profile.airdensity(1:PBL_level).*profile.SO2_vmr(1:PBL_level).*profile.boxheight(1:PBL_level))/DU ...
        + PBL_level_decimal.*profile.airdensity(PBL_level+1).*profile.SO2_vmr(PBL_level+1).*profile.boxheight(PBL_level+1)/DU;
    hcho_vcd_pbl = sum(profile.airdensity(1:PBL_level).*profile.HCHO_vmr(1:PBL_level).*profile.boxheight(1:PBL_level))/DU ...
        + PBL_level_decimal.*profile.airdensity(PBL_level+1).*profile.HCHO_vmr(PBL_level+1).*profile.boxheight(PBL_level+1)/DU;


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

    VCD.o3_surf = profile.O3_vmr(1).*1e9;% all surface mixing ratios are in ppbv
    VCD.no2_surf = profile.NO2_vmr(1).*1e9;
    VCD.no_surf = profile.NO_vmr(1).*1e9;
    VCD.so2_surf = profile.SO2_vmr(1).*1e9;
    VCD.hcho_surf = profile.HCHO_vmr(1).*1e9;

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
%writetable(VCDs_1day,[output_file_path 'VCD_' f_nm(k+3:end-3) '_' site '.csv']);
writetable(VCDs_1day,[output_file_path 'VCD_' datestr(UTC,'yyyymmdd') '_' site '.csv']);
%writetable(profiles_1day,[output_file_path 'profile_' f_nm(k+3:end-3) '_' site '.csv']);
writetable(profiles_1day,[output_file_path 'profile_' datestr(UTC,'yyyymmdd') '_' site '.csv']);






