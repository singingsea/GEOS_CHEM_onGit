function monthly_table = get_GEOS_CHEM_LUT()
DU = 2.6870e+16;
use_total_column = false; % if true, then use total column to calculate conversion ratio
% calculate LUT from GEOS-CHEM summary file
use_time_window = true;
start_time = '2015-01-01';
end_time = '2018-01-01';

site = 'Downsview';timezone_offset = -5;
%site = 'Egbert';timezone_offset = -5;
%site = 'FortMcKay';timezone_offset = -7;
%site = 'Beijing';timezone_offset = 8;
%site = 'LosAngeles';timezone_offset = -7;

save_fig = 1;

if ispc
    addpath('C:\Users\ZhaoX\Documents\MATLAB\matlab\');
    %summary_file_output_file_path = 'C:\Projects\GEOS_CHEM\output\summary_files\';
    summary_file_output_file_path = 'C:\Projects\GEOS_CHEM\output\summary_files\MERRA2_05x0625_47L_NA_V12\';
    if use_total_column
        plot_path = [summary_file_output_file_path site '_LUT_totalcolumn\MERRA2_05x0625_47L_NA_V12\'];  
    else
        plot_path = [summary_file_output_file_path site '_LUT\MERRA2_05x0625_47L_NA_V12\'];
    end
else
    addpath('/export/data/home/xizhao/matlab/');
    summary_file_output_file_path = '/export/data/home/xizhao/GEOS_CHEM/output/summary_files/';
    if use_total_column
        plot_path = [summary_file_output_file_path site '_LUT_totalcolumn/'];  
    else
        plot_path = [summary_file_output_file_path site '_LUT/'];
    end
end
mkdir(plot_path);
matfile_nm = [summary_file_output_file_path 'VCDs_Profile_' site '.mat'];
load(matfile_nm);


%% filter by time window
if use_time_window == true
    TF_time = (VCDs.UTC >= datetime(start_time)) & (VCDs.UTC <= datetime(end_time));
    VCDs(~TF_time,:) = [];
end

%%

VCDs.LST = VCDs.UTC + hours(timezone_offset);% get LST
monthly_table = table;
j=1;
if use_total_column
    VCDs.ratio = VCDs.no2_surf./VCDs.no2;
else
    VCDs.ratio = VCDs.no2_surf./VCDs.no2_pbl;
end
VCDs.doy = day(VCDs.UTC,'dayofyear');

for month = 1:12
    for hour = 0:23
        data = VCDs;
        TF_month = data.UTC.Month == month;
        TF_hour = data.UTC.Hour == hour;
        TF = TF_month & TF_hour;
        data(~TF,:) = [];
        data(isnan(data.ratio),:)=[];
        data(isinf(data.ratio),:)=[];
        
        monthly_table.month(j,:) = month;
        monthly_table.hour(j,:) = hour;
        monthly_table.hour_LST(j,:) = unique(data.LST.Hour);
        monthly_table.ratio(j,:) = mean(data.ratio./DU);
        monthly_table.ratio_std(j,:) = std(data.ratio./DU);
        N = height(data);
        monthly_table.ratio_err(j,:) = std(data.ratio./DU)/(N)^0.5;
        j=j+1;
    end
end

% UTC plot
fig_name = 'Converstion_ratio_UTC';
figure;hold all;
C= jet(12);C(13,:)=C(3,:);
C(14,:)=C(2,:);
C(15,:)=C(1,:);
C(1:3,:)=[];
for month = 1:12
    TF = monthly_table.month == month;
    x = monthly_table.hour(TF,:);
    y = monthly_table.ratio(TF,:);
    b = monthly_table.ratio_std(TF,:);
    plot(x,y.*DU,'.-','Color',C(month,:));
    %[hl, hp] = boundedline(x, y, b,  'transparency', 0.5);
end
xlim([0 23]);
ylabel('Convertion ratio [ppbv/DU]');
xlabel('UTC [Hour]');
month = 1:12;month = month';legend(num2str(month));
print_setting(1/4,save_fig,[plot_path fig_name]);


fig_name = 'Converstion_ratio_LTC';
figure;hold all;
C= jet(12);C(13,:)=C(3,:);
C(14,:)=C(2,:);
C(15,:)=C(1,:);
C(1:3,:)=[];
for month = 1:12
    TF = monthly_table.month == month;
    x = monthly_table.hour(TF,:);
    y = monthly_table.ratio(TF,:);
    start_hr = 24-abs(timezone_offset)+1;
    for i =1:start_hr-1
        y_LST(i,:) = y(i + abs(timezone_offset));% this works for any input timezone
    end
    
    y_LST(start_hr:24,:) = y(1:abs(timezone_offset));
 
    b = monthly_table.ratio_err(TF,:);
    %plot(x,y_LST,'.-','Color',C(month,:));
    [hl, hp] = boundedline(x, y_LST.*DU, b.*DU, 'cmap', C(month,:), 'alpha');
    %[hl, hp] = boundedline(monthly_table.hour_LST(TF,:), y.*DU, b.*DU,'cmap', C(month,:), 'alpha');% this is to test if the LST convertion is correct
end
xlim([0 23]);
ylabel('Convertion ratio [ppbv/DU]');
xlabel('LST [Hour]');
print_setting(1/4,save_fig,[plot_path fig_name]);

if use_total_column
    save([plot_path 'LUT_totalcolumn.mat'],'monthly_table');
else
    save([plot_path 'LUT.mat'],'monthly_table');
end
