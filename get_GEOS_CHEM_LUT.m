function monthly_table = get_GEOS_CHEM_LUT()
% calculate LUT from GEOS-CHEM summary file

site = 'Downsview';timezone_offset = -5;
%site = 'Egbert';timezone_offset = -5;
%site = 'FortMcKay';timezone_offset = -7;
%site = 'Beijing';timezone_offset = 8;
%site = 'LosAngeles';timezone_offset = -7;

save_fig = 1;

if ispc
    addpath('C:\Users\ZhaoX\Documents\MATLAB\matlab\');
    summary_file_output_file_path = 'C:\Projects\GEOS_CHEM\output\summary_files\';
    plot_path = [summary_file_output_file_path site '_LUT\'];
else
    addpath('/export/data/home/xizhao/matlab/');
    summary_file_output_file_path = '/export/data/home/xizhao/GEOS_CHEM/output/summary_files/';
    plot_path = [summary_file_output_file_path site '_LUT/'];
end
mkdir(plot_path);
matfile_nm = [summary_file_output_file_path 'VCDs_Profile_' site '.mat'];
load(matfile_nm);


%%


monthly_table = table;
j=1;
VCDs.ratio = VCDs.no2_surf./VCDs.no2_pbl;
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
        monthly_table.ratio(j,:) = mean(data.ratio);
        monthly_table.ratio_std(j,:) = std(data.ratio);
        N = height(data);
        monthly_table.ratio_err(j,:) = std(data.ratio)/(N)^0.5;
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
    plot(x,y,'.-','Color',C(month,:));
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
    [hl, hp] = boundedline(x, y_LST, b, 'cmap', C(month,:), 'alpha');
end
xlim([0 23]);
ylabel('Convertion ratio [ppbv/DU]');
xlabel('LST [Hour]');
print_setting(1/4,save_fig,[plot_path fig_name]);

save([plot_path 'LUT.mat'],'monthly_table');
