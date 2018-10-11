function plot_all_profiles(profiles)

%data = profilesDownsview;
data = profiles;
timestamps = unique(data.UTC);
figure;hold all;
for i = 1:numel(timestamps)
    TF = data.UTC == timestamps(i,:);
    plot(data.O3_vmr(TF,:),data.Z(TF,:));
    legends{i,:} = datestr(timestamps(i,:));
end
xlabel('O_3 VMR [ppbv]');
ylabel('Height [m]');
legend(legends);

figure;hold all;
for i = 1:numel(timestamps)
    TF = data.UTC == timestamps(i,:);
    plot(data.O3_vmr(TF,:),data.pressure(TF,:));
end
xlabel('O_3 VMR [ppbv]');
ylabel('Pressure [hPa]');
legend(legends);
