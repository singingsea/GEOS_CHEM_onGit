function GEM_MACH = make_dummy_GEM_file()
% this function read in the output of "get_GEOS_CHEM_summary_file.m" and
% make GEM-MACH type dummy file to perform Pandora NO2 surface conversion
load('C:\Projects\GEOS_CHEM\output\summary_files\VCDs_Profile_Downsview.mat');
DU = 2.6870e+16;
GEM_MACH = table;
GEM_MACH.t = VCDs.UTC;
GEM_MACH.gem_vcd = VCDs.no2_pbl.*DU;% this is already pbl VCD!
GEM_MACH.gem_surface = VCDs.no2_surf; % this is surface NO2
N = height(GEM_MACH);
GEM_MACH.gem_freetrop = zeros(N,1); % this is free trop VCD, since the gem_vsd is already pbl VCD, this value is set to 0
GEM_MACH.gem_pbl = VCDs.PBLH; % 
