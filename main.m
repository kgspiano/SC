% =========================================================================================================================================================
% Codes for
% "Periodic clustering of simple and complex cells in visual cortex"
% Gwangsu Kim, Jaeson Jang, and Se-Bum Paik*
%
% *Contact: sbpaik@kaist.ac.kr
%
% =========================================================================================================================================================

%% ----------- Initializing
clear;
clc;

pathtmp = pwd; % set current dir 
addpath(genpath(pathtmp));

%% ----------- Setting parameters

% 1) moire interference parameters
d_off_off = 142; % lattice constant of OFF mosaic
d_on_on = d_off_off*8/7; % lattice constant of ON mosaic
RF_sig_off = d_off_off/2; % Size of OFF RF
RF_sig_on = d_on_on/2; % Size of ON RF

dtheta = 0; % for theta = 0, S = alpha+1/alpha --> for S =8, alpha =1/7
s_factor = (1+1/7)/sqrt((1/7)^2);  % Scaling factor of the interference pattern
% Note : s~8 for monkey. orientation map period : 0.750 mmc ~ 1mm in retina (Paik & Ringach 2011)

RGC_noise_factor = 0.12; %% Noise in units of d_off_off
vec1_off = 1; vec2_off = cos(60*pi/180) + 1i*sin(60*pi/180); % Basis for generating lattices
vec1_on = cos(dtheta) + 1i*sin(dtheta); vec2_on = cos(60*pi/180+dtheta) + 1i*sin(60*pi/180+dtheta);
cell_x_num = 40; % number of off cells in x axis
cell_cut = cell_x_num/2+0.000001;

% 2) Stimulus space parameters (Receptive field)
size_x = d_off_off*6; % size of stimulus space canvus, [um, RS]   ~ 2.78 deg
size_y = d_off_off*6;
x_axis = -size_x:20:size_x; % coordinate axis
y_axis = -size_y:20:size_y;
[xx, yy] = meshgrid(x_axis, y_axis);
rr = xx+1i*yy;

% 3) Retina-cortical mapping parameters
% Dacey, & Petersen,1992, A = 0.1+4.21E+0.038E^2 (from mmr to deg)
% Gauthier et al., 2009 ~9mm eccentricity
% 9mmr eccentric location (~41 degree) --> dA/dE|E=9mmr 4.9 deg/mmr
% Vanessen et al., 1984, cortical magnification factor : 0.146 mmc/deg,
% 0.715 mmc/mmr
RS2VS = 0.715;   % 0.715 [mm, cortex]/[mm, retina]
umr2deg = 0.0049;     % 0.0049 deg/[um, retina]
sig_con = 0.325*d_off_off;  % [um, retina] Std of Gaussian connection strength, ~ 6 ON/OFF RGCs in 2 sigma range, check SC_revision_sigconparametersearch.m


% 4) V1 neurons
num_per_1d = 6; % V1 cell density (number/doff-off); Resolution of the map. In our simulation, num_per_1d = 12 was used

xaxis = -cell_cut*d_off_off+size_x:d_off_off/num_per_1d:cell_cut*d_off_off-size_x; 
yaxis = -cell_cut*d_off_off+size_y:d_off_off/num_per_1d:cell_cut*d_off_off-size_y;
[xtemp, ytemp] = meshgrid(xaxis, yaxis); 
V1_cell_locations = xtemp+1i*ytemp;  % locations of V1 cells, [um, RS]
V1_cell_locations = V1_cell_locations(:);
V1_x = real(V1_cell_locations);
V1_y = imag(V1_cell_locations);
PX2VS = mean(diff(xaxis)*RS2VS);
sig_filter = round(num_per_1d*1.25); % map filtering range

% 5) Colormaps
cmap = flip([[linspace(1, 1, 32) linspace(1,0, 32)]' [linspace(0,1,32) linspace(1,0,32)]' [linspace(0,1,32) linspace(1, 1, 32)]'], 1);
cmap_off = flip(cmap(1:32, :), 1);
cmap_on = cmap(33:64, :);

% 6) Params from Kemkow et al., data (penetration length, interval, etc)
load('Parameters_data') % to compare with data

itermap = 2; % Iteration (~ 90s for 1 iteration)


SI_filtered_crops = cell(1, itermap);
dCOM_filtered_crops = cell(1, itermap);
Ori_filtered_crops = cell(1, itermap);

%% ----------- Main calculation
disp('Simulating cortical map of Orientation, dON-OFF, Simpleness index...')
for itermapind = 1:itermap

    disp(['iteration [' num2str(itermapind) '/' num2str(itermap) ']'])
    %% Generating moire interference between ON and OFF RGC mosaics
    off_cell_locations_ideal = [];
    on_cell_locations_ideal = [];
    for ii = -cell_x_num:cell_x_num
        for jj = -cell_x_num:cell_x_num
            loc_temp_off = d_off_off*(ii*vec1_off + jj*vec2_off);
            off_cell_locations_ideal = [off_cell_locations_ideal, loc_temp_off];
            loc_temp_on = d_on_on*(ii*vec1_on+jj*vec2_on);
            on_cell_locations_ideal = [on_cell_locations_ideal, loc_temp_on];
        end
    end
    
    % Adding noise to the mosaic
    
    sig_noise = d_off_off*RGC_noise_factor;
    noisetmp = sig_noise*randn(size(off_cell_locations_ideal)); % noise applied to on/off mosaics
    noisetmp2 = sig_noise*randn(size(off_cell_locations_ideal));
    noisetmp3 = sig_noise*randn(size(off_cell_locations_ideal));
    noisetmp4 = sig_noise*randn(size(off_cell_locations_ideal));
    
    off_cell_locations = off_cell_locations_ideal+noisetmp+1i*noisetmp2;
    on_cell_locations = on_cell_locations_ideal+noisetmp3+1i*noisetmp4;
    
    % cutting edges

    bool_off = logical((real(off_cell_locations)>cell_cut*d_off_off)+(real(off_cell_locations)<-cell_cut*d_off_off)+(imag(off_cell_locations)>cell_cut*d_off_off)+(imag(off_cell_locations)<-cell_cut*d_off_off));
    bool_on = logical((real(on_cell_locations)>cell_cut*d_off_off)+(real(on_cell_locations)<-cell_cut*d_off_off)+(imag(on_cell_locations)>cell_cut*d_off_off)+(imag(on_cell_locations)<-cell_cut*d_off_off));
    
    off_cell_locations(bool_off) = [];   % [umr]
    on_cell_locations(bool_on) = [];     % [umr]

    
    %% Generating connection to V1 from feedforward RGC inputs
    
    V1_N = length(V1_cell_locations);
    V1_ONsub_locs = []; V1_OFFsub_locs = [];
    V1_ONsub_weights = []; V1_OFFsub_weights = [];
    V1_ind_ON = []; V1_ind_OFF = [];
    
    V1_weights_ON = cell(V1_N,1); % connection weights of ON and OFF RGC
    V1_weights_OFF = cell(V1_N,1);

    
    for ii = 1:V1_N
        V1_loc_temp = V1_cell_locations(ii);
        off_temp = off_cell_locations-V1_loc_temp;
        on_temp = on_cell_locations-V1_loc_temp;
        significant_off_ind = abs(off_temp)<3*sig_con;
        significant_on_ind = abs(on_temp)<3*sig_con;
        significant_cells_off = off_temp(significant_off_ind);
        significant_cells_on = on_temp(significant_on_ind);
        V1_OFFsub_locs = [V1_OFFsub_locs, significant_cells_off];
        V1_ONsub_locs = [V1_ONsub_locs, significant_cells_on];
        
        for jj = 1:length(significant_cells_off)
            r0 = significant_cells_off(jj);
            weight_temp = exp(-(abs(r0))^2/2/sig_con^2);
            V1_OFFsub_weights = [V1_OFFsub_weights, weight_temp];
            V1_ind_OFF = [V1_ind_OFF, ii];
            V1_weights_OFF{ii} = [V1_weights_OFF{ii}, weight_temp];
        end
        for jj = 1:length(significant_cells_on)
            r0 = significant_cells_on(jj);
            weight_temp = exp(-(abs(r0))^2/2/sig_con^2);
            V1_ONsub_weights = [V1_ONsub_weights, weight_temp];
            V1_ind_ON = [V1_ind_ON, ii];
            V1_weights_ON{ii} = [V1_weights_ON{ii},weight_temp];
        end
    end
    
    % find monocontrast cell
    ON_domi_ind = zeros(1,V1_N);
    OFF_domi_ind = zeros(1,V1_N);
    
    thre_mono = 2;  % threshold for ON/OFF dominance
    for ii = 1:V1_N
        ON_weights_sum_temp = sum(V1_weights_ON{ii});
        OFF_weights_sum_temp = sum(V1_weights_OFF{ii});
        if ON_weights_sum_temp>thre_mono*OFF_weights_sum_temp
            ON_domi_ind(ii) = 1;
        elseif OFF_weights_sum_temp>thre_mono*ON_weights_sum_temp
            OFF_domi_ind(ii) = 1;
        end
    end

    
    %% Calculating Ori, dCOM, seg
    SI_V1 = zeros(1, V1_N);
    dCOM_V1 = zeros(1,V1_N);
    Ori_V1 = zeros(1,V1_N);
    dCOM_angle_V1 = zeros(1,V1_N);
   
    for ii = 1:V1_N
        RF_temp_ON = 0;
        RF_temp_OFF = 0;
        ind = ii;
        OFF_temp = V1_OFFsub_locs(V1_ind_OFF == ind);
        ON_temp = V1_ONsub_locs(V1_ind_ON == ind);
        
        for jj = 1:length(OFF_temp)
            r0 = OFF_temp(jj);
            OFF_RF_cen = exp(-0.5*(abs(rr-r0)/RF_sig_off).^2); % center RF
            OFF_RF_cen = OFF_RF_cen/sum(OFF_RF_cen(:)); % normalization
            OFF_RF_sur = exp(-0.5*(abs(rr-r0)/(RF_sig_off*3)).^2); % surround has 3 times larger sigma
            OFF_RF_sur = -OFF_RF_sur/sum(OFF_RF_sur(:)); % normalization
            OFF_RF = -(OFF_RF_cen + OFF_RF_sur); % center + surround ON RF
            weight_temp = exp(-(abs(r0))^2/2/sig_con^2);
            RF_temp_OFF = RF_temp_OFF + weight_temp*OFF_RF;
        end
        for jj = 1:length(ON_temp)
            r0 = ON_temp(jj);
            ON_RF_cen = exp(-0.5*(abs(rr-r0)/RF_sig_on).^2); % center RF
            ON_RF_cen = ON_RF_cen/sum(ON_RF_cen(:)); % normalization
            ON_RF_sur = exp(-0.5*(abs(rr-r0)/(RF_sig_on*3)).^2); % surround has 3 times larger sigma
            ON_RF_sur = -ON_RF_sur/sum(ON_RF_sur(:)); % normalization
            ON_RF = ON_RF_cen + ON_RF_sur; % center + surround ON RF
            
            weight_temp = exp(-(abs(r0))^2/2/sig_con^2);
            RF_temp_ON = RF_temp_ON + weight_temp*ON_RF;
        end
        
        SI_V1(ii) = sum(sum(abs((abs(RF_temp_OFF)-abs(RF_temp_ON)))))/sum(sum((abs(RF_temp_OFF)+abs(RF_temp_ON))));
        
        [xxtmp,yytmp] = meshgrid(1:size(RF_temp_ON,1), 1:size(RF_temp_ON,2));
        RF_temp_ON = RF_temp_ON.*(RF_temp_ON>0);
        RF_temp_OFF = abs(RF_temp_OFF.*(RF_temp_OFF<0));
        xCOMtmp = sum(xxtmp(:).*RF_temp_ON(:))/sum(RF_temp_ON(:));
        yCOMtmp = sum(yytmp(:).*RF_temp_ON(:))/sum(RF_temp_ON(:));
        xCOMtmp2 = sum(xxtmp(:).*RF_temp_OFF(:))/sum(RF_temp_OFF(:));
        yCOMtmp2 = sum(yytmp(:).*RF_temp_OFF(:))/sum(RF_temp_OFF(:));
        dCOM_V1(ii) = sqrt((xCOMtmp-xCOMtmp2)^2+(yCOMtmp-yCOMtmp2)^2)*mean(diff(x_axis))/RF_sig_off; % [sigma_off]
        
        dCOMtmp = (xCOMtmp-xCOMtmp2)+1i*(yCOMtmp-yCOMtmp2);
        dCOM_angle_V1(ii) = angle(dCOMtmp);
        
    end
    
    % N = 100,000 ~ 12 min
    % Reshaping to map
    Ori_V1_save = Ori_V1;
    seg_V1_save = SI_V1;
    dCOM_V1_save = dCOM_V1;
    Ori_V1 = dCOM_angle_V1-pi/2;
    
    SI_V1((ON_domi_ind|OFF_domi_ind)) = NaN;
    dCOM_V1((ON_domi_ind|OFF_domi_ind)) = NaN;
    Ori_V1((ON_domi_ind|OFF_domi_ind)) = NaN;
    SI_map = reshape(SI_V1, [length(xaxis),length(yaxis)]);
    dCOM_map = reshape(dCOM_V1, [length(xaxis),length(yaxis)]);
    Ori_map = reshape(Ori_V1, [length(xaxis),length(yaxis)]);
    
    % Applying 2D spatial filters
    
    [xxx, yyy] = meshgrid(-3*sig_filter:3*sig_filter, -3*sig_filter:3*sig_filter);
    rrr = xxx+1i*yyy;
    Gauss_filter = exp(-abs(rrr).^2/2/sig_filter^2);
    Gauss_filter = Gauss_filter/sum(Gauss_filter(:));
    % seg. map filtering
    SI_map_nan = SI_map;
    SI_pad1 = [zeros(sig_filter*3, length(xaxis))/0;SI_map_nan;zeros(sig_filter*3, length(xaxis))/0];
    SI_pad2 = [zeros((length(yaxis)+sig_filter*6), sig_filter*3)/0, SI_pad1, zeros((length(yaxis)+sig_filter*6), sig_filter*3)/0];
    SI_filtered = zeros(length(xaxis), length(yaxis));
    for ii = 1:length(xaxis)
        for jj = 1:length(yaxis)
            
            map_temp = SI_pad2(jj:jj+sig_filter*6, ii:ii+sig_filter*6);
            linear_temp = Gauss_filter.*map_temp;
            temp = nansum(linear_temp(:))/sum(sum(Gauss_filter(find(~isnan(linear_temp(:))))));
            SI_filtered(jj, ii) = temp;
        end
    end
    SI_filtered = rescale(SI_filtered, min(SI_map(:)), max(SI_map(:)));
    % dCOM map filtering
    dCOM_map_nan = dCOM_map;
    dCOM_pad1 = [zeros(sig_filter*3, length(xaxis))/0;dCOM_map_nan;zeros(sig_filter*3, length(xaxis))/0];
    dCOM_pad2 = [zeros((length(yaxis)+sig_filter*6), sig_filter*3)/0, dCOM_pad1, zeros((length(yaxis)+sig_filter*6), sig_filter*3)/0];
    dCOM_filtered = zeros(length(xaxis), length(yaxis));
    for ii = 1:length(xaxis)
        for jj = 1:length(yaxis)
            
            map_temp = dCOM_pad2(jj:jj+sig_filter*6, ii:ii+sig_filter*6);
            linear_temp = Gauss_filter.*map_temp;
            temp = nansum(linear_temp(:))/sum(sum(Gauss_filter(find(~isnan(linear_temp(:))))));
            dCOM_filtered(jj, ii) = temp;
        end
    end
    dCOM_filtered = rescale(dCOM_filtered, min(dCOM_map(:)), max(dCOM_map(:)));
    
    Ori_pad1 = [zeros(sig_filter*3, length(xaxis))/0;exp(2i*Ori_map);zeros(sig_filter*3, length(xaxis))/0];
    Ori_pad2 = [zeros((length(yaxis)+sig_filter*6), sig_filter*3)/0, Ori_pad1, zeros((length(yaxis)+sig_filter*6), sig_filter*3)/0];
    Ori_filtered = zeros(length(xaxis), length(yaxis));
    
    for ii = 1:length(xaxis)
        for jj = 1:length(yaxis)
            
            map_temp = Ori_pad2(jj:jj+sig_filter*6, ii:ii+sig_filter*6);
            linear_temp = Gauss_filter.*map_temp;
            temp = nansum(linear_temp(:));
            Ori_filtered(jj, ii) = angle(temp)/2;
        end
    end
    
    
    %% Model cortical map of orientation, dON-OFF, and SI

    PX_crop = num_per_1d*10; % 10 d_off_off distance
    RS_crop = PX_crop*PX2VS/RS2VS;
    PX_center = round(size(Ori_filtered, 1)/2);
    PX_range = PX_center-PX_crop:PX_center+PX_crop;
    S_RS = s_factor*d_off_off;  % moire period in retinal space (um)
    S_VS = S_RS*RS2VS; % moire period in cortical space (um)
   
    img_range = (PX_range-PX_center)*PX2VS/S_VS; % image in the unit of moire interference period
    Ori_filtered_added = Ori_filtered+pi/4;
    Ori_filtered_added(Ori_filtered_added > 1/2*pi) = -pi+Ori_filtered_added(Ori_filtered_added > 1/2*pi);
    Ori_filtered_crop = Ori_filtered_added(PX_range, PX_range);

    
    dCOM_filtered_crop = dCOM_filtered(PX_range, PX_range);
    dCOM_filtered_crop = rescale(dCOM_filtered_crop, min(dCOM_map(:)), max(dCOM_map(:)));
    
    
    SI_filtered_crop = SI_filtered(PX_range, PX_range);
    SI_filtered_crop = rescale(SI_filtered_crop, min(SI_map(:)), max(SI_map(:)));

    SI_filtered_crops{itermapind} = SI_filtered_crop;
    dCOM_filtered_crops{itermapind} = dCOM_filtered_crop;
    Ori_filtered_crops{itermapind} = Ori_filtered_crop;
    
end

%% ---------------- Analysis of the simulated cortical maps
disp('Analyzing cortical maps...')

itermapind = 1; % example map
SI_filtered_crop = SI_filtered_crops{itermapind};
dCOM_filtered_crop = dCOM_filtered_crops{itermapind};
Ori_filtered_crop = Ori_filtered_crops{itermapind};
%% Model cortical maps 
%RGC mosaics
figure
subplot(3,4,1)
scatter(real(off_cell_locations)/S_RS, imag(off_cell_locations)/S_RS, 'b.')
axis image xy
xlim([-RS_crop/S_RS RS_crop/S_RS])
ylim([-RS_crop/S_RS RS_crop/S_RS])
hold on
scatter(real(on_cell_locations)/S_RS, imag(on_cell_locations)/S_RS, 'r.')
axis image xy
xlim([-RS_crop/S_RS RS_crop/S_RS])
ylim([-RS_crop/S_RS RS_crop/S_RS])
set(gcf, 'Position', [200, 200, 200, 800])
title('ON-OFF RGC mosaics')
axis off;%colorbar

subplot(3,4,10)
imagesc(SI_filtered_crop)
colormap(gca, cmap);axis image xy;hold on;cb=colorbar;
colormap(gca, 'hot');axis image xy;hold on;cb=colorbar;
caxis([0 0.5]); set(cb, 'YTick', [0 0.5]); set(gcf, 'Position', [200, 200, 300, 800])
title('Map of SI')
hold on;
axis off;
subplot(3,4,6)
imagesc(dCOM_filtered_crop)
colormap(gca, 'gray');axis image xy;hold on;cb=colorbar;
caxis([0 2]); set(cb, 'YTick', [0 2])
title('Map of d_{ON-OFF}')
hold on; axis off;
subplot(3,4,2)
imagesc(Ori_filtered_crop)
colormap(gca, 'hsv');axis image xy;hold on;cb=colorbar;
hold on; axis off;
title('Map of orientation')


%% Analysis of the map structure

S_PX = d_off_off*s_factor*RS2VS/PX2VS;
xlen = size(SI_filtered_crop, 2);
ylen = size(SI_filtered_crop, 1);
dist_thr = 2.6/1.1*S_PX; % 1.1mm = 96 px = 1 lambda thus, data : 2.6mm --> 226.9 px

% Loading maps
SI_filtered_set = zeros(ylen,xlen, itermap);
dCOM_filtered_set = zeros(ylen,xlen, itermap);
ori_filtered_set = zeros(ylen,xlen, itermap);
for itermapind = 1:itermap

    SI_filtered_set(:,:, itermapind) = SI_filtered_crops{itermapind};
    dCOM_filtered_set(:,:, itermapind) = dCOM_filtered_crops{itermapind};
    ori_filtered_set(:,:, itermapind) = Ori_filtered_crops{itermapind};
    
end

%% Calculating Ori-SI (dON-OFF) correlation for each penetration (phase, lag adjusted)
iter = 100;

% Getting cross-section of each map
SI_cs = zeros(iter, NN); % NN : # of data points in each penetration
dCOM_cs = zeros(iter, NN);
ori_cs = zeros(iter, NN);
for ii = 1:iter
    
    % getting random two dots
    dist_tmp = 0;
    while dist_tmp<dist_thr
        xtmp1 = randi(xlen);
        ytmp1 = randi(ylen);
        xtmp2 = randi(xlen);
        ytmp2 = randi(ylen);
        dist_tmp = sqrt((xtmp2-xtmp1)^2+(ytmp2-ytmp1)^2);
    end
    
    xtmp22 = xtmp1+(xtmp2-xtmp1)*dist_thr/dist_tmp;
    ytmp22 = ytmp1+(ytmp2-ytmp1)*dist_thr/dist_tmp;
    xi = [xtmp1, xtmp22];
    yi = [ytmp1, ytmp22];
    % getting profiles for SI, d, ori
    ntmp = randi(itermap);
    SI_filteredtmp = squeeze(SI_filtered_set(:,:,ntmp));
    dCOM_filteredtmp = squeeze(dCOM_filtered_set(:,:,ntmp));
    ori_filteredtmp = squeeze(ori_filtered_set(:,:,ntmp));
    
    SI_lin = improfile(SI_filteredtmp, xi, yi, NN);
    dCOM_lin = improfile(dCOM_filteredtmp, xi, yi, NN);
    ori_lin = improfile(ori_filteredtmp, xi, yi, NN);
    
    SI_cs(ii, :) = SI_lin;
    dCOM_cs(ii, :) = dCOM_lin;
    ori_cs(ii, :) = ori_lin;
    
end


phicand = 0:pi/18:2*pi;
Oris_md = ori_cs; % -pi/2 ~ pi deg
SI_md = SI_cs;
dCOM_md = dCOM_cs;
X = 0:0.1:2.6;
Oris_fitted_md = cell(1,size(ori_cs, 1));
seg_fitted_md = cell(1,size(ori_cs, 1));
dCOM_fitted_md = cell(1,size(ori_cs, 1));
%% calculate transformed signal and correlation
for dd = 1:size(ori_cs, 1)
    Ori_ltmp = Oris_md(dd, :);
    SI_ltmp = SI_md(dd, :);
    dCOM_ltmp = dCOM_md(dd,:);
    
    %% search lag
    lagcand = -5:1:5;
    r_tot = zeros(length(phicand), length(lagcand));
    p_tot = zeros(length(phicand), length(lagcand));
    for ii = 1:length(phicand)
        phitmp = phicand(ii);
        ori_tr = cos(Ori_ltmp*2+phitmp);
        for jj = 1:length(lagcand)
            lagtmp = lagcand(jj);
            if lagtmp<0
                tmp1 = ori_tr(1:end+lagtmp);
                tmp2 = SI_ltmp(1-lagtmp:end);
                [r,p] = corrcoef(tmp1, tmp2, 'Rows', 'complete');
                r_tot(ii,jj) = r(2);
                p_tot(ii,jj) = p(2);
            else
                tmp1 = ori_tr(1+lagtmp:end);
                tmp2 = SI_ltmp(1:end-lagtmp);
                [r,p] = corrcoef(tmp1, tmp2, 'Rows', 'complete');
                r_tot(ii,jj) = r(2);
                p_tot(ii,jj) = p(2);
            end
        end
    end
    
    [itmp, jtmp] = find(r_tot==max(r_tot(:)));
    phitmp = phicand(itmp(1));
    lagtmp = lagcand(jtmp(1));
    ori_tr = cos(Ori_ltmp*2+phitmp);
    if lagtmp<0
        tmp1 = ori_tr(1:end+lagtmp);
        tmp2 = SI_ltmp(1-lagtmp:end);
        tmp3 = dCOM_ltmp(1-lagtmp:end);
        [r,p] = corrcoef(tmp1, tmp2, 'Rows', 'complete');
        
    else
        tmp1 = ori_tr(1+lagtmp:end);
        tmp2 = SI_ltmp(1:end-lagtmp);
        tmp3 = dCOM_ltmp(1:end-lagtmp);
        [r,p] = corrcoef(tmp1, tmp2, 'Rows', 'complete');
        
    end

    Oris_fitted_md{dd} = tmp1;
    seg_fitted_md{dd} = tmp2;
    dCOM_fitted_md{dd} = tmp3;
end

Oris_mdfit = [];
SI_mdfit = [];
dCOM_mdfit = [];
for dd = 1:size(Oris_fitted_md, 2)
    Oris_mdfit = [Oris_mdfit, Oris_fitted_md{dd}];
    SI_mdfit = [SI_mdfit, seg_fitted_md{dd}];
    dCOM_mdfit = [dCOM_mdfit, dCOM_fitted_md{dd}];
end

vis_N = 500;

subplot(3,4,3)
scatter(Oris_mdfit(1:vis_N), SI_mdfit(1:vis_N), 'k', 'markerEdgeAlpha', 0.5)
idx = (isnan(SI_mdfit));
y = polyfit(Oris_mdfit(~idx), SI_mdfit(~idx), 1);
[r,p] = corrcoef(Oris_mdfit, SI_mdfit,'Rows', 'complete');
hold on
plot([-1, 1],[-1*y(1)+y(2), 1*y(1)+y(2)], 'r')
title(['r = ' num2str(r(2)) ', p = ' num2str(p(2))])
xlim([-1.2 1.2]);xticks([-1 0 1]); ylim([0 0.6]);yticks([0 0.1 0.2 0.3 0.4 0.5 0.6])
ylabel('SI');xlabel('cos(2\Theta + \Phi)')
% data line
idx = (isnan(seg_datafit));
tmp = seg_datafit*nanstd(SI_mdfit)/nanstd(seg_datafit);
seg_datafit_rescale = tmp+nanmean(SI_mdfit)-nanmean(tmp);
y = polyfit(Oris_datafit(~idx), seg_datafit_rescale(~idx), 1);
hold on
%plot([-1, 1],[-1*y(1)+y(2), 1*y(1)+y(2)], 'b')
subplot(3,4,4)
scatter(Oris_mdfit(1:vis_N), dCOM_mdfit(1:vis_N), 'k', 'markerEdgeAlpha', 0.5)
idx = (isnan(dCOM_mdfit));
y = polyfit(Oris_mdfit(~idx), dCOM_mdfit(~idx), 1);
[r,p] = corrcoef(Oris_mdfit, dCOM_mdfit, 'Rows', 'complete');
hold on
plot([-1, 1],[-1*y(1)+y(2), 1*y(1)+y(2)], 'r')
title(['r = ' num2str(r(2)) ', p = ' num2str(p(2))])
xlim([-1.2 1.2]);xticks([-1 0 1]); ylim([0 2]);yticks([0 0.5 1 1.5 2])
% data line
idx = (isnan(dCOM_datafit));
tmp = dCOM_datafit*nanstd(dCOM_mdfit)/nanstd(dCOM_datafit);
dCOM_datafit_rescale = tmp+nanmean(dCOM_mdfit)-nanmean(tmp);
y = polyfit(Oris_datafit(~idx), dCOM_datafit_rescale(~idx), 1);
ylabel('dON-OFF');xlabel('cos(2\Theta + \Phi)')
hold on


%% Getting period distribution from model
iter_crosssec = 10000; % 
%% get cross-section of each map
SI_cs = zeros(NN, iter_crosssec); % NN : # of data points in each penetration
dCOM_cs = zeros(NN, iter_crosssec);
ori_cs = zeros(NN, iter_crosssec);
for ii = 1:iter_crosssec
    
    %% 1. get random two dots
    dist_tmp = 0;
    while dist_tmp<dist_thr
        xtmp1 = randi(xlen);
        ytmp1 = randi(ylen);
        xtmp2 = randi(xlen);
        ytmp2 = randi(ylen);
        dist_tmp = sqrt((xtmp2-xtmp1)^2+(ytmp2-ytmp1)^2);
    end
    
    %     thetatmp = atan((ytmp2-ytmp1)/(xtmp2-xtmp1));
    %     xtmp2 = xtmp1+dist_thr*cos(thetatmp);
    %     ytmp2 = ytmp1+dist_thr*sin(thetatmp);
    xtmp22 = xtmp1+(xtmp2-xtmp1)*dist_thr/dist_tmp;
    ytmp22 = ytmp1+(ytmp2-ytmp1)*dist_thr/dist_tmp;
    xi = [xtmp1, xtmp22];
    yi = [ytmp1, ytmp22];
    %% get profiles for seg, dCOM, ori
    %     [cx, cy, c, xii, yii] = improfile(seg_filtered, xi, yi, 27);
    ntmp = randi(itermap);
    SI_filteredtmp = squeeze(SI_filtered_set(:,:,ntmp));
    dCOM_filteredtmp = squeeze(dCOM_filtered_set(:,:,ntmp));
    ori_filteredtmp = squeeze(ori_filtered_set(:,:,ntmp));
    
    SI_lin = improfile(SI_filteredtmp, xi, yi, NN);
    dCOM_lin = improfile(dCOM_filteredtmp, xi, yi, NN);
    ori_lin = improfile(ori_filteredtmp, xi, yi, NN);
    
    SI_cs(:, ii) = SI_lin;
    dCOM_cs(:, ii) = dCOM_lin;
    ori_cs(:, ii) = ori_lin;
    
end


%% Calculating mean pw diff. for two randomly-choosen cross section
max_dist = 26;
pwdf_seg_mean = zeros(max_dist, iter_crosssec);
pwdf_dCOM_mean = zeros(max_dist, iter_crosssec);
pwdf_ori_mean = zeros(max_dist, iter_crosssec);

for ii = 1:iter_crosssec
    
    pw_seg_diff_tmp = 0;
    pw_dCOM_diff_tmp = 0;
    pw_ori_diff_tmp = 0;
    pw_dist = 0;
    
    segtmp = SI_cs(:, ii);
    dCOMtmp = dCOM_cs(:, ii);
    oritmp = ori_cs(:, ii);
    
    for dist_N = 1:max_dist
        for ij = 1:NN-dist_N
            seg_diff_tmp = abs(segtmp(ij)-segtmp(ij+dist_N));
            pw_seg_diff_tmp = [pw_seg_diff_tmp, seg_diff_tmp];
            
            dCOM_diff_tmp = abs(dCOMtmp(ij)-dCOMtmp(ij+dist_N));
            pw_dCOM_diff_tmp = [pw_dCOM_diff_tmp, dCOM_diff_tmp];
            
            ori_diff_tmp = abs(oritmp(ij)-oritmp(ij+dist_N));
            if  ori_diff_tmp>pi/2
                ori_diff_tmp = pi- ori_diff_tmp;
            end
            pw_ori_diff_tmp = [pw_ori_diff_tmp, ori_diff_tmp];
            
            pw_dist = [pw_dist, dist_N];
        end
    end
    
    for dist_N = 1:max_dist
        indtmp = pw_dist==dist_N;
        pwdf_seg_mean(dist_N, ii) = nanmean(pw_seg_diff_tmp(indtmp));
        pwdf_dCOM_mean(dist_N, ii) = nanmean(pw_dCOM_diff_tmp(indtmp));
        pwdf_ori_mean(dist_N, ii) = nanmean(pw_ori_diff_tmp(indtmp));
    end
end

%% Calculating period distribution from two curves (mimic processing in data)
iterr = 100;
pwdf_seg_mean2 = zeros(max_dist, iterr);
pwdf_dCOM_mean2 = zeros(max_dist, iterr);
pwdf_ori_mean2 = zeros(max_dist, iterr);

for ii = 1:iterr
    indtmp = randi(iter_crosssec, [1,2]);
    pwdf_seg_mean2(:,ii) = mean(pwdf_seg_mean(:,indtmp),2);
    pwdf_dCOM_mean2(:,ii) = mean(pwdf_dCOM_mean(:,indtmp),2);
    pwdf_ori_mean2(:,ii) = mean(pwdf_ori_mean(:,indtmp),2);
end

%% Getting period distribution
seg_period_modelsets = zeros(1,iterr);
dCOM_period_modelsets = zeros(1,iterr);
Ori_period_modelsets = zeros(1,iterr);
norm_dist_data = x_dist/Ori_period_data;
xlimtmp = 1:length(norm_dist_data);
pwdf_sets = pwdf_seg_mean(xlimtmp,:);
for ii = 1:iterr
    tmpprof_seg = pwdf_seg_mean2(xlimtmp, ii);
    TF = islocalmin(tmpprof_seg);
    indtmp = find(TF);
    if isempty(indtmp)
        seg_period_modelsets(ii) = NaN;
    else
        indtmp2 = find(tmpprof_seg(indtmp)==min(tmpprof_seg(indtmp)));
        seg_period_modelsets(ii) = indtmp(indtmp2);
    end

    tmpprof_dCOM = pwdf_dCOM_mean2(xlimtmp, ii);
    TF = islocalmin(tmpprof_dCOM);
    indtmp = find(TF);
    if isempty(indtmp)
        dCOM_period_modelsets(ii) = NaN;
    else
        indtmp2 = find(tmpprof_dCOM(indtmp)==min(tmpprof_dCOM(indtmp)));
        dCOM_period_modelsets(ii) = indtmp(indtmp2);
    end
    
    tmpprof_ori = pwdf_ori_mean2(xlimtmp, ii);
    TF = islocalmin(tmpprof_ori);
    indtmp = find(TF);
    if isempty(indtmp)
        Ori_period_modelsets(ii) = NaN;
    else
        indtmp2 = find(tmpprof_ori(indtmp)==min(tmpprof_ori(indtmp)));
        Ori_period_modelsets(ii) = indtmp(indtmp2);
    end
end

%% Distribution of period values obtained from model


subplot(3,4,12)
bar([1,2,3], 1.1*1/11*[nanmean(seg_period_modelsets), nanmean(Ori_period_modelsets), nanmean(dCOM_period_modelsets)])
hold on;box off;ylabel('Period')
errorbar([1,2,3], 1.1*1/11*[nanmean(seg_period_modelsets), nanmean(Ori_period_modelsets), nanmean(dCOM_period_modelsets)], ...
    1.1*1/11*[nanstd(seg_period_modelsets), nanstd(Ori_period_modelsets), nanstd(dCOM_period_modelsets)], 'LineStyle', 'none')
ylim([0 1.5]);xticks([1,2,3]);xticklabels({'\theta', 'd', 'SI'})
hold on
tmp1 = mean(pwdf_seg_mean2, 2);
tmp2 = mean(pwdf_dCOM_mean2, 2);
tmp3 = mean(pwdf_ori_mean2, 2);
stdseg_model = std(tmp1(1:n_x));
stddCOM_model = std(tmp2(1:n_x));
stdori_model = std(tmp3(1:n_x));

%% Visualizing pairwise difference curves

norm_dist_data = x_dist/Ori_period_data;
xlimtmp = 1:length(norm_dist_data);
% figure
tmp = mean(pwdf_ori_mean2, 2)*180/pi;
tmp1 = std(pwdf_ori_mean2, [],2)*180/pi;
subplot(3,4,7)
shadedErrorBar(norm_dist_data, tmp(xlimtmp), tmp1(xlimtmp), 'lineprops','k')
hold on
% set(gcf,'units','points','position',[500,300,600,400])
set(gca, 'YTick', [0 45 90])
set(gca, 'XTick', [0 0.5 1 1.5])
ylim([0 90])
xlabel('Pairwise cortical distance')
ylabel('Orientation difference')
xlim([0 1.5+0.1])
hold on;plot([1 1] ,[0 90], 'k--');
text(0.7,80,'period')
tmp = mean(pwdf_dCOM_mean2, 2);
tmp1 = std(pwdf_dCOM_mean2, [],2);
norm_tmp2 = max(tmp(xlimtmp))/max(pw_dC_diff_mean(1:n_x));
tmp = tmp/norm_tmp2;
tmp1 = tmp1/norm_tmp2;

subplot(3,4,8)
shadedErrorBar(norm_dist_data, tmp(xlimtmp), tmp1(xlimtmp), 'lineprops','k')
xlabel('Pairwise cortical distance')
ylabel('dON-OFF difference')
ylim([0 0.4])
xlim([0 1.5+0.1])
set(gca, 'XTick', [0 0.5 1 1.5])
set(gca, 'YTick', [0 0.2 0.4])
hold on;plot([1 1] ,[0 0.4], 'k--')

tmp = mean(pwdf_seg_mean2, 2);
tmp1 = std(pwdf_seg_mean2, [],2);
norm_tmp = max(tmp(xlimtmp))/max(pw_seg_diff_mean(1:n_x));
tmp = tmp/norm_tmp;
tmp1 = tmp1/norm_tmp;
subplot(3,4,11)
shadedErrorBar(norm_dist_data, tmp(xlimtmp), tmp1(xlimtmp), 'lineprops','k')
xlim([0 1.5+0.1])
ylim([0 0.3])
xlabel('Pairwise cortical distance')
ylabel('SI difference')
set(gca, 'YTick', [0 0.1 0.2 0.3])
set(gca, 'XTick', [0 0.5 1 1.5])
hold on;plot([1 1] ,[0 0.3], 'k--')
set(gcf, 'units','normalized','outerposition',[0.2 0.05 0.8 0.7])
sgtitle('Simulation results for periodic clustering of simple and complex cells in V1')
