function [  ] = fn_qalas_blochsim_map_v3( qalas_nifti_filename, qalas_json_filename, afi_img_nifti_filename, afi_b1_nifti_filename, output_path, num_b1_bins, inv_eff, estimate_pd_map, use_parfor, num_workers, use_bet )
%FN_QALAS_MAP Summary of this function goes here
%   Detailed explanation goes here


% qalas_nifti_filename: full name (including directory) of qalas nifti
% qalas_json_filename: full name (including directory) of qalas .json file
% afi_img_nifti_filename: full name (including directory) of afi image 
%       (not the b1 map, the tfl image itself, either contrast, to be used 
%       in registration to qalas)
% afi_b1_nifti_filename: full name (including directory) of afi b1 map 
% output_path: directory name to save output maps / temp files
%
%       following are optional parameters:
%
% num_b1_bins: number of bins to use when discretizing b1 map, default=25
% inv_eff: range and step size of inversion efficiency values to use in
%       bloch sim. default = 0.75:0.025:1
% estimate_pd_map: set to 0 to avoid estimating proton density -> speeds up
%       the computation
% use_parfor: set to 1 to use parallel processing
% num_workers: number of parfor workers, recommended 1 worker per 64 GB ram
%   available. e.g. use 4 workers on a workstation with 256 GB ram
%
% example usage:
%
%    qalas_nifti_filename = 'path/qalas_R2_QALAS3d_1p15mm_R2_20220823120716_10.nii';
%    qalas_json_filename = 'path/qalas_R2_QALAS3d_1p15mm_R2_20220823120716_10.json';
%    afi_img_nifti_filename = 'path//b1_afi_AFI3D_RS_Spoilx35_R3_20220823120716_18.nii';
%    afi_b1_nifti_filename = 'path//b1_afi_AFI3D_RS_Spoilx35_R3_20220823120716_21.nii';
%    output_path = 'path/output_path';
%
%    fn_qalas_blochsim_map_v3( qalas_nifti_filename, qalas_json_filename, afi_img_nifti_filename, afi_b1_nifti_filename, output_path );


DEBUG_MODE = 0;     % set to 1 to plot images

addpath(genpath('utils/'))


if nargin < 6
    num_b1_bins = 25;
end

if nargin < 7
    % inversion efficiency
    inv_eff = 0.75:0.025:1;
end

if nargin < 8
    estimate_pd_map = 0;     % set to 1 to estimate PD map (makes it slower)
end

if nargin < 9
    use_parfor = 0;
    % set to 1 to use parfor in parameter fitting with 50% of the total
    % cores -> this may lead to memory issues
    % (parfor is used by default in dictionary generation)
end

if nargin < 10
    num_workers = 4;
end

if nargin < 11
    % set to 0 if this is a phantom experiment
    use_bet = 1;
end
 

% load qalas nifti
dt_img = load_nifti(qalas_nifti_filename);
img = dt_img.vol;

N = s(img(:,:,:,1));

img = single(reshape(img, [N, 5]));


%--------------------------------------------------------------------------
% load json header
%--------------------------------------------------------------------------

json = fileread(qalas_json_filename);

res = jsondecode(json);

% The echo spacing time read from the JSON file might be incorrect on VE11C. 
% On XA30, the echo spacing time read from the JSON file is correct.
esp = res.RepetitionTime            % echo spacing in sec
turbo_factor = res.EchoTrainLength  % ETL



%--------------------------------------------------------------------------
% rsos qalas image using fsl 
%--------------------------------------------------------------------------


system(['fslmaths ', qalas_nifti_filename, ' -sqr ', output_path, 'qalas_rsos'])
system(['fslmaths ', output_path, 'qalas_rsos', ' -Tmean ', output_path, 'qalas_rsos'])
system(['fslmaths ', output_path, 'qalas_rsos', ' -sqrt ', output_path, 'qalas_rsos'])

 

%--------------------------------------------------------------------------
% calculate brain mask after neck removal:
%--------------------------------------------------------------------------

if use_bet
    % remove neck region using robustfov
    system(['robustfov -i ', output_path, 'qalas_rsos -r ', output_path, 'qalas_rsos_robustfov'])

    system(['bet ', output_path, 'qalas_rsos_robustfov ', output_path, 'qalas_rsos_robustfov_bet -m '])

    % register bet mask onto the original FOV without neck cropping
    string_flirt = ['flirt -interp nearestneighbour -dof 6 -cost mutualinfo -in ', output_path, 'qalas_rsos_robustfov_bet_mask -ref ', output_path, 'qalas_rsos -out ', output_path, 'qalas_rsos_robustfov_bet_mask_flirt.nii -omat ', output_path, 'omat_qalas_bet'];


    tic
        system(string_flirt)
    toc

    dt_msk = load_nifti([output_path, 'qalas_rsos_robustfov_bet_mask_flirt.nii.gz']);
    msk = dt_msk.vol ~= 0;

else
    msk = ones(N);
end


if DEBUG_MODE
    dt_rsos = load_nifti([output_path, 'qalas_rsos.nii.gz']);

    imagesc3d2( dt_msk.vol, s(dt_msk.vol)/2, 1, [180,0,90], [0,1]),
    imagesc3d2( dt_rsos.vol, s(dt_rsos.vol)/2, 2, [180,0,90], [0,400]),
    imagesc3d2( dt_rsos.vol .* (dt_msk.vol~=0), s(dt_rsos.vol)/2, 3, [180,0,90], [0,400]),
end



%--------------------------------------------------------------------------
% afi niftis appear to be flipped/rotated wrt qalas images-> Flirt fixes this
% register afi onto rsos qalas -> also interpolates to qalas resolution
%--------------------------------------------------------------------------

string_flirt = ['flirt -interp sinc -dof 6 -cost mutualinfo -in ', afi_img_nifti_filename, ' -ref ', output_path, 'qalas_rsos -out ', output_path, 'afi2qalas.nii -omat ', output_path, 'omat_afi2qalas'];

tic
    system(string_flirt)
toc


%--------------------------------------------------------------------------
% apply estimated interp on B1 map
%--------------------------------------------------------------------------

string_flirt = ['flirt -interp sinc -in ', afi_b1_nifti_filename, ' -ref ', output_path, 'qalas_rsos -applyxfm -init ', output_path, 'omat_afi2qalas -out ', output_path, 'b1_afi_flirt'];

tic
    system(string_flirt)
toc

dt_b1_flirt = load_nifti([output_path, 'b1_afi_flirt.nii.gz']);


if DEBUG_MODE
    dt_qalas_rsos = load_nifti([output_path, 'qalas_rsos.nii.gz']);
    dt_afi2qalas = load_nifti([output_path, 'afi2qalas.nii.gz']);

    imagesc3d2( dt_b1_flirt.vol/60 .* dt_msk.vol, s(dt_b1_flirt.vol)/2, 1, [180,-0,90], [0.65,1.35]), colormap jet
    imagesc3d2( dt_qalas_rsos.vol .* dt_msk.vol, s(dt_qalas_rsos.vol)/2, 2, [180,-0,90], [0,400]),
    imagesc3d2( dt_afi2qalas.vol .* dt_msk.vol, s(dt_qalas_rsos.vol)/2, 3, [180,-0,90], [0,1000]),
end

 

%--------------------------------------------------------------------------
% threshold high and low B1 values: use raw b1 map without polyfit
%--------------------------------------------------------------------------

thre_high = 1.35;
thre_low = 0.65;        

img_b1_load = dt_b1_flirt.vol / 60;

temp = img_b1_load .* msk;

temp(temp > thre_high) = thre_high;
temp(temp < thre_low) = thre_low;

temp = temp .* msk;
img_b1 = temp .* msk;    



%--------------------------------------------------------------------------
% create masks for each b1 value
%--------------------------------------------------------------------------

b1_val = linspace( min(img_b1(msk==1)), max(img_b1(msk==1)), num_b1_bins );

sum_msk = sum(msk(:));

if length(b1_val) == 1
    % do not use b1 correction
    msk_b1 = msk;
else
    
    msk_b1 = zeross([N,length(b1_val)]);
    
    for t = 1:length(b1_val)
        if t > 1
            msk_b1(:,:,:,t) = (img_b1 <= b1_val(t)) .* (img_b1 > b1_val(t-1));
        else
            % first bin: contains all pixels less than b1_val(1)
            msk_b1(:,:,:,t) = msk.*(img_b1 <= b1_val(t));
        end
    
        if t == length(b1_val)
            % last bin: contains all pixels greater than b1_val(end-1)
            msk_b1(:,:,:,t) = img_b1 > b1_val(t-1);
        end
        
        msk_b1(:,:,:,t) = msk_b1(:,:,:,t) .* msk;
    
    end
end




%--------------------------------------------------------------------------
% create look up table 
%--------------------------------------------------------------------------

t1_entries = [5:5:3000, 3100:100:5000];
t2_entries = [1:1:700, 710:20:1000, 1100:100:2500];


T1_entries = repmat(t1_entries.', [1,length(t2_entries)]).';
T1_entries = T1_entries(:);
  

T2_entries = repmat(t2_entries.', [1,length(t1_entries)]);
T2_entries = T2_entries(:);


t1t2_lut = cat(2, T1_entries, T2_entries);


% remove cases where T2>T1
idx = 0;
for t = 1:length(t1t2_lut)
    if t1t2_lut(t,1) < t1t2_lut(t,2)
        idx = idx+1;
    end
end

t1t2_lut_prune = zeross( [length(t1t2_lut) - idx, 2] );

idx = 0;
for t = 1:length(t1t2_lut)
    if t1t2_lut(t,1) >= t1t2_lut(t,2)
        idx = idx+1;
        t1t2_lut_prune(idx,:) = t1t2_lut(t,:);
    end
end

length_dict = length(t1t2_lut_prune);


disp(['dictionary entries: ', num2str(length_dict)])


TR = 4500e-3;       % excluding dead time at the end of the sequence 
alpha_deg = 4;

num_reps = 5;       % number of repetitions to simulate to reach steady state
echo2use = 1;

gap_between_readouts = 900e-3;
time2relax_at_the_end = 0;  % actual TR on the console = TR + time2relax_at_the_end

% simulate Mz and Mxy for T1 T2 values in the look up table


length_b1_val = length(b1_val);
length_inv_eff = length(inv_eff);

signal = zeross([length_dict, 5, length_b1_val, length_inv_eff]);

% if nargin < 10

    delete(gcp('nocreate'))
    c = parcluster('local');    % build the 'local' cluster object

    total_cores = c.NumWorkers; 
    parpool( min(ceil(total_cores*.5), length_b1_val) ) 


    parfor b1 = 1:length_b1_val                    
        for ie = 1:length_inv_eff
            [Mz, Mxy] = sim_qalas_pd_b1_eff_T2(TR, alpha_deg, esp, turbo_factor, t1t2_lut_prune(:,1)*1e-3, t1t2_lut_prune(:,2)*1e-3, num_reps, echo2use, gap_between_readouts, time2relax_at_the_end, b1_val(b1), inv_eff(ie));

            temp = abs( Mxy(:,:,end).' );     % transverse magnitude signal

            % normalize every entry of dict.signal
            for n = 1:size(temp,1)
                temp(n,:) = temp(n,:) / sum(abs(temp(n,:)).^2)^0.5;
            end

            signal(:,:,b1,ie) = temp;
        end    
    end

    delete(gcp('nocreate'))


    %--------------------------------------------------------------------------
    % concat dictionary across inv_eff dimension
    %--------------------------------------------------------------------------


    dict = zeross([length_dict * length_inv_eff, 5, length_b1_val]);

    for t = 1:length_inv_eff
        dict(1 + (t-1)*length_dict : t*length_dict, :, :) = signal(:,:,:,t);
    end

% else
%     disp('load pre-existing dictionary, skip bloch sim')
% 
%     tic
%         load(dict_filename)
%     toc
% end


%--------------------------------------------------------------------------
% dictionary fit -> in each slice, bin voxels based on b1 value
%--------------------------------------------------------------------------


T1_map = zeross(N);
T2_map = zeross(N);

PD_map = zeross(N);     % proton density-> won't be estimated if estimate_pd_map is set to 0
IE_map = zeross(N);     % inversion efficiency


if use_parfor
    
    delete(gcp('nocreate'))

    parpool( num_workers ) 

    parfor slc_select = 1:N(3)
        disp(num2str(slc_select))

        for b = 1:length(b1_val)
            msk_slc = msk_b1(:,:,slc_select,b);
            num_vox = sum(msk_slc(:)~=0);

            if num_vox > 0
                img_slc = zeross([5, num_vox]);

                for t = 1:5
                    temp = sq(img(:,:,slc_select,t));
                    img_slc(t,:) = temp(msk_slc~=0);
                end

                tic   
                    res = dict(:,:,b) * img_slc; % dot product    
    
                    % find T1, T2 values    
                    [~, max_idx] = max(abs(res), [], 1);

                    max_idx_t1t2 = mod(max_idx, length_dict);
                    max_idx_t1t2(max_idx_t1t2==0) = length_dict;

                    res_map = t1t2_lut_prune(max_idx_t1t2,:);

                    % find inv_eff                         
                    max_idx_ie = 1 + (max_idx - max_idx_t1t2) / length_dict;

                    ie_to_use = inv_eff(max_idx_ie);

                    if estimate_pd_map
                        [Mz_sim, Mxy_sim] = sim_qalas_pd_b1_eff_T2(TR, alpha_deg, esp, turbo_factor, res_map(:,1)*1e-3, res_map(:,2)*1e-3, num_reps, echo2use, gap_between_readouts, time2relax_at_the_end, b1_val(b), ie_to_use.');
                    end
                toc

                t1_map = zeross(N(1:2));
                t1_map(msk_slc==1) = res_map(:,1);

                t2_map = zeross(N(1:2));
                t2_map(msk_slc==1) = res_map(:,2);

                ie_map = zeross(N(1:2));
                ie_map(msk_slc==1) = ie_to_use;

                if estimate_pd_map
                    Mxy_sim_use = abs(Mxy_sim(:,:,end));
    
                    scl = zeross([num_vox,1]);
    
                    for idx = 1:size(Mxy_sim_use,2)
                        scl(idx) = Mxy_sim_use(:,idx) \ img_slc(:,idx);
                    end
    
                    pd_map = zeross(N(1:2));
                    pd_map(msk_slc~=0) = scl;     
                    PD_map(:,:,slc_select) = PD_map(:,:,slc_select) + pd_map;
                end

                T1_map(:,:,slc_select) = T1_map(:,:,slc_select) + t1_map;
                T2_map(:,:,slc_select) = T2_map(:,:,slc_select) + t2_map;
                IE_map(:,:,slc_select) = IE_map(:,:,slc_select) + ie_map;

            end
        end
    end
    
    delete(gcp('nocreate'))
    
else
    for slc_select = 1:N(3)
        % don't use parfor 
        disp(num2str(slc_select))

        for b = 1:length(b1_val)
            msk_slc = msk_b1(:,:,slc_select,b);
            num_vox = sum(msk_slc(:)~=0);

            if num_vox > 0
                img_slc = zeross([5, num_vox]);

                for t = 1:5
                    temp = sq(img(:,:,slc_select,t));
                    img_slc(t,:) = temp(msk_slc~=0);
                end

                tic   
                    res = dict(:,:,b) * img_slc; % dot product    

                    % find T1, T2 values    
                    [~, max_idx] = max(abs(res), [], 1);

                    max_idx_t1t2 = mod(max_idx, length_dict);
                    max_idx_t1t2(max_idx_t1t2==0) = length_dict;

                    res_map = t1t2_lut_prune(max_idx_t1t2,:);

                    % find inv_eff                         
                    max_idx_ie = 1 + (max_idx - max_idx_t1t2) / length_dict;

                    ie_to_use = inv_eff(max_idx_ie);
    
                    if estimate_pd_map
                        [Mz_sim, Mxy_sim] = sim_qalas_pd_b1_eff_T2(TR, alpha_deg, esp, turbo_factor, res_map(:,1)*1e-3, res_map(:,2)*1e-3, num_reps, echo2use, gap_between_readouts, time2relax_at_the_end, b1_val(b), ie_to_use.');
                    end
                toc

                t1_map = zeross(N(1:2));
                t1_map(msk_slc==1) = res_map(:,1);

                t2_map = zeross(N(1:2));
                t2_map(msk_slc==1) = res_map(:,2);

                ie_map = zeross(N(1:2));
                ie_map(msk_slc==1) = ie_to_use;

                if estimate_pd_map
                    Mxy_sim_use = abs(Mxy_sim(:,:,end));
    
                    scl = zeross([num_vox,1]);
    
                    for idx = 1:size(Mxy_sim_use,2)
                        scl(idx) = Mxy_sim_use(:,idx) \ img_slc(:,idx);
                    end
    
                    pd_map = zeross(N(1:2));
                    pd_map(msk_slc~=0) = scl;     
                    PD_map(:,:,slc_select) = PD_map(:,:,slc_select) + pd_map;
                end

                T1_map(:,:,slc_select) = T1_map(:,:,slc_select) + t1_map;
                T2_map(:,:,slc_select) = T2_map(:,:,slc_select) + t2_map;
                IE_map(:,:,slc_select) = IE_map(:,:,slc_select) + ie_map;
            end
        end
    end
end



%--------------------------------------------------------------------------
% save results
%--------------------------------------------------------------------------

% save as .mat
save([output_path, 'T1_map.mat'], 'T1_map')
save([output_path, 'T2_map.mat'], 'T2_map')
save([output_path, 'IE_map.mat'], 'IE_map')
save([output_path, 'PD_map.mat'], 'PD_map')


% load the qalas nifti to use its header
hdr = load_nifti([output_path, 'qalas_rsos.nii.gz']);

hdr.datatype = 16;
hdr.bitpix = 32;


hdr.vol = T1_map;      
save_nifti(hdr, [output_path, 'T1_map.nii'])

hdr.vol = T2_map;      
save_nifti(hdr, [output_path, 'T2_map.nii'])

hdr.vol = PD_map;      
save_nifti(hdr, [output_path, 'PD_map.nii'])

hdr.vol = IE_map;      
save_nifti(hdr, [output_path, 'IE_map.nii'])


end

