% -------------------------------
% -------------------------------
% -------------------------------
% Code to plot correlation values on caret surfaces
% 
% Gene expression links functional networks across cortex and striatum
% Anderson, K.M., Krienen, F.M., Choi, E.Y., Reinen, J.M., Yeo, B.T., Holmes. A.J.
%
%
% Written by: Kevin M. Anderson
% Contact:    kevin.anderson@yale.edu 
% -------------------------------
% -------------------------------
% -------------------------------

% this code requires:
% spherical demons: https://sites.google.com/site/yeoyeo02/software/sphericaldemonsrelease

% Set up  base directories 
% ----------------
base_dir         = '/Users/kevinanderson/PHD/PROJECTS/2017_CORTICOSTRIATAL_NATCOMM';
fsavg5_label_dir = [base_dir '/data/Yeo_JNeurophysiol11_SplitLabels/fsaverage5/label/'];


% Required: Install freesurfer matlab tools before running this script
% ----------------
% Add your own local freesurfer installation path here
addpath(genpath('/Applications/freesurfer/'))
addpath(genpath([base_dir '/scripts']))
setenv('FREESURFER_HOME', '/Applications/freesurfer')


% Supp Figure S1: Create the 'missing data' caret surfaces
% these surfaces correspond to the 'white' regions lacking experession data
% ------------------
bihemi_in = readtable([base_dir '/output_files/bihemi_regions_59parcels.csv'], 'Delimiter', ',', 'ReadVariableNames', true, 'ReadRowNames', true); 
hemis = {'lh', 'rh'};
for i = 1:2;
    hemi=hemis{i};
    surf_template    = [base_dir '/data/Yeo_JNeurophysiol11_SplitLabels/fsaverage5/label/' hemi '.Yeo2011_7NetworksConfidence_N1000.mgz'];
	annot_file = [fsavg5_label_dir  hemi '.Yeo2011_17Networks_N1000.split_components.annot'];
    name_arr   = {'Limbic', 'Default', 'Cont', 'VentAttn', 'SomMot'};
    
    cors_in    = bihemi_in.V1;
    out_file   = [base_dir '/caret_files/mrna_bihemi_missing_net17_' hemi '_n6.nii.gz'];
    writeCaretFile(annot_file, surf_template, cellfun(@num2str, num2cell(cors_in), 'UniformOutput', false), out_file)
    caret_dir  = [base_dir '/caret_files']; 
    shape_file = [base_dir '/caret_files/mrna_bihemi_missing_net17_' hemi '_n6.surface_shape'];
    WriteFScurvToCaretSurfaceShape(caret_dir, hemi, out_file, 'fsaverage5', shape_file, 0)
end


% Read mrna and rs-fcMRI correlation values (cortex to striatum)
% ----------------
mrna_in = readtable([base_dir '/caret_files/PEARSON_lhrh_n6_mrna_corrvals.csv'], 'Delimiter', ',', 'ReadVariableNames', true, 'ReadRowNames', true); 
rsfc_in = readtable([base_dir '/caret_files/lhrh_n6_net17_rsfcmri_corrvals.csv'], 'Delimiter', ',', 'ReadVariableNames', true, 'ReadRowNames', true); 


% ------------------
% FIGURE 4 - corticostriatal correlations
% ------------------
% For each striatal subregions, write a .nii file with the
% correlation values plugged in, then convert to caret surface file
% ------------------
hemis = {'lh', 'rh'};
for i = 1:2;
    hemi          = hemis{i};
    surf_template = [base_dir '/data/Yeo_JNeurophysiol11_SplitLabels/fsaverage5/label/' hemi '.Yeo2011_7NetworksConfidence_N1000.mgz'];
	annot_file    = [fsavg5_label_dir  hemi '.Yeo2011_17Networks_N1000.split_components.annot'];
    name_arr      = {'Limbic', 'Default', 'Cont', 'VentAttn', 'SomMot'};
    for t = 1:length(name_arr);
        type     = name_arr{t};
        disp([type ' : ' hemi]);
        cors_in  = getfield(mrna_in, type);
        out_file = [base_dir '/caret_files/mrna_' type '_net7_' hemi '_n6.nii.gz'];
        writeCaretFile(annot_file, surf_template, cors_in, out_file)
        
        caret_dir  = [base_dir '/caret_files']; 
        shape_file = [base_dir '/caret_files/PEARSON_mrna_' type '_net7_' hemi '_n6.surface_shape'];

        setenv('FREESURFER_HOME', '/Applications/freesurfer')
        WriteFScurvToCaretSurfaceShape(caret_dir, hemi, out_file, 'fsaverage5', shape_file, 0)
    end
end



% ------------------
% FIGURE 2 - cortico-cortical correlations
% ------------------
regions = {'Limbic_OFC', 'SomMotA'};
for k = 1:2;
    seed_reg = regions{k};
    cort_mrna_in = readtable([base_dir '/caret_files/mRNA_avg_' seed_reg '_cortical_corrs.csv'], 'Delimiter', ',', 'ReadVariableNames', false, 'ReadRowNames', true); 
    %cort_mrna_in = readtable([base_dir '/caret_files/rsfcmri_avg_' seed_reg '_cortical_corrs.csv'], 'Delimiter', ',', 'ReadVariableNames', false, 'ReadRowNames', true); 
    hemis = {'lh', 'rh'};
    for i = 1:2;
        hemi=hemis{i};
        surf_template    = [base_dir '/data/Yeo_JNeurophysiol11_SplitLabels/fsaverage5/label/' hemi '.Yeo2011_7NetworksConfidence_N1000.mgz'];
        annot_file = [fsavg5_label_dir  hemi '.Yeo2011_17Networks_N1000.split_components.annot'];
        cors_in = cort_mrna_in.Var1;  
        out_file = [base_dir '/caret_files/' seed_reg '_' hemi '_cortex_seed_mrna.nii.gz'];
        %out_file = [base_dir '/caret_files/' seed_reg '_' hemi '_cortex_seed_rsfcmri.nii.gz'];
        writeCaretFile(annot_file, surf_template, cors_in, out_file)

        caret_dir  = [base_dir '/caret_files']; 
        shape_file = [base_dir '/caret_files/' seed_reg '_' hemi '_cortex_seed_mrna_n6.surface_shape'];
        %shape_file = [base_dir '/caret_files/' seed_reg '_' hemi '_cortex_seed_rsfcmri_n6.surface_shape'];
        WriteFScurvToCaretSurfaceShape(caret_dir, hemi, out_file, 'fsaverage5', shape_file, 0)
    end
end


regions = {'Limbic_OFC', 'SomMotA'};
for k = 1:2;
    seed_reg = regions{k};
    cor_file = [base_dir '/caret_files/rsfcmri_avg_' seed_reg '_cortical_corrs.csv'];
    cort_cors = readtable(cor_file, 'ReadVariableNames', false);
    hemis = {'lh', 'rh'};
    for i = 1:2;
        hemi=hemis{i};
        surf_template    = [base_dir '/data/Yeo_JNeurophysiol11_SplitLabels/fsaverage5/label/' hemi '.Yeo2011_7NetworksConfidence_N1000.mgz'];
        annot_file = [fsavg5_label_dir  hemi '.Yeo2011_17Networks_N1000.split_components.annot'];
        cors_in = cort_cors.Var2; 
        out_file = [base_dir '/caret_files/' seed_reg '_' hemi 'cortex_seed_rsfcmri.nii.gz'];
        writeCaretFile(annot_file, surf_template, cors_in, out_file)

        caret_dir  = [base_dir '/caret_files']; 
        shape_file = [base_dir '/caret_files/' seed_reg '_' hemi 'cortex_seed_rsfcmri_n6.surface_shape'];
        WriteFScurvToCaretSurfaceShape(caret_dir, hemi, out_file, 'fsaverage5', shape_file, 0)

        idxs_w_values = cellfun(@isempty, strfind(cors_in, 'NA'));
        tmp_array = ones(length(idxs_w_values),1);
        tmp_array(idxs_w_values) = 0;
        tmp_array = num2str(tmp_array);
        tmp_array = num2cell(tmp_array);
        out_file = [base_dir '/caret_files/' seed_reg '_' hemi '_cortex_MISSING_seed_rsfcmri_n6.nii.gz'];
        tmp_array(idxs_w_values) = {'NA'};
        writeCaretFile(annot_file, surf_template, tmp_array, out_file)

        caret_dir  = [base_dir '/caret_files']; 
        shape_file = [base_dir '/caret_files/' seed_reg '_' hemi '_cortex_MISSING_seed_rsfcmri_n6.surface_shape'];
        WriteFScurvToCaretSurfaceShape(caret_dir, hemi, out_file, 'fsaverage5', shape_file, 0)
    end
end


% Do the same as above for the rsfc values
% ------------------
hemis = {'lh', 'rh'};
for i = 1:2;
    hemi=hemis{i};
    surf_template    = [base_dir '/data/Yeo_JNeurophysiol11_SplitLabels/fsaverage5/label/' hemi '.Yeo2011_7NetworksConfidence_N1000.mgz'];
	annot_file = [fsavg5_label_dir  hemi '.Yeo2011_17Networks_N1000.split_components.annot'];
    hemi=hemis{i};
    name_arr   = {'Limbic', 'Default', 'Cont', 'VentAttn', 'SomMot'};
    for t = 1:length(name_arr);
        type     = name_arr{t};
        cors_in  = getfield(rsfc_in, type);
        out_file = [base_dir '/caret_files/rsfc_' type '_net7_' hemi '_n6.nii.gz'];
        writeCaretFile(annot_file, surf_template, cellfun(@num2str, num2cell(cors_in), 'UniformOutput', false), out_file)
                
        caret_dir  = [base_dir '/caret_files']; 
        shape_file = [base_dir '/caret_files/rsfc_' type '_net7_' hemi '_n6.surface_shape'];
        WriteFScurvToCaretSurfaceShape(caret_dir, hemi, out_file, 'fsaverage5', shape_file, 0)
    end
end


% Make a file showing the missing regions. 
% ------------------
hemis = {'lh', 'rh'};
for i = 1:2;
    hemi=hemis{i};
    cors_in  = getfield(mrna_in, type);
    surf_template    = [base_dir '/data/Yeo_JNeurophysiol11_SplitLabels/fsaverage5/label/' hemi '.Yeo2011_7NetworksConfidence_N1000.mgz'];
	annot_file = [fsavg5_label_dir  hemi '.Yeo2011_17Networks_N1000.split_components.annot'];
    idxs_w_values = cellfun(@isempty, strfind(cors_in, 'NA'));
    tmp_array = ones(length(idxs_w_values),1);
    tmp_array(idxs_w_values) = 0;
    tmp_array = num2str(tmp_array);
    tmp_array = num2cell(tmp_array);
    out_file = [base_dir '/caret_files/mrna_missing_net7_' hemi '_n6.nii.gz'];
    tmp_array(idxs_w_values) = {'NA'};
    writeCaretFile(annot_file, surf_template, tmp_array, out_file)
    
    caret_dir  = [base_dir '/caret_files']; 
    shape_file = [base_dir '/caret_files/mrna_missing_net7_' hemi '_n6.surface_shape'];
    WriteFScurvToCaretSurfaceShape(caret_dir, hemi, out_file, 'fsaverage5', shape_file, 0)
end


% For each individual donors, Iterate over the striatal subregions, write a .nii file with the
% correlation values plugged in, then convert to caret .surface_shape file.
% ------------------
donor_nums = {'9861', '10021', '12876', '14380', '15496', '15697'};
for d = 1:6 
    donor   = donor_nums{d};
    mrna_in = readtable([base_dir '/caret_files/lhrh_' donor '_n6_mrna_corrvals.csv'], 'Delimiter', ',', 'ReadVariableNames', true, 'ReadRowNames', true); 

    hemis = {'lh', 'rh'};
    for i = 1;
        hemi = hemis{i};
        surf_template = [base_dir '/data/Yeo_JNeurophysiol11_SplitLabels/fsaverage5/label/' hemi '.Yeo2011_17NetworksConfidence_N1000.mgz'];
        annot_file    = [fsavg5_label_dir  hemi '.Yeo2011_17Networks_N1000.split_components.annot'];
        name_arr      = {'Limbic', 'SomMot'};
        for t = 1:length(name_arr);
            type     = name_arr{t};
            cors_in  = getfield(mrna_in, type);
            out_file = [base_dir '/caret_files/mrna_' donor '_' type '_net7_' hemi '_n6.nii.gz'];
            writeCaretFile(annot_file, surf_template, cors_in, out_file)

            caret_dir  = [base_dir '/caret_files']; 
            shape_file = [base_dir '/caret_files/mrna_' donor '_' type '_net7_' hemi '_n6.surface_shape'];
            WriteFScurvToCaretSurfaceShape(caret_dir, hemi, out_file, 'fsaverage5', shape_file, 0)
        end
    end
end


% For each donor, create surface files showing regions with missing data
% ------------------
donor_nums = {'9861', '10021', '12876', '14380', '15496', '15697'};
for d = 1:6 
    donor = donor_nums{d};
    mrna_in = readtable([base_dir '/caret_files/lhrh_' donor '_n6_mrna_corrvals.csv'], 'Delimiter', ',', 'ReadVariableNames', true, 'ReadRowNames', true); 

    hemis = {'lh', 'rh'};
    for i = 1;
        hemi=hemis{i};
        cors_in  = getfield(mrna_in, type);
        surf_template    = [base_dir '/data/Yeo_JNeurophysiol11_SplitLabels/fsaverage5/label/' hemi '.Yeo2011_7NetworksConfidence_N1000.mgz'];
        annot_file = [fsavg5_label_dir  hemi '.Yeo2011_17Networks_N1000.split_components.annot'];
        idxs_w_values = cellfun(@isempty, strfind(cors_in, 'NA'));
        tmp_array = ones(length(idxs_w_values),1);
        tmp_array(idxs_w_values) = 0;
        tmp_array = num2str(tmp_array);
        tmp_array = num2cell(tmp_array);
        out_file = [base_dir '/caret_files/' donor '_' hemi '_mrna_missing_net7.nii.gz'];
        tmp_array(idxs_w_values) = {'NA'};
        writeCaretFile(annot_file, surf_template, tmp_array, out_file)

        caret_dir  = [base_dir '/caret_files']; 
        shape_file = [base_dir '/caret_files/' donor '_' hemi '_mrna_missing_net7.surface_shape'];
        WriteFScurvToCaretSurfaceShape(caret_dir, hemi, out_file, 'fsaverage5', shape_file, 0)
    end
end








