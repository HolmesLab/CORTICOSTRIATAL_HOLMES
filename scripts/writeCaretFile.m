function writeCaretFile(annot_file, surf_template, cors_in, out_file)

    % Read file that contains information about the network assignment of
    % the Yeo 57 or 114 region parcellation 
    % ------------------
    [v, L, ct] = read_annotation(annot_file); 
    base_curve = MRIread(surf_template);
    
	%data     = readtable(csvin);
    cor_vals = cors_in;%cellfun(@str2num, data.cor_vals, 'UniformOutput', false);
    vals_new = num2cell(zeros(length(v),1));
    
    for row = 1:length(cor_vals)
        cur_label   = row; 
        
        if isempty(strfind(cor_vals{row}, 'NA')) == 0
            continue
        else
        	val_in = str2num(cor_vals{row}); 
        end

        % if it's an empty correlation value, set to zero
        % --------------
        if isempty(val_in) == 1
            val_in = 0;
        end
        if isnan(val_in) == true;
            continue;
        end
        
        % Caret won't display values if they are too close to zero
        % so set the correlation value to the lowest displayable value
        % ----------
        if val_in < .001 && 0 < val_in
            val_in = .0011;
        elseif val_in > -.001 && 0 > val_in
            disp(val_in)
            val_in = -.0011;
        end
        
        replace_val = ct.table(cur_label+1,5);
        idxs = v(L == replace_val)+1;
        vals_new(idxs) = {val_in};
    end
    base_curve.vol = (cell2mat(vals_new)');
    MRIwrite(base_curve, out_file);
end




