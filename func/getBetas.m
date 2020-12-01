function Betas=getBetas(datadir, roi, info)
    parfor bidx=1:numel(info)

        % to show the progress
        if mod(bidx,100)==0
            fprintf([num2str(bidx), '\t']);
        elseif bidx==numel(info)
            fprintf('\n');
        end

        if ~strcmp(info(bidx), ' ')
            data=fullfile(datadir, ['beta_', num2str(bidx, '%04d'), '.nii']);
            [~, menaval]=extract_meanactivity_fromROI(roi,data);
            Betas(bidx,1)=menaval;
        else
            Betas(bidx,1)=nan;
        end
    end
end