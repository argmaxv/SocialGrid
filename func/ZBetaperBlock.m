function z=ZBetaperBlock(bseries, info)
% Z scored a Beta vector within each block

    z=nan(size(bseries));
    for blx=min(info):max(info)
        x=bseries(info==blx);
        z(info==blx)=(x - nanmean(x))/nanstd(x);        
    end
end