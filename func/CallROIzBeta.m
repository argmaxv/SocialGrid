function [ROIname0, curROI1, indexvoi]=CallROIzBeta(ROIname0)
     for ri=1:numel(ROIname0)
    
        ROIname=ROIname0{ri};

        switch ROIname
            
            case 'EC_xB'
                curROI1{ri}='EC_xB';
                indexvoi{ri}='EC xB Fig3B1';
                
            case 'mPFC_xB'
                curROI1{ri}='mPFC_xB'; 
                indexvoi{ri}='mPFC xB Fig3B2';
                
            case 'EC_xD'
                curROI1{ri}='EC_xB';                
                indexvoi{ri}='EC xD Fig3C1'; 
                
            case 'mPFC_xD'
                curROI1{ri}='mPFC_xB'; 
                indexvoi{ri}='mPFC xD Fig3C2';                

            case 'mPFC_GP'
                curROI1{ri}='mPFC_GP';
                indexvoi{ri}='mPFC GP Fig4B1';

            case 'rTPJ_GP'
                curROI1{ri}='rTPJ_GP';
                indexvoi{ri}='rTPJ GP Fig4B2';
            
            case 'lTPJ_GP'
                curROI1{ri}='lTPJ_GP';
                indexvoi{ri}='lTPJ GP Fig4B3';
            
        end %switch
     end % for
end