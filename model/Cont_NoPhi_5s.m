function [spm] = Cont_NoPhi_5s(cont_names, cntmx, spm) %[spm] = Cont_EF01_EF02cF12DecVec_Stck(cont_names, cntmx, spm)
    
% Assigning weights to generate con files and returning 'spm.stats.con.consess'
% spm = Cont_design_name (contrasts names, array of condition and conditionXpmod names, spm (1st level analysis info))

    motion_reg=6;
    
    for n=1:length(cont_names)
        if strcmp(cont_names{n}(1:5),'ftest')
            contype='fcon';
            spm.stats.con.consess{n}.(contype).name =cont_names{n}(7:end);
        else
            contype='tcon';
            spm.stats.con.consess{n}.(contype).name =cont_names{n};
        end
        
        if strcmp(cont_names{n}, 'CvF01F02')==1  
            cntw=zeros(1, length(cntmx)+motion_reg+1);
            mf=1; cntw(mf, find(strcmp(cntmx, 'CvF01F02')))=1;
            spm.stats.con.consess{n}.tcon.weights =cntw;        
            
        elseif strcmp(cont_names{n}, 'SvF01F02')==1
            cntw=zeros(1, length(cntmx)+motion_reg+1);
            mf=1; cntw(mf, find(strcmp(cntmx, 'SvF01F02')))=1;
            spm.stats.con.consess{n}.tcon.weights =cntw; 
        
        elseif strcmp(cont_names{n}, 'ftest_CvSvF01F02')==1
            cntw=zeros(2, length(cntmx)+motion_reg+1);
            mf=1; cntw(mf, find(strcmp(cntmx, 'CvF01F02')))=1;
            mf=mf+1; cntw(mf, find(strcmp(cntmx, 'SvF01F02')))=1;
            spm.stats.con.consess{n}.fcon.weights =cntw; 
        
        else
             cntw=zeros(1, length(cntmx)+motion_reg+1);
             mf=1; cntw(mf, find(strcmp(cntmx, cont_names{n})))=1;
             spm.stats.con.consess{n}.tcon.weights =cntw;
        end
        spm.stats.con.consess{n}.(contype).sessrep = 'repl'; %replicates over sessions
    end
    spm.stats.con.delete = 1;
end