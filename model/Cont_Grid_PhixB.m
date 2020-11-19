function [spm] = Cont_Grid_PhixB(cont_names, cntmx, spm)
    
% Assigning weights to generate con files and returning 'spm.stats.con.consess'
% spm = Cont_design_name (contrasts names, array of condition and conditionXpmod names, spm (1st level analysis info))
    motion_reg=6;
    for n=1:length(cont_names)
        
            spm.stats.con.consess{n}.tcon.name =cont_names{n};
        
        if strcmp(cont_names{n}, 'vF0102')==1                                                    %  Contrast 1
            cntw=zeros(1, length(cntmx)+motion_reg);
            mf=1; cntw(mf,find(strcmp(cntmx, 'vF0102')))=1;
            spm.stats.con.consess{n}.tcon.weights =cntw;
            spm.stats.con.consess{n}.tcon.sessrep = 'repl'; %replicates over sessions

        elseif strcmp(cont_names{n}, 'GP')==1                                                    %  Contrast 2
            cntw=zeros(1, length(cntmx)+motion_reg);
            mf=1; cntw(mf,find(strcmp(cntmx, 'GP')))=1;
            spm.stats.con.consess{n}.tcon.weights =cntw;
            spm.stats.con.consess{n}.tcon.sessrep = 'repl'; %replicates over sessions        

        elseif strcmp(cont_names{n}, 'Easiness')==1                                           %  Contrast 3
            cntw=zeros(1, length(cntmx)+motion_reg);
            mf=1; cntw(mf,find(strcmp(cntmx, 'Easiness')))=1;
            spm.stats.con.consess{n}.tcon.weights =cntw;
            spm.stats.con.consess{n}.tcon.sessrep = 'repl'; %replicates over sessions        
        
        elseif strcmp(cont_names{n}, 'Euc')==1                                                    %  Contrast 4
            cntw=zeros(1, length(cntmx)+motion_reg);
            mf=1; cntw(mf,find(strcmp(cntmx, 'Euc')))=1;
            spm.stats.con.consess{n}.tcon.weights =cntw;
            spm.stats.con.consess{n}.tcon.sessrep = 'repl'; %replicates over sessions        
            
        else
            cntw=zeros(1, length(cntmx)+motion_reg);
            mf=1;
            cntw(mf,find(strcmp(cntmx, cont_names{n})))=1;
            spm.stats.con.consess{n}.tcon.weights =cntw;
        end
    end
    spm.stats.con.delete = 1;
end