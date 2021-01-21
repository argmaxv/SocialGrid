function [spm] = Cont_GPEucF1F2_OnOff(cont_names, cntmx, spm)
    
% Assigning weights to generate con files and returning 'spm.stats.con.consess'
% spm = Cont_design_name (contrasts names, array of condition and conditionXpmod names, spm (1st level analysis info))

    motion_reg=6;
    
    for n=1:length(cont_names)
        spm.stats.con.consess{n}.tcon.name =cont_names{n};
        
        if strcmp(cont_names{n}, 'F12on')==1 
            cntw=zeros(1, length(cntmx)+motion_reg);
            mf=1; cntw(mf,find(strcmp(cntmx, 'F12on')))=1;
            spm.stats.con.consess{n}.tcon.weights =cntw; 
%         
        elseif strcmp(cont_names{n}, 'F12off')==1 
            cntw=zeros(1, length(cntmx)+motion_reg);
            mf=1; cntw(mf,find(strcmp(cntmx, 'F12off')))=1;
            spm.stats.con.consess{n}.tcon.weights =cntw; 
% 
        elseif strcmp(cont_names{n}, 'F12onoff_diff')==1 
            cntw=zeros(1, length(cntmx)+motion_reg);
            mf=1; cntw(mf,find(strcmp(cntmx, 'F12on')))=.5;
            mf=1; cntw(mf,find(strcmp(cntmx, 'F12off')))=-.5;
            spm.stats.con.consess{n}.tcon.weights =cntw; 
%             
        elseif strcmp(cont_names{n}, 'GPoff')==1  
            cntw=zeros(1, length(cntmx)+motion_reg);
            mf=1; cntw(mf,find(strcmp(cntmx, 'GPoff')))=1;
            spm.stats.con.consess{n}.tcon.weights =cntw; 
%         
        elseif strcmp(cont_names{n}, 'GPon')==1 
            cntw=zeros(1, length(cntmx)+motion_reg);
            mf=1; cntw(mf,find(strcmp(cntmx, 'GPon')))=1;
            spm.stats.con.consess{n}.tcon.weights =cntw; 
        
        elseif strcmp(cont_names{n}, 'GPonoff_diff')==1
            cntw=zeros(1, length(cntmx)+motion_reg);
            mf=1; cntw(mf,find(strcmp(cntmx, 'GPon')))=.5;
            mf=1; cntw(mf,find(strcmp(cntmx, 'GPoff')))=-.5;
            spm.stats.con.consess{n}.tcon.weights =cntw; 
            
         else
             cntw=zeros(1, length(cntmx)+motion_reg+1);
             mf=1; cntw(mf,find(strcmp(cntmx, cont_names{n})))=1;
             spm.stats.con.consess{n}.tcon.weights =cntw;
        end
        spm.stats.con.consess{n}.tcon.sessrep = 'repl'; %replicates over sessions
    end
    spm.stats.con.delete = 1;
end