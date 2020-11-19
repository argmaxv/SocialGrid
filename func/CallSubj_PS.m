function [subj, subn] = CallSubj_PS(option)
    
    if nargin<1
        %subj = {'F01','F02','F03','F04','F05','F06','F07','F08','F09','F10','F11','F12','F13','F14','F15','F16','F17','F18','F19','F20','F22','F23','F24','F25'};
        subj = {'F01','F02','F03','F04','F06','F07','F08','F09','F10','F12','F13','F14','F15','F17','F18','F19','F20','F22','F23','F24','F25'};
    elseif option==2
         [~, ~, info, ~, ~]=Call_default_PS;
        prefix{1}=info.prefix.day1;
        prefix{2}=info.prefix.day2;
        subj0 = {'F01','F02','F03','F04','F06','F07','F08','F09','F10','F12','F13','F14','F15','F17','F18','F19','F20','F22','F23','F24','F25'};
        for d=1:2
            for s=1:numel(subj0)
                subj{(d-1)*numel(subj0)+s}=[prefix{d}, subj0{s}];
            end
        end
    end
    subn = numel(subj);
end