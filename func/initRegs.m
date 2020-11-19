function T=initRegs(reg)
    for r=1:numel(reg)
        T.(reg{r}).GP1=[]; % GP_F01
        T.(reg{r}).GP2=[]; % GP_F02
        T.(reg{r}).GrF01=[]; % Gr_F01
        T.(reg{r}).GrF02=[]; % Gr_F02
        T.(reg{r}).EcF01=[]; % Ec_F01
        T.(reg{r}).EcF02=[]; % Ec_F02
    end
end