function IQClist = fct_build_IQClist_public(IQC)

NrNB    = length(IQC.M); % # of the analyzed norm bounds b of the uncertainty
NrIQC   = length(IQC.Psi);

for i = 1 : 1 : NrNB
    
    for j = 1 : 1 : NrIQC
        
        IQClist(i).IQC{j}   = {IQC.Psi{j}, IQC.M{i}};
        
        
    end
    
end

end
