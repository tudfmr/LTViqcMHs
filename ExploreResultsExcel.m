clc
clear all

algo=[{'SODE'}
    {'SOJADE'}
    {'SOSHADE'}
    {'SOLSHADE'}
    {'SOSCA'}    %%%5
    {'SODA'}        
    {'SOGWO'}
    {'SOMFO'}
    {'SOWOA'}
    {'SOISCA'}     %%%%10
    {'SOm_SCA'}
    {'SOSinDE'}
    {'SOLSHADEND'}
    {'SOSPS_LSHADE_EIG'}
    {'SOSCAadbL02'}   %%%% 15
    {'SODEobj3'}   
    {'SOJADEobj3'}
    {'SOSHADEobj3'}
    {'SOLSHADEobj3'}   
    {'SOSCAobj3'}   %%%20
    {'SODAobj3'}        
    {'SOGWOobj3'}
    {'SOMFOobj3'}
    {'SOWOAobj3'}
    {'SOSCAadbL02obj3'}  %%%% 25
    {'SOISCAobj3'}    
    {'SOm_SCAobj3'}
    {'SOLSHADENDobj03'}  
    {'SOSPS_LSHADE_EIGobj03'}
    {'SOSinDEobj03'}    %%%%30
];



al=1:size(algo,1);
b=[0.01 0.03 0.05 0.06 0.07 0.075 0.08 0.085];

nrun=5;

aa=['a';'b';'c';'d';'e';'f';'g';'h';'i';'j';'k';'l';'m';...
    'o';'p';'q';'r';'s';'t';'u';'v';'w';'x';'y';'z'];

for k=1:2%length(b)
    ii=2;
    load(['OrignalbNo' num2str(k)])
    fori=fval;
    xlswrite('ResultsLuncher',{'Original'},['b=' num2str(b(k))],['a' num2str(ii)])
    xlswrite('ResultsLuncher',fori,['b=' num2str(b(k))],['b' num2str(ii)])

    for i=1:length(al)
        ii=ii+1;
        fMHi0=zeros(1,nrun);
        for j=1:nrun
            load(['bNo' num2str(k) 'Al' num2str(al(i)) 'Run' num2str(j)])
            if fpmin<50
                fMHi0(1,j)=fpmin;     
                
            end
            clear fpmin
            xlswrite('ResultsLuncher',{['Run' num2str(j)]},['b=' num2str(b(k))],aa(j+1))
        end
        xlswrite('ResultsLuncher',{'mean'},['b=' num2str(b(k))],aa(nrun+2,:))
        xlswrite('ResultsLuncher',{'std'},['b=' num2str(b(k))],aa(nrun+3,:))
        fMHi=fMHi0;
        [rd,cd]=find(fMHi==0);
        fMHi(cd)=[];
        fmean=mean(fMHi);
        fstd=std(fMHi);
        xlswrite('ResultsLuncher',algo(al(i),:),['b=' num2str(b(k))],['a' num2str(ii)])
        
        xlswrite('ResultsLuncher',[fMHi0,fmean,fstd],['b=' num2str(b(k))],['b' num2str(ii) ':' aa(nrun+3) num2str(ii)])
      
        clear fMHi
    end

end








