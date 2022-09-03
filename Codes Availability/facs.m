clc

%choose where to save the images:                                                      
cd('E:\master\facs_mona\majority figures\colony 4')
                                                                                                    
%enter file loaction:                   
filelocation='E:\master\facs_mona\Data Majority 30112017\Data Majority 30112018\Colony 4 30.11.2018';     

%enter filename (without the well labels):
filenamee='02-Tube-';

ss=strcat(filelocation,'\',filenamee);

%enter well labels:
c={ '1' '2' '3' '4' };                 %columns labels                                 
r=[ 'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'];  %rows labels                                   
                                                                            
%set the wells contents and inducer name:
    %the inducer that changes from row to row column to column                              
        %inducer name
            name_Inducer_C='Arabinose';     
        %inducer concentraion
            Inducer_C=[0.25 0.125 0.0625 0.03125 ]; 
        %concentraions units
            Inducer_CU=' mM';                       

    
    %the inducer that changes from  from row to row
        %inducer name
            name_Inducer_R='ATC';
        %inducer concentration
            Inducer_R={'[0.0.0]' '[0.0.1]' '[0.1.0]' '[0.1.1]' '[1.0.0]' '[1.0.1]' '[1.1.0]' '[1.1.1]'};
        %concentrarion units
            Inducer_RU={''};    
%    
Inducer_C=[num2str(Inducer_C') repmat(Inducer_CU,length(Inducer_C),1)];
Inducer_CC=cellstr((Inducer_C));                                                                                                    
Inducer_R=[cell2mat((Inducer_R')) cell2mat(repmat(Inducer_RU,length(Inducer_R),1))];
Inducer_RR=cellstr((Inducer_R));
 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for z=[1 2]
    if(z==2)
        x=r;
        r=c;
        c=x;
        x=Inducer_R;
        Inducer_R=Inducer_C;
        Inducer_C=x;
        x=Inducer_RR;
        Inducer_RR=Inducer_CC;
        Inducer_CC=x;
        x=name_Inducer_R;
        name_Inducer_R=name_Inducer_C;
        name_Inducer_C=x;
    end
%open the .fcs file for a specific well,... here i used a function that i
    %downloaded from mathworks
    for i=1:length(r)
            for j=1:length(c)
                if(z==1)
%                     s
                s=strcat(ss,r(i),c(j),'.fcs');
                end
                if(z==2)
                 s=strcat(ss,c(j),r(i),'.fcs');
                end
                s=cell2mat(s);
                [fcsdat, fcshdr, fcsdatscaled] = fca_readfcs(s);
                %take histogram values
                %choose which data you wish to plot (a column vector in fcsdat)
%                 [a,b]=hist(fcsdat(:,7),10000);
%                 b=b;
%                 b=b;
                a=sort(fcsdat(:,6));
                b=logspace(0,log10(max(a)),400);
                aa=zeros(1,length(b));
                for h=2:length(b)
                    aa(h)=sum((a>b(h-1)).*(a<b(h)));
                end
                
av=aa;
%                 % remove negative values
%                 a=nonzeros((b>0).*(a+eps));
%                 b=nonzeros((b>0).*b);

%                 smoothing
%                     averaging
%                     N=10;
%                     av=zeros(length(a)-N+1,1);
%                     for k=1:N
%                     av=av+a(k:end-N+k);
%                     end
%                     av=av/N;
%                     av=av/max(av);
%                     hold all
%                     aa=av;
% av=a;
                    %loess smooth
                    aa=smooth(b(1:length(av)),av,0.1,'loess');
                    aa=aa/max(aa);
% aa=a/max(a);
                %plotting
                hf=figure(2);
                hold all
                plot(b(1:length(aa)),aa,'linewidth',2)
            end
            %after this loop has ended, we have multiple plots for a constant
            %column and varying rows (ot vice reversa),
            %now we set the figure's properties 
            set(gca, 'XScale', 'log')
            set(hf,'position',[680 558 760 420])%set figure size
            
            sss=[Inducer_RR(i) '  '  name_Inducer_R];
            sss=cell2mat(sss);
            title(sss,'FontSize',10)
            ylabel('normalized count','FontSize',14)
            xlabel('Flouresence [a.u.]','FontSize',14)
            xlim([100 7000000])
            ylim([-0.1 1.1])
                %these 2 loops here are important for the LEGEND:
            legend(Inducer_CC,'location','northeastoutside','fontsize',12)
            %save figure
            figname=strcat(r(i),'.emf');
            if isequal(class(figname),class({1}))
                figname=cell2mat(figname);
            end
            saveas(hf,figname)
            close all
    end
end