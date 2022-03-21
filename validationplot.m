function validationplot(cb,cf,ttype,fs)

Lf = cf.Phi;          % |DTF|^2 (N x N x freq)
L_sf = cf.Phi_th;     % Significant level
ttype=lower(ttype);
switch ttype
    case 'dtf'
        Lb=cb.dtf;
        L_sb=cb.dtf_th;
    case 'pdc'
        Lb=cb.pdc;
        L_sb=cb.pdc_th;
end


LTraf = cf.Threshold;     % Threshold 
L2vinff = cf.CIupperbound;  % upperbound
L2vsupf = cf.CIlowerbound;  % lower bound
LTrab = cb.th;
L2vinfb = cb.ic1;
L2vsupb = cb.ic2;


metric = cb.metric;   % Valid options: "euc", "diag" or "info"
[N,~,nFreqs]=size(Lf); % N = q = Number of channels/time series;
%                     % nFreqs = number of points on frequency scale
w = 0:fs/(2*nFreqs):fs/2-fs/(2*nFreqs);
w = w/fs;
web=1:5:length(w); % erro bar line
web2=1:3:length(w); % erro bar line

Lf2=Lf;
Lf2(~isnan(L_sf))=nan;

Lvinfb=L2vinfb;
Lvinfb(isnan(L_sb))=nan;
Lvinfb(isnan(Lvinfb))=Lb(isnan(Lvinfb));

Lvsupb=L2vsupb;
Lvsupb(isnan(L_sb))=nan;
Lvsupb(isnan(Lvsupb))=Lb(isnan(Lvsupb));


Lvinff=L2vinff;
Lvinff(isnan(L_sf))=nan;
% Lvinff(isnan(Lvinff))=Lf(isnan(Lvinff));

Lvsupf=L2vsupf;
Lvsupf(isnan(L_sf))=nan;
% Lvsupf(isnan(Lvsupf))=Lf(isnan(Lvsupf));


% % figure()
% % for i=1:N
% %     for j=1:N
% %         if i~=j
% %     subplot(N,N,(i-1)*N+j)
% %     hold on
% % 
% %     h12 = plot(w,squeeze(Lb(i,j,:)),':','Color','m','LineWidth',4);
% %     hsb = plot(w,squeeze(L_sb(i,j,:)),'Color','m','LineWidth',1.8);
% %     
% %     h12b = plot(w(web2),squeeze(Lf2(i,j,web2)),'>','Color','g','LineWidth',2,'MarkerSize',5); 
% %     hsf = plot(w(web2),squeeze(L_sf(i,j,web2)),'d','Color','g','LineWidth',2,'MarkerSize',5);
% % 
% %     CI  = [squeeze(Lvinfb(i,j,:)),squeeze(Lvsupb(i,j,:))];
% %     
% %        
% % %     hEb = boundedline([web(1),web(end)],[Lf(i,j,web(1)),Lf(i,j,web(end))],[CI(web(1),:);CI(web(end),:)],'-k', 'alpha');    
% %     hEb=patch([w,w(end:-1:1)],[CI(:,1);CI(end:-1:1,2)]',[0.68,0.6,0.6],'EdgeColor','none','FaceAlpha',0.1);
% %     LBF=squeeze(Lf(i,j,web))-squeeze(Lvinff(i,j,web));
% %     UBF=squeeze(Lf(i,j,web))-squeeze(Lvsupf(i,j,web));
% %     
% %     hEf = errorbar(w(web),squeeze(Lf(i,j,web)),LBF,UBF,'.','color',[0.93,0.69,.13],'LineWidth',1);
% %     
% %     htf=plot(w,squeeze(LTrab(i,j,:)),'--','Color','k','LineWidth',1);
% %     htb=plot(w(web),squeeze(LTraf(i,j,web)),'^','Color','r','LineWidth',1.8,'MarkerSize',5);
% %     xlim([0,w(end)]);box on;set(gcf,'color','white');
% %         end
% % 
% %     end
% % end
% % 
% % 


figure()
for i=1:N
    for j=1:N
        if i~=j
    subplot2(N,N,(i-1)*N+j)
    hold on

    CI  = [squeeze(Lvinfb(i,j,:)),squeeze(Lvsupb(i,j,:))]; 
    %---------------------------------------------------------------------%
% % %     cp=find(~isnan(CI(:,1)));cp2=find(diff([0;cp;32])>1);
% % %     if ~isempty(cp)
% % %     if length(cp2)==1
% % %         wu=w(cp(1):cp2-1);
% % %         CIu=CI(cp(1):cp2-1,:);
% % %         Scop=1;
% % %         if cp2==1
% % %         wu=w(cp(1):nFreqs);
% % %         CIu=CI(cp(1):nFreqs,:);
% % %         Scop=1;
% % %         end
% % %     elseif length(cp2)==2
% % %         if cp2(1)==1
% % %         wu=w(cp(1):cp(2));
% % %         CIu=CI(cp(1):cp(2),:);
% % %         Scop=1;
% % %         else
% % %         wu1=w(cp(1):cp(cp2(1)-1));
% % %         CIu1=CI(cp(1):cp(cp2(1)-1),:);
% % %         wu2=w(cp(cp2(1)):cp(cp2(2)-1));
% % %         CIu2=CI(cp(cp2(1)):cp(cp2(2)-1),:);
% % %         Scop=2;
% % %         end
% % %     end
% % %     else
% % %         wu=w(cp);
% % %         CIu=CI(cp,:);
% % %         Scop=1;
% % %     end
% % %     
% % %     if Scop==1   
% % %     patch([wu,wu(end:-1:1)],[CIu(:,1);CIu(end:-1:1,2)]',[0.15,0.15,0.15],'EdgeColor',[0.82,0.75,0.75],'FaceAlpha',0.1);
% % %     elseif Scop==2
% % %     patch([wu1,wu1(end:-1:1)],[CIu1(:,1);CIu1(end:-1:1,2)]',[0.15,0.15,0.15],'EdgeColor',[0.82,0.75,0.75],'FaceAlpha',0.1);
% % %     patch([wu2,wu2(end:-1:1)],[CIu2(:,1);CIu2(end:-1:1,2)]',[0.15,0.15,0.15],'EdgeColor',[0.82,0.75,0.75],'FaceAlpha',0.1);  
% % %     end

    
%     for l=2:length(cp2)
%         wu=w(cp(cp2(l-1)):cp(cp2(l)-1));
%         CIu=CI(cp(cp2(l-1)):cp(cp2(l)-1),:);
%         patch([wu,wu(end:-1:1)],[CIu(:,1);CIu(end:-1:1,2)]',[0.15,0.15,0.15],'EdgeColor',[0.82,0.75,0.75],'FaceAlpha',0.1);
%     end
%     if length(cp2)==2 
%         wu1=w(cp(cp2(1)):cp(cp2(2)-1));
%         CIu1=CI(cp(cp2(1)):cp(cp2(2)-1),:);
%         patch([wu1,wu1(end:-1:1)],[CIu1(:,1);CIu1(end:-1:1,2)]',[0.15,0.15,0.15],'EdgeColor',[0.82,0.75,0.75],'FaceAlpha',0.1);
%         patch([wu2,wu2(end:-1:1)],[CIu2(:,1);CIu2(end:-1:1,2)]',[0.15,0.15,0.15],'EdgeColor',[0.82,0.75,0.75],'FaceAlpha',0.1);  
%     end
%     
    %---------------------------------------------------------------------%
    hEb=patch([w,w(end:-1:1)],[CI(:,1);CI(end:-1:1,2)]',[0.43,0.35,0.35],'EdgeColor',[0.82,0.75,0.75],'FaceAlpha',0.1);
    
    LBF=squeeze(Lf(i,j,web))-squeeze(Lvinff(i,j,web));
    UBF=squeeze(Lf(i,j,web))-squeeze(Lvsupf(i,j,web));
   
    hEf = errorbar(w(web),squeeze(Lf(i,j,web)),LBF,UBF,'.','color',[0,0.45,0.74],'LineWidth',1.5);
    
    
    h12 = plot(w,squeeze(Lb(i,j,:)),'Color',[0.93,0.69,0.13],'LineWidth',3);
    hsb = plot(w,squeeze(L_sb(i,j,:)),'Color',[0.85,0.33,0.1],'LineWidth',3);
   
    
    htf=plot(w,squeeze(LTrab(i,j,:)),'--','Color','k','LineWidth',2);
    htb=plot(w(web),squeeze(LTraf(i,j,web)),'^','Color',[0.49,0.18,0.56],'LineWidth',2,'MarkerSize',7);
    xlim([0,w(end)+0.03]);
    box on;set(gcf,'color','white');
    
    if i==N;xlabel([' j = ' ,num2str(j)]);set(gca,'XTick',[0,0.5]);else;set(gca,'XTick',[]);end
        end
     if j==1; ylabel([' i = ' ,num2str(i)]);end
     set(gca,'FontSize',14)
    end
end
supAxes = [.08 .08 .84 .84];
[ax,h1] = suplabel('Normalized frequency','x',supAxes); set(h1, 'FontSize',14)
ttypt='idtf';
% [ax,h1] = suplabel(['{|' upper(ttype) '_i_j (\lambda)|}^{2}'],...
%             'y',supAxes);
% set(h1, 'FontSize',14)
% pos = get(h1,'Position');
% pos(1) = pos(1)+0.05; %0.0545
% set(h1,'Position',pos) % Adjust ylabel position


