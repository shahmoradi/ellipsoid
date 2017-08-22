clear;
clc;
close all;

listing = dir('../out/*.txt');
fvals = split(listing(1).name,'_');
ndMin = str2double(fvals(3));
npMin = str2double(fvals(4));
ndMax = str2double(fvals(3));
npMax = str2double(fvals(4));

for k=2:length(listing)
   fvals = split(listing(k).name,'_');
   
   if str2double(fvals(3)) > ndMax
       ndMax = str2double(fvals(3));
   elseif str2double(fvals(3)) < ndMin
       ndMin = str2double(fvals(3));
   end
       
   if str2double(fvals(4)) > npMax
       npMax = str2double(fvals(4));
   elseif str2double(fvals(4)) < npMin
       npMin = str2double(fvals(4));
   end
end

densitymat = [];

for k = ndMin:ndMax
    len = 14;
    np_scale_data = zeros(2,len);
    neg = zeros(1,len);
    pos = zeros(1,len);
    m = 1;
    
    %legendvals = {};
    %figure
    for i = k+1:k:k*15
        
        dataScale = load(fullfile('../out',['allUniform_ScaleData_',num2str(k),'_',num2str(i),'_.txt']));
        
        np_scale_data(1,m) = log10(i);
        np_scale_data(2,m) = median(log10(log10((dataScale(:,6)))));
        neg(1,m) = abs(np_scale_data(2,m) - quantile(log10(log10(dataScale(:,6))),.25));
        pos(1,m) = abs(np_scale_data(2,m) - quantile(log10(log10(dataScale(:,6))),.75));
        
        dataScale(:,4) = log10(i) - dataScale(:,4);
        dataScale(:,6) = log10(dataScale(:,6));
        
        densitymat = sortrows(dataScale(:,4:2:6));
        densitymat(:,2) = movmean(densitymat(:,2),50);
    
        figure
        scatter(dataScale(:,4), dataScale(:,6),'.','blue');
        xlabel('log_{10}(density)');
        ylabel('log_{10}(scale factor)');
        title(['nd = ',num2str(k),', np = ', num2str(i)]);
        hold on; 
        scatter(densitymat(:,1),densitymat(:,2),'.','red');
        hold off;
        saveas(gcf, fullfile('../out',['density_scale_',num2str(k),'_',num2str(i),'.png']));
        
        %a = [num2str(i)];
        %legendvals = [legendvals a]
        
        m = m + 1;
        
    end
%     densitymat = sortrows(densitymat);
%     legend(legendvals);
%     xlabel('log_{10}(density)');
%     ylabel('log_{10}(scale factor)');
%     title(['nd: ',num2str(k)]);
%     hold off;
%     
%     windowSize = 100; 
%     b = (1/windowSize)*ones(1,windowSize);
%     a = 1;
%     densitymat(:,2) = filter(b,a,densitymat(:,2));
%     
%     figure
%     plot(densitymat(:,1),densitymat(:,2),'.');
%     xlabel('log_{10}(density)');
%     ylabel('log_{10}(scale factor)');
%     title(['nd = ',num2str(k)]);
%     saveas(gcf, fullfile('../out',['density_scale_',num2str(k),'.png']));
    
    figure 
    errorbar(np_scale_data(1,:),np_scale_data(2,:),neg,pos);
    hold on;
    scatter(np_scale_data(1,:),np_scale_data(2,:));
    hold off;
    xlabel('number of points');
    ylabel('log_{10}(scale factor)');
    title(['nd = ',num2str(k)]);
    saveas(gcf, fullfile('../out',['np_scale_',num2str(k),'.png']));
    
end

% for i = 1:10
%     n = num2str(i);
%     m = num2str(np);
%     data1 = load(['_',n,'_allUniform_BoundTrue_1_',nd,'_',np,'.txt']);
%     data2 = load(['_',n,'_allUniform_CenterTrue_1_',nd,'_',np,'.txt']);
%     data3 = load(['_',n,'_allUniform_MembershipTrue_1_',nd,'_',np,'.txt']);
%     data4 = load(['_',n,'_allUniform_AMVETrue_1_',nd,'_',np,'.txt']);
%     
%     figure
%     scatter(data1(:,1),data1(:,2),'.','blue');
%     hold on;
%     scatter(data2(:,1),data2(:,2),'o','red');
%     hold on;
%     scatter(data3(:,1),data3(:,2),200,'.');
%     hold on;
%     scatter(data4(:,1),data4(:,2),'.','green');
%     hold on;
%     scatter(dataScale(i,1),dataScale(i,2),dataScale(i,3),400,'.','red');
%     hold off;
%     legend('Original Ellip Bound','Original Ellip Center','Original Ellip Points','AMVE');
% end

