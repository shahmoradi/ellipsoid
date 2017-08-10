clear;
clc;
close all;

npi = 3;
npf = 20;
nd = 2;
m = 1;

np_scale_data = zeros(2,npf-npi+1);
density_scale_data = zeros(2,npf-npi+1);

for i = npi:npf
    dataScale = load(['_allUniform_ScaleData_',num2str(nd),'_',num2str(i),'.txt']);
    np_scale_data(1,m) = i;
    np_scale_data(2,m) = median(log(dataScale(:,6)));
    dataScale(:,4) = log(i ./ dataScale(:,4));
    figure
    scatter(dataScale(:,4), log(dataScale(:,6)));
    xlabel('log(density)');
    ylabel('log(scale factor)');
    title(['nd: ',num2str(nd), ', np: ', num2str(i)]);
    
    m = m+1;
end

figure 
scatter(np_scale_data(1,:),np_scale_data(2,:));

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

