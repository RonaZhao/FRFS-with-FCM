% y=[0.7,0.46,0.44;0.13,0.1,0.09;1.3,1.03,1.07;
%     0.51,0.35,0.33;0.36,0.34,0.26;
% 0.51,0.47,0.46;2.05,1.68,1.67
% ];
% % y=[21.07,19.32,19.22;
% % 37.19,56.41,55.9;
% % ];
% b=bar(y);
% grid on;
% ch = get(b,'children');
% set(gca,'XTickLabel',{'Glass','Iris','Cleveland','Wine','User','Vertebral','Sonar'})
% % set(gca,'XTickLabel',{'Yeast','Hillvaleey'})
% legend('NEW METHOD','T-FRFS','FRFS');
% xlabel('dataset ');
% ylabel('time(s)');
x=[100,200,300,400,500,600,700,800,900,1000,1100,1200,1500,2000,2500,3000,3500,4000,4500,5000];
 a=[0.08,0.13,0.24,0.36,0.53,0.73,0.99,1.31,1.67,2.14,2.95,3.57,7.18,16.58,31.63,54.17,82.85,122.32,171.75,234.14]; %a数据y值
 b=[0.04,0.09,0.2,0.35,0.53,0.7,0.95,1.27,1.65,2.13,2.8,3.62,7.11,15.2,29.65,53.37,81.8,121.2,173.76,233.54]; %b数据y值
 c=[0.05,0.1,0.2,0.34,0.5,0.71,0.97,1.26,1.64,2.14,2.81,3.53,6.88,15.04,30.21,51.97,81.82,121.67,172.69,233.27];
 plot(x,a,'-*b',x,b,'-or',x,c,'-xk'); %线性，颜色，标记
set(gca,'XTick',[100:200:5000]) %x轴范围1-6，间隔1
set(gca,'YTick',[0:20:250]) %y轴范围0-700，间隔100
legend('NEW METHOD','T-FRFS','FRFS');   %右上角标注
xlabel('objects')  %x轴坐标描述
ylabel('time(s)') %y轴坐标描述