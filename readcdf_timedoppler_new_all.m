clc;
close all;
clear;
%% Read the necessary data %%

file = 'D:\桌面\目标检测\数据\19980204_202525_ANTSTEP.CDF';
%file = 'D:\桌面\目标检测\数据\IPIX\dataset with weak target\#17\19931107_135603_starea.cdf';
ncdisp(file);
finfo = ncinfo(file);
azi = ncread(file,'azimuth_angle');
range = ncread(file,'range');
data = int16(ncread(file,'adc_data'));
N = length(azi); 
n_range = length(range);
data = permute(data,[4 3 2 1]); 
pol = 'hv'; % hh, hv, vh, vv
mode = 'auto' ;% 'raw' (no pre-processing) 
           % 'auto' (automatic pre-processing) 
           % 'dartmouth' (pre-processing for dartmouth files containing land)
is_93_data = false;% 注意这里要根据数据年份改一下！
%% If try to use the data of 93, use this part to correct the data %%          
if is_93_data == true
    data = is_data_93(data);
end
%% 
for rangebin = 1:n_range
    [I(:,rangebin),Q(:,rangebin),meanIQ,stdIQ,inbal] = NewIpixLoad(finfo,data,pol,rangebin,mode);
    R(:,rangebin) = I(:,rangebin) + 1i * Q(:,rangebin);
%    R_abs(:,rangebin) = abs(I(:,rangebin) + 1i * Q(:,rangebin));
%    r(:,rangebin) = abs(R(:,rangebin)) / max(abs(R(:,rangebin)));
end
clear I Q 
%%
t = 0.001 : 0.001 : (N/1000);
figure
imagesc(range,t,abs(R));
colormap('jet');
colorbar;
title('杂波总览俯视-hv');
xlabel('距离（m）');
ylabel('时间（s）');
%% Data handling and image forming %%   
%计算Wigner-Ville分布
%[tfr,t,f]=tfrwv(sig_9_short,1:length(sig_9_short),256);

dis =10; % the range bin you wannt to observe
sig_10 = R(:,dis);
L = length(sig_10);

[tfr,T,f]=tfrstft(sig_10,1:L,1024);
fs=1024;
tfr = fftshift(tfr,1);
% figure,imagesc(T/fs,f*fs,abs(tfr));
figure,imagesc(T/fs,fftshift(f,1)*fs,abs(tfr));
ylim([-500 500]);
xlabel('时间(s)','Fontsize',20);
ylabel('多普勒频率(Hz)','Fontsize',20);
title('IPIX 19980204-202525 RangeCell 10','Fontsize',20);
axis xy;
Path='D:\桌面\目标检测\数据\IPIX！\19980204-202525\referencecell.tif';
print(Path,'-dpng','-r600')
%%
dis = 7; % the range bin you wannt to observe
sig_7 = R(:,dis);

L = length(sig_7);
%计算Wigner-Ville分布
%[tfr,t,f]=tfrwv(sig_9_short,1:length(sig_9_short),256);

[tfr,T,f]=tfrstft(sig_7,1:L,1024);
tfr = fftshift(tfr,1);
fs=1024;
figure,imagesc(T/fs,fftshift(f,1)*fs,abs(tfr));
ylim([-500 500]);
xlabel('时间(s)','Fontsize',20);
ylabel('多普勒频率(Hz)','Fontsize',20);
title('IPIX 19931111-163625 RangeCell 7','Fontsize',20);
axis xy;
Path='D:\桌面\目标检测\数据\IPIX！\19980204-202525\targetcell.tif';
print(Path,'-dpng','-r600')
clear sig_2 sig_5 sig_9 tfr T f
%% SPWVD
L_short = L/4;
[tfr_spwv,T,f]=tfrspwv(sig_7, 1:L_short, 128, kaiser(31), kaiser(63));
fs=1024;
tfr_spwv=fftshift(tfr_spwv,1);
h1=figure;
imagesc(T/fs,(f-0.25)*fs,abs(tfr_spwv));
% set(h1,'position',[100,100,260,200]);
ylim([-200 200]);
set(gca,'fontsize',8);
xlabel('时间(s)');
ylabel('多普勒频率(Hz)');
title('IPIX 19980204-202525 RangeCell 7');
axis xy;
Path='D:\桌面\目标检测\数据\IPIX！\19980204-202525\targetcellspwvd.tif';
print(Path,'-dpng','-r600')

%% SPWVD Normalization
P_c = [1:5,9:28];%REFERENCE CELL
[m,n]=size(tfr_spwv);
temp_c = zeros(m,n,length(P_c));
tfr_spwv_ci = zeros(m,n);
j = 1;
for i = P_c
    [tfr_spwv_ci,T,f]=tfrspwv(R(:,i), 1:L_short, 128, kaiser(31), kaiser(63));
    temp_c(:,:,j) = fftshift(tfr_spwv_ci);
    j=j+1;
end
mu_sum = sum(temp_c(:,:,:),3);
mu = mu_sum/length(P_c); % 海杂波均值函数
sigma_temp=(temp_c-mu).^2;
sigma=sqrt(1/(length(P_c)-1)*sum(sigma_temp(:,:,:),3));% 海杂波标准差函数

% figure,imagesc(T/fs,(f-0.25)*fs,abs(fftshift(mu,1)));
figure,imagesc(T/fs,(f-0.25)*fs,abs(mu));
colormap(jet);
colorbar;
xlabel('时间(s)');
ylabel('多普勒频率(Hz)');
title('IPIX 19980204-202525 参考单元SPWVD均值函数');
axis xy;
Path='D:\桌面\目标检测\数据\IPIX！\19980204-202525\RC-mean.tif';
print(Path,'-dpng','-r600')

% 杂波的归一化SPWVD
Nor_spwv_c = (temp_c - mu)./sigma;
RR_clutter = sum(Nor_spwv_c,3)/length(P_c);

% figure,imagesc(T/fs,(f-0.25)*fs,abs(fftshift(RR_clutter)));
figure,imagesc(T/fs,(f-0.25)*fs,abs(RR_clutter));
colormap;
colorbar;
xlabel('时间(s)');
ylabel('多普勒频率(Hz)');
title('IPIX 19980204-202525 海杂波归一化');
axis xy;
Path='D:\桌面\目标检测\数据\IPIX！\19980204-202525\海杂波回波归一化SPWVD-sum.tif';
print(Path,'-dpng','-r600')

figure,imagesc(T/fs,(f-0.25)*fs,abs(fftshift(sigma)));
colormap(jet);
colorbar;
xlabel('时间(s)');
ylabel('多普勒频率(Hz)');
title('IPIX 19980204-202525 参考单元SPWVD标准差函数');
axis xy;
Path='D:\桌面\目标检测\数据\IPIX！\19980204-202525\RC-standard deviation.tif';
print(Path,'-dpng','-r600')

% 包含目标的归一化SPWVD
P_t = 6:8; % all target CELL
[m,n]=size(tfr_spwv);
tfr_spwv_ti = zeros(m,n);
temp_t = zeros(m,n,length(P_t));
jj=1;
for i = P_t
    [tfr_spwv_ti,T,f]=tfrspwv(R(:,i), 1:L_short, 128, kaiser(31), kaiser(63));
    temp_t(:,:,jj) = fftshift(tfr_spwv_ti);
    jj=jj+1;
end

Nor_spwv = (temp_t - mu)./sigma;
RR_target = sum(Nor_spwv,3)/length(P_t);
figure,imagesc(T/fs,(f-0.25)*fs,abs(RR_target));
% colormap(jet);
colorbar;
xlabel('时间(s)');
ylabel('多普勒频率(Hz)');
title('IPIX 19980204-202525 目标归一化');
axis xy;
Path='D:\桌面\目标检测\数据\IPIX！\19980204-202525\目标回波归一化SPWVD.tif';
print(Path,'-dpng','-r600')

clear temp_c temp_t sigma_temp  RR_target RR_clutter tfr_spwv_i mu_sum mu_sum_tar mu_tar
%% normalization TF ridge
% 生成了具有重要时频点 (Significant Time-Frequency Point, STFP)的二值图像
L_pt=length(P_t);
L_pc=length(P_c);
K = 5;
con_t=8;
con_c=4;
tarcell=2; 
clucell = 1;
[STFP_t,STFP_c,NR_t,L_tt, L_cc,MS_t_m,NR_c,MS_c_m] = ...
NTFD_feature(m, n, L_pt, L_pc, L_short, Nor_spwv, Nor_spwv_c, ...
K, con_t, con_c, tarcell, clucell);

figure,imagesc(T/fs,(f-0.25)*fs,abs(STFP_t(:,:,3)));
colormap;
colorbar;
xlabel('时间(s)');
ylabel('多普勒频率(Hz)');
title('IPIX 19980204-202525 目标回波STFP图');
axis xy;
Path='D:\桌面\目标检测\数据\IPIX！\19980204-202525\STFP-target-L5.tif';
print(Path,'-dpng','-r600')

figure,imagesc(T/fs,(f-0.25)*fs,abs(STFP_c(:,:,1)));
colormap;
colorbar;
xlabel('时间(s)');
ylabel('多普勒频率(Hz)');
title('IPIX 19980204-202525 海杂波STFP图');
axis xy;
Path='D:\桌面\目标检测\数据\IPIX！\19980204-202525\STFP-clutter-l5.tif';
print(Path,'-dpng','-r600')

%% TF三特征：NR MS RI?
% NR MS 从STFP图像中获得的特征，N
%目标
disp(['目标单元NR：',num2str(NR_t),',最大面积MS：',num2str(MS_t_m)]);
%杂波
disp(['杂波单元NR：',num2str(NR_c),',最大面积MS：',num2str(MS_c_m)]);
% 画两个彩色的图
RGB_t = label2rgb(L_t);
figure,imagesc(T/fs,(f-0.25)*fs,abs(fftshift(RGB_t,1)));
axis xy
xlabel('时间(s)');
ylabel('多普勒频率(Hz)');
Path='D:\桌面\目标检测\数据\IPIX！\19980204-202525\STFP-target-rgb-L5.tif';
print(Path,'-dpng','-r600')

RGB_c = label2rgb(L_c);
figure,imagesc(T/fs,(f-0.25)*fs,abs(fftshift(RGB_c,1)));
axis xy
xlabel('时间(s)');
ylabel('多普勒频率(Hz)');
Path='D:\桌面\目标检测\数据\IPIX！\19980204-202525\STFP-clutter-rgb-L5.tif';
print(Path,'-dpng','-r600')

%% RI
ri = -250 : 500/(1024-1) : 250;
RI_t = zeros(size(L_tt));
for i = 1 : length(P_t)
    for j = 1 : L_short
        RI_t(j,i) = ri(L_tt(j,i));
    end
end
RI_c = zeros(size(L_cc));
for i = 1 : length(P_c)
    for j = 1 : L_short
        RI_c(j,i) = ri(L_cc(j,i));
    end
end    

RI_t = sum(RI_t(:,2:3),2)/2;

RI_c = sum(RI_c,2)/length(P_c);
% 绘制"直方图" 
rd_t_mean = mean(RI_t);
rd_t_std = std(RI_t);

section = floor(((RI_t - rd_t_mean) /rd_t_std) *255);
RI_t_x=zeros(1,256);
for k=0:255 
    RI_t_x(k+1)=length(find(section==k));
end 

rd_c_mean = mean(RI_c);
rd_c_std = std(RI_c);
section = floor(((RI_c - rd_c_mean) /rd_c_std) *255);
RI_c_x=zeros(1,256);
for k=0:255 
    RI_c_x(k+1)=length(find(section==k));
end

figure,bar(RI_t_x,'r');
set(gca,'XScale','log')
set(gca,'YScale','log')
legend()
hold on;
bar(RI_c_x,'b');
alpha(0.7);
grid on
legend('target','clutter')
Path='D:\桌面\目标检测\数据\IPIX！\19980204-202525\RI-log.tif';
print(Path,'-dpng','-r600')

%% 暂时不管了

% % NR
% nr_c_max = max(NR_c);
% nr_c_min = min(NR_c);
% nr_c_mean = mean(NR_c);
% section = floor(((NR_c-nr_c_min)/(nr_c_max-nr_c_min))*255);
% nr_c_bar=zeros(1,256);
% for k=0:255 
%     nr_c_bar(k+1)=length(find(section==k));
% end
% nr_t_max = max(NR_t);
% nr_t_min = min(NR_t);
% nr_t_mean = mean(NR_t);
% section = floor(((NR_t-nr_t_min)/(nr_t_max-nr_t_min))*255);
% nr_t_bar=zeros(1,256);
% for k=0:255
%    nr_t_bar(k+1)=length(find(section==k));
% end
% 
% figure,bar(nr_c_bar,'b');
% set(gca,'YScale','log')
% hold on;
% bar(nr_t_bar,'r');
% % MS
% ms_c_max = max(MS_c);
% ms_c_min = min(MS_c);
% ms_c_mean = mean(MS_c);
% section = floor(((MS_c-ms_c_min)/(ms_c_max-ms_c_min))*255);
% ms_c_bar=zeros(1,256);
% for k=0:255 
%     ms_c_bar(k+1)=length(find(section==k));
% end
% ms_t_max = max(MS_t);
% ms_t_min = min(MS_t);
% ms_t_mean = mean(MS_t);
% section = floor(((MS_t-ms_t_min)/(ms_t_max-ms_t_min))*255);
% ms_t_bar=zeros(1,256);
% for k=0:255
%     ms_t_bar(k+1)=length(find(section==k));
% end
% 
% figure,bar(ms_c_bar,'b');
% set(gca,'YScale','log')
% hold on;
% bar(ms_t_bar,'r');