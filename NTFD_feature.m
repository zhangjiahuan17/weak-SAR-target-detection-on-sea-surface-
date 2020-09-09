function [STFP_t,STFP_c,NR_t,L_tt,L_cc,MS_t_m,NR_c,MS_c_m]=NTFD_feature(m, n, L_pt, L_pc, L_short, Nor_spwv, Nor_spwv_c, K, con_t, con_c, tarcell, clucell)
% 生成了具有重要时频点 (Significant Time-Frequency Point, STFP)的二值图像
% input:
%       m,n = size(tfr_spwv)
%       L_pt, L_pc: 分别为只有杂波的参考单元长度和包含目标的回波单元长度
%       L_short: 取了一部分的回波的长度
%       Nor_spwv,Nor_spwv_c：目标所影响单元和杂波单元的归一化SPWVD
%       K：在归一化SPWVD的时间尺度上，取前K个最大像素值进行标记
%       con_t,con_c:4连通or8连通
%       tarcell,clucell: 提取特征所用的目标单元和杂波单元
% output:
%       STFP_t,STFP_c：重要时频点
%       L_tt, L_cc：最大值出现坐标（频率轴的）
%       NR_t,NR_c：连通区域数目，t是目标回波的，c是杂波的
%       MS_t_m,MS_c_m：最大连通区域面积
%% normalization TF ridge
% 生成了具有重要时频点 (Significant Time-Frequency Point, STFP)的二值图像
STFP_c = zeros(m,n,L_pt);
STFP_t = zeros(m,n,L_pc);
% ridge_t=zeros(L_short, L_pt);
% ridge_c=zeros(L_short, L_pc);

L_tt = zeros(L_short,L_pt);
for ii = 1:L_pt
    for jj = 1:L_short
%         [ridge_t(jj,ii), L_tt(jj,ii)] = max(Nor_spwv(:,jj,ii));
        [~, L_tt(jj,ii)] = max(Nor_spwv(:,jj,ii));
        [~,I] = maxk(Nor_spwv(:,jj,ii),K);
        %[B,I] = maxk(___) finds the indices of the largest k values of A and returns them in I.
        STFP_t(I,jj,ii)=1;
    end
end
L_cc = zeros(L_short,L_pc);
for ii = 1:L_pc
    for jj = 1:L_short
        [~, L_cc(jj,ii)] = max(Nor_spwv_c(:,jj,ii));
        [~,I] = maxk(Nor_spwv_c(:,jj,ii),K);
        %[B,I] = maxk(___) finds the indices of the largest k values of A and returns them in I.
        STFP_c(I,jj,ii)=1;
    end
end

%% TF三特征：NR MS RI?
% NR MS 从STFP图像中获得的特征，N
%目标
[L_t,NR_t] =  bwlabel(STFP_t(:,:,tarcell),con_t);%基于4邻域连通（也可为8）
img_reg_t = regionprops(L_t);%计算图片属性
MS_t = [img_reg_t.Area];%获取连通域面积信息
MS_t_m = max(MS_t);
%杂波
[L_c,NR_c] =  bwlabel(STFP_c(:,:,clucell),con_c);%基于4邻域连通（也可为8）
img_reg_c = regionprops(L_c);%计算图片属性
MS_c = [img_reg_c.Area];%获取连通域面积信息
MS_c_m = max(MS_c);


