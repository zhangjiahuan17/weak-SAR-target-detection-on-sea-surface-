function [STFP_t,STFP_c,NR_t,L_tt,L_cc,MS_t_m,NR_c,MS_c_m]=NTFD_feature(m, n, L_pt, L_pc, L_short, Nor_spwv, Nor_spwv_c, K, con_t, con_c, tarcell, clucell)
% �����˾�����ҪʱƵ�� (Significant Time-Frequency Point, STFP)�Ķ�ֵͼ��
% input:
%       m,n = size(tfr_spwv)
%       L_pt, L_pc: �ֱ�Ϊֻ���Ӳ��Ĳο���Ԫ���ȺͰ���Ŀ��Ļز���Ԫ����
%       L_short: ȡ��һ���ֵĻز��ĳ���
%       Nor_spwv,Nor_spwv_c��Ŀ����Ӱ�쵥Ԫ���Ӳ���Ԫ�Ĺ�һ��SPWVD
%       K���ڹ�һ��SPWVD��ʱ��߶��ϣ�ȡǰK���������ֵ���б��
%       con_t,con_c:4��ͨor8��ͨ
%       tarcell,clucell: ��ȡ�������õ�Ŀ�굥Ԫ���Ӳ���Ԫ
% output:
%       STFP_t,STFP_c����ҪʱƵ��
%       L_tt, L_cc�����ֵ�������꣨Ƶ����ģ�
%       NR_t,NR_c����ͨ������Ŀ��t��Ŀ��ز��ģ�c���Ӳ���
%       MS_t_m,MS_c_m�������ͨ�������
%% normalization TF ridge
% �����˾�����ҪʱƵ�� (Significant Time-Frequency Point, STFP)�Ķ�ֵͼ��
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

%% TF��������NR MS RI?
% NR MS ��STFPͼ���л�õ�������N
%Ŀ��
[L_t,NR_t] =  bwlabel(STFP_t(:,:,tarcell),con_t);%����4������ͨ��Ҳ��Ϊ8��
img_reg_t = regionprops(L_t);%����ͼƬ����
MS_t = [img_reg_t.Area];%��ȡ��ͨ�������Ϣ
MS_t_m = max(MS_t);
%�Ӳ�
[L_c,NR_c] =  bwlabel(STFP_c(:,:,clucell),con_c);%����4������ͨ��Ҳ��Ϊ8��
img_reg_c = regionprops(L_c);%����ͼƬ����
MS_c = [img_reg_c.Area];%��ȡ��ͨ�������Ϣ
MS_c_m = max(MS_c);


