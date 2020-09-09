function [I,Q,meanIQ,stdIQ,inbal] = NewIpixLoad(finfo,data,pol,rangebin,mode)

H_txpol     = 1; 
V_txpol     = 2;
Like_adc_I  = 1;  % when use data of 93 set these parameters(L_I,L_Q,C_I,C_Q) to 1 2 3 4;
Like_adc_Q  = 2;  % when use data of 98 set these parameters(L_I,L_Q,C_I,C_Q) to 3 4 1 2;
Cross_adc_I = 3;
Cross_adc_Q = 4;

%% extract correct polarization from cdffile %%


pol = lower(pol);
if length(size(data)) == 3,
    % read global attribute TX_polarization
    txpol = finfo.Attributes(6).Value;
    if pol(1) ~= lower(txpol),
        fname = strrep(finfo.Attributes(7).Value,'\0','');
        disp(['Warning: file ' fname ' does not contain ''' ...
            txpol(1) ''' transmit polarization.']);
    end
    switch pol
        case {'hh','vv'}, xiq = data(:,rangebin,[Like_adc_I Like_adc_Q]);
        case {'hv','vh'}, xiq = data(:,rangebin,[Cross_adc_I Cross_adc_Q]);
    end
    I = permute(xiq(:,:,1),[1 2 3]);
    Q = permute(xiq(:,:,2),[1 2 3]);
else
    switch pol
        case 'hh', xiq = data(:,H_txpol,rangebin,[Like_adc_I Like_adc_Q]);
        case 'hv', xiq = data(:,H_txpol,rangebin,[Cross_adc_I Cross_adc_Q]);
        case 'vv', xiq = data(:,V_txpol,rangebin,[Like_adc_I Like_adc_Q]);
        case 'vh', xiq = data(:,V_txpol,rangebin,[Cross_adc_I Cross_adc_Q]);
    end
    I = permute(xiq(:,1,:,1),[1 3 4 2]);
    Q = permute(xiq(:,1,:,2),[1 3 4 2]);
end

 I = double(I);
 Q = double(Q);

%% apply corrections to I and Q data %%

switch mode,
    case 'raw',
        meanIQ = [0 0];
        stdIQ = [1 1];
        inbal = 0;
    case 'auto',
        % Pre-processing     
        meanIQ = mean([I(:) Q(:)]);
        stdIQ = std([I(:) Q(:)],1);
        I = (I - meanIQ(1)) / stdIQ(1);
        Q = (Q - meanIQ(2)) / stdIQ(2);
        sin_inbal = mean(I(:) .* Q(:));
        inbal = asin(sin_inbal) * 180 / pi;
        I = (I - Q * sin_inbal) / sqrt(1 - sin_inbal^2);
    case 'dartmouth',
        % Define rectangular patches of land in Dartmouth campaign.
        % Format: [azmStart azmEnd  rangeStart rangeEnd]
        landcoord=[
            0  70     0  600;
            305 360     0  600;
            30  55     0 8000;
            210 305     0 4700;
            320 325  2200 2700;
            ];
        % Exclude land from data used to estimate pre-processing parameters
        azm = mod(ipixazm(azi),360);
        range = rangebin;
        nbin = length(rangebin);
        ok = ones(size(I));
        for i=1:size(landcoord,1)
            for r=1:nbin
                if range(r) >= landcoord(i,3) && range(r) <= landcoord(i,4)
                    ok(find(azm >= landcoord(i,1) & azm <= landcoord(i,2)),r)=0;
                end
            end
        end
        ok = find(ok);
        if length(ok) < 100
            disp('Warning: not enough sweeps for land-free pre-processing.');
            ok = ones(size(I));
        end
        % Pre-processing
        meanIQ = mean([I(ok) Q(ok)]);
        stdIQ = std([I(ok) Q(ok)],1);
        I = (I - meanIQ(1)) / stdIQ(1);
        Q = (Q - meanIQ(2)) / stdIQ(2);
        sin_inbal = mean(I(ok) .* Q(ok));
        inbal = asin(sin_inbal) * 180 / pi;
        I = (I - Q * sin_inbal) / sqrt(1 - sin_inbal^2);
end
