% Copyright (C) 2018,2023 Mitsubishi Electric Research Laboratories (MERL)
%
% SPDX-License-Identifier: AGPL-3.0-or-later

% angular-domain mmWave channel estimation with 1-bit ADC: Cramer-Rao bound
%
% SNR is defined on a per receive antenna per pilot basis
%
% Refer to
%
% P. Wang, J. Li, M. Pajovic, P.T. Boufounos, P.V. Orlik,
% "On Angular-Domain Chanel Estimation for One-Bit Massive MIMO Systems with Fixed and Time-Varying Thresholds",
% 2017 Asilomar Conference on Signals, Systems and Computers (ACSSC), Pacific Grove, CA, November, 2017.
%

clc;clear;close all;

%% initialize parameters

nMC  = 100;     % # of Monte-Carlo runs; in the paper we used 100 runs;

N_tx = 8;       % # of Tx antennas
N_rx = 16;      % # of Rx antennas

SNRset = [-10:2.5:30];  % [dB]
snrSet = 10.^(SNRset/10);

M    = 4;    % # of channel paths
K    = 100;  % # of pilots; in the paper we used 100 pilots;
nTQ  = 50;   % # of time-varying thresholds averaged for each Monte-Carlo run; in the paper we used 50;

%% compute angular-domain one-bit CRB with fixed zero-threshold and time-varying thresholds

% channel matrix: load pre-saved channel path gain matrix
load('hTemp.mat');
hMatrix   = diag(hTemp);
hVec      = hMatrix(:);

for lmc = 1:nMC

    % AoA and AoD are uniformly selected from four clusters (see)
    % clusters: [-50,-45], [-10,-5],[30,36], [65,75]
    thetaVec = [-50+5*rand(1,1), -10+5*rand(1,1),30+6*rand(1,1), 65+10*rand(1,1)];  % Rx (AoA)
    % clusters: [15,18], [45,50], [-20,-16], [-72,-70],
    phiVec   = [15+3*rand(1,1), 45+5*rand(1,1), -20+4*rand(1,1), -72+2*rand(1,1)  ];  % Tx (AoD)

    % steering matrices
    A_ms = zeros(N_tx, M);
    A_bs = zeros(N_rx, M);

    indTx = (0:N_tx-1);
    indRx = (0:N_rx-1);
    A_ms  = 1/sqrt(N_tx) * exp(1i*pi*indTx(:)*sind(phiVec));
    A_bs  = 1/sqrt(N_rx) * exp(1i*pi*indRx(:)*sind(thetaVec));

    for lsnr = 1:length(SNRset) % iterate over Monte-Carlo runs

        snr = snrSet(lsnr);
        S   = exp(1i*pi*(randi(4,N_tx,K)-1)/2)/sqrt(N_tx);  % normalized pilots signals; [Nt x K]
        Gamma = kron(S.'*conj(A_ms),A_bs);

        %% CRB
        GammaR = real(Gamma);
        GammaI = imag(Gamma);
        hR     = real(hVec);
        hI     = imag(hVec);

        % I/Q channels
        yR     = GammaR*hR-GammaI*hI;
        yI     = GammaR*hI+GammaI*hR;

        % compute noise variance according to SNR
        nvar = sum(abs(yR).^2+abs(yI).^2)/snr/N_rx/K;

        [lmc,10*log10(sum(abs(yR).^2+abs(yI).^2)/(nvar*N_rx*K))]

        % optimal threshold: known noise variance
        thresR0  = GammaR*hR-GammaI*hI;
        thresI0  = GammaR*hI+GammaI*hR;

        tempFIM1 = zeros(4*M,4*M);
        tempCRB1 = zeros(4*M,4*M);
        tempFIM1 = func_crb_knownNVar(thresR0(:), thresI0(:), thetaVec, phiVec, Gamma, hVec, S, A_ms, A_bs, N_tx, N_rx, M, nvar);
        tempCRB1 = inv(tempFIM1);
        tempcrb1 = diag(tempCRB1);

        crb_Theta1(:,lsnr,lmc)     = tempcrb1(1:M);
        crb_Phi1(:,lsnr,lmc)       = tempcrb1(M+1:2*M);
        crb_AlphaR1(:,lsnr,lmc)    = tempcrb1(2*M+1:3*M);
        crb_AlphaI1(:,lsnr,lmc)    = tempcrb1(3*M+1:4*M);

        % fixed-threshold (zeros): known noise variance
        thresR1  = zeros(N_rx*K,1);
        thresI1  = zeros(N_rx*K,1);

        tempFIM3 = zeros(4*M,4*M);
        tempCRB3 = zeros(4*M,4*M);
        tempFIM3 = func_crb_knownNVar(thresR1(:), thresI1(:), thetaVec, phiVec, Gamma, hVec, S, A_ms, A_bs, N_tx, N_rx, M, nvar);
        tempCRB3 = inv(tempFIM3);
        tempcrb3 = diag(tempCRB3);

        crb_Theta3(:,lsnr,lmc)     = tempcrb3(1:M);
        crb_Phi3(:,lsnr,lmc)       = tempcrb3(M+1:2*M);
        crb_AlphaR3(:,lsnr,lmc)    = tempcrb3(2*M+1:3*M);
        crb_AlphaI3(:,lsnr,lmc)    = tempcrb3(3*M+1:4*M);

        % time-varying thresholds
        nThres    = 8;                % # of levels
        hMaxR     = max(abs(yR));
        thresSetR = linspace(-hMaxR, hMaxR, nThres);
        hMaxI     = max(abs(yI));
        thresSetI = linspace(-hMaxI, hMaxI, nThres);

        for ll = 1:nTQ % iterate over different realization of time-varying thresholds
            indR2    = randi(nThres,N_rx*K,1);   % uniformly select one level as the threshold
            thresR2  = thresSetR(indR2);
            indI2    = randi(nThres,N_rx*K,1);   % uniformly select one level as the threshold
            thresI2  = thresSetI(indI2);

            % unknown noise variance
            tempFIM5 = zeros(4*M+1,4*M+1);
            tempCRB5 = zeros(4*M+1,4*M+1);
            tempFIM5 = func_crb_unknownNVar(thresR2(:), thresI2(:), thetaVec, phiVec, Gamma, hVec, S, A_ms, A_bs, N_tx, N_rx, M, nvar);
            tempCRB5 = inv(tempFIM5);
            tempcrb5 = diag(tempCRB5);

            crb_Theta5(:,lsnr,ll,lmc)     = tempcrb5(1:M);
            crb_Phi5(:,lsnr,ll,lmc)       = tempcrb5(M+1:2*M);
            crb_AlphaR5(:,lsnr,ll,lmc)    = tempcrb5(2*M+1:3*M);
            crb_AlphaI5(:,lsnr,ll,lmc)    = tempcrb5(3*M+1:4*M);
            crb_nVar5(:,lsnr,ll,lmc)      = tempcrb5(4*M+1:4*M+1);

            % known noise variance
            tempFIM6 = zeros(4*M,4*M);
            tempCRB6 = zeros(4*M,4*M);
            tempFIM6 = func_crb_knownNVar(thresR2(:), thresI2(:), thetaVec, phiVec, Gamma, hVec, S, A_ms, A_bs, N_tx, N_rx, M, nvar);

            tempCRB6 = inv(tempFIM6);
            tempcrb6 = diag(tempCRB6);

            crb_Theta6(:,lsnr,ll,lmc)     = tempcrb6(1:M);
            crb_Phi6(:,lsnr,ll,lmc)       = tempcrb6(M+1:2*M);
            crb_AlphaR6(:,lsnr,ll,lmc)    = tempcrb6(2*M+1:3*M);
            crb_AlphaI6(:,lsnr,ll,lmc)    = tempcrb6(3*M+1:4*M);
        end

    end

end

% --------------------- Remarks --------------------------------
%
% a) you can also use 'func_crb_unknownNVar.m' to compute FIMs for
%     1) optimal threshold with unknown noise variance
%     2) zero threshold with unknown noise variance
% and you will find out corresponding FIMs are singular (as we noted in the
% paper) due to the amplitude/noise variance ambiguity.
%
% b) CRBs for time-varying thresholds are non-singular for both known and
% unknown noise variance.

%% plot Fig. 3 of the paper

% AoA (Fig. 3 (a))
figure(100);
subplot(2,2,1)
for lmc =1:nMC
    plotb1 = semilogy(SNRset,crb_Theta1(1,:,lmc),'k-+');  hold on; % OQ
    plotb2 = semilogy(SNRset,crb_Theta3(1,:,lmc),'m--x'); hold on; % FQ
    temp5  = squeeze(mean(crb_Theta5(1,:,:,lmc),3));
    temp6  = squeeze(mean(crb_Theta6(1,:,:,lmc),3));
    plotc1 = semilogy(SNRset,temp5,'b-o'); hold on; % TQ
    plotc2 = semilogy(SNRset,temp6,'r-s'); hold on; % TQ
end
xlim([SNRset(1),SNRset(end)]); grid on;
xlabel('SNR (dB)');
title(sprintf('K=%g, AoA(1)', K));
legend([plotb1(1), plotb2(1), plotc1(1),plotc2(1)],'CQ (known \sigma^2)','FQ (known \sigma^2)','TQ (unknown \sigma^2)','TQ (known \sigma^2)')

subplot(2,2,2)
for lmc =1:nMC
    plotb1 = semilogy(SNRset,crb_Theta1(2,:,lmc),'k-+');  hold on; % OQ
    plotb2 = semilogy(SNRset,crb_Theta3(2,:,lmc),'m--x'); hold on; % FQ
    temp5  = squeeze(mean(crb_Theta5(2,:,:,lmc),3));
    temp6  = squeeze(mean(crb_Theta6(2,:,:,lmc),3));
    plotc1 = semilogy(SNRset,temp5,'b-o'); hold on; % TQ
    plotc2 = semilogy(SNRset,temp6,'r-s'); hold on; % TQ
end
xlim([SNRset(1),SNRset(end)]); grid on;
xlabel('SNR (dB)');
title(sprintf('AoA(2)'));

subplot(2,2,3)
for lmc =1:nMC
    plotb1 = semilogy(SNRset,crb_Theta1(3,:,lmc),'k-+');  hold on; % OQ
    plotb2 = semilogy(SNRset,crb_Theta3(3,:,lmc),'m--x'); hold on; % FQ
    temp5  = squeeze(mean(crb_Theta5(3,:,:,lmc),3));
    temp6  = squeeze(mean(crb_Theta6(3,:,:,lmc),3));
    plotc1 = semilogy(SNRset,temp5,'b-o'); hold on; % TQ
    plotc2 = semilogy(SNRset,temp6,'r-s'); hold on; % TQ
end
xlim([SNRset(1),SNRset(end)]); grid on;
xlabel('SNR (dB)');
title(sprintf('AoA(3)'));

subplot(2,2,4)
for lmc =1:nMC
    plotb1 = semilogy(SNRset,crb_Theta1(4,:,lmc),'k-+');  hold on; % OQ
    plotb2 = semilogy(SNRset,crb_Theta3(4,:,lmc),'m--x'); hold on; % FQ
    temp5  = squeeze(mean(crb_Theta5(4,:,:,lmc),3));
    temp6  = squeeze(mean(crb_Theta6(4,:,:,lmc),3));
    plotc1 = semilogy(SNRset,temp5,'b-o'); hold on; % TQ
    plotc2 = semilogy(SNRset,temp6,'r-s'); hold on; % TQ
end
xlim([SNRset(1),SNRset(end)]); grid on;
xlabel('SNR (dB)');
title(sprintf('AoA(4)'));

% AoD (Fig. 3 (b))
figure(200);
subplot(2,2,1)
for lmc =1:nMC
    plotb1 = semilogy(SNRset,crb_Phi1(1,:,lmc),'k-+');  hold on; % OQ
    plotb2 = semilogy(SNRset,crb_Phi3(1,:,lmc),'m--x'); hold on; % FQ
    temp5  = squeeze(mean(crb_Phi5(1,:,:,lmc),3));
    temp6  = squeeze(mean(crb_Phi6(1,:,:,lmc),3));
    plotc1 = semilogy(SNRset,temp5,'b-o'); hold on; % TQ
    plotc2 = semilogy(SNRset,temp6,'r-s'); hold on; % TQ
end
xlim([SNRset(1),SNRset(end)]); grid on;
xlabel('SNR (dB)');
title(sprintf('K=%g, AoD(1)', K));
legend([plotb1(1), plotb2(1), plotc1(1),plotc2(1)],'CQ (known \sigma^2)','FQ (known \sigma^2)','TQ (unknown \sigma^2)','TQ (known \sigma^2)')

subplot(2,2,2)
for lmc =1:nMC
    plotb1 = semilogy(SNRset,crb_Phi1(2,:,lmc),'k-+');  hold on; % OQ
    plotb2 = semilogy(SNRset,crb_Phi3(2,:,lmc),'m--x'); hold on; % FQ
    temp5  = squeeze(mean(crb_Phi5(2,:,:,lmc),3));
    temp6  = squeeze(mean(crb_Phi6(2,:,:,lmc),3));
    plotc1 = semilogy(SNRset,temp5,'b-o'); hold on; % TQ
    plotc2 = semilogy(SNRset,temp6,'r-s'); hold on; % TQ
end
xlim([SNRset(1),SNRset(end)]); grid on;
xlabel('SNR (dB)');
title(sprintf('AoD(2)'));

subplot(2,2,3)
for lmc =1:nMC
    plotb1 = semilogy(SNRset,crb_Phi1(3,:,lmc),'k-+');  hold on; % OQ
    plotb2 = semilogy(SNRset,crb_Phi3(3,:,lmc),'m--x'); hold on; % FQ
    temp5  = squeeze(mean(crb_Phi5(3,:,:,lmc),3));
    temp6  = squeeze(mean(crb_Phi6(3,:,:,lmc),3));
    plotc1 = semilogy(SNRset,temp5,'b-o'); hold on; % TQ
    plotc2 = semilogy(SNRset,temp6,'r-s'); hold on; % TQ
end
xlim([SNRset(1),SNRset(end)]); grid on;
xlabel('SNR (dB)');
title(sprintf('AoD(3)'));

subplot(2,2,4)
for lmc =1:nMC
    plotb1 = semilogy(SNRset,crb_Phi1(4,:,lmc),'k-+');  hold on; % OQ
    plotb2 = semilogy(SNRset,crb_Phi3(4,:,lmc),'m--x'); hold on; % FQ
    temp5  = squeeze(mean(crb_Phi5(4,:,:,lmc),3));
    temp6  = squeeze(mean(crb_Phi6(4,:,:,lmc),3));
    plotc1 = semilogy(SNRset,temp5,'b-o'); hold on; % TQ
    plotc2 = semilogy(SNRset,temp6,'r-s'); hold on; % TQ
end
xlim([SNRset(1),SNRset(end)]); grid on;
xlabel('SNR (dB)');
title(sprintf('AoD(4)'));

% Real-part Channel Path Gain (Fig. 3 (c))
figure(300);
subplot(2,2,1)
for lmc =1:nMC
    plotb1 = semilogy(SNRset,crb_AlphaR1(1,:,lmc),'k-+');  hold on; % OQ
    plotb2 = semilogy(SNRset,crb_AlphaR3(1,:,lmc),'m--x'); hold on; % FQ
    temp5  = squeeze(mean(crb_AlphaR5(1,:,:,lmc),3));
    temp6  = squeeze(mean(crb_AlphaR6(1,:,:,lmc),3));
    plotc1 = semilogy(SNRset,temp5,'b-o'); hold on; % TQ
    plotc2 = semilogy(SNRset,temp6,'r-s'); hold on; % TQ
end
xlim([SNRset(1),SNRset(end)]); grid on;
xlabel('SNR (dB)');
title(sprintf('K=%g, Real Amp(1)', K));
legend([plotb1(1), plotb2(1), plotc1(1),plotc2(1)],'CQ (known \sigma^2)','FQ (known \sigma^2)','TQ (unknown \sigma^2)','TQ (known \sigma^2)')

subplot(2,2,2)
for lmc =1:nMC
    plotb1 = semilogy(SNRset,crb_AlphaR1(2,:,lmc),'k-+');  hold on; % OQ
    plotb2 = semilogy(SNRset,crb_AlphaR3(2,:,lmc),'m--x'); hold on; % FQ
    temp5  = squeeze(mean(crb_AlphaR5(2,:,:,lmc),3));
    temp6  = squeeze(mean(crb_AlphaR6(2,:,:,lmc),3));
    plotc1 = semilogy(SNRset,temp5,'b-o'); hold on; % TQ
    plotc2 = semilogy(SNRset,temp6,'r-s'); hold on; % TQ
end
xlim([SNRset(1),SNRset(end)]); grid on;
xlabel('SNR (dB)');
title(sprintf('Real Amp(2)'));

subplot(2,2,3)
for lmc =1:nMC
    plotb1 = semilogy(SNRset,crb_AlphaR1(3,:,lmc),'k-+');  hold on; % OQ
    plotb2 = semilogy(SNRset,crb_AlphaR3(3,:,lmc),'m--x'); hold on; % FQ
    temp5  = squeeze(mean(crb_AlphaR5(3,:,:,lmc),3));
    temp6  = squeeze(mean(crb_AlphaR6(3,:,:,lmc),3));
    plotc1 = semilogy(SNRset,temp5,'b-o'); hold on; % TQ
    plotc2 = semilogy(SNRset,temp6,'r-s'); hold on; % TQ
end
xlim([SNRset(1),SNRset(end)]); grid on;
xlabel('SNR (dB)');
title(sprintf('Real Amp(3)'));

subplot(2,2,4)
for lmc =1:nMC
    plotb1 = semilogy(SNRset,crb_AlphaR1(4,:,lmc),'k-+');  hold on; % OQ
    plotb2 = semilogy(SNRset,crb_AlphaR3(4,:,lmc),'m--x'); hold on; % FQ
    temp5  = squeeze(mean(crb_AlphaR5(4,:,:,lmc),3));
    temp6  = squeeze(mean(crb_AlphaR6(4,:,:,lmc),3));
    plotc1 = semilogy(SNRset,temp5,'b-o'); hold on; % TQ
    plotc2 = semilogy(SNRset,temp6,'r-s'); hold on; % TQ
end
xlim([SNRset(1),SNRset(end)]); grid on;
xlabel('SNR (dB)');
title(sprintf('Real Amp(4)'));

% Imaginary-part Channel Path Gain (Fig. 3 (d))
figure(310);
subplot(2,2,1)
for lmc =1:nMC
    plotb1 = semilogy(SNRset,crb_AlphaI1(1,:,lmc),'k-+');  hold on; % OQ
    plotb2 = semilogy(SNRset,crb_AlphaI3(1,:,lmc),'m--x'); hold on; % FQ
    temp5  = squeeze(mean(crb_AlphaI5(1,:,:,lmc),3));
    temp6  = squeeze(mean(crb_AlphaI6(1,:,:,lmc),3));
    plotc1 = semilogy(SNRset,temp5,'b-o'); hold on; % TQ
    plotc2 = semilogy(SNRset,temp6,'r-s'); hold on; % TQ
end
xlim([SNRset(1),SNRset(end)]); grid on;
xlabel('SNR (dB)');
title(sprintf('K=%g, Imag Amp(1)', K));
legend([plotb1(1), plotb2(1), plotc1(1),plotc2(1)],'CQ (known \sigma^2)','FQ (known \sigma^2)','TQ (unknown \sigma^2)','TQ (known \sigma^2)')

subplot(2,2,2)
for lmc =1:nMC
    plotb1 = semilogy(SNRset,crb_AlphaI1(2,:,lmc),'k-+');  hold on; % OQ
    plotb2 = semilogy(SNRset,crb_AlphaI3(2,:,lmc),'m--x'); hold on; % FQ
    temp5  = squeeze(mean(crb_AlphaI5(2,:,:,lmc),3));
    temp6  = squeeze(mean(crb_AlphaI6(2,:,:,lmc),3));
    plotc1 = semilogy(SNRset,temp5,'b-o'); hold on; % TQ
    plotc2 = semilogy(SNRset,temp6,'r-s'); hold on; % TQ
end
xlim([SNRset(1),SNRset(end)]); grid on;
xlabel('SNR (dB)');
title(sprintf('Imag Amp(2)'));

subplot(2,2,3)
for lmc =1:nMC
    plotb1 = semilogy(SNRset,crb_AlphaI1(3,:,lmc),'k-+');  hold on; % OQ
    plotb2 = semilogy(SNRset,crb_AlphaI3(3,:,lmc),'m--x'); hold on; % FQ
    temp5  = squeeze(mean(crb_AlphaI5(3,:,:,lmc),3));
    temp6  = squeeze(mean(crb_AlphaI6(3,:,:,lmc),3));
    plotc1 = semilogy(SNRset,temp5,'b-o'); hold on; % TQ
    plotc2 = semilogy(SNRset,temp6,'r-s'); hold on; % TQ
end
xlim([SNRset(1),SNRset(end)]); grid on;
xlabel('SNR (dB)');
title(sprintf('Imag Amp(3)'));

subplot(2,2,4)
for lmc =1:nMC
    plotb1 = semilogy(SNRset,crb_AlphaI1(4,:,lmc),'k-+');  hold on; % OQ
    plotb2 = semilogy(SNRset,crb_AlphaI3(4,:,lmc),'m--x'); hold on; % FQ
    temp5  = squeeze(mean(crb_AlphaI5(4,:,:,lmc),3));
    temp6  = squeeze(mean(crb_AlphaI6(4,:,:,lmc),3));
    plotc1 = semilogy(SNRset,temp5,'b-o'); hold on; % TQ
    plotc2 = semilogy(SNRset,temp6,'r-s'); hold on; % TQ
end
xlim([SNRset(1),SNRset(end)]); grid on;
xlabel('SNR (dB)');
title(sprintf('Imag Amp(4)'));
