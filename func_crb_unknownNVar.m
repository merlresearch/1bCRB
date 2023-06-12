% Copyright (C) 2018,2023 Mitsubishi Electric Research Laboratories (MERL)
%
% SPDX-License-Identifier: AGPL-3.0-or-later

function tempFIM = func_crb_unknownNVar(thresR, thresI, thetaVec, phiVec, Gamma, hVec, S, A_ms, A_bs, N_tx, N_rx, M, nvar)
% function to compute FIM for angular-domain channel estimation
% (AoA, AoD, and channel path gain) with 1-bit ADC when the noise variance
% is unknown
%
% M: the number of channel paths
%
% Refer to
% P. Wang, J. Li, M. Pajovic, P.T. Boufounos, P.V. Orlik,
% "On Angular-Domain Chanel Estimation for One-Bit Massive MIMO Systems with Fixed and Time-Varying Thresholds",
% 2017 Asilomar Conference on Signals, Systems and Computers (ACSSC), Pacific Grove, CA, November, 2017.

% preliminary
GammaR = real(Gamma);
GammaI = imag(Gamma);
hR     = real(hVec);
hI     = imag(hVec);
SR     = real(S);
SI     = imag(S);
AmsR   = real(A_ms);
AmsI   = imag(A_ms);
AbsR   = real(A_bs);
AbsI   = imag(A_bs);

eta1   = (GammaR*hR-GammaI*hI-thresR)/sqrt(nvar/2);
eta2   = (GammaR*hI+GammaI*hR-thresI)/sqrt(nvar/2);

% replace by erfcx to avoid underflow or overflow errors
tempConstA = ( 2./(erfcx(-eta1/sqrt(2))) + 2./(erfcx(eta1/sqrt(2))) ).*exp(-eta1.^2/2);
tempConstB = ( 2./(erfcx(-eta2/sqrt(2))) + 2./(erfcx(eta2/sqrt(2))) ).*exp(-eta2.^2/2);

for m=1:M

    temp1      = zeros(N_tx, M);
    temp1(:,m) = -sin((0:N_tx-1)*pi*sind(phiVec(m))).*(0:N_tx-1)*pi*cosd(phiVec(m))/sqrt(N_tx);
    partial_AmsR_phi{m} = temp1;

    temp2      = zeros(N_tx, M);
    temp2(:,m) = cos((0:N_tx-1)*pi*sind(phiVec(m))).*(0:N_tx-1)*pi*cosd(phiVec(m))/sqrt(N_tx);
    partial_AmsI_phi{m} = temp2;

    temp1b      = zeros(N_rx, M);
    temp1b(:,m) = -sin((0:N_rx-1)*pi*sind(thetaVec(m))).*(0:N_rx-1)*pi*cosd(thetaVec(m))/sqrt(N_rx);
    partial_AbsR_theta{m} = temp1b;

    temp2b      = zeros(N_rx, M);
    temp2b(:,m) = cos((0:N_rx-1)*pi*sind(thetaVec(m))).*(0:N_rx-1)*pi*cosd(thetaVec(m))/sqrt(N_rx);
    partial_AbsI_theta{m} = temp2b;

    temp3a       = SR.'*partial_AmsR_phi{m};
    temp3b       = SI.'*partial_AmsI_phi{m};
    temp3c       = SR.'*partial_AmsI_phi{m};
    temp3d       = SI.'*partial_AmsR_phi{m};

    temp5a       = SR.'*AmsR;
    temp5b       = SI.'*AmsI;
    temp5c       = SR.'*AmsI;
    temp5d       = SI.'*AmsR;

    % updated for asilomar
    temp_GammaR_phi{m} = kron(temp3a,AbsR)+kron(temp3b,AbsR)+kron(temp3c,AbsI)-kron(temp3d,AbsI);
    temp_GammaI_phi{m} = kron(temp3a,AbsI)+kron(temp3b,AbsI)-kron(temp3c,AbsR)+kron(temp3d,AbsR);

    temp_GammaR_theta{m} = kron(temp5a,partial_AbsR_theta{m})+kron(temp5b,partial_AbsR_theta{m})...
        + kron(temp5c,partial_AbsI_theta{m})-kron(temp5d,partial_AbsI_theta{m});
    temp_GammaI_theta{m} = kron(temp5a,partial_AbsI_theta{m})+kron(temp5b,partial_AbsI_theta{m})...
        - kron(temp5c,partial_AbsR_theta{m})+kron(temp5d,partial_AbsR_theta{m});

end

% FIM for theta: auto-elements
FIM_ThetaTheta = zeros(M,M);

for m1 = 1:M
    tempGammaR1 = temp_GammaR_theta{m1};
    tempGammaI1 = temp_GammaI_theta{m1};

    for m2 = 1:M

        tempGammaR2 = temp_GammaR_theta{m2};
        tempGammaI2 = temp_GammaI_theta{m2};

        tempThetaTheta = tempConstA.*(tempGammaR1*hR - tempGammaI1*hI).*(tempGammaR2*hR-tempGammaI2*hI)+...
                         tempConstB.*(tempGammaR1*hI + tempGammaI1*hR).*(tempGammaR2*hI+tempGammaI2*hR);
        FIM_ThetaTheta(m1,m2) = sum(tempThetaTheta)/pi/nvar;

    end
end


% FIM for theta-phi: cross-elements
FIM_ThetaPhi = zeros(M,M);
FIM_PhiTheta = zeros(M,M);

for m1 = 1:M
    tempGammaR1 = temp_GammaR_theta{m1};
    tempGammaI1 = temp_GammaI_theta{m1};

    for m2 = 1:M

        tempGammaR2 = temp_GammaR_phi{m2};
        tempGammaI2 = temp_GammaI_phi{m2};

        tempThetaPhi = tempConstA.*(tempGammaR1*hR - tempGammaI1*hI).*(tempGammaR2*hR-tempGammaI2*hI)+...
                       tempConstB.*(tempGammaR1*hI + tempGammaI1*hR).*(tempGammaR2*hI+tempGammaI2*hR);
        FIM_ThetaPhi(m1,m2) = sum(tempThetaPhi)/pi/nvar;

    end
end
FIM_PhiTheta = FIM_ThetaPhi.';

% FIM for theta-alphaR: cross-elements
FIM_ThetaAlphaR = zeros(M,M);
FIM_AlphaRTheta = zeros(M,M);

for m1 = 1:M
    tempGammaR1 = temp_GammaR_theta{m1};
    tempGammaI1 = temp_GammaI_theta{m1};

    for m2 = 1:M

        eVec                  = zeros(M*M,1);
        tempVec               = zeros(M,1);
        tempVec(m2)           = 1;
        eVec((m2-1)*M+1:m2*M) = tempVec;

        tempGammaR2 = GammaR;
        tempGammaI2 = GammaI;

        tempThetaAlphaR    = tempConstA.*(tempGammaR1*hR - tempGammaI1*hI).*(tempGammaR2*eVec)+...
                             tempConstB.*(tempGammaR1*hI + tempGammaI1*hR).*(tempGammaI2*eVec);
        FIM_ThetaAlphaR(m1,m2) = sum(tempThetaAlphaR)/pi/nvar;

    end
end

FIM_AlphaRTheta = FIM_ThetaAlphaR.';

% FIM for theta-alphaI: cross-elements
FIM_ThetaAlphaI = zeros(M,M);
FIM_AlphaITheta = zeros(M,M);

for m1 = 1:M
    tempGammaR1 = temp_GammaR_theta{m1};
    tempGammaI1 = temp_GammaI_theta{m1};

    for m2 = 1:M

        eVec                  = zeros(M*M,1);
        tempVec               = zeros(M,1);
        tempVec(m2)           = 1;
        eVec((m2-1)*M+1:m2*M) = tempVec;

        tempGammaR2 = GammaR;
        tempGammaI2 = GammaI;

        tempThetaAlphaI = tempConstA.*(tempGammaR1*hR - tempGammaI1*hI).*(-tempGammaI2*eVec)+...
                          tempConstB.*(tempGammaR1*hI + tempGammaI1*hR).*(tempGammaR2*eVec);
        FIM_ThetaAlphaI(m1,m2) = sum(tempThetaAlphaI)/pi/nvar;

    end
end
FIM_AlphaITheta = FIM_ThetaAlphaI.';

% FIM for phi: auto-elements
FIM_PhiPhi = zeros(M,M);

for m1 = 1:M
    tempGammaR1 = temp_GammaR_phi{m1};
    tempGammaI1 = temp_GammaI_phi{m1};

    for m2 = 1:M

        tempGammaR2 = temp_GammaR_phi{m2};
        tempGammaI2 = temp_GammaI_phi{m2};

        tempPhiPhi = tempConstA.*(tempGammaR1*hR - tempGammaI1*hI).*(tempGammaR2*hR-tempGammaI2*hI)+...
                     tempConstB.*(tempGammaR1*hI + tempGammaI1*hR).*(tempGammaR2*hI+tempGammaI2*hR);
        FIM_PhiPhi(m1,m2) = sum(tempPhiPhi)/pi/nvar;

    end
end

% FIM for phi-alphaR: cross-elements
FIM_PhiAlphaR = zeros(M,M);
FIM_AlphaRPhi = zeros(M,M);

for m1 = 1:M
    tempGammaR1 = temp_GammaR_phi{m1};
    tempGammaI1 = temp_GammaI_phi{m1};

    for m2 = 1:M

        eVec                  = zeros(M*M,1);
        tempVec               = zeros(M,1);
        tempVec(m2)           = 1;
        eVec((m2-1)*M+1:m2*M) = tempVec;

        tempGammaR2 = GammaR;
        tempGammaI2 = GammaI;

        tempPhiAlphaR = tempConstA.*(tempGammaR1*hR - tempGammaI1*hI).*(tempGammaR2*eVec)+...
                        tempConstB.*(tempGammaR1*hI + tempGammaI1*hR).*(tempGammaI2*eVec);
        FIM_PhiAlphaR(m1,m2) = sum(tempPhiAlphaR)/pi/nvar;

    end
end
FIM_AlphaRPhi = FIM_PhiAlphaR.';

% FIM for phi-alphaI: cross-elements
FIM_PhiAlphaI = zeros(M,M);
FIM_AlphaIPhi = zeros(M,M);

for m1 = 1:M
    tempGammaR1 = temp_GammaR_phi{m1};
    tempGammaI1 = temp_GammaI_phi{m1};

    for m2 = 1:M

        eVec                  = zeros(M*M,1);
        tempVec               = zeros(M,1);
        tempVec(m2)           = 1;
        eVec((m2-1)*M+1:m2*M) = tempVec;

        tempGammaR2 = GammaR;
        tempGammaI2 = GammaI;

        tempPhiAlphaI = tempConstA.*(tempGammaR1*hR - tempGammaI1*hI).*(-tempGammaI2*eVec)+...
                        tempConstB.*(tempGammaR1*hI + tempGammaI1*hR).*(tempGammaR2*eVec);
        FIM_PhiAlphaI(m1,m2) = sum(tempPhiAlphaI)/pi/nvar;

    end
end
FIM_AlphaIPhi = FIM_PhiAlphaI.';

% FIM for alphaR-alphaR: auto-elements
FIM_AlphaRAlphaR = zeros(M,M);

for m1 = 1:M

    eVec1                 = zeros(M*M,1);
    tempVec               = zeros(M,1);
    tempVec(m1)           = 1;
    eVec1((m1-1)*M+1:m1*M) = tempVec;

    tempGammaR1 = GammaR;
    tempGammaI1 = GammaI;

    for m2 = 1:M

        eVec2                 = zeros(M*M,1);
        tempVec               = zeros(M,1);
        tempVec(m2)           = 1;
        eVec2((m2-1)*M+1:m2*M) = tempVec;

        tempGammaR2 = GammaR;
        tempGammaI2 = GammaI;

        tempAlphaRAlphaR = tempConstA.*(tempGammaR1*eVec1).*(tempGammaR2*eVec2)+...
                           tempConstB.*(tempGammaI1*eVec1).*(tempGammaI2*eVec2);
        FIM_AlphaRAlphaR(m1,m2) = sum(tempAlphaRAlphaR)/pi/nvar;

    end
end


% FIM for alphaR-alphaI: cross-elements
FIM_AlphaRAlphaI = zeros(M,M);
FIM_AlphaIAlphaR = zeros(M,M);

for m1 = 1:M

    eVec1                 = zeros(M*M,1);
    tempVec               = zeros(M,1);
    tempVec(m1)           = 1;
    eVec1((m1-1)*M+1:m1*M) = tempVec;

    tempGammaR1 = GammaR;
    tempGammaI1 = GammaI;

    for m2 = 1:M

        eVec2                 = zeros(M*M,1);
        tempVec               = zeros(M,1);
        tempVec(m2)           = 1;
        eVec2((m2-1)*M+1:m2*M) = tempVec;

        tempGammaR2 = GammaR;
        tempGammaI2 = GammaI;

        tempAlphaRAlphaI = tempConstA.*(tempGammaR1*eVec1).*(-tempGammaI2*eVec2)+...
                           tempConstB.*(tempGammaI1*eVec1).*(tempGammaR2*eVec2);
        FIM_AlphaRAlphaI(m1,m2) = sum(tempAlphaRAlphaI)/pi/nvar;

    end
end
FIM_AlphaIAlphaR = FIM_AlphaRAlphaI.';

% FIM for alphaI-alphaI: auto-elements
FIM_AlphaIAlphaI = zeros(M,M);

for m1 = 1:M

    eVec1                 = zeros(M*M,1);
    tempVec               = zeros(M,1);
    tempVec(m1)           = 1;
    eVec1((m1-1)*M+1:m1*M) = tempVec;

    tempGammaR1 = GammaR;
    tempGammaI1 = GammaI;

    for m2 = 1:M

        eVec2                 = zeros(M*M,1);
        tempVec               = zeros(M,1);
        tempVec(m2)           = 1;
        eVec2((m2-1)*M+1:m2*M) = tempVec;

        tempGammaR2 = GammaR;
        tempGammaI2 = GammaI;

        tempAlphaIAlphaI = tempConstA.*(- tempGammaI1*eVec1).*(-tempGammaI2*eVec2)+...
                           tempConstB.*(tempGammaR1*eVec1).*(tempGammaR2*eVec2);
        FIM_AlphaIAlphaI(m1,m2) = sum(tempAlphaIAlphaI)/pi/nvar;

    end
end

% FIM for sigma2-sigma2: auto-element
FIM_Sigma2Sigma2 = zeros(1);

tempSigma2Sigma2 = tempConstA.*(GammaR*hR - GammaI*hI-thresR).^2+...
                   tempConstB.*(GammaI*hR + GammaR*hI-thresI).^2;
FIM_Sigma2Sigma2 = sum(tempSigma2Sigma2)/pi/nvar^3/4;


% FIM for theta-sigma2: cross-elements
FIM_ThetaSigma2 = zeros(M,1);
FIM_Sigma2Theta = zeros(1,M);

for m1 = 1:M
    tempGammaR1 = temp_GammaR_theta{m1};
    tempGammaI1 = temp_GammaI_theta{m1};

    tempThetaSigma2 = tempConstA.*(tempGammaR1*hR - tempGammaI1*hI).*(GammaR*hR - GammaI*hI-thresR)+...
                      tempConstB.*(tempGammaR1*hI + tempGammaI1*hR).*(GammaI*hR + GammaR*hI-thresI);
    FIM_ThetaSigma2(m1) = -sum(tempThetaSigma2)/pi/nvar^2/2;

end
FIM_Sigma2Theta = FIM_ThetaSigma2.';

% FIM for phi-sigma2: cross-elements
FIM_PhiSigma2 = zeros(M,1);
FIM_Sigma2Phi = zeros(1,M);

for m1 = 1:M
    tempGammaR1 = temp_GammaR_phi{m1};
    tempGammaI1 = temp_GammaI_phi{m1};

    tempPhiSigma2 = tempConstA.*(tempGammaR1*hR - tempGammaI1*hI).*(GammaR*hR - GammaI*hI-thresR)+...
                    tempConstB.*(tempGammaR1*hI + tempGammaI1*hR).*(GammaI*hR + GammaR*hI-thresI);
    FIM_PhiSigma2(m1) = -sum(tempPhiSigma2)/pi/nvar^2/2;

end
FIM_Sigma2Phi = FIM_PhiSigma2.';

% FIM for alphaR-sigma2: cross-elements
FIM_AlphaRSigma2 = zeros(M,1);
FIM_Sigma2AlphaR = zeros(1,M);

for m1 = 1:M

    eVec1                 = zeros(M*M,1);
    tempVec               = zeros(M,1);
    tempVec(m1)           = 1;
    eVec1((m1-1)*M+1:m1*M) = tempVec;

    tempGammaR1 = GammaR;
    tempGammaI1 = GammaI;

    tempAlphaRSigma2 = tempConstA.*(tempGammaR1*eVec1).*(GammaR*hR - GammaI*hI-thresR)+...
                       tempConstB.*(tempGammaI1*eVec1).*(GammaI*hR + GammaR*hI-thresI);
    FIM_AlphaRSigma2(m1) = -sum(tempAlphaRSigma2)/pi/nvar^2/2;

end
FIM_Sigma2AlphaR = FIM_AlphaRSigma2.';

% FIM for phi-alphaI: cross-elements
FIM_AlphaISigma2 = zeros(M,1);
FIM_Sigma2AlphaI = zeros(1,M);


for m1 = 1:M

    eVec                  = zeros(M*M,1);
    tempVec               = zeros(M,1);
    tempVec(m1)           = 1;
    eVec((m1-1)*M+1:m1*M) = tempVec;

    tempGammaR1 = GammaR;
    tempGammaI1 = GammaI;

    tempAlphaISigma2 = tempConstA.*(-tempGammaI1*eVec).*(GammaR*hR - GammaI*hI-thresR)+...
                       tempConstB.*(tempGammaR1*eVec).*(GammaI*hR + GammaR*hI-thresI);
    FIM_AlphaISigma2(m1) = -sum(tempAlphaISigma2)/pi/nvar^2/2;

end
FIM_Sigma2AlphaI = FIM_AlphaISigma2.';

tempFIM  = zeros(4*M+1,4*M+1);

tempFIM(1:M,1:M)             = FIM_ThetaTheta;
tempFIM(M+1:2*M,M+1:2*M)     = FIM_PhiPhi;
tempFIM(2*M+1:3*M,2*M+1:3*M) = FIM_AlphaRAlphaR;
tempFIM(3*M+1:4*M,3*M+1:4*M) = FIM_AlphaIAlphaI;
tempFIM(4*M+1,4*M+1)         = FIM_Sigma2Sigma2;

tempFIM(1:M,M+1:2*M) = FIM_ThetaPhi;
tempFIM(M+1:2*M,1:M) = FIM_PhiTheta;

tempFIM(1:M,2*M+1:3*M) = FIM_ThetaAlphaR;
tempFIM(2*M+1:3*M,1:M) = FIM_AlphaRTheta;

tempFIM(1:M,3*M+1:4*M) = FIM_ThetaAlphaI;
tempFIM(3*M+1:4*M,1:M) = FIM_AlphaITheta;

tempFIM(1:M,4*M+1) = FIM_ThetaSigma2;
tempFIM(4*M+1,1:M) = FIM_Sigma2Theta;


tempFIM(M+1:2*M,2*M+1:3*M)     = FIM_PhiAlphaR;
tempFIM(2*M+1:3*M,M+1:2*M)     = FIM_AlphaRPhi;

tempFIM(M+1:2*M,3*M+1:4*M)     = FIM_PhiAlphaI;
tempFIM(3*M+1:4*M,M+1:2*M)     = FIM_AlphaIPhi;

tempFIM(M+1:2*M,4*M+1)         = FIM_PhiSigma2;
tempFIM(4*M+1,M+1:2*M)         = FIM_Sigma2Phi;

tempFIM(2*M+1:3*M,3*M+1:4*M)   = FIM_AlphaRAlphaI;
tempFIM(3*M+1:4*M,2*M+1:3*M)   = FIM_AlphaIAlphaR;

tempFIM(2*M+1:3*M,4*M+1)   = FIM_AlphaRSigma2;
tempFIM(4*M+1,2*M+1:3*M)   = FIM_Sigma2AlphaR;

tempFIM(3*M+1:4*M,4*M+1)   = FIM_AlphaISigma2;
tempFIM(4*M+1,3*M+1:4*M)   = FIM_Sigma2AlphaI;
