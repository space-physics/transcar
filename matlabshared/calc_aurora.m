function calc_aurora(datadir,infile,outfile)
%
% Algorithm by Matt Zettergren 2008
% Algorithm description by Michael Hirsch 2014
%
% ALGORITHM:
% -1) run TRANSCAR, which saves emissions.dat
% emissions.dat format:
%  0) run process_transcar(), which reads emissions.dat and converts
%  variables to Matlab format, and resaves the variables into excrates.mat,
%  a Matlab HDF5 data file.
%  1) 

%% LOAD EXCITATION/POPULATION RATES OBTAINED DIRECTLY FROM KINETIC CODE
kinFN = [datadir,'/',infile]; %excrates.mat
display(['calc_aurora: Loading ',kinFN])
load(kinFN,'dipangle','zop','n*','timeop','p*','e*','flux*')

%% DISTANCE ALONG FLUX TUBE IN MAGNETIC ZENITH
% convert altitudes (km) to distance (cm) along flux tube
meanDipAngDeg = mean(dipangle);
display(['Mean dip angle = ',num2str(meanDipAngDeg,'%0.2f'),...
     ' [deg] , std dev dip angle: ',num2str(std(meanDipAngDeg),'%0.3f'),' [deg]'])
%r = zop .* secd(meanDipAngDeg);  %columns of zop are the same!
r = zop(:,1) .* secd(meanDipAngDeg); %vector of r
r = r * 1e5; %km to cm

%% METASTABLE EMISSION RATES
%display('calc_aurora: Computing metastable emissions.')
Ams=[1.215, 0.076, 7.1e-3, 2.2e-3, 0.173, 0.047];     %Einstein coeff. (1/s)
lambdams=[557.7; 297.2; 630.0; 636.4; 732; 247.0];    %wavelength (nm)

for k=1:length(lambdams)
   if k<=2
     ppl=no1s;
   elseif k<=4
     ppl=no1d;
   else
     ppl=noii2p;  
   end
       
   VER=Ams(k)*ppl;
   %[integ,dr]=intrap(VER,r);     %integrate along flux tube
   %[nrow,ncol]=size(integ);
   %blambdams(k,:)=integ(nrow,:); %this array has wavelength (dim 1) and time (dim 2)
   blambdams(k,:) = trapz(r,VER);    %TODO: verify that trapz was an appropriate replacment
end 

%add the metastable emissions to the group
lambda=lambdams;
blambda=blambdams;
%% PROMPT ATOMIC EMISSIONS
%display('calc_aurora: Computing prompt atomic emissions.')
lambdapa=[844.6; 777.4];

for k=1:length(lambdapa)
   if k==1
     ppl=po3p3p;
   else
     ppl=po3p5p;  
   end
       
   VER=ppl;
   %[integ,dr]=intrap(VER,r);       %integrate along flux tube
   %[nrow,ncol]=size(integ);
   %blambdapa(k,:)=integ(nrow,:);   %this array has wavelength (dim 1) and time (dim 2)
   blambdapa(k,:)=trapz(r,VER); 
end 

%add the prompt atomic emissions to the group
lambda=[lambda;lambdapa];
blambda=[blambda;blambdapa];
%% PROMPT BAND EMISSIONS
%N2+ 1NG (from Vallance Jones, 1974)
%display('CALC_AURORA: Computing N2+ 1NG emissions.')
A1NG=[.11e8 , .335e7, .675e6, .111e6, .161e5,  0.00;               %dim 1 is upper state vib. level (nu'); dim 2 is bottom state vib level (nu'')
      .67e7 , .387e7, .381e7, .132e7, .312e6, .590e5;
      .127e7, .916e7, .899e6, .313e7, .171e7, .543e6];
lambda1NG=[391.4, 427.8, 470.9, 522.8, 586.5, 000.0;
           358.2, 388.4, 423.7, 465.2, 514.9, 575.4;
           330.8, 356.3, 385.8, 419.9, 460.0, 507.7];
FCfact1NG=[0.8911; 0.1071; 0.0017];                          %excitation Franck-Condon factors (derived from Vallance Jones, 1974)
tau1NG=1./sum(A1NG,2);
FCfactmat=repmat(FCfact1NG,1,size(A1NG,2));
taumat=repmat(tau1NG,1,size(A1NG,2));
scalemat=A1NG.*taumat.*FCfactmat;

%flatten the matrices for computation
lambda1NGvec=lambda1NG(:);
scalevec=scalemat(:);

%Compute vibrational bands
for k=1:length(scalevec)
    VER=p1ng*scalevec(k);
    %[integ,dr]=intrap(VER,r);                                 %integrate along flux tube
    %[nrow,ncol]=size(integ);
    %b1NGvec(k,:)=integ(nrow,:);
    b1NGvec(k,:) = trapz(r,VER);
end

%Add the 1NG to the rest of the emissions
lambda=[lambda;lambda1NGvec];
blambda=[blambda; b1NGvec];
%% N2+ Meinel band
%display('CALC_AURORA.M: Computing N2+ Meinel emissions.')
AMein=[.455e5, .136e4, .119e4, .279e2, .170e-1,       0,       0,        0,      0,      0,      0;
       .553e5, .352e4, .112e5, .202e4, .675e2 , .216e-1,       0,        0,      0,      0,      0;
       .373e5, .367e5, .105e4, .573e4, .219e4 ,  .992e2, .767e-2,        0,      0,      0,      0;
       .189e5, .517e5, .128e5, .621e4, .186e4 ,  .186e4,  .112e3, .141e-49,      0,      0,      0;
       .817e4, .407e5, .426e5, .140e4, .882e4 ,  .205e3,  .133e4, .107e3  ,      0,      0,      0;
       .320e4, .239e5, .514e5, .246e5, .528e3 ,  .790e4,  .651e2, .806e3  , .885e2,      0,      0;
       .119e4, .119e5, .405e5, .480e5, .927e4 ,  .405e4,  .518e4, .552e3  , .405e3, .649e2,      0;
       .426e3, .538e4, .252e5, .508e5, .353e5 ,  .137e4,  .728e4, .248e4  , .105e4, .156e3, .424e3];
lambdaMein=[1108.7, 1461.0, 2126.5, 3856.1,      0,      0,      0,      0,      0,      0,      0;
             918.2, 1147.2, 1521.1, 2239.6, 4186.0,      0,      0,      0,      0,      0,      0;
             785.3,  947.0, 1188.0, 1585.2, 2362.9, 4568.0,      0,      0,      0,      0,      0;
             687.4,  808.2,  977.4, 1231.2, 1654.1, 2498.3, 5016.8,      0,      0,      0,      0;
             612.3,  706.4,  832.3, 1009.5, 1277.3, 1728.2, 2647.7, 5552.0,      0,      0,      0;
             552.9,  628.5,  726.3,  857.7, 1043.6, 1326.3, 1808.3, 2813.3, 6199.6,      0,      0;
             504.8,  567.1,  645.5,  747.3,  884.5, 1079.6, 1378.6, 1894.7, 2997.3, 6996.3,      0;
             465.1,  517.4,  581.9,  663.3,  769.3,  912.7, 1117.7, 1434.4, 1988.1, 3202.3, 7995.8];
unFCfactMein=[1.28; 1.49; 1.0; 0.52; 0.23; 0.092; 0.035; 0.013];                   %Piper 1986
FCfactMein=unFCfactMein/sum(unFCfactMein,1);                                        %normalize to obtain absolute Franck-Condon factors
tauMein=1./sum(AMein,2);
FCfactmat=repmat(FCfactMein,1,size(AMein,2));
taumat=repmat(tauMein,1,size(AMein,2));
scalemat=AMein.*taumat.*FCfactmat;

%flatten the matrices for computation
lambdaMeinvec=lambdaMein(:);
scalevec=scalemat(:);

%Compute vibrational bands
for k=1:length(scalevec)
    VER=pmein*scalevec(k);
    %[integ,dr]=intrap(VER,r);                                 %integrate along flux tube
    %[nrow,ncol]=size(integ);
    %bMeinvec(k,:)=integ(nrow,:);
    bMeinvec(k,:) = trapz(r,VER);
end

%Add the 2PG to the rest of the emissions
lambda=[lambda;lambdaMeinvec];
blambda=[blambda; bMeinvec];
%% N2 2PG (after Vallance Jones, 1974)
%display('CALC_AURORA.M: Computing N2 2PG emissions.')
A2PG=[.11e8 , .733e7, .294e7, .923e6, .247e6, .606e5,      0,      0,      0,      0;
      .102e8, .528e6, .461e7, .410e7, .249e7, .769e6, .235e6, .674e5,      0,      0;
      .349e7, .844e7, .598e6, .146e7, .377e7, .263e7, .132e7, .526e6, .182e6,      0;
      .507e6, .661e7, .548e7, .218e7, .112e6, .201e7, .248e7, .169e7, .807e6, .333e6;
      .184e5, .133e7, .864e7, .324e7, .201e7, .912e5, .911e6, .195e7, .173e7, .101e7];         %from Vallance Jones, 1974 (Shemansky and Broadfoot, 1971b)
lambda2PG=[337.0, 357.6, 380.4, 405.8, 434.3, 466.5,     0,     0,     0,     0;
           315.8, 333.8, 353.6, 375.4, 399.7, 426.8, 457.3, 491.7,     0,     0;
           297.6, 313.5, 330.9, 349.9, 370.9, 394.2, 420.0, 448.9, 481.3,     0;
           281.8, 296.1, 311.5, 328.4, 346.8, 367.1, 389.4, 414.0, 441.5, 472.2;
           268.4, 281.2, 295.2, 310.2, 326.6, 344.5, 364.1, 385.6, 409.3, 435.5];                   %from Vallance Jones, 1974
FCfact2PG=[.5466; .3050; .1057; .2963e-1; .7573e-2];            %from Benesch et al, 1966a
tau2PG=1./sum(A2PG,2);
FCfactmat=repmat(FCfact2PG,1,size(A2PG,2));
taumat=repmat(tau2PG,1,size(A2PG,2));
scalemat=A2PG.*taumat.*FCfactmat;

%flatten the matrices for computation
lambda2PGvec=lambda2PG(:);
scalevec=scalemat(:);

%Compute vibrational bands
p2PGvec=[];
for k=1:length(scalevec)
    VER=p2pg*scalevec(k);
    %[integ,dr]=intrap(VER,r);                                 %integrate along flux tube
    %[nrow,ncol]=size(integ);
    %b2PGvec(k,:)=integ(nrow,:);
    b2PGvec(k,:) = trapz(r,VER);
end

%Add the 2PG to the rest of the emissions
lambda=[lambda;lambda2PGvec];
blambda=[blambda; b2PGvec];
%% N2 1PG (this one is a bit tricky)
%display('CALC_AURORA: Computing N2 1PG emissions.')
A1PG=[.625e5, .356e5, .112e5, .247e4, .397e3, .424e2,      0,      0,      0,      0,      0,      0,      0;
      .872e5, .412e3, .185e5, .148e5, .569e4, .140e4, .227e3,      0,      0,      0,      0,      0,      0;
      .444e5, .617e5, .125e5, .268e4, .105e5, .729e4, .217e4, .639e3, .887e2,      0,      0,      0,      0;
      .107e5, .773e5, .217e5, .285e5, .797e3, .388e4, .633e4, .368e4, .124e4, .257e3, .254e2,      0,      0;
      .129e4, .302e5, .836e5, .154e4, .294e5, .784e4, .193e3, .152e4, .374e4, .184e4, .544e3, .903e2,      0;
      .742e2, .508e4, .526e5, .686e5, .302e4, .191e5, .152e5, .110e4, .102e4, .281e4, .213e4, .895e3, .224e3;
           0, .371e3, .119e5, .709e5, .438e5, .158e5, .710e4, .176e5, .511e4,      0, .140e4, .193e4, .117e4;
           0,      0, .107e4, .214e5, .810e5, .203e5, .288e5, .509e3, .143e5, .939e4, .112e4, .231e3, .126e4;
           0,      0,      0, .234e4, .325e5, .815e5, .494e4, .350e5, .138e4, .811e4, .115e5, .374e4, .933e2;
           0,      0,      0, .777e2, .429e4, .438e5, .731e5,      0, .328e5, .736e4, .244e4, .104e5, .647e4;
           0,      0,      0,      0, .159e3, .692e4, .538e5, .586e5, .399e4, .244e5, .154e5,      0, .696e4;
           0,      0,      0,      0,      0, .282e3, .101e5, .611e5, .414e5, .133e5, .138e5, .209e5, .158e4;
           0,      0,      0,      0,      0,      0, .449e3, .138e5, .649e5, .248e5, .237e5, .486e4, .219e5];
lambda1PG=[1050.8, 1237.3, 1497.7, 1887.4,      0,      0,      0,      0,      0,      0,      0,0,0;
            891.2,      0, 1192.5, 1426.9, 1768.7, 2306.0,      0,      0,      0,      0,      0,0,0;
            775.4,  872.3,  994.0, 1151.6, 1363.5, 1664.3, 2121.0,      0,      0,      0,      0,0,0;
            687.5,  762.6,  854.2,  968.0,      0, 1306.1, 1571.7, 1961.9,      0,      0,      0,0,0;
            618.7,  679.0,  750.5,  837.0,  943.6, 1078.1,      0, 1489.0,      0,      0,      0,0,0;
                0,  612.7,  670.5,  738.7,  820.6,  920.4, 1044.8,      0,      0, 1706.3,      0,0,0;
                0,      0,  607.0,  662.4,  727.4,  804.8,  898.3, 1013.3, 1159.0,      0,      0,0,0;
                0,      0,      0,  601.4,  654.5,  716.5,  789.7,      0,  984.2, 1117.0,      0,0,0;
                0,      0,      0,      0,  595.9,  646.9,  706.0,  775.2,      0,      0, 1078.0,0,0;
                0,      0,      0,      0,      0,  590.6,  639.8,      0,  761.3,      0,      0,0,0;
                0,      0,      0,      0,      0,      0,  584.5,  632.3,      0,      0,      0,0,0;
                0,      0,      0,      0,      0,      0,      0,  580.4,  625.3,      0,      0,0,0;
                0,      0,      0,      0,      0,      0,      0,      0,  575.5,  618.5,      0,0,0];
tau1PG = 1 ./ sum(A1PG,2);
%% solve for base concentration
%confac=[1.66;1.56;1.31;1.07;.77;.5;.33;.17;.08;.04;.02;.004;.001];             %Cartwright, 1973b, stop at nuprime==12
confac=[1.66;1.86;1.57;1.07;.76;.45;.25;.14;.07;.03;.01;.004;.001];             %Gattinger and Vallance Jones 1974
consfac=confac/sum(confac,1);            %normalize
losscoef=sum(consfac ./ tau1PG,1);
N01pg=p1pg ./ losscoef;
consmat=repmat(consfac,1,size(A1PG,2));
scalemat=A1PG.*consmat;
%% flatten matrices for computation
lambda1PGvec=lambda1PG(:);
scalevec=scalemat(:);

%Compute vibrational bands for 1PG
for k=1:length(scalevec)
    VER=N01pg*scalevec(k);
    %[integ,dr]=intrap(VER,r);   %integrate along flux tube
    %[nrow,ncol]=size(integ);
    %b1PGvec(k,:)=integ(nrow,:);
    b1PGvec(k,:)=trapz(r,VER);
end

%add the 1PG to the rest of the emissions
lambda=[lambda;lambda1PGvec];
blambda=[blambda; b1PGvec];
%% GET RID OF ZERO WAVELENGTH ENTRIES (THAT WERE INCLUDED FOR COMPUTATIONAL EASE)
display('Cleaning up data arrays.')
inds=find(lambda~=0);
blambda=blambda(inds,:);
lambda=lambda(inds);
%% CONVERT TO RAYLEIGH UNITS
blambda=blambda/1e6;
%% SAVE TO FILE
display(['CALC_AURORA: Saving calculations to: ',outfile]);
save([datadir,'/',outfile],'e*','fluxdown*','lambda','lambdams',...
    'lambdapa','lambda1NGvec','lambdaMeinvec','lambda2PGvec',...
    'lambda1PGvec','blambda','timeop','zop')
end