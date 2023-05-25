%Cari's run instructions: altVsC14('[scaling scheme]','[ice-surface
%elevation]')
%Make sure you add "Balculator_v2.2" to the path before running
%The figure'll probably be zoomed out to hell, in which case you'll need to
%use View>Property Editor to change the axes such that your data's visible
function altVsC14(sf,ice_surface)
%syntax for sample data file: locnuc#.txt
%syntax for results data file: loc#res.txt
close all;
%Code to plot a sample or set of samples from a site in concentration
%elevation space, and compare with saturation concentration.

%First load variety of constants
% make_constsTgt;
load consts_LSDPavonTgt.mat;
consts.P14_ref_LS = 12.1173; consts.delP14_ref_LS = 0.85841; %LDEO Greenland
consts.P14_ref_LS = 12.0251; consts.delP14_ref_LS = 0.20489; %TU Cronus A
consts.P14_ref_St = 12.2151; consts.delP14_ref_St = 0.21584; %TU Cronus A

file = input('Please enter the name of the sample data file','s');
FID = fopen(file);
data = textscan(FID,'%s %n %n %n %s %n %n %n %n %n %n %n %n %s %s');
fclose(FID);
dstring='';

% NSstring = input('Nuclide-specific scaling? (y/n): ','s');
% 
% if NSstring == 'y'
%     NucSpec = 1;
% else
%     NucSpec = 0;
% end

%Make the sample structure.
%sample.name = input('Paste vector of sample names','s');

all_sample_name = data{1}; all_lat = data{2}; all_long = data{3};
all_elv = data{4}; all_pressure = data{4}; all_aa = data{5}; 
all_thick = data{6}; all_rho = data{7}; all_shielding = data{8};    
all_E = data{9}; all_N14 = data{10}; all_delN14 = data{11}; 
all_N10 = data{12}; all_delN10 = data{13}; all_be_std_name = data{14};
all_curve = data{15};
for i = 1:length(all_N14)
    N14_sample(i) = all_N14(i)./all_shielding(i)./thickness(all_thick(i),160,all_rho(i)); %Correct concentrations to surface (no shielding, zero thickness)
    delN14_sample(i) = N14_sample(i).*(all_delN14(i)./all_N14(i));
end

%Get local saturation concentrations
alts = linspace(0.75.*min(all_elv),1.25.*max(all_elv),30);

for i = 1:length(alts)
    sample.lat = all_lat(1);
    sample.long = all_long(1);
    sample.elv = alts(i);
    sample.aa = all_aa(1);
    sample.thick = 0;
    sample.rho = mean(all_rho);
    sample.shielding = 1;
    sample.E = 0; %mm/ka
    sample.curve = 'none';
    scalingA = ScalingTgtAtm(sample,consts,14,'q',0,0,'pd','n',sample.curve); 
    sf_LS = scalingA.SF_LS;
    sf_St(:,i) = scalingA.SF_St(1);

    tv = scalingA.tv;
    
switch sf
    case 'LSD'
        P14s = sf_LS.*consts.P14_ref_LS;
    case 'St'
        P14s = sf_St.*consts.P14_ref_St;
end

mconsts.Natoms = consts.NatomsQtzO;
mconsts.sigma0 = consts.sigma014;
mconsts.fstar = consts.fstar14;
mconsts.k_negpartial = consts.k_negpartial_14;
mconsts.k_neg = consts.k_negpartial_14;
mconsts.mfluxRef = consts.mfluxRef;

%   Trajectory-traced dipolar estimate for these purposes
dd = [6.89901,-103.241,522.061,-1152.15,1189.18,-448.004;];
RcEst = (dd(1)*cos(d2r(sample.lat)) + ...
    dd(2)*(cos(d2r(sample.lat))).^2 + ...
    dd(3)*(cos(d2r(sample.lat))).^3 + ...
    dd(4)*(cos(d2r(sample.lat))).^4 + ...
    dd(5)*(cos(d2r(sample.lat))).^5 + ...
    dd(6)*(cos(d2r(sample.lat))).^6);

mu_out = P_mu_total_alpha1(0,ERA40atm(sample.lat,sample.long,sample.elv),mconsts,'no');
% mu_out = P_mu_totalLSD(0,ERA40atm(sample.lat,sample.long,sample.elv),RcEst,consts.SPhiInf,mconsts,'no');

switch sf
    case 'St'
        mu = mu_out;
    case 'LSD'
        mu = mu_out;
end

mt = 70000;
clipindex = find(tv <= mt, 1, 'last' );
tv = tv(1:clipindex);

%First calculate N for all times within defined range in tv. Limit like in
%get age using tv2. clipindex can be 50000 since only for 14C. Then, will
%need to interpolate concentrations between tv time steps.
Lmu = 1500;

P_LS = P14s(1:clipindex); % trim production arrays to agree with age array

dcf = exp(-tv.*consts.l14);
dpfs = exp(-tv.*sample.E.*sample.rho./consts.Lsp); % spallation depth dependence
dpfm = exp(-tv.*sample.E.*sample.rho./Lmu); % muon depth dependence approximation

% sat14 = (P14s+mu)./((sample.rho.*sample.E./160)+consts.l14); %This will increment with each elevation step. And probably need to pick an age that gaurantees saturation (50 kyr).
Nsat(i) = interp1(tv,cumtrapz(tv,(P_LS.*dcf.*dpfs + mu.*dcf.*dpfm)),50000)
Nsat_low(i) = interp1(tv,cumtrapz(tv,(P_LS.*(1-0.056).*dcf.*dpfs + mu.*dcf.*dpfm)),50000);
Nsat_high(i) = interp1(tv,cumtrapz(tv,(P_LS.*(1+0.056).*dcf.*dpfs + mu.*dcf.*dpfm)),50000);
%Calculate concentrations for a given set of exposure lengths
time = 1e3.*[2:2:24];

    for j = 1:length(time)%This might be able to come out of the for loop
        N_time(j,i) = interp1(tv,cumtrapz(tv,(P_LS.*dcf.*dpfs + mu.*dcf.*dpfm)),time(j));
    %     N14(:,i) = (P14s+mu)./((sample.rho.*sample.E./160)+consts.l14).*(1-exp(-(((sample.rho.*sample.E./160)+consts.l14).*time(i))));
    end

end %End elevation loop here. Do each alt one at a time and log.
tiledlayout(1,2)
% figure('PaperSize',[3.5 5],'PaperUnits','inches');
nexttile
p1=patch([Nsat_low fliplr(Nsat_high)], [alts-ice_surface fliplr(alts-ice_surface)], [.8 .8 .8]);%shades saturation errorwindow (YOU'LL NEED TO MANUALLY LAY THE BOX OVER THIS PATCH IF YOU ZOOM IN)
hold on;
plot(Nsat,alts-ice_surface,'k',Nsat_low,alts-ice_surface,'k',Nsat_high,alts-ice_surface,'k'); %this plots up the saturation and its errorbounds
box on;
hold on;
% if exist('rep_unc') == 1
%     h = errorbar_x(N14_sample,all_elv-ice_surface,rep_unc.*N14_sample,rep_unc.*N14_sample,'bo');
% else
    h = errorbar_x(N14_sample,all_elv-ice_surface,delN14_sample,delN14_sample,'bo');%plots samples' data
    labelpoints(all_N14,all_elv-ice_surface,all_sample_name,'SE')%adds datapoint labels
% end

%Plot time curves
hold on;
for i = 1:size(N_time,1)
    plot(N_time(i,:),alts-ice_surface,'k--'); %this loop plots up the 2ka-saturation dashed lines
    hold on;
end
xlabel('C-14 Concentration (10^{5} atoms g^{-1})');
grid on;

if max(N14_sample) > max(alts-ice_surface)
    xmax = max(N14_sample);
else
    xmax = max(sat14);
end

if nargin == 1
    ymin = min(alts);
    ymax = max(alts);
else
    ymin = 0.9*min(alts);
    ymax = 1.1*max(alts);
end

axis([1e5 6.5e5 ymin ymax]);

if ice_surface > 0
    ylabel('Elevation Above Modern Glacier (m)');
else
    ylabel('Elevation (m asl)');
end

% foo = [alts' Nsat' Nsat_low' Nsat_high' N_time];
% save foo.txt foo -ASCII;
% 
% figure(2);
% %Make plot of Age vs % of saturation concentration. Just do for one
% %elevation for now.
% 
% pct_sat = N_time./repmat(Nsat14',1,size(N_time,2));
% 
% plot(time,pct_sat);
%  xlabel('Age (ka)'); ylabel('Frac. Saturation');

%makes an elevation-age plot.  you'll need to take the headings off your
%results .txt and add an elevation column.  you'll need to add xaxis lines and axis titles yourself in the editor 
% clear all
file = input('Please enter the name of the results data file','s');
FID = fopen(file);
data = textscan(FID,'%s %n %n %n %s %n %n %n %n %n %n %n %n %s %s');
fclose(FID);
dstring='';
all_sample_name = data{1}; all_t = data{2}; all_int = data{3};
all_ext = data{4};
all_elv = data{8};
nexttile
x=[all_t];
y=[all_elv];
sz=42;
s=scatter(x,y,sz,'s','r','filled');
s.LineWidth=0.6;
s.MarkerEdgeColor='k';
labelpoints(x,y,all_sample_name,'SE');
xlabel('Age (ka)');
ylabel('Elevation (m a.s.l.)');
hold on;
yneg=[0 0 0 0 0 0 0 0];
ypos=[yneg];
errorbar(x,y,yneg,ypos,all_int,all_ext,'r','LineStyle','none');
hold on;
chex=[13500 14700 14700 13500];
chey=[1800 1800 1950 1950];
p2=patch(chex,chey,'c','FaceAlpha',.5,'LineStyle','none');%shades mwp1a