addpath(genpath('/MATLAB Drive/MUST'))

% We want a 2.7-MHz 64-element cardiac phased array.
param = getparam('P4-2v');

% Add transmit apodization & set Delay 
% The left ventricle will be insonified with 
% seven diverging waves 60 degrees wide 
% and tilted at -20 to +20 degrees.

tilt = deg2rad(linspace(-20,20,7)); % tilt angles in rad
txdel = cell(7,1); % this cell will contain the transmit delays

% to calculate the transmit delays for the 7 diverging waves.
for k = 1:7
    txdel{k} = txdelay(param,tilt(k),deg2rad(60));
end

% Plot the transmit delays to obtain a 60 degrees wide circular waves steered at -20 degrees
% -------------------------------------------------------
figure(1);
stem(txdel{1}*1e6)
xlabel('Element number')
ylabel('Delays (\mus)')
title('TX delays for a 60{\circ}-wide -20{\circ}-tilted wave')
axis tight square

% Simulate an acoustic pressure field with PFIELD
[xi,zi] = impolgrid([100 100],15e-2,deg2rad(120),param); % 100x100 polar grid using IMPOLGRID

% RMS Pressure Field
option.WaitBar = false;
P = pfield(xi,0*xi,zi,txdel{1},param,option);

% Plot the Field 
% -----------------------------
figure(2);
pcolor(xi*1e2,zi*1e2,20*log10(P/max(P(:))))
shading interp
xlabel('x (cm)')
ylabel('z (cm)')
title('RMS pressure field for a 60{\circ} wide -20{\circ} tilted wave')
axis equal ij tight
caxis([-20 0]) % dynamic range = [-20,0] dB
cb = colorbar;
cb.YTickLabel{end} = '0 dB';
colormap(hot)


I = rgb2gray(imread('heart.png'));
% Pseudorandom distribution of scatterers (depth is 15 cm)
[x,y,z,RC] = genscat([NaN 15e-2],1540/param.fc,I);

figure(3);
scatter(x*1e2,z*1e2,2,abs(RC).^.25,'filled')
colormap([1-hot;hot])
axis equal ij tight
set(gca,'XColor','none','box','off')
title('Scatterers for a cardiac 5-chamber view')
ylabel('[cm]')

% Simulate the seven series of RF signals with SIMUS.
% The RF signals will be sampled at 4 $\times$ center frequency.

RF = cell(7,1); % this cell will contain the RF series
param.fs = 4*param.fc; % sampling frequency in Hz

option.WaitBar = false; % remove the wait bar of SIMUS
h = waitbar(0,'');
for k = 1:7
    waitbar(k/7,h,['SIMUS: RF series #' int2str(k) ' of 7'])
    RF{k} = simus(x,y,z,RC,txdel{k},param,option);
end
close(h)

% Plot RF of 1st series - 32th element
% rf = RF{1}(:,32);
% t = (0:numel(rf)-1)/param.fs*1e6; % time (ms)
% plot(t,rf)
% set(gca,'YColor','none','box','off')
% xlabel('time (\mus)')
% title('RF signal of the 32^{th} element (1^{st} series, tilt = -20{\circ})')
% axis tight


IQ = cell(7,1);  % this cell will contain the I/Q series

for k = 1:7
    IQ{k} = rf2iq(RF{k},param.fs,param.fc);
end

% iq = IQ{1}(:,32);
% plot(t,real(iq),t,imag(iq))
% set(gca,'YColor','none','box','off')
% xlabel('time (\mus)')
% title('I/Q signal of the 32^{th} element (1^{st} series, tilt = -20{\circ})')
% legend({'in-phase','quadrature'})
% axis tight


[xi,zi] = impolgrid([256 128],15e-2,deg2rad(80),param);
bIQ = zeros(256,128,7);  % this array will contain the 7 I/Q images

h = waitbar(0,'');
for k = 1:7
    waitbar(k/7,h,['DAS: I/Q series #' int2str(k) ' of 7'])
    bIQ(:,:,k) = das(IQ{k},xi,zi,txdel{k},param);
end
close(h)

bIQ = tgc(bIQ);


I = bmode(bIQ(:,:,1),50); % log-compressed image
figure(4);
pcolor(xi*1e2,zi*1e2,I)
shading interp, colormap gray
title('DW-based echo image with a tilt angle of -20{\circ}')

axis equal ij
set(gca,'XColor','none','box','off')
c = colorbar;
c.YTick = [0 255];
c.YTickLabel = {'-50 dB','0 dB'};
ylabel('[cm]')


cIQ = sum(bIQ,3); % this is the compound beamformed I/Q
I = bmode(cIQ,50); % log-compressed image
figure(5);
pcolor(xi*1e2,zi*1e2,I)
shading interp, colormap gray
title('Compound DW-based cardiac echo image')

axis equal ij
set(gca,'XColor','none','box','off')
c = colorbar;
c.YTick = [0 255];
c.YTickLabel = {'-50 dB','0 dB'};
ylabel('[cm]')
