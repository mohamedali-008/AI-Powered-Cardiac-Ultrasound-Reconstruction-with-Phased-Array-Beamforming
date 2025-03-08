function [RP,param,SPECT,IDX] = pfield3(varargin)

%PFIELD3   3-D RMS acoustic pressure field of a planar 2-D array
%   RP = PFIELD3(X,Y,Z,DELAYS,PARAM) returns the three-dimensional
%   radiation pattern of a planar 2-D array whose elements are excited at
%   different time delays (given by the vector DELAYS). The radiation
%   pattern RP is given in terms of the root-mean-square (RMS) of acoustic
%   pressure. The characteristics of the array and transmission must be
%   given in the structure PARAM. The x- and y-coordinates of the elements
%   must be in PARAM.elements (see below for details). The radiation
%   pattern is calculated at the points specified by (X,Y,Z).
%
%   Use PFIELD for uniform linear or convex arrays.
%
%   >--- Try it: enter "pfield3" in the command window for an example ---<
%
%   Units: X,Y,Z must be in m; DELAYS must be in s.
%
%   DELAYS can be a matrix. This syntax can be used to simulate MPT
%   (multi-plane transmit) sequences, for example. In this case, each ROW
%   represents a delay series. For example, to create a 3-MPT sequence with
%   a 1024-element matrix array, the DELAYS matrix must have 3 rows and
%   1024 columns (size = [3 1024]).
%
%   Note: Use TXDELAY3 to create standard delays (focus point, focus line,
%         plane waves, diverging waves) with a matrix array.
%
%   PFIELD3 is called by SIMUS3 to simulate ultrasound RF radio-frequency
%   signals generated by a planar 2-D array.
%
%   ---
%   NOTE #1: X-, Y-, and Z-axes
%   Conventional axes are used:
%   The X-axis is PARALLEL to the transducer and points from the first
%   (leftmost) element to the last (rightmost) element (X = 0 at the CENTER
%   of the transducer). The Z-axis is PERPENDICULAR to the transducer and
%   points downward (Z = 0 at the level of the transducer, Z increases as
%   depth increases). The Y-axis is such that the coordinates are
%   right-handed.
%   ---
%   NOTE #2: Simplified method: Directivity
%   By default, the calculation is made faster by assuming that the
%   directivity of the elements is dependent only on the central frequency.
%   This simplification very little affects the pressure field in most
%   situations (except near the array). To turn off this option, use
%   OPTIONS.FullFrequencyDirectivity = true.
%   (see ADVANCED OPTIONS below).
%   ---
%   NOTE #3: Present version of PFIELD3
%   The present version of PFIELD3 only considers planar 2-D array whose
%   z-coordinates of the elements are 0. All the elements are rectangular
%   and have equal width and height.
%   ---
%
%   PARAM is a structure that contains the following fields:
%   -------------------------------------------------------
%       *** TRANSDUCER PROPERTIES ***
%   1)  PARAM.fc: central frequency (in Hz, REQUIRED)
%   2)  PARAM.elements: x- and y-coordinates of the element centers
%            (in m, REQUIRED). It MUST be a two-row matrix, with the 1st
%            and 2nd rows containing the x and y coordinates, respectively. 
%   3)  PARAM.width: element width, in the x-direction (in m, REQUIRED)
%   4)  PARAM.height: element height, in the y-direction (in m, REQUIRED)
%   5)  PARAM.bandwidth: pulse-echo 6dB fractional bandwidth (in %)
%            The default is 75%.
%   6)  PARAM.baffle: property of the baffle:
%            'soft' (default), 'rigid', or a scalar > 0.
%            See "Note on BAFFLE properties" below for details
%
%       *** MEDIUM PARAMETERS ***
%   7)  PARAM.c: longitudinal velocity (in m/s, default = 1540 m/s)
%   8)  PARAM.attenuation: attenuation coefficient (dB/cm/MHz, default: 0)
%            Notes: A linear frequency-dependence is assumed.
%                   A typical value for soft tissues is ~0.5 dB/cm/MHz.
%
%       *** TRANSMIT PARAMETERS ***
%   9)  PARAM.TXapodization: transmit apodization (default: no apodization)
%   10) PARAM.TXnow: number of wavelengths of the TX pulse (default: 1)
%   11) PARAM.TXfreqsweep: frequency sweep for a linear chirp (default: [])
%                          To be used to simulate a linear TX down-chirp.
%
%   Other syntaxes:
%   --------------
%   i}  [RP,PARAM] = PFIELD3(...) also returns the complete list of
%       parameters including the default values.
%   ii} [...] = PFIELD3 without any input argument runs an example of a
%       focused acoustic field generated by a 3.5 MHz matrix array.
%
%
%   Note on CHIRP signals:
%   ---------------------
%   Linear chirps are characterized by PARAM.TXnow, PARAM.fc and
%   PARAM.TXfreqsweep. The transmitted pulse has a duration of
%   approximately T (= PARAM.TXnow/PARAM.fc), with the amplitude and phase
%   defined over the time interval -T/2 to +T/2. The total frequency sweep
%   is DeltaF (= PARAM.TXfreqsweep): the frequencies changes linearly from
%   (PARAM.fc + DeltaF/2) to (PARAM.fc - DeltaF/2) in the defined time
%   interval.
%
%
%   Note on BAFFLE property:
%   -----------------------
%   In PFIELD3, it is assumed by default that the array elements are
%   embedded in an infinite SOFT baffle. To modify the property of the
%   baffle, modify the field PARAM.baffle:
%       1) 'rigid'
%       2) 'soft' (this is the default)
%       3) a nonnegative scalar Alpha,
%          with Alpha = (medium impedance)/(baffle impedance)
%          Note: Alpha = 0 => 'rigid'; Alpha >> 1 => 'soft'
%
%   The baffle property affects the obliquity factor included in the
%   directivity of the elements. This obliquity factor is not 1 if the
%   baffle is not rigid. A general case (see case #3 below) can be chosen
%   by specifying an impedance ratio. For details, refer to the
%   corresponding papers.
%   1) For a rigid baffle => obliquity factor = 1.
%   2) For a soft baffle => obliquity factor = cos(Theta).
%      Selfridge et al. Appl Phys Lett 37(1), 35-36 (1980)
%      "A theory for the radiation pattern of a narrow-strip acoustic
%      transducer." <a
%      href="matlab:web('http://scitation.aip.org/content/aip/journal/apl/37/1/10.1063/1.91692')">Paper here</a>
%   3) General baffle => obliquity factor = cos(Theta)/(cos(Theta)+Alpha)
%      with Alpha = (medium impedance)/(baffle impedance).
%      Pesqu� et al. IEEE Ultrasonics Symposium, (1984)
%      "Effect of the planar baffle impedance in acoustic radiation of a
%      phased array element theory and experimentation." <a
%      href="matlab:web('http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=1535402')">Paper here</a>
%      Example: For a baffle of impedance 2.8 MRayl (epoxy) adjacent to
%               soft tissues of impedance 1.6 MRayls, alpha = 0.57.
%
%
%   ADVANCED OPTIONS:
%   ----------------
%       %-- FREQUENCY SAMPLES --%
%   1)  Only frequency components of the transmitted signal in the range
%       [0,2fc] with significant amplitude are considered. The default
%       relative amplitude is -60 dB in PFIELD. You can change this value
%       by using the following:
%           [...] = PFIELD(...,OPTIONS),
%       where OPTIONS.dBThresh is the threshold in dB (default = -60).
%   ---
%       %-- FULL-FREQUENCY DIRECTIVITY --%   
%   2)  By default, the directivity of the elements depends on the center
%       frequency only. This makes the algorithm faster. To make the
%       directivities fully frequency-dependent, use: 
%           [...] = PFIELD(...,OPTIONS),
%       with OPTIONS.FullFrequencyDirectivity = true (default = false).
%   ---
%       %-- ELEMENT SPLITTING --%   
%   3)  Each transducer element of the array is split into small rectangles.
%       The width and height and of these small rectangles must be small
%       enough to ensure that the far-field model is accurate. By default,
%       the elements are split into M-by-N rectangles, with M and N being
%       defined by:
%           M = ceil(element_width/smallest_wavelength);
%           N = ceil(element_height/smallest_wavelength);
%       To modify the number MN of subelements by splitting, you may adjust
%       OPTIONS.ElementSplitting, which must contain two elements. For
%       example, OPTIONS.ElementSplitting = [1 3].
%   ---
%       %-- WAIT BAR --%   
%   4)  If OPTIONS.WaitBar is true, a wait bar appears (only if the number
%       of frequency samples >10). Default is true.
%   ---
%
%
%   Notes regarding the model & REFERENCES:
%   --------------------------------------
%   1) PFIELD3 works for planar 2-D arrays. It considers arrays that have
%      identical rectangular elements on the z=0 plane. Each element is
%      split into small rectangles (if required). As the sub-elements are
%      small enough, the three-dimensional radiation patterns are derived
%      by using Fraunhofer (far-field) equations.
%   2) The paper that describes the first 2-D version of PFIELD is:
%      SHAHRIARI S, GARCIA D. Meshfree simulations of ultrasound vector
%      flow imaging using smoothed particle hydrodynamics. Phys Med Biol,
%      2018;63:205011. <a
%      href="matlab:web('https://www.biomecardio.com/publis/physmedbio18.pdf')">PDF here</a>
%   3) The papers that describe the theory and validation of the 2-D + 3-D
%      versions of PFIELD and SIMUS are:
%      i)  GARCIA D. SIMUS: an open-source simulator for medical ultrasound
%          imaging. Part I: theory & examples. Comput Methods Programs
%          Biomed, 2022;218:106726. <a
%          href="matlab:web('https://www.biomecardio.com/publis/cmpb22.pdf')">PDF here</a>
%      ii) CIGIER A, VARRAY F, GARCIA D. SIMUS: an open-source simulator
%          for medical ultrasound imaging. Part II:comparison with four
%          simulators. Comput Methods Programs Biomed, 2022;220:106774.
%          <a href="matlab:web('https://www.biomecardio.com/publis/cmpb22a.pdf')">PDF here</a>
%   4) Use the fonction <a href="matlab:cite")>CITE</a> to guide you in citations.
%   5) There is yet no publication for PFIELD3 (it is planned for 2023-24).
%
%
%   Example #1:
%   ----------
%   %-- Generate a focused pressure field with a matrix array
%   % 3-MHz matrix array with 32x32 elements
%   param = [];
%   param.fc = 3e6;
%   param.bandwidth = 70;
%   param.width = 250e-6;
%   param.height = 250e-6;
%   % position of the elements (pitch = 300 microns)
%   pitch = 300e-6;
%   [xe,ye] = meshgrid(((1:32)-16.5)*pitch);
%   param.elements = [xe(:).'; ye(:).'];
%   % Focus position
%   x0 = 0; y0 = -2e-3; z0 = 30e-3;
%   % Transmit time delays: 
%   dels = txdelay3(x0,y0,z0,param);
%   % 3-D grid
%   n = 32;
%   [xi,yi,zi] = meshgrid(linspace(-5e-3,5e-3,n),linspace(-5e-3,5e-3,n),...
%       linspace(0,6e-2,4*n));
%   % RMS pressure field
%   RP = pfield3(xi,yi,zi,dels,param);
%   % Display the pressure field
%   slice(xi*1e2,yi*1e2,zi*1e2,20*log10(RP/max(RP(:))),...
%       x0*1e2,y0*1e2,z0*1e2)
%   shading flat
%   colormap(hot), caxis([-20 0])
%   set(gca,'zdir','reverse'), axis equal
%   alpha color % some transparency
%   c = colorbar; c.YTickLabel{end} = '0 dB';
%   zlabel('[cm]')
%
%   Example #2:
%   ----------
%   %-- Generate a pressure field focused on a line
%   % 3-MHz matrix array with 32x32 elements
%   param = [];
%   param.fc = 3e6;
%   param.bandwidth = 70;
%   param.width = 250e-6;
%   param.height = 250e-6;
%   % position of the elements (pitch = 300 microns)
%   pitch = 300e-6;
%   [xe,ye] = meshgrid(((1:32)-16.5)*pitch);
%   param.elements = [xe(:).'; ye(:).'];
%   % Oblique focus-line @ z = 2.5cm
%   x0 = [-1e-2 1e-2]; y0 = [-1e-2 1e-2]; z0 = [2.5e-2 2.5e-2];
%   % Transmit time delays: 
%   dels = txdelay3(x0,y0,z0,param);
%   % 3-D grid
%   n = 32;
%   [xi,yi,zi] = meshgrid(linspace(-5e-3,5e-3,n),linspace(-5e-3,5e-3,n),...
%       linspace(0,6e-2,4*n));
%   % RMS pressure field
%   RP = pfield3(xi,yi,zi,dels,param);
%   % Display the elements
%   figure, plot3(xe*1e2,ye*1e2,0*xe,'b.')
%   % Display the pressure field
%   contourslice(xi*1e2,yi*1e2,zi*1e2,RP,[],[],.5:.5:6,15)
%   set(gca,'zdir','reverse'), axis equal
%   colormap(hot)
%   zlabel('[cm]')
%   view(-35,20), box on
%
%
%   This function is part of <a
%   href="matlab:web('https://www.biomecardio.com/MUST')">MUST</a> (Matlab UltraSound Toolbox).
%   MUST (c) 2020 Damien Garcia, LGPL-3.0-or-later
%
%   See also PFIELD, TXDELAY3, SIMUS3, GETPARAM, GETPULSE, CITE.
%
%   -- Damien Garcia -- 2022/10, last update 2023/08/11
%   website: <a
%   href="matlab:web('https://www.biomecardio.com')">www.BiomeCardio.com</a>


MyFolder = 'C:\Users\garcia\Documents'; % PLZ do not modify,
% unless you are getting tired of the advertising message...
% In this case, use any folder that you own.
% 
if ~exist(MyFolder,'dir')
    mlock % NOTE: use MUNLOCK if you need to clear PFIELD from memory
    persistent pfieldLastUse %#ok
end

if nargin==0
    if nargout>0
        [RP,param] = RunTheExample;
    else
        RunTheExample;
    end
    return
end
narginchk(5,6)

%-- Input variables: X,Y,Z,DELAYS,PARAM,OPTIONS
x = varargin{1};
y = varargin{2};
z = varargin{3};
assert(isequal(size(x),size(y),size(z)),...
    'X, Y, and Z must be of same size.')
delaysTX = varargin{4};
param = varargin{5};
if nargin==6
    options = varargin{6};
else
    options = [];
end

%-- Check the transmit delays
assert(isnumeric(delaysTX) && all(delaysTX(~isnan(delaysTX))>=0),...
    'DELAYS must be a nonnegative array.')
if isvector(delaysTX) % we need a row vector
    delaysTX = delaysTX(:).'; 
end
NumberOfElements = size(delaysTX,2);

%--
% Note: delaysTX can be a matrix. This option can be used for MPT
% (multi-plane transmit) for example. In this case, each row represents a
% delay series. For example, for a 4-MPT sequence with a 1024-element
% matrix array, delaysTX has 4 rows and 1024 columns, i.e. size(delaysTX) =
% [4 1024].
%--
delaysTX = delaysTX.';

%-- Check if PFIELD3 is called by SIMUS3
isSIMUS3 = false;
if isfield(options,'CallFun')
    isSIMUS3 = strcmpi(options.CallFun,'simus3');
end

%-- Check if SIMUS3 is running on a parallel pool
try
    onppool = isSIMUS3 && options.ParPool;
catch
    onppool = false;
end

%-- Advertising message & Statistics
if ~onppool && ~exist(MyFolder,'dir')
    if isempty(pfieldLastUse) || (now-pfieldLastUse)>1
        pfieldLastUse = now;
        
        % An ad appears once a day (once PFIELD is used) in a small message
        % box. Please do not remove it!
        AdMessage

        % A warning message appears if one of the functions is outdated
        try checkUpdate, catch, end
        
        % -- Statistics on daily use -- The number of uses is incremented
        % daily, and the information is sent in an FTP file for statistical
        % purpose. Don't worry, there is no spam! The MUSTstat file is
        % encrypted (p code) because it contains a password. Please do not
        % remove it!
        p = fileparts(mfilename('fullpath'));
        try MUSTstat(p), catch, end

    end
end


%---------------------------%
% Check the PARAM structure %
%---------------------------%

param = IgnoreCaseInFieldNames(param);

%-- 1) Center frequency (in Hz)
assert(isfield(param,'fc'),...
    'A center frequency value (PARAM.fc) is required.')
fc = param.fc; % central frequency (Hz)

%-- 2) Coordinates of the transducer elements (xe,ye)
% note: ze = 0 in this PFIELD3 version.
assert(isfield(param,'elements'),...
    ['PARAM.elements must contain the x- and y-locations ',...
    'of the transducer elements.'])
assert(size(param.elements,1)==2,...
    ['PARAM.elements must have two rows that contain the ',...
    'x (1st row) and y (2nd row) coordinates of the transducer elements.'])
xe = param.elements(1,:);
ye = param.elements(2,:);
assert(numel(xe)==NumberOfElements,...
    'The number of elements must match the number of transmit delays.')

%-- 3) Element width (in m)
assert(isfield(param,'width'),...
    'An element width (PARAM.width) is required.')
ElementWidth = param.width;
assert(isnumeric(ElementWidth) && isscalar(ElementWidth) &&...
    ElementWidth>0,'The element width must be positive.')

%-- 4) Element height (in m)
assert(isfield(param,'height'),...
    'An element height (PARAM.height) is required with PFIELD3 and SIMUS3.')
ElementHeight = param.height;
assert(isnumeric(ElementHeight) && isscalar(ElementHeight) &&...
    ElementHeight>0,'The element height must be positive.')

%-- 5) Fractional bandwidth at -6dB (in %)
if ~isfield(param,'bandwidth')
    param.bandwidth = 75;
end
assert(param.bandwidth>0 && param.bandwidth<200,...
    'The fractional bandwidth at -6 dB (PARAM.bandwidth, in %) must be in ]0,200[')

%-- 6) Baffle
%   An obliquity factor will be used if the baffle is not rigid
%   (default = SOFT baffle)
if ~isfield(param,'baffle')
    param.baffle = 'soft'; % default
end
if strcmpi(param.baffle,'rigid')
    NonRigidBaffle = false;
elseif strcmpi(param.baffle,'soft')
    NonRigidBaffle = true;
elseif isscalar(param.baffle)
    assert(param.baffle>0,...
        'The ''baffle'' field scalar must be positive')
    NonRigidBaffle = true;
else
    error('The ''baffle'' field must be ''rigid'',''soft'' or a positive scalar')
end

%-- 7) Longitudinal velocity (in m/s)
if ~isfield(param,'c')
    param.c = 1540; % default value
end
c = param.c; % speed of sound (m/s)

%-- 8) Attenuation coefficient (in dB/cm/MHz)
if ~isfield(param,'attenuation') % no attenuation, alpha_dB = 0
    param.attenuation = 0;
    alpha_dB = 0;
else
    alpha_dB = param.attenuation;
    assert(isscalar(alpha_dB) && isnumeric(alpha_dB) && alpha_dB>=0,...
        'PARAM.attenuation must be a nonnegative scalar')
end

%-- 9) Transmit apodization (no unit)
if ~isfield(param,'TXapodization')
    param.TXapodization = ones(1,NumberOfElements);
else
    assert(isvector(param.TXapodization) && isnumeric(param.TXapodization),...
        'PARAM.TXapodization must be a vector')
    assert(numel(param.TXapodization)==NumberOfElements,...
        'PARAM.TXapodization must be of length = (number of elements)')
end
% apodization is 0 where TX delays are NaN:
idx = isnan(delaysTX);
param.TXapodization(any(idx,2)) = 0;
delaysTX(idx) = 0;

%-- 10) TX pulse: Number of wavelengths
if ~isfield(param,'TXnow')
    param.TXnow = 1;
end
NoW = param.TXnow;
assert(isscalar(NoW) && isnumeric(NoW) && NoW>0,...
    'PARAM.TXnow must be a positive scalar.')

%-- 11) TX pulse: Frequency sweep for a linear chirp
if ~isfield(param,'TXfreqsweep') || isinf(NoW)
    param.TXfreqsweep = [];
end
FreqSweep = param.TXfreqsweep;
assert(isempty(FreqSweep) ||...
    (isscalar(FreqSweep) && isnumeric(FreqSweep) && FreqSweep>0),...
    'PARAM.TXfreqsweep must be empty (windowed sine) or a positive scalar (linear chirp).')

%----------------------------------%
% END of Check the PARAM structure %
%----------------------------------%



%-----------------------------%
% Check the OPTIONS structure %
%-----------------------------%

options = IgnoreCaseInFieldNames(options);

%-- 1) dB threshold
%     (in dB: faster computation if lower value, but less accurate)
if ~isfield(options,'dBThresh')
    options.dBThresh = -60; % default is -60dB in PFIELD3
end
assert(isscalar(options.dBThresh) && isnumeric(options.dBThresh) &&...
    options.dBThresh<=0,'OPTIONS.dBThresh must be a nonpositive scalar.')

%-- 2) Frequency-dependent directivity?
if isfield(options,'FullFrequencyDirectivity')
    isFFD = options.FullFrequencyDirectivity;
else
    isFFD = false; % default
    % By default, the directivity of the elements depends on the center
    % frequency only. This makes the algorithm faster. 
end
assert(isscalar(isFFD) && islogical(isFFD),...
    'OPTIONS.FullFrequencyDirectivity must be a logical scalar (true or false).')

%-- 3) Element splitting
%
% --- A short note about the algorithm:
% Far-field equations are used in PFIELD3. Each transducer element of the
% array is split into M-by-N small rectangles, so that these MN rectangles
% have size smaller than one wavelength by one wavelength. The far-field
% condition is acceptable for these small rectangles.
%---
if isfield(options,'ElementSplitting') && ~isempty(options.ElementSplitting)
    assert(numel(options.ElementSplitting)==2,...
        'OPTIONS.ElementSplitting must be a two-element vector.')
    M = options.ElementSplitting(1);
    N = options.ElementSplitting(2);
    assert(isscalar(M) & M==round(M) & M>0 &...
        isscalar(N) & N==round(N) & N>0,...
        'OPTIONS.ElementSplitting must contain two positive integers.')
else
    LambdaMin = c/(fc*(1+param.bandwidth/200));
    M = ceil(ElementWidth/LambdaMin);
    N = ceil(ElementHeight/LambdaMin);
end

%-- 4) Wait bar
if ~isfield(options,'WaitBar')
    options.WaitBar = true;
end
assert(isscalar(options.WaitBar) && islogical(options.WaitBar),...
    'OPTIONS.WaitBar must be a logical scalar (true or false).')

%-- Advanced (masked) options: Frequency step (scaling factor)
% The frequency step is determined automatically. It is tuned to avoid
% significant interferences due to unadapted discretization. The frequency
% step can also be adjusted by using a scaling factor. For a fast check,
% you may use a scaling factor>1. For a smoother result, you may use a
% scaling factor<1.
if ~isfield(options,'FrequencyStep')
    options.FrequencyStep = 1;
end
assert(isscalar(options.FrequencyStep) &&...
    isnumeric(options.FrequencyStep) && options.FrequencyStep>0,...
    'OPTIONS.FrequencyStep must be a positive scalar.')

%------------------------------------%
% END of Check the OPTIONS structure %
%------------------------------------%


%-----
% SIMUS3 first runs PFIELD3 with empty X,Y,Z to detect syntax errors.
if isSIMUS3 && isempty(x), RP = []; return, end
%-----


%------------------------------------%
% POINT LOCATIONS, DISTANCES & GRIDS %
%------------------------------------%

siz0 = size(x);
nx = numel(x);

%-- Coordinates of the points where pressure is needed
x = x(:); y = y(:); z = z(:);

%-- Cast x, y, and z to single class
x = cast(x,'single');
y = cast(y,'single');
z = cast(z,'single');

%-- Centroids of the sub-elements
%-- note: Each elements is split into M-by-N sub-elements.
% X-position (xi) and Y-position (yi) of the centroids of the sub-elements
% (relative to the centers of the transducer elements).
% The values in xi are in the range ]-ElementWidth/2 ElementWidth/2[.
% The values in yi are in the range ]-ElementHeight/2 ElementHeight/2[.
% (if M = 1 and N = 1, then xi = yi = 0).
SegWidth = ElementWidth/M;
xi = -ElementWidth/2 + SegWidth/2 + (0:M-1)*SegWidth;
SegHeight = ElementHeight/N;
yi = -ElementHeight/2 + SegHeight/2 + (0:N-1)*SegHeight;
[xi,yi] = meshgrid(xi,yi);
xi = reshape(xi,[1 1 M*N]);
yi = reshape(yi,[1 1 M*N]);

%-- Out-of-field points
% Null pressure will be assigned to out-of-field points.
isOUT = z<0;

%-- Variables that we need:
%
% Note: We work in an ISO spherical system for each sub-element
%       r = distance between the segment centroid and the point of interest
%       sinT = sine theta: theta is the polar angle.
%       cosT = cosine theta.
%       sinP = sine phi: phi is the azimuthal angle.
%       cosP = cosine phi.
%       They are of size [numel(x) NumberOfElements M*N].
%
dxi = x-xi-xe;
dyi = y-yi-ye;
d2 = dxi.^2+dyi.^2;
r = sqrt(d2+z.^2);
%---
epss = eps('single');
cosT = (z+epss)./(r+epss);
sinT = (sqrt(d2)+epss)./(r+epss);
cosP = (dxi+epss)./(sqrt(d2)+epss);
sinP = (dyi+epss)./(sqrt(d2)+epss);
clear dxi dyi d2
%---
% The term 1/r is present in the equations (problems if r is very small!):
% small r values are replaced by lambda/2
lambda = c/fc;
r(r<lambda/2) = lambda/2;

%-------------------------------------------%
% end of POINT LOCATIONS, DISTANCES & GRIDS %
%-------------------------------------------%


mysinc = @(x) sin(abs(x)+epss)./(abs(x)+epss); % cardinal sine
% [note: In MATLAB, sinc is sin(pi*x)/(pi*x)]

%-------------------%
% FREQUENCY SPECTRA %
%-------------------%

%-- FREQUENCY SPECTRUM of the transmitted pulse
if isempty(FreqSweep)
    % We want a windowed sine of width PARAM.TXnow
    T = NoW/fc; % temporal pulse width
    if isinf(T); T = 1e6; end
    wc = 2*pi*fc;
    pulseSpectrum = @(w) 1i*(mysinc(T*(w-wc)/2)-mysinc(T*(w+wc)/2));
else
    % We want a linear chirp of width PARAM.TXnow
    % (https://en.wikipedia.org/wiki/Chirp_spectrum#Linear_chirp)
    T = NoW/fc; % temporal pulse width
    if isinf(T); T = 1e6; end
    wc = 2*pi*fc;
    dw = 2*pi*FreqSweep;
    s2 = @(w) sqrt(pi*T/dw)*exp(-1i*(w-wc).^2*T/2/dw).*...
        (fresnelint((dw/2+w-wc)/sqrt(pi*dw/T)) + fresnelint((dw/2-w+wc)/sqrt(pi*dw/T)));
    pulseSpectrum = @(w) (1i*s2(w)-1i*s2(-w))/T;
end

%-- FREQUENCY RESPONSE of the ensemble PZT + probe
% We want a generalized normal window (6dB-bandwidth = PARAM.bandwidth)
% (https://en.wikipedia.org/wiki/Window_function#Generalized_normal_window)
wB = param.bandwidth*wc/100; % angular frequency bandwidth
p = log(126)/log(2*wc/wB); % p adjusts the shape
probeSpectrum = @(w) exp(-(abs(w-wc)/(wB/2/log(2)^(1/p))).^p);
% The frequency response is a pulse-echo (transmit + receive) response. A
% square root is thus required when calculating the pressure field:
probeSpectrum = @(w) sqrt(probeSpectrum(w));
% Note: The spectrum of the pulse (pulseSpectrum) will be then multiplied
% by the frequency-domain tapering window of the transducer (probeSpectrum)

%-- FREQUENCY STEP
if isSIMUS3 % PFIELD3 has been called by SIMUS3
    df = options.FrequencyStep;
else % We are in PFIELD3 only (i.e. not called by SIMUS3)
    % The frequency step df is chosen to avoid interferences due to
    % inadequate discretization.
    % -- df = frequency step (must be sufficiently small):
    % One has exp[-i(k r + w delay)] = exp[-2i pi(f r/c + f delay)] in the Eq.
    % One wants: the phase increment 2pi(df r/c + df delay) be < 2pi.
    % Therefore: df < 1/(r/c + delay).
    df = 1/(max(r(:)/c) + max(delaysTX(:)));
    df = options.FrequencyStep*df;
    % note: df is here an upper bound; it will be recalculated below
end

%-- FREQUENCY SAMPLES
Nf = 2*ceil(param.fc/df)+1; % number of frequency samples
f = linspace(0,2*param.fc,Nf); % frequency samples
df = f(2); % update the frequency step
%- we keep the significant components only by using options.dBThresh
S = abs(pulseSpectrum(2*pi*f).*probeSpectrum(2*pi*f));
GdB = 20*log10(S/max(S)); % gain in dB
IDX = GdB>options.dBThresh;
IDX(find(IDX,1):find(IDX,1,'last')) = true;
f = f(IDX);
nSampling = length(f);

%-- we need VECTORS
pulseSPECT = pulseSpectrum(2*pi*f); % pulse spectrum
probeSPECT = probeSpectrum(2*pi*f); % probe response

%--------------------------%
% end of FREQUENCY SPECTRA %
%--------------------------%



%-- Wait bar
options.WaitBar = options.WaitBar & (nSampling>10);
if options.WaitBar
    if isSIMUS3
        wbtitle = 'Let SIMUS3 do the work for you...';
        wbname = 'SIMUS3 / www.biomecardio.com';
    else
        wbtitle = 'Let PFIELD3 do the work for you...';
        wbname = 'PFIELD3 / www.biomecardio.com';
    end
    hwb = waitbar(0,wbtitle,'Name',wbname);
end

%-- Initialization
RP = 0; % RP = Radiation Pattern
if isSIMUS3
    %- For SIMUS3 only (we need the full spectrum of RX signals):
    SPECT = zeros([nSampling NumberOfElements],'like',single(1i));
elseif nargout==3
    SPECT = zeros([nSampling nx],'like',single(1i));
end

%-- Obliquity factor (baffle property)
%   An obliquity factor is required if the baffle is not rigid.
%   [Th = angle relative to the element normal axis]
if NonRigidBaffle
    if strcmpi(param.baffle,'soft')
        ObliFac = cosT;
    else % param.baffle is a scalar
        ObliFac = cosT./(cosT+param.baffle);
    end
else % 1 if rigid baffle
    ObliFac = ones(size(cosT));
end

%-- Note on Attenuation
% Reference: Diagnostic ultrasound imaging - inside out (T.L. Szabo)
%            Chapter 4: Attenuation
% Key reference: Acoustics for ultrasound imaging (Ben Cox, 2013)
%                Chapter 5: Acoustic attenuation and absorption
% We will use this attenuation-based wavenumber:
%   kwa = alpha_dB/8.69*f(k)/1e6*1e2; % P(z,f) = P0 exp(-alpha*f*z/8.69)
%   note: 20/log(10) ~ 8.69

%-- EXPONENTIAL arrays of size [numel(x) NumberOfElements M]
kw = 2*pi*f(1)/c; % wavenumber
kwa = alpha_dB/8.69*f(1)/1e6*1e2; % attenuation-based wavenumber
EXP = exp(-kwa*r + 1i*mod(kw*r,2*pi)); % faster than exp(-kwa*r+1i*kw*r)
%-- Exponential array for the increment wavenumber dk
dkw = 2*pi*df/c;
dkwa = alpha_dB/8.69*df/1e6*1e2;
EXPdf = exp((-dkwa + 1i*dkw)*r);

%-- We replace EXP by EXP.*ObliFac./r
EXP = EXP.*ObliFac./r;
clear ObliFac r

%-- TX apodization
APOD = param.TXapodization(:);

%-- Simplified directivity (if not dependent on frequency)
% In the "simplified directivity" version, the directivities of the
% sub-elements depend on the center frequency ONLY. It is thus not needed
% to calculate the directivity arrays (DIRx and DIRy) in the following
% for-loop. These directivities DIRx and DIRy are included in the variable
% EXP to reduce storage.
if ~isFFD
    kc = 2*pi*fc/c; % center wavenumber
    DIRx = mysinc(kc*SegWidth/2*cosP.*sinT); % x-directivity of each segment
    DIRy = mysinc(kc*SegHeight/2*sinP.*sinT); % y-directivity of each segment
    EXP = EXP.*DIRx.*DIRy;
    clear DIRx DIRy
end



%-----------------------------%
% SUMMATION OVER THE SPECTRUM %
%-----------------------------%

tstart = tic;
for k = 1:nSampling

    kw = 2*pi*f(k)/c; % wavenumber
    
    %-- Exponential array of size [numel(x) NumberOfElements MxN]
    % For all k, we need: EXP = exp((-kwa+1i*kw)*r)
    %                         with kw = 2*pi*f(k)/c;
    %                     and with kwa = alpha_dB/8.7*f(k)/1e6*1e2;
    % Since f(k) = f(1)+(k-1)df, we use the following recursive product:
    if k>1
        EXP = EXP.*EXPdf;
    end
        
    %-- Directivity (if frequency-dependent)
    if isFFD % isFFD = true -> frequency-dependent directivity
        DIRx = mysinc(kw*SegWidth/2*cosP.*sinT); % x-directivity
        DIRy = mysinc(kw*SegHeight/2*sinP.*sinT); % y-directivity
        DIR = DIRx.*DIRy;
    end
        
    %-- Radiation patterns of the single elements
    % They are the combination of the far-field patterns of the M small
    % segments that make up the single elements
    %--
    if isFFD % isFFD = true -> frequency-dependent directivity
        RPmono = mean(DIR.*EXP,3); % summation over the M*N small segments
    else % isFFD = false: the directivity depends on center frequency only
         % note: the directivity (DIR) has already been included in EXP
        if M*N>1           
            RPmono = mean(EXP,3); % summation over the M*N small segments
        else
            RPmono = EXP;
        end
    end
        
    %-- Transmit delays + Transmit apodization
    % use of SUM: summation over the number of delay series (e.g. MLT)
    DELAPOD = sum(exp(1i*kw*c*delaysTX),2).*APOD;
    
    %-- Summing the radiation patterns generating by all the elements
    RPk = RPmono*DELAPOD;
    %- include spectrum responses:
    RPk = pulseSPECT(k)*probeSPECT(k)*RPk;
    RPk(isOUT) = 0;
       
    %-- Output
    if isSIMUS3 % Receive: for SIMUS3 only (spectra of the RF signals)
        SPECT(k,:) = probeSPECT(k) *... % the array bandwidth is considered
            (RPk.*options.RC(:)).'*RPmono ... % pressure received by the elements
            ; % *f(k)^2/fc^2; % Rayleigh scattering (OPTIONAL)
        if any(param.RXdelay) % reception delays, if any
            SPECT(k,:) = SPECT(k,:).*exp(1i*kw*c*param.RXdelay);
        end
    else % using PFIELD3 alone
        RP = RP + abs(RPk).^2; % acoustic intensity
        if nargout==3
            SPECT(k,:) = RPk;
        end
    end
    
    %- update the wait bar
    if options.WaitBar && rem(k,10)==0
        tstep = toc(tstart);
        trem = tstep*(nSampling/k-1);
        waitbar(k/nSampling,hwb,...
            {['Elapsed: ',int2str(floor(tstep/60)) ,' min ',...
            int2str(floor(rem(tstep,60))),' s'],...
            ['Remaining: ',int2str(floor(trem/60)) ,' min ',...
            int2str(floor(rem(trem,60))),' s']})
    end
    
end

%------------------------------------%
% end of SUMMATION OVER THE SPECTRUM %
%------------------------------------%



% Close the wait bar
if options.WaitBar, close(hwb), end

% Correcting factor (including integration step, df)
if isinf(NoW)
    CorFac = 1;
else
    CorFac = df;
end
% if exist('SPECT','var'), SPECT = SPECT*CorFac; end
RP = RP*CorFac;

% RMS acoustic pressure (if we are in PFIELD3 only)
if ~isSIMUS3
    RP = reshape(sqrt(RP),siz0);
    if nargout>2
        SPECT = reshape(SPECT.',[siz0 nSampling]);
    end
end

end



function [RP,param] = RunTheExample

% Parameters of a 3.5 MHz matrix array
param.fc = 3.5e6;
param.c = 1540;
param.bandwidth = 70;
param.width = 250e-6;
param.height = 250e-6;

% Positions of the 1024 elements
spacing = 300e-6;
[xe,ye] = meshgrid(((1:32)-16.5)*spacing);
param.elements = [xe(:).'; ye(:).'];

% Transmit time delays: focus at (-2,0,30) mm
dels = txdelay3(-2e-3,0,30e-3,param);

% 3-D grid
n = 32;
[xi,yi,zi] = meshgrid(linspace(-5e-3,5e-3,n),linspace(-5e-3,5e-3,n),...
    linspace(0,6e-2,4*n));

% RMS pressure field
RP = pfield3(xi,yi,zi,dels,param);

% 3-D figure
figure
%--
% y-z slice @ x = -2 mm
% x-z slice @ y = 0 mm
% x-y slice @ z = 30 mm
slice(xi*1e3,yi*1e3,zi*1e3,RP,-2,0,30)
%--
set(gca,'zdir','reverse')
zlabel('[mm]')
shading flat
axis equal
colormap([1-hot;hot])
hold on
plot3(xe(:)*1e3,ye(:)*1e3,0*xe(:),'.')
title('A focused field by a 32-by-32 matrix array')

% Animate
for i = 1:720
    camorbit(.5,0,'data',[0 0 1])
    drawnow
end

end



function f = fresnelint(x)

% FRESNELINT Fresnel integral.
%
% J = FRESNELINT(X) returns the Fresnel integral J = C + 1i*S.
%
% We use the approximation introduced by Mielenz in
%       Klaus D. Mielenz, Computation of Fresnel Integrals. II
%       J. Res. Natl. Inst. Stand. Technol. 105, 589 (2000), pp 589-590
%

siz0 = size(x);
x = x(:);

issmall = abs(x)<=1.6;
c = zeros(size(x));
s = zeros(size(x));

% When |x| < 1.6, a Taylor series is used (see Mielenz's paper)
if any(issmall)
    n = 0:10;
    cn = [1 cumprod(-pi^2*(4*n+1)./(4*(2*n+1).*(2*n+2).*(4*n+5)))];
    sn = [1 cumprod(-pi^2*(4*n+3)./(4*(2*n+2).*(2*n+3).*(4*n+7)))]*pi/6;
    n = [n 11];
    c(issmall) = sum(cn.*x(issmall).^(4*n+1),2);
    s(issmall) = sum(sn.*x(issmall).^(4*n+3),2);    
end

% When |x| > 1.6, we use the following:
if any(~issmall)
    n = 0:11;
    fn = [0.318309844, 9.34626e-8, -0.09676631, 0.000606222, ...
        0.325539361, 0.325206461, -7.450551455, 32.20380908, ...
        -78.8035274, 118.5343352, -102.4339798, 39.06207702];
    gn = [0, 0.101321519, -4.07292e-5, -0.152068115, -0.046292605, ...
        1.622793598, -5.199186089, 7.477942354, -0.695291507, ...
        -15.10996796, 22.28401942, -10.89968491];
    fx = sum(fn.*x(~issmall).^(-2*n-1),2);
    gx = sum(gn.*x(~issmall).^(-2*n-1),2);    
    c(~issmall) = 0.5*sign(x(~issmall)) + ...
        fx.*sin(pi/2*x(~issmall).^2) - gx.*cos(pi/2*x(~issmall).^2);
    s(~issmall) = 0.5*sign(x(~issmall)) - ...
        fx.*cos(pi/2*x(~issmall).^2) - gx.*sin(pi/2*x(~issmall).^2);
end

f = reshape(c,siz0) + 1i*reshape(s,siz0);

end


function structArray = IgnoreCaseInFieldNames(structArray)

switch inputname(1)
    case 'param'
        fieldLIST = {'attenuation','baffle','bandwidth','c','elements','fc',...
            'fnumber','focus','fs','height','kerf','movie','Nelements',...
            'passive','pitch','radius','RXangle','RXdelay',...
            'TXapodization','TXfreqsweep','TXnow','t0','width'};
    case 'options'
        if isstruct(structArray)
            fieldLIST = {'dBThresh','ElementSplitting',...
                'FullFrequencyDirectivity','FrequencyStep','ParPool',...
                'WaitBar'};
        else
            return
        end
end

OldFieldNames = fieldnames(structArray);
tmp = lower(OldFieldNames);
assert(length(tmp)==length(unique(tmp)),...
    ['The structure ' upper(inputname(1)),...
    ' contains duplicate field names (when ignoring case).'])

[idx,loc] = ismember(lower(fieldLIST),tmp);
idx = find(idx); loc = loc(idx);
for k = 1:length(idx)
    tmp = eval(['structArray.' OldFieldNames{loc(k)}]); %#ok
    structArray = rmfield(structArray,OldFieldNames{loc(k)});
    eval(['structArray.' fieldLIST{idx(k)} ' = tmp;']) %#ok
end

end


function AdMessage
msgboxStruct.Interpreter = 'tex';
msgboxStruct.WindowStyle = 'modal';
try
[icondata,iconcmap] =...
    imread('https://www.biomecardio.com/images/MUSTicon.jpg');
catch
    icondata = ind2rgb(round(255*rescale(peaks(100))+1),...
        [1-hot(128); hot(128)]);
    iconcmap = [];
end
uiwait(msgbox(...
    {'\fontsize{9}MUST is an open-source toolbox for research purposes.',...
    'Please cite the corresponding articles.',...
    '',...
    'MUST \copyright 2020 Damien Garcia',...
    'LGPL-3.0-or-later','',...
    '\bf{garcia.damien@gmail.com}',...
    '\color[rgb]{0,0.5,0.5}\bullet \bf{www.biomecardio.com} \bullet'},...
    'MUST','custom',icondata,iconcmap,msgboxStruct));
end


function checkUpdate

% update according to the website
tmp = webread('https://www.biomecardio.com/MUST/index.html');
last_update = regexp(tmp,'\d\d\d\d/\d\d/\d\d','match','once');
last_update = datenum(last_update,'yyyy/mm/dd');

% update according to the .m files
dirname = fileparts(mfilename('fullpath'));
mlist = cellstr(ls([dirname '\*.m']));
for k = 1:length(mlist)
    tmpdate = regexp(fileread(mlist{k}),...
        '\d\d\d\d/\d\d/\d\d','match','once');
    if isempty(tmpdate), continue, end
    tmpdate = datenum(tmpdate);
    isupdated = tmpdate>=last_update;
    if isupdated, return, end
end

opts = struct('WindowStyle','modal','Interpreter','tex');
uiwait(warndlg({'\bf{}Your MUST toolbox is outdated.',...
    '\rm\color{red}\bullet Consider an upgrade \bullet'},...
    'Out-Of-Date',opts));

end