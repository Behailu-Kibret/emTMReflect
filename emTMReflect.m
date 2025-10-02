function out = emTMReflect(varargin)
%EMTMREFLECT  TM (parallel) plane-wave reflection/transmission at a planar interface.
%
% Medium 1: z <= 0   (incident + reflected)
% Medium 2: z >= 0   (transmitted or evanescent)
% Instantaneous fields use Re{· e^{+jωt}}.
%
% Usage (all inputs are optional name–value pairs):
%   emTMReflect;                                  % quick demo (defaults)
%   emTMReflect('Eps1',1,'Eps2',4,'ThetaInc',30); % change media & angle
%   out = emTMReflect('ComputeOnly',true);        % compute Fresnel only
%
% Key name–value inputs:
%   'Eps1',1, 'Mu1',1, 'Eps2',4, 'Mu2',1, 'Freq',1e9, 'ThetaInc',30, 'E0',1
%   'AnimateMagnitude',true, 'SlowMo',0.5, 'NPeriods',12
%   'GridSize',[500 500], 'DomainLambda',[4 4]   % extents in λ1
%   'RotateZRight',true                          % +z → right, +x ↑
%   'ShowRays',true, 'ShowAngleArcs',true, 'ShowAngleValues',true
%   'ShowMediumLabels',true, 'Colormap','turbo', 'FreezeCLim',true
%   'FigSizePx',[650 600], 'AxesSizePx',500, 'Renderer','opengl'
%   'SaveGIF',false,'GIFFile','tm_wave.gif'
%   'SaveVideo',false,'VideoFile','tm_wave.avi','VideoProfile','Motion JPEG AVI','VideoFPS',30
%   'MaxFrames',300, 'UpdateEvery',4, 'Verbose',true
%   'ComputeOnly',false, 'ReturnFields',false
%
% Output (struct):
%   Fresnel:  Gamma_b, T_b, R, Tpow, Regime ('prop'|'crit'|'evan')
%   Angles:   ThetaBrew, ThetaCrit, ThetaInc, ThetaTrans (if real)
%   Waves:    beta1, beta2, kx, kz1, kx2, kz2
%   Grid:     x, z, X, Z
%   Handles:  fig, ax (if plotted)
%   Fields:   Ex0, Ez0 (if 'ReturnFields',true)

%% ---------------- Parse inputs ----------------
ip = inputParser; ip.FunctionName = mfilename;
addParameter(ip,'Eps1',1,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'Mu1',1,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'Eps2',4,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'Mu2',1,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'Freq',1e9,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'ThetaInc',30,@(x)isnumeric(x)&&isscalar(x));
addParameter(ip,'E0',1,@(x)isnumeric(x)&&isscalar(x));

addParameter(ip,'AnimateMagnitude',true,@islogical);
addParameter(ip,'SlowMo',0.5,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'NPeriods',12,@(x)isnumeric(x)&&isscalar(x)&&x>0);

addParameter(ip,'GridSize',[500 500],@(v)isnumeric(v)&&numel(v)==2);
addParameter(ip,'DomainLambda',[4 4],@(v)isnumeric(v)&&numel(v)==2);
addParameter(ip,'RotateZRight',true,@islogical);

addParameter(ip,'ShowRays',true,@islogical);
addParameter(ip,'ShowAngleArcs',true,@islogical);
addParameter(ip,'ShowAngleValues',true,@islogical);
addParameter(ip,'ShowMediumLabels',true,@islogical);
addParameter(ip,'Colormap','turbo',@(s)ischar(s)||isstring(s));
addParameter(ip,'FreezeCLim',true,@islogical);

addParameter(ip,'FigSizePx',[650 600],@(v)isnumeric(v)&&numel(v)==2);
addParameter(ip,'AxesSizePx',500,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'Renderer','opengl',@(s)ischar(s)||isstring(s));

addParameter(ip,'SaveGIF',false,@islogical);
addParameter(ip,'GIFFile','tm_wave.gif',@(s)ischar(s)||isstring(s));
addParameter(ip,'SaveVideo',false,@islogical);
addParameter(ip,'VideoFile','tm_wave.avi',@(s)ischar(s)||isstring(s));
addParameter(ip,'VideoProfile','Motion JPEG AVI',@(s)ischar(s)||isstring(s));
addParameter(ip,'VideoFPS',30,@(x)isnumeric(x)&&isscalar(x)&&x>0);

addParameter(ip,'MaxFrames',300,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(ip,'UpdateEvery',4,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(ip,'Verbose',true,@islogical);

addParameter(ip,'ComputeOnly',false,@islogical);
addParameter(ip,'ReturnFields',false,@islogical);

parse(ip,varargin{:});
S = ip.Results;

%% ---------------- Constants ----------------
c0   = 299792458;
eta0 = 376.730313668;
w    = 2*pi*S.Freq;
k0   = w/c0;

%% ---------------- Material wave params ----------------
eta1  = eta0*sqrt(S.Mu1/S.Eps1);
eta2  = eta0*sqrt(S.Mu2/S.Eps2);
beta1 = k0*sqrt(S.Eps1*S.Mu1);
beta2 = k0*sqrt(S.Eps2*S.Mu2);

ThetaBrew = atand(sqrt(S.Eps2/S.Eps1));
ThetaCrit = NaN;
if S.Eps1 > S.Eps2
    ThetaCrit = asind(sqrt(S.Eps2/S.Eps1));
end

%% ---------------- Snell & branch selection (no clamp; correct evanescence) ----------------
si = sind(S.ThetaInc); 
ci = cosd(S.ThetaInc);

% Snell (allow st>1 beyond critical so cosθt becomes imaginary)
st  = (beta1/beta2)*si;          % sin(theta_t)
ct2 = 1 - st.^2;                 % cos^2(theta_t)
tol = 1e-12;

if ct2 > tol
    ct = sqrt(ct2);              % propagating in medium 2
    Regime = "prop";
elseif abs(ct2) <= tol
    ct = 0;                      % critical (grazing)
    Regime = "crit";
else
    % Evanescent: choose branch so kz2 = beta2*ct has NEGATIVE imaginary part,
    % giving exp(-alpha z) decay for z>0 with exp(-j·) spatial phasors.
    ct = -1i*sqrt(-ct2);         
    Regime = "evan";
end

%% ---------------- Fresnel TM (parallel) ----------------
Gamma_b = ( -eta1*ci + eta2*ct ) ./ ( eta1*ci + eta2*ct );
T_b     = ( 2*eta2*ci ) ./ ( eta1*ci + eta2*ct );

%% ---------------- Grid ----------------
E0 = S.E0;
lambda1 = 2*pi/beta1;
Nx = S.GridSize(1); Nz = S.GridSize(2);
Lx = S.DomainLambda(1)*lambda1;  Lz = S.DomainLambda(2)*lambda1;
x  = linspace(-Lx, Lx, Nx);
z  = linspace(-Lz, Lz, Nz);
[X,Z] = meshgrid(x, z);

%% ---------------- Phases & z-dependent Γ(z) ----------------
kx  = beta1*si;     kz1 = beta1*ci;
kx2 = beta2*st;     kz2 = beta2*ct;
phi1 = kx.*X + kz1.*Z;                          % β1(x sinθi + z cosθi)
phi2 = kx2.*X + kz2.*Z;                         % β2(x sinθt + z cosθt)
Gamma_z = Gamma_b .* exp(1i*2*kz1.*Z);          % standing-wave factor in medium 1

%% ---------------- Fields (phasors) ----------------
% Medium 1 (z<=0): total = incident + reflected
Ex1 =  E0*ci .* exp(-1i*phi1) .* (1 + Gamma_z);
Ez1 = -E0*si .* exp(-1i*phi1) .* (1 - Gamma_z);

% Medium 2 (z>=0): transmitted (or evanescent)
Ex2 =  (T_b*E0).*ct .* exp(-1i*phi2);
Ez2 = -(T_b*E0).*st .* exp(-1i*phi2);

% Assemble by region
m1 = (Z<=0); m2 = (Z>=0);
Ex0 = Ex1.*m1 + Ex2.*m2;
Ez0 = Ez1.*m1 + Ez2.*m2;

%% ---------------- Power sanity check ----------------
R = abs(Gamma_b).^2;
if Regime=="prop"
    Tpow = (eta1/eta2) * (ct/ci) * abs(T_b).^2;
else
    Tpow = 0;
end
if S.Verbose
    fprintf('R = %.4f,  T = %.4f,  R+T = %.4f   (regime: %s)\n', R, Tpow, R+Tpow, Regime);
end

%% ---------------- Populate outputs ----------------
out.Gamma_b = Gamma_b; out.T_b = T_b; out.R = R; out.Tpow = Tpow;
out.Regime = char(Regime);
out.ThetaBrew = ThetaBrew; out.ThetaCrit = ThetaCrit;
out.ThetaInc  = S.ThetaInc;
out.ThetaTrans = (Regime=="prop")*asind(st);
out.beta1=beta1; out.beta2=beta2; out.kx=kx; out.kz1=kz1; out.kx2=kx2; out.kz2=kz2;
out.x=x; out.z=z; out.X=X; out.Z=Z;

if S.ReturnFields
    out.Ex0 = Ex0; out.Ez0 = Ez0;
end

if S.ComputeOnly
    return;
end

%% ---------------- Plot setup (ROTATED: +z right, +x up) ----------------
frame0_inst = hypot(real(Ex0), real(Ez0));
maxAmpPhasor = max(hypot(abs(Ex0(:)), abs(Ez0(:))));

% Figure & axes (fixed pixels so video frames are consistent)
hFig = figure('Color','w','Units','pixels');
hFig.Position(1:2) = [100 100];
hFig.Position(3:4) = [S.FigSizePx];
ax = axes('Parent',hFig,'Units','pixels','ActivePositionProperty','position', ...
          'Position',[70 70 S.AxesSizePx S.AxesSizePx]);

% Plot orientation
if S.RotateZRight
    hImg = imagesc(ax, z*1e3, x*1e3, frame0_inst.');  % z → horizontal, x → vertical
    xlabel(ax,'z (mm)'); ylabel(ax,'x (mm)');
else
    hImg = imagesc(ax, x*1e3, z*1e3, frame0_inst);    % standard x→h, z→v
    xlabel(ax,'x (mm)'); ylabel(ax,'z (mm)');
end
axis(ax,'image'); set(ax,'YDir','normal'); colormap(ax, S.Colormap);
hCb = colorbar(ax);
title(ax, sprintf('TM: f=%.2f GHz, \\theta_i=%d^\\circ, t=0', S.Freq/1e9, S.ThetaInc));
if S.FreezeCLim, clim(ax,[0, maxAmpPhasor]); end

% Park colorbar to the right without shrinking axes:
drawnow;
if ~isempty(hCb) && isgraphics(hCb)
    set(hCb,'Units','pixels');
    cbw = 20; gap = 86;
    cbpos = get(ax,'Position');
    set(hCb,'Position',[cbpos(1)+cbpos(3)+gap, cbpos(2), cbw, cbpos(4)]);
end

hold(ax,'on');
% Axes lines through origin (in the chosen view)
% Axes lines through the origin (valid for both orientations)
xline(ax,0,'k-','LineWidth',1);  % vertical: z=0 (RotateZRight=true) or x=0 (false)
yline(ax,0,'k-','LineWidth',1);  % horizontal: x=0 (RotateZRight=true) or z=0 (false)


% Rays (mapped to plotting coordinates)
ui_phys = [sind(S.ThetaInc),  cosd(S.ThetaInc)];   % (x,z) incident (+z)
ur_phys = [sind(S.ThetaInc), -cosd(S.ThetaInc)];   % (x,z) reflected (-z)
hasTrans = (Regime=="prop") || (Regime=="crit");
if hasTrans
    %theta_t = asind(st);
    % Clamp ONLY for display to avoid asin domain warnings from tiny roundoff.
    theta_t = asind( min(1, max(-1, real(st))) );
    ut_phys = [sind(theta_t), cosd(theta_t)];
end
if S.RotateZRight
    % (z,x) on axes
    ui = [ui_phys(2), ui_phys(1)];
    ur = [ur_phys(2), ur_phys(1)];
    if hasTrans, ut = [ut_phys(2), ut_phys(1)]; end
else
    % (x,z) on axes
    ui = ui_phys; ur = ur_phys; if hasTrans, ut = ut_phys; end
end

% Set ray lengths & angle-arc radius from current limits
xl = xlim(ax); yl = ylim(ax);
L  = 0.22*min(diff(xl), diff(yl));
r  = 0.12*min(diff(xl), diff(yl));

if S.ShowRays
    quiver(ax, -ui(1)*L, -ui(2)*L, ui(1)*L, ui(2)*L, 0, 'LineWidth',1.5,'MaxHeadSize',0.8,'Color','k'); % incident
    quiver(ax, 0, 0, ur(1)*L, ur(2)*L, 0, 'LineWidth',1.5,'MaxHeadSize',0.8,'Color','k');               % reflected
    if hasTrans
        quiver(ax, 0, 0, ut(1)*L, ut(2)*L, 0, 'LineWidth',1.5,'MaxHeadSize',0.8,'Color','k');           % transmitted
    end
end

% Angle arcs & labels
a = deg2rad(S.ThetaInc);
if S.ShowAngleArcs
    if S.RotateZRight
        % parameterize circle in (z,x): z=r*cosψ, x=r*sinψ; ψ=0 → +z (right), ψ=π → -z (left)
        % INCIDENT arc (left, below): ψ in [π, π+a]
        psi_i = linspace(pi, pi + a, 60);
        plot(ax, r*cos(psi_i), r*sin(psi_i), 'k', 'LineWidth',1.2);
        text(ax, r*cos(pi+a/2), r*sin(pi+a/2), '\theta_i', 'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold','BackgroundColor','w');

        % REFLECTION arc (left, above): ψ in [π-a, π]
        psi_r = linspace(pi - a, pi, 60);
        plot(ax, r*cos(psi_r), r*sin(psi_r), 'k', 'LineWidth',1.2);
        text(ax, r*cos(pi-a/2), r*sin(pi-a/2), '\theta_r', 'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold','BackgroundColor','w');

        % TRANSMISSION arc (right, above): ψ in [0, a]
        if hasTrans
            at = deg2rad(theta_t);
            psi_t = linspace(0, at, 60);
            plot(ax, r*cos(psi_t), r*sin(psi_t), 'k', 'LineWidth',1.2);
            text(ax, r*cos(at/2), r*sin(at/2), '\theta_t', 'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold','BackgroundColor','w');
        end
    else
        % Standard (x,z) orientation; arcs mirrored accordingly
        % INCIDENT (below, right of -z axis)
        psi_i = linspace(pi, pi + a, 60);
        plot(ax, r*sin(psi_i), r*cos(psi_i), 'k', 'LineWidth', 1.2);
        text(ax, r*sin(pi + a/2), r*cos(pi + a/2), '\theta_i', ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'FontWeight','bold','BackgroundColor','w');

        % REFLECTION (above, right of -z axis)
        psi_r = linspace(pi - a, pi, 60);
        plot(ax, r*sin(psi_r), r*cos(psi_r), 'k', 'LineWidth', 1.2);
        text(ax, r*sin(pi - a/2), r*cos(pi - a/2), '\theta_r', ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'FontWeight','bold','BackgroundColor','w');

        % TRANSMISSION (above, left of +z axis)
        if hasTrans
            at = deg2rad(theta_t);
            psi_t = linspace(0, at, 60);   % z>0
            plot(ax, r*sin(psi_t), r*cos(psi_t), 'k', 'LineWidth', 1.2);
            text(ax, r*sin(at/2), r*cos(at/2), '\theta_t', ...
                'HorizontalAlignment','center','VerticalAlignment','middle', ...
                'FontWeight','bold','BackgroundColor','w');
        end
    end
end

% Angle values in boxes at bottom centers
if S.ShowAngleValues
    zr = xlim(ax); xr = ylim(ax);
    if S.RotateZRight
        z_mid_left  = mean([zr(1), 0]);
        z_mid_right = mean([0,    zr(2)]);
        x_bot = xr(1) + 0.04*diff(xr);
        txt_left  = sprintf('\\theta_i = %.1f^\\circ   |   \\theta_r = %.1f^\\circ', S.ThetaInc, S.ThetaInc);
        if hasTrans
            txt_right = sprintf('\\theta_t = %.1f^\\circ', theta_t);
        else
            txt_right = '\theta_t: evanescent';
        end
        text(ax, z_mid_left,  x_bot, txt_left,  'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold','BackgroundColor','w','EdgeColor','k','Margin',2,'Interpreter','tex');
        text(ax, z_mid_right, x_bot, txt_right, 'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold','BackgroundColor','w','EdgeColor','k','Margin',2,'Interpreter','tex');
    else
        % similar placement for (x,z) orientation (optional)
    end
end

% Medium labels at top centers
if S.ShowMediumLabels
    zl = xlim(ax); xl_ = ylim(ax);
    if S.RotateZRight
        z_mid_left  = mean([zl(1), 0]);
        z_mid_right = mean([0,    zl(2)]);
        x_top = xl_(2) - 0.04*diff(xl_);
        text(ax, z_mid_left,  x_top, sprintf('\\epsilon_1 = %.3g,  \\mu_1 = %.3g', S.Eps1, S.Mu1), 'HorizontalAlignment','center','VerticalAlignment','top','FontWeight','bold','BackgroundColor','w');
        text(ax, z_mid_right, x_top, sprintf('\\epsilon_2 = %.3g,  \\mu_2 = %.3g', S.Eps2, S.Mu2), 'HorizontalAlignment','center','VerticalAlignment','top','FontWeight','bold','BackgroundColor','w');
    else
        % similar placement for (x,z) orientation (optional)
    end
end
hold(ax,'off');

out.fig = hFig; out.ax = ax; out.img = hImg;

%% ---------------- Animation & Recording ----------------
T   = 1/S.Freq; 
dt  = T/40; 
steps = round(S.NPeriods*T/dt);

% Decide how many frames to write (cap to keep files small)
framesToWrite = min(S.MaxFrames, steps+1);
captureEvery  = max(1, floor((steps+1)/framesToWrite));

% Prepare writers
if S.SaveVideo
    v = VideoWriter(S.VideoFile, S.VideoProfile);
    v.FrameRate = S.VideoFPS;
    v.Quality   = 95;
    open(v);
end
if S.SaveGIF && exist(S.GIFFile,'file')
    delete(S.GIFFile);
end
isFirstGifFrame = true;
gifDelay = max(0.01, 1/S.VideoFPS);

% Loop
for n = 0:steps
    tEff = S.SlowMo * n*dt;
    tf   = exp(1i*w*tEff);                        % instantaneous: Re{· e^{+jωt}}

    if S.AnimateMagnitude
        frame = hypot(real(Ex0.*tf), real(Ez0.*tf));
    else
        frame = real(Ex0.*tf);
    end

    % Update UI every UpdateEvery steps
    if mod(n, S.UpdateEvery) == 0
        if S.RotateZRight, set(hImg,'CData', frame.'); else, set(hImg,'CData', frame); end
        title(ax, sprintf('TM: \\theta_i=%d^\\circ, t=%.2f ns', S.ThetaInc, tEff*1e9));
        drawnow;
    end

    % Capture every captureEvery steps
    if mod(n, captureEvery) == 0
        fr = getframe(ax);   % fixed-size frames (AxesSizePx × AxesSizePx)

        if S.SaveVideo
            writeVideo(v, fr);
        end

        if S.SaveGIF
            [imind, cmap] = rgb2ind(fr.cdata, 256, 'nodither');
            if isFirstGifFrame
                imwrite(imind, cmap, S.GIFFile, 'gif', 'LoopCount', inf, 'DelayTime', gifDelay);
                isFirstGifFrame = false;
            else
                imwrite(imind, cmap, S.GIFFile, 'gif', 'WriteMode', 'append', 'DelayTime', gifDelay);
            end
        end
    end
end

if S.SaveVideo
    close(v);
    if S.Verbose, fprintf('Saved video: %s (%.0f fps)\n', S.VideoFile, S.VideoFPS); end
end
if S.SaveGIF && S.Verbose
    fprintf('Saved GIF:   %s (~%.0f fps)\n', S.GIFFile, 1/gifDelay);
end
end
