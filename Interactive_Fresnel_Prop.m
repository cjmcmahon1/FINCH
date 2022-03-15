% More Fresnel Propagation Tests
% Same code as Short_Distance_Fresnel but with tests in response to
% comments
%{
TODO:
-Plot FT of the reconstructed hologram, look at line spread for different
propagation distances
-Check that the -1mm plot properly has both images displayed when fresnel
propagating
-Make the intensity axis consistent for different images for easier
comparison
-add images to README for reference
%}
addpath('./MATLAB_functions/'); %include helper functions

% Parameters; units mm
PARAMS = struct;
PARAMS.Lx = 250e-3;      %x side length of input image
PARAMS.Ly = 250e-3;      %y side length of input image
PARAMS.lambda = 490e-6; %wavelength
PARAMS.Mx = 2048;        %x samples
PARAMS.My = 2048;        %y samples
PARAMS.NA = 0.1;        %numerical aperture

%Generate fields by Fresnel propagating constant amplitude,
%circular aperture fields two different distances z1 & z2. 
%using propagate(z, parameters)
%the Brooker papers have z1~-10mm, z2~10mm
z1 = -1; %mm
z2 = 1; %mm
z_prop = -1; %mm
p1 = propagate_init(z1, PARAMS);
p2 = propagate_init(z2, PARAMS);
%create phase shifted interference patterns
shifted1 = shifted_hologram(p1, p2, 0 * pi / 3);
shifted2 = shifted_hologram(p1, p2, 2 * pi / 3);
shifted3 = shifted_hologram(p1, p2, 4 * pi / 3);
%generate the complex-valued hologram
%This will internally generate and use the phase shifted holograms above
hol = complex_hologram(p1, p2, 3);
%fresnel propagate the complex hologram backwards/forwards
hol_plane = fresnel_prop(hol.intensity, z_prop, PARAMS);
hol_struct = struct('intensity', hol_plane, 'x', hol.x, 'y', hol.y);
%The backwards propagation should look like two images, z1
%propagated 0mm, and z2 propagated -2mm
% back_comp_plane1 = propagate_init(0, PARAMS);
% back_comp_plane2 = propagate_init(-2, PARAMS);
% back_comparison = struct('field', ...
%     back_comp_plane1.field + back_comp_plane2.field, ...
%     'x', back_comp_plane1.x, 'y', back_comp_plane1.y);
%Forward propagation should look like z1=0, z2=2mm
% fwd_comp_plane1 = propagate_init(0, PARAMS);
% fwd_comp_plane2 = propagate_init(2, PARAMS);
% fwd_comparison = struct('intensity', ...
%     abs(fwd_comp_plane1.field^2) + abs(fwd_comp_plane2.field^2), ...
%     'x', fwd_comp_plane1.x, 'y', fwd_comp_plane1.y);

%plot comparision of the two intensities with the reconstructed hologram at
%values z = -1mm and z = +1mm
hfig = figure('Name', 'Ground Truth Field vs Propagated Hologram');
ax = axes('Parent',hfig,'position',[0.13 0.39  0.77 0.54]);
pos = get(hfig,'position');
set(hfig,'position',pos.*[0.25 0.25 2.2 1.9]); %make plot window wider
b = uicontrol('Parent',hfig,'Style','slider','Position',[375,220,475,24],...
              'value', z_prop, 'min', -1.5, 'max',1.5);
bgcolor = hfig.Color;
bl1 = uicontrol('Parent',hfig,'Style','text','Position',[368,190,30,20],...
                'String','-1.5','BackgroundColor',bgcolor);
bl2 = uicontrol('Parent',hfig,'Style','text','Position',[829,190,30,20],...
                'String','1.5','BackgroundColor',bgcolor);
bl3 = uicontrol('Parent', hfig,'Style','text','Position',[530,180,150,30],...
                'String','Propagation Distance z (mm)','BackgroundColor',bgcolor);
sys = gen_hol_im(hol, z_prop, PARAMS);
callback_fxn = imagesc(sys.x, sys.y, abs(sys.intensity));
axis('square');
colormap('gray');
%b.Callback = @(es,ed) updateSystem(callback_fxn,gen_hol_im(hol, es.Value, PARAMS)); 
b.Callback = {@myUpdateCB, callback_fxn, hol, PARAMS};

    function myUpdateCB(hObject, EvenData, s, hol, bench_params)
    updateSystem(s, gen_hol_im(hol, hObject.Value, bench_params));
end
function hol_struct = gen_hol_im(hol, z, bench_params)
    hol_plane = fresnel_prop(hol.intensity, z, bench_params);
    hol_struct = struct('intensity', hol_plane, 'x', hol.x, 'y', hol.y);
end