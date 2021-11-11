%% SWIFT (semantic wavelet-induced frequency-tagging) version 1.0
%
% This function generates image sequences by modulating the semantic
% information of an input image at the temporal frequency f0 using cyclic
% wavelet scrambling.
%
% sequence = swift(f0,frame_rate,in_img)
%
% Inputs:
% f0            - Temporal frequency (in Hertz) of the cyclic wavelet scrambling
% frame_rate    - Speed (in frames per second) at which the sequence will be presented
% in_img        - Input image matrix in 2D (grayscale) or 3D (RGB)
%
% Output:
% sequence      - 3D matrix corresponding to the cyclic wavelet scrambled
%                 sequence. The first dimension correspond to the individual
%                 frames in the scrambling cycle and the two last dimensions
%                 correspond to the size in pixels of each frame.
%
% In order to conserve the tag frequency (f0) in the experiment, make sure
% that the sequences are presented at the frame rate stated above.
% You need have the wavelet toolbox instaled in order to run this
% function.
%
% Authors: Roger Koenig-Robert and Rufin VanRullen
%
% Please refer to the following article:
% Koenig-Robert, R., VanRullen, R., SWIFT: A novel method to track the neural correlates of recognition, NeuroImage
% (2013), http://dx.doi.org/10.1016/j.neuroimage.2013.04.116
 

function sequence = swift(f0,frame_rate,in_img)
%% Input parameters:
freqdist = [0.5 0.5]; % distribution of harmonics in each scrambling cycle (default: [0.2 0.2 0.2 0.2 0.2])
levels = 6; % number of decomposition levels (default: 6)
filter = 'dmey'; % wavelet family used in the decomposition (default: dmey)
normPix = 1; % local luminance spectra normalization (default: 1)
%%
nbframes = frame_rate/f0; % number of frames per cycle of wavelet scrambling
if rem(nbframes,1)~=0
    nbframes = round(nbframes);
    f0 = frame_rate/nbframes;
    display(['f0 changed to ' num2str(f0) ' Hz to match an entire number of frames per cycle']);
end
harms = 1:length(freqdist);
%%
clear('image','C','S', 'hor','vert','diag','h1','v1','d1','h2','v2','d2','hc','vc','dc','sequence','circle')
%% transform to grayscale if RGB
if size(in_img,3)==3
    in_img = rgb2gray(in_img);
end
%% wavelet decomposition
[C,S] = wavedec2(double(in_img),levels,filter);
%% extract coefficients
oldindex  = 1;
index = prod(S(1,:));
app = reshape(C(oldindex:index),S(1,:));
for level = 1:levels
    %%extract horizontal
    oldindex=index+1;
    index = index+prod(S(level+1,:));
    hor(level).coeffs=reshape(C(oldindex:index),S(level+1,:));
    %%extract vertical
    oldindex=index+1;
    index = index+prod(S(level+1,:));
    vert(level).coeffs=reshape(C(oldindex:index),S(level+1,:));
    %%extract diagonal
    oldindex=index+1;
    index = index+prod(S(level+1,:));
    diag(level).coeffs=reshape(C(oldindex:index),S(level+1,:));
end
%% energy extraction in each dimension
for level = 1:levels
    energy = sqrt(diag(level).coeffs.^2 + vert(level).coeffs.^2 + hor(level).coeffs.^2); %total energy per location & level
    newvec1 = 2*rand(size(energy,1),size(energy,2),3)-1; %defines random vector1
    newvec2 = 2*rand(size(energy,1),size(energy,2),3)-1; %defines random vector2
    newenergy1 = squeeze(sqrt(sum(newvec1.^2,3))); %calculates energy of vector1
    newenergy2 = squeeze(sqrt(sum(newvec2.^2,3))); %calculates energy of vector2
    h1(level).coeffs = energy.*(squeeze(newvec1(:,:,1)))./newenergy1; %vector1 energy normalization horizontal
    v1(level).coeffs = energy.*(squeeze(newvec1(:,:,2)))./newenergy1; %vector1 energy normalization vertical
    d1(level).coeffs = energy.*(squeeze(newvec1(:,:,3)))./newenergy1; %vector1 energy normalization diagonal
    h2(level).coeffs = energy.*(squeeze(newvec2(:,:,1)))./newenergy2; %vector2 energy normalization horizontal
    v2(level).coeffs = energy.*(squeeze(newvec2(:,:,2)))./newenergy2; %vector2 energy normalization vertical
    d2(level).coeffs = energy.*(squeeze(newvec2(:,:,3)))./newenergy2; %vector2 energy normalization diagonal
end
%% definition of the circular path per level and location
for level = 1:levels
    p1=[];
    p2=[];
    p3=[];
    c =[];
    pprime = [];
    xprime = [];
    %% defining final points of the original vector (p1) and the new random vectors (p2 and p3) which define a sphere
    p1(:,:,1) = hor(level).coeffs;
    p1(:,:,2) = vert(level).coeffs;
    p1(:,:,3) = diag(level).coeffs;
    p2(:,:,1) = h1(level).coeffs;
    p2(:,:,2) = v1(level).coeffs;
    p2(:,:,3) = d1(level).coeffs;
    p3(:,:,1) = h2(level).coeffs;
    p3(:,:,2) = v2(level).coeffs;
    p3(:,:,3) = d2(level).coeffs;
    %% resolving the the circular path defined by the 3 vectors
    t = p2-p1; u = p3-p1; v = p3-p2;
    w = cross(t,u,3);
    dottt = squeeze(dot(t,t,3));
    dotuv = squeeze(dot(u,v,3));
    dotuu = squeeze(dot(u,u,3));
    dottv = squeeze(dot(t,v,3));
    dotww = squeeze(dot(w,w,3));
    dotvv = squeeze(dot(v,v,3));
    c(:,:,1) = squeeze(p1(:,:,1)) + (dottt.*dotuv.* (squeeze(u(:,:,1))) - dotuu.*dottv.* (squeeze(t(:,:,1)))) ./ (2*dotww);
    c(:,:,2) = squeeze(p1(:,:,2)) + (dottt.*dotuv.* (squeeze(u(:,:,2))) - dotuu.*dottv.* (squeeze(t(:,:,2)))) ./ (2*dotww);
    c(:,:,3) = squeeze(p1(:,:,3)) + (dottt.*dotuv.* (squeeze(u(:,:,3))) - dotuu.*dottv.* (squeeze(t(:,:,3)))) ./ (2*dotww);
    r = (1/2).*sqrt(dottt).*sqrt(dotuu).*sqrt(dotvv)./sqrt(dotww);
    x1 = p1-c;
    x2 = p2-c;
    x3 = p3-c;
    alphasurbeta = -dot(x1,x3,3)./dot(x1,x2,3);
    pprime(:,:,1)=squeeze(c(:,:,1))+squeeze(x3(:,:,1))+alphasurbeta.*(squeeze(x2(:,:,1)));
    pprime(:,:,2)=squeeze(c(:,:,2))+squeeze(x3(:,:,2))+alphasurbeta.*(squeeze(x2(:,:,2)));
    pprime(:,:,3)=squeeze(c(:,:,3))+squeeze(x3(:,:,3))+alphasurbeta.*(squeeze(x2(:,:,3)));
    xprime(:,:,1)=(r.*(squeeze((pprime(:,:,1)-c(:,:,1))))./(sqrt(sum((pprime-c).^2,3))));
    xprime(:,:,2)=(r.*(squeeze((pprime(:,:,2)-c(:,:,2))))./(sqrt(sum((pprime-c).^2,3))));
    xprime(:,:,3)=(r.*(squeeze((pprime(:,:,3)-c(:,:,3))))./(sqrt(sum((pprime-c).^2,3))));
    freq = randsample(harms,numel(hor(level).coeffs),true,freqdist); % applies different harmonic modulations at each coefficient
    freq = reshape(freq, size(hor(level).coeffs));
    for frame = 1:nbframes %discrete steps going around the circle
        coswave = cos(freq.*2*pi*frame/nbframes);
        coswave = repmat(coswave,[1 1 3]);
        sinwave = sin(freq.*2*pi*frame/nbframes);
        sinwave = repmat(sinwave,[1 1 3]);
        circle.level(frame,level).coord = c + coswave .* x1 + sinwave .* xprime;
    end %end frame
end %end level
%% coef reconstruction
for frame = 1:nbframes/4
    for level = 1:levels % nb of levels to reconstruct
        hc(level).coeffs = circle.level(frame,level).coord(:,:,1);
        vc(level).coeffs = circle.level(frame,level).coord(:,:,2);
        dc(level).coeffs = circle.level(frame,level).coord(:,:,3);
    end
    %% reconstruct the wavelet multi-scale pyramid
    C = [];
    C = [C app(:)'];
    for level = 1:levels
        C = [C hc(level).coeffs(:)'];
        C = [C vc(level).coeffs(:)'];
        C = [C dc(level).coeffs(:)'];
    end;
    %% image domian sequence
    sequence(frame,:,:) = waverec2(C,S,filter);
end
clear('circle', 'C')
%%
sequence(sequence>255) = 255;
sequence(sequence<0) = 0;
%% Local luminance modulation normalization
if normPix == 1
    %% cutoff computing
    fpixel = linspace(0,frame_rate,nbframes);
    lastHarm = f0*harms(end);
    penultHarm = f0*harms(end-1);
    [~, penultHarmCoef] = min(abs(fpixel-penultHarm));
    [~, lastHarmCoef] = min(abs(fpixel-lastHarm));
    harmDelta = lastHarmCoef - penultHarmCoef;
    cutoff = lastHarmCoef + round(harmDelta/2);
    %% luminance spectra normalization
    Y = fft(sequence,[],1);
    Ynorm = Y;
    Ynorm(2:end,:,:) = 0.8*squeeze(mean(mean(mean(abs(Y(2:cutoff,:,:)),3),2),1))*Y(2:end,:,:)./(abs(Y(2:end,:,:)));
    Ynorm(cutoff:end-cutoff+2,:,:) = 0;
    sequence = abs(ifft(Ynorm,[],1));
end
%% Frame mean luminance normalization
refMean = 85;%mean(mean(sequence(size(sequence,1),:,:)));
for frame = 1:nbframes/4
    frameMean = squeeze(mean(mean(sequence(frame,:,:),2),3));
    sequence(frame,:,:) = sequence(frame,:,:).*(refMean/frameMean);
end
