function bilateraldehaze(I)
% Authors: Shay Karadi  & Tamir Karamani
%{
The following function applies guided joint bilateral filter to hazed
images in order to achieve a better dehazed image.
%}
%tic is used to calculate the run time of the function
tic
Original = im2double(imread(I));
%converts the original RGB image to grayscale
gray = rgb2gray(Original);
[row,col,~] = size(Original);
W = zeros(row,col);
%calculates the minimum channel image W ,Eq 2.1
for i = 1:row
    for j = 1:col
        W(i,j) = min([Original(i,j,1),Original(i,j,2),Original(i,j,3)]);
    end
end
%chosen parameters taken from the article 
p = 0.7;
omega = 0.95;    
Neighborhood = floor(2 * max(row,col) / 50 + 1);
%applies median filtering to get the atmospheric scattering light V(x),Eq 2.2-2.4
B = medfilt2(W,[Neighborhood Neighborhood]);
C = B - medfilt2(abs(W - B),[Neighborhood Neighborhood]);
V = max(min(p * C,W),0);
%applies bilateral filter to get the referance image R(x)
R = imbilatfilt(W,0.01,0.03 * max(row,col));
%applies guided joint bilateral filter to get the corrected atmospheric veil VR(x)
VR = imguidedfilter(V,R,'DegreeOfSmoothing',0.01, 'neighborhood',[3 3]);    
%calculates the transmission t(x)
vec = maxk(unique(W),ceil(0.002 * row * col));
MaxIm = zeros(1,length(vec));
for i = 1:length(vec)
  MaxIm(i) = max(max(gray(W == vec(i))));
end
A = max(MaxIm);
t = 1 - omega * VR / A;
t0 = 1.15*min(t(:));
%calculates the recovered image J(x)
J = zeros(row,col,3);
for i=1:row
    for j=1:col
        for k=1:3
            J(i,j,k) = (Original(i,j,k) - A) ./ max(t(i,j),t0) + A;
        end
    end
end
figure('Name','Dehazing using guided joint bilateral filter','NumberTitle','off');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5, 0.1, 0.5, 0.9]);
subplot(3,3,1),imshow(Original,[]),title("Haze");
subplot(3,3,2),imshow(W,[]);title("W-Minimum Channel");
subplot(3,3,3),imshow(V,[]);title("V-Median Filtering");
subplot(3,3,4),imshow(R,[]);title("R-Reference");
subplot(3,3,5),imshow(VR,[]);title("VR-Corrected Atmospheric Veil");
subplot(3,3,6),imshow(t,[]);title("t-Transmission");
subplot(3,3,7),imshow(J,[]);title("J-Recoverd");
figure('Name','Dehazing using guided joint bilateral filter','NumberTitle','off');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.55, 0.5, 0.45]);
sgtitle("Original hazed & Dehazed"+newline+"Algorithm time: "+round(1000*toc)+" ms");
subplot(1,2,1),imshow(Original,[]);
subplot(1,2,2),imshow(J,[]);
end
