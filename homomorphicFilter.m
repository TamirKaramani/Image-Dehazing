function homomorphicFilter(I,sigma,gammaH,gammaL)
% Authors: Shay Karadi  & Tamir Karamani
%{
The following function applies Homomorphic filter to hazed
images in order to achieve a better dehazed image.
%}
%tic is used to calculate the run time of the function
tic
Original = im2double(imread(I));
[r,c,k] = size(Original);
%Gauss HPF
M = 2 * r + 1;
N = 2 * c + 1;
[X, Y] = meshgrid(1:N,1:M);
Arg = (X - ceil(N / 2)) .^ 2 + (Y - ceil(M / 2)) .^ 2;
GHPF = (gammaH - gammaL) * (1 - exp(- Arg./(2 * sigma .^ 2))) + gammaL;
%Logarithm transform
LOG = log(1 + Original);
%Fourier transform
FFT = fftshift(fft2(LOG,M,N));
%apply HPF on fourier transform
GHPF_FFT = repmat(GHPF,1,1,k) .* FFT;
%Inverse fourier transform
IFFT = real(ifft2(ifftshift(GHPF_FFT),M,N));
IFFT = IFFT(1:r,1:c,:);
%Exponent transform
RGB = exp(IFFT) - 1;
figure('Name','Homomorphic Filter Dehazing','NumberTitle','off');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5, 0.1, 0.5, 0.9]);
subplot(2,4,1),imshow(Original,[]);title("Haze Image");
subplot(2,4,2),imshow(abs(LOG),[]);title("LOG");
subplot(2,4,3),imshow(log(abs(FFT)+1),[]);title("FFT");
subplot(2,4,4),imshow(log(GHPF+1),[]);title("Gaussian High Pass Filter");
subplot(2,4,5),imshow(log(abs(FFT)+1),[]);title("GHPF*FFT");
subplot(2,4,6),imshow(IFFT,[]);title("IFFT");
subplot(2,4,7),imshow(RGB,[]);title("EXP");
figure('Name','Homomorphic Filter Dehazing','NumberTitle','off');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.55, 0.5, 0.45]);
sgtitle("Original hazed & Dehazed"+newline+"Algorithm time: "+round(1000*toc)+" ms");
subplot(1,2,1),imshow(Original,[]);
subplot(1,2,2),imshow(RGB,[]);
end