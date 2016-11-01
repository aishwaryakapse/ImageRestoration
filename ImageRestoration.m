clc;
close all;
clear all;

tic;

I=imread('lena.jpg');

subplot(2,3,1);
imshow(I);  %Original Image
title('Original Image');

Is=fftshift(I);
subplot(2,3,2);
imshow(Is); %Center Shifted Image 
title('Center Shifted Image');

Is=double(Is);

%Fourier Transform of Center Shifted Image
Ifft=fft2(Is);
subplot(2,3,3);
imshow(uint8(Ifft));    %Fourier Transform of Center Shifted Image
title('FFT of Center Shifted Image');

[p,q]=size(Ifft);   %Size of Fourier Transformed Image

a=0.001; b=0.01; T=1;  %Parameters

H=zeros(p,q);
Ifft_deg=Ifft;

%Degradation of the Original Image in Frequency Domain
for u=1:p 
    for v=1:q
        H(u,v)=(T/(pi*(u*a+v*b)))*(sin(pi*(u*a+v*b)))*exp(-1i*pi*(u*a+v*b));
        Ifft_deg(u,v)=H(u,v)*Ifft(u,v);
    end
end

subplot(2,3,4);
imshow(uint8(abs(Ifft_deg)));    %Fourier Transform of degraded image
title('FFT of Degraded Image');

I_deg=ifft2(Ifft_deg);  %Inverse Fourier Transform of degraded image

subplot(2,3,5);
imshow(uint8(I_deg));   %Degraded Image before cancelling center shift
title('Degraded Image before cancelling center shift');

Is_deg=ifftshift(I_deg);    %Cancel center Shift
Is_deg=uint8(Is_deg);
subplot(2,3,6);
imshow(Is_deg);  %Degraded Image after cancelling center shift
title('Degraded Image after cancelling center shift');

%Adding noise to the degraded image
%noisy_Is_deg=imnoise(Is_deg,'gaussian',0,0.01);
noisy_Is_deg=imnoise(Is_deg,'salt & pepper',0.0001);

output_image_qvn=double(noisy_Is_deg);
output_image_qvn1=double(noisy_Is_deg);
output_image_qvn2=double(noisy_Is_deg);

output_image_qwmn=double(noisy_Is_deg);
output_image_qwmn1=double(noisy_Is_deg);
output_image_qwmn2=double(noisy_Is_deg);

[m,n]=size(Is_deg);

K1=2000;
K2=0.015;

%-------- Quadratic Volterra Filter Equation---------
for i=2:m-1
    for j=2:n-1  
        output_image_qvn1(i,j)=(3*(output_image_qvn(i,j)^2))-(0.5*output_image_qvn(i+1,j+1)*output_image_qvn(i-1,j-1))...
            -(0.5*output_image_qvn(i+1,j-1)*output_image_qvn(i-1,j+1))-(output_image_qvn(i+1,j)*output_image_qvn(i-1,j))...
        -(output_image_qvn(i,j+1)*output_image_qvn(i,j-1));
        output_image_qvn2(i,j)=output_image_qvn1(i,j)+K1*double(noisy_Is_deg(i,j));
    end
end

%---------Quadratic Weighted Median Filter------------
count=0;
while count<=2
    for i=2:m-1
        for j=2:n-1
            A=[(3*(output_image_qwmn(i,j)^2)),(0.5*output_image_qwmn(i+1,j+1)*output_image_qwmn(i-1,j-1)),...
            (0.5*output_image_qwmn(i+1,j-1)*output_image_qwmn(i-1,j+1)),(1*output_image_qwmn(i+1,j)*output_image_qwmn(i-1,j)),...
        (1*output_image_qwmn(i,j+1)*output_image_qwmn(i,j-1))];
        output_image_qwmn1(i,j)=median(A);
        output_image_qwmn2(i,j)=output_image_qwmn1(i,j)+K2*double(noisy_Is_deg(i,j));
        end
    end
    count=count+1;
end


figure;

subplot(2,2,1);
imshow(uint8(Is_deg));
title('Degraded Image');

subplot(2,2,2);
imshow(uint8(noisy_Is_deg));
title('Noise added to degraded image');

subplot(2,2,3);
imagesc(output_image_qvn2);colormap(gray);
title('QV Filtered Degraded Image');

subplot(2,2,4);
imagesc(output_image_qwmn2);colormap(gray);
title('QWM Filtered Degraded Image');

toc;