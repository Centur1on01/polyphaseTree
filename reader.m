% a = fopen('result.txt');
% 
% b = textscan(a, '%f\t\t%f');
% fs = 100000;
% figure();
% plot(pwelch((b{1})));
% hold on;
% %figure();
% plot(pwelch((b{2})));
% hold off;
% 
% 
% figure();
% as = fopen('sine.txt');
% 
% %bs = textscan(as, '%f\t\t%f');
% bs = textscan(as, '%f');
% %plot(pwelch((bs)));
% [pxx,fss] = pwelch(bs);
% plot(pwelch((bs)));
% %plot(fss,10*log10(pxx))

% figure();
% sines = dlmread('sine.txt');
% resss = dlmread('result.txt');
% ii =1:1000;
jj =1:499;
% plot(ii, sines);
% hold on;
% 
% plot(jj, resss);
% hold off;


%a = dlmread('result_back.txt');
%aa = dlmread('result_front.txt');
aaa = fopen('out.txt');
%aaaa = textscan(aaa, '%f\t\t%f');
%aaaa = textscan(aaa, '%f\t\t%f\t\t%f\t\t%f');
aaaa = textscan(aaa, '%f\t\t%f\t\t%f\t\t%f\t\t%f\t\t%f\t\t%f\t\t%f');
%aaaa = textscan(aaa, '%f\t\t\t%f\t\t\t%f\t\t\t%f\t\t\t%f\t\t\t%f\t\t\t%f\t\t\t%f');
b = dlmread('sum.txt');

figure();
plot(pwelch(b));

figure();
hold on
plot(aaaa{1});
plot(aaaa{8});
hold off;



% figure();
% hold on;
% plot(jj, a);
% plot(jj, aa);
% hold off;

figure();
plot(0:396800-1, b);


figure();
hold on;
subplot(2,1,1);
plot(pwelch(aaaa{1}));
subplot(2,1,2);
plot(pwelch(aaaa{2}));
hold off;

% 
% figure();
% plot(0:49599-1,aaaa{1});
% hold off;
% figure();
% plot(0:49599-1,aaaa{2});
% figure();
% plot(0:49599-1,aaaa{3});
% figure();
% plot(0:49599-1,aaaa{4});
% figure();
% plot(0:49599-1,aaaa{5});
% figure();
% plot(0:49599-1,aaaa{6});
% figure();
% plot(0:49599-1,aaaa{7});
figure();
plot(0:49599-1,aaaa{8});

figure();
hold on;
subplot(2,1,1);
plot(0:396800-1, b);
subplot(2,1,2);
plot(0:49599-1,aaaa{8});
grid on;
hold off;


figure();
hold on;
grid on;
subplot(2,1,1);
plot(pwelch(b));
subplot(2,1,2);
plot(pwelch(aaaa{8}));
hold off;

% 
% figure();
% [H1,f] = freqz(b{1},1,512,2);
% plot(f,angle(H1)*180/pi);
% hold on;
% [H2,f] = freqz(b{2},1,512,2);
% plot(f,angle(H2)*180/pi);
% hold off;
% 
% 








