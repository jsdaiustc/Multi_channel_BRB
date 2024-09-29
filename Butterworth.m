function y=Butterworth(y,fp,fs)
[N,fw]=buttord(fp/500,fs/500,3,20);
[b,a]=butter(N,fw,'low');
y=filter(b,a,y);