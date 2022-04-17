fs = 1000;
x=load('s0064lrem');
y=x.val/2000;
y=y';
y1=y(1:5*fs,1);
plot(y1);

figure;

fresult=fft(y1);
fresult(1 : round(length(fresult)*5/1000))=0;
fresult(end - round(length(fresult)*5/1000) : end)=0;
yfil=real(ifft(fresult));
plot(yfil);

WinSize = floor(1000 * 300 / 1000);
if rem(WinSize,2)==0
    WinSize = WinSize+1;
end
filtered1=ecgdemowinmax(yfil, WinSize);

figure();
plot(filtered1);

peaks1=filtered1/(max(filtered1)/7);

figure();
subplot(2,1,1);
plot(peaks1);
    %   Filter by threshold filter
for data = 1:1:length(peaks1)
    if peaks1(data) < 4
        peaks1(data) = 0;
    else
        peaks1(data)=1;
    end
end

subplot(2,1,2);
plot(peaks1);

positions=find(peaks1);
figure();
plot(positions);
distance=positions(2)-positions(1);
    %   Returns minimum distance between two peaks
for data=1:1:length(positions)-1
    if positions(data+1)-positions(data)<distance 
       distance=positions(data+1)-positions(data);
    end
end

QRdistance=floor(0.04*1000);
if rem(QRdistance,2)==0
   QRdistance=QRdistance+1;
end
WinSize=2*distance-QRdistance;

filtered2=ecgdemowinmax(yfil, WinSize);
figure();
plot(filtered2);
peaks2=filtered2/(max(filtered2)/7);
figure();
plot(peaks2);
for data=1:1:length(peaks2)
    if peaks2(data)<4
       peaks2(data)=0;
    else
       peaks2(data)=1;
    end
end

positions2=find(peaks2);
distanceBetweenFirstAndLastPeaks = positions2(length(positions2))-positions2(1);
 
averageDistanceBetweenPeaks = distanceBetweenFirstAndLastPeaks/(length(positions2)-1);
    
averageHeartRate = 60 * 1000/averageDistanceBetweenPeaks;
    
disp('Average Heart Rate = ');
disp(averageHeartRate);

p = averageHeartRate;
last=0;
upflag=0;
pulse = zeros(length(filtered2),1);
len = length(positions2);
op = 0;
for i=1:1:length(filtered2)
    if(upflag == 0)
        if(op < len && op ~= 0)
            p = (1/(positions2(op+1)-positions2(op)))*1000*60;
            upflag = 1;
        end
    end
    
    pulse(i) = p;
    if(op < len)
        if(i == positions(op+1))
            op = op+1;
            upflag = 0;
        end
    end
end

figure();
plot(pulse);