s ='C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\JC\';
load(fullfile(s,'trcResults\TRCs.mat'))
for t=1:length(TRCs)% which trial
    figure
n = TRCs(t).name;
a=TRCs(t).data.gamepad1-TRCs(t).data.gamepad1(1,1:3);

time = TRCs(t).data.time;
% plot(a)

x=(TRCs(t).data.gamepad1(:,1)-TRCs(t).data.gamepad1(1,1)).^2;
y =(TRCs(t).data.gamepad1(:,2)-TRCs(t).data.gamepad1(1,2)).^2; 
z=(TRCs(t).data.gamepad1(:,3)-TRCs(t).data.gamepad1(1,3)).^2;
d = sqrt(x+y+z);
subplot(2,1,1);plot(time,d);title(n);xlabel('Second');ylabel('Reaction Disp(mm)')
subplot(2,1,2);plot(time,gradient(d));title('Speed of Reaction');xlabel('Second');ylabel('Reaction Speed(mm/s)')
end