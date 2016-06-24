function shape = shape3(x,x_ori)
%%

part1= le(x,3*x_ori(end)/6);
part2= ge(x,3*x_ori(end)/6)& le(x,4*x_ori(end)/6);
part3= ge(x,4*x_ori(end)/6)& le(x,5*x_ori(end)/6);
ind=find(part3,1,'first');
part3(ind-1)=1;
part4= ge(x,5*x_ori(end)/6);

shape(part1)=1;
line=1:sum(part2);
adapt=1-(line.^4)/(2*max(line.^4));
shape(part2)=adapt;

line=1:sum(part3);
line=max(line)-line;
adapt=(line.^4)/(2*max(line.^4));
shape(part3)=adapt;
shape(part4)=0;
shape=shape*1000;
% shape=shape*.5;
% close all 
% plot(x,shape)
return;

%Varaible convexity
% part1= le(x,1*x_ori(end)/8);
% part2= ge(x,1*x_ori(end)/8)& le(x,2*x_ori(end)/4);
% part3= ge(x,2*x_ori(end)/4)& le(x,7*x_ori(end)/8);
% ind=find(part3,1,'first');
% part3(ind-1)=1;
% part4= ge(x,7*x_ori(end)/8);
% 
% shape(part1)=1;
% line=1:sum(part2);
% adapt=1-(line.^4)/(2*max(line.^4));
% shape(part2)=adapt;
% 
% line=1:sum(part3);
% line=max(line)-line;
% adapt=(line.^4)/(2*max(line.^4));
% shape(part3)=adapt;
% shape(part4)=0;
% shape=shape*.5;
% 
% return;