function [ trigger ] = generate_trigger(int,start,stop)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

t1=start:int:stop;
t2=(start:int:stop)+1;

trigger=zeros([1,2*numel(t1)]);

for qq=1:2:2*numel(t1)
    
    trigger(qq)=t1((qq+1)/2);
    trigger(qq+1)=t2((qq+1)/2);  
    
end

end

