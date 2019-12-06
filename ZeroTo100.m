function [F R]=ZeroTo100(PF,PR,s)

%Author: Hossein Mokhtarzadeh
% to make 0-100% For any any two Parameters PF and PR. 
% you can simply use PF and PR the same ro F and R would be similar.
% Example
% [F R]=ZeroTo100(PF,PR,100);
S1=size(PF);
S2=size(PR);
for i=1:S1(2)
   
    x=1:S1(1);
     xx = linspace(x(1),x(end),s+1);
        yy(:,i) = spline(x,PF(:,i));
        F(:,i)=(ppval(yy(:,i),xx))';
        
      
        
end
    
    
    for i=1:S2(2)
   
    x=1:S2(1);
     xx = linspace(x(1),x(end),s+1);
              
         yyR(:,i) = spline(x,PR(:,i));
        R(:,i)=(ppval(yyR(:,i),xx))';
        
        
end

    