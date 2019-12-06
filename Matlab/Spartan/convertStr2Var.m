%% convert a string to a variable name
% https://www.mathworks.com/matlabcentral/answers/57818-using-a-matlab-string-as-a-valid-variable-name

X = 'omega';
eval([X '= [2 3 4 5]']);