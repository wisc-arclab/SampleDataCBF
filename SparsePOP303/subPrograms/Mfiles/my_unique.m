function [C, ia, ic] = my_unique(A, msg)
if ~strcmp(msg, 'rows');
	error('msg should be rows in my_unique.');
end
if verLessThan('matlab', '8.0.1')
	% This part is the same as unique in R2012b or earlier version
	[C, ia, ic] = unique(A, msg);
else
	% This part is for unique in R2013a or later version
	% We use legacy mode of unique.
	[C, ia, ic] = unique(A, msg, 'last', 'legacy');
end
return
