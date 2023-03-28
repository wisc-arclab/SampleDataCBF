% Compiling Mex Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is a component of SparsePOP 
% Copyright (C) 2007 SparsePOP Project
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

if exist('verLessThan') ~= 2
	error('Get verLessThan.m');
end


if verLessThan('matlab', '7.3') 
	MexFlags = ' -O -Dlinux=0 ';
elseif strcmp(computer, 'GLNXA64') || strcmp(computer, 'MACI64')
	MexFlags = ' -O -Dlinux=1 -largeArrayDims ';
elseif strcmp(computer, 'GLNX86')  || strcmp(computer, 'MACI')
	MexFlags = ' -O -Dlinux=0 ';
elseif strcmp(computer, 'PCWIN64') 
	MexFlags = ' -O -Dlinux=1 -largeArrayDims ';
elseif strcmp(computer, 'PXWIN') 
	MexFlags = ' -O -Dlinux=0 ';
else 
	MexFlags = ' -O -Dlinux=0 ';
end

LIBfiles = ' conversion.cpp spvec.cpp polynomials.cpp sup.cpp clique.cpp mysdp.cpp Parameters.cpp ';
if ispc % Windows family create .obj files
        OBJfiles = strrep(LIBfiles,'.cpp','.obj');
else
        OBJfiles = strrep(LIBfiles,'.cpp','.o');
end

eval('cd subPrograms/Mex');
fprintf('Compiling Libraries...');
command = ['mex -c ' MexFlags LIBfiles];
eval(command);
fprintf('done\n');
fprintf('Generating mexconv1...');
command = ['mex ' MexFlags ' mexconv1.cpp'  OBJfiles ];
eval(command);
fprintf('done\n');
fprintf('Generating mexconv2...');
command = ['mex ' MexFlags ' mexconv2.cpp ' OBJfiles ];
eval(command);
fprintf('done\n');
eval('cd ../../');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mex files of SparseCoLO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if verLessThan('matlab', '7.3') 
	MexFlags = ' -O -Dlinux=0 ';
elseif strcmp(computer, 'GLNXA64') || strcmp(computer, 'MACI64')  
	MexFlags = ' -O -Dlinux=1 -largeArrayDims ';
elseif strcmp(computer, 'GLNX86')  || strcmp(computer, 'MACI')
	MexFlags = ' -O -Dlinux=0 ';
elseif strcmp(computer, 'PCWIN64')
	MexFlags = ' -O -Dlinux=1 -largeArrayDims ';
elseif strcmp(computer, 'PCWIN')
	MexFlags = ' -O -Dlinux=0 ';
else
	MexFlags = ' -O -Dlinux=0 ';
end

LIBfiles = [' ccputime.cpp'];
if ispc % Windows family create .obj files
        OBJfiles = strrep(LIBfiles,'.cpp','.obj');
else
        OBJfiles = strrep(LIBfiles,'.cpp','.o');
end

eval('cd V260SubPrograms/SparseCoLO/mex'); 
fprintf('Compiling Libraries...');
command = ['mex -c ' MexFlags LIBfiles];
eval(command);
fprintf('done\n');

clear mexFiles

mexFiles{1} = 'mexForestConvert.cpp';
mexFiles{2} = 'mexMaxSpanningTree2.cpp';
mexFiles{3} = 'mexPrimalOneSDP2.cpp';
% mexFiles{4} = 'mexArrowTriDQOP.cpp';
% mexFiles{5} = 'mexDiagTriDQOP.cpp';
for i=1:length(mexFiles)
    mexFileName = mexFiles{i};
    fprintf('Compiling %s...',mexFileName);
    command = ['mex ' MexFlags mexFileName OBJfiles];
    eval(command);
    fprintf('done\n');
end

eval('cd ../../../');

fprintf('Compilation finished successfully.\n');
