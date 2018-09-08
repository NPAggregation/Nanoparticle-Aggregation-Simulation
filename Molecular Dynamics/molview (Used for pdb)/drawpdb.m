function drawpdb(flag)
%DRAWPDB Make a simple GUI and render simple PDB molecule files.

% Author: Joe Hicklin
% September 2001

% See http://extremevr.mirarco.org/wiki/index.php/Molecules or
% http://elchem.kaist.ac.kr/jhkwak/OkanaganPdb97/molecule/molecule.html
% for some molecule files.

persistent lst; % the list box

if(nargin == 0)
    files = dir('*.pdb');
    if(isempty(files))
        error('%s\n%s\n%s\n%s','There are no .pdb files in the current directory.',...
            'drawpdb.m requires some files to view.',...
            'I sometimes get mine from:',...
            'http://extremevr.mirarco.org/wiki/index.php/Molecules');
    end
    figure('Color','black')
    set(gca,'Position',[0 0 0.8 1],'visible','off','DataAspectRatio',[1 1 1])
    cameratoolbar;
    lst = uicontrol('Units','normalized', ...
        'Position',[.8, .05,.19,.9],...
        'String',{files.name},'Style','listbox','Callback','drawpdb(1)');
    drawpdb(1)
else
    cla
    light;
    light('Position',[-1 -1 -2]);
    [x,y,z] = sphere(40);
    nm = get(lst,'String');
    fid = fopen(nm{get(lst,'Value')},'r');
    line = fgetl(fid);
    while ischar(line)
        if(strncmp('HETATM',line,6) || strncmp('ATOM',line,4))
            switch(line(14))
            case  'H', color = [0.7 0.7 0.7]; r = 0.6;
            case  'C', color = [0.3 0.3 1.0]; r = 1.0;
            case  'O', color = [0.3 1.0 0.3]; r = 1.0;
            case  'N', color = [1.0 0.3 1.0]; r = 0.8;
            otherwise, color = [1.0 0.0 0.0]; r = 1.0;
            end
            c =  sscanf(line(31:54),'%f %f %f');
            surface('XData',c(1) + r*x,'YData',c(2) + r*y,...
                'ZData',c(3) + r*z,'FaceColor',color,...
                'EdgeColor','none','FaceLighting','gouraud')
        end
        line = fgetl(fid);
    end
    fclose(fid);
end
