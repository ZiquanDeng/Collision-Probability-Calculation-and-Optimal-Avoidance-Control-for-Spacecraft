function gpopsPrint(setup);
%------------------------------------------------------------------%
% Prints the information about a multiple-phase optimal control    %
% problem to a filename contained in the string SETUP.NAME         %
%------------------------------------------------------------------%
% GPOPS Copyright (c) Anil V. Rao and David Benson                 %
%------------------------------------------------------------------%

sizes = setup.sizes;
limits = setup.limits;
numphases = setup.numphases;
filename = strcat([setup.name,'.txt']);
fid = fopen(filename,'w');
fprintf(fid,' _______________________________________________________________________________________________________\n');
fprintf(fid,'|                                                                                                       |\n');
fprintf(fid,'| GPOPS Version 5.1 Running with the hp-Adaptive Radau Pseudospectral Method and SNOPT 7.1              |\n');
fprintf(fid,'|_______________________________________________________________________________________________________|\n');
fprintf(fid,'|                                                                                                       |\n');
fprintf(fid,'| GPOPS Copyright (c) 2008-2011 Anil V. Rao and David Benson.                                           |\n');
fprintf(fid,'|_______________________________________________________________________________________________________|\n');
fprintf(fid,'|                                                                                                       |\n');
fprintf(fid,'| SNOPT (c) Copyright Regents of the University of California and Stanford University.                  |\n');
fprintf(fid,'|_______________________________________________________________________________________________________|\n');
fprintf(fid,'|                                                                                                       |\n');
fprintf(fid,'| Downloading, using, copying, or modifying the GPOPS code constitutes an agreement to ALL of the terms |\n');
fprintf(fid,'| of the GPOPS license. Please see the LICENSE page on the GPOPS website at http://www.gpops.org.       |\n');
fprintf(fid,'|_______________________________________________________________________________________________________|\n');
fprintf(fid,'                                                                                                        \n');
ssfilename = strcat(['| >>>>>>>> Summary of Problem Written to File: ',filename,' <<<<<<< |']);
strdashed1 = '';
for k=1:length(ssfilename)-2
    strdashed1 = strcat([strdashed1,'_']);
end;
strdashed1 = strcat([' ',strdashed1,' ']);
strdashed3 = '';
for k=1:length(ssfilename)-2
    strdashed3 = strcat([strdashed3,'_']);
end;
strdashed3 = strcat(['|',strdashed3,'|']);
disp(strdashed1)
disp(['|' blanks(length(ssfilename)-2),'|']);
disp(ssfilename);
disp(strdashed3);
for i=1:numphases
    linestr1 =          ' __________________________________________\n';
    linestr2 =          '|                                          |\n';
    linestr3 =          '|__________________________________________|\n';
    phasestr = '|>>>>>>>>> Information in Phase %d <<<<<<<<<|\n';
    fprintf(fid,'\n');
    fprintf(fid,linestr1);
    fprintf(fid,linestr2);
    fprintf(fid,phasestr,i);
    fprintf(fid,linestr3);
    fprintf(fid,'\n');
    nstates = sizes(i,1);
    ncontrols = sizes(i,2);
    nparameters = sizes(i,3);
    npaths = sizes(i,4);
    nevents = sizes(i,5);
    if nstates>0,
        xlow1 = limits(i).state.min(:,1);
        xlow2 = limits(i).state.min(:,2);
        xlow3 = limits(i).state.min(:,3);
        xupp1 = limits(i).state.max(:,1);
        xupp2 = limits(i).state.max(:,2);
        xupp3 = limits(i).state.max(:,3);
        for j=1:nstates
            sxlow1 = num2str(xlow1(j));
            sxupp1 = num2str(xupp1(j));
            sxlow2 = num2str(xlow2(j));
            sxupp2 = num2str(xupp2(j));
            sxlow3 = num2str(xlow3(j));
            sxupp3 = num2str(xupp3(j));
            ee = ' <= ';
            sname = strcat(['State ',num2str(j)]);
            s1 = strcat(['\t Start of Phase:    \t',sxlow1,ee,sname,ee,sxupp1,'\n']);
            s2 = strcat(['\t During Phase:      \t',sxlow2,ee,sname,ee,sxupp2,'\n']);
            s3 = strcat(['\t Terminus of Phase: \t',sxlow3,ee,sname,ee,sxupp3,'\n']);
            fprintf(fid,strcat(['State \t',num2str(j),'\n']));
            fprintf(fid,s1);
            fprintf(fid,s2);
            fprintf(fid,s3);
        end;
    else
    	s1 = strcat(['No States in Phase ',num2str(i),'\n']);
        fprintf(fid,s1);
    end;
    fprintf(fid,'\n');
    if ncontrols>0,
        ulow = limits(i).control.min;
        uupp = limits(i).control.max;
        for j=1:ncontrols
          sulow = num2str(ulow(j));
            suupp = num2str(uupp(j));
            ee = ' <= ';
            sname = strcat(['Control ',num2str(j)]);
            s1 = strcat(['\t During Phase:    \t',sulow,ee,sname,ee,suupp,'\n']);
            fprintf(fid,strcat(['Control \t',num2str(j),'\n']));
            fprintf(fid,s1);
        end;
    else
    	s1 = strcat(['No Controls in Phase ',num2str(i),'\n']);
        fprintf(fid,s1);
    end;
    fprintf(fid,'\n');
    if nparameters>0,
        plow = limits(i).parameter.min;
        pupp = limits(i).parameter.max;
        for j=1:nparameters
            splow = num2str(plow(j));
            spupp = num2str(pupp(j));
            ee = ' <= ';
            sname = strcat(['Parameter ',num2str(j)]);
            s1 = strcat(['\t During Phase:    \t',splow,ee,sname,ee,spupp,'\n']);
            fprintf(fid,strcat(['Parameter \t',num2str(j),'\n']));
            fprintf(fid,s1);
        end;
    else
        s1 = strcat(['No Parameters in Phase ',num2str(i),'\n']);
        fprintf(fid,s1);
    end
    fprintf(fid,'\n');
    if npaths>0,
        pathlow = limits(i).path.min;
        pathupp = limits(i).path.max;
        for j=1:npaths
            spathlow = num2str(pathlow(j));
            spathupp = num2str(pathupp(j));
            ee = ' <= ';
            sname = strcat(['Path ',num2str(j)]);
            s1 = strcat(['\t During Phase:    \t',spathlow,ee,sname,ee,spathupp,'\n']);
            fprintf(fid,strcat(['Path Constraint \t',num2str(j),'\n']));
            fprintf(fid,s1);
        end;
    else
        s1 = strcat(['No Path Constraints in Phase ',num2str(i),'\n']);
        fprintf(fid,s1);
    end
    fprintf(fid,'\n');
    if nevents>0,
        eventlow = limits(i).event.min;
        eventupp = limits(i).event.max;
        for j=1:nevents
            seventlow = num2str(eventlow(j));
            seventupp = num2str(eventupp(j));
            ee = ' <= ';
            sname = strcat(['Event ',num2str(j)]);
            s1 = strcat(['\t\t\t\t',seventlow,ee,sname,ee,seventupp,'\n']);
            fprintf(fid,strcat(['Event Constraint \t',num2str(j),'\n']));
            fprintf(fid,s1);
        end;
    else
        s1 = strcat(['No Event Constraints in Phase ',num2str(i),'\n']);
        fprintf(fid,s1);
    end    
    fprintf(fid,'\n');
    if isfield(limits(i),'duration'),
        durationlow = limits(i).duration.min;
        durationupp = limits(i).duration.max;
        sdurationlow = num2str(durationlow);
        sdurationupp = num2str(durationupp);
        ee = ' <= ';
        s1 = strcat(['Phase Duration: ',sdurationlow,ee,'duration',ee,sdurationupp,'\n']);
        fprintf(fid,s1);
    else
        s1 = strcat(['No Limits on Phase Duration in Phase ',num2str(i),'\n']);
        fprintf(fid,s1);
    end
    fprintf(fid,'\n');
end;
fclose(fid);
