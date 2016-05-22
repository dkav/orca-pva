% Orca population viability analysis using projection matrix
% Filename: whpvacc.m
% Author: Darren Kavanagh
% Date: 11.12.2008

clear all
load swhlfhs % load life history


startyr = 1992; % start year
endyr = 2001;   % end year
range = (startyr:1:endyr); % time range
noyr = (endyr-startyr);    % no years
lfhstbl = zeros(noyr,6);   % life history table
morttbl = zeros(noyr,6);   % mortality table 
movouttbl = zeros(noyr,6); % individuals moving in table
movintbl = zeros(noyr,6);  % individuals moving out table

% loop through time range
for i=1:(noyr)
    curyr = range(i);
    
    % loop through and place in age classes
    % and identify where they are in the following year
    for j=1:(length(swhlfhs))        
        birthyr = swhlfhs(j,3);
        deathyr = swhlfhs(j,4);
        
        % check to see if individual is alive at current year
        if (birthyr <= curyr) && (deathyr > curyr)
            whage = curyr-birthyr;
            gender = swhlfhs(j,2);
            
            % Calf
            if whage == 0
                lfhstbl(i,1) = lfhstbl(i,1)+1;
                if deathyr == curyr+1
                   morttbl(i,1) = morttbl(i,1)+1;
                else
                   movouttbl(i,1) = movouttbl(i,1)+1; 
                   movintbl(i,2) = movintbl(i,2)+1;
                end  
            % Juvenile    
            elseif (whage > 0) && (whage < 11)
                lfhstbl(i,2) = lfhstbl(i,2)+1;
                if deathyr == curyr+1
                   morttbl(i,2) = morttbl(i,2)+1;
                elseif whage == 10 
                   movouttbl(i,2) = movouttbl(i,2)+1;
                   if gender ==2
                    movintbl(i,3) = movintbl(i,3)+1; 
                   else
                    movintbl(i,5) = movintbl(i,5)+1; 
                   end
                end
            % Female reproductive    
            elseif (whage > 10) && (whage < 42) && (gender == 2)
                lfhstbl(i,3) = lfhstbl(i,3)+1;
                if deathyr == curyr+1
                   morttbl(i,3) = morttbl(i,3)+1;
                elseif whage == 41 
                   movouttbl(i,3) = movouttbl(i,3)+1;
                   movintbl(i,4) = movintbl(i,4)+1;
                end
            % Female post reproductive    
            elseif (whage > 41) && (gender == 2)
                lfhstbl(i,4) = lfhstbl(i,4)+1;
                if deathyr == curyr+1
                   morttbl(i,4) = morttbl(i,4)+1;
                end
            % Male young    
            elseif (whage > 10) && (whage < 22) && (gender ~= 2)
                lfhstbl(i,5) = lfhstbl(i,5)+1;
                if deathyr == curyr+1
                   morttbl(i,5) = morttbl(i,5)+1;
                elseif whage == 21 
                   movouttbl(i,5) = movouttbl(i,5)+1;
                   movintbl(i,6) = movintbl(i,6)+1;      
                end
            % Male old    
            elseif (whage > 21) && (gender ~= 2)
                lfhstbl(i,6) = lfhstbl(i,6)+1;
                if deathyr == curyr+1
                   morttbl(i,6) = morttbl(i,6)+1;
                end
            end
        end
    end
end

projectcell = cell(noyr,1);

% create projection matrix from each year
for i=1:(noyr-1)
  projmatx = zeros(6, 6);  
% Calf to Juvenile
  projmatx(2,1) = movintbl(i,2) / lfhstbl(i,1);  

% Juvenile to Juvenile
  projmatx(2,2) = (lfhstbl(i,2)-movouttbl(i,2)-morttbl(i,2)) / ... 
                    lfhstbl(i,2);

% Junvenile to Reproductive Female
  projmatx(3,2) = movintbl(i,3) / lfhstbl(i,2);

% Junvenile to Young Male
  projmatx(5,2) = movintbl(i,5) / lfhstbl(i,2);
  
% Reproductive Female to Calf
  projmatx(1,3) = lfhstbl(i+1,1) / lfhstbl(i,3);

% Reproductive Female to Reproductive Female
  projmatx(3,3) = (lfhstbl(i,3)-movouttbl(i,3)-morttbl(i,3)) / ... 
                    lfhstbl(i,3);
                
% Reproductive Female to Post-Reproductive Female
  projmatx(4,3) = movintbl(i,4) / lfhstbl(i,3);

% Post-Reproductive Female to Post-Reproductive Female
  projmatx(4,4) = (lfhstbl(i,4)-morttbl(i,4)) / lfhstbl(i,4);

% Young Male to Young Male
  projmatx(5,5) = (lfhstbl(i,5)-movouttbl(i,5)-morttbl(i,5)) / ... 
                    lfhstbl(i,5);
                
% Young Male to Old Male
  projmatx(6,5) = movintbl(i,6) / lfhstbl(i,5);
  
% Old Male to Old Male
  projmatx(6,6) = (lfhstbl(i,6)-morttbl(i,6)) / lfhstbl(i,6);

  projmatx(isnan(projmatx)) = 0;
 
  projectcell(i)= {projmatx};
end

% Population trajectory
% Multiple matrices approach
% --------------------------
tfutpop = zeros(50,100);

for i=1:100
   futpop = zeros(50,6);
   futpop(1,:) = lfhstbl(end,:); 
   tfutpop(1,i) = sum(futpop(1,:));
   
    for j=1:49
      rcellno = floor(((noyr-1)-1+1)*rand+1);   
      futpop(j+1,:) = (projectcell{rcellno,1} * (futpop(j,:))')'; 
      tfutpop(j+1,i) = sum(futpop(j+1,:));
    end          
end  

mtfutpop = mean(tfutpop,2);
stdpop = std(tfutpop,[],2)';
se=  (stdpop')/ (sqrt(50));
utfutpop = (mtfutpop)+ 1.96 * se;
ltfutpop = (mtfutpop)- 1.96 * se;

figure (1)
plot(1:1:50,mtfutpop,1:1:50, utfutpop,1:1:50,ltfutpop )
    xlabel('Years')
    ylabel('Population Size')
    box off
    legend('Average','Upper 95% CL', 'Upper 95% CL', 'Location', 'SouthOutside', ...
             'Orientation','horizontal');
    legend('boxoff')
    tstr = {'Multiple matrices approach: ',strcat(num2str(startyr), ' -',num2str(endyr))};
    title(tstr)
    

% Population trajectory
% Variable entries approach
% --------------------------

tfutpop = zeros(50,100);

for i=1:1000
   futpop = zeros(50,6);
   futpop(1,:) = lfhstbl(end,:); 
   tfutpop(1,i) = sum(futpop(1,:));
   
   
    varmatrix = zeros(6,6);
    rcellno = floor(((noyr-1)-1+1)*rand+1);
    varmatrix(2,1) = projectcell{floor(((noyr-1)-1+1)*rand+1),1}(2,1);
    varmatrix(2,2) = projectcell{floor(((noyr-1)-1+1)*rand+1),1}(2,2);
    varmatrix(1,3) = projectcell{floor(((noyr-1)-1+1)*rand+1),1}(1,3);
    varmatrix(3,2) = projectcell{floor(((noyr-1)-1+1)*rand+1),1}(3,2);
    varmatrix(5,2) = projectcell{floor(((noyr-1)-1+1)*rand+1),1}(5,2);
    varmatrix(3,3) = projectcell{floor(((noyr-1)-1+1)*rand+1),1}(3,3);
    varmatrix(4,3) = projectcell{floor(((noyr-1)-1+1)*rand+1),1}(4,3);
    varmatrix(4,4) = projectcell{floor(((noyr-1)-1+1)*rand+1),1}(4,4);
    varmatrix(5,5) = projectcell{floor(((noyr-1)-1+1)*rand+1),1}(5,5);
    varmatrix(6,5) = projectcell{floor(((noyr-1)-1+1)*rand+1),1}(6,5);
    varmatrix(6,6) = projectcell{floor(((noyr-1)-1+1)*rand+1),1}(6,6);
    
    for j=1:49  
      futpop(j+1,:) =  (varmatrix * (futpop(j,:))')'; 
      tfutpop(j+1,i) = sum(futpop(j+1,:));
    end          
end 

mtfutpop = mean(tfutpop,2);
stdpop = std(tfutpop,[],2)';
se=  (stdpop')/ (sqrt(50));
utfutpop = (mean(tfutpop,2))+ 1.96 * se;
ltfutpop = (mean(tfutpop,2))- 1.96 * se;

figure (2)
plot(1:1:50,mean(tfutpop,2),1:1:50, utfutpop,1:1:50,ltfutpop )
    xlabel('Years')
    ylabel('Population Size')
    box off
    legend('Average','Upper 95% CL', 'Upper 95% CL', 'Location', 'SouthOutside', ...
             'Orientation','horizontal');
    legend('boxoff')
    tstr = {'Variable entries approach: ',strcat(num2str(startyr), ' -',num2str(endyr))};
    title(tstr)
    axis ([0 50 0 1000])

