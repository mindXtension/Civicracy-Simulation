function [ output_args ] = TestOpenDemocracy_1( nUnits ) 
%TESTOPENDEMOCRACY  
%what shall be done: 
% set up a structures variable called 'Unit' that stands for a person 
% each unit has: 
%   - an overal IQ, used to better equealize knowlege 
%   - 'Field' fields of knowledge 
%       - abilities ranging from 1 to 100 in each field of knowledge 
%       ** values to implement to make it more human: 
%           - additionally a variation 'V' of how much they believe in their own ability in each specific field of knowledge 
%           - a value how much they can convince others about their abilities 
% 
%       - a pool o trusts which is a minimum of 1 for each field of knowledge 
%   - the ability to gain trust in each F 
%   - the ability to send their pool of trusts to another unit 
%   - a connection 'Con' to a random number 'nCon' of Units they 'know' 
%       ** values to implement to make it more human: 
%           - a factor of how much they know the knowlegde of the connected 
%           unit in each field of knowledge 'depth' 
%   dispThrusts 
nUnits=1000; %set ammount of units to be simulated 
nFields=1; %set amount of fields of knowledge 
nConMax=100; %set maximum connection of a Unit 
nConMin=25; %set minimum connection of a Unit 
nReCon=25; %set number of possible bidirectional connenctions 
Alpha=2.5; %set shape of statistical skilldistribution, 
IQ_Variation=1; %set minor variations of skills (value 1-3) to blur discrete steps 
Unit=generateUnits(nUnits,nFields,nConMax,nConMin,nReCon, Alpha, IQ_Variation); %generate population of 'Units' 
%DispConnections(Unit,1,nUnits,'nconnections') 
%dispThrusts %TODO 
% Lets go for 5 iteration to set up the System, usually it is quite stable 
% after that 
i=0; 
sumTrusts=nUnits; 
while sum(sumTrusts)>0 %while still trusts have to be given, otherwise sytem is define 
%for i=1:5 
    i=i+1; 
    Unit=setTrusts(Unit); 
    [IterationResults{i} Responsibles nSmartest sumTrusts sumRes(i,1:11)]=analyzeUnits(Unit); %get Result for every iteration 
    TrustsPerIteration(i,1:11)=sumTrusts; %stores how many Trusts have to be gioven per Iteration 
    ResponsiblesPerIteration(i,:)=Responsibles; 
end 
disp('number of iterations:') 
disp(i); 
disp('System is set up'); 
DispResults1 (Unit,IterationResults,ResponsiblesPerIteration, nUnits) 
Unit=modifyUnits(Unit,5); 
for i=1:5 
    Unit=setTrusts(Unit); 
    [IterationResults{i} Responsibles]=analyzeUnits(Unit); %get Result for every iteration 
    ResponsiblesPerIteration(i,:)=Responsibles; 
end 
% set collection of Units 
  
  
disp('done') 
end 
  
function [Level Responsibles nSmartest sumTrusts sumRes]=analyzeUnits(Unit) 
%generate a structure that holds the indexed of Units that are above a 
%certain threshold level of Trusts 
%INPUT: unit 
%OUTPUT:  
Levels=[1 2 5 10 20 50 100 200 500 1000 2000 5000]; 
    for i=1:length(Levels)-1 
         
        Level(i).Responsibility=Levels(i); %add value of needed Responibility for each level 
        [Level(i).Responsibles nSmartest]=getLevels(Unit, Levels(i), Levels(i+1)); %get IDs of units with sufficient Responsibility 
        Responsibles(i)=length(find(Level(i).Responsibles)); %counts how many units have been selected for given level and places this value into  variable Responsibles for quick overview 
        sumRes(i)=0; 
        sumTrusts(i)=0; 
        if Level(i).Responsibles~=0 
            for r=Level(i).Responsibles' 
               sumRes(i)=sumRes(i)+Unit(r).Field(1).Responsibility2; 
               sumTrusts(i)=sumTrusts(i)+Unit(r).Field(1).Trust; 
            end 
        end 
    end 
  
end 
  
function [Responsibles nSmartest]=getLevels(Unit, minResponsibility, maxResponsibility) 
%identify units with more than 'minTrust' in all Fields of knowledge 
nUnits=length(Unit); 
nFields=length(Unit(1).Field); 
nSmartest=0; 
for f=1:nFields %for all Fields 
    i=1; 
    for u=1:nUnits 
        if Unit(u).Field(f).Responsibility>=minResponsibility&&Unit(u).Field(f).Responsibility<maxResponsibility %if a Unit has one more than 'minResponsibility' and less than Responsibility of next level Responsibility 
            Responsibles(i,f)=u; %add its index to the Responsibles list) 
            i=i+1; % 
        end 
        if i==1 
            Responsibles(i,f)=0; 
        end 
        if Unit(u).Field(f).smartest==1 
            nSmartest=nSmartest+1; %counts how many units dont know any smarter units 
        end 
  
    end 
end 
end %function 
  
function [Unit]=setTrusts(Unit) 
nUnits=length(Unit); 
nFields=length(Unit(1).Field); 
x=0; 
for u=1:nUnits %for every Unit 
    for f=1:nFields %for all fields of knowledge 
        %all friend shall be checked: abilities can change 
        Unit(u).Field(f).FriendAbilities=zeros(Unit(u).nCon,1); % empty previous ability-results of friens 
        %get abilities of all connections in all fields: 
        %if ~exist(Unit(u).Field(f).FriendAbilities,'var') %if 
        %friendAbilities have not been set (line not used: all friend shall be checked: abilities can change) 
        for c=1:length(Unit(u).Con) %for all connections 
            Friend=Unit(u).Con(c); %get Friend ID 
            Unit(u).Field(f).FriendAbilities(c)=Unit(Friend).Field(f).Ability; % Add friends current ability to ability listing of all friends 
        end 
        %end 
        %TODO 
        Unit(u).Field(f).smartest=0; 
        MaxFriendsAbility=max(Unit(u).Field(f).FriendAbilities); 
        Unit(u).Field(f).Responsibility=Unit(u).Field(f).Responsibility+Unit(u).Field(f).Trust; %add gathered Trusts (or deduct if Trust have been taken) to responsibility to sum up all gathered trust votes (if variable Trust has been set to a value higher than 0, someone gave this unit its trust 
        if MaxFriendsAbility>Unit(u).Field(f).Ability %if any friend has higher ability 
            iFriendToTrust=find(Unit(u).Field(f).FriendAbilities==MaxFriendsAbility); %find the friend with the highest ability value (in real life quite complicated process - so here is work for future simulation studies) 
            if size(iFriendToTrust,1)>1 %if on has friends with similar highest abilitues 
                iFriendToTrust=iFriendToTrust(randi(size(iFriendToTrust,2))); % randomly choose one of them 
            end 
            FriendToTrust=Unit(u).Con(iFriendToTrust); %get ID of highest ranked friend 
            %Unit(u).Field(f).FriendToTrust~=FriendToTrust 
            if Unit(u).Field(f).FriendToTrust~=FriendToTrust&&Unit(u).Field(f).FriendToTrust~=u %if someone else is the new smartest Friend 
                Unit(u).Field(f).Trust=Unit(u).Field(f).Responsibility; %place all Responsibility into own Trusts to pass on Unit(FriendToTrust).Field(f).Trust+Unit(u).Field(f).Trust; 
                Unit(Unit(u).Field(f).FriendToTrust).Field(f).Trust=Unit(Unit(u).Field(f).FriendToTrust).Field(f).Trust-Unit(u).Field(f).Responsibility; %deduct given Trusts from former FriendToTrust, old Best Friend can even have a negative Trusts-value to be taken into account the step it is beeing calculated 
            end 
            %pass gathered Trusts to the friend who has higher ability       
            if (Unit(u).Field(f).Trust>0)&&(Unit(FriendToTrust).Field(f).Ability>Unit(u).Field(f).Ability) % check if ability of connection is above own abiltiy, and trust has not been passed on 
  
                Unit(u).Field(f).FriendToTrust=FriendToTrust;              
                Unit(FriendToTrust).Field(f).Trust=Unit(FriendToTrust).Field(f).Trust+Unit(u).Field(f).Trust; % if so,pass on own Trusts to connected smart Friend-unit 
                Unit(u).Field(f).Responsibility2=0; 
                %Unit(u).Field(f).FriendToTrust=FriendToTrust; 
  
            end %pass on trusts 
        else %if own ability is higher than at all Connections  
            Unit(u).Field(f).smartest=1; 
            Unit(u).Field(f).Responsibility2=Unit(u).Field(f).Responsibility; %Responsibility2 holds the value of units one is beeing Responisible for and cannot pass on his Trusts due to not beeing connected to someone smarter  
        end 
        Unit(u).Field(f).Trust=0; %set a Units Trusts to zero -> vote given, trusts passed on - or added to own responsibility 
        if Unit(u).Field(f).Responsibility>1 
            x=x+1; 
        end 
    end% for all fields 
end% for all units 
disp(x); 
end %function 
  
function [Unit]=modifyUnits(Unit,Variation) 
%change a units ability in a certain field by a random offset 
nUnits=length(Unit); 
nFields=length(Unit(1).Field); 
for u=1:nUnits %for every Unit 
    for f=1:nFields %for all fields of knowledge 
        Unit(u).Field(f).Ability=Unit(u).Field(f).Ability+randn(1)*Variation;    
    end% for all fields 
end% for all units 
end %function 
  
function [Unit]=generateUnits(nUnits,nFields,nConMax,nConMin,nReCon, Alpha, IQ_Variation) 
% Generates a set of Units, compareabe to a group of people 
% INPUT: 
% - nUnits: set number of units 
% - nFields: set number of fields that can be voted for  
% - nConMax: set maximum possible number of connections of a unit (maximum ammount of people  
%            one knows (connections) that are possible) - is currently  
%            fixed, but may slightly change (grow) by the time 
% - nConMin :set minimum number of connections of a unit - (minimum ammount of people  
%            one knows (connections) that are at least possible) is currently fixed, 
%            but may slightly change (grow) by the time 
% - nReCon: Number possible of bidirectional connections (lower than nConMax) 
% - Alpha: vlaue between 2 and 3 to define the characteristic of the skill distribution 
% - IQ_Variation: value that states how much the skill might vary in a certain field (around 1-3)  
SkillDistribution=GenerateSkillDistribution(100,2.5); % set Alpha (2nd Value between 2 and 3) 
Unit=struct; %define a structure 'Unit' 
for u=1:nUnits %go trough all units 
    Unit(u).IQ=SkillDistribution(randi(length(SkillDistribution)))+randn(1); %set the basic IQ of a Unit 
    for f=1:nFields 
        Unit(u).Field(f).Ability=Unit(u).IQ+randn(1)*IQ_Variation; % set ablitities according to maximum IQ-value 
        Unit(u).Field(f).Trust=1; %set the Trust for each ability to 1 
        Unit(u).Field(f).Responsibility=0; 
        Unit(u).Field(1).Responsibility2=0; 
        Unit(u).Field(f).ChangeCredits=1; 
        Unit(u).Field(f).FriendToTrust=u; %first trust goes to each unit itself 
    end 
    % for visuaisation: set up random coordinates 
    window=800; %set size of rectangular window 
    Unit(u).Position= randi(window,3,1); 
end 
% after creation of Units, set up connections between them  
    % 1. initilize connections  
    for u=1:nUnits 
        Unit(u).nCon=nConMin+randi(nConMax-nConMin); %setup number 
        Unit(u).nReCon=[randi(nReCon)]; %setup number of bidiretional connections 
        Unit(u).Con=zeros(Unit(u).nCon,1); %setup empty connections 
    end 
    % 2. set connetions 
    for u=1:nUnits % go through all units 
        iCon=find(Unit(u).Con==0); %get indizes of not set connections 
        if ~isempty(iCon) %if any connections need to be set (some to all connections could have been set via bidirectional connections) 
            for i=1:length(iCon) %go through all empty connections of a unit 
                con=u; 
                while con==u % make sure that con points to another unit than unit:'u' 
                    con=randi(nUnits); %set value for a connection to a random Unit 
                end 
                Unit(u).Con(iCon(i))=con; %set connection to the random unit for the processed unit 
                if iCon(i)<=nReCon %if index in iCon points to an element in the range of bidirectional connections and           
                    iConTarget=find(Unit(con).Con==0); 
                    if ~isempty(iConTarget) %if other unit has still free spots 
                        Unit(con).Con(iConTarget(1))=u; %set procesed units index:'u' as a connection for the randum unit:'con' at the first empty entry 
                    end %otherwise, so when all connections are full nothing happens 
                end 
            end 
        end 
         
        % what if connection allready exists? (more than one connection to 
        % same unit? 
           % #HIER WEITERMACHEN!!! 
        %Unit(u).ReCon=randi(nUnits,[nCon,1]); 
        % 
        % schau ReCons an (die ersten bis nReCon) 
        %   wenn  
  
        %Unit(u).Con=randi(nUnits,[nCon,1]); % set connections to other units MIND: how about closed system which only know themselves (clusters), and its now possible to also 'know yourself'  
    end % go through all units 
end %function generateUnits 
  
function [SkillDistribution]=GenerateSkillDistribution(nSkillSteps,Alpha) 
% generate a function that represents the skillevel in a group of units 
% is based on an gaussian function 
% the function is discretisized by nSkillSteps 
% to add more variation add randn(1) when assigning skill 
N=nSkillSteps; % Number of discrete skill levels 
%Alpha=2.5; % set Alpha between 2 and 3 
x = gausswin(N,Alpha); %generates a gaussian distribution of values 
x1=1./x; %inverts gaussian values 
% hold on 
% for i=1:length(x) 
%     plot(i,x1(i)) 
% end 
% hold off 
m=min(x1); 
% Make SkillDistribution-Shape  
SkillDistribution(1:N/2)=100-1./x(1:N/2)+m; 
SkillDistribution(N/2+1:N)=100+1./x(N/2+1:N)-m; 
end %GenerateSkillDistribution 
  
function DispConnections(Unit,minUnits,maxUnits,what) 
%function to display units and connections in the window between minUnits 
%and maxUnits 
if strcmp(what,'connections')==1 
    hold on 
    for u=minUnits:maxUnits %go through all units 
       %scatter(Unit(u).Position(1),Unit(u).Position(2),'*','r'); %plot unit at unit-own coordinates 
       for i=1:Unit(u).nCon %go through all connetions 
           if ~isempty(find(Unit(Unit(u).Con(i)).Con==i)) 
               C=[1 0 1]; %bidirectional connection 
           else 
               C=[0 0 1]; %unidirectional connection 
           end 
           line([Unit(u).Position(1) Unit(Unit(u).Con(i)).Position(1)],... 
                [Unit(u).Position(2) Unit(Unit(u).Con(i)).Position(2)],... 
                [Unit(u).Position(3) Unit(Unit(u).Con(i)).Position(3)],'color',C); 
       end %for goint through ll connections 
       scatter3(Unit(u).Position(1),Unit(u).Position(2),Unit(u).Position(3),'*','r'); %plot unit as marker at unit-own coordinates 
    end %for going through all units 
    hold off 
end %if connections shall be displayed 
if strcmp(what,'trust')==1 
    nFields=length(Unit(1).Field); 
    hold on 
    for u=minUnits:maxUnits %go through all units 
        for f=1:nFields %for all Field of knowledge 
        end %for filds of knowledge to be displayed 
         
    end %for going through all units 
end %if trusted connections shall be displayed 
end %function DispConnections 
  
function DispResults1 (Unit,IterationResults,ResponsiblesPerIteration, nUnits) 
%create a polar-coordinate graph of all units 
%angle: spread of all units at a certain level alon a circle 
%   - sorted:first those who trust in ones of next level, then the following level and so on 
%currently only for one field 
  
maxR=400; 
cmap=colormap(jet); 
maxIterations=size(IterationResults,2); 
%get maximum of Levles used 
maxLevel=size(IterationResults{maxIterations},2); 
while IterationResults{maxIterations}(1,maxLevel).Responsibles==0 
    maxLevel=maxLevel-1; 
end 
% go throug all Levels, check which units connct to each other higher 
% levels and plot the units sorted according to levels they are connected 
% to 
for iLevel=1:maxLevel-1 %for all Levels 
    r=ResponsiblesPerIteration(maxIterations,iLevel); %ammount of units that stay at this level  
    deltaAlpha=2*pi/r; %increment of angle to display all units in a certain level 
    Alpha=0; 
    LevelCount=zeros(1,maxLevel); 
    for i=1:r %go through all levels of the given level 
        u=IterationResults{maxIterations}(1,iLevel).Responsibles(i);%select unit 
        %check into which level unit is connected to 
        FriendLevel=iLevel+1; 
        while ~ismember(Unit(u).Field(1).FriendToTrust,IterationResults{maxIterations}(1,FriendLevel).Responsibles) 
            FriendLevel=FriendLevel+1; 
        end 
        LevelCount(FriendLevel)=LevelCount(FriendLevel)+1; %increase index of last Value for 'ConnectedToLevel' 
        ConnectedToLevel{FriendLevel}(LevelCount(FriendLevel))=u; 
    end 
  
    for iLevel2=iLevel:maxLevel %go from actual to last level 
        iLevel2 
        size(ConnectedToLevel{iLevel2},2) 
        for i=1:size(ConnectedToLevel{iLevel2},2) %go through all units in ConnectedToLevel  
            iUnit=ConnectedToLevel{iLevel2}(i); %select unit 
            Unit(iUnit).Field(1).R=maxR/iLevel; %set Radius 
            Unit(iUnit).Field(1).Alpha=Alpha;%set Angle 
            Alpha=Alpha+deltaAlpha; %increase Angle 
            plotdata(i,1:2)=[Unit(iUnit).Field(1).Alpha Unit(iUnit).Field(1).R]; 
        end 
         
    end 
    hold off 
end %for all Levels 
  
  
end 
  
function getFriendToTrust(Unit) 
%makes a list of units and trusted connections 
  
for iUnit=1:length(Unit) 
    A(iUnit,1)=Unit(iUnit).Field(1).FriendToTrust; 
    A(iUnit,2)=Unit(iUnit).Field(1).Responsibility; 
    A(iUnit,3)=Unit(iUnit).Field(1).Responsibility2; 
end 
xlswrite('TrustedFriends.xls',A); 
end 
