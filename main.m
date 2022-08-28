clear; close all; clc

%% Import and Define constants
file2process = ["NoSurfactantFl.csv" "WithSurfactantFl.csv" "NoSurfactantPhoto.csv" "WithSurfactantPhoto.csv"];
yAxisNameVector = ["Fluorescent Emission (Normalised) / a.u." "Fluorescent Emission (Normalised) / a.u." "Absorbance / a.u." "Absorbance / a.u."];

global figureLetterNumber
figureLetterNumber = 65;

%% Make the Plots
maximumMatrix = zeros(12,4);
nameMatrix = zeros(12,4);

figure(1);
for i = 1:2
    subplot(2,2,i)
    [maximumMatrix(:,i),nameMatrix(:,i)] = plotter(file2process(i),[0,0],yAxisNameVector(i));
    figureLetterNumber = figureLetterNumber + 1;
end

for i = 3:4
    subplot(2,2,i)
    [maximumMatrix(:,4),nameMatrix(:,4)] = plotter(file2process(i),[16,31],yAxisNameVector(i));
    % 16 31
    figureLetterNumber = figureLetterNumber + 1;
end

maximumMatrix(maximumMatrix==0) = NaN;
nameMatrix(nameMatrix==0) = NaN;

figure(2)
logMaximumMatrix = log(maximumMatrix);

for i = 1:2
    if isnan(logMaximumMatrix(end,i))
        myFit = polyfit(nameMatrix(1:end-1,i),logMaximumMatrix(1:end-1,i),1);
    else
        myFit = polyfit(nameMatrix(:,i),logMaximumMatrix(:,i),1);
    end
    x = linspace(min(nameMatrix(:,i)),max(nameMatrix(:,i)));
    y = polyval(myFit,x);
    subplot(1,2,i)
    hold on
    plot(nameMatrix(:,i),logMaximumMatrix(:,i),'DisplayName','')
    plot(x,y)
    hold off
    xlabel('Growth Time / s') ; ylabel('ln(\lambda_{max}) / nm')
    legend('',strcat([num2str(myFit(1)*10^4) 'e-04x+' num2str(myFit(2))]))
end

radMatrix = zeros(size(maximumMatrix));
for i = 1:4
    for j = 1:numel(maximumMatrix(:,i))
        radMatrix(j,i) = wav2radius(maximumMatrix(j,i));
    end
end

volMatrix = zeros(size(radMatrix));
for i = 1:4
    for j = 1:numel(radMatrix(:,i))
        volMatrix(j,i) = radius2volume(radMatrix(j,i));
    end
end

massMatrix = zeros(size(volMatrix));
for i = 1:4
    for j = 1:numel(volMatrix(:,i))
        massMatrix(j,i) = volume2mass(volMatrix(j,i));
    end
end

%% Function Definitions
function [maximumList,nameList] = plotter(filename,cutoffs,yAxisName)
    global figureLetterNumber
    table2plot = readtable(filename);

    colours = flip(jet(13));

    if cutoffs == [0,0]
        cutoffs = [1,height(table2plot)-1];
    end

    maximumList = zeros(12,1);
    nameList = zeros(22,1);
    wav = table2array(table2plot(2:end,1));
    hold on

    for i = 2:width(table2plot)
        intensity = table2array(table2plot(2:end,i));

        if table2array(table2plot(1,1)) == 1
            intensity = normalize(intensity,'range');
        end

        wavMod = wav(cutoffs(1):cutoffs(2));
        intensityMod = intensity(cutoffs(1):cutoffs(2));

        [~,maximum] = max(intensityMod);
        maximumList(i-1) = wavMod(maximum);

        %intensityMod = smooth(intensityMod);
        SeriesName = string(table2array(table2plot(1,i)));
        nameList(2*i-3) = SeriesName;
        nameList(2*i-2) = NaN;

        
        plot(wavMod,intensityMod,'Color',colours(i,:))
        xline(wavMod(maximum),'--',SeriesName,'Color',colours(i,:));
        title(char(figureLetterNumber))
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        ax.TitleFontSizeMultiplier = 2;
        
        xlabel('Wavelength / nm'); ylabel(yAxisName);
        grid on
    end
    %legendList = [reshape(num2str(nameList),1,[]);nan(1,numel(num2str(nameList)))];
    %legendList = legendList(:);
    names = cell(12, 1);
%     nameList(2:2:end) = NaN;
    for i = 1:length(nameList)
        names{i} = nameList(i);
    end

    for i = 1:numel(names)
        if ~isnan(names{i})
            names{i} = num2str(names{i});
        else
            names{i} = '';
        end
    end
    legend(names)
    hold off
    
    tempNameList = zeros(12,1);
    j = 1;
    for i = 1:numel(nameList)
        if ~isnan(nameList(i))
            tempNameList(j) = nameList(i);
            j = j + 1;
        end
    end
    nameList = tempNameList;
end

% function outputMatrix = calculateMatrix(function2apply,matrix2apply)
%     outputMatrix = zeros(size(matrix2apply));
%     for i = 1:size(matrix2apply,2)
%         for j = 1:numel(matrix2apply(:,i))
%             outputMatrix(j,i) = function2apply(matrix2apply(j,i));
%         end
%     end
% end

function radius = wav2radius(wav)
    % radius in m
    h = 6.63*10^-34;
    c = 3.00*10^8;
    me = 1.18*10^-31;
    mh = 4.10*10^-31;
    Egap = 2.79*10^-19;

    radius = sqrt((h^2*((1/me)+(1/mh)))/(8*(((h*c)/(wav*10^-9)-Egap))));
end

function volume = radius2volume(radius)
    % volume is in m^3
    volume = (4/3)*pi*radius^3;
end

function mass = volume2mass(volume)
    %volume is in g
    rho = 5.81*10^6;
    mass = rho*volume;
end