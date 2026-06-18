function plotSummaryCandReconstr_OneC(whichData,whichPhi,whichDynCoeffs,whichLatentStates,inf_opts,trialIdentifier,varExplPerNeuron,samplingRate, varargin)
behavior1 = [];
behavior2 = [];
listOfCs  = [];
if nargin > 11
    dataIdentifier = varargin{1};
    behavior1 = varargin{2};
    behavior2 = varargin{3};
    listOfCs  = varargin{4};
    ifBehavior1 = 1;
    ifBehavior2 = 1;
elseif nargin > 10
    dataIdentifier = varargin{1};
    behavior1 = varargin{2};
    behavior2 = varargin{3};
    ifBehavior1 = 1;
    ifBehavior2 = 1;
elseif nargin > 9
    behavior1 = varargin{2};
    dataIdentifier = varargin{1};
    ifBehavior1 = 1;
    ifBehavior2 = 0;
elseif nargin > 8
    dataIdentifier = varargin{1};
    ifBehavior1 = 0;
    ifBehavior2 = 0;
else
    dataIdentifier = "zebgen";
    ifBehavior1 = 0;
    ifBehavior2 = 0;
end

if isempty(behavior1)
    ifBehavior1 = 0;
end

if isempty(behavior2)
    ifBehavior2 = 0;
end

if inf_opts.AcrossIndividuals
    for ii = 1:size(whichPhi,1)

        disp('Fixme: Not updated')
        fig1 = figure();
        subplot(3,1,1)
        c_values = cell2mat(whichDynCoeffs{ii});
        plot(c_values.')
        hold on
        ylabel("c")
        hold off
 
        [~,whichTraces] = maxk(varExplPerNeuron(:,ii),10);
        subplot(3,1,2)
        x_values = cell2mat(whichLatentStates{ii});
        reconstrY = (whichPhi{ii} * x_values).';
        plot(reconstrY(:,whichTraces))
        hold on
        ylabel("reconstr. fluor.")
        hold off
  
        subplot(3,1,3)
        origY = whichData{ii}.';
        plot(origY(whichTraces))
        hold on
        ylabel("scaled fluor.")
        hold off

        set(gcf,"Color", "white")
        %filename = input("figure filename:","s");
        filename = sprintf("%s_%s_%d.fig", dataIdentifier, trialIdentifier,ii);
        saveas(fig1,filename);

    end
else
    for ii = 1:size(whichLatentStates,1) %in the event of multiple trials

        c_values_all = whichDynCoeffs{ii};

        if isempty(listOfCs)
            listOfCs = 1:size(c_values_all,1);
        end

        for whichC = listOfCs %1:size(c_values_all,1) 
            c_values        = c_values_all(whichC,:);

            
            x_values        = whichLatentStates{ii};
            reconstrY       = (whichPhi * x_values).';
            origY           = whichData{ii}.';
            nTimepoints     = size(x_values,2);

            corrC       = corr(c_values.',reconstrY); % one for each neuron

            [~,whichTraces] =  maxk(corrC.^2,5);

    
            tVals = linspace(0,nTimepoints/samplingRate,nTimepoints);
    
            % fig1 = figure();
            % subplot(3,1,1)
            % c_values = whichDynCoeffs{ii};
            % plot(c_values.')
            % hold on
            % ylabel("c")
            % hold off
            % 
            % [~,whichTraces] = maxk(varExplPerNeuron,5); 
            % subplot(3,1,2)
            % x_values = whichLatentStates{ii};
            % reconstrY = (whichPhi * x_values).';
            % plot(reconstrY(:,whichTraces))
            % hold on
            % ylabel("reconstr. fluor.")
            % hold off
            % 
            % subplot(3,1,3)
            % origY = whichData{ii}.';
            % plot(origY(:,whichTraces))
            % hold on
            % ylabel("scaled fluor.")
            % set(gcf,"Color", "white")
            % hold off
    
            fig1 = figure();
            fig1.WindowState = 'maximized';
            
            % subplot: reconstruction    
            maxvalDFF = max(origY,[],'all');
            stdzdOn = true;
            if stdzdOn
                dataReconstructedRescaled = maxvalDFF .*reconstrY;
                dFFRescaled = maxvalDFF .* origY;
            else
                dataReconstructedRescaled = reconstrY;
                dFFRescaled = origY;
            end
            % dFF = dFFRescaled;
            % dataReconstructed = dataReconstructedRescaled;
        
                
            dFF = checkAndFlip(origY,tVals);
            dFFRescaled = checkAndFlip(dFFRescaled,tVals);
            dataReconstructed = checkAndFlip(reconstrY,tVals);
            dataReconstructedRescaled = checkAndFlip(dataReconstructedRescaled,tVals);
        
            neurSel = whichTraces;
    
            if ifBehavior2
                tiledlayout(6,1)
            elseif ifBehavior1
                tiledlayout(5,1)
            else
                tiledlayout(4,1)
            end
            % subplot(5,1,1),cla;
            nexttile
            hold on;
            bottom = 0;
            top = max(dFFRescaled,[],'all');
            plot(tVals, dFFRescaled(:, neurSel)); % best reconstructed
        %     xlabel('Time (s)');
            axis([0 max(tVals) bottom top]);
            ylabel('Fluor. (AU)');
        %     title({'';'Original Data Y'});
            set(gca, 'TickDir', 'out');
            set(gca,'XTickLabel',[]);
            box off;
            hold off;
            
            % subplot(5,1,2),cla;
            nexttile
            hold on;
            bottom = 0;
            if stdzdOn
                top = max(dFFRescaled,[],'all');
                    plot(tVals, dataReconstructedRescaled(:, neurSel));
            end
            axis([0 max(tVals) bottom top]);
            ylabel('Est. Fluor.');
        %     title({'';'Estimated Y'});
            set(gca, 'TickDir', 'out');
            set(gca,'XTickLabel',[]);
            box off;
            hold off;
        
            % x traces
            % subplot(5,1,3),cla;
            nexttile
            hold on;
            bottom = min(x_values,[],'all');
            top = max(x_values,[],'all');
            plot(tVals, x_values); 
        %     xlabel('Time (s)');
            axis([0 max(tVals) bottom top]); 
            ylabel('Latent states (x)');
            set(gca, 'TickDir', 'out');
            set(gca,'XTickLabel',[]);
            box off;
            hold off;
        
            % c traces
       
            cmapMaxDistinct = distinguishable_colors(size(c_values,2));
            % subplot(5,1,4), cla;
            nexttile
            hold on;
            % bottom = min(c_values,[],'all');
            % top = max(c_values,[],'all');
            [row0, ~] = find(abs(c_values) < 0.0001*max(abs(c_values(:)),[],'all')); % find small values
            GC = groupcounts(row0); % tally up small values in each row
            % for jj=size(c_values,1)
            Y = c_values;
            if isempty(GC) | GC~=size(c_values,2) % check: are all the values small? if not, allow to plot
                plot(tVals,Y,'color',cmapMaxDistinct(whichC,:),'DisplayName',num2str(whichC));
            end
            hold on;
            % end
            axis([0 max(tVals) -2 2]);
            set(gca, 'TickDir', 'out');
            if ~ifBehavior1
                xlabel('Time (s)');
            else
                set(gca,'XTickLabel',[]);
            end
            box off;
            % leg = legend([],'Orientation','horizontal');
            leg = legend('NumColumns',10,'Location','southOutside');
            leg.ItemTokenSize = [2,2];
            
            % legend([],'Orientation','vertical')
            ylabel('Dyn. coeffs (c)');
            hold off;
    
            
            if ifBehavior1
                nexttile
                hold on;
                plot(tVals, behavior1); 
                if ~ifBehavior2
                    xlabel('Time (s)');
                else
                    set(gca,'XTickLabel',[]);
                end
                axis([0 max(tVals) -Inf Inf]);
                ylabel('Behavior (motor signal)');
                set(gca, 'TickDir', 'out');
                box off;
                hold off;
            end
    
            if ifBehavior2
                nexttile
                hold on;
                plot(tVals, behavior2); 
                xlabel('Time (s)');
                axis([0 max(tVals) -Inf Inf]);
                ylabel('Visual velocity');
                set(gca, 'TickDir', 'out');
                box off;
                hold off;
            end
    
            fontsize(fig1, scale=0.8) 
            filename = sprintf("%s_%s_DynCoef%d.fig", dataIdentifier, trialIdentifier,whichC);
            saveas(fig1,filename);

           
            filename2 = sprintf("%s_%s_DynCoef%d.pdf", dataIdentifier, trialIdentifier,whichC);
            print('-painters','-dpdf','-fillpage',filename2)
            % % filename2 = sprintf("%s_%s.svg", dataIdentifier, trialIdentifier);
            % % print('-painters','-dsvg',filename2)
        end
    end

end
end
