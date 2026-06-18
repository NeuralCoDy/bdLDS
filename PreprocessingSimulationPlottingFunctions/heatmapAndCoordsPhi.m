function heatmapAndCoordsPhi(whichData,whichCoords,loadingMatrixPhi,inf_opts,trialIdentifier)

if inf_opts.AcrossIndividuals
    for i = 1:size(whichData,1) % individual
        myPhi = loadingMatrixPhi{i};

        threshold = median(abs(myPhi));

        fig1 = figure();
        cla
        hold on
        colormap(redbluecmap)
        imagesc(myPhi)
        colorbar
        title('Phi')
        hold off
        filename = sprintf("zebgen_%s_Phi_%d.fig", trialIdentifier,i);
        saveas(fig1,filename);

        fig1 = figure();
        cla
        hold on
        colormap("parula")
        imagesc(abs(myPhi))
        colorbar
        title('Abs Phi')
        hold off
        filename = sprintf("zebgen_%s_absPhi_%d.fig", trialIdentifier,i);
        saveas(fig1,filename);
        close all

        fig1 = figure();
        cla
        hold on
        groupNums = kmeans(myPhi,int8(round(size(myPhi,2)/5)));
        [~,idx] = sort(groupNums);
        colormap(redbluecmap)
        %myPhiAugmented = [idx, myPhi];
        %myPhiSorted = sortrows(myPhiAugmented);
        %myPhiSorted = myPhiSorted(:,2:end);
        imagesc(myPhi(idx,:))
        %imagesc(myPhiSorted)
        colorbar
        title('Phi, k means clustered')
        hold off
        filename = sprintf("zebgen_%s_kmeansPhi.fig", trialIdentifier);
        saveas(fig1,filename);
        
        % disp(size(myPhi,2))
        %myPhiToDisplay = myPhi(abs(myPhi)>threshold);
        myPhiToDisplay = myPhi;
        myPhiToDisplay(abs(myPhiToDisplay)<threshold) = 0;
    
        for j = 1:size(myPhi,2) % each latent dim
            
            fig1 = figure();
            cla
            hold on
            colormap(redbluecmap)
            scatter3(whichCoords(:,3),whichCoords(:,2),whichCoords(:,1),50,myPhiToDisplay(:,j))
            colorbar
            %caxis([-2 2])
            %caxis([0 1])
            title('Phi > median of abs(Phi)')
            hold off
            filename = sprintf("zebgen_%s_mapPhiThresholded_latent%d.fig", trialIdentifier,j);
            saveas(fig1,filename);
            close all
        end
    end
    % close all

else

    myPhi = loadingMatrixPhi;

    threshold = 3*std(myPhi,0,"all");

    fig1 = figure();
    cla
    hold on
    colormap(redbluecmap)
    imagesc(myPhi)
    colorbar
    title('Phi')
    hold off
    filename = sprintf("zebgen_%s_Phi.fig", trialIdentifier);
    saveas(fig1,filename);

    fig1 = figure();
    cla
    hold on
    colormap("parula")
    imagesc(abs(myPhi))
    colorbar
    title('Abs Phi')
    hold off
    filename = sprintf("zebgen_%s_absPhi.fig", trialIdentifier);
    saveas(fig1,filename);
    close all

    fig1 = figure();
    cla
    hold on

    groupNums = kmeans(myPhi,int8(round(size(myPhi,2)/5)));
    [~,idx] = sort(groupNums);
    colormap(redbluecmap)
    %myPhiAugmented = [idx, myPhi];
    %myPhiSorted = sortrows(myPhiAugmented);
    %myPhiSorted = myPhiSorted(:,2:end);
    imagesc(myPhi(idx,:))
    %imagesc(myPhiSorted)
    colorbar
    title('Phi, k means clustered')
    hold off
    filename = sprintf("zebgen_%s_kmeansPhi.fig", trialIdentifier);
    saveas(fig1,filename);
    
    % disp(size(myPhi,2))
    %myPhiToDisplay = myPhi(abs(myPhi)>threshold);
    myPhiToDisplay = myPhi;
    myPhiToDisplay(abs(myPhiToDisplay)<threshold) = 0;

    for j = 1:size(myPhi,2) % each latent dim
        
        fig1 = figure();
        cla
        hold on
        colormap(redbluecmap)
        scatter3(whichCoords(:,3),whichCoords(:,2),whichCoords(:,1),50,myPhiToDisplay(:,j))
        colorbar
        %caxis([-2 2])
        %caxis([0 1])
        title('abs(Phi) > 3 std(Phi)')
        hold off
        filename = sprintf("zebgen_%s_mapPhiThresholded_latent%d.fig", trialIdentifier,j);
        saveas(fig1,filename);
        close all
    end
    % close all

end

end