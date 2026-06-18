function batch_PSID_Sim(whichFile)
% Based on a script by Bryan Tseng (btseng2@jh.edu) and the PSID repository, Sani et al. 2021. Modeling behaviorally relevant neural dynamics enabled by preferential subspace identification. Nature Neuroscience. https://doi.org/10.1038/s41593-020-00733-0, https://github.com/ShanechiLab/PSID
% To run:
% batch_PSID_Sim('saveFMT_1bhvgraded_forCLDS_GOOD_251113.mat')
%% Add PSID utilities
addpath(genpath('.'));
addpath(genpath('C:/Users/eyeze/Documents/GradSchool/JohnsHopkinsStudent/CharlesLab/PSID/'))
base_path = '.';
save_path = './PSIDcomparison/';
save_dir = fullfile(save_path, 'Sim');
if ~exist(save_dir, 'dir'); mkdir(save_dir); end

load(whichFile)
%%
% Conditions
% condition_files = {
%     'BuH_yaw_active_like', 'BuH_yaw_sine', ...
%     'HoB_yaw_active', 'HoB_yaw_active_like', 'HoB_yaw_active_multiple_perturbations', ...
%     'HoB_yaw_sine', 'WB_yaw_active_like', 'WB_yaw_active_like_multilevel', 'WB_yaw_sine'
%     };

% cluster_numbers = [213 205 219 45 67 214 212 217 142 151 158 165];
n_ahead = 50;  
% nx = 50; n1 = 5;
nx = 8; n1 = 5;

bhv_steps = [1, 50];
% hsv_steps = [10, 100];


dFF_transposed = cellfun(@transpose, dFF, 'UniformOutput', false);
bhv_transposed = cellfun(@transpose, behaviorData, 'UniformOutput', false);

% because CLDS sees 80/20 test/train split
dFF_transposed = dFF_transposed(1:40);
bhv_transposed = bhv_transposed(1:40);


id_bhv = PSID(dFF_transposed, bhv_transposed, nx, n1, 12); % works with cell arrays
whichSample = 1;
[Z_recon, X_recon, r2bhv] = reconstructPSID(id_bhv.A, id_bhv.Cy, id_bhv.Cz, dFF_transposed{whichSample}, bhv_transposed{whichSample});
%%
figure()
hold on
plot(behaviorData{whichSample},'DisplayName','True bhv')
plot(X_recon.','DisplayName','PSID reconstr. bhv')
behaviorReco = Psi*B_cell{whichSample}; % b-dLDS
plot(behaviorReco.','DisplayName','b-dLDS reconstr. bhv')
legend
hold off

%%
figure()
tiledlayout('vertical')

nexttile()
plot(dFF{whichSample}(1,:).')
ylabel('True x')

nexttile()
plot(X_recon(:,1).')
ylabel('PSID reconstr. x')

nexttile()
behaviorReco = Psi*B_cell{whichSample}; % b-dLDS
plot(A_cell{whichSample}(1,:).')
ylabel('b-dLDS reconstr. x')

%%
figure();
tiledlayout(3,3)
% nexttile()
% imagesc(trueFC(:,:,750,1))
% colorbar
% clim([-1.5 1.5])
% grid off
% colormap redbluecmap
% box off
% title('Frame 750')
% ylabel('Ground truth')
% nexttile()
% imagesc(trueFC(:,:,1250,1))
% colorbar
% clim([-1.5 1.5])
% grid off
% colormap redbluecmap
% box off
% title('Frame 1250')
% nexttile()
% imagesc(trueFC(:,:,2250,1))
% colorbar
% clim([-1.5 1.5])
% grid off
% colormap redbluecmap
% box off
% title('Frame 2250')


% nexttile()
% imagesc(fcMatDLDS(:,:,750,1))
% colorbar
% clim([-1.5 1.5])
% grid off
% colormap redbluecmap
% box off
% ylabel('b-dLDS')
% nexttile()
% imagesc(fcMatDLDS(:,:,1250,1))
% colorbar
% clim([-1.5 1.5])
% grid off
% colormap redbluecmap
% box off
% nexttile()
% imagesc(fcMatDLDS(:,:,2250,1))
% colorbar
% clim([-1.5 1.5])
% grid off
% colormap redbluecmap
% box off

nexttile()
imagesc(id_bhv.A(:,:))
colorbar
clim([-1.5 1.5])
grid off
colormap redbluecmap
box off
ylabel('PSID')
nexttile()
imagesc(id_bhv.A(:,:))
colorbar
clim([-1.5 1.5])
grid off
colormap redbluecmap
box off
nexttile()
imagesc(id_bhv.A(:,:))
colorbar
clim([-1.5 1.5])
grid off
colormap redbluecmap
box off
%%
% for whichSample = 1:size(dFF,1)
%     % bhv = behaviorData{whichSample};
%     % === Train/test split ===
%     % neuralData       = dFF{whichSample};
%     % neuralDataSubset = neuralData(1:10000,:); %check if PSID works on smaller data
% 
%     % T = size(neuralData, 2);
%     % T_train = floor(0.8 * T);
%     % Y_train = neuralData(:,1:T_train).';
%     % Z_train_bhv = bhv(1,1:T_train).';
%     % Y_test = neuralData(:,T_train+1:end).';
%     % Z_test_bhv = bhv(1,T_train+1:end).';
%     % Y_test = neuralData;
%     % Z_test_bhv = bhv;
% 
%     %% === Train PSID ===
%     % id_bhv = PSID(Y_train, Z_train_bhv, nx, n1, 12);
%     id_bhv = PSID(dFF_transposed, bhv_transposed, nx, n1, 12); % works with cell arrays
% 
%     % id_hsv = PSID(Y_train, Z_train_hsv, nx, n1, 12);
% 
%     %% === Predict ===
%     % doesn't work with cell arrays - select one sample
%     [Zpb, Ztb, mse_bhv, r_bhv] = predictPSID(id_bhv.A, id_bhv.Cy, id_bhv.Cz, dFF_transposed{whichSample}, bhv_transposed{whichSample}, n_ahead);
%     % [Zph, Zth, mse_hsv, r_hsv] = predictPSID(id_hsv.A, id_hsv.Cy, id_hsv.Cz, Y_test, Z_test_hsv, n_ahead);
%     %%
%     % fig = figure('Visible', 'off', ...
%     %     'Units', 'inches', ...
%     %     'Position', [0 0 8.5 15], ... % standard portrait layout
%     %     'PaperPositionMode', 'auto');
%     fig = figure();
% 
%     % t = tiledlayout(6, 2, ...
%     %     'TileSpacing', 'tight', ...
%     %     'Padding', 'tight');
%     t = tiledlayout("vertical");
% 
%     font_size = 9;
%     lw = 0.5;
% 
%     % bhv
%     nexttile(t, 1); plot(1:n_ahead, mse_bhv, 'k-.', 'LineWidth', lw); title('bhv Rel. MSE');
%     xlabel('Time steps ahead'); ylabel('Relative MSE'); set(gca,'TickDir','out');  box off
% 
%     nexttile(t, 2); plot(1:n_ahead, r_bhv, 'b-.', 'LineWidth', lw); title('bhv Corr');
%     xlabel('Time steps ahead'); ylabel('Correlation'); set(gca,'TickDir','out'); box off
% 
%     bhv_tiles = [3, 5];
% 
%     for idx = 1:length(bhv_steps)
%         s = bhv_steps(idx);
%         if s > n_ahead
%             continue;  % skip if step exceeds prediction range
%         end
%         nexttile(t, bhv_tiles(idx), [1 2]);
%         hold on
%         box off
%         t_ax = (1:size(Zpb,1)) + s;
%         plot(t_ax, Ztb(:,s), 'b', 'LineWidth', lw); hold on;
%         plot(t_ax, Zpb(:,s), 'r-.', 'LineWidth', lw);
%         title(sprintf('bhv @%d time steps', s));
%         xlabel('Time steps'); ylabel('Behavior (a.u.)');
%         hold off
%     end
% 
%     % Add legend directly to last plot (e.g., hsv @30ms)
%     legend({'True', 'Predicted'}, ...
%         'Orientation', 'horizontal', ...
%         'Location', 'southoutside', ...
%         'Box', 'off', ...
%         'FontSize', font_size);
%      box off
% 
% 
%     % sgtitle(sprintf('Behavior predictions\nTrain: %.1f seconds | Test: %.1f seconds', ...
%     %      T_train/1000, (T - T_train)/1000), ...
%     %     'FontSize', font_size + 2);
%     set(findall(fig, '-property', 'FontSize'), 'FontSize', font_size);
% 
% 
%     % Save figure in PDF, SVG, and JPG to separate folders
%     pdf_dir = fullfile(save_dir, 'pdfFigs');
%     svg_dir = fullfile(save_dir, 'svgFigs');
%     jpg_dir = fullfile(save_dir, 'jpgFigs');
%     if ~exist(pdf_dir, 'dir'); mkdir(pdf_dir); end
%     if ~exist(svg_dir, 'dir'); mkdir(svg_dir); end
%     if ~exist(jpg_dir, 'dir'); mkdir(jpg_dir); end
% 
%     fname_pdf = fullfile(pdf_dir, sprintf('Behaviorprediction_%d', whichSample));
%     fname_svg = fullfile(svg_dir, sprintf('Behaviorprediction_%d', whichSample));
%     fname_jpg = fullfile(jpg_dir, sprintf('Behaviorprediction_%d', whichSample));
% 
%     print([fname_pdf '.pdf'], '-dpdf', '-bestfit');
%     print([fname_svg '.svg'], '-dsvg');
%     print([fname_jpg '.jpg'], '-djpeg', '-r300');  % 300 dpi for good quality
% 
%     close(fig);
% 
% end

% for cond_i = 1:length(condition_files)
%     condition_name = condition_files{cond_i};
%     disp(['Processing: ' condition_name]);
%     file_path = fullfile(base_path, [condition_name '.mat']);
%     load(file_path);  % loads Data
    % 
    % % === Neural data ===
    % % spike_rates_mat = [];
    % % for i = 1:length(cluster_numbers)
    % %     field_name = ['fr_' num2str(cluster_numbers(i))];
    % %     if isfield(Data, field_name)
    % %         spike_rates_mat = [spike_rates_mat; Data.(field_name)'];
    % %     end
    % % end
    % % dead_neurons = all(spike_rates_mat == 0, 2);
    % % spike_rates = spike_rates_mat(~dead_neurons, :)';
    % 
    % % === Behavioral signals ===
    % % bsv = Data.bsv(:);
    % % hsv = Data.hsv(:);
    % bhv = behavior_data{1};
    % 
    % % === Train/test split ===
    % T = size(dFF{1}, 1);
    % T_train = floor(0.8 * T);
    % Y_train = dFF{1}(1:T_train, :);
    % Z_train_bhv = bhv(1:T_train);
    % % Z_train_hsv = hsv(1:T_train);
    % Y_test = dFF{1}(T_train+1:end, :);
    % Z_test_bhv = bhv(T_train+1:end);
    % % Z_test_hsv = hsv(T_train+1:end);
    % 
    % % === Train PSID ===
    % id_bhv = PSID(Y_train, Z_train_bhv, nx, n1, 12);
    % % id_hsv = PSID(Y_train, Z_train_hsv, nx, n1, 12);
    % 
    % % === Predict ===
    % [Zpb, Ztb, mse_bhv, r_bhv] = predictPSID(id_bhv.A, id_bhv.Cy, id_bhv.Cz, Y_test, Z_test_bhv, n_ahead);
    % % [Zph, Zth, mse_hsv, r_hsv] = predictPSID(id_hsv.A, id_hsv.Cy, id_hsv.Cz, Y_test, Z_test_hsv, n_ahead);
    % 
    % fig = figure('Visible', 'off', ...
    %     'Units', 'inches', ...
    %     'Position', [0 0 8.5 15], ... % standard portrait layout
    %     'PaperPositionMode', 'auto');
    % 
    % t = tiledlayout(6, 2, ...
    %     'TileSpacing', 'tight', ...
    %     'Padding', 'tight');
    % 
    % font_size = 9;
    % lw = 2;
    % 
    % % bhv
    % nexttile(t, 1); plot(1:n_ahead, mse_bhv, 'k-x', 'LineWidth', lw); title('bhv Rel. MSE');
    % xlabel('Step Ahead (ms)'); ylabel('Relative MSE'); set(gca,'TickDir','out');
    % 
    % nexttile(t, 2); plot(1:n_ahead, r_bhv, 'b-x', 'LineWidth', lw); title('bhv Corr');
    % xlabel('Step Ahead (ms)'); ylabel('Correlation'); set(gca,'TickDir','out');
    % 
    % bhv_tiles = [3, 5];
    % 
    % for idx = 1:length(bhv_steps)
    %     s = bhv_steps(idx);
    %     if s > n_ahead
    %         continue;  % skip if step exceeds prediction range
    %     end
    %     nexttile(t, bhv_tiles(idx), [1 2]);
    %     t_ax = (1:size(Zpb,1)) + s;
    %     plot(t_ax/1000, Ztb(:,s), 'b', 'LineWidth', lw); hold on;
    %     plot(t_ax/1000, Zpb(:,s), 'r--', 'LineWidth', lw);
    %     title(sprintf('bhv @%dms', s));
    %     xlabel('Time (s)'); ylabel('Body velocity (deg/s)');
    % end
    % 
    % % % hsv
    % % nexttile(t, 7); plot(1:n_ahead, mse_hsv, 'k-x', 'LineWidth', lw); title('hsv Rel. MSE');
    % % xlabel('Step Ahead (ms)'); ylabel('Relative MSE'); set(gca,'TickDir','out');
    % % 
    % % nexttile(t, 8); plot(1:n_ahead, r_hsv, 'b-x', 'LineWidth', lw); title('hsv Corr');
    % % xlabel('Step Ahead (ms)'); ylabel('Correlation'); set(gca,'TickDir','out');
    % % 
    % % hsv_tiles = [9, 11];
    % % for idx = 1:length(hsv_steps)
    % %     s = hsv_steps(idx);
    % %     if s > n_ahead
    % %         continue;  % skip if step exceeds prediction range
    % %     end
    % %     nexttile(t, hsv_tiles(idx), [1 2]);
    % %     t_ax = (1:size(Zph,1)) + s;
    % %     plot(t_ax/1000, Zth(:,s), 'b', 'LineWidth', lw); hold on;
    % %     plot(t_ax/1000, Zph(:,s), 'r--', 'LineWidth', lw);
    % %     title(sprintf('hsv @%dms', s));
    % %     xlabel('Time (s)'); ylabel('Head velocity (deg/s)');
    % % end
    % 
    % 
    % % Add legend directly to last plot (e.g., hsv @30ms)
    % legend({'True', 'Predicted'}, ...
    %     'Orientation', 'horizontal', ...
    %     'Location', 'southoutside', ...
    %     'Box', 'off', ...
    %     'FontSize', font_size);
    % 
    % 
    % sgtitle(sprintf('%s - head and body velocity predictions\nTrain: %.1f seconds | Test: %.1f seconds (fs = 1kHz)', ...
    %     strrep(condition_name, '_', ' '), T_train/1000, (T - T_train)/1000), ...
    %     'FontSize', font_size + 2);
    % set(findall(fig, '-property', 'FontSize'), 'FontSize', font_size);
    % % Save figure in PDF, SVG, and JPG to separate folders
    % pdf_dir = fullfile(save_dir, 'pdfFigs');
    % svg_dir = fullfile(save_dir, 'svgFigs');
    % jpg_dir = fullfile(save_dir, 'jpgFigs');
    % if ~exist(pdf_dir, 'dir'); mkdir(pdf_dir); end
    % if ~exist(svg_dir, 'dir'); mkdir(svg_dir); end
    % if ~exist(jpg_dir, 'dir'); mkdir(jpg_dir); end
    % 
    % fname_pdf = fullfile(pdf_dir, [condition_name '_prediction']);
    % fname_svg = fullfile(svg_dir, [condition_name '_prediction']);
    % fname_jpg = fullfile(jpg_dir, [condition_name '_prediction']);
    % 
    % print([fname_pdf '.pdf'], '-dpdf', '-bestfit');
    % print([fname_svg '.svg'], '-dsvg');
    % print([fname_jpg '.jpg'], '-djpeg', '-r300');  % 300 dpi for good quality
    % 
    % close(fig);
% end
end
%% === Reconstruction Function ===
function [Z_recon, X_recon, r2] = reconstructPSID(A, C_y, C_z, Y, Z)

T_test = size(Y, 1);
valid_T = T_test - 1;
Z_recon = nan(valid_T, size(Z,2));
X_recon = nan(valid_T, size(Y,2));
for k = 1:valid_T
    y_k     = Y(k,:)';
    x_t     = pinv(C_y) * y_k;
    x_t     = A * x_t;
    X_recon(k,:) = x_t;
    Z_recon(k,:) = C_z * x_t;
end
corrRecon = corrcoef(Z(2:end,:),Z_recon(:));
r2        = corrRecon(1,2)^2;
end

%% === Predict Function ===
function [Z_pred, Z_true, rel_mse, corrs] = predictPSID(A, C_y, C_z, Y, Z, n_ahead)

T_test = size(Y, 1);
valid_T = T_test - n_ahead;
Z_pred = nan(valid_T, n_ahead);
Z_true = nan(valid_T, n_ahead);
for k = 1:valid_T
    y_k = Y(k,:)';
    x_t = pinv(C_y) * y_k;
    for i = 1:n_ahead
        x_t = A * x_t;
        Z_pred(k, i) = C_z * x_t;
    end
    Z_true(k, :) = Z(k+1:k+n_ahead)';
end
rel_mse = sum((Z_pred - Z_true).^2, 1) ./ sum(Z_true.^2, 1);
corrs = arrayfun(@(i) corr(Z_pred(:,i), Z_true(:,i), 'rows','complete'), 1:n_ahead);
end
