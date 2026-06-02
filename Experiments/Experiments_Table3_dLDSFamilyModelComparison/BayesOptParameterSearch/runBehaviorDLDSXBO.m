function objectiveToMinimize = runBehaviorDLDSXBO(dFF,behavior_data,...
    lambda_history,lambda_behavior,inf_opts)

disp('optimize for behavior R2')
disp('minimize 1-R^2')

% inf_opts.lambda_D = 0; % in main script
inf_opts.lambda_behavior = lambda_behavior;
inf_opts.lambda_history  = lambda_history;
%%
tic

% https://www.mathworks.com/matlabcentral/fileexchange/51945-memorylinux
[~,meminfo] = system('cat /proc/meminfo'); 
tokens = regexpi(meminfo,'^MemTotal:\s*(\d+)\s', 'tokens'); 
totalmem = str2double(tokens{1}{1});  
 % get available memory                                                  
tokens = regexpi(meminfo,'^*MemFree:\s*(\d+)\s','tokens');    
freemem = str2double(tokens{1}{1});                 
startmem = totalmem-freemem;    

currentState = rng('shuffle');

if inf_opts.behaviordLDS
    [Phi, F, Psi] = bpdndf_dynamics_learning_behavior_x(dFF, [], [], [], behavior_data, inf_opts);              % Run the dictionary learning algorithm
else
    % [Phi, F] = bpdndf_dynamics_learning(dFF, [], [], inf_opts);              % Run the dictionary learning algorithm
end

if inf_opts.behaviordLDS
    if inf_opts.AcrossIndividuals % EY added 08/27/2023
        for ii = 1:size(data_obj,1)
            [A_cell{ii},B_cell{ii}] = parallel_bilinear_dynamic_inference_behavior(dFF(ii), Phi{ii}, F, Psi{ii}, behaviorData{ii}, ...
                                       @bpdndf_bilinear_handle_behavior_x, inf_opts); % Infer sparse coefficients
        end
    else
        [A_cell,B_cell] = parallel_bilinear_dynamic_inference_behavior(dFF, Phi, F, Psi, behavior_data, ...
                                       @bpdndf_bilinear_handle_behavior_x, inf_opts); % Infer sparse coefficients
    end
else
    % if inf_opts.AcrossIndividuals
    %     for ii = 1:size(Phi,1)
    %         [A_cell{ii},B_cell{ii}] = parallel_bilinear_dynamic_inference(dFF(ii), Phi{ii}, F, ...
    %                                        @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients
    %     end
    % else
    %     [A_cell,B_cell] = parallel_bilinear_dynamic_inference(dFF, Phi, F, ...
    %                                        @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients
    % 
    % end
end


dLDSruntime = toc;

fprintf('\nElapsed time: %.2f seconds\n', dLDSruntime);

% https://www.mathworks.com/matlabcentral/fileexchange/51945-memorylinux
[~,meminfo] = system('cat /proc/meminfo'); 
tokens = regexpi(meminfo,'^MemTotal:\s*(\d+)\s', 'tokens'); 
totalmem = str2double(tokens{1}{1});  
 % get available memory                                                  
tokens = regexpi(meminfo,'^*MemFree:\s*(\d+)\s','tokens');    
freemem = str2double(tokens{1}{1});                 
endmem = totalmem-freemem; 

memused = endmem-startmem;
fprintf('\nUsed memory: %0.2f\n', memused)

reconstructedBehavior = Psi*A_cell{1};
allbhv = behavior_data{1};
rval = corrcoef(allbhv(:), reconstructedBehavior(:));
varExplBhv = (rval(1,2))^2;
objectiveToMinimize = 1 - varExplBhv;

% varExpl = overallVarExpl(dFF,A_cell,Phi,inf_opts);
%%
% varExplPerNeuron = perNeuronVarExpl(dFF,A_cell,Phi,inf_opts);
%%