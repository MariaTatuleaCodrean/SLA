function [x0,t0,tend,tsim,vel,len,way,active,alpha,nbact] = bactchainSimulation_lengthdistrib(n0,X0,vel_pdf,len_pdf,outputfile,Nsim,Nevent,Pright)
%BACTCHAINSIMULATION Simulation of single lane swimming to predict chain
%length and velocity evolution for bacterial aggregation in liquid crystals
%
% Input arg     Type            Description
%   n0          int             initial number of bacteria
%   X0          double array    initial positions of bacteria
%   vel_pdf     distrib         bacterial speed distribution
%   len_pdf     distrib         bacterial length distribution
%
% Optional arguments
%   outputfile  str             name of outputfile
%   Nsim        int             number of simulation repeats
%   Nevent      int             number of merging events
%   Pright      double          probability of moving right
%   
% Output arg
%   x0      double array        start position ...
%   t0      double array        start time of chain trajectory
%   tend    double array        end time ...
%   tsim    double array        time of merger events
%   vel     double array        speed (constant for duration of trajectory)
%   len     double array        absolute chain length / drag
%   way     double array        direction: +1 for right, -1 for left
%   active  logical array       flag for whether chain is active or not
%   alpha   double array        active fraction of agents
%   nbact   double array        chain length as number of bacteria
%
% Note: each column in output variables is an independent trajectory, so
% the index along the second dimension of the array is the unique ID of
% each bacterium or bacterial chain observed during the simulation

% Choose default input values for optional arguments, if needed
if ~exist('outputfile','var') || isempty(outputfile)
    outputfile='output_autosave.mat';
    save(outputfile,'n0','X0','vel_pdf')
end
if ~exist('Nsim','var') || isempty(Nsim)
    Nsim=1;
end
if ~exist('Nevent','var') || isempty (Nevent)
    Nevent=ceil(.8*n0); % default to 80% aggregation
end
if ~exist('Pright','var') || isempty (Pright)
    Pright=0.5; % default to equal left/right probability
end

% Define chain variables
t0     = nan(Nsim,n0+Nevent);  % start time of chain trajectory
tend   = nan(Nsim,n0+Nevent);  % end time ...
x0     = nan(Nsim,n0+Nevent);  % start position ...
vel    = nan(Nsim,n0+Nevent);  % speed (constant for duration of trajectory)
len    = nan(Nsim,n0+Nevent);  % chain length/drag 
way    = nan(Nsim,n0+Nevent);  % direction: +1 for right, -1 for left
active = false(Nsim,n0+Nevent); % flag for whether chain is active or not
nbact  = nan(Nsim,n0+Nevent);  % chain length as number of bacteria

% Define time vector
tsim  = nan(Nsim,Nevent+1);
alpha = nan(Nsim,Nevent+1);

% Start multiple independent simulations
for jj = 1:Nsim
    way(jj,1:n0) = 2*(rand(1, n0) > (1-Pright)) - 1; % swimming direction

    v = random(vel_pdf, n0, 1);  % experimentally measured distribution
    vel(jj,1:n0) = v;            % Initialize velocities

    L = random(len_pdf, n0, 1);  % cell length distribution
    len(jj,1:n0)   = L;          % Initialize lengths
    nbact(jj,1:n0) = 1;          % Initialize number of bacteria

    x0(jj,1:n0) = X0;            % Initialize positions
    active(jj,1:n0) = true;      % Activate swimming bacteria

    % Start time
    t0(jj,1:n0) = 0;
    tsim(jj,1) = 0;

    % Iterate over finite number of merging events
    for kk=1:Nevent
        % Current position of all bacteria
        x = x0(jj,:) + (tsim(jj,kk)-t0(jj,:)).*vel(jj,:).*way(jj,:);

        % Sort from left to right (NaN values always remain at the end)
        [x,I]   = sort(x);
        x0(jj,:)     = x0(jj,I);
        t0(jj,:)     = t0(jj,I);
        tend(jj,:)   = tend(jj,I);
        vel(jj,:)    = vel(jj,I);
        len(jj,:)    = len(jj,I);
        way(jj,:)    = way(jj,I);
        active(jj,:) = active(jj,I);
        nbact(jj,:)  = nbact(jj,I);

        % Find active chains
        idx = find(active(jj,:));

        % Compute relative velocities and positions of active chains only
        [x_selected, v_selected] = deal(x(idx), vel(jj,idx).*way(jj,idx));
        dx = diff(x_selected);  % dx = x2 - x1, should be all positive
        dv = -diff(v_selected); % minus means dv = v1 - v2, can #1 catch up?
        if min(dx)<0
            error('Error. Ordering went wrong.')
        end

        % Find next merging time
        nextdt = dx./dv;
        nextdt(nextdt<0) = Inf; % discard negative timestepping values
        [dt, Istar] = min(nextdt);

        % If all times are infinite, no more merging events
        if dt==Inf
            disp(['Simulation ' num2str(jj) '/' num2str(Nsim) ' finished without complete merging'])
            tsim(jj,kk+1) = tsim(jj,kk);
        else
            % Identify chains involved in merger event
            i1 = idx(Istar); % from active ones, select the one with lowest merging time
            i2 = idx(Istar+1);

            % Deactivate bacteria that are merging
            [active(jj,i1), active(jj,i2)] = deal(false, false);
            [tend(jj,i1), tend(jj,i2)]     = deal(tsim(jj,kk) + dt, tsim(jj,kk) + dt);

            % Activate new index
            active(jj,n0+kk) = true;

            % Determine new parameters
            x0(jj,n0+kk) = x(i1) + dt*vel(jj,i1)*way(jj,i1);
            if ((x(i1) + dt*vel(jj,i1)*way(jj,i1))-(x(i2) + dt*vel(jj,i2)*way(jj,i2)))>1e-14
                disp((x(i1) + dt*vel(jj,i1)*way(jj,i1))-(x(i2) + dt*vel(jj,i2)*way(jj,i2)))
                warning('Warning. Merger may be faulty.')
            end
            t0(jj,n0+kk)    = tsim(jj,kk)  + dt;
            len(jj,n0+kk)   = len(jj,i1)   + len(jj,i2);
            nbact(jj,n0+kk) = nbact(jj,i1) + nbact(jj,i2);

            % Determine force and new direction
            idx = [i1 i2];
            F = [len(jj,i1)*vel(jj,i1)   len(jj,i2)*vel(jj,i2)];
            [~,Istar] = max(F); % select agent with larger thrust force
            way(jj,n0+kk) = way(jj,idx(Istar));
            vel(jj,n0+kk) = (vel(jj,i1)*len(jj,i1)+vel(jj,i2)*len(jj,i2))/len(jj,n0+kk);

            % Update simulation time for next time step
            tsim(jj,kk+1) = tsim(jj,kk) + dt;
        end

        % Caculate active fraction at each time point (trivial) 
        alpha(jj,:) = 1 - (0:1:Nevent)/n0;
    end

    % Re-order chains by start time and start position
    x0temp = x0(jj,:);
    [~,I]  = sort(x0temp);
    x0(jj,:)     = x0(jj,I);
    t0(jj,:)     = t0(jj,I);
    tend(jj,:)   = tend(jj,I);
    vel(jj,:)    = vel(jj,I);
    len(jj,:)    = len(jj,I);
    way(jj,:)    = way(jj,I);
    active(jj,:) = active(jj,I);
    nbact(jj,:)  = nbact(jj,I);

    t0temp = t0(jj,:);
    [~,I]  = sort(t0temp);
    x0(jj,:)     = x0(jj,I);
    t0(jj,:)     = t0(jj,I);
    tend(jj,:)   = tend(jj,I);
    vel(jj,:)    = vel(jj,I);
    len(jj,:)    = len(jj,I);
    way(jj,:)    = way(jj,I);
    active(jj,:) = active(jj,I);
    nbact(jj,:)  = nbact(jj,I);

end

% Save input and output variables
save(outputfile,'n0','X0','vel_pdf','len_pdf','Nsim','Nevent','x0','t0',...
    'tend','tsim','vel','len','way','active','alpha','nbact','-append')
end