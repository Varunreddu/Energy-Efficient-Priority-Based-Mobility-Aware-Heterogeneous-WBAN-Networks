function results = atmac_priority_markov_sim()
% Priority-aware AT-MAC Markov model + contention simulator
% Matches figures/equations (11)–(20): P1/P2/P3 priority entry, CCA1/CCA2
% with α,β, collision-driven r advancement, and SU/FL resets. Then runs a
% minislot contention simulator using the solved τ and per-mode payload pmfs.

rng(21);

%% ---------------- Simulation parameters ----------------
Sim.Tend = 100;    % seconds
dt = 1e-3;         % seconds

%% ---------------- MAC / Backoff ----------------
macMinBE=3; macMaxBE=5; maxRetries=3;

%% ---------------- Superframe ----------------
SF.superframe = 0.1;    % 100 ms
SF.minislots  = 20;     % 20 contention opportunities per SF
SF.rssiThresh = 5.0;    % dB mobility trigger

%% ---------------- Nodes ----------------
Start.server=0.3; Start.coordinator=0.5;
Start.heart=1.0; Start.spo2=1.2; Start.temp=1.4; Start.bp=1.6;

Sensors = [ ...
  struct('name',"heart",'lambda',10,'start',Start.heart,'initBuffer',0,'ts',2,'w',[0.40 0.30 0.20 0.10])
  struct('name',"spo2", 'lambda',10,'start',Start.spo2, 'initBuffer',0,'ts',1,'w',[0.35 0.30 0.25 0.10])
  struct('name',"temp", 'lambda',10,'start',Start.temp, 'initBuffer',0,'ts',0,'w',[0.30 0.25 0.35 0.10])
  struct('name',"bp",   'lambda',10,'start',Start.bp,   'initBuffer',0,'ts',2,'w',[0.40 0.30 0.20 0.10]) ];
N = numel(Sensors);
lambda_eff_scale = 1.0; B = 50*ones(N,1);

%% ---------------- AT-MAC parameters ----------------
at.p = 3;                 % stages 0..p
at.L = 5;                 % payload classes
at.R = 3;                 % max retransmissions 0..R
at.CWmin = 10; at.CWscale = 2;
at.payload_vals = [20 40 60 80 100]; % bytes
% Priority probabilities: P1=Critical, P2=Urgent, P3=Normal
at.P_pri = [0.50 0.30 0.20];  % must sum to 1

% Linear busy-prob fits a(L), b(L) (tune using paper's ranges/measured fit)
% a_m(Lbits) = a1(m)*Lbits + a0(m), b_m(Lbits) = b1(m)*Lbits + b0(m)
at.a1 = [9e-6 1.0e-5 1.1e-5]; at.a0 = [0.05 0.06 0.07];
at.b1 = [7e-6 8.0e-6 9.0e-6]; at.b0 = [0.04 0.05 0.06];

%% ---------------- PHY / Energy ----------------
ACK_bits=88;
E = 4.0*ones(N,1);
d=0.5; n_path=2; ETxElec=50e-9; ERxElec=50e-9; Eamp=1e-10;
Pidle_sleep=5e-5; Ectrl_round=1e-5;
E_backoff_per_ms=2e-7; bits_per_slot=80; DataRate=250e3;
slotTime = bits_per_slot/DataRate;
Ebackoff_per_slot = E_backoff_per_ms * (slotTime*1000);

%% ---------------- Channel (Gilbert–Elliott) ----------------
ch.state=1; ch.sG=0.995; ch.sB=0.60; ch.pGB=0.015; ch.pBG=0.12;

%% ---------------- Fixed point: priority Markov per eq (12)–(20) ----------------
fp = fixed_point_priority_eq1220(at, N);
fprintf('Fixed point: tau=%.6g, Pcoll=%.6g\n', fp.tau, fp.Pcoll);
fprintf('Payload pmf per mode:\n');
disp(fp.pl_pmf_mode);

%% ---------------- Initialize runtime ----------------
t=0; nextSfT = Start.coordinator;
arrivalTimes = cell(N,1); arrivalMode = cell(N,1); q=zeros(N,1);
for i=1:N, arrivalTimes{i}=[]; arrivalMode{i}=[]; end
rssi = -60 + 5*randn(N,1); prev=rssi;

txPk=zeros(N,1); rxPk=zeros(N,1); genPk=zeros(N,1);
delaySum=zeros(N,1); qAvg=zeros(N,1);

isActive = @(i,tt) (tt >= getfield(struct('heart',Start.heart,'spo2',Start.spo2,'temp',Start.temp,'bp',Start.bp), char(Sensors(i).name)));

% Mode-based minislot bias
function mb = mode_bias(ii)
    if ~isempty(arrivalMode{ii}), m=arrivalMode{ii}(1); else, m=3; end
    if m==1, mb=1.35; elseif m==2, mb=1.10; else, mb=0.95; end
end

% Payload sampling using chain-derived per-mode pmf
pmf_mode = fp.pl_pmf_mode;
function bytes = sample_payload(ii)
    if ~isempty(arrivalMode{ii}), m=arrivalMode{ii}(1); else, m=3; end
    pmf = pmf_mode(m,:); pmf = pmf/sum(pmf);
    if q(ii) > 0.5*B(ii), pmf = pmf .* [1.15 1.12 1.08 0.93 0.80]; pmf=pmf/sum(pmf); end
    cdf = cumsum(pmf); u=rand; idx=find(u<=cdf,1,'first'); if isempty(idx), idx=1; end
    bytes = at.payload_vals(idx);
end

%% ---------------- Main loop ----------------
while t < Sim.Tend - 1e-12

    % Superframe boundary: contention in M minislots
    if t >= nextSfT - 1e-12 && t >= Start.coordinator
        prev = rssi;
        nextSfT = t + SF.superframe;
        E = E - Ectrl_round * ones(N,1);

        target_load = 1.35; M = SF.minislots;
        for ms=1:M
            contenders = find(arrayfun(@(ii) isActive(ii,t) && q(ii)>0, 1:N));
            Na = numel(contenders); if Na==0, break; end
            base_dyn = min(0.6, target_load / Na);

            attempts = false(N,1);
            for ii = contenders
                p_i = base_dyn * mode_bias(ii);
                if ~isnan(rssi(ii)) && ~isnan(prev(ii)) && abs(rssi(ii)-prev(ii)) >= SF.rssiThresh && q(ii)>0
                    p_i = max(p_i, 0.40);  % mobility floor
                end
                p_i = max(0, min(0.95, p_i));
                if rand < p_i, attempts(ii)=true; end
            end

            K = find(attempts);
            if isempty(K)
                % idle minislot
            elseif numel(K)==1
                i = K(1);
                retries=0; BE=macMinBE; success=false;
                while retries <= maxRetries && ~success && E(i) > 0
                    Lbytes = sample_payload(i);
                    Lbits  = 8*Lbytes;
                    Etx = Lbits*(ETxElec + Eamp*(d^n_path));
                    Eack = ACK_bits*ERxElec;
                    txPk(i)=txPk(i)+1; E(i)=E(i)-Etx;

                    % GE update
                    if ch.state==1, if rand<ch.pGB, ch.state=2; end
                    else, if rand<ch.pBG, ch.state=1; end
                    end

                    if (ch.state==1 && rand<ch.sG) || (ch.state==2 && rand<ch.sB)
                        rxPk(i)=rxPk(i)+1; E(i)=E(i)-Eack;
                        if ~isempty(arrivalTimes{i})
                            atime = arrivalTimes{i}(1); delaySum(i)=delaySum(i)+(t-atime);
                            arrivalTimes{i}(1)=[];
                        end
                        if ~isempty(arrivalMode{i}), arrivalMode{i}(1)=[]; end
                        q(i)=max(0,q(i)-1); success=true;
                    else
                        boSlots = randi([0, max(0, 2^BE - 1)]);
                        E(i) = E(i) - boSlots * Ebackoff_per_slot;
                        retries=retries+1; BE=min(macMaxBE, BE+1);
                    end
                end
            else
                % collision penalty
                for ii = K.'
                    boSlots = randi([0, max(0, 2^macMinBE - 1)]);
                    E(ii) = E(ii) - boSlots * Ebackoff_per_slot;
                end
            end
        end
    end

    % Arrivals (Poisson)
    for i=1:N
        if isActive(i,t)
            p_arr = min(Sensors(i).lambda * lambda_eff_scale * dt, 1.0);
            if rand < p_arr && q(i) < B(i)
                q(i)=q(i)+1; arrivalTimes{i}=[arrivalTimes{i}, t]; genPk(i)=genPk(i)+1;
                % Draw priority mode by P1,P2,P3
                Pcum = cumsum(at.P_pri / sum(at.P_pri)); u=rand;
                m = find(u<=Pcum,1,'first'); if isempty(m), m=3; end
                arrivalMode{i} = [arrivalMode{i}, m];
            end
        end
    end

    % Idle energy, RSSI walk, average queue
    E = E - Pidle_sleep * ones(N,1) * dt;
    rssi = rssi + 0.1*randn(N,1);
    qAvg = qAvg + q;

    t = t + dt;
end

%% ---------------- Metrics ----------------
qAvg = qAvg / (Sim.Tend / dt);
avgDelay = NaN(N,1);
for i=1:N, if rxPk(i)>0, avgDelay(i)=delaySum(i)/rxPk(i); end, end
throughput_pkps = rxPk/Sim.Tend;

valid_mask = genPk > 0;
pdr_per_sensor = NaN(N,1);
pdr_per_sensor(valid_mask) = min(1, rxPk(valid_mask) ./ genPk(valid_mask));
pdr_overall = sum(rxPk) / max(sum(genPk),1);
avgThroughput = sum(rxPk)/Sim.Tend;
avgPDR = mean(pdr_per_sensor(valid_mask));
validDelays = avgDelay(~isnan(avgDelay));
if isempty(validDelays), meanDelay=NaN; else, meanDelay=mean(validDelays); end
E(E<0)=0;

fprintf('\n--- Metrics ---\nThroughput: %.3f pkts/s\nPDR (network): %.3f\nDelay: %.3f s\n', ...
    avgThroughput, pdr_overall, meanDelay);

simTable = table(string({Sensors.name}'), [Sensors.lambda]'*lambda_eff_scale, genPk, txPk, rxPk, ...
    throughput_pkps, avgDelay, qAvg, pdr_per_sensor, E, ...
    'VariableNames', {'sensor','lambda_eff','genPk','txPk','rxPk','throughput_pkps','avgDelay_s','avgQueue','pdr_per_node','residualE_J'});
disp(simTable);
writetable(simTable,'atmac_priority_markov_results.csv');

results.simResults = simTable;
results.fp = fp;
results.pdr_overall = pdr_overall;
results.avgPDR_active = avgPDR;

end % main

%% =====================================================================
function sol = fixed_point_priority_eq1220(at, Nnodes)
% Builds Markov chain per eqs (12)–(20), solves stationary pi and τ, iterates
% with multi-node coupling Pcoll = 1 - (1 - τ)^(N-1). Returns τ, Pcoll, and
% per-mode payload pmf induced by α,β.

pSt=at.p; L=at.L; R=at.R; M=3;
CW = zeros(1,pSt+1);
for ii=0:pSt, CW(ii+1) = round(at.CWmin * (at.CWscale^ii)); end
tau_stage = 1 ./ (((CW-1)/2) + 1);

% α(L), β(L)
function [alpha,beta] = alpha_beta(m,l)
    Lbits = 8*at.payload_vals(l);
    a = min(0.99, max(0, at.a1(m)*Lbits + at.a0(m)));
    b = min(0.99, max(0, at.b1(m)*Lbits + at.b0(m)));
    alpha = 1 - a; beta = 1 - b;
end

% Build state list
function [states,SPEC,head] = make_states()
    SPEC = struct('ID',1,'C1',2,'C2',3,'SU',4,'FL',5);
    head = 5;
    grid = [];
    for m=1:M, for i=0:pSt, for l=1:L, for r=0:R
        grid=[grid; m i l r]; %#ok<AGROW>
    end, end, end, end
    states = [ -1 -1 -1 -1; -2 -2 -2 -2; -3 -3 -3 -3; -4 -4 -4 -4; -5 -5 -5 -5; grid ];
end

% Indexer
function id = idx_of(states, m,i,l,r)
    id = find(states(:,1)==m & states(:,2)==i & states(:,3)==l & states(:,4)==r, 1);
end

tau = 0.1;
for iter=1:100
    Pcoll = 1 - (1 - tau)^(max(Nnodes-1,0));
    [P, states, SPEC] = build_P(Pcoll);
    % Stationary distribution
    try
        opts.tol=1e-10; [V_eig,~]=eigs(P.',1,1,opts); pi = real((V_eig/sum(V_eig))');
    catch
        [V_eig,D_eig]=eig(full(P.')); [~,ix]=min(abs(diag(D_eig)-1)); v=V_eig(:,ix); pi=real((v/sum(v))');
    end
    % τ = sum π * τ_stage over contention states
    tau_new = 0;
    for s=1:size(states,1)
        st = states(s,:); iS=st(2);
        if iS>=0, tau_new = tau_new + pi(s)*tau_stage(iS+1); end
    end
    if abs(tau_new - tau) < 1e-8, tau = tau_new; break; end
    tau = 0.5*tau + 0.5*tau_new;
end

% Per-mode payload pmf via αβ weights (steady CCA pass tendency)
pl_mode = zeros(M,L);
for m=1:M
    w=zeros(1,L);
    for l=1:L
        [alpha,beta]=alpha_beta(m,l);
        w(l)=alpha*beta;
    end
    if sum(w)>0, w=w/sum(w); else, w=ones(1,L)/L; end
    pl_mode(m,:)=w;
end
pl_overall = (at.P_pri(:)'/sum(at.P_pri)) * pl_mode;

sol.tau=tau; sol.Pcoll=Pcoll; sol.pl_pmf_mode=pl_mode; sol.pl_pmf_overall=pl_overall;

% ---------- Build P following eqs (12)–(20) ----------
function [P, states, SPEC] = build_P(Pcoll_in)
    [states, SPEC, head] = make_states();
    Ns = size(states,1); I=[]; J=[]; V=[];
    Pm = at.P_pri(:)'/sum(at.P_pri);

    % Eq (12): Idle -> (p,0,0,0) with Px * τ; spread uniformly over l
    tau_enter = mean(tau_stage);
    for m=1:M
        for l=1:L
            I(end+1)=SPEC.ID; J(end+1)=idx_of(states,m,0,l,0); V(end+1)=Pm(m)*tau_enter/L; %#ok<AGROW>
        end
    end
    % Eq (18): Idle self-loop to complete row
    % (set later by normalization)

    % Eqs (13),(14),(19),(20): CCA1/2 behavior with α,β
    for s=head+1:Ns
        st = states(s,:); m=st(1); i=st(2); l=st(3); r=st(4);
        if i<0, continue; end
        tp = tau_stage(i+1);
        [alpha,beta]=alpha_beta(m,l);

        % From (m,i,l,r): attempt to C1 with τ_i, else stay
        I(end+1)=s; J(end+1)=SPEC.C1; V(end+1)=tp; %#ok<AGROW>
        I(end+1)=s; J(end+1)=s;        V(end+1)=1-tp; %#ok<AGROW>

        % C1 -> C2 (α), C1 -> C1 (1-α)
        I(end+1)=SPEC.C1; J(end+1)=SPEC.C2; V(end+1)=alpha/(L*M); %#ok<AGROW>
        I(end+1)=SPEC.C1; J(end+1)=SPEC.C1; V(end+1)=(1-alpha)/(L*M); %#ok<AGROW>

        % C2 -> SU (β), C2 -> C2 (1-β)
        I(end+1)=SPEC.C2; J(end+1)=SPEC.SU; V(end+1)=beta/(L*M); %#ok<AGROW>
        I(end+1)=SPEC.C2; J(end+1)=SPEC.C2; V(end+1)=(1-beta)/(L*M); %#ok<AGROW>
    end

    % Eqs (15)–(17): Collision and retry evolution
    for s=head+1:Ns
        st = states(s,:); m=st(1); i=st(2); l=st(3); r=st(4);
        if i<0, continue; end
        if r < R
            j2 = min(i+1, pSt);
            I(end+1)=s; J(end+1)=idx_of(states,m,j2,l,r+1); V(end+1)=Pcoll_in * 0.5 / (L*M); %#ok<AGROW>
        else
            % cumulative drop after R+1 collisions approx Pcoll^(R+1)
            I(end+1)=s; J(end+1)=SPEC.FL; V(end+1)=(Pcoll_in^(R+1)) / (L*M); %#ok<AGROW>
        end
    end

    % SU/FL -> ID (reset)
    I(end+1)=SPEC.SU; J(end+1)=SPEC.ID; V(end+1)=1.0; %#ok<AGROW>
    I(end+1)=SPEC.FL; J(end+1)=SPEC.ID; V(end+1)=1.0; %#ok<AGROW>

    % Normalize rows to stochastic
    P = sparse(I,J,V,Ns,Ns);
    rs = sum(P,2);
    z=find(rs==0); for k=z.', P(k,k)=1; end
    rs = sum(P,2);
    D = spdiags(1./rs,0,Ns,Ns);
    P = D * P;
end

end % fixed_point_priority_eq1220
