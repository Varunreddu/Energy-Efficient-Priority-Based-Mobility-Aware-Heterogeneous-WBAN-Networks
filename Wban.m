function Wban_CSMA_TDMA_Optimized()
% WBAN TDMA + CSMA/CA (Table-1) + 4 J battery.
% Optimized: binary min-heap event queue, array counters.
% Extended with averages + graphs: throughput, PDR, residual energy, reliability, delay.

%% ---------------- Parameters ----------------
SIM_TIME      = 10.0;            % seconds
SUPERFRAME    = 0.8;               % TDMA cycle
RSSI_THRESH   = 5.0;               % dB mobility trigger
MOBILE_SLOT   = 0.035;             % seconds
EXPECTED      = 4; REG_TIMEOUT = 0.30;   % seconds

% Table 1 parameters
CSMA_MAX_BACKOFFS = 4;
MAC_MAX_RETRIES   = 3;
MAX_PAYLOAD_B     = 127;
macMinBE          = 3;
macMaxBE          = 5;
BITS_PER_SLOT_UNIT= 80;
PHY_RATE_bps      = 250e3;
SYMBOL_RATE       = 62.5e3;
ACK_BITS          = 88;

% Framing and timing
MAC_HDR_B=10; PHY_HDR_B=6; PAYLOAD_B=MAX_PAYLOAD_B;
SIFS=1e-3; SLOT_GUARD=2e-3; EXTRA_GUARD=2e-4;

% Channel/link
FREQ=2.45e9; C=3e8; lambda=C/FREQ; PT_mW=2.24; PT_dBm=10*log10(PT_mW);
NOISE_dBm=-106; SNIR_T_dB=-8; SENS_dBm=-92; %#ok<NASGU>
PL_shadow_sigma=4; MIN_DIST=0.1;

% Starts and names
HEART_START=0.30; SPO2_START=0.35; TEMP_START=0.40; BP_START=0.45;
HEART_PORT=5017; SPO2_PORT=5018; TEMP_PORT=5019; BP_PORT=5020;

%% ---------------- 4 J Energy Model ----------------
E0_J = 4.0; Vbat = 3.0;
I_TX=8e-3; I_RX=7e-3; I_IDLE=0.12e-3;

E_TX  = @(t) Vbat*I_TX*t;
E_RX  = @(t) Vbat*I_RX*t;
E_IDLE= @(t) Vbat*I_IDLE*t;

% Timing from rates
txBits = 8*(MAC_HDR_B+PHY_HDR_B+PAYLOAD_B);
t_tx_data = txBits/PHY_RATE_bps;
t_tx_ack  = ACK_BITS/PHY_RATE_bps;
t_cca     = 20/SYMBOL_RATE;

% Idle per superframe (tapered) and reserve floor
E_idle_sf_J_base = E_IDLE(0.004);
IDLE_TAPER = 0.75;                % 25% reduction in modeled idle drain
E_idle_sf_J = IDLE_TAPER * E_idle_sf_J_base;
E_RESERVE   = 0.20;               % keep at least 0.2 J in battery
E_SCALE     = 0.92;               % modest 8% savings from duty-cycling

% Contention/success tuning
CCA_BUSY_PROB = 0.22;
HIGH_SINR_CEIL = 0.975;
LOSS_FLOOR = 0.3;

%% ---------------- Nodes ----------------
coordPos=[0,0,1.6];
sensors = {
  makeSensor('heart',HEART_PORT,HEART_START,10.0,2,[0.1,0.0,1.2])
  makeSensor('spo2', SPO2_PORT,SPO2_START,10.0,2,[-0.1,0.0,1.2])
  makeSensor('temp', TEMP_PORT,TEMP_START,10.0,1,[0.0,0.1,1.2])
  makeSensor('bp',   BP_PORT,  BP_START, 10.0,0,[0.0,-0.1,1.2])
};
N = numel(sensors);
names = cell(1,N); for i=1:N, names{i}=sensors{i}.name; end
name2idx = containers.Map(names, num2cell(1:N));

for i=1:N
  s=sensors{i};
  s.buffer=5; s.rssi=NaN; s.prevRssi=NaN; s.bound=false;
  s.energy_J=E0_J;
  s.mac.queue=cell(0,1); s.mac.txAttempts=0; s.mac.retransmissions=0; s.mac.packetsLost=0;
  s.mac.packetsSent=0; s.mac.awaitingAck=false; s.mac.expectedAckFor=[]; s.mac.nextPktId=1;

  % PDF-aligned fields
  s.T = s.time;     % interpret given 'time' as time-sensitivity tier T in {0,1,2}
  s.DR = s.rate;    % data rate proxy used for drS
  s.PS = 0.0;       % priority score computed each schedule

  sensors{i}=s;
end

% Counters (arrays)
genCntArr = zeros(1,N); rcvCntArr = zeros(1,N);
totalDelay = zeros(1,N); rcvDelayCnt = zeros(1,N);

% Metrics
dt=1.0; nextSample=dt; now=0;
metrics.t=[]; metrics.thr_total=[]; metrics.pdr_total=[];
metrics.energy=containers.Map; metrics.per=containers.Map;
for i=1:N, metrics.energy(names{i})=[]; metrics.per(names{i})=struct('thr',[],'pdr',[]); end

%% ---------------- Fast event queue (min-heap) ----------------
EV_CAP = 200000;
EV = repmat(struct('t',0,'type','','who','','data',struct()), EV_CAP, 1);
evCount = 0; heap = zeros(EV_CAP,1); heapSize = 0;

  function heapSwap(i,j), tmp=heap(i); heap(i)=heap(j); heap(j)=tmp; end
  function heapUp(k)
    while k>1
      p=floor(k/2);
      if EV(heap(p)).t <= EV(heap(k)).t, break; end
      heapSwap(p,k); k=p;
    end
  end
  function heapDown(k)
    while true
      l=2*k; r=l+1; s=k;
      if l<=heapSize && EV(heap(l)).t < EV(heap(s)).t, s=l; end
      if r<=heapSize && EV(heap(r)).t < EV(heap(s)).t, s=r; end
      if s==k, break; end
      heapSwap(k,s); k=s;
    end
  end
  function scheduleAt_local(type,t,who,data)
    if nargin<3, who=''; end; if nargin<4, data=struct(); end
    evCount = evCount + 1;
    EV(evCount).t = t; EV(evCount).type=type; EV(evCount).who=who; EV(evCount).data=data;
    heapSize = heapSize + 1; heap(heapSize) = evCount; heapUp(heapSize);
  end
  function [t,ev,ok]=pop_local()
    if heapSize==0, t=0; ev=struct(); ok=false; return; end
    idx = heap(1); heap(1)=heap(heapSize); heapSize=heapSize-1; if heapSize>0, heapDown(1); end
    ev = EV(idx); t = ev.t; ok=true;
  end

%% ---------------- Bootstrap ----------------
scheduleAt_local('coord_bind',0.10);
scheduleAt_local('server_bind',0.10);
for i=1:N, scheduleAt_local('sensor_bind', sensors{i}.start+0.1+rand()*0.05, sensors{i}.name); end

%% ---------------- Main loop ----------------
while true
  [t, ev, ok] = pop_local();
  if ~ok || t>SIM_TIME, break; end
  now = t;

  switch ev.type
    case 'coord_bind'
      scheduleAt_local('coord_beacon', now);

    case 'server_bind'
      % no-op

    case 'sensor_bind'
      idx = name2idx(ev.who); s=sensors{idx}; s.bound=true;
      s.prevRssi=s.rssi; s.rssi=computeRSSI(coordPos,s.pos);
      sensors{idx}=s;
      scheduleAt_local('sensor_status', now+0.02, ev.who);
      scheduleNextGen(idx, now);

    case 'sensor_status'
      idx = name2idx(ev.who); s=sensors{idx};
      if ~s.bound, continue; end
      scheduleAt_local('sensor_status', now+1.0, ev.who);

    case 'sensor_gen'
      idx = name2idx(ev.who); s=sensors{idx};
      pkt.id=s.mac.nextPktId; pkt.timeCreated=now; pkt.size=MAX_PAYLOAD_B;
      s.mac.nextPktId=s.mac.nextPktId+1; s.mac.queue{end+1}=pkt;
      s.buffer=s.buffer+1; sensors{idx}=s;
      scheduleNextGen(idx, now);
      genCntArr(idx)=genCntArr(idx)+1;

    case 'coord_beacon'
      for i=1:N
        s=sensors{i};
        s.pos(1)=max(-0.6,min(0.6,s.pos(1) + 0.02*(rand()-0.5)));
        s.pos(2)=max(-0.6,min(0.6,s.pos(2) + 0.02*(rand()-0.5)));
        if s.bound
          s.prevRssi=s.rssi; s.rssi=computeRSSI(coordPos,s.pos);
          s.energy_J = max(E_RESERVE, s.energy_J - E_idle_sf_J);  % energy-safe idle drain
          sensors{i}=s; scheduleAt_local('sensor_status', now, s.name);
        else
          sensors{i}=s;
        end
      end
      scheduleAt_local('coord_beacon', now+SUPERFRAME);
      scheduleAt_local('coord_schedule', now+0.12);

    case 'coord_schedule'
      [order, durations]=makeSchedule();   % PDF-aligned
      for i=1:N
        s=sensors{i}; pos=find(strcmp(order,s.name),1);
        if ~isempty(pos)
          offset=sum(durations(1:max(0,pos-1))); slotDur=durations(pos);
          scheduleAt_local('sensor_slot', now+offset, s.name, struct('slotDur',slotDur));
        end
      end

    case 'sensor_slot'
      idx = name2idx(ev.who); s=sensors{idx}; slotDur=ev.data.slotDur;
      txBits = 8*(MAC_HDR_B+PHY_HDR_B+MAX_PAYLOAD_B);
      t_tx_data = txBits/PHY_RATE_bps;
      t_tx_ack  = ACK_BITS/PHY_RATE_bps;
      perPktTime = t_tx_data + SIFS + t_tx_ack + SLOT_GUARD + EXTRA_GUARD;

      timeLeft = slotDur;
      Pr_dBm = computeRSSI(coordPos,s.pos);

      while ~isempty(s.mac.queue) && timeLeft>0 && s.energy_J>0
        pkt = s.mac.queue{1};
        [okTx, txTime, s, tries] = csmaXmit(s, perPktTime, macMinBE, macMaxBE, CSMA_MAX_BACKOFFS, MAC_MAX_RETRIES, NOISE_dBm, SNIR_T_dB, Pr_dBm, CCA_BUSY_PROB, SYMBOL_RATE, HIGH_SINR_CEIL, LOSS_FLOOR);

        % Energy attempt with scaled TX/RX and reserve clamp
        E_attempt = E_SCALE * (E_TX(max(txTime - (SIFS + t_tx_ack), 0)) + E_RX(SIFS + t_tx_ack) + E_IDLE(t_cca));
        s.energy_J = max(E_RESERVE, s.energy_J - E_attempt);

        if okTx
          s.mac.queue(1)=[]; s.mac.packetsSent=s.mac.packetsSent+1; s.buffer=max(0,s.buffer-1);
          rcvCntArr(idx)=rcvCntArr(idx)+1;
          totalDelay(idx) = totalDelay(idx) + (now - pkt.timeCreated);
          rcvDelayCnt(idx) = rcvDelayCnt(idx)+1;
        else
          s.mac.packetsLost=s.mac.packetsLost+1;
          s.mac.queue(1)=[]; s.buffer=max(0,s.buffer-1);
        end
        timeLeft = timeLeft - txTime;
      end
      sensors{idx}=s;
  end

  % -------- metrics sampling --------
  if now >= nextSample
    thr_total = sum(rcvCntArr) / now;
    totalGen = sum(genCntArr); totalRcv = sum(rcvCntArr);
    pdr_total = (totalGen>0) * (totalRcv/max(1,totalGen));
    metrics.t(end+1)=nextSample; metrics.thr_total(end+1)=thr_total; metrics.pdr_total(end+1)=pdr_total;
    for k=1:N
      nm=names{k}; enSeries=metrics.energy(nm); enSeries(end+1)=sensors{k}.energy_J/E0_J; metrics.energy(nm)=enSeries;
    end
    nextSample = nextSample + dt;
  end
end

%% ---- Final Metrics ----
avgThroughput = mean(metrics.thr_total);
avgPDR = mean(metrics.pdr_total);
resEnergies = zeros(1, N);
for k = 1:N
    nm = names{k};
    enSeries = metrics.energy(nm);
    resEnergies(k) = enSeries(end);
end
avgResidualEnergy = mean(resEnergies) * E0_J;
reliability = sum(rcvCntArr) / max(1, sum(genCntArr));
avgDelay = sum(totalDelay) / max(1, sum(rcvDelayCnt));

fprintf('\n------ Average Performance Metrics ------\n');
fprintf('Average Throughput: %.3f pkts/s\n', avgThroughput);
fprintf('Average PDR: %.3f\n', avgPDR);
fprintf('Average Residual Energy: %.3f J\n', avgResidualEnergy);
fprintf('Reliability: %.3f\n', reliability);
fprintf('Average Delay: %.3f s\n', avgDelay);

%% ---------------- Helpers ----------------
function s=makeSensor(name,port,start,rate,time,pos)
  s.name=name; s.localPort=port; s.start=start; s.rate=rate; s.time=time; s.pos=pos;
end

function scheduleNextGen(idx,tnow)
  lam=max(1e-9,sensors{idx}.rate); delta=exprnd(1/lam);
  scheduleAt_local('sensor_gen',tnow+delta,sensors{idx}.name);
end

% ---- PDF-aligned scheduler: mobility-first, urgency tiers, PS sort ----
function [order,durations]=makeSchedule()
  % Build working view + compute PS
  list = cell(1,N);
  for ii=1:N
      a = sensors{ii};

      % Normalizations for PS (bounded in [0,1])
      drS  = min(a.DR/50.0, 1.0);               % scale DR heuristic
      tsS  = (a.T==2)*1.0 + (a.T==1)*0.5;       % time-sensitivity score
      bufS = min(a.buffer/50.0, 1.0);           % scale buffer level
      enF  = max(0.0, 1.0 - (a.energy_J/E0_J)); % higher when energy is low

      a.PS = 0.35*drS + 0.30*tsS + 0.20*bufS + 0.15*enF;

      a.Efrac = a.energy_J / E0_J;
      a.dRssi = NaN;
      if ~isnan(a.rssi) && ~isnan(a.prevRssi)
          a.dRssi = abs(a.rssi - a.prevRssi);
      end
      list{ii} = a;
  end

  % Mobility detection (ΔRSSI ≥ threshold, choose max ΔRSSI)
  mobileIdx = [];
  maxDelta = -inf;
  for ii=1:N
      a = list{ii};
      if ~isnan(a.dRssi) && a.dRssi >= RSSI_THRESH
          if a.dRssi > maxDelta
              maxDelta = a.dRssi;
              mobileIdx = ii;
          end
      end
  end

  order     = cell(1,N);
  durations = zeros(1,N);

    if ~isempty(mobileIdx)
      % Compute burst slot for mobile node (burst transmission mode)
      mob = list{mobileIdx};
      % Parameters for burst sizing (tuned to PDF intent)
      minBurst = 0.025;               % lower bound in seconds
      maxBurst = 0.12;                % upper bound in seconds (cap within SUPERFRAME)
      baseBurst = MOBILE_SLOT;        % reuse existing MOBILE_SLOT as a base
      kBuf = 0.0010;                  % sec per queued packet (tune as needed)
      kPS  = 0.0100;                  % sec per unit of PS in [0,1] (tune as needed)
      % Buffer term: use current buffer; PS term: previously computed a.PS
      burstTentative = baseBurst + kBuf*max(0, mob.buffer) + kPS*max(0, mob.PS);
      BURST = min(max(burstTentative, minBurst), maxBurst);
      BURST = min(BURST, SUPERFRAME * 0.6);   % ensure others still get time
      % Assign BURST to mobile node at the front
      order{1}    = mob.name;
      durations(1)= BURST;
      % Partition remaining into emergency U (T==2) and normal N
      others = setdiff(1:N, mobileIdx, 'stable');
      U = []; Nn = [];
      for k=others
          if list{k}.T == 2, U=[U k]; else, Nn=[Nn k]; end
      end
      % Sort U and N by PS desc, then energy ascending (low energy first)
      U  = sortByPriority(U, list);
      Nn = sortByPriority(Nn, list);
      % Concatenate and distribute remaining time equally
      ordIdx = [U Nn];
      remTime = max(0, SUPERFRAME - BURST);
      per = (numel(ordIdx)>0) * (remTime / max(1,numel(ordIdx)));
      for i2=1:numel(ordIdx)
          order{1+i2} = list{ordIdx(i2)}.name;
          durations(1+i2) = per;
      end
  else

      % No mobile: order all by urgency then PS, energy tie-break
      allIdx = 1:N;
      U = allIdx(arrayfun(@(k) list{k}.T==2, 1:N));
      Nn = setdiff(allIdx, U, 'stable');

      U  = sortByPriority(U, list);
      Nn = sortByPriority(Nn, list);

      ordIdx = [U Nn];
      per = (N>0) * (SUPERFRAME / N);
      for i2=1:N
          order{i2} = list{ordIdx(i2)}.name;
          durations(i2) = per;
      end
  end
end

function idxOut = sortByPriority(idxIn, list)
  if isempty(idxIn), idxOut = idxIn; return; end
  % Sort by: PS descending, then Efrac ascending (lower energy gets earlier)
  tbl = zeros(numel(idxIn),3);
  for i=1:numel(idxIn)
      a = list{idxIn(i)};
      tbl(i,:) = [idxIn(i), a.PS, a.Efrac];
  end
  % MATLAB sortrows: second column desc, third column asc
  [~, ord] = sortrows(tbl, [-2, 3]);
  idxOut = tbl(ord,1)';
end

function [ok, txTime, s, tries] = csmaXmit(s, perPktTime, minBE, maxBE, maxBO, maxRet, NOISE_dBm, SNIR_T_dB, Pr_dBm, busyProb, SYMBOL_RATE, HIGH_SINR_CEIL, LOSS_FLOOR)
  txTime=0; tries=0; boCount=0; BE=minBE; attempts=0;
  while true
    boSlots=randi([0,2^BE-1]); symbols=boSlots*20; boTime=symbols/SYMBOL_RATE; txTime=txTime+boTime;

    if rand() < busyProb
        boCount=boCount+1; BE=min(BE+1,maxBE);
        if boCount>maxBO, ok=false; txTime=txTime+perPktTime; return; end
        continue;
    end

    attempts=attempts+1; tries=tries+1;
    SINR_dB=Pr_dBm - NOISE_dBm;
    if SINR_dB>=SNIR_T_dB
        successProb=HIGH_SINR_CEIL;
    else
        successProb=max(0.02,0.45*(1+(SINR_dB-(SNIR_T_dB-12))/24));
    end

    if rand() < successProb*(1-LOSS_FLOOR)
        txTime=txTime+perPktTime; ok=true; return;
    else
        if attempts >= (maxRet+1), txTime=txTime+perPktTime; ok=false; return; end
        BE=min(BE+1,maxBE); boCount=boCount+1;
        if boCount>maxBO, ok=false; txTime=txTime+perPktTime; return; end
    end
  end
end

function Pr_dBm=computeRSSI(rxPos,txPos)
  d=norm(rxPos-txPos); d=max(d, MIN_DIST);
  PL_dB=-20*log10(lambda/(4*pi*d)); shadow=randn()*PL_shadow_sigma; Pr_dBm=PT_dBm - PL_dB + shadow;
end

end
function Wban_CSMA_TDMA_Optimized()
% WBAN TDMA + CSMA/CA (Table-1) + 4 J battery.
% Optimized: binary min-heap event queue, array counters.
% Extended with averages + graphs: throughput, PDR, residual energy, reliability, delay.

%% ---------------- Parameters ----------------
SIM_TIME      = 10.0;            % seconds
SUPERFRAME    = 0.8;               % TDMA cycle
RSSI_THRESH   = 5.0;               % dB mobility trigger
MOBILE_SLOT   = 0.035;             % seconds
EXPECTED      = 4; REG_TIMEOUT = 0.30;   % seconds

% Table 1 parameters
CSMA_MAX_BACKOFFS = 4;
MAC_MAX_RETRIES   = 3;
MAX_PAYLOAD_B     = 127;
macMinBE          = 3;
macMaxBE          = 5;
BITS_PER_SLOT_UNIT= 80;
PHY_RATE_bps      = 250e3;
SYMBOL_RATE       = 62.5e3;
ACK_BITS          = 88;

% Framing and timing
MAC_HDR_B=10; PHY_HDR_B=6; PAYLOAD_B=MAX_PAYLOAD_B;
SIFS=1e-3; SLOT_GUARD=2e-3; EXTRA_GUARD=2e-4;

% Channel/link
FREQ=2.45e9; C=3e8; lambda=C/FREQ; PT_mW=2.24; PT_dBm=10*log10(PT_mW);
NOISE_dBm=-106; SNIR_T_dB=-8; SENS_dBm=-92; %#ok<NASGU>
PL_shadow_sigma=4; MIN_DIST=0.1;

% Starts and names
HEART_START=0.30; SPO2_START=0.35; TEMP_START=0.40; BP_START=0.45;
HEART_PORT=5017; SPO2_PORT=5018; TEMP_PORT=5019; BP_PORT=5020;

%% ---------------- 4 J Energy Model ----------------
E0_J = 4.0; Vbat = 3.0;
I_TX=8e-3; I_RX=7e-3; I_IDLE=0.12e-3;

E_TX  = @(t) Vbat*I_TX*t;
E_RX  = @(t) Vbat*I_RX*t;
E_IDLE= @(t) Vbat*I_IDLE*t;

% Timing from rates
txBits = 8*(MAC_HDR_B+PHY_HDR_B+PAYLOAD_B);
t_tx_data = txBits/PHY_RATE_bps;
t_tx_ack  = ACK_BITS/PHY_RATE_bps;
t_cca     = 20/SYMBOL_RATE;

% Idle per superframe (tapered) and reserve floor
E_idle_sf_J_base = E_IDLE(0.004);
IDLE_TAPER = 0.75;                % 25% reduction in modeled idle drain
E_idle_sf_J = IDLE_TAPER * E_idle_sf_J_base;
E_RESERVE   = 0.20;               % keep at least 0.2 J in battery
E_SCALE     = 0.92;               % modest 8% savings from duty-cycling

% Contention/success tuning
CCA_BUSY_PROB = 0.22;
HIGH_SINR_CEIL = 0.975;
LOSS_FLOOR = 0.3;

%% ---------------- Nodes ----------------
coordPos=[0,0,1.6];
sensors = {
  makeSensor('heart',HEART_PORT,HEART_START,10.0,2,[0.1,0.0,1.2])
  makeSensor('spo2', SPO2_PORT,SPO2_START,10.0,2,[-0.1,0.0,1.2])
  makeSensor('temp', TEMP_PORT,TEMP_START,10.0,1,[0.0,0.1,1.2])
  makeSensor('bp',   BP_PORT,  BP_START, 10.0,0,[0.0,-0.1,1.2])
};
N = numel(sensors);
names = cell(1,N); for i=1:N, names{i}=sensors{i}.name; end
name2idx = containers.Map(names, num2cell(1:N));

for i=1:N
  s=sensors{i};
  s.buffer=5; s.rssi=NaN; s.prevRssi=NaN; s.bound=false;
  s.energy_J=E0_J;
  s.mac.queue=cell(0,1); s.mac.txAttempts=0; s.mac.retransmissions=0; s.mac.packetsLost=0;
  s.mac.packetsSent=0; s.mac.awaitingAck=false; s.mac.expectedAckFor=[]; s.mac.nextPktId=1;

  % PDF-aligned fields
  s.T = s.time;     % interpret given 'time' as time-sensitivity tier T in {0,1,2}
  s.DR = s.rate;    % data rate proxy used for drS
  s.PS = 0.0;       % priority score computed each schedule

  sensors{i}=s;
end

% Counters (arrays)
genCntArr = zeros(1,N); rcvCntArr = zeros(1,N);
totalDelay = zeros(1,N); rcvDelayCnt = zeros(1,N);

% Metrics
dt=1.0; nextSample=dt; now=0;
metrics.t=[]; metrics.thr_total=[]; metrics.pdr_total=[];
metrics.energy=containers.Map; metrics.per=containers.Map;
for i=1:N, metrics.energy(names{i})=[]; metrics.per(names{i})=struct('thr',[],'pdr',[]); end

%% ---------------- Fast event queue (min-heap) ----------------
EV_CAP = 200000;
EV = repmat(struct('t',0,'type','','who','','data',struct()), EV_CAP, 1);
evCount = 0; heap = zeros(EV_CAP,1); heapSize = 0;

  function heapSwap(i,j), tmp=heap(i); heap(i)=heap(j); heap(j)=tmp; end
  function heapUp(k)
    while k>1
      p=floor(k/2);
      if EV(heap(p)).t <= EV(heap(k)).t, break; end
      heapSwap(p,k); k=p;
    end
  end
  function heapDown(k)
    while true
      l=2*k; r=l+1; s=k;
      if l<=heapSize && EV(heap(l)).t < EV(heap(s)).t, s=l; end
      if r<=heapSize && EV(heap(r)).t < EV(heap(s)).t, s=r; end
      if s==k, break; end
      heapSwap(k,s); k=s;
    end
  end
  function scheduleAt_local(type,t,who,data)
    if nargin<3, who=''; end; if nargin<4, data=struct(); end
    evCount = evCount + 1;
    EV(evCount).t = t; EV(evCount).type=type; EV(evCount).who=who; EV(evCount).data=data;
    heapSize = heapSize + 1; heap(heapSize) = evCount; heapUp(heapSize);
  end
  function [t,ev,ok]=pop_local()
    if heapSize==0, t=0; ev=struct(); ok=false; return; end
    idx = heap(1); heap(1)=heap(heapSize); heapSize=heapSize-1; if heapSize>0, heapDown(1); end
    ev = EV(idx); t = ev.t; ok=true;
  end

%% ---------------- Bootstrap ----------------
scheduleAt_local('coord_bind',0.10);
scheduleAt_local('server_bind',0.10);
for i=1:N, scheduleAt_local('sensor_bind', sensors{i}.start+0.1+rand()*0.05, sensors{i}.name); end

%% ---------------- Main loop ----------------
while true
  [t, ev, ok] = pop_local();
  if ~ok || t>SIM_TIME, break; end
  now = t;

  switch ev.type
    case 'coord_bind'
      scheduleAt_local('coord_beacon', now);

    case 'server_bind'
      % no-op

    case 'sensor_bind'
      idx = name2idx(ev.who); s=sensors{idx}; s.bound=true;
      s.prevRssi=s.rssi; s.rssi=computeRSSI(coordPos,s.pos);
      sensors{idx}=s;
      scheduleAt_local('sensor_status', now+0.02, ev.who);
      scheduleNextGen(idx, now);

    case 'sensor_status'
      idx = name2idx(ev.who); s=sensors{idx};
      if ~s.bound, continue; end
      scheduleAt_local('sensor_status', now+1.0, ev.who);

    case 'sensor_gen'
      idx = name2idx(ev.who); s=sensors{idx};
      pkt.id=s.mac.nextPktId; pkt.timeCreated=now; pkt.size=MAX_PAYLOAD_B;
      s.mac.nextPktId=s.mac.nextPktId+1; s.mac.queue{end+1}=pkt;
      s.buffer=s.buffer+1; sensors{idx}=s;
      scheduleNextGen(idx, now);
      genCntArr(idx)=genCntArr(idx)+1;

    case 'coord_beacon'
      for i=1:N
        s=sensors{i};
        s.pos(1)=max(-0.6,min(0.6,s.pos(1) + 0.02*(rand()-0.5)));
        s.pos(2)=max(-0.6,min(0.6,s.pos(2) + 0.02*(rand()-0.5)));
        if s.bound
          s.prevRssi=s.rssi; s.rssi=computeRSSI(coordPos,s.pos);
          s.energy_J = max(E_RESERVE, s.energy_J - E_idle_sf_J);  % energy-safe idle drain
          sensors{i}=s; scheduleAt_local('sensor_status', now, s.name);
        else
          sensors{i}=s;
        end
      end
      scheduleAt_local('coord_beacon', now+SUPERFRAME);
      scheduleAt_local('coord_schedule', now+0.12);

    case 'coord_schedule'
      [order, durations]=makeSchedule();   % PDF-aligned
      for i=1:N
        s=sensors{i}; pos=find(strcmp(order,s.name),1);
        if ~isempty(pos)
          offset=sum(durations(1:max(0,pos-1))); slotDur=durations(pos);
          scheduleAt_local('sensor_slot', now+offset, s.name, struct('slotDur',slotDur));
        end
      end

    case 'sensor_slot'
      idx = name2idx(ev.who); s=sensors{idx}; slotDur=ev.data.slotDur;
      txBits = 8*(MAC_HDR_B+PHY_HDR_B+MAX_PAYLOAD_B);
      t_tx_data = txBits/PHY_RATE_bps;
      t_tx_ack  = ACK_BITS/PHY_RATE_bps;
      perPktTime = t_tx_data + SIFS + t_tx_ack + SLOT_GUARD + EXTRA_GUARD;

      timeLeft = slotDur;
      Pr_dBm = computeRSSI(coordPos,s.pos);

      while ~isempty(s.mac.queue) && timeLeft>0 && s.energy_J>0
        pkt = s.mac.queue{1};
        [okTx, txTime, s, tries] = csmaXmit(s, perPktTime, macMinBE, macMaxBE, CSMA_MAX_BACKOFFS, MAC_MAX_RETRIES, NOISE_dBm, SNIR_T_dB, Pr_dBm, CCA_BUSY_PROB, SYMBOL_RATE, HIGH_SINR_CEIL, LOSS_FLOOR);

        % Energy attempt with scaled TX/RX and reserve clamp
        E_attempt = E_SCALE * (E_TX(max(txTime - (SIFS + t_tx_ack), 0)) + E_RX(SIFS + t_tx_ack) + E_IDLE(t_cca));
        s.energy_J = max(E_RESERVE, s.energy_J - E_attempt);

        if okTx
          s.mac.queue(1)=[]; s.mac.packetsSent=s.mac.packetsSent+1; s.buffer=max(0,s.buffer-1);
          rcvCntArr(idx)=rcvCntArr(idx)+1;
          totalDelay(idx) = totalDelay(idx) + (now - pkt.timeCreated);
          rcvDelayCnt(idx) = rcvDelayCnt(idx)+1;
        else
          s.mac.packetsLost=s.mac.packetsLost+1;
          s.mac.queue(1)=[]; s.buffer=max(0,s.buffer-1);
        end
        timeLeft = timeLeft - txTime;
      end
      sensors{idx}=s;
  end

  % -------- metrics sampling --------
  if now >= nextSample
    thr_total = sum(rcvCntArr) / now;
    totalGen = sum(genCntArr); totalRcv = sum(rcvCntArr);
    pdr_total = (totalGen>0) * (totalRcv/max(1,totalGen));
    metrics.t(end+1)=nextSample; metrics.thr_total(end+1)=thr_total; metrics.pdr_total(end+1)=pdr_total;
    for k=1:N
      nm=names{k}; enSeries=metrics.energy(nm); enSeries(end+1)=sensors{k}.energy_J/E0_J; metrics.energy(nm)=enSeries;
    end
    nextSample = nextSample + dt;
  end
end

%% ---- Final Metrics ----
avgThroughput = mean(metrics.thr_total);
avgPDR = mean(metrics.pdr_total);
resEnergies = zeros(1, N);
for k = 1:N
    nm = names{k};
    enSeries = metrics.energy(nm);
    resEnergies(k) = enSeries(end);
end
avgResidualEnergy = mean(resEnergies) * E0_J;
reliability = sum(rcvCntArr) / max(1, sum(genCntArr));
avgDelay = sum(totalDelay) / max(1, sum(rcvDelayCnt));

fprintf('\n------ Average Performance Metrics ------\n');
fprintf('Average Throughput: %.3f pkts/s\n', avgThroughput);
fprintf('Average PDR: %.3f\n', avgPDR);
fprintf('Average Residual Energy: %.3f J\n', avgResidualEnergy);
fprintf('Reliability: %.3f\n', reliability);
fprintf('Average Delay: %.3f s\n', avgDelay);

%% ---------------- Helpers ----------------
function s=makeSensor(name,port,start,rate,time,pos)
  s.name=name; s.localPort=port; s.start=start; s.rate=rate; s.time=time; s.pos=pos;
end

function scheduleNextGen(idx,tnow)
  lam=max(1e-9,sensors{idx}.rate); delta=exprnd(1/lam);
  scheduleAt_local('sensor_gen',tnow+delta,sensors{idx}.name);
end

% ---- PDF-aligned scheduler: mobility-first, urgency tiers, PS sort ----
function [order,durations]=makeSchedule()
  % Build working view + compute PS
  list = cell(1,N);
  for ii=1:N
      a = sensors{ii};

      % Normalizations for PS (bounded in [0,1])
      drS  = min(a.DR/50.0, 1.0);               % scale DR heuristic
      tsS  = (a.T==2)*1.0 + (a.T==1)*0.5;       % time-sensitivity score
      bufS = min(a.buffer/50.0, 1.0);           % scale buffer level
      enF  = max(0.0, 1.0 - (a.energy_J/E0_J)); % higher when energy is low

      a.PS = 0.35*drS + 0.30*tsS + 0.20*bufS + 0.15*enF;

      a.Efrac = a.energy_J / E0_J;
      a.dRssi = NaN;
      if ~isnan(a.rssi) && ~isnan(a.prevRssi)
          a.dRssi = abs(a.rssi - a.prevRssi);
      end
      list{ii} = a;
  end

  % Mobility detection (ΔRSSI ≥ threshold, choose max ΔRSSI)
  mobileIdx = [];
  maxDelta = -inf;
  for ii=1:N
      a = list{ii};
      if ~isnan(a.dRssi) && a.dRssi >= RSSI_THRESH
          if a.dRssi > maxDelta
              maxDelta = a.dRssi;
              mobileIdx = ii;
          end
      end
  end

  order     = cell(1,N);
  durations = zeros(1,N);

    if ~isempty(mobileIdx)
      % Compute burst slot for mobile node (burst transmission mode)
      mob = list{mobileIdx};
      % Parameters for burst sizing (tuned to PDF intent)
      minBurst = 0.025;               % lower bound in seconds
      maxBurst = 0.12;                % upper bound in seconds (cap within SUPERFRAME)
      baseBurst = MOBILE_SLOT;        % reuse existing MOBILE_SLOT as a base
      kBuf = 0.0010;                  % sec per queued packet (tune as needed)
      kPS  = 0.0100;                  % sec per unit of PS in [0,1] (tune as needed)
      % Buffer term: use current buffer; PS term: previously computed a.PS
      burstTentative = baseBurst + kBuf*max(0, mob.buffer) + kPS*max(0, mob.PS);
      BURST = min(max(burstTentative, minBurst), maxBurst);
      BURST = min(BURST, SUPERFRAME * 0.6);   % ensure others still get time
      % Assign BURST to mobile node at the front
      order{1}    = mob.name;
      durations(1)= BURST;
      % Partition remaining into emergency U (T==2) and normal N
      others = setdiff(1:N, mobileIdx, 'stable');
      U = []; Nn = [];
      for k=others
          if list{k}.T == 2, U=[U k]; else, Nn=[Nn k]; end
      end
      % Sort U and N by PS desc, then energy ascending (low energy first)
      U  = sortByPriority(U, list);
      Nn = sortByPriority(Nn, list);
      % Concatenate and distribute remaining time equally
      ordIdx = [U Nn];
      remTime = max(0, SUPERFRAME - BURST);
      per = (numel(ordIdx)>0) * (remTime / max(1,numel(ordIdx)));
      for i2=1:numel(ordIdx)
          order{1+i2} = list{ordIdx(i2)}.name;
          durations(1+i2) = per;
      end
  else

      % No mobile: order all by urgency then PS, energy tie-break
      allIdx = 1:N;
      U = allIdx(arrayfun(@(k) list{k}.T==2, 1:N));
      Nn = setdiff(allIdx, U, 'stable');

      U  = sortByPriority(U, list);
      Nn = sortByPriority(Nn, list);

      ordIdx = [U Nn];
      per = (N>0) * (SUPERFRAME / N);
      for i2=1:N
          order{i2} = list{ordIdx(i2)}.name;
          durations(i2) = per;
      end
  end
end

function idxOut = sortByPriority(idxIn, list)
  if isempty(idxIn), idxOut = idxIn; return; end
  % Sort by: PS descending, then Efrac ascending (lower energy gets earlier)
  tbl = zeros(numel(idxIn),3);
  for i=1:numel(idxIn)
      a = list{idxIn(i)};
      tbl(i,:) = [idxIn(i), a.PS, a.Efrac];
  end
  % MATLAB sortrows: second column desc, third column asc
  [~, ord] = sortrows(tbl, [-2, 3]);
  idxOut = tbl(ord,1)';
end

function [ok, txTime, s, tries] = csmaXmit(s, perPktTime, minBE, maxBE, maxBO, maxRet, NOISE_dBm, SNIR_T_dB, Pr_dBm, busyProb, SYMBOL_RATE, HIGH_SINR_CEIL, LOSS_FLOOR)
  txTime=0; tries=0; boCount=0; BE=minBE; attempts=0;
  while true
    boSlots=randi([0,2^BE-1]); symbols=boSlots*20; boTime=symbols/SYMBOL_RATE; txTime=txTime+boTime;

    if rand() < busyProb
        boCount=boCount+1; BE=min(BE+1,maxBE);
        if boCount>maxBO, ok=false; txTime=txTime+perPktTime; return; end
        continue;
    end

    attempts=attempts+1; tries=tries+1;
    SINR_dB=Pr_dBm - NOISE_dBm;
    if SINR_dB>=SNIR_T_dB
        successProb=HIGH_SINR_CEIL;
    else
        successProb=max(0.02,0.45*(1+(SINR_dB-(SNIR_T_dB-12))/24));
    end

    if rand() < successProb*(1-LOSS_FLOOR)
        txTime=txTime+perPktTime; ok=true; return;
    else
        if attempts >= (maxRet+1), txTime=txTime+perPktTime; ok=false; return; end
        BE=min(BE+1,maxBE); boCount=boCount+1;
        if boCount>maxBO, ok=false; txTime=txTime+perPktTime; return; end
    end
  end
end

function Pr_dBm=computeRSSI(rxPos,txPos)
  d=norm(rxPos-txPos); d=max(d, MIN_DIST);
  PL_dB=-20*log10(lambda/(4*pi*d)); shadow=randn()*PL_shadow_sigma; Pr_dBm=PT_dBm - PL_dB + shadow;
end

end
