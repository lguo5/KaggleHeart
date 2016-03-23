%% load and compare Call vs True
path_csv_call='/Users/liheng/Dropbox/KaggleBowl/LG/regGroTest_train/vol.csv';
path_csv_true='/Users/liheng/Datasets/DSB16/train.csv';

vol_call_raw=dlmread(path_csv_call,',',1,0); % skip row 1??
vol_call_sys=vol_call_raw(:,2);
vol_call_dia=vol_call_raw(:,3);
nsubjcsv=length(vol_call_sys);

vol_true_raw=dlmread(path_csv_true,',',1,0); % skip row 1
vol_true_sys=vol_true_raw(1:nsubjcsv,2);
vol_true_dia=vol_true_raw(1:nsubjcsv,3);

% plot
figure; plot(vol_call_dia,vol_true_dia,'s'); xlabel('call (dia)'); ylabel('true (dia)');
figure; plot(vol_call_sys,vol_true_sys,'s'); xlabel('call (sys)'); ylabel('true (sys)');

% fit
tbl = table(vol_call_dia,vol_true_dia,vol_call_sys,vol_true_sys);
fit_dia=fitlm(tbl,'vol_true_dia~vol_call_dia')
fit_sys=fitlm(tbl,'vol_true_sys~vol_call_sys')

%% For submission: "Call" CSV only: load csv
path_csv_call='/Users/liheng/Dropbox/KaggleBowl/LG/regGroTest_test/vol.csv';
vol_call_raw=dlmread(path_csv_call,',',0,0); % <<<<<<< skip row 1??
% vol_call_ids=vol_call_raw(:,1);
vol_call_sys=vol_call_raw(:,2);
vol_call_dia=vol_call_raw(:,3);
nsubjcsv=length(vol_call_sys)

% vol_call_sys(nsubjcsv+1:1140)=mean(vol_call_sys);
% vol_call_dia(nsubjcsv+1:1140)=mean(vol_call_dia);
% vol_call_sys=vol_call_raw(:,2)+20;
% vol_call_dia=vol_call_raw(:,3)+20;
% nsubjcsv=length(vol_call_sys)

%% write submission
vol_call_export=[vol_call_dia,vol_call_sys]';
vol_call_export=vol_call_export(:);e

CDFs=softstep(repmat(0:599,[nsubjcsv*2 1]),vol_call_export,5);  % debug plot:   CDFs'    
writeCSV(CDFs,701:1140);
% writeCSV(CDFs,701:701+nsubjcsv-1);