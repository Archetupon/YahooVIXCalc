function [MFIV] = YahooMFIVCalc(ticker)

%FUNCTION
%   Function calculate the model-free implied variance (VIX style) for any
%   stock, ETF, or index with option chains available on Yahoo Finance
%   ****Note that you need to have cell2str.m in your directory as well. 
%   This is available for download on the FileExchange****
%
%INPUT arguments: 
%   Ticker: string format! i.e. 'SPY'. Yahoo Finance compatible ticker, for index there should be
%   '^...', i.e. '^SPX' for the S&P500 index.
%
%OUTPUT arguments:
%   MFIV - the model-free implied volatility (square root of model-free implied variance), 
%   calculated via the VIX methodology.
%
%EXAMPLE
%   The following would calculate the MFIV for the SPY ETF
%   [MFIV_SPY] = YahooMFIVCalc('SPY');

%++ Contrust url ++%   
Basename1 = 'https://finance.yahoo.com/quote/';
Basename2 = '/options?p=';
url1 = [Basename1,ticker,Basename2,ticker];

%++ Grab Yahoo Finance base url and html code ++%
html1 = webread(url1);
%***************************************************%

%++ Grab html code which contains active maturities ++%
pattern1 = 'OptionContractsStore';% start of list
pattern2 = 'hasMiniOptions';     % end of list
list = extractBetween(html1, pattern1, pattern2); 
list = cell2str(list);
%***************************************************%

%++ Grab active maturities and choose closest +/- 30 days per VIX method ++%
pattern1 = '"expirationDates":[';% start of list
pattern2 = '],"';     % end of list
maturities = extractBetween(list, pattern1, pattern2);
maturities = cell2str(maturities);
maturities = maturities(3:end-3);
maturities = strsplit(maturities,',');
maturities_store = strings(length(maturities),1);
for i=1:length(maturities_store)
    maturities_store(i) = maturities{i};
end
maturities_store = str2double(maturities_store);
%++ Convert from UNIX/Epoch time ++%
epoch = '01-Jan-1970';
formatIn = 'dd-mmm-yyyy';
time_matlab = ones(length(maturities_store),1).*datenum(epoch,formatIn) + maturities_store./86400;
t_30 = datenum(datetime(date) + days(30)); % 30 calendar days from now
%++ Find maturities to use, ensure they are Fridays (other days can have
% few strikes and liquidity problems ++%
mat_diff_pos = daysact(t_30,time_matlab);
mat_diff_pos(mat_diff_pos<=0) = nan;
mat_diff_neg = daysact(t_30,time_matlab);
mat_diff_neg(mat_diff_neg>0) = nan;
mat_diff_neg = abs(mat_diff_neg);
mat_index = zeros(2,1);
[~,mat_index(1)] = min(mat_diff_neg);
[~,mat_index(2)] = min(mat_diff_pos);
for i=1:100 %100 just a shill number necessarily large
   if  weekday(time_matlab(mat_index(1))) ~= 6
       mat_index(1) = mat_index(1)-1;
   else
       break
   end
end
for i=1:100 %100 just a shill number necessarily large
   if  weekday(time_matlab(mat_index(2))) ~= 6
       mat_index(2) = mat_index(2)+1;
   else
       break
   end
end
maturities_hooman_read = datestr(time_matlab(mat_index),'yyyy-mm-dd');
user_maturities = maturities_hooman_read; % Remnant of YahooOptionScraper
%***************************************************%

%++ Pull option chains from Yahoo Finance ++% 
call_cell = cell(size(user_maturities,1),1);
put_cell = cell(size(user_maturities,1),1);
join_cell = cell(size(user_maturities,1),1);
for k=1:size(user_maturities,1)

 %++ Change expiration date back to UNIX time per Yahoo url format ++%
 d = datetime(user_maturities(k,:));
 d = num2str(posixtime(d));
 %***************************************************%

 %++ Contruct url ++%
 Basename3 = 'https://finance.yahoo.com/quote/';
 Basename4 = '/options?p=';
 Basename5 = '&date=';
 url2 = [Basename3,ticker,Basename4,ticker,Basename5,d];
 
 %++ Grab Yahoo Finance base url and html code ++%
 html2 = webread(url2);
 %***************************************************%

 %++ Start with call options first ++%
 pattern1 = '{"calls":[{';% start of list
 pattern2 = '"puts":[{';     % end of list
 list = extractBetween(html2, pattern1, pattern2); 
 list = cell2str(list);
 %***************************************************%

 %++ First extract the formatted values ++%
 %++ This applies to all values expect the percentages ++%
 pattern1 = '"fmt":"';
 pattern2 = '"';
 floats = extractBetween(list, pattern1, pattern2); 
 %***************************************************%

 %++ Now extract the raw values ++%
 %++ This is needed for the percetages ++%
 pattern1 = '"raw":';
 pattern2 = ',';
 decimals = extractBetween(list, pattern1, pattern2); 
 %***************************************************%

 %++ Initialize string vectors to store values ++%
 implied_volatilities = strings(length(decimals)/11,1);
 expirations = strings(length(floats)/11,1);
 changes = strings(length(floats)/11,1);
 strikes = strings(length(floats)/11,1);
 last_prices = strings(length(floats)/11,1);
 open_interests = strings(length(floats)/11,1);
 percentage_changes = strings(length(decimals)/11,1);
 ask_prices = strings(length(floats)/11,1);
 volumes = strings(length(floats)/11,1);
 last_trade_dates = strings(length(floats)/11,1);
 bid_prices = strings(length(floats)/11,1);
 %***************************************************%

 %++ Loop through, extracting option data ++%
 count = 1;
 for i=1:11:length(floats)
     implied_volatilities(count) = decimals{i};
     expirations(count) = floats{i+1};
     changes(count) = floats{i+2};
     strikes(count) = floats{i+3};
     last_prices(count) = floats{i+4};
     open_interests(count) = floats{i+5};
     percentage_changes(count) = decimals{i+6};
     ask_prices(count) = floats{i+7};
     volumes(count) = floats{i+8};
     last_trade_dates(count) = floats{i+9};
     bid_prices(count) = floats{i+10};
     count = count + 1;
 end
 %***************************************************%

 %+++ Convert data points to appropriate formats +++%
 implied_volatilities = str2double(implied_volatilities);
 expirations = datetime(expirations,'InputFormat','yyyy-MM-dd');
 changes = str2double(changes);
 strikes = str2double(strikes);
 last_prices = str2double(last_prices);
 open_interests = str2double(open_interests);
 percentage_changes = str2double(percentage_changes);
 ask_prices = str2double(ask_prices);
 volumes = str2double(volumes);
 last_trade_dates = datetime(last_trade_dates,'InputFormat','yyyy-MM-dd');
 bid_prices = str2double(bid_prices);
 %***************************************************%

 %++ Create table of option data ++%
 call_table = table(last_trade_dates,expirations,strikes,last_prices,bid_prices,ask_prices,...
     changes,percentage_changes,volumes,open_interests,implied_volatilities);
 call_table.Properties.VariableNames = {'Last_Trade_Date','Expiration','Strike','Last_Price','Bid',...
     'Ask','Change','Per_Change','Volume','Open_Interest','Implied_Volatility'};
 call_cell{k} = call_table;
 %++ Calculate bid/ask mid as market proxy
 b_a_mid = (table2array(call_cell{k}(:,5)) + table2array(call_cell{k}(:,6)))./2;
 call_cell{k} = [call_cell{k}(:,1:6) array2table(b_a_mid) call_cell{k}(:,7:end)];
 %***************************************************%

 %++ Now onto the put options ++%
 pattern1 = '"puts":[{';% start of list
 pattern2 = '}}],"sortColumn"';     % end of list
 list = extractBetween(html2, pattern1, pattern2); 
 list = cell2str(list);
 %***************************************************%

 %++ First extract the formatted values ++%
 %++ This applies to all values expect the percentages ++%
 pattern1 = '"fmt":"';
 pattern2 = '"';
 floats = extractBetween(list, pattern1, pattern2); 
 %***************************************************%

 %++ Now extract the raw values ++%
 %++ This is needed for the percetages ++%
 pattern1 = '"raw":';
 pattern2 = ',';
 decimals = extractBetween(list, pattern1, pattern2); 
 %***************************************************%

 %++ Initialize string vectors to store values ++%
 implied_volatilities = strings(length(decimals)/11,1);
 expirations = strings(length(floats)/11,1);
 changes = strings(length(floats)/11,1);
 strikes = strings(length(floats)/11,1);
 last_prices = strings(length(floats)/11,1);
 open_interests = strings(length(floats)/11,1);
 percentage_changes = strings(length(decimals)/11,1);
 ask_prices = strings(length(floats)/11,1);
 volumes = strings(length(floats)/11,1);
 last_trade_dates = strings(length(floats)/11,1);
 bid_prices = strings(length(floats)/11,1);
 %***************************************************%

 %++ Loop through, extracting option data ++%
 count = 1;
 for i=1:11:length(floats)
     implied_volatilities(count) = decimals{i};
     expirations(count) = floats{i+1};
     changes(count) = floats{i+2};
     strikes(count) = floats{i+3};
     last_prices(count) = floats{i+4};
     open_interests(count) = floats{i+5};
     percentage_changes(count) = decimals{i+6};
     ask_prices(count) = floats{i+7};
     volumes(count) = floats{i+8};
     last_trade_dates(count) = floats{i+9};
     bid_prices(count) = floats{i+10};
     count = count + 1;
 end
 %***************************************************%

 %+++ Convert data points to appropriate formats +++%
 implied_volatilities = str2double(implied_volatilities);
 expirations = datetime(expirations,'InputFormat','yyyy-MM-dd');
 changes = str2double(changes);
 strikes = str2double(strikes);
 last_prices = str2double(last_prices);
 open_interests = str2double(open_interests);
 percentage_changes = str2double(percentage_changes);
 ask_prices = str2double(ask_prices);
 volumes = str2double(volumes);
 last_trade_dates = datetime(last_trade_dates,'InputFormat','yyyy-MM-dd');
 bid_prices = str2double(bid_prices);
 %***************************************************%

 %++ Create table of option data ++%
 put_table = table(last_trade_dates,expirations,strikes,last_prices,bid_prices,ask_prices,...
     changes,percentage_changes,volumes,open_interests,implied_volatilities);
 put_table = put_table(1:(size(put_table,1)-size(call_table,1)),:);
 put_table.Properties.VariableNames = {'Last_Trade_Date','Expiration','Strike','Last_Price','Bid',...
     'Ask','Change','Per_Change','Volume','Open_Interest','Implied_Volatility'};
 put_cell{k} = put_table;
 %++ Calculate bid/ask mid as market proxy
 b_a_mid = (table2array(put_cell{k}(:,5)) + table2array(put_cell{k}(:,6)))./2;
 put_cell{k} = [put_cell{k}(:,1:6) array2table(b_a_mid) put_cell{k}(:,7:end)];
 %***************************************************%
 
 %++ Create join table for the "Forward SPX level" calc ++%
 join_cell{k} = innerjoin(call_cell{k},put_cell{k},'LeftKeys','Strike','RightKeys','Strike');
 
 %***************************************************%

%  %++ Save tables as Excel files ++%
%  call_file_name =[date,ticker,'CallOptionChain',user_maturities(k,:),'.xlsx'];
%  writetable(call_table,call_file_name)
%  put_file_name =[date,ticker,'PutOptionChain',user_maturities(k,:),'.xlsx'];
%  writetable(put_table,put_file_name)
%  %***************************************************%
end

%++ Pull some additional data ++%
% Spot Price%
pattern1 = '"regularMarketPrice":{"raw":';% start of list
pattern2 = ',"fmt":"';     % end of list
list = extractBetween(html2, pattern1, pattern2); 
spot_price = str2double(list{1});

% 1M Interest Rate
ratename = 'http://www.global-rates.com/interest-rates/libor/american-dollar/usd-libor-interest-rate-1-month.aspx';
html3 = webread(ratename);
pattern1 = 'class="tabledata1">';% start of list
pattern2 = 'class="tabledata2">';     % end of list
list = extractBetween(html3, pattern1, pattern2); 
list = cell2str(list);
pattern1 = '"center">';
pattern2 = '&nbsp;';
rates = extractBetween(list, pattern1, pattern2);
interest_rate = str2double(rates{1})/100;
%***************************************************%

%++ Find strikes for minimum call-put for "forward" ++% 
[~,f_strike_near_index] = min(abs(table2array(join_cell{1}(:,7))-table2array(join_cell{1}(:,18))));
f_strike_near = table2array(join_cell{1}(f_strike_near_index,3));
[~,f_strike_next_index] = min(abs(table2array(join_cell{2}(:,7))-table2array(join_cell{2}(:,18))));
f_strike_next = table2array(join_cell{2}(f_strike_next_index,3));
days_to_near = datenum(table2array(join_cell{1}(1,2)))-datenum(date);
days_to_next = datenum(table2array(join_cell{2}(1,2)))-datenum(date);
F_1 = f_strike_near + exp(interest_rate*(days_to_near/365))*(...
    table2array(join_cell{1}(f_strike_near_index,7))-table2array(join_cell{1}(f_strike_near_index,18)));
F_2 = f_strike_next + exp(interest_rate*(days_to_next/365))*(...
    table2array(join_cell{2}(f_strike_next_index,7))-table2array(join_cell{2}(f_strike_next_index,18)));
F_1_test = (table2array(join_cell{1}(:,3)) <= F_1).*table2array(join_cell{1}(:,3));
F_2_test = (table2array(join_cell{2}(:,3)) <= F_2).*table2array(join_cell{2}(:,3));
K_vec = zeros(size(user_maturities,1),1);
K_vec(1) = max(F_1_test); 
K_vec(2) = max(F_2_test);
%***************************************************%

%++ Cut off put and call matrices based on K values ++%
put_call_averages = zeros(size(user_maturities,1),1);
for i=1:size(user_maturities,1)
 %++ Calculate first the put/call averages of the Kth strikes ++%
 put_avg_index = table2array(put_cell{i}(:,3)) == K_vec(i);
 call_avg_index = table2array(call_cell{i}(:,3)) == K_vec(i);
 put_call_averages(i) = (table2array(call_cell{i}(call_avg_index,7))+...
     table2array(put_cell{i}(put_avg_index,7)))/2;
 put_logical = table2array(put_cell{i}(:,3)) < K_vec(i);
 put_indices = find(put_logical);
 put_cell{i} = put_cell{i}(put_indices,:);
 call_logical = table2array(call_cell{i}(:,3)) > K_vec(i);
 call_indices = find(call_logical);
 %++ Now cut off the tables ++%
 put_cell{i} = put_cell{i}(put_indices,:);
 call_cell{i} = call_cell{i}(call_indices,:);
end
%***************************************************%

%++ Remove zero bid price options, and also cut off when 2 zero bids ++%
for i=1:size(user_maturities,1)
    put_zero_bid_logical = find(table2array(put_cell{i}(:,5)));
    put_zero_bid_index = zeros(size(put_zero_bid_logical,1),1);
    put_zero_bid_index(end) = size(put_zero_bid_logical,1);
    for j=(size(put_zero_bid_logical,1)-1):-1:1
        if put_zero_bid_logical(j) ~= 0
            put_zero_bid_index(j) = j;
        elseif put_zero_bid_logical(j) == 0 && put_zero_bid_logical(j+1) == 0
            break
        end
    end
    call_zero_bid_logical = find(table2array(call_cell{i}(:,5)));
    call_zero_bid_index = zeros(size(call_zero_bid_logical,1),1);
    call_zero_bid_index(1) = 1;
    for j=2:size(call_zero_bid_logical,1)
        if call_zero_bid_logical(j) ~= 0
            call_zero_bid_index(j) = j;
        elseif call_zero_bid_logical(j) == 0 && call_zero_bid_logical(j-1) == 0
            break
        end
    end
    call_cell{i} = call_cell{i}(call_zero_bid_index,:);
    put_cell{i} = put_cell{i}(put_zero_bid_index,:);
end
%***************************************************%

%++ New tables of relevant info, including atm half call/put ++%
near_temp = array2table([K_vec(1) put_call_averages(1)]);
near_temp.Properties.VariableNames = {'Strike','b_a_mid'};
next_temp = array2table([K_vec(2) put_call_averages(2)]);
next_temp.Properties.VariableNames = {'Strike','b_a_mid'};
table_indices = [3;7];
near_term_table = [put_cell{1}(:,table_indices);near_temp;call_cell{1}(:,table_indices)];
next_term_table = [put_cell{2}(:,table_indices);next_temp;call_cell{2}(:,table_indices)];
%***************************************************%

%++ Calculate delta Ks ++%
delta_k_1 = zeros(size(near_term_table,1),1);
delta_k_1(1) = table2array(near_term_table(2,1)) - table2array(near_term_table(1,1));
delta_k_1(end) = table2array(near_term_table(end,1)) - table2array(near_term_table(end-1,1));
delta_k_2 = zeros(size(next_term_table,1),1);
delta_k_2(1) = table2array(next_term_table(2,1)) - table2array(next_term_table(1,1));
delta_k_2(end) = table2array(next_term_table(end,1)) - table2array(next_term_table(end-1,1));
for i =2:(size(delta_k_1,1)-1)
    delta_k_1(i) = (table2array(near_term_table(i+1,1)) - table2array(near_term_table(i-1,1)))/2;
end
for i =2:(size(delta_k_2,1)-1)
    delta_k_2(i) = (table2array(next_term_table(i+1,1)) - table2array(next_term_table(i-1,1)))/2;
end
%***************************************************%

%++ Calculate implied variances ++%
var_1 = zeros(size(delta_k_1,1),1);
var_2 = zeros(size(delta_k_2,1),1);
for i =1:size(delta_k_1,1)
    var_1(i) = (delta_k_1(i)/(table2array(near_term_table(i,1))^2))*exp(interest_rate*(days_to_near/365))...
        *table2array(near_term_table(i,2));
end
imp_var_1 = sum(var_1);
imp_var_1 = imp_var_1*2/(days_to_near/365) - 1/(days_to_near/365)*(F_1/K_vec(1)-1)^2;
for i =1:size(delta_k_2,1)
    var_2(i) = (delta_k_2(i)/(table2array(next_term_table(i,1))^2))*exp(interest_rate*(days_to_next/365))...
        *table2array(next_term_table(i,2));
end
imp_var_2 = sum(var_2);
imp_var_2 = imp_var_2*2/(days_to_next/365) - 1/(days_to_next/365)*(F_2/K_vec(2)-1)^2;
%***************************************************%

%++ Calculate the volatility index via interpolation ++%
int_factor_1 = (days_to_next-30)/(days_to_next-days_to_near);
int_factor_2 = (30-days_to_near)/(days_to_near-days_to_next);
MFIV = 100*sqrt(((days_to_near/365)*imp_var_1*int_factor_1+(days_to_next/365)*imp_var_2*int_factor_2)*(365/30));
%***************************************************%

end
