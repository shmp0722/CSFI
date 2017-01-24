function [ID1,ID2,Type,Gender,Birth,Date,Eye,RGC_HFA,...
    SUP,INF,RGC_OCT,RGCNORM,SUP1,INF1,Area_OD,cpRNFL,...
    wRGC,CSFI_rate,InvCSFI_rate,MD_302,MD_242,...
    Mdsup_24,Mdinf_24,age,Axial_length,VFI_rate,...
    Corneal_thickness,SE,VA] = DataRead(filename, 2, endRow)

%IMPORTFILE テキスト ファイルから数値データを列ベクトルとしてインポートします。
%   [ID1,ID2,TYPE,GENDER,BIRTH,DATE,EYE,RGC_HFA,SUP,INF,RGC_OCT,RGCNORM,SUP1,INF1,AREA_OD,CPRNFL,WRGC,CSFI_RATE,INVCSFI_RATE,MD_302,MD_242,MDSUP_24,MDINF_24,AGE,AXIAL_LENGTH,VFI_RATE,CORNEAL_THICKNESS,SE,VA]
%   = IMPORTFILE(FILENAME) 既定の選択については テキスト ファイル FILENAME からデータを読み取ります。
%
%   [ID1,ID2,TYPE,GENDER,BIRTH,DATE,EYE,RGC_HFA,SUP,INF,RGC_OCT,RGCNORM,SUP1,INF1,AREA_OD,CPRNFL,WRGC,CSFI_RATE,INVCSFI_RATE,MD_302,MD_242,MDSUP_24,MDINF_24,AGE,AXIAL_LENGTH,VFI_RATE,CORNEAL_THICKNESS,SE,VA]
%   = IMPORTFILE(FILENAME, STARTROW, ENDROW) テキスト ファイル FILENAME の STARTROW
%   行から ENDROW 行までのデータを読み取ります。
%
% Example:
%   [ID1,ID2,Type,Gender,Birth,Date,Eye,RGC_HFA,SUP,INF,RGC_OCT,RGCNORM,SUP1,INF1,Area_OD,cpRNFL,wRGC,CSFI_rate,InvCSFI_rate,MD_302,MD_242,Mdsup_24,Mdinf_24,age,Axial_length,VFI_rate,Corneal_thickness,SE,VA] = importfile('CSFI_data.csv',2, 745);
%
%    TEXTSCAN も参照してください。

% MATLAB による自動生成 2017/01/24 14:29:53

%% 変数を初期化します。
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% データの列を文字列として読み取る:
% 詳細は TEXTSCAN のドキュメンテーションを参照してください。
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% テキスト ファイルを開きます。
fileID = fopen(filename,'r');

%% データの列を書式文字列に従って読み取ります。
% この呼び出しは、このコードの生成に使用されたファイルの構造に基づいています。別のファイルでエラーが発生する場合は、インポート
% ツールからコードの再生成を試みてください。
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% テキスト ファイルを閉じます。
fclose(fileID);

%% 数値文字列を含む列の内容を数値に変換します。
% 非数値文字列を NaN で置き換えます。
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29]
    % 入力セル配列の文字列を数値に変換します。非数値文字列が NaN で置き換えられました。
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % 数値でない接頭辞と接尾辞を検出して削除する正規表現を作成します。
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % 桁区切り以外の場所でコンマが検出されました。
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % 数値文字列を数値に変換します。
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end


%% データを数値列とセル列に分割します。
rawNumericColumns = raw(:, [1,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29]);
rawCellColumns = raw(:, [2,3,4,5,6,7]);


%% 非数値セルを次の値で置き換え NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % 非数値セルを検索
rawNumericColumns(R) = {NaN}; % 非数値セルを置き換え

%% インポートした配列を列変数名に割り当てます
ID1 = cell2mat(rawNumericColumns(:, 1));
ID2 = rawCellColumns(:, 1);
Type = rawCellColumns(:, 2);
Gender = rawCellColumns(:, 3);
Birth = rawCellColumns(:, 4);
Date = rawCellColumns(:, 5);
Eye = rawCellColumns(:, 6);
RGC_HFA = cell2mat(rawNumericColumns(:, 2));
SUP = cell2mat(rawNumericColumns(:, 3));
INF = cell2mat(rawNumericColumns(:, 4));
RGC_OCT = cell2mat(rawNumericColumns(:, 5));
RGCNORM = cell2mat(rawNumericColumns(:, 6));
SUP1 = cell2mat(rawNumericColumns(:, 7));
INF1 = cell2mat(rawNumericColumns(:, 8));
Area_OD = cell2mat(rawNumericColumns(:, 9));
cpRNFL = cell2mat(rawNumericColumns(:, 10));
wRGC = cell2mat(rawNumericColumns(:, 11));
CSFI_rate = cell2mat(rawNumericColumns(:, 12));
InvCSFI_rate = cell2mat(rawNumericColumns(:, 13));
MD_302 = cell2mat(rawNumericColumns(:, 14));
MD_242 = cell2mat(rawNumericColumns(:, 15));
Mdsup_24 = cell2mat(rawNumericColumns(:, 16));
Mdinf_24 = cell2mat(rawNumericColumns(:, 17));
age = cell2mat(rawNumericColumns(:, 18));
Axial_length = cell2mat(rawNumericColumns(:, 19));
VFI_rate = cell2mat(rawNumericColumns(:, 20));
Corneal_thickness = cell2mat(rawNumericColumns(:, 21));
SE = cell2mat(rawNumericColumns(:, 22));
VA = cell2mat(rawNumericColumns(:, 23));


