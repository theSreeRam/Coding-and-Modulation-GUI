function varargout = Sreeram_2017B2AA0803G(varargin)
% SREERAM_2017B2AA0803G MATLAB code for Sreeram_2017B2AA0803G.fig
%      SREERAM_2017B2AA0803G, by itself, creates a new SREERAM_2017B2AA0803G or raises the existing
%      singleton*.
%
%      H = SREERAM_2017B2AA0803G returns the handle to a new SREERAM_2017B2AA0803G or the handle to
%      the existing singleton*.
%
%      SREERAM_2017B2AA0803G('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SREERAM_2017B2AA0803G.M with the given input arguments.
%
%      SREERAM_2017B2AA0803G('Property','Value',...) creates a new SREERAM_2017B2AA0803G or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Sreeram_2017B2AA0803G_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Sreeram_2017B2AA0803G_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Sreeram_2017B2AA0803G

% Last Modified by GUIDE v2.5 21-Nov-2020 00:05:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Sreeram_2017B2AA0803G_OpeningFcn, ...
                   'gui_OutputFcn',  @Sreeram_2017B2AA0803G_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Sreeram_2017B2AA0803G is made visible.
function Sreeram_2017B2AA0803G_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Sreeram_2017B2AA0803G (see VARARGIN)

% Choose default command line output for Sreeram_2017B2AA0803G
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Sreeram_2017B2AA0803G wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Sreeram_2017B2AA0803G_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
if hObject.Value == 2
    set(handles.pushbutton2, 'Visible','Off');
    set(handles.pushbutton4, 'Visible','Off');
    set(handles.pushbutton5, 'Visible','Off');
    set(handles.pushbutton7, 'Visible','Off');
    set(handles.pushbutton8, 'Visible','Off');
    set(handles.pushbutton9, 'Visible','Off');
    set(handles.pushbutton10, 'Visible','Off');
    set(handles.pushbutton13, 'Visible', 'Off');
    set(handles.pushbutton14, 'Visible', 'Off');
    
    set(handles.text10, 'Visible', 'On');
    set(handles.text11, 'Visible', 'On');
    set(handles.text12, 'Visible', 'Off');
    set(handles.text13, 'Visible', 'On');
    
    set(handles.axes2, 'Visible', 'Off');
    set(handles.axes4, 'Visible', 'Off');
    set(handles.axes7, 'Visible', 'Off');
    
    N = 256 ;
    x = zeros(1,N);
    a = rand(1,N);
    for i = 1:N
        if a(i)>0.5
            x(i) = 1;
        end
    end
    disp(x);
    handles.source = x;
    set(handles.text10, 'String', num2str(x));
    fprintf("The Random Bits Generated are: \n");
    disp(x);
    handles.source = x;
    handles.source_choice = 2;

elseif hObject.Value == 3
    set(handles.pushbutton2, 'Visible','Off');
    set(handles.pushbutton4, 'Visible','Off');
    set(handles.pushbutton5, 'Visible','Off');
    set(handles.pushbutton7, 'Visible','Off');
    set(handles.pushbutton8, 'Visible','Off');
    set(handles.pushbutton9, 'Visible','Off');
    set(handles.pushbutton10, 'Visible','On');
    set(handles.pushbutton13, 'Visible', 'Off');
    set(handles.pushbutton14, 'Visible', 'Off');
    
    set(handles.text10, 'Visible', 'On');
    set(handles.text11, 'Visible', 'On');
    set(handles.text12, 'Visible', 'Off');
    set(handles.text13, 'Visible', 'On');
    
    set(handles.axes2, 'Visible', 'Off');
    set(handles.axes4, 'Visible', 'Off');
    set(handles.axes7, 'Visible', 'Off');
    
    handles.source_choice = 3;
    

elseif hObject.Value == 4
    set(handles.pushbutton2, 'Visible','On');
    set(handles.pushbutton4, 'Visible','Off');
    set(handles.pushbutton5, 'Visible','Off');
    set(handles.pushbutton7, 'Visible','Off');
    set(handles.pushbutton8, 'Visible','Off');
    set(handles.pushbutton9, 'Visible','Off');
    set(handles.pushbutton10, 'Visible','Off');
    set(handles.pushbutton13, 'Visible', 'Off');
    set(handles.pushbutton14, 'Visible', 'Off');

    set(handles.text10, 'Visible', 'Off');
    set(handles.text11, 'Visible', 'Off');
    set(handles.text12, 'Visible', 'Off');
    set(handles.text13, 'Visible', 'Off');
    
    set(handles.axes2, 'Visible', 'On');
    set(handles.axes4, 'Visible', 'On');
    set(handles.axes7, 'Visible', 'On');
    
    handles.source_choice = 4;

elseif hObject.Value == 5
    set(handles.pushbutton2, 'Visible','Off');
    set(handles.pushbutton4, 'Visible','On');
    set(handles.pushbutton5, 'Visible','On');
    set(handles.pushbutton7, 'Visible','On');
    set(handles.pushbutton8, 'Visible','On');
    set(handles.pushbutton9, 'Visible','On');
    set(handles.pushbutton10, 'Visible','Off');
    set(handles.pushbutton13, 'Visible', 'On');
    set(handles.pushbutton14, 'Visible', 'On');
    
    set(handles.text10, 'Visible', 'Off');
    set(handles.text11, 'Visible', 'Off');
    set(handles.text12, 'Visible', 'On');
    set(handles.text13, 'Visible', 'Off');
    
    set(handles.axes2, 'Visible', 'Off');
    set(handles.axes4, 'Visible', 'Off');
    set(handles.axes7, 'Visible', 'Off');
    
    handles.source_choice = 5;

    
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
slider_value = get(handles.slider1,'Value');
handles.awgn_SNR = slider_value;
set(handles.text7, 'String', num2str(slider_value));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
handles.awgn_SNR = 10;
guidata(hObject, handles);


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.mod_choice == 1 & handles.code_choice == 1

disp("16 FSK with Hamming Code Started, Please Wait");    

M = 16;         
k = log2(M);   
EbNo = 5;    
Fs = 340;      
nsamp = 128;     
freqsep = 20;  
ebnoVec = -5:7;

z_match = handles.awgn_SNR + 6;

data1 = handles.coded;
data2 = handles.source;

original_length_data2 = length(data2);
original_length_data1 = length(data1);
original_source_fsk = data2';

fskMod = comm.FSKModulator;
fskMod.ModulationOrder = M;
fskMod.FrequencySeparation = freqsep;  
fskMod.SamplesPerSymbol = nsamp;  
fskMod.SymbolMapping = 'Gray';  
fskMod.SymbolRate = Fs; 

fskDemod = comm.FSKDemodulator;
fskDemod.ModulationOrder = M;
fskDemod.FrequencySeparation = freqsep;  
fskDemod.SamplesPerSymbol = nsamp;  
fskDemod.SymbolMapping = 'Gray';  
fskDemod.SymbolRate = Fs; 

err = comm.ErrorRate;
 
if mod(length(data1),4) ~= 0
   for a = 3:-1:mod(length(data1),4)
      data1 = [data1,0];
   end
end
data1 = data1';

if mod(length(data2),4) ~= 0
   for a = 3:-1:mod(length(data2),4)
      data2 = [data2,0] ;
   end
end
data2 = data2';

release(fskMod);
modData1 = fskMod(data1);
release(fskMod);
modData2 = fskMod(data2);

awgnChannel = comm.AWGNChannel('BitsPerSymbol',log2(16));
errorRate = comm.ErrorRate;

ber1 = zeros(size(ebnoVec));
ber2 = zeros(size(ebnoVec));

 n = 7; %size of codeword
 k = 4; %size of message
 A = [ 1 1 1; 1 1 0; 1 0 1; 0 1 1 ]; %Parity submatrix            
 G = [ eye(k) A ]; %Generator matrix
 H = [ A' eye(n-k) ]; %Parity-check matrix

s = zeros([1 original_length_data2]);
et = zeros([1 original_length_data1]);

for z = 1:length(ebnoVec)
    reset(errorRate)
    awgnChannel.EbNo = ebnoVec(z);
        
    rxSig1 = awgnChannel(modData1);
    rxSig2 = awgnChannel(modData2);
        
    release(fskDemod);
    rxData1 = fskDemod(rxSig1);
    release(fskDemod);
    rxData2 = fskDemod(rxSig2);
    
     for h=1:1:original_length_data2
          s(h)= rxData2(h);
          if z == z_match
              handles.dem2 = s';
              guidata(hObject,handles);

          end
     end
     
     for h=1:1:original_length_data1
          et(h)= rxData1(h);
     end
     
    errVec2 = errorRate(original_source_fsk,s');
    ber2(z) = errVec2(1);

    y = et;
    len = length(y)-6;
    x = zeros([1 original_length_data2]);
    const = 1;
    for i=1:7:len
        codewrd = y(1,i:i+6);
        syndrome = mod(codewrd * H',2);

        correctedcode = codewrd;
        msg_decoded=correctedcode;
        msg_decoded=msg_decoded(1:4);
        x(1, const:const + 3) = msg_decoded;
        const = const + 4;
     end
    t=zeros([1 original_length_data2]);
    for h=1:1:original_length_data2
         t(h)= x(h);
    end
    
    [numerr1(z),ber1(z)] = biterr(original_source_fsk,t');
    if z == z_match
        handles.dem1 = t';
        guidata(hObject,handles);
    end
    
end

x = modData1;
Fs = 1000;
t = 0:1/Fs:1-1/Fs;
figure(2)
plot(psd(spectrum.periodogram,x,'Fs',Fs,'NFFT',length(x)))
    

t1=0:length(x)-1;
figure(3)
plot(t1,real(x));
title('Time domain');

for i = 1:length(ebnoVec)
    BER_theory(i) = berawgn(ebnoVec(i),'fsk',M,'noncoherent');
end
    
figure(1)
semilogy(ebnoVec,[ber1; ber2; BER_theory])
legend('Coded Simulation BER', 'UnCoded Simulation BER','Theoretical BER','location','ne')
xlabel('Eb/No (dB)')
ylabel('BER')
title('16 FSK Modulation - BER Curves')
grid on

guidata(hObject, handles);
disp("16 FSK Modulation with Hamming Code done"); 
elseif handles.mod_choice == 1 & handles.code_choice == 2
%-------------------------16-FSK WITH BCH--------------------------------%
disp("16 FSK with BCH Started, Please Wait");
M = 16;
k = log2(M);   
EbNo = 5;      
Fs = 340;       
nsamp = 128;     
freqsep = 20;  
ebnoVec = -4:12;

z_match = handles.awgn_SNR + 5;

data1 = handles.coded;
data2 = handles.source;

original_length_data2 = length(data2);
original_length_data1 = length(data1);
original_source_fsk = data2';

fskMod = comm.FSKModulator;
fskMod.ModulationOrder = M;
fskMod.FrequencySeparation = freqsep;  
fskMod.SamplesPerSymbol = nsamp;  
fskMod.SymbolMapping = 'Gray';  
fskMod.SymbolRate = Fs; 

fskDemod = comm.FSKDemodulator;
fskDemod.ModulationOrder = M;
fskDemod.FrequencySeparation = freqsep;  
fskDemod.SamplesPerSymbol = nsamp;  
fskDemod.SymbolMapping = 'Gray';  
fskDemod.SymbolRate = Fs; 

err = comm.ErrorRate;
 
if mod(length(data1),4) ~= 0
   for a = 3:-1:mod(length(data1),4)
      data1 = [data1,0];
   end
end
data1 = data1';

if mod(length(data2),4) ~= 0
   for a = 3:-1:mod(length(data2),4)
      data2 = [data2,0] ;
   end
end
data2 = data2';
    release(fskMod);
    modData1 = fskMod(data1);
    release(fskMod);
    modData2 = fskMod(data2);
awgnChannel = comm.AWGNChannel('BitsPerSymbol',log2(16));
errorRate = comm.ErrorRate;

ber1 = zeros(size(ebnoVec));
ber2 = zeros(size(ebnoVec));

s = zeros([1 original_length_data2]);
et = zeros([1 original_length_data1]);
     
for z = 1:length(ebnoVec)
    reset(errorRate)
    awgnChannel.EbNo = ebnoVec(z);
    rxSig1 = awgnChannel(modData1);
    rxSig2 = awgnChannel(modData2);
    release(fskDemod);
    rxData1 = fskDemod(rxSig1);
    release(fskDemod);
    rxData2 = fskDemod(rxSig2);
    
    for h=1:1:original_length_data2
          s(h)= rxData2(h);
    end
     if z == z_match;
         handles.dem2 = s';
         guidata(hObject,handles);
     end
     for h=1:1:original_length_data1
          
          if (rxData1(h) ~= 1 && rxData1(h) ~= 0)
              et(h)=1;
          else
              et(h)= rxData1(h);
          end
     end
    errVec2 = errorRate(original_source_fsk,s');
    ber2(z) = errVec2(1);
   
    binary_db = reshape(et,127,[]);%15*127
    binary_db = binary_db';
    binary_db = gf(binary_db);
    [newmsg,err,ccode] = bchdec(binary_db,127,64);
    x = newmsg.x;
    DecOutput = reshape(x,[],1);
      
     t=zeros([1 original_length_data2]);
     for h=1:1:original_length_data2
          t(h)= DecOutput(h);
     end

     [numerr1(z),ber1(z)] = biterr(original_source_fsk,t');
     
     if z == z_match
         handles.dem1 = t';
         guidata(hObject,handles);
     end
    
end
x = modData1;
Fs = 1000;
t = 0:1/Fs:1-1/Fs;
figure(2)
plot(psd(spectrum.periodogram,x,'Fs',Fs,'NFFT',length(x)))
    
t1=0:length(x)-1;
figure(3)
plot(t1,real(x));
title('Time domain');

for i = 1:length(ebnoVec)
    BER_theory(i) = berawgn(ebnoVec(i),'fsk',M,'noncoherent');
end
    
figure(1)
semilogy(ebnoVec,[ber1; ber2; BER_theory])
legend('Coded Simulation BER', 'Uncoded Simulation BER','Theoretical BER','location','ne')
xlabel('Eb/No (dB)')
ylabel('BER')
title('16 FSK Modulation - BER Curves')
grid on

disp("16 FSK with BCH Done");
guidata(hObject, handles);

elseif handles.mod_choice==2 & handles.code_choice == 1
%-----------------------16-QAM WITH HAMMING------------------------------%%
disp("16 QAM with Hamming Started, Please Wait");
%%
x1 = handles.coded;
x2 = handles.source;
x1 = x1';
x2 = x2';

Modulator = comm.GeneralQAMModulator; 
Demodulator = comm.GeneralQAMDemodulator;  
modulatedata1 = Modulator(x1);  
modulatedata2 = Modulator(x2);  

E_bN_0= -5:10; 

z_match = 6 + handles.awgn_SNR;
awgnChannel = comm.AWGNChannel('BitsPerSymbol',log2(16));

for j=1:length(E_bN_0)

    awgnChannel.EbNo = E_bN_0(j);
    
    receivedSig2(:,j)= awgnChannel(modulatedata2);       
    Demodulatedata2 = Demodulator(receivedSig2(:,j)); 
    [numerr2(j),errrate2(j)] = biterr(x2,Demodulatedata2); 
    
    if j == z_match
        handles.dem2 = Demodulatedata2;
        guidata(hObject,handles);
    end
    
end

 n = 7; 
 k = 4; 
 A = [ 1 1 1; 1 1 0; 1 0 1; 0 1 1 ];           
 G = [ eye(k) A ]; 
 H = [ A' eye(n-k) ]; 
for j=1:length(E_bN_0)     
    
    awgnChannel.EbNo = E_bN_0(j);
    
    receivedSig1(:,j) = awgnChannel(modulatedata1);
    Demodulatedata1 = step(Demodulator,receivedSig1(:,j)); 
    y = Demodulatedata1';
    len = length(y)-6;
    x = zeros([1 length(x2)]);
    const = 1;
    for i=1:7:len
        codewrd = y(1,i:i+6);
        syndrome = mod(codewrd * H',2);

        correctedcode = codewrd;
        msg_decoded=correctedcode;
        msg_decoded=msg_decoded(1:4);
        x(1, const:const + 3) = msg_decoded;
        const = const + 4;
    end
    
    t=zeros([1 length(x2)]);
    for h=1:1:length(x2)
         t(h)= x(h);
    end
    
    [numerr1(j),errrate1(j)] = biterr(x2,t'); 

    if j == z_match
        handles.dem1 = t';
        guidata(hObject,handles);
    end

end

tber = berawgn(E_bN_0,'qam',16,'nondiff');   
figure(1)
semilogy(E_bN_0,tber,'mx-','linewidth',2) 
hold on
semilogy(E_bN_0,errrate1,'k-','linewidth',2)  
semilogy(E_bN_0,errrate2,'b-','linewidth',2)  
axis([-5 11 10^-5 1])
legend ('Theoretical BER','Coded Simulation BER', 'Uncoded Simulation BER');
xlabel('E_b/N_0 ')
ylabel('Bit Error Rate')
title('BER Plots for 16-QAM')
grid on

Fs = 1000;
t = 0:1/Fs:1-1/Fs;
x = modulatedata1;
figure(2)
plot(psd(spectrum.periodogram,x,'Fs',Fs,'NFFT',length(x)))

t1=0:length(x)-1;
figure(3)
plot(t1,real(x));
title('Time domain');

M = 16;                       
x = (0:15)';                   
y1 = qammod(x,16,'bin');       
scatterplot(y1)
text(real(y1)+0.1, imag(y1), dec2bin(x))
title('16-QAM : Constellation Diagram')
axis([-4 4 -4 4]);

disp("16 QAM with Hamming Code Done");
guidata(hObject, handles);

elseif handles.mod_choice==2 & handles.mod_choice==2
 %%
%-----------------------16-QAM WITH BCH----------------------------------%
disp("16 QAM with BCH Started, Please Wait");
%%
x1 = handles.coded;
x2 = handles.source;

original_length_data2 = length(x2);
original_length_data1 = length(x1);
original_source_fsk = x2';

x1 = x1';
x2 = x2';

Modulator = comm.GeneralQAMModulator; 
Demodulator = comm.GeneralQAMDemodulator;   

ebnoVec = -5:10; 

z_match = 6 + handles.awgn_SNR;

ber1 = zeros(size(ebnoVec));
ber2 = zeros(size(ebnoVec));

awgnChannel = comm.AWGNChannel('BitsPerSymbol',log2(16));
errorRate = comm.ErrorRate;

s = zeros([1 original_length_data2]);
et = zeros([1 original_length_data1]);

for j=1:length(ebnoVec)
    disp("iteration :");
    j
    reset(errorRate)
    errVec1 = [0 0 0];
    errVec2 = [0 0 0];
    awgnChannel.EbNo = ebnoVec(j);
    
    while errVec2(2) < 200 && errVec2(3) < 1e7
    
    modData1 = Modulator(x1);
    modData2 = Modulator(x2);
        
    rxSig1 = awgnChannel(modData1);
    rxSig2 = awgnChannel(modData2);
       
    rxData1 = Demodulator(rxSig1);
    rxData2 = Demodulator(rxSig2);
    
    for h=1:1:original_length_data2
          s(h)= rxData2(h);
          if j == 8
              handles.dem2 = s';
              guidata(hObject,handles);
          end
    end
     
     for h=1:1:original_length_data1
          
          if (rxData1(h) ~= 1 && rxData1(h) ~= 0)
              et(h)=1;
          else
              et(h)= rxData1(h);
          end
     end
     
     errVec2 = errorRate(original_source_fsk,s');
     ber2(j) = errVec2(1);
        
   
    binary_db = reshape(et,127,[]);%15*127
    binary_db = binary_db';
    binary_db = gf(binary_db);
    [newmsg,err,ccode] = bchdec(binary_db,127,64);
    x = newmsg.x;
    DecOutput = reshape(x,[],1);
    
    t=zeros([1 original_length_data2]);
     for h=1:1:original_length_data2
          t(h)= DecOutput(h);
     end
    
    errVec1 = errorRate(original_source_fsk,t');
    ber1(j) = errVec1(1);
    
    if j == z_match
        handles.dem1 = t';
        guidata(hObject,handles);
    end

    end

end


tber = berawgn(ebnoVec,'qam',16,'nondiff');   % Theoretical BER of 16 QAM in AWGN Channel
figure(1)
semilogy(ebnoVec,tber,'mx-','linewidth',2) %Plot Theoretical BER in AWGN
hold on
semilogy(ebnoVec,ber1,'k-','linewidth',2)  %Plot BER of Coded Simulated Data
semilogy(ebnoVec,ber2,'b-','linewidth',2)  %Plot BER of Uncoded Simulated Data
axis([-7 14 10^-5 1])
legend ('Theoretical BER','Coded Simulation BER', 'Uncoded Simulation BER');
xlabel('E_b/N_0 ')
ylabel('Bit Error Rate')
title('BER Plots for 16-QAM')
grid on

Fs = 1000;
t = 0:1/Fs:1-1/Fs;
x = modData1;
figure(2)
plot(psd(spectrum.periodogram,x,'Fs',Fs,'NFFT',length(x)))


%-----Plot Time Domain Graph of modulated signal, 16-QAM-------------%

t1=0:length(x)-1;
figure(3)
plot(t1,real(x));
title('Time domain plot');

%-------------------16-QAM Constellation Diagram-------------------%

  M = 16;                      
  x = (0:15)';                 
 y1 = qammod(x,16,'bin');       
scatterplot(y1)
text(real(y1)+0.1, imag(y1), dec2bin(x))
title('16-QAM : Constellation Diagram')
axis([-4 4 -4 4]);

disp("16 QAM with BCH Done.");
guidata(hObject, handles);

elseif handles.mod_choice==3 & handles.code_choice==1
%%
%---------------------------16-PSK WITH HAMMING--------------------------%
disp("16 PSK with HAMMING Started, Please Wait");

custMap = [0 2 4 6 8 10 12 14 15 13 11 9 7 5 3 1];

pskModulator = comm.PSKModulator(16,'BitInput',true,'SymbolMapping','Custom', 'CustomSymbolMapping',custMap);
pskDemodulator = comm.PSKDemodulator(16,'BitOutput',true,'SymbolMapping','Custom','CustomSymbolMapping',custMap);

awgnChannel = comm.AWGNChannel('BitsPerSymbol',log2(16));
errorRate = comm.ErrorRate;

ebnoVec = 1:15;
ber1 = zeros(size(ebnoVec));
ber2 = zeros(size(ebnoVec));

z_match = handles.awgn_SNR;

data1 = handles.coded;
    if mod(length(data1),4) ~= 0
        for a = 3:-1:mod(length(data1),4)
            data1 = [data1,0];
        end
    end
data1 = data1';
      
data2 = handles.source;
    if mod(length(data2),4) ~= 0
        for a = 3:-1:mod(length(data2),4)
            data2 = [data2,0] ;
        end
    end
data2 = data2';

n = 7; 
k = 4;
A = [ 1 1 1; 1 1 0; 1 0 1; 0 1 1 ];        
G = [ eye(k) A ]; 
H = [ A' eye(n-k) ]; 

for z = 1:length(ebnoVec)
    
    reset(errorRate)
    errVec1 = [0 0 0];
    errVec2 = [0 0 0];
    awgnChannel.EbNo = ebnoVec(z);
    
    while errVec2(2) < 200 && errVec2(3) < 1e7
               
        % Modulate the binary data
        modData1 = pskModulator(data1);
        modData2 = pskModulator(data2);
        
        % Pass the modulated data through the AWGN channel
        rxSig1 = awgnChannel(modData1);
        rxSig2 = awgnChannel(modData2);
        
        % Demodulate the received signal
        rxData1 = pskDemodulator(rxSig1);
        rxData2 = pskDemodulator(rxSig2);
        
        % Collect the error statistics
        errVec2 = errorRate(data2,rxData2);
        ber2(z) = errVec2(1);
        
   
    y = rxData1';
    len = length(y)-6;
    x = zeros([1 length(data2)]);
    const = 1;
    for i=1:7:len
        codewrd = y(1,i:i+6);
        syndrome = mod(codewrd * H',2);

        correctedcode = codewrd;
        msg_decoded=correctedcode;
        msg_decoded=msg_decoded(1:4);
        x(1, const:const + 3) = msg_decoded;
        const = const + 4;
    end
    
    t=zeros([1 length(data2)]);
    for h=1:1:length(data2)
         t(h)= x(h);
    end
    
    [numerr1(z),ber1(z)] = biterr(data2,t'); 
    
    
    if z == z_match
        handles.dem1 = t';
        handles.dem2 = rxData2;
        guidata(hObject, handles);
    end
    
    end
    
end

Fs = 10000;
t = 0:1/Fs:1-1/Fs;
x = modData1;
figure(2)
plot(psd(spectrum.periodogram,x,'Fs',Fs,'NFFT',length(x)))


t1=0:length(x)-1;
figure(3)
plot(t1,real(x));
title('Time domain');

%----------------BER PLOTS : 16-PSK-------------------------------%

berTheory = berawgn(ebnoVec,'psk',16,'nondiff');

figure(1)
semilogy(ebnoVec,[ber1; ber2; berTheory])
legend('Coded Simulation BER','Uncoded Simulation BER','Theoretical BER')
xlabel('Eb/No (dB)')
ylabel('BER')
title('16 PSK Modulation - BER curves') 
grid on

M = 16;             
phOffset = 0;       %
symMap = 'binary';  % Symbol mapping (either 'binary' or 'gray')
%Construct the modulator object.

pskModulator = comm.PSKModulator(M,phOffset,'SymbolMapping',symMap);

constellation(pskModulator);
title('Constellation of 16-PSK modulation')
guidata(hObject, handles);


elseif handles.mod_choice == 3 & handles.code_choice == 2
%%
%---------------------------16-PSK WITH BCH------------------------------%
disp("16 PSK with BCH Started, Please Wait");

data1 = handles.coded;
data2 = handles.source;

original_length_data2 = length(data2);
original_length_data1 = length(data1);
original_source_fsk = data2';

custMap = [0 2 4 6 8 10 12 14 15 13 11 9 7 5 3 1];

pskModulator = comm.PSKModulator(16,'BitInput',true,'SymbolMapping','Custom', 'CustomSymbolMapping',custMap);
pskDemodulator = comm.PSKDemodulator(16,'BitOutput',true,'SymbolMapping','Custom','CustomSymbolMapping',custMap);

awgnChannel = comm.AWGNChannel('BitsPerSymbol',log2(16));
errorRate = comm.ErrorRate;

ebnoVec = 1:15;
z_match = handles.awgn_SNR;

ber1 = zeros(size(ebnoVec));
ber2 = zeros(size(ebnoVec));

    if mod(length(data1),4) ~= 0
        for a = 3:-1:mod(length(data1),4)
            data1 = [data1,0];
        end
    end
data1 = data1';
      
    if mod(length(data2),4) ~= 0
        for a = 3:-1:mod(length(data2),4)
            data2 = [data2,0] ;
        end
    end
data2 = data2';

s = zeros([1 original_length_data2]);
et = zeros([1 original_length_data1]);

for z = 1:length(ebnoVec)
    % Reset the error counter for each Eb/No value
    reset(errorRate)
    % Reset the array used to collect the error statistics
    errVec1 = [0 0 0];
    errVec2 = [0 0 0];
    % Set the channel Eb/No
    awgnChannel.EbNo = ebnoVec(z);
    
    while errVec2(2) < 200 && errVec2(3) < 1e7
               
        % Modulate the binary data
        modData1 = pskModulator(data1);
        modData2 = pskModulator(data2);
        
        % Pass the modulated data through the AWGN channel
        rxSig1 = awgnChannel(modData1);
        rxSig2 = awgnChannel(modData2);
        
        % Demodulate the received signal
        rxData1 = pskDemodulator(rxSig1);
        rxData2 = pskDemodulator(rxSig2);
                
        for h=1:1:original_length_data2
          s(h)= rxData2(h);
          if z==8
              handles.dem2 = s';
              guidata(hObject,handles);
          end
        end
     
     for h=1:1:original_length_data1
          
          if (rxData1(h) ~= 1 && rxData1(h) ~= 0)
              et(h)=1;
          else
              et(h)= rxData1(h);
          end
     end
     
        % Collect the error statistics
        errVec2 = errorRate(original_source_fsk,s');
        ber2(z) = errVec2(1);
        
    %----------------------Decoding------------------------------%
   
    binary_db = reshape(et,127,[]);%15*127
    binary_db = binary_db';
    binary_db = gf(binary_db);
    [newmsg,err,ccode] = bchdec(binary_db,127,64);
    x = newmsg.x;
    DecOutput = reshape(x,[],1);
      
     t=zeros([1 original_length_data2]);
     for h=1:1:original_length_data2
          t(h)= DecOutput(h);
     end

     [numerr1(z),ber1(z)] = biterr(original_source_fsk,t'); %For Coded case -> Find number of Error Bits and Error Rate
    
    if z == 8
        handles.dem1 = t';  
        guidata(hObject,handles);
    end
    
    end
    
end

%-----------------Plot PSD of modulated signal, 16-PSK-----------------%
Fs = 10000;
t = 0:1/Fs:1-1/Fs;
x = modData1;
figure(3)
plot(psd(spectrum.periodogram,x,'Fs',Fs,'NFFT',length(x)))

%-----Plot Time Domain Graph of modulated signal, 16-PSK-------------%

t1=0:length(x)-1;
figure(2)
plot(t1,real(x));
title('Time domain plot');

%----------------BER PLOTS : 16-PSK-------------------------------%

berTheory = berawgn(ebnoVec,'psk',16,'nondiff');

figure(1)
semilogy(ebnoVec,[ber1; ber2; berTheory])
legend('Coded Simulation BER','Uncoded Simulation BER','Theoretical BER')
xlabel('Eb/No (dB)')
ylabel('BER')
title('16 PSK Modulation - BER curves') 
grid on

%-------------------16-PSK CONSTELLATION---------------------%
M = 16;             % Modulation alphabet size
phOffset = 0;       % Phase offset
symMap = 'binary';  % Symbol mapping (either 'binary' or 'gray')
%Construct the modulator object.

pskModulator = comm.PSKModulator(M,phOffset,'SymbolMapping',symMap);

constellation(pskModulator);
title('Constellation of 16-PSK modulation')
disp("16 PSK with BCH Done");

guidata(hObject, handles);


elseif handles.mod_choice==4 & handles.code_choice == 1
%%
%--------------------32-QAM WITH HAMMING---------------------------------%
disp("32 QAM with Hamming Started, Please Wait");

x1 = handles.coded;
x2 = handles.source;
x1 = x1';
x2 = x2';

%------------------Modulation-------------------------%

Modulator = comm.GeneralQAMModulator; 
Demodulator = comm.GeneralQAMDemodulator;   
modulatedata1 = Modulator(x1);  
modulatedata2 = Modulator(x2);  

E_bN_0= 1:15;  
z_match =  handles.awgn_SNR;
awgnChannel = comm.AWGNChannel('BitsPerSymbol',log2(32));


for j=1:length(E_bN_0)

    awgnChannel.EbNo = E_bN_0(j);
    
    receivedSig2(:,j)= awgnChannel(modulatedata2);      
    Demodulatedata2 = Demodulator(receivedSig2(:,j)); 
    [numerr2(j),errrate2(j)] = biterr(x2,Demodulatedata2); 
    
    if j == z_match
        handles.dem2 = Demodulatedata2;
        guidata(hObject,handles);
    end
    
end

 n = 7; 
 k = 4; 
 A = [ 1 1 1; 1 1 0; 1 0 1; 0 1 1 ];        
 G = [ eye(k) A ]; 
 H = [ A' eye(n-k) ];

for j=1:length(E_bN_0)    
    
    awgnChannel.EbNo = E_bN_0(j);
    
    receivedSig1(:,j) = awgnChannel(modulatedata1);  
    Demodulatedata1 = step(Demodulator,receivedSig1(:,j)); 
    y = Demodulatedata1';
    len = length(y)-6;
    x = zeros([1 length(x2)]);
    const = 1;
    for i=1:7:len
        codewrd = y(1,i:i+6);
        syndrome = mod(codewrd * H',2);

        correctedcode = codewrd;
        msg_decoded=correctedcode;
        msg_decoded=msg_decoded(1:4);
        x(1, const:const + 3) = msg_decoded;
        const = const + 4;
    end
    
    t=zeros([1 length(x2)]);
    for h=1:1:length(x2)
         t(h)= x(h);
    end
    
    [numerr1(j),errrate1(j)] = biterr(x2,t'); 

    if j == z_match
        handles.dem1 = t';
        guidata(hObject,handles);
    end

end

tber = berawgn(E_bN_0,'qam',32,'nondiff');   % Theoretical BER of 32 QAM in AWGN Channel 
figure(1)
semilogy(E_bN_0,tber,'mx-','linewidth',2) %Plot Theoretical BER in AWGN
hold on
semilogy(E_bN_0,errrate1,'k-','linewidth',2)  %Plot BER of Coded Simulated Data
semilogy(E_bN_0,errrate2,'b-','linewidth',2)  %Plot BER of Uncoded Simulated Data
axis([-7 15 10^-8 1])
legend ('Theoretical BER','Coded Simulation BER', 'Uncoded Simulation BER');
xlabel('E_b/N_0 ')
ylabel('Bit Error Rate')
title('32 QAM Modulation - BER curves')
grid on

Fs = 1000;
t = 0:1/Fs:1-1/Fs;
x = modulatedata1;
figure(2)
plot(psd(spectrum.periodogram,x,'Fs',Fs,'NFFT',length(x)))

t1=0:length(x)-1;
figure(3)
plot(t1,real(x));
title('Time domain plot');
disp("32 QAM with Hamming DONE");
guidata(hObject, handles);

elseif handles.mod_choice==4 & handles.code_choice ==2
%%
%-----------------------32-QAM WITH BCH----------------------------------%
disp("32 QAM with BCH Started, Please Wait");
x1 = handles.coded;
x2 = handles.source;

original_length_data2 = length(x2);
original_length_data1 = length(x1);
original_source_fsk = x2';

x1 = x1';
x2 = x2';

Modulator = comm.GeneralQAMModulator; 
Demodulator = comm.GeneralQAMDemodulator;  

ebnoVec = 1:15;

z_match = handles.awgn_SNR;

ber1 = zeros(size(ebnoVec));
ber2 = zeros(size(ebnoVec));

awgnChannel = comm.AWGNChannel('BitsPerSymbol',log2(32));
errorRate = comm.ErrorRate;

s = zeros([1 original_length_data2]);
et = zeros([1 original_length_data1]);

for j=1:length(ebnoVec)
    reset(errorRate)
    errVec1 = [0 0 0];
    errVec2 = [0 0 0];
    awgnChannel.EbNo = ebnoVec(j);
    
    while errVec2(2) < 200 && errVec2(3) < 1e7
    
    modData1 = Modulator(x1);
    modData2 = Modulator(x2);
        
    rxSig1 = awgnChannel(modData1);
    rxSig2 = awgnChannel(modData2);
    rxData1 = Demodulator(rxSig1);
    rxData2 = Demodulator(rxSig2);
    
    for h=1:1:original_length_data2
          s(h)= rxData2(h);
          if j == z_match
              handles.dem2 = s';
              guidata(hObject,handles);
          end
    end
     
     for h=1:1:original_length_data1
          
          if (rxData1(h) ~= 1 && rxData1(h) ~= 0)
              et(h)=1;
          else
              et(h)= rxData1(h);
          end
     end
     
     errVec2 = errorRate(original_source_fsk,s');
     ber2(j) = errVec2(1);
        
   
    binary_db = reshape(et,127,[]);%15*127
    binary_db = binary_db';
    binary_db = gf(binary_db);
    [newmsg,err,ccode] = bchdec(binary_db,127,64);
    x = newmsg.x;
    DecOutput = reshape(x,[],1);
    
    t=zeros([1 original_length_data2]);
     for h=1:1:original_length_data2
          t(h)= DecOutput(h);
     end
    
    errVec1 = errorRate(original_source_fsk,t');
    ber1(j) = errVec1(1);
    
    if j == z_match
        handles.dem1 = t';
        guidata(hObject,handles);
    end

    end

end


tber = berawgn(ebnoVec,'qam',32,'nondiff');  
figure(1)
semilogy(ebnoVec,tber,'mx-','linewidth',2) 
hold on
semilogy(ebnoVec,ber1,'k-','linewidth',2)  
semilogy(ebnoVec,ber2,'b-','linewidth',2)  
axis([-7 14 10^-5 1])
legend ('Theoretical BER','Coded Simulation BER', 'Uncoded Simulation BER');
xlabel('E_b/N_0 ')
ylabel('Bit Error Rate')
title('BER Plots for 32-QAM')
grid on


Fs = 1000;
t = 0:1/Fs:1-1/Fs;
x = modData1;
figure(2)
plot(psd(spectrum.periodogram,x,'Fs',Fs,'NFFT',length(x)))

t1=0:length(x)-1;
figure(3)
plot(t1,real(x));
title('Time domain');

disp("32 QAM with BCH Done");
guidata(hObject, handles);


end







%%
%-----------------------Reconfiguration Code---------------------%

if handles.source_choice == 2;
%     Random Binary Bits
    disp("The decoded message is (Coded)");
    disp((handles.dem1'));
    disp("The decoded message is (Uncoded)");
    disp(handles.dem2');
    set(handles.text11, 'String', num2str(handles.dem1'));
    set(handles.text13, 'String', num2str(handles.dem2'));

elseif handles.source_choice == 3
%     Text File input
reconstruct_uncode = handles.dem2';
bitstream_uncode = (reconstruct_uncode);
length_uncode = size(bitstream_uncode,2)
nchar_uncode = floor(length_uncode/7);
x_uncode = char(zeros(1,nchar_uncode));

k = 1;
for i = 1:7:(nchar_uncode*7)
    c = native2unicode([64 32 16 8 4 2 1]*bitstream_uncode(1,i:i+6)');
    x_uncode(k) = c;
    k = k+1;
end
disp("The text in Uncoded :");
set(handles.text13, 'String', x_uncode);
disp(x_uncode);

reconstruct_code = handles.dem1';
bitstream_code = (reconstruct_code);
length_code = size(bitstream_code,2)
nchar_code = floor(length_code/7);
x_code = char(zeros(1,nchar_code));

k = 1;
for i = 1:7:(nchar_code*7)
    c = native2unicode([64 32 16 8 4 2 1]*bitstream_code(1,i:i+6)');
    x_code(k) = c;
    k = k+1;
end
disp("The text in Coded :");
set(handles.text11, 'String', x_code);
disp(x_code);

elseif handles.source_choice == 4 
%     Image file input                             
x = handles.image_x;              
y = handles.image_y;              
z = 3;

bits_code = handles.dem1';
% size(bits_code);
% size(bits_code,2)/8;
a_code = zeros(1,size(bits_code,2)/8);
k=1;
for i = 1:8:size(a_code,2)*8
    a_code(1,k) = [128 64 32 16 8 4 2 1]*(bits_code(1,i:i+7))';
    k=k+1;
end
image_code = reshape(a_code',[x y 3]);
image_code = image_code/255;

axes(handles.axes4);
imshow(image_code);
set(handles.axes4,'Units','normalized');

bits_uncode  = handles.dem2';
size(bits_uncode)
size(bits_uncode,2)/8
a_uncode = zeros(1,size(bits_uncode,2)/8);
k=1;
for i = 1:8:size(a_uncode,2)*8
    a_uncode(1,k) = [128 64 32 16 8 4 2 1]*(bits_uncode(1,i:i+7))';
    k=k+1;
end
image_uncode = reshape(a_uncode',[x y 3]);
image_uncode = image_uncode/255;

figure('name','Image Coded');
imshow(image_code);

axes(handles.axes7);
imshow(image_uncode);
set(handles.axes4,'Units','normalized');

figure('name','Image Uncoded');
imshow(image_uncode);

elseif handles.source_choice == 5
    uncoded_aud = handles.dem2;
    original_size = handles.originalsize;
    Fs = handles.fs;
    uncoded_aud2 = reshape(uncoded_aud, [8 length(uncoded_aud)/8])';
    uncoded_aud3 = num2str(uncoded_aud2(:,1:8));
    uncoded_aud3 = uncoded_aud3(~isspace(uncoded_aud3));

    for i=1:length(uncoded_aud3)        % AWGN might have corrupted symbols to 
         if (uncoded_aud3(i) ~= '0')    % become >1 eventhough we modulated 
                 uncoded_aud3(i) = '1'; % binary data. These symbols are 
         else                           % changed back to 1.
         uncoded_aud3(i) = '0';
         end
    end

    uncoded_aud3 = reshape(uncoded_aud3, [length(uncoded_aud)/8 8]);
    uncoded_aud4 = uint8(bin2dec( char(uncoded_aud3) ));

    uncoded_aud5 = reshape( typecast( uncoded_aud4, 'single' ), original_size );
%     sound(uncoded_aud5,Fs);
    handles.uncoded_audio = uncoded_aud5;

    aud = handles.demodr1;
    original_size = handles.originalsize;
    Fs = handles.fs;
    aud2 = reshape(aud, [8 length(aud)/8])';
    aud3 = num2str(aud2(:,1:8));
    aud3 = aud3(~isspace(aud3));

    for i=1:length(aud3)        % AWGN might have corrupted symbols to 
         if (aud3(i) ~= '0')    % become >1 eventhough we modulated 
                 aud3(i) = '1'; % binary data. These symbols are 
         else                   % changed back to 1.
                aud3(i) = '0';
         end
    end

    aud3 = reshape(aud3, [length(aud)/8 8]);
    aud4 = uint8(bin2dec( char(aud3) ));

    aud5 = reshape( typecast( aud4, 'single' ), original_size );
%     sound(aud5,Fs);   %The Reconstructed Coded Audio File is played
    handles.coded_audio = aud5;

    guidata(hObject, handles);
    disp("Reconstruction Done");
end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, filepath] = uigetfile({'*.*';'*.jpeg';'*.png';'*.bmp';'*.jpg'},'File Selector');
if ~ischar(filename)
    return;  % User aborted the file selection
end

fullname = [filepath filename];

ImageFile = imread(fullname);
axes(handles.axes2)
imagesc(ImageFile);
axis off

[x y z] = size(ImageFile);
ImageFileBits = dec2bin(ImageFile);
ImageFileBits = ImageFileBits - '0';
ImageFileBits = ImageFileBits';
ImageFileBits = reshape(ImageFileBits,1,[]);
disp("Image coded to Binary successfully");
disp("Binary file size = ");
disp(size(ImageFileBits,2));
handles.image_x = x;
handles.image_y = y;
handles.source = ImageFileBits;
guidata(hObject, handles);  % Store updated handles struct in the GUI

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, filepath] = uigetfile({'*.*';'*.mp3';'*.wav'}, 'Audio Selector');
if ~ischar(filename)
    return;  % User aborted the file selection
end

fullname = [filepath filename];
disp(fullname);
set(handles.text12, 'String', fullname);
[y, Fs] = audioread(fullname);
% Fs_new = 4000;
% [Numer, Denom] = rat(Fs_new/Fs);
% y = resample(y, Numer, Denom);

 y = downsample(y, 2);

handles.audioFile = y;
handles.audioFs = Fs;

y = (y(:,1)+y(:,2))/2;

wave_binary = dec2bin( typecast( single(y(:)), 'uint8'), 8 ) - '0';
original_size = size(y);
x = uint8(bin2dec( char(wave_binary + '0') ));

temp_1 = dec2bin(x,8);
binary_temp = transpose(temp_1);
temp = binary_temp(:)-'0';
handles.source = temp';
handles.originalsize = original_size;
handles.fs = Fs;
disp("Audiofile converted to bits successfully.");
disp("Size of the binary file : ");
disp(handles.originalsize);
guidata(hObject, handles);


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, ~, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sound(handles.audioFile, handles.audioFs);


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear sound;


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sound(handles.uncoded_audio , handles.audioFs);

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*.txt'},'File Selector');
if ~ischar(filename)
    return;  % User aborted the file selection
end
file = fullfile(pathname, filename);
[fid, msg] = fopen(file, 'r');
if fid == -1
    error(msg);
end
Data = fscanf(fid, '%c');  % Or how your file is formatted
fclose(fid);
ascii = unicode2native(Data);
binary = dec2bin(ascii);
[i,l] = size(binary);
bin_trans = binary';
bits = reshape(bin_trans,1,[]);
bits = bits-'0';
disp(bits);
handles.source = bits;
set(handles.text10, 'String', Data);
disp("The Input text file is :");
disp(Data);
guidata(hObject, handles);  % Store updated handles struct in the GUI


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1

% coded = encode(handles.source, 7, 4, 'hamming/binary');
handles.code_choice = 1;
% handles.coded = coded;
% disp(coded);
% disp("Hamming code encoding Done");
% guidata(hObject, handles);


x = handles.source;
n = 7; %size of codeword
k = 4; %size of message
A = [ 1 1 1; 1 1 0; 1 0 1; 0 1 1 ]; %Parity submatrix            
G = [ eye(k) A ]; %Generator matrix
H = [ A' eye(n-k) ]; %Parity-check matrix

%appending zeroes to make input a multiple of 4
remainder = mod(length(x),4);
x_1 = x;
if (remainder ~= 0)
    for i=1:1:4-remainder
        x_1 = [x_1,0];
    end
end  

j=1;
len=length(x_1);

%ENCODING%
for i=1:4:len
    inp = x_1(1,i:i+3);
    msg = inp; %Message block vector-change to any 4 bit sequence
    code = mod(msg*G,2);%Encode message
    y(1,j:j+6) = code;
    j = j+7;
end

disp("data Encoded in hamming scheme");%coded bit stream
handles.coded = y;
guidata(hObject, handles);


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3

% HEADINGS
set(handles.text25, 'Visible', 'On');
set(handles.text26, 'Visible', 'Off');
set(handles.text4, 'Visible', 'Off');
set(handles.text27, 'Visible', 'Off');

% SCALE VALUES
set(handles.text19, 'Visible', 'On');
set(handles.text21, 'Visible', 'On');
set(handles.text20, 'Visible', 'On');

set(handles.text5, 'Visible', 'Off');
set(handles.text7, 'Visible', 'Off');
set(handles.text6, 'Visible', 'Off');

set(handles.text22, 'Visible', 'Off');
set(handles.text24, 'Visible', 'Off');
set(handles.text23, 'Visible', 'Off');

set(handles.text16, 'Visible', 'Off');
set(handles.text18, 'Visible', 'Off');
set(handles.text17, 'Visible', 'Off');

set(handles.slider3, 'Visible', 'On');
set(handles.slider1, 'Visible', 'Off');
set(handles.slider4, 'Visible', 'Off');
set(handles.slider2, 'Visible', 'Off');


% coded = handles.coded;
% y = fskmod(coded, 16, 10, 8, 340);
% handles.modulated = y;
%disp(y);
handles.mod_choice = 1;
guidata(hObject, handles);


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4

% HEADINGS
set(handles.text25, 'Visible', 'Off');
set(handles.text26, 'Visible', 'On');
set(handles.text4, 'Visible', 'Off');
set(handles.text27, 'Visible', 'Off');

% SCALE VALUES
set(handles.text19, 'Visible', 'Off');
set(handles.text21, 'Visible', 'Off');
set(handles.text20, 'Visible', 'Off');

set(handles.text5, 'Visible', 'Off');
set(handles.text7, 'Visible', 'Off');
set(handles.text6, 'Visible', 'Off');

set(handles.text22, 'Visible', 'On');
set(handles.text24, 'Visible', 'On');
set(handles.text23, 'Visible', 'On');

set(handles.text16, 'Visible', 'Off');
set(handles.text18, 'Visible', 'Off');
set(handles.text17, 'Visible', 'Off');

set(handles.slider3, 'Visible', 'Off');
set(handles.slider1, 'Visible', 'Off');
set(handles.slider4, 'Visible', 'On');
set(handles.slider2, 'Visible', 'Off');

% coded = handles.coded;
% y = qammod(coded, 16);
% handles.modulated = y;
handles.mod_choice = 2;
guidata(hObject, handles);


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
y = handles.source;
handles.code_choice = 2;
m = 7;
n = (2^m)- 1;
k=64;
data = y;
if mod(length(data),64) ~= 0
            for a = 63:-1:mod(length(data),64)
            data = [data,0] ;
            end
end
tdata = reshape(data,[ ],64);
msg = gf(tdata);
code = bchenc(msg,n,k);
d=code.x;                       %Converts GF to uint32
d_vect = dec2bin(d');
dataout = d_vect - '0';
BCHout = dataout';
handles.coded = BCHout;
handles.rearrange = d;
guidata(hObject, handles);  
disp("BCH Encoding done");

% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% HEADINGS
set(handles.text25, 'Visible', 'Off');
set(handles.text26, 'Visible', 'Off');
set(handles.text4, 'Visible', 'On');
set(handles.text27, 'Visible', 'Off');

% SCALE VALUES
set(handles.text19, 'Visible', 'Off');
set(handles.text21, 'Visible', 'Off');
set(handles.text20, 'Visible', 'Off');

set(handles.text5, 'Visible', 'On');
set(handles.text7, 'Visible', 'On');
set(handles.text6, 'Visible', 'On');

set(handles.text22, 'Visible', 'Off');
set(handles.text24, 'Visible', 'Off');
set(handles.text23, 'Visible', 'Off');

set(handles.text16, 'Visible', 'Off');
set(handles.text18, 'Visible', 'Off');
set(handles.text17, 'Visible', 'Off');

set(handles.slider3, 'Visible', 'Off');
set(handles.slider1, 'Visible', 'On');
set(handles.slider4, 'Visible', 'Off');
set(handles.slider2, 'Visible', 'Off');

% Hint: get(hObject,'Value') returns toggle state of radiobutton5
% coded = handles.coded;
% y = pskmod(coded, 16);
% handles.modulated = y;
handles.mod_choice = 3;
guidata(hObject, handles);


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% HEADINGS
set(handles.text25, 'Visible', 'Off');
set(handles.text26, 'Visible', 'Off');
set(handles.text4, 'Visible', 'Off');
set(handles.text27, 'Visible', 'On');

% SCALE VALUES
set(handles.text19, 'Visible', 'Off');
set(handles.text21, 'Visible', 'Off');
set(handles.text20, 'Visible', 'Off');

set(handles.text5, 'Visible', 'Off');
set(handles.text7, 'Visible', 'Off');
set(handles.text6, 'Visible', 'Off');

set(handles.text22, 'Visible', 'Off');
set(handles.text24, 'Visible', 'Off');
set(handles.text23, 'Visible', 'Off');

set(handles.text16, 'Visible', 'On');
set(handles.text18, 'Visible', 'On');
set(handles.text17, 'Visible', 'On');

set(handles.slider3, 'Visible', 'Off');
set(handles.slider1, 'Visible', 'Off');
set(handles.slider4, 'Visible', 'Off');
set(handles.slider2, 'Visible', 'On');

% % Hint: get(hObject,'Value') returns toggle state of radiobutton6
% coded = handles.coded;
% y = qammod(coded, 32);
% handles.modulated = y;
handles.mod_choice = 4;
guidata(hObject, handles);


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sound(handles.uncoded_audio, handles.audioFs);

% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
slider_value = get(handles.slider2,'Value');
handles.awgn_SNR = slider_value;
set(handles.text18, 'String', num2str(slider_value));
guidata(hObject, handles);

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
handles.awgn_SNR = 10;
guidata(hObject, handles);


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
slider_value = get(handles.slider3,'Value');
handles.awgn_SNR = slider_value;
set(handles.text21, 'String', num2str(slider_value));
guidata(hObject, handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
handles.awgn_SNR = 5;
guidata(hObject, handles);


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
slider_value = get(handles.slider4,'Value');
handles.awgn_SNR = slider_value;
set(handles.text24, 'String', num2str(slider_value));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

handles.awgn_SNR = 7;
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function radiobutton5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2


% --- Executes during object creation, after setting all properties.
function axes4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes4
axis off
