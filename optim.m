
function varargout = optim(varargin)
% OPTIM MATLAB code for optim.fig
%      OPTIM, by itself, creates text_a new OPTIM or raises the existing
%      singleton*.
%
%      TEXT_H = OPTIM returns the handle to text_a new OPTIM or the handle
%      to the existing singleton*.
%
%      OPTIM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OPTIM.M with the given input arguments.
%
%      OPTIM('Property','Value',...) creates text_a new OPTIM or raises the
%      existing singleton*.  Starting from the left, property value pairs
%      are applied to the GUI before optim_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property
%      application stop.  All inputs are passed to optim_OpeningFcn via
%      varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help optim

% Last Modified by GUIDE v2.5 13-Dec-2016 11:06:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @optim_OpeningFcn, ...
                   'gui_OutputFcn',  @optim_OutputFcn, ...
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
function optim_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn. hObject    handle to
% figure eventdata  reserved - to be defined in text_a future version of
% MATLAB handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to optim (see VARARGIN)

% Choose default command line output for optim
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
function varargout = optim_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT); hObject
% handle to figure eventdata  reserved - to be defined in text_a future
% version of MATLAB handles    structure with handles and user data (see
% GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_grad_Callback(hObject, eventdata, handles)
function checkbox_hessian_Callback(hObject, eventdata, handles)
function checkbox_maxfval_Callback(hObject, eventdata, handles)
function text_A_CreateFcn(hObject, eventdata, handles)
function popupmenu_f_Callback(hObject, eventdata, handles)
function checkbox_A_b_CreateFcn(hObject, eventdata, handles)
function checkbox_Aeq_beq_CreateFcn(hObject, eventdata, handles)
function checkbox_lb_ub_CreateFcn(hObject, eventdata, handles)
function text_b_CreateFcn(hObject, eventdata, handles)
function text_Aeq_CreateFcn(hObject, eventdata, handles)
function text_beq_CreateFcn(hObject, eventdata, handles)
function text_lb_CreateFcn(hObject, eventdata, handles)
function text_ub_CreateFcn(hObject, eventdata, handles)
function radiobutton_on_CreateFcn(hObject, eventdata, handles)
function radiobutton_off_CreateFcn(hObject, eventdata, handles)
function text_x0_Callback(hObject, eventdata, handles)
function checkbox_A_b_Callback(hObject, eventdata, handles)
function text_b_Callback(hObject, eventdata, handles)
function text_A_Callback(hObject, eventdata, handles)
function text_lb_Callback(hObject, eventdata, handles)
function text_ub_Callback(hObject, eventdata, handles)
function text_beq_Callback(hObject, eventdata, handles)
function text_Aeq_Callback(hObject, eventdata, handles)
function text_ceq_Callback(hObject, eventdata, handles)
function text_c_Callback(hObject, eventdata, handles)
function text_h_Callback(hObject, eventdata, handles)
function text_g_Callback(hObject, eventdata, handles)
function text_w_Callback(hObject, eventdata, handles)
function checkbox_c_ceq_Callback(hObject, eventdata, handles)
function checkbox_exitflag_Callback(hObject, eventdata, handles)
function text_f_Callback(hObject, eventdata, handles)
function x_y_text_CreateFcn(hObject, eventdata, handles)
function checkbox_exitflag_CreateFcn(hObject, eventdata, handles)
function checkbox_output_CreateFcn(hObject, eventdata, handles)
function checkbox_output_Callback(hObject, eventdata, handles)
function checkbox_lambda_CreateFcn(hObject, eventdata, handles)
function checkbox_lambda_Callback(hObject, eventdata, handles)
function checkbox_h_Callback(hObject, eventdata, handles)
function checkbox_x0_Callback(hObject, eventdata, handles)
function edit23_Callback(hObject, eventdata, handles)
function edit23_CreateFcn(hObject, eventdata, handles)
function edit25_Callback(hObject, eventdata, handles)
function edit25_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function text_ceq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function text_c_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu_solver_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu_f_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function text_f_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function text_x0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function listbox_options_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu_opt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton_optimise_Callback(hObject, eventdata, handles)

if(get(handles.func,'Value'))
    Ffun = fopen(strcat(get(handles.func_name,'String'),'.m'),'w');
    write(Ffun,'',get(handles.text_f,'String'));
    fclose(Ffun);
    fun = str2func(get(handles.func_name,'String'));
else
    fun = get(handles.text_f,'String');
    if(get(handles.popupmenu_solver,'Value') >= 3 )
        fun = str2num(fun);
    end;
end;

if(get( handles.checkbox_x0,'Value'))
x0txt = strsplit( get(handles.text_x0,'String'),' ');
for k = 1:length(x0txt)
    x0(k) = str2double(char(x0txt(k)));
end
else 
    x0 = [];
end;

if(get( handles.checkbox_A_b,'Value'))
    Atxt = get(handles.text_A,'String');
    for k = 1:size((get(handles.text_A,'String')),1)
        tmp = char((strsplit(char(Atxt(k,:)),'  ')));
        for l = 1:length(tmp)
            A(k,l) = str2double(tmp(l,:));
        end
    end
    btxt = get(handles.text_b,'String');
    for k = 1:length(get(handles.text_b,'String'))
        b(k) = str2double(char(btxt(k,:)));
    end
else
    A = [];
    b = [];
end    

if(get( handles.checkbox_Aeq_beq,'Value'))
    Aeqtxt = get(handles.text_Aeq,'String');
    for k = 1:size(get(handles.text_Aeq,'String'),1)
        tmp = char((strsplit(char(Aeqtxt(k)),' ')));
        for l = 1:length(t)
            Aeq(k,l) = str2double(tmp(l,:));
        end
    end
    beqtxt = get(handles.text_beq,'String');
    for k = 1:length(get(handles.text_beq,'String'))
        beq(k) = str2double(char(beqtxt(k)));
    end
else
    Aeq = [];
    beq = [];
end 

if(get( handles.checkbox_lb,'Value'))
    lbtxt = get(handles.text_lb,'String');
    for k = 1:length(get(handles.text_lb,'String'))
        lb(k) = str2double(char(lbtxt(k)));
    end
else
      lb = [];
end
if(get( handles.checkbox_ub,'Value'))
    ubtxt = get(handles.text_ub,'String');
    for k = 1:length(get(handles.text_ub,'String'))
        ub(k) = str2double(char(ubtxt(k)));
    end
else
  
    ub = [];
end

if(get( handles.checkbox_h,'Value'))
    Htxt = get(handles.text_h,'String');
    for k = 1:size((get(handles.text_h,'String')),1)
        tmp = char((strsplit(char(Htxt(k,:)),'  ')));
        for l = 1:length(tmp)
            H(k,l) = str2double(tmp(l,:));
        end
    end
else
    H = [];
end 

if(get( handles.checkbox_c,'Value'))
  if(get(handles.nonl,'Value'))
    Fcon = fopen(strcat(get(handles.nonlcon_name,'String'),'.m'),'w');
    write(Fcon,'',get(handles.text_c,'String'));
    fclose(Fcon);
    nonlcon = str2func(get(handles.nonlcon_name,'String'));
    else
        nonlcon = get(handles.text_c,'String');
    end; 
else
    nonlcon = [];
end
try
switch get(handles.popupmenu_solver,'Value')
    case 1 %'fmincon'
       [x,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub)
    case 2%'fminimax'
       [x, fval, maxfval, exitflag, output, lambda] = fminimax(fun, x0, A, b, Aeq, beq, lb, ub, nonlcon)
    case 3%'linprog'
       [x, fval, exitflag, output, lambda] = linprog(fun, A, b, Aeq, beq, lb, ub, x0)
    case 4%'quadprog
       [x, fval, exitflag, output, lambda] = quadprog(H, fun, A, b, Aeq, beq, lb, ub, x0)
    otherwise
end
% if(get(handles.checkbox_plot,'Value'))
%    try
%        figure;
%        fplot(fun,[x-[5,5],x+[5,5]]);
%    catch
%       errordlg('неможливо побудувати графік','plot error'); 
%    end
% end;
set(handles.x_y_text,'String',strvcat( '~~~x~~~',num2str(transpose(x)),'~~~f(x)~~~',num2str(fval)));
str = '';
if(get(handles.checkbox_output,'Value'))
    str = strvcat(str,'~~~output~~~');
    str = strvcat(str,fwrite_struct(output));
    %str = strvcat(str,strcat('output.iterations  :',num2str( output.iterations)), strcat('output.funcCount  :',num2str(output.funcCount)), strcat('output.constrviolation  :',num2str(output.constrviolation)),strcat('output.stepsize  :', num2str(output.stepsize)), strcat('output.algorithm  :',output.algorithm),strcat('output.firstorderopt  :', num2str(output.firstorderopt)), strcat('output.cgiterations  :',num2str(output.cgiterations)),' ');
end;
if(get(handles.checkbox_exitflag,'Value'))
    str = strvcat(str,'~~~evitflag~~~');
    str = strvcat(str,num2str(exitflag),' ');
end
if(get(handles.checkbox_lambda,'Value'))
    str = strvcat(str,'~~~lambda~~~');
    str = strvcat(str,fwrite_struct(lambda));
end;
if(get(handles.checkbox_maxfval,'Value'))
    str = strvcat(str,'~~~maxfval~~~');
    str = strvcat(str,num2str(maxfval),' ');
end;
if(get(handles.checkbox_grad,'Value'))
    str = strvcat(str,'~~~grad~~~');
    str = strvcat(str,num2str(grad),' ');
end;
if(get(handles.checkbox_hessian,'Value'))
    str = strvcat(str,'~~~hessian~~~');
    str = strvcat(str,num2str(hessian),' ');
end;
set(handles.res_text,'String',str);
catch
    errordlg('перевірте вхідні умови','function error');
end;
function str = fwrite_struct(s)
str = '';
names = fieldnames(s);
for k = 1:length(names)-1
    elem = s.(char(names(k)));
    if(isnumeric(elem))
        elem = num2str(elem);
    end
    str = strvcat(str,strcat(char(names(k)),'   :',elem));
end 
str = strvcat(str,' ');

function popupmenu_solver_Callback(hObject, eventdata, handles)
switch get(hObject,'Value')
    case 1 %'fmincon'
        set(handles.h,'Visible','off');
        set(handles.c_ceq,'Visible','on');

        set(handles.checkbox_maxfval,'Visible','off');
        set(handles.checkbox_grad,'Visible','on');
        set(handles.checkbox_hessian,'Visible','on');
    case 2%'fminimax'
        set(handles.h,'Visible','off');
        set(handles.c_ceq,'Visible','on');

        set(handles.checkbox_maxfval,'Visible','on');
        set(handles.checkbox_grad,'Visible','off');
        set(handles.checkbox_hessian,'Visible','off');
    case 3%'linprog'
        set(handles.h,'Visible','off');
        set(handles.c_ceq,'Visible','off');

        set(handles.checkbox_maxfval,'Visible','off');
        set(handles.checkbox_grad,'Visible','off');
        set(handles.checkbox_hessian,'Visible','off');
    case 4%'quadprog'
        set(handles.h,'Visible','on');
        set(handles.c_ceq,'Visible','off');

        set(handles.checkbox_maxfval,'Visible','off');
        set(handles.checkbox_grad,'Visible','off');
        set(handles.checkbox_hessian,'Visible','off');
    otherwise
        set(handles.pushbutton_optimise,'Visible','off');
end

function open_pushbutton_Callback(hObject, eventdata, handles)
fname = uigetfile('*.txt','name dialog window'); % получение имени и пути к файлу
if fname ~= 0
F = fopen(fname, 'r');
set(handles.text_f,'String',read(F));
set(handles.text_x0,'String',read(F));
set(handles.text_Aeq,'String',read(F));
set(handles.text_beq,'String',read(F));
set(handles.text_A,'String',read(F));
set(handles.text_b,'String',read(F));
set(handles.text_lb,'String',read(F));
set(handles.text_ub,'String',read(F));
if(get(handles.popupmenu_solver,'Value') == 4)
    set(handles.text_h,'String',read(F));
end;
% if(get(handles.popupmenu_solver,'Value') == 1)
%     set(handles.text_c,'String',read(F));
%     set(handles.text_ceq,'String',read(F));
% end;
fclose(F);
%errordlg();
end
function txt = read(F)
if ~feof(F)
txt = '';
line = fgetl(F);
while ~strcmp( line ,'#')
        txt = strvcat( txt, line);
        line = fgetl(F);
end
else
    errordlg('недостатньо даних для відображення','file error');
end

function save_pushbutton_Callback(hObject, eventdata, handles)
fname = uiputfile('*.txt','name dialog window'); % получение имени и пути к файлу
if fname ~= 0
F = fopen(fname, 'w');
write(F,'~~~f~~~',get(handles.text_f,'String'));
fprintf(F,'%-5s\n %-20s\n\n','~~~x0~~~',get(handles.text_x0,'String'));
if(get( handles.checkbox_A_b,'Value'))
    write(F,'~~~A~~~',get(handles.text_A,'String'));
    write(F,'~~~b~~~',get(handles.text_b,'String'));
end
if(get( handles.checkbox_Aeq_beq,'Value'))
    write(F,'~~~Aeq~~~',get(handles.text_Aeq,'String'));
    write(F,'~~~Beq~~~',get(handles.text_beq,'String'));
end
if(get( handles.checkbox_lb,'Value'))
    write(F,'~~~lb~~~',get(handles.text_lb,'String'));
end;
if(get( handles.checkbox_ub,'Value'))
    write(F,'~~~ub~~~',get(handles.text_ub,'String'));
end;
if(get( handles.checkbox_h,'Value'))
    write(F,'~~~h~~~',get(handles.text_h,'String'));
end;
if(get( handles.checkbox_c,'Value'))
    write(F,'~~~c~~~ceq~~~',get(handles.text_c,'String'));
end;
end
write(F,'******RESULT******',get(handles.x_y_text,'String'));
write(F,'******DETAILS******',get(handles.res_text,'String'));
fclose(F);

function  write(F,title,txt)
fprintf(F,'\n %-20s\n',title);
for k = 1:size(txt,1)
 fprintf(F,'%-20s\n',char(txt(k,:)));
end

function pushbutton_func_Callback(hObject, eventdata, handles)
if(get(hObject,'Value'))
fname = uigetfile('*.m','name dialog window');
F = fopen(fname, 'r');
str = '';
while (~feof(F))
    str = strvcat(str,fgetl(F));
end
if(strcmp( get(hObject,'String') , 'use nonlcon'))
    set(handles.nonlcon_name ,'String', fname(1:length(fname)-2));
    set(handles.text_c,'String',str);
else
    set(handles.func_name ,'String', fname(1:length(fname)-2));
    set( handles.text_f,'String',str);    
end;
fclose(F);
end


% --------------------------------------------------------------------
function mn_file_Callback(hObject, eventdata, handles)
% hObject    handle to mn_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function help_mn_Callback(hObject, eventdata, handles)
% hObject    handle to help_mn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
help

% --------------------------------------------------------------------
function open_mn_Callback(hObject, eventdata, handles)
% hObject    handle to open_mn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function save_mn_Callback(hObject, eventdata, handles)
% hObject    handle to save_mn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_plot.
function checkbox_plot_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plot
