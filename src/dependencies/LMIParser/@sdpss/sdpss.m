    classdef sdpss
    % Objects sdpss are state-space systems which can be used for optimization with Yalmip
    % Any sdpss object is a structure which consists of four fields:
    %
    % sdpss.A 
    % sdpss.B 
    % sdpss.C
    % sdpss.D
    %
    % In operations, if the first argument is of type sdpss, then all
    % others are converted to spdss for compatibility.
    %
    %   s=sdpss(A);
    %   A is sdpss then output is just A
    %   A is double or spdvar then s.D=A and all other empty
    %   A is ss then output is the sdpss version of ss
    %
    %   s=sdpss(A,B,C,D);
    %   Output is sdpss system with A,B,C,D
            
        
    % Public, tunable properties
    properties
    A=[]; 
    B=[];
    C=[];
    D=[];
    end

    % Pre-computed constants or internal states
    properties (Access = private)   

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Definition
    methods (Access = protected)
    end
    methods
        function s=sdpss(A,B,C,D);            
            if nargin==1;
                % Accept in first argument
                %   {double,sdpvar,ss}->sdpss
                %   sdpss: leave untouched
                if isa(A,'sdpss');
                    s=A;
                end;
                if isa(A,'double') | isa(A,'sdpvar') | isa(A,'lmim');                    
                    s.D=A;
                    [k,m]=size(s.D);
                    s.A=[];
                    s.B=double.empty(0,m);
                    s.C=double.empty(k,0);
                end;
                if isa(A,'ss');
                    [A,B,C,D]=ssdata(A);
                    s.A=A;
                    s.B=B;
                    s.C=C;
                    s.D=D;
                end;
            end
            if nargin==4;
                % Generates sdpss object from A,B,C,D
                if size(A,1)~=size(A,2);error('A must be square.');end;
                s.A=A;
                s.B=B;
                s.C=C;
                s.D=D;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [A,B,C,D]=ssdata(s1)
            s1=sdpss(s1);
            % Extracts A,B,C,D from sdpss object.
            A=s1.A;B=s1.B;C=s1.C;D=s1.D;
        end
        function [M,n]=ssdatam(s1)
            % Extracts M=[A B;C D] and dimension n of A from sdpss object.
            s1=sdpss(s1);
            n=size(s1.A,1);
            M=[s1.A s1.B;s1.C s1.D];
        end
        function s=sdpss2ss(s1);
            % Converts sdpss to ss object.
            % Only possible for type double.
            if ~isa(s1,'sdpss');
                error('Input must be sdpss system.')
            else
                [A,B,C,D]=ssdata(s1);
                s=ss(A,B,C,D);
            end
        end
        function s=ssvalue(s1)
            % Transforms sdpss object with sdpvar entries
            % using yalmip command 'value' into sdpss object with doubles
            s=sdpss(s1);
            s.A=value(s1.A);
            s.B=value(s1.B);
            s.C=value(s1.C);
            s.D=value(s1.D);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         function ssvar(s1)
%             s1=sdpss(s1);
%             mysee(s1.A,'A');
%             mysee(s1.B,'B');
%             mysee(s1.C,'C');
%             mysee(s1.D,'D');
%         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function see(s1)
            disp('Matrix A -------------------------------------------------------------------------------')
            mysee(s1.A);
            disp('Matrix B -------------------------------------------------------------------------------')
            mysee(s1.B);
            disp('Matrix C -------------------------------------------------------------------------------')
            mysee(s1.C);
            disp('Matrix D -------------------------------------------------------------------------------')
            mysee(s1.D);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end 