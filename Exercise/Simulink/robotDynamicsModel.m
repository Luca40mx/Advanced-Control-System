function robotDynamicsModel(block)

    % The setup method is used to set up the basic attributes of the
    % S-function such as ports, parameters, etc. Do not add any other
    % calls to the main body of the function.
    %
    setup(block);

    %endfunction

    % Function: setup ===================================================
    % Abstract:
    %   Set up the basic characteristics of the S-function block such as:
    %   - Input ports
    %   - Output ports
    %   - Dialog parameters
    %   - Options

    function setup(block)

        % Register number of ports
        block.NumInputPorts = 2; % tau, he
        block.NumOutputPorts = 2; % q, dq

        % Setup port properties to be inherited or dynamic
        block.SetPreCompInpPortInfoToDynamic;
        block.SetPreCompOutPortInfoToDynamic;

        %% INPUT:
        % tau
        block.InputPort(1).Dimensions = 3;
        block.InputPort(1).DatatypeID = 0; % double
        block.InputPort(1).Complexity = 'Real';
        block.InputPort(1).DirectFeedthrough = true;

        % external wrench, he
        block.InputPort(2).Dimensions = 6;
        block.InputPort(2).DatatypeID = 0; % double
        block.InputPort(2).Complexity = 'Real';
        block.InputPort(2).DirectFeedthrough = true;

        %% OUTPUT:
        % q
        block.OutputPort(1).Dimensions = 3;
        block.OutputPort(1).DatatypeID = 0; % double
        block.OutputPort(1).Complexity = 'Real';

        % dq
        block.OutputPort(2).Dimensions = 3;
        block.OutputPort(2).DatatypeID = 0; % double
        block.OutputPort(2).Complexity = 'Real';

        % Register parameters
        block.NumDialogPrms = 3; % robotStructure

        block.NumContStates = 6; % first 3 element for position, the other for velocity

        % Register sample times
        block.SampleTimes = [0 0]; % Continuous sample time

        block.SimStateCompliance = 'DefaultSimState';

        % Register methods
        block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup);
        block.RegBlockMethod('InitializeConditions', @InitializeConditions);
        block.RegBlockMethod('Start', @Start);
        block.RegBlockMethod('Outputs', @Outputs); % Required
        % block.RegBlockMethod('Update', @Update);
        block.RegBlockMethod('Derivatives', @Derivatives);
        block.RegBlockMethod('Terminate', @Terminate);

        %end setup

        % PostPropagationSetup:
        %   Functionality    : Setup work areas and state variables. Can
        %                      also register run-time methods here

        function DoPostPropSetup(block)
            block.NumDworks = 1;

            block.Dwork(1).Name = 'F';
            block.Dwork(1).Dimensions = 3;
            block.Dwork(1).DatatypeID = 0; % double
            block.Dwork(1).Complexity = 'Real'; % real
            block.Dwork(1).UsedAsDiscState = true;

            % InitializeConditions:
            %   Functionality    : Called at the start of simulation and if it is
            %                      present in an enabled subsystem configured to reset
            %                      states, it will be called when the enabled subsystem
            %                      restarts execution to reset the states.

            function InitializeConditions(block)

                %end InitializeConditions

                % Start:
                %   Functionality    : Called once at start of model execution. If you
                %                      have states that should be initialized once, this
                %                      is the place to do it.

                q_init = block.DialogPrm(2).Data';
                dq_init = block.DialogPrm(3).Data';
                block.ContStates.Data = [q_init; dq_init];

                function Start(~)

                    %end Start

                    % Outputs:
                    %   Functionality    : Called to generate block outputs in
                    %                      simulation step

                    function Outputs(block)

                        states = block.ContStates.Data;
                        q = states(1:3); % joint position
                        dq = states(4:6); % joint velocity

                        block.OutputPort(1).Data = q;
                        block.OutputPort(2).Data = dq;

                        %end Outputs

                        % Update:
                        %   Functionality    : Called to update discrete states
                        %                      during simulation step

                        function Derivatives(block)
                            
                            simrobot = block.DialogPrm(1).Data;
                            states = block.ContStates.Data;

                            x1 = states(1:3);
                            x2 = states(4:6);
                            tau = block.InputPort(1).Data;
                            he = block.InputPort(2).Data;


                            B = simrobot.func.Bmatrix(x1(1), x1(2), x1(3));
                            C = simrobot.func.Cmatrix(x1(1), x1(2), x1(3), x2(1), x2(2), x2(3));
                            g = simrobot.func.Gmatrix(x1(1), x1(2), x1(3));
                            J_geometric = simrobot.func.GeometricJacobian(x1(1), x1(2), x1(3));

                            dx1 = x2;
                            dx2 = B \ ((tau - J_geometric' * he) - C * x2 - g);

                            block.Derivatives.Data = [dx1; dx2];

                            %end Derivatives

                            % Terminate:
                            %   Functionality    : Called at the end of simulation for cleanup

                            function Terminate(~)

                                %end Terminate
