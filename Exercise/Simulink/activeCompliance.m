function activeCompliance(block)
    %%
    %% The setup method is used to set up the basic attributes of the
    %% S-function such as ports, parameters, etc. Do not add any other
    %% calls to the main body of the function.
    %%
    setup(block);

    %endfunction

    %% Function: setup ===================================================
    %% Abstract:
    %%   Set up the basic characteristics of the S-function block such as:
    %%   - Input ports
    %%   - Output ports
    %%   - Dialog parameters
    %%   - Options
    %%
    %%   Required         : Yes
    %%   C MEX counterpart: mdlInitializeSizes
    %%
    function setup(block)

        % Register number of ports
        block.NumInputPorts = 2;
        block.NumOutputPorts = 1;

        % Setup port properties to be inherited or dynamic
        block.SetPreCompInpPortInfoToDynamic;
        block.SetPreCompOutPortInfoToDynamic;

        % Override input port properties

        % desired position
        block.InputPort(1).Dimensions = 3;
        block.InputPort(1).DatatypeID = 0; % double
        block.InputPort(1).Complexity = 'Real';
        block.InputPort(1).DirectFeedthrough = true;

        % actual position
        block.InputPort(2).Dimensions = 3;
        block.InputPort(2).DatatypeID = 0; % double
        block.InputPort(2).Complexity = 'Real';
        block.InputPort(2).DirectFeedthrough = true;

        % Override output port properties
        block.OutputPort(1).Dimensions = [3, 3];
        block.OutputPort(1).DatatypeID = 0; % double
        block.OutputPort(1).Complexity = 'Real';

        % Register parameters
        block.NumDialogPrms = 1;

        % Register sample times
        %  [0 offset]            : Continuous sample time
        %  [positive_num offset] : Discrete sample time
        %
        %  [-1, 0]               : Inherited sample time
        %  [-2, 0]               : Variable sample time
        block.SampleTimes = [0 0];

        % Specify the block simStateCompliance. The allowed values are:
        %    'UnknownSimState', < The default setting; warn and assume DefaultSimState
        %    'DefaultSimState', < Same sim state as a built-in block
        %    'HasNoSimState',   < No sim state
        %    'CustomSimState',  < Has GetSimState and SetSimState methods
        %    'DisallowSimState' < Error out when saving or restoring the model sim state
        block.SimStateCompliance = 'DefaultSimState';

        %% -----------------------------------------------------------------
        %% The MATLAB S-function uses an internal registry for all
        %% block methods. You should register all relevant methods
        %% (optional and required) as illustrated below. You may choose
        %% any suitable name for the methods and implement these methods
        %% as local functions within the same file. See comments
        %% provided for each function for more information.
        %% -----------------------------------------------------------------

        block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup);
        block.RegBlockMethod('InitializeConditions', @InitializeConditions);
        block.RegBlockMethod('Start', @Start);
        block.RegBlockMethod('Outputs', @Outputs); % Required
        block.RegBlockMethod('Update', @Update);
        block.RegBlockMethod('Derivatives', @Derivatives);
        block.RegBlockMethod('Terminate', @Terminate);

        %end setup

        %%
        %% PostPropagationSetup:
        %%   Functionality    : Setup work areas and state variables. Can
        %%                      also register run-time methods here
        %%   Required         : No
        %%   C MEX counterpart: mdlSetWorkWidths
        %%
        function DoPostPropSetup(block)
            block.NumDworks = 1;

            block.Dwork(1).Name = 'x1';
            block.Dwork(1).Dimensions = 1;
            block.Dwork(1).DatatypeID = 0; % double
            block.Dwork(1).Complexity = 'Real'; % real
            block.Dwork(1).UsedAsDiscState = true;

            %%
            %% InitializeConditions:
            %%   Functionality    : Called at the start of simulation and if it is
            %%                      present in an enabled subsystem configured to reset
            %%                      states, it will be called when the enabled subsystem
            %%                      restarts execution to reset the states.
            %%   Required         : No
            %%   C MEX counterpart: mdlInitializeConditions
            %%
            function InitializeConditions(block)

                %end InitializeConditions

                %%
                %% Start:
                %%   Functionality    : Called once at start of model execution. If you
                %%                      have states that should be initialized once, this
                %%                      is the place to do it.
                %%   Required         : No
                %%   C MEX counterpart: mdlStart
                %%
                function Start(block)

                    block.Dwork(1).Data = 0;

                    %end Start

                    %%
                    %% Outputs:
                    %%   Functionality    : Called to generate block outputs in
                    %%                      simulation step
                    %%   Required         : Yes
                    %%   C MEX counterpart: mdlOutputs
                    %%
                    function Outputs(block)

                        % block.OutputPort(1).Data = block.Dwork(1).Data + block.InputPort(1).Data;

                        robotStructure = block.DialogPrm(1).Data;
                        jointVar = block.InputPort(2).Data;
                        A_b_ee = robotStructure.func.rotPosEE(jointVar(1), jointVar(2), jointVar(3));
                        J_analytical = robotStructure.func.AnalitycalJacobianComplete(jointVar(1), jointVar(2), jointVar(3));
                        desired_position = block.InputPort(1).Data;


                        Td = [eye(3), desired_position; 
                              0 0 0 1];

                        % end effector frame
                        Te = A_b_ee;
                        
                        % desired position in the end effector frame
                        T_D_e = [Td(1:3, 1:3)' * Te(1:3, 1:3), Td(1:3, 1:3)' * (Te(1:3, 4) - Td(1:3, 4));
                                 0 0 0 1];

                        phi = atan2(T_D_e(2, 3), T_D_e(1, 3));
                        gamma = acos(T_D_e(3, 3));
                        % psi = atan2(A_b_ee(3, 2), -A_b_ee(3, 1));

                        T_zyz = [0, -sin(phi), cos(phi) * sin(gamma);
                                 0, cos(phi), sin(phi) * sin(gamma);
                                 1, 0, cos(gamma)];

                        Ta = blkdiag(eye(3), T_zyz);

                        Jad = inv(Ta) * blkdiag(Td(1:3, 1:3)', Td(1:3, 1:3)') * J_analytical; %#ok<*MINV>

                        block.OutputPort(1).Data = Jad(1:3, 1:3);

                        %end Outputs

                        %%
                        %% Update:
                        %%   Functionality    : Called to update discrete states
                        %%                      during simulation step
                        %%   Required         : No
                        %%   C MEX counterpart: mdlUpdate
                        %%
                        function Update(block)

                            % block.Dwork(1).Data = block.InputPort(1).Data;

                            %end Update

                            %%
                            %% Derivatives:
                            %%   Functionality    : Called to update derivatives of
                            %%                      continuous states during simulation step
                            %%   Required         : No
                            %%   C MEX counterpart: mdlDerivatives
                            %%
                            function Derivatives(block)

                                %end Derivatives

                                %%
                                %% Terminate:
                                %%   Functionality    : Called at the end of simulation for cleanup
                                %%   Required         : No
                                %%   C MEX counterpart: mdlTerminate
                                %%
                                function Terminate(block)

                                    %end Terminate
