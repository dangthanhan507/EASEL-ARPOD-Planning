classdef MPCMHE_properties
    properties
        %objective
        P1objective
        P2objective

        %decision variables
        P1optimizationVariablesCell
        P2optimizationVariablesCell
        latentVariablesCell

        %constraints
        P1constraintsCell
        P2constraintsCell
        latentConstraintsCell

        %output cell
        outputExpressionsCell

        %parameters
        parametersCell

        %optimization hyperparams
        muAggressive
        muConservative
        targetDualityGap
        maxIter
        gradientTolerance
        delta
        compilerFlags
        allowSave
        solverVerboseLevel
        verboseLevel
    end
    methods 
        function obj = init(obj, objective, p1cell, p2cell, latentcell, p1constraint, p2constraint, latentconstraint, outputcell, paramcell)

            %set optimization problem by hand

            %objective functions
            obj.P1objective = objective;
            obj.P2objective = -objective;

            %decision variable cells
            obj.P1optimizationVariablesCell = p1cell;
            obj.P2optimizationVariablesCell = p2cell;
            obj.latentVariablesCell = latentcell;

            %constraints
            obj.P1constraintsCell = p1constraint;
            obj.P2constraintsCell = p2constraint;
            obj.latentConstraintsCell = latentconstraint;

            %output and params
            obj.outputExpressionsCell = outputcell;
            obj.parametersCell = paramcell;

            %set hyperparams by default
            obj.muAggressive = 0.5;
            obj.muConservative = 0.99;
            obj.targetDualityGap = 0.1;
            obj.maxIter = 500;
            obj.gradientTolerance = 1e-3;
            obj.delta = 2;
            obj.compilerFlags = '-O1';
            obj.allowSave = true;
            obj.solverVerboseLevel = 1;
            obj.verboseLevel = 1;
        end

        function opt = setupOptimization(obj)
            classname = class2equilibriumLatentCS(...
                'classname','tmpC_target_chaser_main',...
                'P1objective',obj.P1objective,...
                'P2objective',obj.P2objective,...
                'P1optimizationVariables',obj.P1optimizationVariablesCell,...
                'P2optimizationVariables',obj.P2optimizationVariablesCell,...
                'latentVariables',obj.latentVariablesCell,...
                'P1constraints',obj.P1constraintsCell,...
                'P2constraints',obj.P2constraintsCell,...
                'latentConstraints',obj.latentConstraintsCell,...
                'outputExpressions',obj.outputExpressionsCell,...
                'parameters',obj.parametersCell,...
                'muFactorAggressive',obj.muAggressive,...
                'muFactorConservative',obj.muConservative,...
                'desiredDualityGap',obj.targetDualityGap,...
                'maxIter',obj.maxIter,...
                'gradTolerance',obj.gradientTolerance,...
                'skipAffine',true,...
                'delta',obj.delta,...
                'compilerOptimization',obj.compilerFlags,...
                'allowSave',obj.allowSave,...
                'profiling',false,...
                'solverVerboseLevel',obj.solverVerboseLevel,... 
                'verboseLevel',obj.verboseLevel); 
            opt = feval(classname);
        end

    end
end