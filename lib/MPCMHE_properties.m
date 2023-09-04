classdef MPCMHE_properties
    properties
        %objective
        objective

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
            obj.objective = objective;

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
        function obj = concat(obj, other)
            %{
                concatenate two constraints together.
                mostly used for concatenating attitude and translationsal formulations when needed
            %}
            obj.objective = obj.objective + other.objective;

            obj.P1optimizationVariablesCell = [obj.P1optimizationVariablesCell, other.P1optimizationVariablesCell];
            obj.P2optimizationVariablesCell = [obj.P2optimizationVariablesCell, other.P2optimizationVariablesCell];
            obj.latentVariablesCell = [obj.latentVariablesCell, other.latentVariablesCell];

            obj.P1constraintsCell = [obj.P1constraintsCell, other.P1constraintsCell];
            obj.P2constraintsCell = [obj.P2constraintsCell, other.P2constraintsCell];
            obj.latentConstraintsCell = [obj.latentConstraintsCell, other.latentConstraintsCell];

            obj.outputExpressionsCell = [obj.outputExpressionsCell, other.outputExpressionsCell];
            obj.parametersCell = [obj.parametersCell, other.parametersCell];
        end
        function opt = setupOptimization(obj)
            classname = cmex2equilibriumLatentCS(...
                'classname','tmpC_target_chaser_main',...
                'P1objective',obj.objective,...
                'P2objective',-obj.objective,...
                'P1optimizationVariables',obj.P1optimizationVariablesCell,...
                'P2optimizationVariables',obj.P2optimizationVariablesCell,...
                'latentVariables',obj.latentVariablesCell,...
                'P1constraints',obj.P1constraintsCell,...
                'P2constraints',obj.P2constraintsCell,...
                'latentConstraints',obj.latentConstraintsCell,...
                'outputExpressions',[{obj.objective}, obj.outputExpressionsCell],...
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