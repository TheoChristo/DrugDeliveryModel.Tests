/*using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Numerics.Interpolation.Jacobians;
using MGroup.DrugDeliveryModel.Tests.EquationModels;
using MGroup.NumericalAnalyzers.Dynamic;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.NumericalAnalyzers.Staggered;
using MGroup.NumericalAnalyzers.Discretization.NonLinear;
using MGroup.Constitutive.Structural;
using MGroup.DrugDeliveryModel.Tests.Commons;
using MGroup.NumericalAnalyzers;
using MGroup.Solvers.Direct;
using Xunit;
using MGroup.Constitutive.ConvectionDiffusion;
using MGroup.Constitutive.ConvectionDiffusion.BoundaryConditions;
using MGroup.Constitutive.ConvectionDiffusion.InitialConditions;

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
    public class EquationsTests13DistributedModelBuilderCox
    {
        //---------------------------------------Equation Cox---------------------------------------

        //---------------------------------------Variables------------------------------------------
        /// <summary>
        /// The average value of the three components of the fluid velocity vector  [m/s]
        /// </summary>
        private const double FluidSpeed = 2.32E-4; // [m/s]

        /// <summary>
        /// Diffusivity of oxygen [m2/s]
        /// </summary>
        private const double Dox = 1.79E-5; // [m2/s]

        /// <summary>
        /// Oxygen uptake [mol/(m3*s)]
        /// </summary>
        private const double Aox = 2.55E-2; // [mol/(m3*s)]

        /// <summary>
        /// Oxygen uptake [mol/m3]
        /// </summary>
        private const double Kox = 4.64E-3; // [mol / m3]

        /// <summary>
        /// Oxygen permeability across tumor vessel walls [m/s]
        /// </summary>
        private const double PerOx = 3.55E-4; // [m/s]

        /// <summary>
        /// Vascular Density [1/m]
        /// </summary>
        private const double Sv = 7E3; // [1/m]

        /// <summary>
        /// Initial Oxygen Concentration [mol/m3]
        /// </summary>
        private const double CInitOx = 0.2; // [mol/m3]

        /// <summary>
        /// Cancer cell density [1]
        /// </summary>
        private const double T = 500d; // [cells]

        //---------------------------------------Logging----------------------------------
        /// <summary>
        /// The degree of freedom that will be monitored for equation cox
        /// </summary>
        private ConvectionDiffusionDof coxMonitorDOF = ConvectionDiffusionDof.UnknownVariable;

        /// <summary>
        /// The coordinates of the monitored node
        /// </summary>
        private double[] monitorNodeCoords = { 0.05, 0.05, 0.05 };


        //---------------------------------------Time Discretization Specs------------------------------
        private const double TotalTime = 2.5E-3; //0.0025

        /// <summary>
        /// For increased accuracy use time-step of order 1E-5
        /// </summary>
        private const double TimeStep = 1E-5;// sec

        /// <summary>
        /// Simplified version of the production term without non-linear term
        /// </summary>
        private readonly Func<double> simpleDependentLinearSource =() => -PerOx * Sv;

        private readonly Func<double> independentLinearSource =() => PerOx * Sv * 0.2d;

        public void EquationsTests13DistributedModelBuilder()
        {
            IsoparametricJacobian3D.DeterminantTolerance = 1e-20;
        }

        //Non - Linear Change here for linear problem
        public static double ProductionFuncWithoutConstantTerm(double Cox)
        {
            return -PerOx * Sv * Cox - Aox * T * Cox / (Cox + Kox);
        }

        //Non - Linear Derivative
        public static double ProductionFuncWithoutConstantTermDDerivative(double Cox)
        {
            return -PerOx * Sv - Aox * T / (Cox + Kox) + Aox * T * Cox * Math.Pow(Cox + Kox, -2);
        }

        public static double LinearProductionFuncWithoutConstantTerm(double Cox)
        {
            return -PerOx*Sv*Cox;
        }

        public static double LinearProductionFuncWithoutConstantTermDDerivative(double Cox)
        {
            return -PerOx * Sv;
        }

        [Theory]
        [InlineData("../../../DataFiles/workingTetMesh2185_1Domain.mphtxt")]
        public void SolveEquationCOxLinearProduction(string fileName)
        {
            return;
            var capacity = 1;
            var diffusionCoefficient = Dox;
            var convectionCoefficient = FluidSpeed;
            var dependentProductionCoefficient = simpleDependentLinearSource();
            var independentSourceCoefficient = independentLinearSource();

            //Read Mesh From comsol file
            var mesh = new ComsolMeshReader(fileName);

            //Assign equation properties to the domain elements
            var convectionDomainCoefficients = new Dictionary<int, double[]>();
            var dependentProductionCoefficients = new Dictionary<int, double>();
            var independentProductionCoefficients = new Dictionary<int, double>();

            foreach (var elementConnectivity in mesh.ElementConnectivity)
            {
                convectionDomainCoefficients[elementConnectivity.Key] = new double[] { convectionCoefficient, convectionCoefficient, convectionCoefficient };
                dependentProductionCoefficients[elementConnectivity.Key] = dependentProductionCoefficient;
                independentProductionCoefficients[elementConnectivity.Key] = independentSourceCoefficient;
            }

            //Create Model
            var modelProvider = new GenericComsol3DConvectionDiffusionProductionModelProviderDistributedSpace(mesh);
            var model = modelProvider.CreateModelFromComsolFile(convectionDomainCoefficients, diffusionCoefficient,
                dependentProductionCoefficients, independentProductionCoefficients, capacity);

            //Assign Boundary Conditions
            AddTopRightBackNodesBC(model, 0d, 0.1, 0, 0.1, 0, 0.1, CInitOx);
            //AddInitialConditions(model, 0d, 0.1, 0, 0.1, 0, 0.1, CInitOx);


            var solverFactory = new DenseMatrixSolver.Factory() { IsMatrixPositiveDefinite = false };
            var algebraicModel = solverFactory.BuildAlgebraicModel(model);
            var solver = solverFactory.BuildSolver(algebraicModel);
            var problem = new ProblemConvectionDiffusion(model, algebraicModel);

            var linearAnalyzer = new LinearAnalyzer(algebraicModel, solver, problem);

            // var dynamicAnalyzerBuilder = new NewmarkDynamicAnalyzer.Builder(algebraicModel, problem, linearAnalyzer, timeStep: TimeStep, totalTime: TotalTime);
            var dynamicAnalyzerBuilder = new BDFDynamicAnalyzer.Builder(algebraicModel, problem, linearAnalyzer, timeStep: TimeStep, totalTime: TotalTime, bdfOrder: 5);
            var dynamicAnalyzer = dynamicAnalyzerBuilder.Build();

            // Create a log for the desired dof
            var nodeIdToMonitor =Utilities.FindNodeIdFromNodalCoordinates(mesh.NodesDictionary, monitorNodeCoords, 1e-3);
            var watchDofs = new List<(INode node, IDofType dof)>()
            {
                (model.NodesDictionary[nodeIdToMonitor], coxMonitorDOF),
            };

            linearAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs, algebraicModel);


            dynamicAnalyzer.ResultStorage = new ImplicitIntegrationAnalyzerLog();
            dynamicAnalyzer.Initialize();
            Console.WriteLine("Solving Cox Linear prod");
            dynamicAnalyzer.Solve();

            int totalNewmarkstepsNum = (int)Math.Truncate(TotalTime / TimeStep);
            var cox = new double[totalNewmarkstepsNum];
            for (int i1 = 0; i1 < totalNewmarkstepsNum; i1++)
            {
                var timeStepResultsLog = dynamicAnalyzer.ResultStorage.Logs[i1];
                cox[i1] = ((DOFSLog)timeStepResultsLog).DOFValues[model.GetNode(nodeIdToMonitor), coxMonitorDOF];
            }

            CSVExporter.ExportVectorToCSV(cox, "../../../Integration/cox_linear_nodes_mslv.csv");
            Console.WriteLine("FINISHED solving Cox Linear prod");
        }

        [Theory]
        [InlineData("../../../DataFiles/workingTetMesh2185_1Domain.mphtxt")]
        public void SolveEquationCOxNonLinearProduction(string fileName)
        {
            var capacity = 1;
            var diffusionCoefficient = Dox;
            var convectionCoefficient = FluidSpeed;
            var independentSourceCoefficient = independentLinearSource();
            var dependentSourceCoefficient = 0;//simpleDependentLinearSource();

            //Read Mesh From comsol file
            var mesh = new ComsolMeshReader(fileName);

            //Assign equation properties to the domain elements
            var convectionDomainCoefficients = new Dictionary<int, double[]>();
            var dependentProductionCoefficients = new Dictionary<int, double>();
            var independentProductionCoefficients = new Dictionary<int, double>();

            foreach (var elementConnectivity in mesh.ElementConnectivity)
            {
                convectionDomainCoefficients[elementConnectivity.Key] = new double[] { convectionCoefficient, convectionCoefficient, convectionCoefficient };
                dependentProductionCoefficients[elementConnectivity.Key] = dependentSourceCoefficient;
                independentProductionCoefficients[elementConnectivity.Key] = independentSourceCoefficient;
            }

            //Create Model
            var modelProvider = new GenericComsol3DConvectionDiffusionProductionModelProviderDistributedSpace(mesh);
            var model = modelProvider.CreateModelFromComsolFile(convectionDomainCoefficients, diffusionCoefficient,
                dependentProductionCoefficients, independentProductionCoefficients, capacity,
                LinearProductionFuncWithoutConstantTerm, LinearProductionFuncWithoutConstantTermDDerivative);
            //var model = modelProvider.CreateModelFromComsolFile(convectionDomainCoefficients, diffusionCoefficient,
            //    dependentProductionCoefficients, independentProductionCoefficients, capacity);

            //Assign Boundary Conditions
            AddTopRightBackNodesBC(model, 0d, 0.1, 0, 0.1, 0, 0.1, CInitOx);
            AddInitialConditions(model, 0d, 0.1, 0, 0.1, 0, 0.1, CInitOx);


            var solverFactory = new DenseMatrixSolver.Factory() { IsMatrixPositiveDefinite = false };
            var algebraicModel = solverFactory.BuildAlgebraicModel(model);
            var solver = solverFactory.BuildSolver(algebraicModel);
            var problem = new ProblemConvectionDiffusion(model, algebraicModel);

            var loadControlAnalyzerBuilder = new LoadControlAnalyzer.Builder(algebraicModel, solver, problem, numIncrements: 2)
            {
                ResidualTolerance = 1E-3,
                MaxIterationsPerIncrement = 100,
                NumIterationsForMatrixRebuild = 1
            };

            var linearAnalyzer = loadControlAnalyzerBuilder.Build();//new LinearAnalyzer(algebraicModel, solver, problem);

            var dynamicAnalyzerBuilder = new NewmarkDynamicAnalyzer.Builder(algebraicModel, problem, linearAnalyzer, timeStep: TimeStep, totalTime: TotalTime);
            //var dynamicAnalyzerBuilder = new BDFDynamicAnalyzer.Builder(algebraicModel, problem, linearAnalyzer, timeStep: TimeStep, totalTime: TotalTime, bdfOrder: 5);
            var dynamicAnalyzer = dynamicAnalyzerBuilder.Build();

            // Create a log for the desired dof
            var nodeIdToMonitor =Utilities.FindNodeIdFromNodalCoordinates(mesh.NodesDictionary, monitorNodeCoords, 1e-2);
            var watchDofs = new List<(INode node, IDofType dof)>()
            {
                (model.NodesDictionary[nodeIdToMonitor], coxMonitorDOF),
            };

            linearAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs, algebraicModel);


            dynamicAnalyzer.ResultStorage = new ImplicitIntegrationAnalyzerLog();
            dynamicAnalyzer.Initialize();
            Console.WriteLine("Solving Cox Non-Linear prod");
            dynamicAnalyzer.Solve();

            int totalNewmarkstepsNum = (int)Math.Truncate(TotalTime / TimeStep);
            var cox = new double[totalNewmarkstepsNum];
            for (int i1 = 0; i1 < totalNewmarkstepsNum; i1++)
            {
                var timeStepResultsLog = dynamicAnalyzer.ResultStorage.Logs[i1];
                cox[i1] = ((DOFSLog)timeStepResultsLog).DOFValues[model.GetNode(nodeIdToMonitor), coxMonitorDOF];
            }

            CSVExporter.ExportVectorToCSV(cox, "../../../Integration/cox_non_linear_nodes_mslv.csv");
            Console.WriteLine("FINISHED solving Cox Non-Linear prod");

        }


        private void AddTopRightBackNodesBC(Model model,double modelMinX,double modelMaxX,
                                                        double modelMinY,double modelMaxY,
                                                        double modelMinZ,double modelMaxZ, double boundaryCondition)
        {
            var topNodes = new List<INode>();
            var rightNodes = new List<INode>();
            var backNodes = new List<INode>();

            var tol = 1E-5;

            foreach (var node in model.NodesDictionary.Values)
            {
                if  (Math.Abs(modelMaxZ - node.Z) < tol) topNodes.Add(node);
                if (Math.Abs(modelMaxX - node.X) < tol) rightNodes.Add(node);
                if (Math.Abs(modelMaxY - node.Y) < tol) backNodes.Add(node);
            }
            //Union of all boundary nodes in a single enumerable.
            var peripheralNodes = topNodes.Union(backNodes).Union(rightNodes);

            var dirichletBCs = new List<NodalUnknownVariable>();

            //Add the prescribed value to all boundary nodes
            foreach (var node in peripheralNodes)
            {
                dirichletBCs.Add(new NodalUnknownVariable(node, ConvectionDiffusionDof.UnknownVariable, boundaryCondition));
            }

            model.BoundaryConditions.Add(new ConvectionDiffusionBoundaryConditionSet(dirichletBCs, new INodalConvectionDiffusionNeumannBoundaryCondition[] {}));
        }

        private void AddInitialConditions(Model model,
                                                  double modelMinX,double modelMaxX,
                                                  double modelMinY,double modelMaxY,
                                                  double modelMinZ,double modelMaxZ, double initialCondition)
        {
            var tol = 1E-5;
            var initialConditions = new List<INodalConvectionDiffusionInitialCondition>();
            foreach (var node in model.NodesDictionary.Values)
            {
                if  ((Math.Abs(modelMaxZ - node.Z) >= tol) && (Math.Abs(modelMaxX - node.X) >= tol) && (Math.Abs(modelMaxY - node.Y) >= tol))
                    initialConditions.Add(new NodalInitialUnknownVariable(node, ConvectionDiffusionDof.UnknownVariable, initialCondition));
            }
            model.InitialConditions.Add(new ConvectionDiffusionInitialConditionSet(initialConditions, new DomainInitialUnknownVariable[]{ }));
        }


    }
}
*/