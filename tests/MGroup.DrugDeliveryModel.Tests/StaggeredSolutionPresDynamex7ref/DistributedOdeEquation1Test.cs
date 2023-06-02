/*using System;
using System.Collections.Generic;
using System.Linq;
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

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
	public class DistributedOdeEquation1Test
    {
        const double Sc = 0.1;
        const double timeStep = 1; // in days
        const double totalTime = 10; // in days

        static double miNormal = 5; //KPa
        static double kappaNormal = 6.667; //Kpa
        static double miTumor = 22.44; //Kpa
        static double kappaTumor = 201.74; //Kpa
        static int currentTimeStep = 0;
        static double lambda0 = 1;
        static Dictionary<double, double[]> Solution = new Dictionary<double, double[]>(); private static List<(INode node, IDofType dof)> watchDofs = new List<(INode node, IDofType dof)>();

		public DistributedOdeEquation1Test()
		{
            IsoparametricJacobian3D.DeterminantTolerance = 1e-20;
        }

        [Theory]
        [InlineData("../../../DataFiles/workingTetMesh155.mphtxt", 1, 1, 6.94e-6, 0.0083)]
        public void TestDistributedOdeApproach(string fileName, double coxCommmonValue, double TcommonValue, double k1, double k2)
        {
            double timeStep = 1;
            double totalTime = 10;
            ConvectionDiffusionDof lamdaMonitorDOF = ConvectionDiffusionDof.UnknownVariable;
            int nodeIdToMonitor = 0;

            Dictionary<int, double> T = new Dictionary<int, double>();// 500 [cells]
            Dictionary<int, double> cox = new Dictionary<int, double>();// 500 [cells]
            //Read geometry
            var comsolReader = new ComsolMeshReader(fileName);

            // initialize Shared quantities of Coupled model
            Dictionary<int, double> lambda = new Dictionary<int, double>(comsolReader.ElementConnectivity.Count());
            foreach (var elem in comsolReader.ElementConnectivity)
            {
                T.Add(elem.Key, TcommonValue);
                cox.Add(elem.Key, coxCommmonValue);
            }

            var modelBuilder = new DistributedOdeModelBuilder(comsolReader, k1, k2, cox, T, 1d,lamdaMonitorDOF,nodeIdToMonitor);

            var model = modelBuilder.GetModel();

            //Add the spatially distributed velocity field




            modelBuilder.AddBoundaryConditions(model);
            (var analyzer, var solver, var nlAnalyzers) =
                modelBuilder.GetAppropriateSolverAnalyzerAndLog(model, timeStep, totalTime, 0,1);
            modelBuilder.AddInitialConditionsForTheRestOfBulkNodes(model);
            ((NewmarkDynamicAnalyzer)analyzer).ResultStorage = new ImplicitIntegrationAnalyzerLog();

            analyzer.Initialize(true);
            analyzer.Solve();

            int totalNewmarkstepsNum = (int)Math.Truncate(totalTime / timeStep);
            var tCells = new double[totalNewmarkstepsNum];
            for (int i1 = 0; i1 < totalNewmarkstepsNum; i1++)
            {
                var timeStepResultsLog = ((NewmarkDynamicAnalyzer)analyzer).ResultStorage.Logs[i1];
                tCells[i1] = ((DOFSLog)timeStepResultsLog).DOFValues[model.GetNode(nodeIdToMonitor), lamdaMonitorDOF];
            }
        }

            [Theory]
        [InlineData("../../../DataFiles/workingTetMesh155.mphtxt", 1, 1, 6.94e-6, 0.0083)]
        public void SolveEquation1(string fileName, double cox, double T, double k1, double k2)
        {
            double capacity = 1;
            double dependentProductionCoeff = ((double)1 / 3) * k1 * cox * T / (k2 + cox);
            double convectionCoeff = 0;
            double independentSource = 0;
            double diffusion = 0;
            var modelProvider = new GenericComsol3DConvectionDiffusionProductionModelProvider(new double[] { convectionCoeff, convectionCoeff, convectionCoeff },
                diffusion, dependentProductionCoeff, independentSource, capacity);

            var model = modelProvider.CreateModelFromComsolFile(fileName);
            modelProvider.AddTopAndBottomBCs(model, 2, 1, 0, 1);
            modelProvider.AddInitialConditionsForTheRestOfBulkNodes(model, 2, 0, 1);


            var solverFactory = new DenseMatrixSolver.Factory() { IsMatrixPositiveDefinite = false}; //Dense Matrix Solver solves with zero matrices!
            var algebraicModel = solverFactory.BuildAlgebraicModel(model);
            var solver = solverFactory.BuildSolver(algebraicModel);
            var problem = new ProblemConvectionDiffusion(model, algebraicModel);

            var linearAnalyzer = new LinearAnalyzer(algebraicModel, solver, problem);

            var dynamicAnalyzerBuilder = new NewmarkDynamicAnalyzer.Builder( algebraicModel, problem, linearAnalyzer, timeStep: 1, totalTime: 10);
            var dynamicAnalyzer = dynamicAnalyzerBuilder.Build();

            var watchDofs = new List<(INode node, IDofType dof)>()
            {
                (model.NodesDictionary[36], ConvectionDiffusionDof.UnknownVariable),
            };

            linearAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs, algebraicModel);

            dynamicAnalyzer.ResultStorage = new ImplicitIntegrationAnalyzerLog();
            dynamicAnalyzer.Initialize();
            dynamicAnalyzer.Solve();

            DOFSLog log = (DOFSLog)linearAnalyzer.Logs[0];
            double computedValue = log.DOFValues[watchDofs[0].node, watchDofs[0].dof];

        }

        [Theory]
        [InlineData("../../../DataFiles/workingTetMesh155.mphtxt", 1, 6.94e-6, 0.0083,1,500)]
        public void SolveEquation2(string fileName, double cox, double k1, double k2, double vs, double T_initial)
        {
            double capacity = 1;
            double dependentProductionCoeff = k1 * cox / (k2 + cox);
            double convectionCoeff = vs;
            double independentSource = 0;
            double diffusion = 0;
            var modelProvider = new GenericComsol3DConvectionDiffusionProductionModelProvider(new double[] { convectionCoeff, convectionCoeff, convectionCoeff },
               diffusion, dependentProductionCoeff, independentSource, capacity);

            var model = modelProvider.CreateModelFromComsolFile(fileName);
            modelProvider.AddTopAndBottomBCs(model, 2, 1, 0, 1);
            modelProvider.AddInitialConditionsForTheRestOfBulkNodes(model, 2, 0, T_initial);


            var solverFactory = new DenseMatrixSolver.Factory() { IsMatrixPositiveDefinite = false }; 
            var algebraicModel = solverFactory.BuildAlgebraicModel(model);
            var solver = solverFactory.BuildSolver(algebraicModel);
            var problem = new ProblemConvectionDiffusion(model, algebraicModel);

            var linearAnalyzer = new LinearAnalyzer(algebraicModel, solver, problem);

            //var dynamicAnalyzerBuilder = new BDFDynamicAnalyzer.Builder( algebraicModel, problem, linearAnalyzer, timeStep: 1, totalTime: 10, bdfOrder: 5);
            //var dynamicAnalyzer = dynamicAnalyzerBuilder.Build();
            var dynamicAnalyzerBuilder = new NewmarkDynamicAnalyzer.Builder(algebraicModel, problem, linearAnalyzer, timeStep: 1, totalTime: 10);
            var dynamicAnalyzer = dynamicAnalyzerBuilder.Build();


            var watchDofs = new List<(INode node, IDofType dof)>()
            {
                (model.NodesDictionary[36], ConvectionDiffusionDof.UnknownVariable),
            };

            linearAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs, algebraicModel);

            dynamicAnalyzer.Initialize();
            dynamicAnalyzer.Solve();

        }

        

        /// <summary>
        /// 
        /// </summary>
        /// <param name="fileName"></param>
        /// <param name="kth"> gia to normal part 7.52×10-13 m2/kPasec  </param>
        /// <param name="Lp"> gia to normal part 2.7x10-9 m/kPa?sec </param>
        /// <param name="Sv">    7×103 m-1 for normal tissue</param>
        /// <param name="pv"> 4kPa </param>
        /// <param name="LplSvl">3.75 x10-1 [1/(kPa?sec)]</param>
        /// <param name="pl">0 Kpa</param>
        /// <param name="div_vs">scalar. In next development step it will be distributed in space from structural solution</param>
        [Theory]
        [InlineData("../../../DataFiles/workingTetMesh155.mphtxt", 7.52e-13, *//*2.7e-9*//*1, 7e3, *//*4*//*1,
                                                                  *//*3.75e-1*//* 1, 0, 0)]
        public void SolveEquation7and8ofUpdatedReportStaticZero(string fileName, double kth, double Lp, double Sv, double pv,
                                                                double LplSvl, double pl, double div_vs)
        {
            double convectionCoeff = 1;

            double capacity = 0;
            double dependentProductionCoeff = -(Lp * Sv + LplSvl);
            double independentSource = Lp * Sv * pv + LplSvl * pl - div_vs;
            double diffusion = kth;

            var modelProvider = new GenericComsol3DConvectionDiffusionProductionModelProvider(new double[] { convectionCoeff, convectionCoeff, convectionCoeff },
                diffusion, dependentProductionCoeff, independentSource, capacity);

            var model = modelProvider.CreateModelFromComsolFile(fileName);
            modelProvider.AddTopAndBottomBCs(model, 2, 100, 0, 50);
            var solverFactory = new DenseMatrixSolver.Factory() { IsMatrixPositiveDefinite = false}; //Dense Matrix Solver solves with zero matrices!
            var algebraicModel = solverFactory.BuildAlgebraicModel(model);
            var solver = solverFactory.BuildSolver(algebraicModel);
            var problem = new ProblemConvectionDiffusion(model, algebraicModel);

            var linearAnalyzer = new LinearAnalyzer(algebraicModel, solver, problem);

            //var dynamicAnalyzerBuilder = new StaticAnalyzer.Builder(model, algebraicModel, problem, linearAnalyzer, timeStep: 1, totalTime: 10, bdfOrder: 5);
            //var cAnalyzer = dynamicAnalyzerBuilder.Build();
            var staticAnalyzer = new StaticAnalyzer( algebraicModel, problem, linearAnalyzer);


            //dynamicAnalyzer.ResultStorage = new ImplicitIntegrationAnalyzerLog();

            //var watchDofs = new List<(INode node, IDofType dof)>()
            //{
            //    (model.NodesDictionary[36], ConvectionDiffusionDof.UnknownVariable),
            //};

            //linearAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs, algebraicModel);

            staticAnalyzer.Initialize();
            staticAnalyzer.Solve();

            DOFSLog log = (DOFSLog)linearAnalyzer.Logs[0];
            double computedValue = log.DOFValues[watchDofs[0].node, watchDofs[0].dof];

        }

        
        [InlineData("../../../DataFiles/3d8Hexa.mphtxt", 1, *//*2.7e-9*//*1, -1, *//*4*//*-1,
                                                                  *//*3.75e-1*//* 0, 0, 0)]
        public void SolveMSolveTestofEquation7and8(string fileName, double kth, double Lp, double Sv, double pv,
                                                                double LplSvl, double pl, double div_vs)
        {
            double convectionCoeff = 1;

            double capacity = 0;
            double dependentProductionCoeff = -(Lp * Sv + LplSvl);
            double independentSource = Lp * Sv * pv + LplSvl * pl - div_vs;
            double diffusion = kth;

            var modelProvider = new GenericComsol3DConvectionDiffusionProductionModelProvider(new double[] { convectionCoeff, convectionCoeff, convectionCoeff },
                diffusion, dependentProductionCoeff, independentSource, capacity);

            var model = modelProvider.CreateModelFromComsolFile(fileName);
            modelProvider.AddTopAndBottomBCs(model, 2, 100, 0, 50);
            var solverFactory = new DenseMatrixSolver.Factory() { IsMatrixPositiveDefinite = false }; //Dense Matrix Solver solves with zero matrices!
            var algebraicModel = solverFactory.BuildAlgebraicModel(model);
            var solver = solverFactory.BuildSolver(algebraicModel);
            var problem = new ProblemConvectionDiffusion(model, algebraicModel);

            var linearAnalyzer = new LinearAnalyzer(algebraicModel, solver, problem);

            //var dynamicAnalyzerBuilder = new StaticAnalyzer.Builder(model, algebraicModel, problem, linearAnalyzer, timeStep: 1, totalTime: 10, bdfOrder: 5);
            //var cAnalyzer = dynamicAnalyzerBuilder.Build();
            var staticAnalyzer = new StaticAnalyzer(algebraicModel, problem, linearAnalyzer);


            //dynamicAnalyzer.ResultStorage = new ImplicitIntegrationAnalyzerLog();

            //var watchDofs = new List<(INode node, IDofType dof)>()
            //{
            //    (model.NodesDictionary[36], ConvectionDiffusionDof.UnknownVariable),
            //};

            //linearAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs, algebraicModel);

            staticAnalyzer.Initialize();
            staticAnalyzer.Solve();

            DOFSLog log = (DOFSLog)linearAnalyzer.Logs[0];
            double computedValue = log.DOFValues[watchDofs[0].node, watchDofs[0].dof];

        }

    }
}
*/