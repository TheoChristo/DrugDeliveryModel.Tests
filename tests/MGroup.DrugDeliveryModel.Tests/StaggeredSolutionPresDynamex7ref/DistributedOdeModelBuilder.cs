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
using MGroup.FEM.ConvectionDiffusion.Tests.Commons;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.Solution;
using MGroup.Constitutive.ConvectionDiffusion.BoundaryConditions;
using MGroup.Constitutive.ConvectionDiffusion.InitialConditions;
using System.Security.AccessControl;
using MGroup.LinearAlgebra.Matrices;
using MGroup.Solvers.AlgebraicModel;
using MGroup.FEM.ConvectionDiffusion.Isoparametric;
using System.Reflection.PortableExecutable;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Numerics.Integration.Quadratures;
using MGroup.MSolve.Numerics.Interpolation;
using System.Xml.Linq;

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
    public class DistributedOdeModelBuilder
    {
        private double K1 { get; }
        private double K2 { get; }

        private double InitialCondition { get; }

        private Dictionary<int, double> DomainCOx { get; }

        private readonly Dictionary<int, double> domainT; // [cells]

        private ComsolMeshReader modelReader;
        public GlobalAlgebraicModel<Matrix> algebraicModel;

        public Dictionary<int, INode> DummyNodesDictionary { get; } = new Dictionary<int, INode>();
        public Dictionary<int, Tuple<CellType, Node[], int>> ElementDummyNodesConnectivity { get; set; }

        private ConvectionDiffusionDof MonitorDOFType { get; }
        private int MonitorNodeId { get; }

        public DistributedOdeModelBuilder(ComsolMeshReader modelReader, double k1, double k2, Dictionary<int, double> domainCOx, Dictionary<int, double> T, double initialCondition,
            ConvectionDiffusionDof monitorDOFType, int monitorNodeId)
        {
            this.MonitorDOFType = monitorDOFType;
            this.MonitorNodeId = monitorNodeId;

            K1 = k1;
            K2 = k2;
            DomainCOx = domainCOx;
            this.domainT = T;
            this.InitialCondition = initialCondition;

            this.modelReader = modelReader;
            IsoparametricJacobian3D.DeterminantTolerance = 1e-30;


            ElementDummyNodesConnectivity = new Dictionary<int, Tuple<CellType, Node[], int>>();

            int dummyNodeId = 0;    
            foreach (var elementConnectivity in modelReader.ElementConnectivity)
            {
                if(elementConnectivity.Value.Item1== CellType.Tet4)
                {
                    var interpolation = InterpolationTet4.UniqueInstance;
                    var quadrature = TetrahedronQuadrature.Order1Point1;

                    var numOfNecessaryDymmies = quadrature.IntegrationPoints.Count();
                    var elemDummyNodes = new Node[numOfNecessaryDymmies];
                    for (int i1 = 0; i1 < numOfNecessaryDymmies; i1++)
                    {
                        var node = new Node(id: dummyNodeId, x: 0, y: 0, z: 0);
                        elemDummyNodes[i1] = node;
                        DummyNodesDictionary.Add(key: dummyNodeId, node);
                       dummyNodeId++;
                    }

                    ElementDummyNodesConnectivity.Add(key: elementConnectivity.Key, value: new Tuple<CellType, Node[], int>(CellType.Tet4, elemDummyNodes, elementConnectivity.Value.Item3));
                }
                else
                {
                    throw new NotImplementedException();
                }

            }

        }

        public Model GetModel()
        {
            var model = new Model();
            model.SubdomainsDictionary[0] = new Subdomain(id: 0);


            foreach (var node in DummyNodesDictionary.Values)
            {
                model.NodesDictionary.Add(node.ID, node);
            }



            foreach (var elementConnectivity in ElementDummyNodesConnectivity)
            {
                if (elementConnectivity.Value.Item1 == CellType.Tet4)
                {
                    var quadrature = TetrahedronQuadrature.Order1Point1;
                    var numOfNecessaryDymmies = quadrature.IntegrationPoints.Count();

                    var coeffs = new double[numOfNecessaryDymmies];
                    var cox = DomainCOx[elementConnectivity.Key];
                    var t = domainT[elementConnectivity.Key];
                    coeffs[0] = ((double)1 / 3) * K1 * cox * t / (K2 + cox);
                    var element1 = new DistributedFirstOrderODEElement3D(elementConnectivity.Value.Item2, quadrature, coeffs);
                    element1.ID = elementConnectivity.Key; ;
                    model.ElementsDictionary.Add(elementConnectivity.Key, element1);
                    model.SubdomainsDictionary[0].Elements.Add(element1);
                }
                else
                {
                    throw new NotImplementedException();
                }

            }


            return model;
        }


        public void AddInitialConditionsForTheRestOfBulkNodes(Model model)
        {
            var intitalConditions = new List<INodalConvectionDiffusionInitialCondition>();
            foreach (var node in DummyNodesDictionary.Values)
            {
                
                    intitalConditions.Add(new NodalInitialUnknownVariable(node, ConvectionDiffusionDof.UnknownVariable, InitialCondition));
               

            }
            
            model.InitialConditions.Add(new ConvectionDiffusionInitialConditionSet(intitalConditions,
                new DomainInitialUnknownVariable[]
                { }));


        }

        public void AddBoundaryConditions(Model model)
        {
            //BoundaryAndInitialConditionsUtility.AssignConvectionDiffusionDirichletBCsToModel(model, convectionDiffusionDirichletBC, 1e-3);
        }

        //public (IParentAnalyzer analyzer, ISolver solver, IChildAnalyzer loadcontrolAnalyzer) GetAppropriateSolverAnalyzerAndLog
        //(Model model, double pseudoTimeStep, double pseudoTotalTime, int currentStep, int nIncrements)
        //{
        //    var solverFactory = new DenseMatrixSolver.Factory() { IsMatrixPositiveDefinite = false }; //Dense Matrix Solver solves with zero matrices!
        //    //var solverFactory = new SkylineSolver.Factory() { FactorizationPivotTolerance = 1e-8 };
        //    algebraicModel = solverFactory.BuildAlgebraicModel(model);
        //    var solver = solverFactory.BuildSolver(algebraicModel);
        //    var provider = new ProblemConvectionDiffusion(model, algebraicModel);

        //    var loadControlAnalyzerBuilder = new LoadControlAnalyzer.Builder(algebraicModel, solver, provider, numIncrements: 1)
        //    {
        //        ResidualTolerance = 1E-8,
        //        MaxIterationsPerIncrement = 100,
        //        NumIterationsForMatrixRebuild = 1
        //    };
        //    var loadControlAnalyzer = loadControlAnalyzerBuilder.Build();

        //    loadControlAnalyzer.TotalDisplacementsPerIterationLog = new TotalDisplacementsPerIterationLog(new List<(INode node, IDofType dof)>()
        //    {(model.NodesDictionary[MonitorNodeId], MonitorDOFType)}, algebraicModel);

        //    var analyzerBuilder = new NewmarkDynamicAnalyzer.Builder(algebraicModel, provider, loadControlAnalyzer, timeStep: pseudoTimeStep, totalTime: pseudoTotalTime, false, currentStep: currentStep);
        //    analyzerBuilder.SetNewmarkParametersForConstantAcceleration();
        //    /*var analyzerBuilder = new BDFDynamicAnalyzer.Builder(algebraicModel, provider, loadControlAnalyzer,
        //        timeStep: pseudoTimeStep, totalTime: pseudoTotalTime, 5);*/
        //    var analyzer = analyzerBuilder.Build();

        //    var watchDofs = new[]
        //    {
        //        new List<(INode node, IDofType dof)>()
        //        {
        //            (model.NodesDictionary[MonitorNodeId], MonitorDOFType),
        //        }
        //    };
        //    loadControlAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs[0], algebraicModel);

        //    return (analyzer, solver, loadControlAnalyzer);
        //}
/*
        public (IParentAnalyzer analyzer, ISolver solver, IChildAnalyzer loadcontrolAnalyzer) GetAppropriateSolverAnalyzerAndLog(Model model, double TimeStep, double TotalTime, int currentStep, int nIncrements)
        {
            var solverFactory = new DenseMatrixSolver.Factory() { IsMatrixPositiveDefinite = false };
            algebraicModel = solverFactory.BuildAlgebraicModel(model);
            var solver = solverFactory.BuildSolver(algebraicModel);
            var problem = new ProblemConvectionDiffusion(model, algebraicModel);

            var linearAnalyzer = new LinearAnalyzer(algebraicModel, solver, problem);


            var analyzerBuilder = new NewmarkDynamicAnalyzer.Builder(algebraicModel, problem, linearAnalyzer, timeStep: TimeStep, totalTime: TotalTime, false, currentStep: currentStep);
            analyzerBuilder.SetNewmarkParametersForConstantAcceleration();

            var analyzer = analyzerBuilder.Build();
            var watchDofs = new[]
            {
                new List<(INode node, IDofType dof)>()
                {
                    (model.NodesDictionary[MonitorNodeId], MonitorDOFType),
                }
            };
            linearAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs[0], algebraicModel);

            return (analyzer, solver, linearAnalyzer);
        }

        public void RetrieveLambdaSolution(ISolver solver, IChildAnalyzer childAnalyzer, Model model, Dictionary<int, double> lambda)
        {
            //List<double> p = new List<double>();

            var u = childAnalyzer.CurrentAnalysisResult;

            var lamdaSolution = algebraicModel.ExtractAllResults(u);
            var lamdaNodals = lamdaSolution.Data;

            foreach (var elementConnectivity in modelReader.ElementConnectivity)
            {
                var dummyNodes = ElementDummyNodesConnectivity[elementConnectivity.Key].Item2;

                if (elementConnectivity.Value.Item1 == CellType.Tet4)
                {
                    var nodalLambda = lamdaNodals[dummyNodes[0].ID, 0];
                    lambda[elementConnectivity.Key] = nodalLambda;
                }
                else
                {
                    throw new NotImplementedException();
                }
            }


            
        }
    }
}
*/