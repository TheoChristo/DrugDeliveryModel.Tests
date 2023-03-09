using System;
using System.Collections.Generic;
using MGroup.Constitutive.ConvectionDiffusion;
using MGroup.DrugDeliveryModel.Tests.Commons;
using MGroup.DrugDeliveryModel.Tests.EquationModels;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Numerics.Interpolation.Jacobians;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.NumericalAnalyzers.Dynamic;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.NumericalAnalyzers.Discretization.NonLinear;
using MGroup.Solvers.Direct;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.Solution;
using MGroup.NumericalAnalyzers;
using MGroup.Constitutive.ConvectionDiffusion.BoundaryConditions;
using System.Linq;
using System.Xml.Linq;
using BC = MGroup.DrugDeliveryModel.Tests.Commons.BoundaryAndInitialConditionsUtility.BoundaryConditionCase;

namespace MGroup.DrugDeliveryModel.Tests.PreliminaryModels;

public class TCellModelProvider
{
    private double K1 { get; }
    private double K2 { get; }

    private Dictionary<int, double> DomainCOx { get; }

    private Dictionary<int, double[][]> SolidVelocity { get; }
    private ComsolMeshReader Mesh { get; }
    private ConvectionDiffusionDof MonitorDOFType { get; }
    private int MonitorNodeId { get; }
    private List<(BC, ConvectionDiffusionDof[], double[][], double[])> DirichletBCs { get; }
    private List<(BC, ConvectionDiffusionDof[], double[][], double[])> NeumannBCs { get; }

    private double InitialCondition { get; }

    public TCellModelProvider(double k1, double k2, Dictionary<int, double> domainCOx, Dictionary<int, double[][]> solidVelocity,
        ComsolMeshReader mesh,
        ConvectionDiffusionDof tCellMonitorDOFType, int monitorNodeId,
        List<(BC, ConvectionDiffusionDof[], double[][], double[])> dirichletBCs,
        List<(BC, ConvectionDiffusionDof[], double[][], double[])> neumannBCs,
        double initialCondition)
    {
        K1 = k1;
        K2 = k2;
        DomainCOx = domainCOx;
        SolidVelocity = solidVelocity;
        Mesh = mesh;
        MonitorDOFType = tCellMonitorDOFType;
        MonitorNodeId = monitorNodeId;
        DirichletBCs = dirichletBCs;
        NeumannBCs = neumannBCs;
        InitialCondition = initialCondition;

        IsoparametricJacobian3D.DeterminantTolerance = 1E-20;
    }

    public Model GetModel()
    {
        var capacity = 1;
        var diffusionCoefficient = 0d;


        //Assign equation properties to the domain elements
        var convectionDomainCoefficients = new Dictionary<int, double[]>();
        var dependentProductionCoefficients = new Dictionary<int, double>();
        var independentProductionCoefficients = new Dictionary<int, double>();

        foreach (var elementConnectivity in Mesh.ElementConnectivity)
        {
            var vs = SolidVelocity[elementConnectivity.Key];
            convectionDomainCoefficients[elementConnectivity.Key] = new double[] { vs[0][0], vs[0][1], vs[0][2]};
            
            var elementCOx = DomainCOx[elementConnectivity.Key];
            var dependentProductionCoefficient = (K1 * elementCOx) / (K2 + elementCOx);
            dependentProductionCoefficients[elementConnectivity.Key] = dependentProductionCoefficient;

            independentProductionCoefficients[elementConnectivity.Key] = 0d;
        }

        //Create Model
        var modelProvider = new GenericComsol3DConvectionDiffusionProductionModelProviderDistributedSpace(Mesh);
        var model = modelProvider.CreateModelFromComsolFile(convectionDomainCoefficients, diffusionCoefficient,
            dependentProductionCoefficients, independentProductionCoefficients, capacity);

        return model;
    }

    //public void AddBoundaryConditions2(Model model)
    //{
    //    BoundaryAndInitialConditionsUtility.AssignConvectionDiffusionDirichletBCsToModel(model, DirichletBCs, 1E-5);

    //}

    public void AddBoundaryConditions(Model model)
    {
        foreach (var BCdata in DirichletBCs)
        {
            var RegionType = BCdata.Item1;
            var Bcstype = BCdata.Item2;
            double[][] CaracteristicCoords = BCdata.Item3;
            double[] prescrVal = BCdata.Item4;

            switch (RegionType)
            {
                case BC.TopRightBackDiriclet:
                    {
                        double modelMinX = CaracteristicCoords[0][0]; double modelMinY = CaracteristicCoords[0][1]; double modelMinZ = CaracteristicCoords[0][2];
                        double modelMaxX = CaracteristicCoords[1][0]; double modelMaxY = CaracteristicCoords[1][1]; double modelMaxZ = CaracteristicCoords[1][2];
                        AddTopRightBackNodesBC(model, prescrVal[0], modelMinX, modelMaxX, modelMinY, modelMaxY, modelMinZ, modelMaxZ);
                        break;
                    }
                case BC.LeftDirichlet or BC.RightDirichlet:
                    {
                        double xCoordin = CaracteristicCoords[0][0];
                        AddXFaceBcs(model, xCoordin, prescrVal[0], ConvectionDiffusionDof.UnknownVariable);
                        break;
                    }
                case BC.FrontDirichlet or BC.BackDirichlet:
                    {
                        double yCoordin = CaracteristicCoords[0][1];
                        AddYFaceBcs(model, yCoordin, prescrVal[0], ConvectionDiffusionDof.UnknownVariable);
                        break;
                    }
                case BC.TopDirichlet or BC.BottomDirichlet:
                    {
                        // TODO: pithanws to ConvectionDiffusionDof.UnknownVariable antikathistatai me duo periptwseis analoga to Bcstype
                        double zCoordin = CaracteristicCoords[0][2];
                        AddZFaceBcs(model, zCoordin, prescrVal[0], ConvectionDiffusionDof.UnknownVariable);
                        break;
                    }

            }
        }

    }

    public void AddXFaceBcs(Model model, double xCoordOfFace, double dirichleValue, ConvectionDiffusionDof dofTypeToPrescribeDircihle)
    {

        var xFaceNodes = new List<INode>();
        foreach (var node in model.NodesDictionary.Values)
        {
            if (Math.Abs(xCoordOfFace - node.X) < 1E-9) xFaceNodes.Add(node);
        }

        var dirichletBCs = new List<NodalUnknownVariable>();
        foreach (var node in xFaceNodes)
        {
            dirichletBCs.Add(new NodalUnknownVariable(node, dofTypeToPrescribeDircihle, dirichleValue));
        }


        model.BoundaryConditions.Add(new ConvectionDiffusionBoundaryConditionSet(
            dirichletBCs,
            new INodalConvectionDiffusionNeumannBoundaryCondition[] { }
        ));

    }

    public void AddYFaceBcs(Model model, double yCoordOfFace, double dirichleValue, ConvectionDiffusionDof dofTypeToPrescribeDircihle)
    {

        var yFaceNodes = new List<INode>();
        foreach (var node in model.NodesDictionary.Values)
        {
            if (Math.Abs(yCoordOfFace - node.Y) < 1E-9) yFaceNodes.Add(node);
        }

        var dirichletBCs = new List<NodalUnknownVariable>();
        foreach (var node in yFaceNodes)
        {
            dirichletBCs.Add(new NodalUnknownVariable(node, dofTypeToPrescribeDircihle, dirichleValue));
        }


        model.BoundaryConditions.Add(new ConvectionDiffusionBoundaryConditionSet(
            dirichletBCs,
            new INodalConvectionDiffusionNeumannBoundaryCondition[] { }
        ));

    }


    public void AddZFaceBcs(Model model, double zCoordOfFace, double dirichleValue, ConvectionDiffusionDof dofTypeToPrescribeDircihle)
    {

        var zFaceNodes = new List<INode>();
        foreach (var node in model.NodesDictionary.Values)
        {
            if (Math.Abs(zCoordOfFace - node.Z) < 1E-9) zFaceNodes.Add(node);
        }

        var dirichletBCs = new List<NodalUnknownVariable>();
        foreach (var node in zFaceNodes)
        {
            dirichletBCs.Add(new NodalUnknownVariable(node, dofTypeToPrescribeDircihle, dirichleValue));
        }


        model.BoundaryConditions.Add(new ConvectionDiffusionBoundaryConditionSet(
            dirichletBCs,
            new INodalConvectionDiffusionNeumannBoundaryCondition[] { }
        ));

    }

    private void AddTopRightBackNodesBC(Model model, double boundaryCondition,
                                                        double modelMinX, double modelMaxX,
                                                        double modelMinY, double modelMaxY,
                                                        double modelMinZ, double modelMaxZ)
    {
        var topNodes = new List<INode>();
        var rightNodes = new List<INode>();
        var backNodes = new List<INode>();

        var tol = 1E-5;

        foreach (var node in model.NodesDictionary.Values)
        {
            if (Math.Abs(modelMaxZ - node.Z) < tol) topNodes.Add(node);
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

        model.BoundaryConditions.Add(new ConvectionDiffusionBoundaryConditionSet(dirichletBCs, new INodalConvectionDiffusionNeumannBoundaryCondition[] { }));
    }

    public void AddInitialConditions(Model model)
    {
        BoundaryAndInitialConditionsUtility.AssignConvectionDiffusionICToModel(model, InitialCondition);
    }

    public (IParentAnalyzer analyzer, ISolver solver, IChildAnalyzer loadcontrolAnalyzer) GetAppropriateSolverAnalyzerAndLog
    (Model model, double pseudoTimeStep, double pseudoTotalTime, int currentStep)
    {
        var solverFactory = new DenseMatrixSolver.Factory() { IsMatrixPositiveDefinite = false }; //Dense Matrix Solver solves with zero matrices!
                                                                                                  //var solverFactory = new SkylineSolver.Factory() { FactorizationPivotTolerance = 1e-8 };
        var algebraicModel = solverFactory.BuildAlgebraicModel(model);
        var solver = solverFactory.BuildSolver(algebraicModel);
        var provider = new ProblemConvectionDiffusion(model, algebraicModel);


        var linearAnalyzer = new LinearAnalyzer(algebraicModel, solver, provider);

        //var loadControlAnalyzerBuilder = new LoadControlAnalyzer.Builder(algebraicModel, solver, provider, numIncrements: 1)
        //{
        //    ResidualTolerance = 1E-8,
        //    MaxIterationsPerIncrement = 100,
        //    NumIterationsForMatrixRebuild = 1
        //};
        //var loadControlAnalyzer = loadControlAnalyzerBuilder.Build();

        //loadControlAnalyzer.TotalDisplacementsPerIterationLog = new TotalDisplacementsPerIterationLog(new List<(INode node, IDofType dof)>()
        //    {(model.NodesDictionary[MonitorNodeId], MonitorDOFType)}, algebraicModel);

        var analyzerBuilder = new NewmarkDynamicAnalyzer.Builder(algebraicModel, provider, linearAnalyzer, timeStep: pseudoTimeStep, totalTime: pseudoTotalTime, true, currentStep: currentStep);
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
}