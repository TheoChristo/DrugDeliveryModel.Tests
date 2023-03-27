using System;
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

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
    public class CoxModelBuilder
    {
        /// <summary>
        /// The average value of the three components of the fluid velocity vector  [m/s]
        /// </summary>
        public readonly Dictionary<int, double[]> FluidSpeed; // [m/s]

        /// <summary>
        /// Diffusivity of oxygen [m2/s]
        /// </summary>
        private readonly double Dox; // [m2/s]

        /// <summary>
        /// Oxygen uptake [mol/(m3*s)]
        /// </summary>
        private readonly double Aox; // [mol/(m3*s)]

        /// <summary>
        /// Oxygen uptake [mol/m3]
        /// </summary>
        private readonly double Kox; // [mol / m3]

        /// <summary>
        /// Oxygen permeability across tumor vessel walls [m/s]
        /// </summary>
        private readonly double PerOx; // [m/s]

        /// <summary>
        /// Vascular Density [1/m]
        /// </summary>
        private readonly double Sv; // [1/m]

        /// <summary>
        /// Initial Oxygen Concentration [mol/m3]
        /// </summary>
        private readonly double CInitOx; // [mol/m3]

        /// <summary>
        /// Cancer cell density [1]
        /// </summary>
        private readonly Dictionary<int, double> T; // [cells]

        /// <summary>
        /// Term without non-linear term
        /// </summary>
        private readonly Func<double> independentLinearSource;

        private readonly Dictionary<int, Func<double, double>> ProductionFuncWithoutConstantTerm;

        private readonly Dictionary<int, Func<double, double>> ProductionFuncWithoutConstantTermDDerivative;

        private readonly ComsolMeshReader mesh;

        /// <summary>
        /// List containing the DIRICHLET boundary conditions for the Convection Diffusion problem.
        /// Item1 : Boundary condition case with respect to the face of the domain (LeftDirichlet, TopDirichlet etc).
        /// Item2 : An StructuralDof array containing the DOFs that are constrained.
        /// Item3 : A jagged array of arrays that contain all the coordinate sets of the bounded dofs. The lenght of
        ///         array is equal to the number of constrained dofs. Each array contains the coordinates of the constrained
        ///         dofs.
        /// Item3 : A double array containing the values of the constrained dofs (1-1 correspondence with the dofs in Item2).
        /// </summary>
        private List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> convectionDiffusionDirichletBC;

        /// <summary>
        /// List containing the NEUMANN boundary conditions for the Convection Diffusion problem.
        /// Item1 : Boundary condition case with respect to the face of the domain (RightPointFlux, TopDistributedFlux etc)
        /// Item2 : An StructuralDof array containing the information about the direction of the dofs where the force is
        ///         applied.
        /// Item3 : A jagged array of arrays that contain all the coordinate sets of the bounded dofs. The lenght of
        ///         array is equal to the number of constrained dofs. Each array contains the coordinates of the constrained
        ///         dofs
        /// Item3 : A double array containing the values of the constrained dofs (1-1 correspondence with the dofs in Item2)
        /// </summary>
        private List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> convectionDiffusionNeumannBC;

        private double initialCondition;

        public Dictionary<int, double[]> div_vs { get; set; }

        private int nodeIdToMonitor; //TODO put it where it belongs (coupled7and9eqsSolution.cs)

        private ConvectionDiffusionDof dofTypeToMonitor = ConvectionDiffusionDof.UnknownVariable;


        public CoxModelBuilder(ComsolMeshReader modelReader,
            Dictionary<int, double[]> FluidSpeed, double Dox, double Aox, double Kox, double PerOx, double Sv, double CInitOx, Dictionary<int, double> T, double initialCondition,
            Func<double> independentLinearSource, Dictionary<int, Func<double, double>> ProductionFuncWithoutConstantTerm, Dictionary<int, Func<double, double>> ProductionFuncWithoutConstantTermDDerivative,
            int nodeIdToMonitor, ConvectionDiffusionDof dofTypeToMonitor,
            List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> convectionDiffusionDirichletBC,
            List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> convectionDiffusionNeumannBC )
        {
            this.mesh = modelReader;
            this.FluidSpeed = FluidSpeed;
            this.Dox = Dox;
            this.Aox = Aox;
            this.Kox = Kox;
            this.PerOx = PerOx;
            this.Sv = Sv;
            this.CInitOx = CInitOx;
            this.T = T;
            this.initialCondition = initialCondition;

            this.independentLinearSource = independentLinearSource;
            this.ProductionFuncWithoutConstantTerm = ProductionFuncWithoutConstantTerm;
            this.ProductionFuncWithoutConstantTermDDerivative = ProductionFuncWithoutConstantTermDDerivative;

            this.mesh = modelReader;
            IsoparametricJacobian3D.DeterminantTolerance = 1e-30;

            this.convectionDiffusionDirichletBC = convectionDiffusionDirichletBC;
            this.convectionDiffusionNeumannBC = convectionDiffusionNeumannBC;

            this.initialCondition = initialCondition;

            //log
            this.nodeIdToMonitor = nodeIdToMonitor;
            this.dofTypeToMonitor = dofTypeToMonitor;

        }

        public Model GetModel()
        {
            var capacity = 1;
            var diffusionCoefficient = Dox;
            var independentSourceCoefficient = independentLinearSource();
            var dependentSourceCoefficient = 0;

            //Assign equation properties to the domain elements
            var convectionDomainCoefficients = new Dictionary<int, double[]>();
            var dependentProductionCoefficients = new Dictionary<int, double>();
            var independentProductionCoefficients = new Dictionary<int, double>();

            foreach (var elementConnectivity in mesh.ElementConnectivity)
            {
                convectionDomainCoefficients[elementConnectivity.Key] = FluidSpeed[elementConnectivity.Key];
                dependentProductionCoefficients[elementConnectivity.Key] = dependentSourceCoefficient;
                independentProductionCoefficients[elementConnectivity.Key] = independentSourceCoefficient;
            }

            //Create Model
            var modelProvider = new GenericComsol3DConvectionDiffusionProductionModelProviderDistributedSpace(mesh);
            var model = modelProvider.CreateModelFromComsolFile(convectionDomainCoefficients, diffusionCoefficient,
                dependentProductionCoefficients, independentProductionCoefficients, capacity,
                ProductionFuncWithoutConstantTerm, ProductionFuncWithoutConstantTermDDerivative);
            return model;
        }

        public void AddBoundaryConditions(Model model)
        {
            BoundaryAndInitialConditionsUtility.AssignConvectionDiffusionDirichletBCsToModel(model, convectionDiffusionDirichletBC, 1e-3);
            BoundaryAndInitialConditionsUtility.AssignConvectionDiffusionICToModel(model, initialCondition);
        }

        public (IParentAnalyzer analyzer, ISolver solver, IChildAnalyzer loadcontrolAnalyzer) GetAppropriateSolverAnalyzerAndLog(Model model, double TimeStep, double TotalTime, int currentStep, int nIncrements)
        {
            var solverFactory = new DenseMatrixSolver.Factory() { IsMatrixPositiveDefinite = false };
            var algebraicModel = solverFactory.BuildAlgebraicModel(model);
            var solver = solverFactory.BuildSolver(algebraicModel);
            var problem = new ProblemConvectionDiffusion(model, algebraicModel);

            var loadControlAnalyzerBuilder = new LoadControlAnalyzer.Builder(algebraicModel, solver, problem, numIncrements: 1)
            {
                ResidualTolerance = 1E-3,
                MaxIterationsPerIncrement = 100,
                NumIterationsForMatrixRebuild = 1
            };

            var loadControlAnalyzer = loadControlAnalyzerBuilder.Build();

            var dynamicAnalyzerBuilder = new NewmarkDynamicAnalyzer.Builder(algebraicModel, problem, loadControlAnalyzer, timeStep: TimeStep, totalTime: TotalTime, false, currentStep: currentStep);
            var dynamicAnalyzer = dynamicAnalyzerBuilder.Build();

            // Create a log for the desired dof
            var watchDofs = new List<(INode node, IDofType dof)>()
            {
                (model.NodesDictionary[nodeIdToMonitor], dofTypeToMonitor),
            };

            loadControlAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs, algebraicModel);
            dynamicAnalyzer.ResultStorage = new ImplicitIntegrationAnalyzerLog();

            return (dynamicAnalyzer, solver, loadControlAnalyzer);
        }
    }
}
