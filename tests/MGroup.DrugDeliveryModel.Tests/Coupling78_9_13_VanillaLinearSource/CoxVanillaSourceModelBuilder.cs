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
using BC = MGroup.DrugDeliveryModel.Tests.Commons.BoundaryAndInitialConditionsUtility.BoundaryConditionCase;
using MGroup.LinearAlgebra.Matrices;
using MGroup.Solvers.AlgebraicModel;
using MGroup.MSolve.Numerics.Integration.Quadratures;
using MGroup.MSolve.Numerics.Interpolation;
using MGroup.MSolve.Discretization;

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
    public class CoxVanillaSourceModelBuilder
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
        //private readonly double CInitOx; // [mol/m3]

        /// <summary>
        /// Cancer cell density [1]
        /// </summary>
        private readonly Dictionary<int, double> T; // [cells]

        /// <summary>
        /// Term without non-linear term
        /// </summary>
        private readonly Func<double> independentLinearSource;

        private readonly Func<double> dependentLinearSource;


        private readonly ComsolMeshReader mesh;

        public GlobalAlgebraicModel<Matrix> algebraicModel;

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


        public CoxVanillaSourceModelBuilder(ComsolMeshReader modelReader,
            Dictionary<int, double[]> FluidSpeed, Func<double> independentLinearSource, Func<double> dependentLinearSource,
            double Dox, double Aox, double Kox, double PerOx, double Sv, double initialCondition,
            int nodeIdToMonitor, ConvectionDiffusionDof dofTypeToMonitor,
            List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> convectionDiffusionDirichletBC,
            List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> convectionDiffusionNeumannBC)
        {
            this.mesh = modelReader;
            this.FluidSpeed = FluidSpeed;
            this.Dox = Dox;
            this.Aox = Aox;
            this.Kox = Kox;
            this.PerOx = PerOx;
            this.Sv = Sv;
            this.initialCondition = initialCondition;

            this.independentLinearSource = independentLinearSource;
            this.dependentLinearSource = dependentLinearSource;

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
            var dependentSourceCoefficient = dependentLinearSource();

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
                dependentProductionCoefficients, independentProductionCoefficients, capacity);
            return model;
        }

        //public void AddBoundaryConditions(Model model)
        //{
        //    BoundaryAndInitialConditionsUtility.AssignConvectionDiffusionDirichletBCsToModel(model, convectionDiffusionDirichletBC, 1e-3);
        //    BoundaryAndInitialConditionsUtility.AssignConvectionDiffusionICToModel(model, initialCondition);
        //}

        public void AddBoundaryConditions(Model model)
        {
            foreach (var BCdata in convectionDiffusionDirichletBC)
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
                    (model.NodesDictionary[nodeIdToMonitor], dofTypeToMonitor),
                }
            };
            linearAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs[0], algebraicModel);

            return (analyzer, solver, linearAnalyzer);
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

            var analyzerBuilder = new NewmarkDynamicAnalyzer.Builder(algebraicModel, provider, linearAnalyzer, timeStep: pseudoTimeStep, totalTime: pseudoTotalTime, false, currentStep: currentStep);
            analyzerBuilder.SetNewmarkParametersForConstantAcceleration();
            var analyzer = analyzerBuilder.Build();
            var watchDofs = new[]
            {
                    new List<(INode node, IDofType dof)>()
                    {
                        (model.NodesDictionary[nodeIdToMonitor], dofTypeToMonitor),
                    }
                };
            linearAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs[0], algebraicModel);

            return (analyzer, solver, linearAnalyzer);
        }

        public void UpdateGausspointValuesOfElements(Dictionary<int, double> domainCox, ISolver solver, IChildAnalyzer childAnalyzer, Model model, GlobalAlgebraicModel<Matrix> algebraicModel)
        {
            var p = childAnalyzer.CurrentAnalysisResult;
            IIsoparametricInterpolation3D interpolation;
            IQuadrature3D quadrature;
            foreach (var elem in model.ElementsDictionary.Values)
            {
                if (mesh.ElementConnectivity[elem.ID].Item1 == CellType.Tet4)
                {
                    interpolation = InterpolationTet4.UniqueInstance;
                    quadrature = TetrahedronQuadrature.Order1Point1;
                }
                else if (mesh.ElementConnectivity[elem.ID].Item1 == CellType.Hexa8)
                {
                    interpolation = InterpolationHexa8.UniqueInstance;
                    quadrature = GaussLegendre3D.GetQuadratureWithOrder(2, 2, 2);
                }
                else if (mesh.ElementConnectivity[elem.ID].Item1 == CellType.Wedge6)
                {
                    interpolation = InterpolationWedge6.UniqueInstance;
                    quadrature = WedgeQuadrature.Points8;
                }
                else throw new ArgumentException("Wrong Cell Type");
                var elemSoultion = algebraicModel.ExtractElementVector(p, elem);
                domainCox[elem.ID] = UpdateGauspointQuantitiesOfElement(elemSoultion, interpolation, quadrature, elem);
            }

        }

        private double UpdateGauspointQuantitiesOfElement(double[] localDisplacements, IIsoparametricInterpolation3D Interpolation,
         IQuadrature3D quadratureForMass, IElementType element)
        {
            var shapeFunctions = Interpolation.EvaluateFunctionsAtGaussPoints(quadratureForMass);


            //double[][] pressureTensorDivergenceAtGaussPoints = new double[quadratureForMass.IntegrationPoints.Count][];
            //IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives;
            //shapeFunctionNaturalDerivatives = Interpolation.EvaluateNaturalGradientsAtGaussPoints(quadratureForMass);
            //var jacobians = shapeFunctionNaturalDerivatives.Select(x => new IsoparametricJacobian3D(element.Nodes, x));
            //Matrix[] jacobianInverse = jacobians.Select(x => x.InverseMatrix.Transpose()).ToArray();
            //for (int gp = 0; gp < quadratureForMass.IntegrationPoints.Count; ++gp)
            //{
            //    double[] dphi_dnatural = new double[3]; //{ dphi_dksi, dphi_dheta, dphi_dzeta}
            //    for (int i1 = 0; i1 < shapeFunctionNaturalDerivatives[gp].NumRows; i1++)
            //    {
            //        dphi_dnatural[0] += shapeFunctionNaturalDerivatives[gp][i1, 0] * localDisplacements[i1];
            //        dphi_dnatural[1] += shapeFunctionNaturalDerivatives[gp][i1, 1] * localDisplacements[i1];
            //        dphi_dnatural[2] += shapeFunctionNaturalDerivatives[gp][i1, 2] * localDisplacements[i1];
            //    }

            //    var dphi_dnaturalMAT = Matrix.CreateFromArray(dphi_dnatural, 1, 3);

            //    var dphi_dcartesian = dphi_dnaturalMAT * jacobianInverse[gp].Transpose();

            //    pressureTensorDivergenceAtGaussPoints[gp] = new double[3] { dphi_dcartesian[0, 0], dphi_dcartesian[0, 1], dphi_dcartesian[0, 2] };

            //}

            double elementCox = 0;

            for (int i1 = 0; i1 < shapeFunctions[0].Length; i1++)
            {
                elementCox += shapeFunctions[0][i1] * localDisplacements[i1];
            }

            //an exoume perissotera GP vazoume for (int gp = 0; gp < quadratureForMass.IntegrationPoints.Count; ++gp)

            return elementCox;
        }
    }
}
