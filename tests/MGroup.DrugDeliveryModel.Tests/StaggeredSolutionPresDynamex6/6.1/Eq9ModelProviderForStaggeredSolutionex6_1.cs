using System;
using System.Collections.Generic;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Solution;
using MGroup.MSolve.Discretization.Entities;
using MGroup.Constitutive.Structural;
using MGroup.Constitutive.Structural.BoundaryConditions;
using MGroup.Constitutive.Structural.Continuum;
using MGroup.Constitutive.Structural.Transient;
using MGroup.DrugDeliveryModel.Tests.Commons;
using MGroup.DrugDeliveryModel.Tests.Materials;
using MGroup.FEM.Structural.Continuum;
using MGroup.NumericalAnalyzers.Discretization.NonLinear;
using MGroup.NumericalAnalyzers.Dynamic;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.Solvers.Direct;
using MGroup.LinearAlgebra.Matrices;
using MGroup.Solvers.AlgebraicModel;
using MGroup.MSolve.Numerics.Integration.Quadratures;
using MGroup.MSolve.Numerics.Interpolation.Jacobians;
using MGroup.MSolve.Numerics.Interpolation;
using MGroup.MSolve.Solution.LinearSystem;
using System.Linq;
using MGroup.MSolve.Discretization;

namespace MGroup.DrugDeliveryModel.Tests.EquationModels
{
    public class Eq9ModelProviderForStaggeredSolutionex6_1
    {
        private GlobalAlgebraicModel<SkylineMatrix> algebraicModel;
        private IParentAnalyzer analyzer;
        //private double sc = 0.1;
        private double miNormal;// = 5;//KPa
        private double kappaNormal;// = 6.667; //Kpa
        private double miTumor;// = 22.44; //Kpa
        private double kappaTumor;// = 216.7; //Kpa
        private StructuralDof loadedDof;
        public double load_value { get; private set; }

        private double modelMinX;
        private double modelMaxX;
        private double modelMinY;
        private double modelMaxY;
        private double modelMinZ;
        private double modelMaxZ;

        private ComsolMeshReader reader;

        private Dictionary<int, double> lambda;
        Dictionary<int, double[][]> pressureTensorDivergenceAtElementGaussPoints;

        public int nodeIdToMonitor { get; private set; } //TODO put it where it belongs (coupled7and9eqsSolution.cs)
        private StructuralDof dofTypeToMonitor;



        public int loadedNode_Id { get; private set; }

        public Eq9ModelProviderForStaggeredSolutionex6_1(ComsolMeshReader comsolReader, double sc, double miNormal, double kappaNormal, double miTumor,
            double kappaTumor, double timeStep, double totalTime,
            Dictionary<int, double> lambda, Dictionary<int, double[][]> pressureTensorDivergenceAtElementGaussPoints,
            int nodeIdToMonitor, StructuralDof dofTypeToMonitor, StructuralDof loadedDof,
            double load_value, double modelMinX, double modelMaxX, double modelMinY, double modelMaxY, double modelMinZ, double modelMaxZ)
        {
            //this.sc = sc;
            this.miNormal = miNormal;
            this.kappaNormal = kappaNormal;
            this.miTumor = miTumor;
            this.kappaTumor = kappaTumor;
            this.pressureTensorDivergenceAtElementGaussPoints = pressureTensorDivergenceAtElementGaussPoints;
            this.lambda = lambda;
            reader = comsolReader;

            this.modelMinX = modelMinX;
            this.modelMaxX = modelMaxX;
            this.modelMinY = modelMinY;
            this.modelMaxY = modelMaxY;
            this.modelMinZ = modelMinZ;
            this.modelMaxZ = modelMaxZ;

            //log
            this.nodeIdToMonitor = nodeIdToMonitor;
            this.dofTypeToMonitor = dofTypeToMonitor;

            //Load
            this.loadedDof = loadedDof;
            this.load_value = load_value;

        }

        public Model GetModel()
        {
            var nodes = reader.NodesDictionary;
            var model = new Model();
            model.SubdomainsDictionary[0] = new Subdomain(id: 0);

            foreach (var node in nodes.Values)
            {
                model.NodesDictionary.Add(node.ID, node);
            }

            var materialNormal = new NeoHookeanMaterial3dJ3Isochoric(miNormal, kappaNormal);
            var materialTumor = new NeoHookeanMaterial3dJ3Isochoric(miTumor, kappaTumor);

            var elasticMaterial = new ElasticMaterial3D(youngModulus: 1, poissonRatio: 0.3);
            var DynamicMaterial = new TransientAnalysisProperties(density: 1, rayleighCoeffMass: 0, rayleighCoeffStiffness: 0);
            var elementFactory = new ContinuumElement3DFactory(elasticMaterial, DynamicMaterial);

            var volumeLoads = new List<IElementStructuralNeumannBoundaryCondition>();
            //var domains = new Dictionary<int, double[]>(2);
            foreach (var elementConnectivity in reader.ElementConnectivity)
            {
                var domainId = elementConnectivity.Value.Item3;
                var element = elementFactory.CreateNonLinearElementGrowt(elementConnectivity.Value.Item1, elementConnectivity.Value.Item2, domainId == 0 ? materialTumor : materialNormal, DynamicMaterial, lambda[elementConnectivity.Key]);
                //element.volumeForce = pressureTensorDivergenceAtElementGaussPoints[elementConnectivity.Key][0];
                var volumeForceX = new ElementDistributedLoad(element, StructuralDof.TranslationX, pressureTensorDivergenceAtElementGaussPoints[elementConnectivity.Key][0][0]);
                var volumeForceY = new ElementDistributedLoad(element, StructuralDof.TranslationY, pressureTensorDivergenceAtElementGaussPoints[elementConnectivity.Key][0][1]);
                var volumeForceZ = new ElementDistributedLoad(element, StructuralDof.TranslationZ, pressureTensorDivergenceAtElementGaussPoints[elementConnectivity.Key][0][2]);
                volumeLoads.AddRange(new[] { volumeForceX, volumeForceY, volumeForceZ });
                element.ID = elementConnectivity.Key;
                model.ElementsDictionary.Add(elementConnectivity.Key, element);
                model.SubdomainsDictionary[0].Elements.Add(element);
            }
            var modelNeumannConditions = new StructuralBoundaryConditionSet(null, null, null, volumeLoads, null, null);
            model.BoundaryConditions.Add(modelNeumannConditions);
            return model;
        }

        public void AddEq9ModelLoads(Model model)
        {
            INode maxDistanceNode = null;
            double currentMaxDistance = 0;
            foreach (INode node in model.NodesDictionary.Values)
            {
                double distance = Math.Sqrt(Math.Pow(node.X, 2) + Math.Pow(node.Y, 2) + Math.Pow(node.Z, 2));
                if (distance > currentMaxDistance)
                {
                    currentMaxDistance = distance;
                    maxDistanceNode = node;
                }
            }


            loadedNode_Id = maxDistanceNode.ID;
            //nodeIdToMonitor = loadedNode_Id;
            var loads = new List<INodalLoadBoundaryCondition>();

            loads.Add(new NodalLoad
            (
                maxDistanceNode,
                loadedDof,
                amount: load_value
            ));

            var emptyConstraints = new List<INodalDisplacementBoundaryCondition>();
            model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(emptyConstraints, loads));
        }




        public void AddBottomLeftRightFrontBackBCs(Model model)
        {
            var bottomNodes = new List<INode>();
            var leftNodes = new List<INode>();
            var rightNodes = new List<INode>();
            var frontNodes = new List<INode>();
            var backNodes = new List<INode>();

            var innerBulkNodes = new List<INode>();
            var tol = 1E-5;
            //Search for all boundary nodes
            foreach (var node in model.NodesDictionary.Values)
            {
                if (Math.Abs(modelMinZ - node.Z) < tol) bottomNodes.Add(node);
                if (Math.Abs(modelMinX - node.X) < tol) leftNodes.Add(node);
                if (Math.Abs(modelMaxX - node.X) < tol) rightNodes.Add(node);
                if (Math.Abs(modelMinY - node.Y) < tol) frontNodes.Add(node);
                if (Math.Abs(modelMaxY - node.Y) < tol) backNodes.Add(node);
            }

            foreach (var node in model.NodesDictionary.Values)
            {
                if (Math.Abs(modelMinZ - node.Z) < tol) { }

                else if (Math.Abs(modelMinX - node.X) < tol) { }
                else if (Math.Abs(modelMaxX - node.X) < tol) { }

                else if (Math.Abs(modelMinY - node.Y) < tol) { }
                else if (Math.Abs(modelMaxY - node.Y) < tol) { }
                else innerBulkNodes.Add(node);
            }

            var constraints = new List<INodalDisplacementBoundaryCondition>();

            //Apply roller constraint to bottom nodes. (constrained movement in z direction)
            foreach (var node in bottomNodes)
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationZ, amount: 0d));

            //Apply roller constraint to left nodes. (constrained movement in x direction)
            foreach (var node in leftNodes)
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationX, amount: 0d));

            //Apply roller constraint to right nodes. (constrained movement in x direction)
            foreach (var node in rightNodes)
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationX, amount: 0d));

            //Apply roller constraint to front nodes. (constrained movement in y direction)
            foreach (var node in frontNodes)
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationY, amount: 0d));

            //Apply roller constraint to back nodes. (constrained movement in y direction)
            foreach (var node in backNodes)
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationY, amount: 0d));

            var emptyloads = new List<INodalLoadBoundaryCondition>();
            model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(constraints, emptyloads));
        }

        public void AddBottomLeftFrontBackBCs(Model model)
        {
            var bottomNodes = new List<INode>();
            var leftNodes = new List<INode>();
            var frontNodes = new List<INode>();
            var backNodes = new List<INode>();

            var innerBulkNodes = new List<INode>();
            var tol = 1E-5;
            //Search for all boundary nodes
            foreach (var node in model.NodesDictionary.Values)
            {
                if (Math.Abs(modelMinZ - node.Z) < tol) bottomNodes.Add(node);
                if (Math.Abs(modelMinX - node.X) < tol) leftNodes.Add(node);
                if (Math.Abs(modelMinY - node.Y) < tol) frontNodes.Add(node);
                if (Math.Abs(modelMaxY - node.Y) < tol) backNodes.Add(node);
            }

            foreach (var node in model.NodesDictionary.Values)
            {
                if (Math.Abs(modelMinZ - node.Z) < tol) { }

                else if (Math.Abs(modelMinX - node.X) < tol) { }
                else if (Math.Abs(modelMaxX - node.X) < tol) { }

                else if (Math.Abs(modelMinY - node.Y) < tol) { }
                else if (Math.Abs(modelMaxY - node.Y) < tol) { }
                else innerBulkNodes.Add(node);
            }

            var constraints = new List<INodalDisplacementBoundaryCondition>();

            //Apply roller constraint to bottom nodes. (constrained movement in z direction)
            foreach (var node in bottomNodes)
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationZ, amount: 0d));

            //Apply roller constraint to left nodes. (constrained movement in x direction)
            foreach (var node in leftNodes)
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationX, amount: 0d));

            //Apply roller constraint to front nodes. (constrained movement in y direction)
            foreach (var node in frontNodes)
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationY, amount: 0d));

            //Apply roller constraint to back nodes. (constrained movement in y direction)
            foreach (var node in backNodes)
                constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationY, amount: 0d));

            var emptyloads = new List<INodalLoadBoundaryCondition>();
            model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(constraints, emptyloads));
        }



        public (IParentAnalyzer analyzer, ISolver solver, IChildAnalyzer loadcontrolAnalyzer) GetAppropriateSolverAnalyzerAndLog(Model model, double pseudoTimeStep, double pseudoTotalTime, int currentStep, int nIncrements)
        {
            var solverFactory = new SkylineSolver.Factory() { FactorizationPivotTolerance = 1e-8 };
            algebraicModel = solverFactory.BuildAlgebraicModel(model);
            var solver = solverFactory.BuildSolver(algebraicModel);
            var provider = new ProblemStructural(model, algebraicModel);
            var loadControlAnalyzerBuilder = new LoadControlAnalyzer.Builder(algebraicModel, solver, provider, numIncrements: nIncrements)
            {
                ResidualTolerance = 1E-8,
                MaxIterationsPerIncrement = 100,
                NumIterationsForMatrixRebuild = 1
            };
            var nlAnalyzer = loadControlAnalyzerBuilder.Build();
            var loadControlAnalyzer = (LoadControlAnalyzer)nlAnalyzer;
            loadControlAnalyzer.TotalDisplacementsPerIterationLog = new TotalDisplacementsPerIterationLog(
                new List<(INode node, IDofType dof)>()
                {
                    (model.NodesDictionary[nodeIdToMonitor], dofTypeToMonitor),


                }, algebraicModel
            );

            //var analyzer = (new PseudoTransientAnalyzer.Builder(algebraicModel, provider, loadControlAnalyzer, timeStep: pseudoTimeStep, totalTime: pseudoTotalTime, currentStep: currentStep)).Build();
            var analyzerBuilder = new NewmarkDynamicAnalyzer.Builder(algebraicModel, provider, loadControlAnalyzer, timeStep: pseudoTimeStep, totalTime: pseudoTotalTime, false, currentStep: currentStep);
            analyzerBuilder.SetNewmarkParametersForConstantAcceleration();
            analyzer = analyzerBuilder.Build();


            //Sparse tet Mesh
            var watchDofs = new[]
            {
                new List<(INode node, IDofType dof)>()
                {
                    (model.NodesDictionary[nodeIdToMonitor], dofTypeToMonitor),
                    (model.NodesDictionary[nodeIdToMonitor], StructuralDof.TranslationY),
                    (model.NodesDictionary[nodeIdToMonitor], StructuralDof.TranslationZ),
                }
            };

            loadControlAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs[0], algebraicModel);

            return (analyzer, solver, loadControlAnalyzer);
        }
        public Dictionary<int, double[][]> GetVelocities()
        {
            var nodalDofs = new StructuralDof[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
            var modelVelocities = analyzer.CreateState().StateVectors["First order derivative of solution"].Copy();
            modelVelocities.CheckForCompatibility = false;


            IIsoparametricInterpolation3D interpolation;
            IQuadrature3D quadrature;

            int nGaussPoints = 1;
            var velcocities = new Dictionary<int, double[][]>();

            foreach (var elem in reader.ElementConnectivity)
            {
                var nodes = elem.Value.Item2;
                var nodalVelocities = new double[3 * nodes.Count()];


                for (int i = 0; i < nodes.Length; i++)
                {
                    for (int i1 = 0; i1 < 3; i1++)
                    {
                        var dof = nodalDofs[i1];
                        double dofVelocity = 0;
                        try
                        {
                            dofVelocity = algebraicModel.ExtractSingleValue(modelVelocities, nodes[i], dof);
                        }
                        catch (KeyNotFoundException e)
                        {
                            // recover from exception
                        }

                        nodalVelocities[3 * i + i1] = dofVelocity;
                    }
                }


                var velocity = new double[nGaussPoints][];

                int gausspoinNo = 0; // anti gia loop to ntegation.Gausspoint1

                if (elem.Value.Item1 == CellType.Tet4)
                {
                    interpolation = InterpolationTet4.UniqueInstance;
                    quadrature = TetrahedronQuadrature.Order1Point1;
                }
                else if (elem.Value.Item1 == CellType.Hexa8)
                {
                    interpolation = InterpolationHexa8.UniqueInstance;
                    quadrature = GaussLegendre3D.GetQuadratureWithOrder(2, 2, 2);
                }
                else if (elem.Value.Item1 == CellType.Wedge6)
                {
                    interpolation = InterpolationWedge6.UniqueInstance;
                    quadrature = WedgeQuadrature.Points8;
                }
                else throw new ArgumentException("Wrong Cell Type");

                var shapeFunctionValues = interpolation.EvaluateFunctionsAt(quadrature.IntegrationPoints[gausspoinNo]);


                velocity[0] = new double[] { 0, 0, 0 };

                for (int i1 = 0; i1 < shapeFunctionValues.Length; i1++)
                {
                    velocity[0][0] += shapeFunctionValues[i1] * nodalVelocities[3 * i1 + 0];
                    velocity[0][1] += shapeFunctionValues[i1] * nodalVelocities[3 * i1 + 1];
                    velocity[0][2] += shapeFunctionValues[i1] * nodalVelocities[3 * i1 + 2];

                }

                velcocities[elem.Key] = velocity;



            }

            modelVelocities.CheckForCompatibility = false;
            return velcocities;
        }

        private IGlobalVector modelVelocities;
        public Dictionary<int, double[][]> GetVelocities2()
        {
            var nodalDofs = new StructuralDof[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
            modelVelocities = analyzer.CreateState().StateVectors["First order derivative of solution"].Copy();
            modelVelocities.CheckForCompatibility = false;
            var velocityNodalResults = algebraicModel.ExtractAllResults(modelVelocities);
            var velocityNodalData = velocityNodalResults.Data;

            IIsoparametricInterpolation3D interpolation;
            IQuadrature3D quadrature;

            int nGaussPoints = 1;
            var velcocities = new Dictionary<int, double[][]>();

            foreach (var elem in reader.ElementConnectivity)
            {
                var nodes = elem.Value.Item2;
                var nodalVelocities = new double[3 * nodes.Count()];

                if (elem.Value.Item1 == CellType.Tet4)
                {
                    interpolation = InterpolationTet4.UniqueInstance;
                    quadrature = TetrahedronQuadrature.Order1Point1;
                }
                else if (elem.Value.Item1 == CellType.Hexa8)
                {
                    interpolation = InterpolationHexa8.UniqueInstance;
                    quadrature = GaussLegendre3D.GetQuadratureWithOrder(2, 2, 2);
                }
                else if (elem.Value.Item1 == CellType.Wedge6)
                {
                    interpolation = InterpolationWedge6.UniqueInstance;
                    quadrature = WedgeQuadrature.Points8;
                }
                else throw new ArgumentException("Wrong Cell Type");

                for (int i = 0; i < nodes.Length; i++)
                {
                    for (int i1 = 0; i1 < 3; i1++)
                    {
                        var dof = nodalDofs[i1];
                        double dofVelocity = 0;

                        bool foundVelocity = velocityNodalData.TryGetValue(nodes[i].ID, i1, out double velocityVal);
                        if (foundVelocity) { dofVelocity = velocityVal; }


                        nodalVelocities[3 * i + i1] = dofVelocity;
                    }
                }


                var velocity = new double[nGaussPoints][];

                int gausspoinNo = 0; // anti gia loop to ntegation.Gausspoint1
                var shapeFunctionValues = interpolation.EvaluateFunctionsAt(quadrature.IntegrationPoints[gausspoinNo]);


                velocity[gausspoinNo] = new double[] { 0, 0, 0 };

                for (int i1 = 0; i1 < shapeFunctionValues.Length; i1++)
                {
                    velocity[gausspoinNo][0] += shapeFunctionValues[i1] * nodalVelocities[3 * i1 + 0];
                    velocity[gausspoinNo][1] += shapeFunctionValues[i1] * nodalVelocities[3 * i1 + 1];
                    velocity[gausspoinNo][2] += shapeFunctionValues[i1] * nodalVelocities[3 * i1 + 2];

                }

                velcocities[elem.Key] = velocity;



            }

            modelVelocities.CheckForCompatibility = false;
            return velcocities;
        }

        public Dictionary<int, double[]> GetVelocityDIV()
        {
            var nodalDofs = new StructuralDof[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
            //var modelVelocities = analyzer.CreateState().StateVectors["First order derivative of solution"].Copy();
            modelVelocities.CheckForCompatibility = false;
            var velocityNodalResults = algebraicModel.ExtractAllResults(modelVelocities);
            var velocityNodalData = velocityNodalResults.Data;

            IIsoparametricInterpolation3D interpolation;
            IQuadrature3D quadrature;

            int nGaussPoints = 1;
            var velcocities = new Dictionary<int, double[]>();

            foreach (var elem in reader.ElementConnectivity)
            {
                var nodes = elem.Value.Item2;
                var nodalVelocities = new double[3 * nodes.Count()];


                for (int i = 0; i < nodes.Length; i++)
                {
                    for (int i1 = 0; i1 < 3; i1++)
                    {
                        var dof = nodalDofs[i1];
                        double dofVelocity = 0;

                        bool foundVelocity = velocityNodalData.TryGetValue(nodes[i].ID, i1, out double velocityVal);
                        if (foundVelocity) { dofVelocity = velocityVal; }


                        nodalVelocities[3 * i + i1] = dofVelocity;
                    }
                }


                var velocityDiv = new double[nGaussPoints];

                int gausspoinNo = 0; // anti gia loop to ntegation.Gausspoint1
                                     //IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives;

                if (elem.Value.Item1 == CellType.Tet4)
                {
                    interpolation = InterpolationTet4.UniqueInstance;
                    quadrature = TetrahedronQuadrature.Order1Point1;
                }
                else if (elem.Value.Item1 == CellType.Hexa8)
                {
                    interpolation = InterpolationHexa8.UniqueInstance;
                    quadrature = GaussLegendre3D.GetQuadratureWithOrder(2, 2, 2);
                }
                else if (elem.Value.Item1 == CellType.Wedge6)
                {
                    interpolation = InterpolationWedge6.UniqueInstance;
                    quadrature = WedgeQuadrature.Points8;
                }
                else throw new ArgumentException("Wrong Cell Type");

                var shapeFunctionNaturalDerivatives = interpolation.EvaluateNaturalGradientsAt(quadrature.IntegrationPoints[gausspoinNo]);
                var jacobian = new IsoparametricJacobian3D(nodes, shapeFunctionNaturalDerivatives);
                var jacobianInverse = jacobian.InverseMatrix.Transpose();


                double[,] dvi_dnaturalj = new double[3, 3]; //{ dphi_dksi, dphi_dheta, dphi_dzeta}
                for (int i1 = 0; i1 < shapeFunctionNaturalDerivatives.NumRows; i1++)
                {
                    dvi_dnaturalj[0, 0] += shapeFunctionNaturalDerivatives[i1, 0] * nodalVelocities[3 * i1 + 0];
                    dvi_dnaturalj[0, 1] += shapeFunctionNaturalDerivatives[i1, 1] * nodalVelocities[3 * i1 + 0];
                    dvi_dnaturalj[0, 2] += shapeFunctionNaturalDerivatives[i1, 2] * nodalVelocities[3 * i1 + 0];

                    dvi_dnaturalj[1, 0] += shapeFunctionNaturalDerivatives[i1, 0] * nodalVelocities[3 * i1 + 1];
                    dvi_dnaturalj[1, 1] += shapeFunctionNaturalDerivatives[i1, 1] * nodalVelocities[3 * i1 + 1];
                    dvi_dnaturalj[1, 2] += shapeFunctionNaturalDerivatives[i1, 2] * nodalVelocities[3 * i1 + 1];

                    dvi_dnaturalj[2, 0] += shapeFunctionNaturalDerivatives[i1, 0] * nodalVelocities[3 * i1 + 2];
                    dvi_dnaturalj[2, 1] += shapeFunctionNaturalDerivatives[i1, 1] * nodalVelocities[3 * i1 + 2];
                    dvi_dnaturalj[2, 2] += shapeFunctionNaturalDerivatives[i1, 2] * nodalVelocities[3 * i1 + 2];
                }

                var dvi_dnaturaljMAT = Matrix.CreateFromArray(dvi_dnaturalj);

                var dvi_dcartesianj = dvi_dnaturaljMAT * jacobianInverse.Transpose();
                velocityDiv[gausspoinNo] = dvi_dcartesianj[0, 0] + dvi_dcartesianj[1, 1] + dvi_dcartesianj[2, 2];

                velcocities[elem.Key] = velocityDiv;



            }

            modelVelocities.CheckForCompatibility = false;
            return velcocities;
        }
    }
}
