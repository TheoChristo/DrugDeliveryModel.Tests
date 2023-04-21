using System;
using System.Collections.Generic;
using System.Linq;
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
using MGroup.MSolve.Numerics.Integration;
using MGroup.LinearAlgebra.Matrices;
using MGroup.Constitutive.ConvectionDiffusion.BoundaryConditions;
using MGroup.Constitutive.ConvectionDiffusion;
using MGroup.Solvers.AlgebraicModel;
using MGroup.MSolve.Numerics.Interpolation;
using MGroup.MSolve.Numerics.Integration.Quadratures;
using MGroup.MSolve.Solution.LinearSystem;
using TriangleNet.Tools;
using MGroup.MSolve.Numerics.Interpolation.Jacobians;
using MGroup.MSolve.Discretization.BoundaryConditions;

namespace MGroup.DrugDeliveryModel.Tests.EquationModels
{
	public class Eq9ModelProvider
    {
        private GlobalAlgebraicModel<SkylineMatrix> algebraicModel;
        private IParentAnalyzer analyzer;

        //TODO Orestis :if AddLoads9BCs() is implemented in a right way these will not be necessary and be deleted.
        //TODO Orestis :if AddLoads9BCs() is implemented in a right way these will not be necessary and be deleted.
        //TODO Orestis :if AddEquation9BCs() is implemented in a right way these will not be necessary and be deleted.
        // Model Min,Max(X,Y,Z) Deleted and removed from constructor
        // load_value is deleted and removed from constructor
        // loadedDof is deleted and removed from constructor

        //private double sc = 0.1;
        private double miNormal;// = 5;//KPa
        private double kappaNormal;// = 6.667; //Kpa
        private double miTumor;// = 22.44; //Kpa
        private double kappaTumor;// = 216.7; //Kpa
        
        /// <summary>
        /// List containing the DIRICHLET boundary conditions for the structural problem
        /// Item1 : Boundary condition case with respect to the face of the domain (LeftDirichlet, TopDirichlet etc)
        /// Item2 : An StructuralDof array containing the DOFs that are constrained
        /// Item3 : A jagged array of arrays that contain all the coordinate sets of the bounded dofs. The lenght of
        ///         array is equal to the number of constrained dofs. Each array contains the coordinates of the constrained
        ///         dofs
        /// Item3 : A double array containing the values of the constrained dofs (1-1 correspondence with the dofs in Item2)
        /// </summary>
        private List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, StructuralDof[], double[][], double[])> structuralDirichletBC;
        
        /// <summary>
        /// List containing the NEUMANN boundary conditions for the structural problem
        /// Item1 : Boundary condition case with respect to the face of the domain (RightPointFlux, TopDistributedFlux etc)
        /// Item2 : An StructuralDof array containing the information about the direction of the dofs where the force is
        ///         applied.
        /// Item3 : A jagged array of arrays that contain all the coordinate sets of the bounded dofs. The lenght of
        ///         array is equal to the number of constrained dofs. Each array contains the coordinates of the constrained
        ///         dofs
        /// Item3 : A double array containing the values of the constrained dofs (1-1 correspondence with the dofs in Item2)
        /// </summary>
        private List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, StructuralDof[], double[][], double[])> structuralNeumannBC;
        
        private ComsolMeshReader reader;
       
        private Dictionary<int, double> lambda;
        
        Dictionary<int, double[][]> pressureTensorDivergenceAtElementGaussPoints;

        public int nodeIdToMonitor { get; private set; } //TODO put it where it belongs (coupled7and9eqsSolution.cs)
        
        private StructuralDof dofTypeToMonitor;
        
        Dictionary<int, double[]> elementslastConvergedDisplacements;

        private bool elementSavedDisplacementsIsInitialized = false;
        
        private double density;

        
        //TODO Orestis : OLd Constructor is deleted.
 
        public Eq9ModelProvider(
            ComsolMeshReader comsolReader, double sc, double miNormal, double kappaNormal, double miTumor,
            double kappaTumor, double density,
            double timeStep, double totalTime,
            Dictionary<int, double> lambda, Dictionary<int, double[][]> pressureTensorDivergenceAtElementGaussPoints,
            int nodeIdToMonitor, StructuralDof dofTypeToMonitor,
            List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, StructuralDof[], double[][], double[])> structuralNeumannBC,
            List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, StructuralDof[], double[][], double[])> structuralDirichletBC
            )
        {
            //this.sc = sc;
            this.miNormal = miNormal;
            this.kappaNormal = kappaNormal;
            this.miTumor = miTumor;
            this.kappaTumor = kappaTumor;
            this.pressureTensorDivergenceAtElementGaussPoints = pressureTensorDivergenceAtElementGaussPoints;
            this.lambda = lambda;
            reader = comsolReader;

            //log
            this.nodeIdToMonitor = nodeIdToMonitor;
            this.dofTypeToMonitor = dofTypeToMonitor;

            this.structuralNeumannBC = structuralNeumannBC;
            this.structuralDirichletBC = structuralDirichletBC;

            this.density = density;
            
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
            var dynamicMaterial = new TransientAnalysisProperties(density: density, rayleighCoeffMass: 0, rayleighCoeffStiffness: 0);
            var elementFactory = new ContinuumElement3DFactory(elasticMaterial, dynamicMaterial);

            var volumeLoads = new List<IElementStructuralNeumannBoundaryCondition>();

            //var domains = new Dictionary<int, double[]>(2);
            foreach (var elementConnectivity in reader.ElementConnectivity)
            {
                var domainId = elementConnectivity.Value.Item3;
                var element = elementFactory.CreateNonLinearElementGrowt(elementConnectivity.Value.Item1, elementConnectivity.Value.Item2, domainId == 0 ? materialTumor : materialNormal, dynamicMaterial, lambda[elementConnectivity.Key]);
                //element.volumeForce = pressureTensorDivergenceAtElementGaussPoints[elementConnectivity.Key][0];
                var volumeForceX = new ElementDistributedLoad(element, StructuralDof.TranslationX, pressureTensorDivergenceAtElementGaussPoints[elementConnectivity.Key][0][0]);
                var volumeForceY = new ElementDistributedLoad(element, StructuralDof.TranslationY, pressureTensorDivergenceAtElementGaussPoints[elementConnectivity.Key][0][1]);
                var volumeForceZ = new ElementDistributedLoad(element, StructuralDof.TranslationZ, pressureTensorDivergenceAtElementGaussPoints[elementConnectivity.Key][0][2]);
                volumeLoads.AddRange(new [] {volumeForceX, volumeForceY, volumeForceZ});
                element.ID = elementConnectivity.Key;
                //if (elementSavedDisplacementsIsInitialized) { element.lastConvergedDisplacements = elementslastConvergedDisplacements[element.ID]; }
                model.ElementsDictionary.Add(elementConnectivity.Key, element);
                model.SubdomainsDictionary[0].Elements.Add(element);
            }

            var modelNeumannConditions = new StructuralBoundaryConditionSet(null,null,null,volumeLoads,null,null);
            model.BoundaryConditions.Add(modelNeumannConditions);

            return model;
        }
        
        public void AddBoundaryConditions(Model model)
        {
            BoundaryAndInitialConditionsUtility.AssignStructuralDirichletBCsToModel(model, structuralDirichletBC, 1e-3);
            BoundaryAndInitialConditionsUtility.AssignStructuralNeumannBCsToModel(model, structuralNeumannBC, 1E-3);
        }
        
        //TODO Gerasimos add if for dynamic or peudostatic analyzer

        public (IParentAnalyzer analyzer, ISolver solver, IChildAnalyzer loadcontrolAnalyzer) GetAppropriateSolverAnalyzerAndLog(Model model, double pseudoTimeStep, double pseudoTotalTime, int currentStep, int nIncrements)
        {
            var solverFactory = new SkylineSolver.Factory() { FactorizationPivotTolerance = 1e-8 };
            algebraicModel =  solverFactory.BuildAlgebraicModel(model);
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
                },
                algebraicModel
            );

            //var analyzer = (new PseudoTransientAnalyzer.Builder(algebraicModel, provider, loadControlAnalyzer, timeStep: pseudoTimeStep, totalTime: pseudoTotalTime, currentStep: currentStep)).Build();
            var analyzerBuilder = new NewmarkDynamicAnalyzer.Builder(algebraicModel, provider, loadControlAnalyzer, timeStep: pseudoTimeStep, totalTime: pseudoTotalTime,false, currentStep: currentStep);
            analyzerBuilder.SetNewmarkParametersForConstantAcceleration();
            //var analyzerBuilder = new BDFDynamicAnalyzer.Builder(algebraicModel, provider, loadControlAnalyzer, timeStep: pseudoTimeStep, totalTime: pseudoTotalTime, currentTimeStep: currentStep, bdfOrder: 5);
            analyzer = analyzerBuilder.Build();


            //Sparse tet Mesh
            var watchDofs = new[]
            {
                new List<(INode node, IDofType dof)>()
                {
                    (model.NodesDictionary[nodeIdToMonitor], dofTypeToMonitor),
                    //(model.NodesDictionary[nodeIdToMonitor], StructuralDof.TranslationY),
                    //(model.NodesDictionary[nodeIdToMonitor], StructuralDof.TranslationZ),
                }
            };

            loadControlAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs[0], algebraicModel);

            return (analyzer, solver, loadControlAnalyzer);
        }

        public void SaveStateFromElements(Model model)
        {
            //elementslastConvergedDisplacements = new Dictionary<int, double[]>();
            //foreach (var elem in reader.ElementConnectivity)
            //{
            //    elementslastConvergedDisplacements[elem.Key] = ((ContinuumElement3DGrowth)model.ElementsDictionary[elem.Key]).localDisplacements.Copy();
            //}

            //elementSavedDisplacementsIsInitialized = true;
        }

        public Dictionary<int, double[][]> GetVelocities()
        {
            var nodalDofs = new StructuralDof[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
            var modelVelocities = analyzer.CreateState().StateVectors["First order derivative of solution"].Copy();
            modelVelocities.CheckForCompatibility=false;


            var interpolation = InterpolationTet4.UniqueInstance;
            var quadrature = TetrahedronQuadrature.Order1Point1;
            
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

                        nodalVelocities[3*i+i1]= dofVelocity;
                    }
                }


                var velocity = new double[nGaussPoints][];

                int gausspoinNo = 0; // anti gia loop to ntegation.Gausspoint1
                var shapeFunctionValues = interpolation.EvaluateFunctionsAt(quadrature.IntegrationPoints[gausspoinNo]);


                velocity[0]= new double[] {0,0,0};

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

            var interpolation = InterpolationTet4.UniqueInstance;
            var quadrature = TetrahedronQuadrature.Order1Point1;

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

                        bool foundVelocity = velocityNodalData.TryGetValue(nodes[i].ID, i1, out double velocityVal);
                        if (foundVelocity) { dofVelocity= velocityVal; }


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

            var interpolation = InterpolationTet4.UniqueInstance;
            var quadrature = TetrahedronQuadrature.Order1Point1;

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
