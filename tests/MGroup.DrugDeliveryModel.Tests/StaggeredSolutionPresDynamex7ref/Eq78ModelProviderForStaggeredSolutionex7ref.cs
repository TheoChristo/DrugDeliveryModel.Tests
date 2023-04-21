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
using MGroup.Solvers.AlgebraicModel;
using MGroup.LinearAlgebra.Matrices;
using MGroup.FEM.ConvectionDiffusion.Isoparametric;
using System.Reflection.PortableExecutable;
using System.Xml.Linq;
using TriangleNet.Tools;
using MGroup.MSolve.Numerics.Integration.Quadratures;
using MGroup.MSolve.Numerics.Interpolation;
using MGroup.MSolve.Discretization;

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
    public class Eq78ModelProviderForStaggeredSolutionex7ref
    {
        private double Sv;
        private double k_th_tumor;
        private double k_th_host;
        private double pv;
        private double pl;
        private double Lp;
        private double LplSvl_tumor;
        private double LplSvl_host;

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
        
        private ConvectionDiffusionDof dofTypeToMonitor;

        private ComsolMeshReader modelReader;

        public GlobalAlgebraicModel<Matrix> algebraicModel;


        public Eq78ModelProviderForStaggeredSolutionex7ref(ComsolMeshReader modelReader,
            double k_th_tumor, double k_th_host, double Lp, double Sv, double pv, double LplSvl_tumor, double LplSvl_host,
            double pl, Dictionary<int, double[]> div_vs,
            int nodeIdToMonitor, ConvectionDiffusionDof dofTypeToMonitor,
            List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> convectionDiffusionDirichletBC,
            List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> convectionDiffusionNeumannBC )
        {
            this.Sv = Sv;
            this.k_th_tumor = k_th_tumor;
            this.k_th_host = k_th_host;
            this.pv = pv;
            this.pl = pl;
            this.Lp = Lp;
            this.LplSvl_host = LplSvl_host;
            this.LplSvl_tumor = LplSvl_tumor;
            this.div_vs = div_vs;

            this.modelReader = modelReader;
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
            double capacity = 0;
            double convectionCoeff = 0;
            //double diffusion = k_th;
            

            Dictionary<int, double> diffusionCoeffs = new Dictionary<int, double>();
            Dictionary<int, double[]>ConvectionCoeffs = new Dictionary<int, double[]>(); //=> new[]  {1d, 1d, 1d};               
            Dictionary<int, double> DependentProductionCoeffs = new Dictionary<int, double>();
            Dictionary<int, double> IndependentProductionCoeffs = new Dictionary<int, double>();
            foreach (var elementConnectivity in modelReader.ElementConnectivity)
            {
                var domainId = elementConnectivity.Value.Item3;

                double LplSvl = domainId == 0 ? LplSvl_tumor : LplSvl_host;

                double dependentProductionCoeff = -(Lp * Sv + LplSvl);

                ConvectionCoeffs[elementConnectivity.Key] = new double[]
                    { convectionCoeff, convectionCoeff, convectionCoeff };
                DependentProductionCoeffs[elementConnectivity.Key] = dependentProductionCoeff;
                
                var nodes = elementConnectivity.Value.Item2;

                diffusionCoeffs[elementConnectivity.Key] = domainId == 0 ? k_th_tumor : k_th_host;

                //var independentSource = Lp * Sv * pv  - div_vs[elementConnectivity.Key][0]; 
                //var independentSource = Lp * Sv * pv; 
                var independentSource = Lp * Sv * pv + LplSvl * pl - div_vs[elementConnectivity.Key][0]; //TODO [0] is the first gauss point Make it more genreal for all guss paints
                IndependentProductionCoeffs[elementConnectivity.Key] = independentSource;
            }

            //initialize mpdel provider solution
            var modelProvider = new GenericComsol3DConvectionDiffusionProductionModelProviderDistributedSpace(modelReader);

            var model = modelProvider.CreateModelFromComsolFile(ConvectionCoeffs, diffusionCoeffs,
                DependentProductionCoeffs, IndependentProductionCoeffs, capacity);
            return model;
        }

        public void AddBoundaryConditions(Model model)
        {
            BoundaryAndInitialConditionsUtility.AssignConvectionDiffusionDirichletBCsToModel(model, convectionDiffusionDirichletBC, 1e-3);
        }

        public (IParentAnalyzer analyzer, ISolver solver, IChildAnalyzer loadcontrolAnalyzer) GetAppropriateSolverAnalyzerAndLog
        (Model model, double pseudoTimeStep, double pseudoTotalTime, int currentStep, int nIncrements)
        {
            var solverFactory = new DenseMatrixSolver.Factory() { IsMatrixPositiveDefinite = false }; //Dense Matrix Solver solves with zero matrices!
            //var solverFactory = new SkylineSolver.Factory() { FactorizationPivotTolerance = 1e-8 };
            algebraicModel = solverFactory.BuildAlgebraicModel(model);
            var solver = solverFactory.BuildSolver(algebraicModel);
            var provider = new ProblemConvectionDiffusion(model, algebraicModel);

            var loadControlAnalyzerBuilder = new LoadControlAnalyzer.Builder(algebraicModel, solver, provider, numIncrements: 1)
            {
                ResidualTolerance = 1E-8,
                MaxIterationsPerIncrement = 100,
                NumIterationsForMatrixRebuild = 1
            };
            var loadControlAnalyzer = loadControlAnalyzerBuilder.Build();

            loadControlAnalyzer.TotalDisplacementsPerIterationLog = new TotalDisplacementsPerIterationLog(new List<(INode node, IDofType dof)>()
            {(model.NodesDictionary[nodeIdToMonitor], dofTypeToMonitor)}, algebraicModel);

            var analyzerBuilder = new NewmarkDynamicAnalyzer.Builder(algebraicModel, provider, loadControlAnalyzer, timeStep: pseudoTimeStep, totalTime: pseudoTotalTime, false, currentStep: currentStep);
            analyzerBuilder.SetNewmarkParametersForConstantAcceleration();
            var analyzer = analyzerBuilder.Build();
            var watchDofs = new[]
            {
                new List<(INode node, IDofType dof)>()
                {
                    (model.NodesDictionary[nodeIdToMonitor], dofTypeToMonitor),
                }
            };
            loadControlAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs[0], algebraicModel);

            return (analyzer, solver, loadControlAnalyzer);
        }

        public List<double> RetrievePressureSolution(ISolver solver, IChildAnalyzer childAnalyzer, Model model, GlobalAlgebraicModel<Matrix> algebraicModel)
        {
            List<double> p = new List<double>();

            var u = childAnalyzer.CurrentAnalysisResult;

            var currentSolution = algebraicModel.ExtractAllResults(u);
            //var freedofs = algebraicModel.SubdomainFreeDofOrdering;

            foreach (var node in model.NodesDictionary.Values)
            {
                p.Add(currentSolution.Data[node.ID, 0]);
            }

            return p;
        }

        internal void UpdatePressureDivergenceDictionary(Dictionary<int, double[][]> pressureTensorDivergenceAtElementGaussPoints, ISolver solver, IChildAnalyzer childAnalyzer, Model model, GlobalAlgebraicModel<Matrix> algebraicModel)
        {
            var p = childAnalyzer.CurrentAnalysisResult;
            var interpolation = InterpolationTet4.UniqueInstance;
            var quadrature = TetrahedronQuadrature.Order1Point1;
            foreach (var elem in model.ElementsDictionary.Values)
            {
                var elemSoultion = algebraicModel.ExtractElementVector(p, elem);
                pressureTensorDivergenceAtElementGaussPoints[elem.ID] = UpdatePressureAndGradientsOfElement(elemSoultion, interpolation, quadrature, elem);
            }
           
        }

        private double[][] UpdatePressureAndGradientsOfElement(double[] localDisplacements, IIsoparametricInterpolation3D Interpolation,
         IQuadrature3D quadratureForMass, IElementType element)
        {
            double[][] pressureTensorDivergenceAtGaussPoints = new double[quadratureForMass.IntegrationPoints.Count][];
            IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives;
            shapeFunctionNaturalDerivatives = Interpolation.EvaluateNaturalGradientsAtGaussPoints(quadratureForMass);
            var jacobians = shapeFunctionNaturalDerivatives.Select(x => new IsoparametricJacobian3D(element.Nodes, x));
            Matrix[] jacobianInverse = jacobians.Select(x => x.InverseMatrix.Transpose()).ToArray();
            for (int gp = 0; gp < quadratureForMass.IntegrationPoints.Count; ++gp)
            {
                double[] dphi_dnatural = new double[3]; //{ dphi_dksi, dphi_dheta, dphi_dzeta}
                for (int i1 = 0; i1 < shapeFunctionNaturalDerivatives[gp].NumRows; i1++)
                {
                    dphi_dnatural[0] += shapeFunctionNaturalDerivatives[gp][i1, 0] * localDisplacements[i1];
                    dphi_dnatural[1] += shapeFunctionNaturalDerivatives[gp][i1, 1] * localDisplacements[i1];
                    dphi_dnatural[2] += shapeFunctionNaturalDerivatives[gp][i1, 2] * localDisplacements[i1];
                }

                var dphi_dnaturalMAT = Matrix.CreateFromArray(dphi_dnatural, 1, 3);

                var dphi_dcartesian = dphi_dnaturalMAT * jacobianInverse[gp].Transpose();

                pressureTensorDivergenceAtGaussPoints[gp] = new double[3] { dphi_dcartesian[0, 0], dphi_dcartesian[0, 1], dphi_dcartesian[0, 2] };

            }

            return pressureTensorDivergenceAtGaussPoints;
        }
    }
}
